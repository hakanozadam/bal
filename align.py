#!/bin/env python3

# AUTHORS:
#        Hakan Ozadam
#        Rachel Brown
#
#        Moore Laboratory
#        UMASS Medical School / HHMI
#        RNA Therapeutics Institute
#        Albert Sherman Center, ASC4-1009
#        368 Plantation Street
#        Worcester, MA 01605
#        USA
#
#################################################################


import os
import subprocess
from functools import wraps
from time import gmtime, strftime , localtime
import datetime

from bal.core.prepare import get_arguments, get_executables, get_bowtie2_5p_arguments
from bal.core.hisat import Hisat
from bal.core.engines import SequentialEngine
from bal.core.bowtie2 import Bowtie2
from bal.core.trim_reads import TrimReads
from bal.core.candidate_branchpoints import CandidateBranchpoints
from bal.core.make_bp_reference import BpReference
from bal.settings import *
from bal.core.bowtie2_align_trimmed import Bowtie2AlignTrimmed
from bal.core.bp_bed_from_sam import BpBedFromSam
from bal.genomic_io.functions import determine_5p_sequence_length

##################################################################

def time_this(func):
    @wraps(func)
    def wrapper(*args, **kwargs):
        start_time = datetime.datetime.now()
        result = func(*args, **kwargs)
        finish_time = datetime.datetime.now()
        print('Running time is ' +\
            str( finish_time - start_time ))
        return result
    return wrapper

##################################################################

def prep():
    try:
        import pysam
    except ImportError:
        print("Couldn't find the module pysam. "
              "Install the python3 package pysam before running bal.")
        exit(1)

##################################################################


##################################################################
@time_this
def main():
    prep()
    bal_directory    = os.path.dirname(os.path.realpath(__file__))
    executables      = get_executables(bal_directory)
    arguments        = get_arguments()
    single_end       = True if arguments.U else False
    run_success_file = os.path.join(arguments.o, ".success")

    five_prime_intron_fasta_file = os.path.join(arguments.x, settings['five_prime_intron_directory'],
                                                settings['five_prime_intron_name'] + '.fa')

    five_prime_sequence_length = determine_5p_sequence_length(five_prime_intron_fasta_file)

    print("five_prime_sequence_length = " , five_prime_sequence_length)

    bowtie2_5p_arguments = get_bowtie2_5p_arguments(five_prime_sequence_length,
                                                            arguments.five_p_allowed_mismatch)

    print('five_prime_bowtie2_arguments =', bowtie2_5p_arguments)

    #exit(0) # delete this line later

    if os.path.isfile(run_success_file):
        os.remove(run_success_file)

    bowtie2_mate_1_strand_parameter = " --norc "
    bowtie2_mate_2_strand_parameter = " --nofw "

    if arguments.rna_strandness == 'R':
        bowtie2_mate_1_strand_parameter = " --nofw "
        bowtie2_mate_2_strand_parameter = " --norc "

    job_runner = SequentialEngine()

    ###############################################################
    #### HISAT STEP ###############################################
    ###############################################################

    hisat_job = Hisat('Hisat' , [], os.path.join(arguments.o, 'hisat'), executable = executables['hisat'],
                      executable_arguments = ' -N 1 ',
                      genome_reference = os.path.join(arguments.x, settings['hisat_genome_directory'],
                                                      settings['hisat_genome_name']),
                      known_splice_sites_file = os.path.join(arguments.x, settings['hisat_genome_directory'],
                                                             settings['hisat_known_splice_sites_file']),
                      single_end_fastq = arguments.U,
                      paired_end_fastqs = (arguments.mate_1 , arguments.mate_2), p = arguments.p )
    job_runner[hisat_job.name] = hisat_job

    ###################################################################
    ###### Bowtie2 Jobs ###############################################
    ###### Align Reads against 5p end of the introns ##################
    ###################################################################


    bowtie2_5p_refrence_base = os.path.join(arguments.x, settings['five_prime_intron_reference_directory'],
                                            settings['five_prime_intron_reference_name'])
    bowtie2_gene_reference_base = os.path.join( arguments.x, settings['bowtie2_gene_directory'] ,
                                                settings['bowtie2_gene_sub_directory'] )

    if single_end:
        bowtie2_5p_single_end_job = Bowtie2( name = "bowtie2_5p_single_end",
                                         input_files = list(),
                                         output_directory = os.path.join(arguments.o, 'bowtie2_5p_single_end'),
                                         executable = executables['bowtie2'],
                                         executable_arguments = bowtie2_5p_arguments + " " +\
                                                                bowtie2_mate_1_strand_parameter  ,\
                                         genome_reference  = bowtie2_5p_refrence_base ,
                                         single_end_fastq = os.path.join(arguments.o, 'hisat', 'unaligned.fastq'),
                                         paired_end_fastqs = list() ,
                                         p = arguments.p)
        job_runner[bowtie2_5p_single_end_job.name] = bowtie2_5p_single_end_job
    else:
        bowtie2_5p_mate_1_job = Bowtie2( name = "bowtie2_5p_mate_1",
                                         input_files = list(),
                                         output_directory = os.path.join(arguments.o,
                                                                         'bowtie2_5p' ,'bowtie2_5p_mate_1'),
                                         executable = executables['bowtie2'],
                                         executable_arguments = bowtie2_5p_arguments + " " +\
                                                                bowtie2_mate_1_strand_parameter ,\
                                         genome_reference  = bowtie2_5p_refrence_base ,
                                         single_end_fastq = os.path.join(arguments.o,
                                                                         'hisat', 'unaligned.1.fastq'),
                                         paired_end_fastqs = list(),
                                         p = arguments.p)
        job_runner[bowtie2_5p_mate_1_job.name] = bowtie2_5p_mate_1_job

        bowtie2_5p_mate_2_job = Bowtie2( name = "bowtie2_5p_mate_2",
                                         input_files = list(),
                                         output_directory = os.path.join(arguments.o,
                                                                         'bowtie2_5p', 'bowtie2_5p_mate_2'),
                                         executable = executables['bowtie2'],
                                         executable_arguments = bowtie2_5p_arguments + " " +\
                                                                bowtie2_mate_2_strand_parameter ,\
                                         genome_reference  = bowtie2_5p_refrence_base ,
                                         single_end_fastq = os.path.join(arguments.o,
                                                                         'hisat', 'unaligned.2.fastq'),
                                         paired_end_fastqs = list(),
                                         p = arguments.p)
        job_runner[bowtie2_5p_mate_2_job.name] = bowtie2_5p_mate_2_job

    ########################################################################
    #################     Trim Reads     ###################################
    ########################################################################

    if single_end:
        trim_reads_single_end_job = TrimReads(name = "trim_reads",
                                              input_files = (os.path.join(bowtie2_5p_single_end_job.output_directory,
                                                                         "alignment.sam"),),
                                              output_directory = os.path.join(arguments.o, "trimmed_reads"),
                                              executable = "", executable_arguments = "",
                                              length_threshold = settings["trim_length_threshold"] )
        job_runner[trim_reads_single_end_job.name] = trim_reads_single_end_job
    else:
        trim_reads_single_mate_1_job = TrimReads(name = "trim_reads_mate_1",
                                              input_files = (os.path.join(bowtie2_5p_mate_1_job.output_directory,
                                                                         "alignment.sam"),),
                                              output_directory = os.path.join(arguments.o, "trimmed_reads", "mate_1"),
                                              executable = "", executable_arguments = "",
                                              length_threshold = settings["trim_length_threshold"] )
        job_runner[trim_reads_single_mate_1_job.name] = trim_reads_single_mate_1_job

        trim_reads_single_mate_2_job = TrimReads(name = "trim_reads_mate_2",
                                              input_files = (os.path.join(bowtie2_5p_mate_2_job.output_directory,
                                                                         "alignment.sam"),),
                                              output_directory = os.path.join(arguments.o, "trimmed_reads", "mate_2"),
                                              executable = "", executable_arguments = "",
                                              length_threshold = settings["trim_length_threshold"] )
        job_runner[trim_reads_single_mate_2_job.name] = trim_reads_single_mate_2_job



    ##############################################################################
    ################# Align Trimmed Reads   ######################################
    ##############################################################################

    if single_end:
        bowtie2_trimmed_alignment_single_end_job = \
            Bowtie2AlignTrimmed(name                     = "trimmed_alignment_single_end",
                                input_files              = list(),
                                output_directory         = os.path.join(arguments.o, 'alignment_of_trimmed'),
                                executable               = executables['bowtie2'],
                                executable_arguments     = settings['bowtie2_trimmed_read_parameters']  ,
                                bt2_reference_directory  = bowtie2_gene_reference_base,
                                fastq_directory          = trim_reads_single_end_job.fastq_dir
                                )
        job_runner[bowtie2_trimmed_alignment_single_end_job.name] = bowtie2_trimmed_alignment_single_end_job
    else:
        bowtie2_trimmed_alignment_mate_1_job = \
            Bowtie2AlignTrimmed(name                     = "trimmed_alignment_mate_1",
                                input_files              = list(),
                                output_directory         = os.path.join(arguments.o, 'alignment_of_trimmed', 'mate_1'),
                                executable               = executables['bowtie2'],
                                executable_arguments     = settings['bowtie2_trimmed_read_parameters']  ,
                                bt2_reference_directory  = bowtie2_gene_reference_base,
                                fastq_directory          = trim_reads_single_mate_1_job.fastq_dir
                                )
        job_runner[bowtie2_trimmed_alignment_mate_1_job.name] = bowtie2_trimmed_alignment_mate_1_job

        bowtie2_trimmed_alignment_mate_2_job = \
            Bowtie2AlignTrimmed(name                     = "trimmed_alignment_mate_2",
                                input_files              = list(),
                                output_directory         = os.path.join(arguments.o, 'alignment_of_trimmed', 'mate_2'),
                                executable               = executables['bowtie2'],
                                executable_arguments     = settings['bowtie2_trimmed_read_parameters']  ,
                                bt2_reference_directory  = bowtie2_gene_reference_base,
                                fastq_directory          = trim_reads_single_mate_2_job.fastq_dir
                                )
        job_runner[bowtie2_trimmed_alignment_mate_2_job.name] = bowtie2_trimmed_alignment_mate_2_job

    ###############################################################################
    ######## Determine BP Candidates  #############################################

    if single_end:
        candidate_branchpoints_singe_end_job = \
            CandidateBranchpoints('candidate_branchpoints_single_end',\
                                  input_files = ( bowtie2_trimmed_alignment_single_end_job.alignment_directory , ),
                                  output_directory = os.path.join(arguments.o, 'candidate_branchpoints'))
        job_runner[candidate_branchpoints_singe_end_job.name] = candidate_branchpoints_singe_end_job
    else:
        candidate_branchpoints_paired_end_job = \
            CandidateBranchpoints('candidate_branchpoints_paired_end',\
                                  input_files = (bowtie2_trimmed_alignment_mate_1_job.alignment_directory ,
                                                 bowtie2_trimmed_alignment_mate_2_job.alignment_directory) ,
                                  output_directory = os.path.join(arguments.o,
                                                                  'candidate_branchpoints'))
        job_runner[candidate_branchpoints_paired_end_job.name] = candidate_branchpoints_paired_end_job

    job_runner.run()
    job_runner = SequentialEngine()

    if arguments.hpc == True:
        subprocess.call('touch ' + run_success_file , shell=True )
        return 0

    ####################################################################################
    ######### Make Branchpoint Reference  ##############################################

    make_bp_job = BpReference('branchpoint_reference',
                              input_files          = [ os.path.join(arguments.x, settings['genome_fasta_file']) ,
                                                       os.path.join(arguments.o, 'candidate_branchpoints',
                                                                    settings['bp_candidate_file_name'])],
                              output_directory     = os.path.join(arguments.o, 'branchpoint_reference'),
                              executable           = executables['bowtie2-build'],
                              number_of_nucleotides = settings['number_of_nucleotides'],
                              executable_arguments = ""
                              )
    job_runner[make_bp_job.name] = make_bp_job

    if single_end:
        bp_candidate_bed_file = candidate_branchpoints_singe_end_job.bed_file
    else:
        bp_candidate_bed_file = candidate_branchpoints_paired_end_job.bed_file

    bp_candidate_bed_file_line_count = 0
    with open(bp_candidate_bed_file, 'r') as input_bed:
        for line in input_bed:
            if len(line.strip()) > 0:
                bp_candidate_bed_file_line_count += 1

    if bp_candidate_bed_file_line_count == 0:
        with open(run_success_file, 'w') as run_success_out:
            print("No branchpoints found. Exiting.")
            print("No branchpoints found", file = run_success_out)
        return

    ####################################################################################
    ######### Align Reads Against Branchpoints #########################################

    if single_end:
        align_bp_single_end_job = Bowtie2( name = "align_bp_single_end",
                                           input_files = list(),
                                           output_directory = os.path.join(arguments.o, 'align_bp_single_end'),
                                           executable = executables['bowtie2'],
                                           executable_arguments = settings["align_bp_bt2_arguments"] +
                                                                  bowtie2_mate_1_strand_parameter  ,
                                           genome_reference  = make_bp_job.reference_base ,
                                           single_end_fastq = os.path.join(arguments.o, 'hisat', 'unaligned.fastq'),
                                           paired_end_fastqs = list(),
                                           p = arguments.p )
        job_runner[align_bp_single_end_job.name] = align_bp_single_end_job
    else:
        align_bp_mate_1_job = Bowtie2( name = "align_bp_mate_1",
                                           input_files = list(),
                                           output_directory = os.path.join(arguments.o, 'align_bp_mate_1'),
                                           executable = executables['bowtie2'],
                                           executable_arguments = settings["align_bp_bt2_arguments"] +
                                                                  bowtie2_mate_1_strand_parameter ,\
                                           genome_reference  = make_bp_job.reference_base ,
                                           single_end_fastq = os.path.join(arguments.o, 'hisat', 'unaligned.1.fastq'),
                                           paired_end_fastqs = list(),
                                           p = arguments.p)
        job_runner[align_bp_mate_1_job.name] = align_bp_mate_1_job

        align_bp_mate_2_job = Bowtie2( name = "align_bp_mate_2",
                                           input_files = list(),
                                           output_directory = os.path.join(arguments.o, 'align_bp_mate_2'),
                                           executable = executables['bowtie2'],
                                           executable_arguments = settings["align_bp_bt2_arguments"] +
                                                                  bowtie2_mate_2_strand_parameter ,\
                                           genome_reference  = make_bp_job.reference_base ,
                                           single_end_fastq = os.path.join(arguments.o, 'hisat', 'unaligned.2.fastq'),
                                           paired_end_fastqs = list(),
                                           p = arguments.p )
        job_runner[align_bp_mate_2_job.name] = align_bp_mate_2_job

    job_runner.run()
    job_runner = SequentialEngine()

    #############################################################################################
    ######### BP From Sam Ref  ##################################################################
    ######### Get the branchpoints that have a hit from the reads not aligned to the genome
    #############################################################################################

    bp_from_sam_ref_job_input_files = list()

    if single_end:
        bp_from_sam_ref_job_input_files.append( os.path.join(align_bp_single_end_job.output_directory,
                                                             "alignment.sam") )
    else:
        bp_from_sam_ref_job_input_files.append( os.path.join(align_bp_mate_1_job.output_directory,
                                                             "alignment.sam") )
        bp_from_sam_ref_job_input_files.append( os.path.join(align_bp_mate_2_job.output_directory,
                                                             "alignment.sam") )

    bp_from_sam_ref_job = BpBedFromSam(name             = 'BP_From_Sam_Ref',
                                       input_files      = bp_from_sam_ref_job_input_files,
                                       output_directory = os.path.join(arguments.o, 'bp_candidates_genome_level_1'))

    job_runner[bp_from_sam_ref_job.name] = bp_from_sam_ref_job

    ####################################################################################

    job_runner.run()
    subprocess.call('touch ' + run_success_file , shell=True )
    return 0



if __name__ == '__main__':
    main()