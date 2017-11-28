#!/bin/env python3

# AUTHORS:
#        Hakan Ozadam
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

from bal.umass_cluster.prepare import get_arguments, get_executables,\
                                arrange_input_files
from bal.core.engines import SequentialEngine
from bal.settings import *
from bal.umass_cluster.lsf.PartitionFastqJobGroup import PartitionFastqJobGroup
from bal.umass_cluster.lsf.BalJobGroup import BalJobGroup
from bal.umass_cluster.lsf.MergeBedJobGroup import MergeBedJobGroup
from bal.umass_cluster.lsf.BranchpointRefJobGroup import BranchpointRefJobGroup
from bal.umass_cluster.lsf.AlignBpJobGroup import AlignBpJobGroup
from bal.umass_cluster.lsf.ExtractBpFromSamRefJobGroup import ExtractBpFromSamRefJobGroup

###################################################################

###################################################################

def main():
    bal_directory = os.path.dirname(os.path.realpath(__file__))
    executables   = get_executables(bal_directory)

    arguments           = get_arguments()
    input_directory     = os.path.abspath(arguments.i)
    output_directory    = os.path.abspath(arguments.o)
    reference_directory = os.path.abspath(arguments.x)
    alignment_mode      = arguments.m
    reads_per_file      = arguments.n

    cluster_output_directory     = os.path.join(output_directory, ".lsf_out")
    partitioned_input_directory  = os.path.join(output_directory, "partitioned_input")
    partitioned_output_directory = os.path.join(output_directory, "partitioned_output")
    merged_bed_directory         = os.path.join(output_directory, "merged_bed_files")
    candidate_bp_directory       = os.path.join(output_directory, "candidate_bp_directory")
    bp_reference_directory       = os.path.join(output_directory, "bp_bt2_reference")
    bp_alignments_directory      = os.path.join(output_directory, "bp_alignments")
    main_output_directory        = os.path.join(output_directory, "output")
    genome_fasta_file            = os.path.join(reference_directory, "genome.fa")

    candidate_bp_genome_level_1_directory     = os.path.join(output_directory, "candidate_bp_genome_level_1")

    #job_runner = SequentialEngine()
    #print("Umass pipeline is working!")
    arranged_input_files = arrange_input_files(input_directory, alignment_mode )

    ####################################
    ###   Partition Input Files    #####
    ####################################
    jobGroup = PartitionFastqJobGroup(
                 input_directory  = input_directory ,
                 output_directory = partitioned_input_directory ,
                 alignment_mode   = alignment_mode ,
                 run_time         = 3 ,
                 reads_per_file   = reads_per_file ,
                 memory           = 4096 ,
                 executable       = executables['partition_fastq']
                  )
    jobGroup.run_in_main()

    #####################################
    ###   BAL on Partitioned Data   #####
    #####################################

    bal_main_arguments = arguments.a + " --hpc "
    if arguments.rna_strandness == "R":
        bal_main_arguments += " --rna-strandness R "
    else:
        bal_main_arguments += " --rna-strandness F "

    bal_threads = arguments.p
    if bal_threads < 10:
        bal_threads = 10

    jobGroup = BalJobGroup(
                 input_directory  = partitioned_input_directory ,
                 output_directory = partitioned_output_directory ,
                 arguments        = bal_main_arguments ,
                 reference        = reference_directory ,
                 alignment_mode   = alignment_mode ,
                 threads          = bal_threads ,
                 run_time         = arguments.t ,
                 memory           = arguments.memory ,
                 executable       = executables['bal']
                  )
    jobGroup.run_in_main()

    #############################################
    #### Merge Piece Bed Files   ################
    #############################################
    # Merge the bed files of each library
    # Bed files come in pieces as we partition the lib input file
    # to speed up the alignment
    # First merge the bed files for each library

    jobGroup = MergeBedJobGroup(
                 input_directory  = partitioned_output_directory ,
                 output_directory = merged_bed_directory ,
                 executable       = executables['merge_bed_files']
                  )
    jobGroup.run_in_main()

    ############################################################
    #### Merge Library Candidate BP Bed Files   ################
    ############################################################
    # Now we have the merged bed files from each library.
    # So we merge ll library files to get the whole list of
    # branchpoints coming from all libraries.

    merge_bed_jobGroup = MergeBedJobGroup(
                 input_directory  = partitioned_output_directory ,
                 output_directory = candidate_bp_directory ,
                 single_file_list = os.path.join(merged_bed_directory, 'candidate_bp_files_list.txt'),
                 executable       = executables['merge_bed_files']
                  )
    merge_bed_jobGroup.run_in_main()

    #############################################################
    ##### Make Bowtie2 Reference  ###############################
    #############################################################

    bp_ref_jobGroup = BranchpointRefJobGroup(
        input_directory       = os.path.dirname( os.path.abspath(merge_bed_jobGroup.bp_candidates_bed_file) ), # Note that input direcory has no importance in this class
        output_directory      = bp_reference_directory ,
        bed_file              = merge_bed_jobGroup.bp_candidates_bed_file ,
        fasta_file            = genome_fasta_file,
        executable            = executables['make_bp_ref'],
        number_of_nucleotides = settings['number_of_nucleotides'],
        name                  = "Ref_Branchpoint" ,
        run_time              = 3, # hours
        memory                = 8000)
    bp_ref_jobGroup.run_in_main()

    #############################################################
    ##### Align Reads Against The Reference  ####################
    #############################################################

    align_bp_jobGroup = AlignBpJobGroup(
        input_directory = partitioned_output_directory,
        output_directory = bp_alignments_directory,
        arguments = ' ',
        reference = bp_reference_directory,
        alignment_mode = alignment_mode,
        threads = 8,
        name = "bp_align" ,
        run_time = 3,
        executable = executables['align_bp'],
        memory = 8192 ,
        rna_strandness = arguments.rna_strandness)

    align_bp_jobGroup.run_in_main()

    ###############################################################
    ###### Extract BP from SAM Ref   ##############################
    ###############################################################

    extract_bp_from_sam_ref_jobGroup = ExtractBpFromSamRefJobGroup(
        input_directory  = bp_alignments_directory,
        output_directory = candidate_bp_genome_level_1_directory,
        individual_lib_directory = align_bp_jobGroup.individual_sam_files_directory,
        executable = executables['bp_bed_from_sam_ref'] )

    extract_bp_from_sam_ref_jobGroup.run_in_main()

#####################################################################

if __name__ == '__main__':
    main()