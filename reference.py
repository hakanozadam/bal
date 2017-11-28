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

from bal.reference.prepare import get_arguments, get_executables
from bal.reference.bowtie2_reference import Bowtie2Reference
from bal.reference.hisat_reference import HISATReference
from bal.reference.extract_intron_sequences import ExtractIntronSequences
from bal.reference.extract_gene_sequences import ExtractGeneSequences
from bal.reference.bowtie2_gene_reference import Bowtie2GeneReference
from bal.core.engines import SequentialEngine
from bal.genomic_io.gtf import check_gtf_validity
from bal.settings import *

##################################################################

###################################################################

def main():
    bal_directory = os.path.dirname(os.path.realpath(__file__))
    executables = get_executables(bal_directory)
    arguments = get_arguments()

    job_runner = SequentialEngine()

    bowtie2_intron_sequence_directory   = os.path.join(arguments.o, settings['five_prime_intron_directory'] )
    bowtie2_intron_reference_directory  = os.path.join(arguments.o, settings['five_prime_intron_reference_directory'])
    bowtie2_genome_directory            = os.path.join(arguments.o, settings['bowtie2_genome_directory'])
    gene_sequence_directory             = os.path.join(arguments.o, settings['gene_sequence_directory'])
    bowtie2_gene_directory              = os.path.join(arguments.o, settings['bowtie2_gene_directory'])
    HISAT_genome_directory              = os.path.join(arguments.o, settings['hisat_genome_directory'])
    genome_fasta_file                   = os.path.join(arguments.o, settings['genome_fasta_file'])

    gtf_file = os.path.abspath(arguments.g)
    gtf_result = check_gtf_validity(gtf_file)
    if gtf_result:
        print("Error! The given gtf file %s can not be used with this software"%gtf_file)
        print(gtf_result)
        exit(1)



    os.makedirs(arguments.o, exist_ok=True)
    print("Copying the genome fasta file ", genome_fasta_file)
    command = " ".join( ("cp", arguments.f, genome_fasta_file) )
    p = subprocess.Popen([command], stdout=subprocess.PIPE,
                                                  stderr=subprocess.PIPE,
                                                  shell = True)
    std_out , std_err = p.communicate()
    cp_return_code    = p.returncode

    if cp_return_code:
        raise IOError("The copying of the file", genome_fasta_file, "failed.")

    print("Preparing references...")

    extract_job = ExtractIntronSequences('Five_Prime_Intron_Sequence' , [arguments.g, arguments.f],
                                         bowtie2_intron_sequence_directory,
                                         five_prime_length = arguments.N )
    job_runner[extract_job.name] = extract_job

    ####################

    bowtie2_intron_ref_job = Bowtie2Reference('Five_Prime_Intron_Reference' ,
                                            [ os.path.join(bowtie2_intron_sequence_directory, "five_prime_introns.fa") ] ,
                                            bowtie2_intron_reference_directory , executables["bowtie2-build"] , "",
                                            genome_reference_name = "intron")
    job_runner[bowtie2_intron_ref_job.name] = bowtie2_intron_ref_job

    #######################################################################

    extract_gene_sequences_job = ExtractGeneSequences(
                 name = "Extract_Gene_Sequences",
                 input_files  = [arguments.g, arguments.f] ,
                 output_directory = gene_sequence_directory
                  )

    job_runner[extract_gene_sequences_job.name] = extract_gene_sequences_job

    #######################################################################

    bowtie2_gene_reference_job = Bowtie2GeneReference(
                 name             = "Bowtie2_Gene_Reference",
                 input_files      = [ extract_gene_sequences_job.gene_directory ] ,
                 output_directory = bowtie2_gene_directory ,
                 executable       = executables['bowtie2-build']
                  )

    job_runner[bowtie2_gene_reference_job.name] = bowtie2_gene_reference_job

    ######################

    hisat_genome_ref_job = HISATReference(executables['hisat_extract_splice_sites'], 'HISAT_Genome_Reference' ,
                                            [ arguments.f , arguments.g] ,
                                            HISAT_genome_directory , executables["hisat-build"] , "",
                                            genome_reference_name = "genome")
    job_runner[hisat_genome_ref_job.name] = hisat_genome_ref_job

    job_runner.run()
    print("You can run bal aligner using the reference parameter as -x ", arguments.o)

#####################################################################

if __name__ == '__main__':
    main()