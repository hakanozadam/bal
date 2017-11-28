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

from shutil import which
import os
import subprocess

from .step import Step
from .exceptions import *
# from ..genomic_io.fastq import FastqFile
from ..genomic_io.fasta import FastaFile, FastaEntry
# from ..genomic_io.functions import make_fasta_from_fastq

#################################################################

class ReferenceFromIntrons(Step):
    '''Adapted this from ref_from_reads module - only need current intron fasta file'''
    def __init__(self, name, input_files, output_directory, executable, executable_arguments=""):
        super().__init__(name, list(), output_directory, executable, executable_arguments)
        if len(input_files) != 1:
            raise( StepError("In ReferenceFromReads Step, you need to provide one intron fastq files."
                             "The given file list was " + "\n".join(input_files)) )

        self.pre_input_files        = input_files # we postpone file existence check till prepare step
        self.reference_prefix = "reference"
        self.reference_base   = os.path.join( self.output_directory, self.reference_prefix )
#        self.paired_reference_basis  = ( os.path.join( self.output_directory, self.reference_prefix + "_1" ),
#                                          os.path.join( self.output_directory, self.reference_prefix + "_2") )

    ##############################################################

    def prepare(self):
        '''
        Bowtie2-build needs fasta files for input so convert fastq files to fasta
        Also for paired-end case, take the reverse complement of the second mate before making a reference
        The output (as a side effect) of this function is (are) bowtie2 reference(s).
        '''

        self.input_files        = self.pre_input_files
#        self.mate_1_fastq_file  = self.input_files[0]
#       self.mate_2_fastq_file  = os.path.join(self.output_directory, "mate_2.fastq")
        self.intron_fasta_file  = os.path.join(self.output_directory, "intron.fasta")
#        self.paired_fasta_files = ( os.path.join(self.output_directory, "mate_1.fasta") ,
#                                        os.path.join(self.output_directory, "mate_2.fasta"))


        self.intron_fasta_file = self.input_files[0]
        self.command = " ".join((self.executable, self.intron_fasta_file, self.reference_base))
        print(self.command)
            
    ###############################################################

    def post_run(self):
        output_extensions = ('.1.bt2', '.2.bt2' , '.rev.1.bt2', '.rev.2.bt2')
        missing_files = list()

        if len(self.input_files) == 1:
            required_files = map( lambda x: self.single_reference_base + x , output_extensions )
        else:
            required_files = list(map( lambda x: self.paired_reference_basis[0] + x , output_extensions )) +\
                             list(map( lambda x: self.paired_reference_basis[1] + x , output_extensions ) )

        for file in required_files:
            if not os.path.isfile(file):
                missing_files.append(file)

        if len(missing_files) > 0:
            missing_file_message = "The following files are missing:\n" +\
                "\n".join(missing_files)
            self.error_messages.append(missing_file_message)
            subprocess.call('touch ' + self.failure_file , shell=True )
        else:
            subprocess.call('touch ' + self.success_file , shell=True )
