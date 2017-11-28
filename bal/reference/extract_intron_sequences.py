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
from ..core.step import Step
from ..core.exceptions import *
from ..genomic_io.fasta import FastaFile, FastaEntry
from ..annotation.intron import get_intron_five_prime_sequences

#################################################################

class ExtractIntronSequences(Step):
    '''Adapted this from ref_from_reads module - only need current intron fasta file'''
    def __init__(self, name, input_files, output_directory, five_prime_length):
        '''
        N is the number of nucleotides to be taken from the five prime end
        '''

        super().__init__(name, list(), output_directory)

        self.gtf_file           = input_files[0]
        self.fasta_file         = input_files[1]
        self.five_prime_length  = five_prime_length
        self.output_file        = os.path.join(self.output_directory, "five_prime_introns.fa")
        self.intron_bed_file    = os.path.join(self.output_directory, "five_prime_introns.bed")

    ##############################################################

    def prepare(self):
        pass
            
    ###############################################################

    def _module_run(self):
        get_intron_five_prime_sequences(self.gtf_file, self.fasta_file, self.five_prime_length,
                                    self.output_file, self.intron_bed_file)


    ################################################################

    def post_run(self):
        missing_files = list()
        required_files = [self.output_file, self.intron_bed_file]

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


########################################################################
