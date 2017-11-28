
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

from .step import Step
from .exceptions import *
from ..genomic_io.fastq import FastqFile
from ..genomic_io.fasta import FastaFile, FastaEntry
from ..genomic_io.functions import make_fasta_from_fastq
from ..annotation.intron import get_intron_sequences

#################################################################

class IntronSequences(Step):
    ''' The first input file is the GTF file that annotates the whole genome file (the second input)
    The second input file is the genome fasta file
    The five_prime length is the number of nucleotides to be included in the five prime intron file
    '''
    def __init__(self, name, input_files, output_directory,
                 executable='', executable_arguments = '', five_prime_length = 20):
        super().__init__(name, input_files, output_directory, executable, executable_arguments)
        if len(input_files) != 2:
            raise StepError("IntronSequences object requires two input files: gtf_file and fasta_file")

        if five_prime_length < 1:
            raise StepError("Five prime length must be greater than zero. {} given".format(five_prime_length))

        self.five_prime_length = five_prime_length
        self.gtf_file          = input_files[0]
        self.fasta_file        = input_files[1]

        self.intron_sequence_file   = os.path.join(self.output_directory, "introns.fasta")
        self.five_prime_intron_file = os.path.join(self.output_directory, "introns_five_prime.fasta")

    ###################################################################
    def prepare(self):
        self.command = ""

    ##############################################################################
    def _module_run(self):
        get_intron_sequences(gtf_file              = self.gtf_file,
                             fasta_file            = self.fasta_file,
                             five_prime_slice_len  = self.five_prime_length,
                             intron_sequence_file  = self.intron_sequence_file,
                             five_prime_slice_file = self.five_prime_intron_file)

    ###############################################################################

    def post_run(self):
        if not os.path.isfile(self.intron_sequence_file):
            self.error_messages.append("Couldn't find the intron sequence file {}".\
                                       format(self.intron_sequence_file))
        if not os.path.isfile(self.five_prime_intron_file):
            self.error_messages.append("Couldn't find the five prime sequence file {}".\
                                       format(self.five_prime_intron_file))

        if len(self.error_messages) > 0:
           subprocess.call('touch ' + self.failure_file , shell=True )
        else:
           subprocess.call('touch ' + self.success_file , shell=True )




