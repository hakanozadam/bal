#
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
from collections import OrderedDict, defaultdict

from .step import Step
from .exceptions import *
from ..genomic_io.fastq import FastqFile
from ..genomic_io.fasta import FastaFile, FastaEntry
from ..genomic_io.functions import make_fasta_from_fastq
from ..annotation.intron import get_intron_sequences
from ..settings import *
import pysam

#################################################################

class BpReference(Step):
    '''
    To be completed
    '''
    def __init__(self, name, input_files, output_directory,
                 executable='', executable_arguments = '' , number_of_nucleotides = 150):
        super().__init__(name, [], output_directory, executable, executable_arguments)

        self.genome_fasta_file     = input_files[0]
        self.bp_bed_file           = input_files[1]
        self.bp_fasta_file         = os.path.join(self.output_directory, "bp_sequences.fa")
        self.number_of_nucleotides = number_of_nucleotides
        self.reference_base        = os.path.join(self.output_directory, settings['bp_reference_base'])
        
    ###################################################################
    def prepare(self):

        self.get_bp_sequences()
        if not os.path.isfile(self.bp_fasta_file):
            raise StepError("There was a problem in getting the bp sequences."
                                       "BP reference file %s doesn't exist."%self.bp_fasta_file)

        self.command = " ".join( [ self.executable, self.executable_arguments,
                                   self.bp_fasta_file, self.reference_base ] )

    ##############################################################################

    def get_bp_sequences(self):
        with open(self.bp_bed_file, "r") as bed_input,\
             FastaFile(self.genome_fasta_file) as genome_input,\
             open(self.bp_fasta_file, "w") as bp_output:

             # first read the bed file into a dict grouped by chromosome
             bps_by_chr = defaultdict(list)
             for bp_entry in bed_input:
                 bp_contents = bp_entry.rstrip().split("\t")
                 bp_chr      = bp_contents[0]
                 bps_by_chr[bp_chr].append(bp_entry)
                 
             # Then go through the fasta file and get the sequences
             for chr_entry in genome_input:
                 this_chr = chr_entry.header
                 for branchpoint in bps_by_chr[this_chr]:
                     bp_contents = branchpoint.rstrip().split("\t")
                     bp_location = int(bp_contents[1])
                     bp_header_contents  = bp_contents[3].split(settings['field_separator'])
                     five_prime_location = int(bp_header_contents[3])
                     this_sequence = ''
                     if bp_contents[5] == '+':
                         bp_fragment_start = bp_location - self.number_of_nucleotides
                         if bp_fragment_start < 0 :
                             bp_fragment_start = 0
                         this_sequence = chr_entry.sequence[ bp_fragment_start : bp_location + 1 ] +\
                             chr_entry.sequence[ five_prime_location : five_prime_location +\
                                                                       self.number_of_nucleotides ]
                     elif bp_contents[5] == '-':
                         bp_fragment_raw = chr_entry.sequence[ bp_location :\
                                                               bp_location + self.number_of_nucleotides + 1 ]
                         five_p_fragment_start = five_prime_location - self.number_of_nucleotides
                         five_p_fragment_raw   = chr_entry.sequence[ five_p_fragment_start + 1 :\
                                                                     five_prime_location + 1 ]
                         bp_fragment_raw_fasta = FastaEntry('bp' , bp_fragment_raw)
                         bp_fragment_raw_fasta.reverse_complement()
                         five_p_fragment_raw_fasta = FastaEntry('five_p' , five_p_fragment_raw )
                         five_p_fragment_raw_fasta.reverse_complement()
                         this_sequence = bp_fragment_raw_fasta.sequence + five_p_fragment_raw_fasta.sequence
                     else:
                         raise(StepError("Invalid strand type:", bp_contents[5]))
                     this_bp_sequence_entry = FastaEntry(bp_contents[3], this_sequence)
                     print(this_bp_sequence_entry, file = bp_output)

    ###############################################################################

    ###############################################################################

    ###############################################################################

    def post_run(self):
        missing_references = list()
        suffixes = ('.1.bt2', '.2.bt2', '.3.bt2', '.4.bt2', '.rev.1.bt2',  '.rev.2.bt2')
        error_messages = list()

        for suffix in suffixes:
          if (not os.path.isfile(self.reference_base + suffix) ) :
              missing_references.append("Couldn't find the bowtie2 reference: " + self.reference_base + suffix)
        if len(missing_references) > 0:
           error_messages.append("Couldn't find the following bowtie2 reference(s):\n" +\
                  "\n".join(missing_references))

        if len(error_messages) > 0:
           subprocess.call('touch ' + self.failure_file , shell=True )
        else:
           subprocess.call('touch ' + self.success_file , shell=True )

        self.error_messages = error_messages

