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
import pysam
from collections import defaultdict
from .step import Step
from .exceptions import *
from ..genomic_io.fastq import FastqFile, FastqEntry
from ..genomic_io.fasta import FastaFile, FastaEntry
from ..settings import *

#################################################################

class TrimReads(Step):
    '''Complete this part !!!!!!!'''
    def __init__(self, name, input_files, output_directory, executable, executable_arguments = "",
                 length_threshold = 18 ):
        super().__init__(name, list(), output_directory, executable, executable_arguments)
        self.sam_file         = input_files[0]
        self.length_threshold = length_threshold
        self.fastq_dir     = os.path.join(self.output_directory, "by_genes_fastq")
        os.makedirs(self.fastq_dir, exist_ok=True)

    ##############################################################

    def prepare(self):
        pass

    ##############################################################

    def _module_run(self):
        samfile = pysam.AlignmentFile(self.sam_file, "r")
        # Keys are genes and each value is a list of fastq entries in string format
        reads_by_gene = defaultdict(list)

        # get the trimmed reads separated by genes
        for read in samfile.fetch():
            #Make sure that there are at least soft clip and the rest part
            if len(read.cigar) < 2 or read.cigar[0][0] != 4 :
                continue
            if read.cigar[0][1] < self.length_threshold:
                continue
            # get the soft clipped part only
            intron_info_contents = samfile.getrname(read.tid).split(settings['field_separator'])
            this_gene = intron_info_contents[2]
            this_header = read.query_name + settings['fastq_field_separator'] + samfile.getrname(read.tid)
            soft_clipped_part_sequence = read.seq[0:read.cigar[0][1]]
            soft_clipped_part_quality  = read.qual[0:read.cigar[0][1]]
            this_fastq_entry = FastqEntry(this_header, soft_clipped_part_sequence, "+",
                                          soft_clipped_part_quality)
            reads_by_gene[this_gene].append( str(this_fastq_entry) )

        missing_files = list()
        for key in reads_by_gene.keys():
            this_fastq_file = os.path.join(self.fastq_dir, key + ".fastq")
            with open(this_fastq_file, 'w') as fastq_output:
                for this_read in reads_by_gene[key]:
                    print(this_read, file = fastq_output)
            if not os.path.isfile(this_fastq_file):
                missing_files.append(this_fastq_file)

        self.error_messages = list()
        if len(missing_files)> 0 :
            self.error_messages = ['The following files are missing:'] + missing_files


    ###############################################################

    def post_run(self):
        if len(self.error_messages) > 0:
           print("\n".join(self.error_messages))
           subprocess.call('touch ' + self.failure_file , shell=True )
        else:
           subprocess.call('touch ' + self.success_file , shell=True )

    #################################################################


