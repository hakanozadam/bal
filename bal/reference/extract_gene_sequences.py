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

from shutil import which
import os
import subprocess
from ..core.step import Step
from ..core.exceptions import *
from ..genomic_io.fasta import FastaFile, FastaEntry
from ..genomic_io.gtf import GtfFile, GtfEntry
from ..settings import *

from collections import OrderedDict, defaultdict

#################################################################

class ExtractGeneSequences(Step):
    '''Extract gene Codrdinates and produce fasta files'''
    def __init__(self, name, input_files, output_directory):
        '''
        N is the number of nucleotides to be taken from the five prime end
        '''

        super().__init__(name, list(), output_directory)

        self.gtf_file           = os.path.abspath(input_files[0])
        self.fasta_file         = os.path.abspath(input_files[1])
        self.output_directory   = os.path.abspath(self.output_directory)
        self.output_file        = os.path.join(self.output_directory, "gene_sequences")
        self.gene_bed_file      = os.path.join(self.output_directory, "genes.bed")

        self.gene_directory = os.path.join(self.output_directory, 'genes')
        os.makedirs(self.gene_directory, exist_ok=True)


        os.makedirs(self.output_directory, exist_ok=True)

        missing_files = list()
        for f in (self.gtf_file, self.fasta_file):
           if not os.path.isfile(f):
              missing_files.append(f)

        if len(missing_files) > 0:
           raise FileNotFoundError("The following files are missing:\n" +
                                   "\n".join(missing_files))

        self.genes              = OrderedDict()

        # keys are chr and values are genes of that chr in a list form
        self.chr_genes          = defaultdict(list) # keys are chr and values are genes of that chr in a list form

    ##############################################################

    def prepare(self):
        pass
            
    ###############################################################

    def make_gene_bed_file(self):
        print("Making gene bed file and getting coordinates...")
        with GtfFile(self.gtf_file) as gtf_input:
           for entry in gtf_input:
              this_gene_id = entry.attribute_contents["gene_id"]
              this_gene    = self.genes.get(this_gene_id, None)
              # print(this_gene , "---", entry)
              if not this_gene:
                 self.genes[this_gene_id] = [entry.seqname, entry.start, entry.end, this_gene_id, "0", entry.strand]
                 self.chr_genes[entry.seqname].append(this_gene_id)
              else:
                 if int(this_gene[1]) > int(entry.start):
                    this_gene[1] = entry.start
                 if int(this_gene[2]) < int(entry.end):
                    this_gene[2] = entry.end

        fasta_genes = list()
        with FastaFile(self.fasta_file) as fasta_input:
           for fasta_entry in fasta_input:
              fasta_genes.append(fasta_entry.header)


        with open(self.gene_bed_file, 'w') as bed_output:
           for gene_id, gene_entry in self.genes.items():
              if not gene_entry[0] in fasta_genes:
                 continue
              gene_entry[1] = str(gene_entry[1] - 1)
              gene_entry[2] = str(gene_entry[2])
              print("\t".join(gene_entry), file=bed_output)


    ###############################################################

    def make_fasta_files(self):
        print("Getting gene sequences...")
        i=0
        with FastaFile(self.fasta_file) as fasta_input:
           for fasta_entry in fasta_input:
              if not self.chr_genes.get(fasta_entry.header):
                 continue
              for gene_id in self.chr_genes[fasta_entry.header]:
                 gene_bed_entry = self.genes[gene_id]
                 gene_sequence = fasta_entry.sequence[int(gene_bed_entry[1]) : int(gene_bed_entry[2]) ]

                 N_counter = 0

                 for nucleotide in gene_sequence:
                    if nucleotide not in ['A', 'C', 'G', 'T', 'a', 'c', 'g', 't' ]:
                       N_counter += 1

                 if N_counter >= len(gene_sequence) -1 :
                    continue

                 gene_header   = settings['field_separator'].join(gene_bed_entry)
                 gene_fasta_entry = FastaEntry(header=gene_header, sequence=gene_sequence)
                 if gene_bed_entry[5] == "-":
                    gene_fasta_entry.reverse_complement()
                 gene_file = os.path.join(self.gene_directory, gene_id + ".fa")
                 with open(gene_file, 'w') as gene_output:
                    print(gene_fasta_entry, file = gene_output)

                 i += 1
                 if i % 1000 == 0:
                    print('\r' + str(i) + " genes have been extracted....", end= '\r')

    ###############################################################

    def _module_run(self):
        self.make_gene_bed_file()
        self.make_fasta_files()

    ################################################################

    def post_run(self):
        missing_files = list()
        required_files = [self.gene_bed_file]



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


