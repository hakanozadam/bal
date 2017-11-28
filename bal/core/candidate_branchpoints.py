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
from collections import OrderedDict
from glob import glob

from .step import Step
from .exceptions import *
from ..genomic_io.fastq import FastqFile
from ..genomic_io.fasta import FastaFile, FastaEntry
from ..genomic_io.functions import make_fasta_from_fastq
from ..annotation.intron import get_intron_sequences
from ..settings import *
import pysam

#################################################################

class CandidateBranchpoints(Step):
    '''
        The input_files variable contains the directories containing alignment sam files
    '''
    def __init__(self, name, input_files, output_directory,
                 executable='', executable_arguments = ''):
        self.alignment_directories = map(os.path.abspath, input_files)
        super().__init__(name, [], output_directory, executable, executable_arguments)
        self.bed_file = os.path.join(self.output_directory, settings['bp_candidate_file_name'])
        self.bed_entries = OrderedDict()

    ###################################################################
    def prepare(self):
        self.command = ""

    ##############################################################################
    def _module_run(self):

        for directory in self.alignment_directories:
            if not os.path.isdir(directory):
                raise(FileNotFoundError("Could not find the directory ", directory))
            gene_directories = glob(directory + "/*")
            for gene_directory in gene_directories:
                if not os.path.isdir(gene_directory):
                    continue
                sam_file = os.path.join(gene_directory, "alignments.sam")
                if not os.path.isfile(sam_file):
                    raise(FileNotFoundError("Could not find the sam file" + sam_file))
                self.process_sam_file(sam_file)
        self.write_bed_file()

    ###############################################################################

    def write_bed_file(self):
        with open(self.bed_file, "w") as out_bed_file:
            for this_entry in self.bed_entries.values():
                print( "\t".join(map(str,this_entry))  , file = out_bed_file)

    ###############################################################################

    def process_sam_file(self, input_samfile):
        # Branchpoint candidates are reported in bed file
        # The entry name is in the following format
        # chromosome__strand__gene__bpLocation__3'distance___5'end__3'end1__3'end2__...3'endN

        samfile = pysam.AlignmentFile(input_samfile, "r")

        for read in samfile.fetch():
            read_id_and_intron = read.qname.split(settings['fastq_field_separator'])
            if len(read_id_and_intron) < 2 :
                raise(StepError("There is a problem with the sam header " + read.qname +\
                      " in the file " + self.sam_file ))
            intron_contents = read_id_and_intron[1].split(settings['field_separator'])

            if len(intron_contents) < 5:
                print("Warning: There is a problem with the entry " , read.qname, "\nIt is missing fields." )
                continue

            #ref_chr = samfile.getrname(read.reference_id)
            aligned_ref_label = samfile.getrname(read.reference_id)
            aligned_ref_contents = aligned_ref_label.split(settings['field_separator'])
            if len(aligned_ref_contents) < 6:
                print("Warning: There is a problem with the entry", aligned_ref_label,
                      "in the sam file", input_samfile)
                continue
            ref_chr    = aligned_ref_contents[0]
            gene_start = int(aligned_ref_contents[1])
            gene_end   = int(aligned_ref_contents[2])

            five_prime_end     = int(intron_contents[3])
            three_prime_ends   = tuple(map(int,intron_contents[4:]))

            distance_to_three_prime = -1
            branchpoint_location =    -1

            # make sure trimmed read is mapped to the right chr
            if intron_contents[0] != ref_chr:
                continue

            # make sure the mapped strand and the intron strand are the same
            if intron_contents[1] == '+' :
                branchpoint_location = read.reference_end + gene_start - 1
                if branchpoint_location < five_prime_end:
                    continue
                # The read shouldn't be too close to the five prime end.
                if read.pos + gene_start - five_prime_end < 30:
                    continue

                distance_to_three_prime = three_prime_ends[0] - branchpoint_location
                for tp in three_prime_ends[1:]:
                    this_distance = tp - branchpoint_location
                    if this_distance > 0 and this_distance < distance_to_three_prime:
                        distance_to_three_prime = this_distance
            elif intron_contents[1] == '-':
                branchpoint_location = gene_end - read.reference_end
                if branchpoint_location >= five_prime_end:
                    continue
                distance_to_three_prime = branchpoint_location - three_prime_ends[0]
                for tp in three_prime_ends[1:]:
                    this_distance = branchpoint_location - tp
                    if this_distance > 0 and this_distance < distance_to_three_prime:
                        distance_to_three_prime = this_distance
            else:
                print("Warning: Invalid strand " + intron_contents[1] + " in entry " + read.qname +\
                      " in file " + self.sam_file)

            if distance_to_three_prime < 0 :
                continue

            bed_array = [ref_chr, intron_contents[1], intron_contents[2],
                         str(five_prime_end - 1), str(branchpoint_location),
                         str(distance_to_three_prime) ] +\
                        list( map(str, list(map(lambda x:x-1, three_prime_ends)) ) )

            bed_entry_name = settings['field_separator'].join( bed_array )

            if self.bed_entries.get(bed_entry_name) == None:
                bed_entry = [ref_chr, branchpoint_location  ,
                                         branchpoint_location + 1, bed_entry_name,
                                         1, intron_contents[1]]
                self.bed_entries[bed_entry_name] = bed_entry
            else:
                self.bed_entries[bed_entry_name][4] += 1

        samfile.close()

    ###############################################################################

    def post_run(self):
        if not os.path.isfile(self.bed_file):
            self.error_messages.append("Couldn't find the bed file", self.bed_file)

        if len(self.error_messages) > 0:
           subprocess.call('touch ' + self.failure_file , shell=True )
        else:
           subprocess.call('touch ' + self.success_file , shell=True )
