#
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
import subprocess
from glob import glob
from collections import OrderedDict, defaultdict

from .step import Step
from .exceptions import *
from ..settings import *
from ..genomic_io.functions import bp_coverage_filter_sam, merge_sam_files
import pysam

#################################################################

class BpBedFromSam(Step):
    '''
        From the sam alignment file, extracts branchpoints and reports them in bed file.
        It also filters sam file so that only the alignments that cover the branchpoint
        are taken. The branchpoints are reported according to this filtered data.
    '''
    def __init__(self, name, input_files, output_directory,
                 executable='', executable_arguments = ''):
        super().__init__(name, input_files, output_directory, executable, executable_arguments)
        self.bed_file          = os.path.join(self.output_directory, settings['bp_candidate_file_name'])
        self.bed_entries       = OrderedDict()
        self.filtered_sam_file = os.path.join(self.output_directory, "bp_coverage_filtered_alignment.sam")
        self.filtered_files    = list()

    ###################################################################
    def prepare(self):
        self.command = ""

    ##############################################################################
    def _module_run(self):

        for file in self.input_files:
            file_directory = os.path.dirname(file)
            file_name = os.path.basename(file)
            filtered_file = os.path.join(file_directory, "bp_coverage_filtered_" + file_name )
            bp_coverage_filter_sam(file, filtered_file)

            if not os.path.isfile(filtered_file):
                raise FileNotFoundError("Could not create the filtered sam file" + filtered_file )

            self.process_sam_file(filtered_file)
            self.filtered_files.append(filtered_file)
        self.write_bed_file()
        merge_sam_files(self.filtered_files, self.filtered_sam_file)

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
            if read.is_unmapped:
                continue
            aligned_ref_label    = samfile.getrname(read.reference_id)
            aligned_ref_contents = aligned_ref_label.split(settings['field_separator'])

            if self.bed_entries.get(aligned_ref_label) == None:
                self.bed_entries[aligned_ref_label] = [aligned_ref_contents[0],
                                                       aligned_ref_contents[4],
                                                       str(int(aligned_ref_contents[4]) + 1),
                                                       aligned_ref_label,
                                                       1,
                                                       aligned_ref_contents[1]]
            else:
                self.bed_entries[aligned_ref_label][4] += 1

        samfile.close()

    ###############################################################################

    def post_run(self):
        if not os.path.isfile(self.bed_file):
            self.error_message = ("Couldn't find the bed file", self.bed_file)
        super().post_run()
