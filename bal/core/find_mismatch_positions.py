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
from collections import OrderedDict, defaultdict
import re

from .step import Step
from .exceptions import *
from ..settings import *
import pysam

#################################################################

class FindMismatchPositions(Step):
    '''
    To be completed
    '''
    def __init__(self, name, fasta_file, sam_file, output_directory,
                 executable='', executable_arguments = '' ):
        super().__init__(name, [], output_directory, executable, executable_arguments)

        self.fasta_file = os.path.abspath(fasta_file)
        self.sam_file   = os.path.abspath(sam_file)

        self.sam_file_with_annotated_mismatches = os.path.join(self.output_directory, "annotated_mismatches.sam")

        self.mismatch_counts_per_branchpoint       = defaultdict(int)
        self.mismatch_counts_at_bp_per_branchpoint = defaultdict(int)
        self.read_counts_per_branchpoint           = defaultdict(int)
        self.mismatch_records_nuc_based            = list()

        self.mismatch_by_position_file = os.path.join( self.output_directory , 'mismatch_by_position.txt')
        self.bed_output_file           = os.path.join( self.output_directory , 'bp_candidates_with_mismatch_ratio.bed')

        self.bed_entries = dict()
        
    ###################################################################
    def prepare(self):
        self.command = ""
        os.makedirs(self.output_directory, exist_ok=True)

    ##############################################################################

    def _get_reference_read_tag(self,read,tag):
        for elt in read.tags:
            if elt[0] == tag:
                return elt[1]
        raise(KeyError('The key' + tag + ' not found in '+ str(read.tags)))

    ##############################################################################

    def _parse_MD_tag(self, MD_string):

        result = list()
        my_re_main = re.compile('(?P<first>[0-9]+)(?P<rest>.*)')
        my_re_internal = re.compile('(?P<first>([A-Z]|\^[A-Z]+))(?P<rest>.*)')
        search_result = my_re_main.search(MD_string)
        while search_result:
            result.append(search_result.group('first'))
            internal_result = my_re_internal.search( search_result.group('rest')  )
            if internal_result:
               result.append( internal_result.group('first'))
               search_result = my_re_main.search( internal_result.group('rest')  )
            else:
               return result

        return result

    ############################################################################

    def _find_index_at_read(self, ref_index, read):
        result = -1
        for pair in read.aligned_pairs:
            if ref_index == pair[1]:
                if pair[0] == None:
                    print('Warning: Unexpected read ref found: None for' + str(read))
                    print('<' , ref_index , '>')
                    print('------\n', read.aligned_pairs  , '---------\n')
                    continue
                return(pair[0])

        if result == -1:
            print('Warning: Could not find the relative index in the read' + str(read))
        return(0)

    ############################################################################

    def _extract_mismatch_locations(self, read, MD_array, this_bp):
        # this is the active index of the sequenced read
        ref_offset = 0
        for elt in MD_array:
            if len(elt) < 1:
                continue
            if elt[0] in list( map(str, range(0,10)) ):
                ref_offset += int(elt)
            elif elt[0] == '^':
                ref_offset += (len(elt) - 1)

            elif elt[0] in ('A', 'C', 'G', 'T'):
                self.mismatch_counts_per_branchpoint[ this_bp ] += 1
                index_at_read = self._find_index_at_read(read.pos + ref_offset, read)
                misincorporated_nucleotide = read.query[index_at_read].upper()
                if not misincorporated_nucleotide in ('A', 'C', 'G', 'T'):
                    continue
                self.mismatch_records_nuc_based[read.pos + ref_offset][misincorporated_nucleotide] += 1
                self.mismatch_records_nuc_based[read.pos + ref_offset]['total'] += 1
                if (read.pos + ref_offset) == (( int(self.ref_length) - 1) / 2):
                    self.mismatch_counts_at_bp_per_branchpoint[this_bp] += 1
                ref_offset += 1
            else:
                print('Warning: Uexpected MD field:', elt)


    #############################################################################

    def _process_MD_tag(self,read, this_bp):
        this_MD_tag = self._get_reference_read_tag(read, 'MD')
        MD_contents = self._parse_MD_tag(this_MD_tag)
        self._extract_mismatch_locations(read, MD_contents, this_bp)
        return( MD_contents )


    ##############################################################################

    def process_annotated_sam_file(self, sam_file):
        sam_file = pysam.AlignmentFile(self.sam_file, "r")

        ref_length = -1
        for header_elt in sam_file.header.items():
            if header_elt[0] == "SQ":
                for seq_info in  header_elt[1]:
                    # this is the length
                    if seq_info.get('LN', -2)  > ref_length:
                        ref_length = seq_info.get('LN', -2)

        if ref_length < 0 :
            raise(Exception('There is a problem with getting the ref sequence error from the sam file header. '
                            'Check the SQ and LN fields in the header.'))

        self.ref_length = ref_length

        for i in range(0, self.ref_length):
            self.mismatch_records_nuc_based.append({'A' : 0, 'C' : 0, 'G' : 0, 'T' : 0, 'total' : 0})

        for read in sam_file.fetch():
            #print('----------------------------')
            #print(read)
            #print(read.aligned_pairs)
            #print(dir(read))
            #print('MD tag:' + self._get_reference_read_tag(read, 'MD'))
            this_bp = sam_file.getrname(read.reference_id)
            self.read_counts_per_branchpoint[this_bp] += 1
            #print( self._process_MD_tag( read, this_bp ) )
            self._process_MD_tag( read, this_bp )
            #print('-----------------------------')


            #print(dir(read))

    #############################################################################

    def _write_output(self):

        header = 'nucleotide\t' + '\t'.join(map(str, range(1,len(self.mismatch_records_nuc_based) + 1)) )

        with open(self.mismatch_by_position_file, 'w') as output_stream:
            print(header, file = output_stream)

            for nuc in ['A', 'C', 'G', 'T', 'total']:
                this_line = nuc + '\t'
                counts = list()
                for i in range(0, len( self.mismatch_records_nuc_based ) ):
                    counts.append( self.mismatch_records_nuc_based[i][nuc] )
                this_line += '\t'.join( map(str, counts) )
                print(this_line, file = output_stream)

        samfile = pysam.AlignmentFile(self.sam_file, 'r')

        # the fifth column is for the number of read counts
        # the seventh column is the the number of mismatches at the branchpoint
        # the eigth column is the percentage of mismatches at the branchpoint

        with open(self.bed_output_file , 'w') as output_stream:

            for read in samfile.fetch():
                if read.is_unmapped:
                    continue
                aligned_ref_label    = samfile.getrname(read.reference_id)
                aligned_ref_contents = aligned_ref_label.split(settings['field_separator'])

                if self.bed_entries.get(aligned_ref_label) == None:
                    percentage = "%.2f"%( 100 * (self.mismatch_counts_at_bp_per_branchpoint[aligned_ref_label] /\
                                                 self.read_counts_per_branchpoint[aligned_ref_label]) )
                    self.bed_entries[aligned_ref_label] = [aligned_ref_contents[0],
                                                           aligned_ref_contents[4],
                                                           str(int(aligned_ref_contents[4]) + 1),
                                                           aligned_ref_label,
                                                           str(self.read_counts_per_branchpoint[aligned_ref_label]),
                                                           aligned_ref_contents[1],
                                                           str( self.mismatch_counts_at_bp_per_branchpoint[aligned_ref_label], ),
                                                           percentage
                                                           ]
                    print("\t".join(self.bed_entries[aligned_ref_label]), file = output_stream)

        samfile.close()





    ##############################################################################

    def _module_run(self):

        #pysam.calmd("-e", "-S", self.sam_file, self.fasta_file )
        annotated_file_contents = pysam.calmd("-e", "-S", self.sam_file, self.fasta_file)

        with open(self.sam_file_with_annotated_mismatches, "w") as output_sam_stream:
            for elt in annotated_file_contents:
                print(elt.strip(), file=output_sam_stream)

        self.process_annotated_sam_file(self.sam_file_with_annotated_mismatches)
        self._write_output()

    ###############################################################################


    ###############################################################################

    def post_run(self):
        # TO BE COMPLETED!!!
        error_messages = list()

        if len(error_messages) > 0:
           subprocess.call('touch ' + self.failure_file , shell=True )
        else:
           subprocess.call('touch ' + self.success_file , shell=True )

        self.error_messages = error_messages

