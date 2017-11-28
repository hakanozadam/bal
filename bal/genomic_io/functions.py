
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
from .fasta import FastaEntry, FastaFile
from .fastq import FastqFile, FastqEntry
import pysam
from ..settings import *

def make_fasta_from_fastq(input_file, output_file):
    '''
    Given a fastq file (input_file), converts it to a fasta file (output_file)
    '''
    with FastqFile(input_file) as input_stream,\
        open(output_file, 'w') as fasta_stream:
        for fastq_entry in input_stream:
            fasta_entry = FastaEntry(header = fastq_entry.header, 
                                    sequence = fastq_entry.sequence)
            print(fasta_entry, file = fasta_stream)

###################################################################

def bp_coverage_filter_sam(input_sam_file, output_sam_file,
                           reference_half_length = settings['number_of_nucleotides'],
                           min_required_nucs_on_coverage = settings['bp_side_coverage_threshold'] ):
    '''
    Given a sam file, keeps only those enties that have coverage of at least min_required_nucs_on_coverage nucleotides
    aroiund the center of the reference.
    Note that the center of the reference is determined by reference_half_length

    :param input_sam_file:
    :param output_sam_file:
    :param reference_half_length:
    :param min_required_nucs_on_coverage:
    '''

    input_samfile_stream   = pysam.AlignmentFile(input_sam_file,  "r")
    output_samfile_stream  = pysam.AlignmentFile(output_sam_file, "wh", header = input_samfile_stream.header)

    for read in input_samfile_stream.fetch():
        if read.is_unmapped:
            continue
        # make sure the aligned read covers at least min_required_nucs_on_coverage
        # nucleotides on either side of the coverage!
        if (reference_half_length - read.reference_start) < min_required_nucs_on_coverage or\
                (read.reference_end - reference_half_length) < min_required_nucs_on_coverage:
            print( read.reference_start ,"   " ,  read.reference_end , "   ", str(read))
            continue

        output_samfile_stream.write(read)

    input_samfile_stream.close()
    output_samfile_stream.close()

#########################################################################

def merge_sam_files(input_sam_list, merged_sam_file):
    '''
    Merges the given input sam files into the target merged_sam_file
    :param input_sam_list:
    :param merged_sam_file:
    '''
    if len(input_sam_list) < 1:
        return

    input_samfile_stream = pysam.AlignmentFile(input_sam_list[0],  "r")
    output_samfile_stream  = pysam.AlignmentFile(merged_sam_file, "wh", header = input_samfile_stream.header)
    input_samfile_stream.close()

    for file in input_sam_list:
        input_samfile_stream = pysam.AlignmentFile(file,  "r")
        for read in input_samfile_stream.fetch():
            if read.is_unmapped:
                continue
            output_samfile_stream.write(read)
        input_samfile_stream.close()
    output_samfile_stream.close()

#########################################################################

def determine_5p_sequence_length(five_p_sequence_fasta_file):
    result = 0

    if not os.path.isfile(five_p_sequence_fasta_file):
        raise (FileNotFoundError(five_p_sequence_fasta_file + " not found!"))

    with FastaFile(five_p_sequence_fasta_file) as fasta_input:
        for fasta_entry in fasta_input:
            this_length = len(fasta_entry.sequence)
            if this_length > result:
                result = this_length

    if result == 0:
        error_message = "There is a problem with the five prime intron sequences." \
                        "The maximum length was determined as 0 in the file" +\
                        five_p_sequence_fasta_file
        raise(IOError(error_message))

    return result


