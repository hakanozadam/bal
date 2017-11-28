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
import argparse



from ..genomic_io.fastq import FastqFile
from ..genomic_io.fasta import FastaFile, FastaEntry

import pysam

#################################################################

def get_commandline_arguments():

   parser = argparse.ArgumentParser(description=
   '''
   A Modified Version of Lasso Algorithm

   The user provides the coordinates, gene name, strand and chrsomosome of the intron.
   The provided fastq file is mapped to all possible branchpoints from that intron.
   ''')
   parser.add_argument("-f" ,
                       help = "Input fasta file" ,
                       required = True ,
                       metavar = "input_fasta_file" ,
                       type = str)
   parser.add_argument("-U" ,
                       help = "Input fastq file single end" ,
                       required = False ,
                       metavar = "input_fastq_file" ,
                       type = str)
   parser.add_argument("-1" ,
                       help = "Input fastq file, mate 1, paired end" ,
                       required = False ,
                       metavar = "input_fastq_file_mate_1" ,
                       dest = 'mate_1',
                       type = str)
   parser.add_argument("-2" ,
                       help = "Input fastq file, mate 2, paired end" ,
                       required = False ,
                       metavar = "input_fastq_file_mate_2" ,
                       dest = 'mate_2',
                       type = str)
   parser.add_argument("-o" ,
                       help = "Output directory" ,
                       required = True ,
                       metavar = "output_fastq_file" ,
                       type = str)
   parser.add_argument("--start" ,
                       help = "Intron start coordinate. 0-based inclusive." ,
                       required = True ,
                       metavar = "intron_start" ,
                       type = int)
   parser.add_argument("--end" ,
                       help = "Intron end coordinate. 0-based inclusive." ,
                       required = True ,
                       metavar = "intron_end" ,
                       type = int)
   parser.add_argument("--chr" ,
                       help = "Chromosome" ,
                       required = True ,
                       metavar = "chromosome of the intron" ,
                       type = str)
   parser.add_argument("--gene" ,
                       help = "Gene Name" ,
                       required = True ,
                       metavar = "gene" ,
                       type = str)
   parser.add_argument("--strand" ,
                       help = "Strand" ,
                       required = True ,
                       choices = ['+', '-'],
                       metavar = "Strand" ,
                       type = str)
   parser.add_argument("--radius" ,
                       help = "Radius: Number of nucleotides to the left and to the right of the bp" ,
                       required = True ,
                       metavar = "radius" ,
                       type = int)
   parser.add_argument("--min-bp-coverage" ,
                       help = "Minimum number of nucleotids needed to cover either side of the candidate branchpoint" ,
                       required = False ,
                       default = 5,
                       metavar = "min_bp_coverage" ,
                       dest = 'min_bp_coverage',
                       type = int)

   return parser.parse_args()

###########################################################################

def get_intron_sequence(fasta_file, chromosome, intron_start, intron_end, strand):
    found_chr = False

    with FastaFile(fasta_file) as opened_fasta_file:
        for entry in opened_fasta_file:
            if entry.header == chromosome:
                found_chr = True
                if intron_start >= len(entry.sequence) or\
                    intron_end >= len(entry.sequence):
                    print('The intron boundaries are out of range.'
                          'Intron start = ', intron_start, 'intron_end = ', intron_end,
                          'chromosome has', len(entery.sequence), 'nucleotides.')
                    raise(IndexError('Out of chr range'))

                intron_fasta_entry = FastaEntry(header = 'intron_identifier',
                                                sequence = entry.sequence[intron_start:(intron_end + 1)])
                if strand == '-':
                    intron_fasta_entry.reverse_complement()

                print(intron_fasta_entry)

        if found_chr == False:
            raise(Exception('Could not find the chromosome ', chromosome))

    return intron_fasta_entry.sequence

##############################################################################

def make_reference_sequences(arguments, output_fasta_file):

    intron_sequence = get_intron_sequence(arguments.f, arguments.chr, arguments.start, arguments.end, arguments.strand)

    this_five_prime_end = arguments.start
    this_three_prime_end = arguments.end

    radius = arguments.radius

    if arguments.strand == '-':
        this_five_prime_end = arguments.end
        this_three_prime_end = arguments.start

    with open(output_fasta_file, 'w') as output_stream:

        for i in range(arguments.min_bp_coverage, len(intron_sequence)):

            this_distance_to_three_prime = len(intron_sequence) - i - 1
            this_bp_location = arguments.start + i

            piece_A_start = max( i - radius, 0 )
            piece_B_start = 0
            piece_A_end = i + 1
            piece_B_end = min( radius, len(intron_sequence) )

            if arguments.strand == "-":
                this_bp_location = len(intron_sequence) - i - 1 + arguments.start

            bp_label = "__".join( (
                str(i), arguments.chr, arguments.strand, arguments.gene, str(this_five_prime_end),
                str(this_bp_location), str(this_distance_to_three_prime), str(this_three_prime_end) )
            )

            this_sequence = intron_sequence[  piece_A_start : piece_A_end  ] +\
                            intron_sequence[piece_B_start : piece_B_end]

            this_fasta_entry = FastaEntry( header = bp_label, sequence =  this_sequence)
            print(this_fasta_entry, file = output_stream)

#####################################################################################

def make_bt2_reference(reference_base, fasta_file , executable):
    command = " ".join((executable , fasta_file, reference_base))
    p = subprocess.Popen([command], stdout=subprocess.PIPE,
                                                  stderr=subprocess.PIPE,
                                                  shell = True)
    std_out , std_err = p.communicate()
    std_out = std_out.decode("utf-8").rstrip()
    std_err = std_err.decode("utf-8").rstrip()

    returncode = p.returncode
    if returncode:
        print(std_out, std_err)
    return returncode

########################################################################################

def align_reads(input_files, reference_base, output_sam, executable, parameters):
    if len(input_files) == 1:
        command_arg = " -U " + input_files[0]
    else:
        command_arg = " -1 % -2 % "%input_files
    command_arg +=  " -x " + reference_base +  " -S " + output_sam + " --no-unal " + parameters

    command = " ".join( (executable, command_arg) )
    print(command)
    p = subprocess.Popen([command], stdout=subprocess.PIPE,
                                                  stderr=subprocess.PIPE,
                                                  shell = True)
    std_out , std_err = p.communicate()
    std_out = std_out.decode("utf-8").rstrip()
    std_err = std_err.decode("utf-8").rstrip()

    returncode = p.returncode
    print(std_out, std_err)
    return returncode


##########################################################################################

def find_bp_alignments(input_sam, output_sam, radius, coverage_threshold):

    input_file        = pysam.AlignmentFile(input_sam, 'r')
    output_sam_stream = pysam.AlignmentFile(output_sam, "wh", header = input_file.header)
    lengths = dict()
    for ref in input_file.header['SQ']:
        print(ref['SN'], ' --- ', ref['LN'])
        lengths[ref['SN']] = int(ref['LN'])


    for read in input_file:
        if read.is_unmapped:
            continue
        output_sam_stream.write(read)
        aligned_ref_label = input_file.getrname(read.reference_id)
        ref_contents = aligned_ref_label.split("__")
        bp_position = int(ref_contents[0])
        if bp_position - read.pos >= coverage_threshold and\
            read.aend - bp_position >= coverage_threshold:
            output_sam_stream.write(read)
