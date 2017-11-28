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
####################################################################

import argparse
import os
import sys
from collections import defaultdict

#####################################################################

def handle_arguments():
   parser = argparse.ArgumentParser(description=
   '''
   ??????
   ''')
   parser.add_argument("-f" ,
                       help = "Input fasta file" ,
                       required = True ,
                       metavar = "input_fasta_file" ,
                       type = str)
   parser.add_argument("-b" ,
                       help = "Input bed file" ,
                       required = True ,
                       metavar = "input_bed_file" ,
                       type = str)
   parser.add_argument("-n" ,
                       help = "Radius: Number of nucleotides to the left and to the right of the bp" ,
                       required = True ,
                       metavar = "radius" ,
                       type = int)
   parser.add_argument("-o" ,
                       help = "Output directory" ,
                       required = True ,
                       metavar = "output_fastq_file" ,
                       type = str)

   arguments = parser.parse_args()
   return arguments

#####################################################################

def main():
    arguments = handle_arguments()
    os.makedirs(arguments.o, exist_ok=True)
    position_index = 1
    strand_index = 5
    label_index  = 3
    chr_index = 0

    #contains parsed bed entries
    # they are grouped by chromosomes
    branchpoints_by_chr = defaultdict(dict)

    #keys are positions relative to the branchpoint
    # values are total number of nucleotides in the given position
    nucleotide_counts = dict()

    nucleotides = ("A", "C", "G", "T")
    radius = arguments.n

    for nuc in nucleotides:
        nucleotide_counts[nuc] = dict()
        for i in range(0,radius*2 + 1):
            nucleotide_counts[nuc][i] = 0

    with open(arguments.b, 'r') as bed_input_stream:
        for line in bed_input_stream:
            contents = line.split()
            if len(contents) < 6:
                continue
            branchpoints_by_chr[contents[chr_index]][contents[label_index]] = contents

    motif_bed_file = os.path.join(arguments.o, "sequences_around_branchpoints.bed")

    with FastaFile(arguments.f) as fasta_input_stream,\
        open(motif_bed_file, "w") as bed_output_stream:
        for fasta_entry in fasta_input_stream:
            current_chr = fasta_entry.header
            if branchpoints_by_chr.get(current_chr, None) == None:
                continue
            for bp in branchpoints_by_chr[current_chr].values():
                bp_position = int( bp[position_index] )
                sequence_around_bp = (fasta_entry.sequence[bp_position - radius : bp_position + radius + 1]).upper()

                if bp[strand_index] == "-":
                    dummy_entry = FastaEntry(header = "dummy", sequence = sequence_around_bp)
                    dummy_entry.reverse_complement()
                    sequence_around_bp = dummy_entry.sequence

                for i in range(0,radius*2 + 1):
                    current_nucleotide = sequence_around_bp[i]
                    if current_nucleotide in nucleotides:
                        nucleotide_counts[current_nucleotide][i] += 1

                print( "\t".join( bp )  + "\t" +  sequence_around_bp, file = bed_output_stream)

    nucleotide_distribution_file = os.path.join(arguments.o, "nucleotide_distribution.txt")

    with open(nucleotide_distribution_file, "w") as nucleotide_distribution_output_stream:
        header = "Poisition\t" + "\t".join( list( map( str, map(lambda x: x+1, range(0, radius*2 + 1) ) ) ) )
        print(header, file = nucleotide_distribution_output_stream)

        for nuc in nucleotides:
            this_line =  nuc +  "\t" + "\t".join( map( str, map(lambda i: nucleotide_counts[nuc][i], range(0,radius*2 + 1) ) ) )
            print(this_line, file = nucleotide_distribution_output_stream)




#####################################################################ß

script_directory = os.path.abspath(os.path.dirname(os.path.realpath(__file__)) )
bal_dir = os.path.split( os.path.split( script_directory)[0] )[0]

if __name__ == '__main__':

    sys.path.append( bal_dir  )
    from genomic_io.fasta import FastaFile, FastaEntry
    from genomic_io.fastq import FastqFile, FastqEntry
    
    main()
else:
    exit(1)

####################################################################
