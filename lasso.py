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


import os
import subprocess
from functools import wraps
from time import gmtime, strftime , localtime
import datetime

from bal.lasso.lasso import *
from bal.core.prepare import get_executables

arguments = get_commandline_arguments()

if arguments.U != None:
    input_files = [arguments.U]
elif arguments.mate_1 and arguments.mate_2:
    input_files = [arguments.mate_1, arguments.mate_2]
else:
    print("Error in input files. Either Provide fastq files in the -U option"
          ", or provide paired end files in the -1 -2 options.")
    exit()

for f in input_files:
    if not os.path.isfile(f):
        raise(FileNotFoundError(f))

os.makedirs(arguments.o, exist_ok = True)
reference_directory = os.path.join(arguments.o, 'reference')
os.makedirs(reference_directory, exist_ok = True)
sequences_file = os.path.join(reference_directory, 'sequences.fa')
make_reference_sequences(arguments, sequences_file)
bal_directory    = os.path.dirname(os.path.realpath(__file__))
executables      = get_executables(bal_directory)
ref_base = os.path.join(reference_directory, 'reference')
ref_result = make_bt2_reference(reference_base = ref_base,
                                fasta_file = sequences_file ,
                                executable = executables['bowtie2-build'])

alignment_folder = os.path.join(arguments.o, 'alignment')
os.makedirs(alignment_folder, exist_ok = True)
output_sam = os.path.join(alignment_folder, "alignment.sam")

if not os.path.isfile(output_sam):
    print("Aligning Reads...")
    alignment_result = align_reads(input_files = input_files,
                                   reference_base = ref_base,
                                   output_sam = output_sam,
                                   executable = executables['bowtie2'],
                                   parameters = " -D 20 -R 3 -N 1 -L 10 -i S,1,0.50 ")
else:
    print(output_sam +  " exists, skipping alignment step.")

print("Picking bp alignments.")

filtered_sam = os.path.join(arguments.o, "alignments_through_bp.sam")

find_bp_alignments(input_sam = output_sam, output_sam = filtered_sam,
                   radius = arguments.radius, coverage_threshold = arguments.min_bp_coverage)

