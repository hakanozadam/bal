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

import argparse
import os
from shutil import which
import sys
from collections import namedtuple
import re
import glob

#################################################################

def correct_command_line_arguments(parameter):
    corrected_arguments = list()
    correct_next = False
    for elt in sys.argv[1:]:
        if correct_next:
            correct_next = False
            corrected_arguments.append(" " + elt)
        else:
            corrected_arguments.append(elt)

        if elt == parameter:
            correct_next = True

    return corrected_arguments

##################################################################

def get_commandline_arguments():
    ''' Parse and return the command line arguments'''

    parser = argparse.ArgumentParser(description=
    '''
    BAL Pipeline For Umass Cluster

    This pipeline partitions the given data into smaller files
    and processes them using the main bal pipeline.
    It has been particularly implemented for umass hpc cluster.
    ''')

    parser.add_argument("-i" ,
                    metavar  = 'input_directory' ,
                    help     = "Input directory containing the fastq files." ,
                    required = True ,
                    type     = str)
    parser.add_argument("-o" ,
                    metavar  = 'output_directory' ,
                    help     = "Output directory" ,
                    required = True ,
                    type     = str)
    parser.add_argument("-x" ,
                    metavar  = 'reference_directory' ,
                    help     = "BAL reference directory." ,
                    required = True ,
                    type     = str)
    parser.add_argument("-n" ,
                    metavar  = 'reads_per_file' ,
                    help     = "Each big fasta file is partitioned into smaller files"
                               "each smaller file has n reads, possibly the last file "
                               "will have less than n reads" ,
                    default  = 100000,
                    required = False ,
                    type     = int)
    parser.add_argument("-a" ,
                    metavar  = 'additional_bal_arguments' ,
                    help     = "Additional commandline arguments to be passed to BAL." ,
                    required = False ,
                    default  = "",
                    type     = str)
    parser.add_argument("-m" ,
                    metavar  = 'alignment_mode' ,
                    help     = "Alignment mode: single, paired ."
                               "Are these Single-end libraries or paired end libraries?" ,
                    required = True ,
                    choices  = ['single' , 'paired'],
                    default  = "",
                    type     = str)
    parser.add_argument("-p" ,
                    metavar  = 'cores' ,
                    help     = "Number of cores for the aligners" ,
                    default  = 1,
                    required = False ,
                    type     = int)
    parser.add_argument("-t" ,
                    metavar  = 'time_limit' ,
                    help     = "time limit for BAL in hours" ,
                    default  = 3,
                    required = False ,
                    type     = int)
    parser.add_argument("-memory" ,
                    metavar  = 'memory_limit' ,
                    help     = "emory limit for BAL jobs" ,
                    default  = 8182 ,
                    required = False ,
                    type     = int)
    parser.add_argument("--rna-strandness" ,
                    metavar = 'strand of the library' ,
                    help = "This is the strand specificity of the RNA library. "
                           "F: forward, R: Reverse. The default is F\n"
                           "For single end reads:"
                           "F: A read corresponds to a transcript\n"
                           "R: A Read corresponds to the reverse complement of a transcript\n"
                           "For paired-end reads, the above holds for the read coming from mate-1\n\n" ,
                    choices = ["F", "R"],
                    default = "F",
                    required = False ,
                    dest = "rna_strandness",
                    type = str)
    return parser.parse_args( correct_command_line_arguments("-a") )

#################################################################################


##################################################################################

def get_arguments():
    command_line_arguments = get_commandline_arguments()

    error_messages = list()
    if not os.path.isdir(command_line_arguments.x):
        error_messages.append("Couldn't find the Bal reference directory %s"\
                              %command_line_arguments.x)
    if not os.path.isdir(command_line_arguments.i):
        error_messages.append("Couldn't find the input directory %s"\
                              %str(command_line_arguments.i) )

    if error_messages:
      print("Error!\nThere are error(s) in the command line arguments:")
      for error in enumerate(error_messages):
         print("{n}) {e}".format(n = error[0] +  1, e = error[1]))
      exit(1)

    return command_line_arguments

##################################################################################

def get_file_2(input_files, lib_name):
   ext2 = re.compile(r"(?P<ext2>.*)[._]2\.f(ast)?q$" , flags = re.IGNORECASE)
   for file in input_files:
       file_name = os.path.basename(file)
       ext2_search = ext2.search(file_name)
       if ext2_search and ext2_search.group("ext2") == lib_name:
           return file
   return ""

##################################################################################

def arrange_input_files(input_directory, alignment_mode):

       input_files = glob.glob(input_directory + "/*")
       Library = namedtuple("Library" , ["name" , "strand_1_file" , "strand_2_file"])
       arranged_input_files = list()
       ext  = re.compile(r"(?P<ext>.*)\.f(ast)?q$" , flags = re.IGNORECASE)
       ext1 = re.compile(r"(?P<ext1>.*)[._]1\.f(ast)?q$" , flags = re.IGNORECASE)
       ext2 = re.compile(r"(?P<ext2>.*)[._]2\.f(ast)?q$" , flags = re.IGNORECASE)

       for file in sorted(input_files):
           file_name = os.path.basename(file)

           if alignment_mode == "single":
               ext_search = ext.search(file_name)
               if not ext_search:
                   continue
               lib_name = ext_search.group("ext")
               arranged_input_files.append( Library(lib_name, file, "") )
           elif alignment_mode == "paired":
               ext1_search = ext1.search(file_name)
               if not ext1_search:
                   continue
               lib_name = ext1_search.group("ext1")
               file_2   = get_file_2(input_files, lib_name)
               if not file_2:
                   print("Warning, strand 2 file couldn't be found for the file %s"%file)
                   continue
               arranged_input_files.append(Library(lib_name, file, file_2))
           else:
               print("Unknown alignemnt mode %s"%alignment_mode)
               exit(1)
       return arranged_input_files

###################################################################################

def get_executables(main_directory):
   executables = dict()
   executables['partition_fastq'] = os.path.join(main_directory, 'bal', 'umass_cluster',
                                                 'partition_fastq_file.py' )
   executables['merge_bed_files'] = os.path.join(main_directory, 'bal', 'umass_cluster',
                                                 'merge_bed_files.py' )
   executables['make_bp_ref'] = os.path.join(main_directory, 'bal', 'umass_cluster',
                                                 'make_bp_ref.py' )
   executables['align_bp'] = os.path.join(main_directory, 'bal', 'umass_cluster',
                                                 'align_bp.py' )
   executables['bal']             = os.path.join(main_directory, 'align.py')
   executables['bowtie2-build']             = os.path.join(main_directory, 'bal', 'bin',
                                                           'bowtie2','linux_x86_64','bowtie2-build' )
   executables['bowtie2']             = os.path.join(main_directory, 'bal', 'bin',
                                                           'bowtie2','linux_x86_64','bowtie2' )
   executables['bp_bed_from_sam_ref'] = os.path.join(main_directory, 'bal', 'umass_cluster',
                                                 'extract_bp_from_sam_ref.py' )

   for key, value in executables.items():
       if not os.path.isfile(value):
           raise( FileNotFoundError("The file %s, which is the executable for %s, couldn't be found"\
               %(value, key)) )

   return executables

####################################################################################