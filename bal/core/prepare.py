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
from sys import platform as _os
from ..settings import *

#################################################################

def get_commandline_arguments():
    ''' Parse and return the command line arguments'''

    parser = argparse.ArgumentParser(description=
    '''
    TODO: Write a detailed explanation of the pipeline here!

    We need to revise this function. Make sure we need all the arguments and that
    these are exactly the arguments we need.

    Warning: In the parameters -1,-2,-U, we do not support comma separated files yet.
    Please give a single file for now. Multiple file support will be separated later.
    ''')

    parser.add_argument("-1" ,
                    metavar = 'mate_1_file' ,
                    help = "For paired-end data:The fastq file coming from strand (mate) 1."
                           " The reads in the matching files must match, i.e.,"
                           "they must come from the same RNA fragment (but from opposite ends and opposite strands)" ,
                    required = False ,
                    dest = 'mate_1',
                    type = str)
    parser.add_argument("-2" ,
                    metavar = 'mate_2_file' ,
                    help = "For paired-end data:The fastq file coming from strand (mate) 2."
                           " The reads in the matching files must match, i.e.,"
                           "they must come from the same RNA fragment (but from opposite ends and opposite strands)" ,
                    required = False ,
                    dest = 'mate_2' ,
                    type = str)
    parser.add_argument("-U" ,
                    metavar = 'Unpaired fastq file' ,
                    help = "For single-end data: Singe-end fastq file." ,
                    required = False ,
                    type = str)
    parser.add_argument("-x" ,
                    metavar = 'Reference Directory' ,
                    help = "This directory contains genomic and intronic references needed by HISAT and bowtie2." ,
                    required = True ,
                    type = str)
    parser.add_argument("-o" ,
                    metavar = 'Output Directory' ,
                    help = "Output directory" ,
                    required = True ,
                    type = str)
    parser.add_argument("-p" ,
                    metavar = 'number of threads' ,
                    help = "Output directory" ,
                    required = False ,
                    default = 1,
                    type = int)
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
    parser.add_argument('--hpc',
                        dest='hpc',
                        action='store_true',
                        help = 'Special mode in umass hpc cluster. '
                               'Causes Bal not to go beyond a certain step.')

    parser.add_argument('--5p-allowed-mismatch',
                        dest='five_p_allowed_mismatch',
                        default = 1,
                        required = False,
                        type = int,
                        help = 'The number of allowed nucleotide mismatches in 5p intron local alignment')

    parser.add_argument('--branchpoint-diameter',
                        dest='branchpoint_diameter',
                        default = 100,
                        required = False,
                        type = int,
                        help = 'The number of nucleotides on either side of the bp in bp reference sequence')

    parser.add_argument('--bp-side-min-coverage',
                        dest='bp_side_min_coverage',
                        default = 8,
                        required = False,
                        type = int,
                        help = 'The alignments agains the branchpoint are required to span at least this many'
                               'nucleotides on either side of the branchpoint.')

    parser.add_argument('--trim-length-threshold',
                        dest='trim_length_threshold',
                        default = 8,
                        required = False,
                        type = int,
                        help = 'The minimum length of the reads after the trimming step. '
                               'This piece is the other half of the candidate branchpoint read.')

    parser.set_defaults(hpc=False)
    return parser.parse_args()

#################################################################################

def check_reference_files(ref_dir):
   ''' TODO: Check for the existence of other files as well'''
   results = list()
   suffixes = ('.1.bt2', '.2.bt2', '.3.bt2', '.4.bt2', '.rev.1.bt2',  '.rev.2.bt2')
   hisat_ref_base = os.path.join(ref_dir , settings['hisat_genome_directory'] ,
                                 settings['hisat_genome_name'])

   intron_ref_base = os.path.join(ref_dir, settings['five_prime_intron_reference_directory'],
                                  settings['five_prime_intron_reference_name'])

   for suffix in suffixes:
      if (not os.path.isfile(hisat_ref_base + suffix) ) and\
         (not os.path.isfile(hisat_ref_base + suffix + "l")):
          results.append("Couldn't find the HISAT reference: " + hisat_ref_base + suffix + " or " +
          hisat_ref_base + suffix + "l")
      if not os.path.isfile(intron_ref_base + suffix) :
          results.append("Couldn't find the Bowtie2 Intron reference: " +\
              intron_ref_base + suffix)
   return results

#################################################################################

def process_commandline_arguments(cmd_args):
   ''' Check if the input files exist or not and do some consistency checks'''
   error_messages = list()

   if cmd_args.mate_1:
       setattr(cmd_args, 'paired' , True)
       setattr(cmd_args, 'U' , None)
       if(not os.path.isfile(cmd_args.mate_1)):
           error_messages.append("Couldn't find the input file " + cmd_args.mate_1)
       if(not os.path.isfile(cmd_args.mate_2)):
           error_messages.append("Couldn't find the input file " + cmd_args.mate_2)
   elif cmd_args.U:
       setattr(cmd_args, 'paired' , False)
       if(not os.path.isfile(cmd_args.U)):
           error_messages.append("Couldn't find the input file " + cmd_args.U)
   else:
       error_messages.append("No input deep-sequencing files were provided."
                             "Provide a single-end fastq file in the -U parameter or paired-end "
                             "fastq-files in the -1 and -2 parameters.")

   reference_result = check_reference_files(cmd_args.x)
   if reference_result:
       error_messages = error_messages + reference_result

   if error_messages:
      print("Error!\nThe following error(s) occurred:")
      for error in enumerate(error_messages):
         print("{n}) {e}".format(n = error[0] +  1, e = error[1]))
      exit(1)

   return cmd_args

##################################################################################
def get_arguments():
   return process_commandline_arguments(get_commandline_arguments())

###################################################################################



###################################################################################

def get_executables(bin_directory):
   ''' Check the existence of executables: hisat, bowtie2, samtools
   Put their paths in a dictionary and return it'''
   
   
   #check the os and define bin variables for executables accordingly
   
   if _os == "linux" or _os == "linux2":
       hisat_relative_path    = 'bal/bin/hisat/linux_x86_64'
       bowtie2_relative_path  = 'bal/bin/bowtie2/linux_x86_64'
       bowtie2_build_relative_path  = 'bal/bin/bowtie2/linux_x86_64/bowtie2-build'

       
   elif _os == "darwin":
       hisat_relative_path    = 'bal/bin/hisat/mac_os_x_x86_64'
       bowtie2_relative_path  = 'bal/bin/bowtie2/mac_os_x_x86_64'
       bowtie2_build_relative_path  = 'bal/bin/bowtie2/mac_os_x_x86_64/bowtie2-build'
       print(bowtie2_build_relative_path)

   executables = dict()
   error_messages = list()

   executables['hisat']   = os.path.join(bin_directory, hisat_relative_path, 'hisat')
   executables['hisat_extract_splice_sites'] = os.path.join(bin_directory, hisat_relative_path,\
                                                            'extract_splice_sites.py')

   executables['bowtie2'] = os.path.join(bin_directory, bowtie2_relative_path,'bowtie2')
   executables['bowtie2-build'] = os.path.join(bin_directory, bowtie2_build_relative_path)

   for executable, path in executables.items():
      if not which(path):
         error_messages.append("Couldn't find the {executable} executable at {path}"\
            .format(executable = executable, path = path))

   if(error_messages):
      print('The following executable(s) are missing. If you have the files in the indicated path,'
            'make sure that the files are executable.')
      print("\n".join(error_messages))
      exit(1)

   return executables

#######################################################################################

def get_bowtie2_5p_arguments(five_p_seq_length, allowed_mismatches):
    import math

    bowtie2_L  = "-L " + str( min(10, (math.ceil(int(five_p_seq_length) / 3)) ) )
    mismatch_penalty = 2
    bowtie2_mp = "--mp {pmax},{pmin}".format(pmax = str(mismatch_penalty), pmin = str(mismatch_penalty)  )
    match_bonus = 2
    bowtie2_ma = "--ma " + str(match_bonus)

    min_required_match_score = (int(five_p_seq_length) - allowed_mismatches) * match_bonus \
                               - (mismatch_penalty * allowed_mismatches)
    bowtie2_score_min = "--score-min L,{m},0".format(m = str(min_required_match_score) )

    bowtie2_5p_arguments = ' -a --local -N 1 -R 20 -i S,1,0.40 {L} --rdg 0,2 --rfg 0,2 {ma} {mp} {score_min} '.\
                            format(L = bowtie2_L,
                                   ma = bowtie2_ma,
                                   mp = bowtie2_mp,
                                   score_min = bowtie2_score_min)

    return bowtie2_5p_arguments