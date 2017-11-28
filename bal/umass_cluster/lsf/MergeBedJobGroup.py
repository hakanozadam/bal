#!/bin/env python3

# AUTHOR:
#        Hakan Ozadam
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
from time import gmtime, strftime , localtime
import re
import shutil
from collections import namedtuple
import glob
import sys

###############################################################################

if __name__ == '__main__':
    script_directory = os.path.dirname(os.path.realpath(__file__))
    sys.path.append( os.path.split( script_directory)[0]  )
    import JobGroup , Job
else:
    from . import JobGroup, Job

class MergeBedJobGroupException(Exception):
   pass

class MergeBedJobGroup(JobGroup.JobGroup):

   def __init__(self , input_directory, output_directory ,
                  name = "BAL" , run_time = 2,
                  executable = "", memory = 4096   ,
                  single_file_list = "", # if single file is given, only the bed files in this file is merged
                  # Otherwise the default behavior is looking for the file lists
                  # in the directories inside input_directory
                  file_selection_mode = "directories" , **kwargs ):

      self.input_extensions = ['txt']
      self.single_file_list = single_file_list

      super().__init__( name = name, input_directory = input_directory , output_directory = output_directory ,
                        executable = executable, input_extensions = self.input_extensions ,
                        file_selection_mode = "directories",
                        **kwargs)

      self.run_time        = run_time * 60 # Convert hours to minutes

      if self.run_time < 60:
          self.run_time = 60

      self.memory = 4096 if memory else memory

      self.exclude_system_files()
      self.jobs = list()
      self.merged_files = os.path.join(self.output_directory, "candidate_bp_files_list.txt")
      self.list_of_merged_files = list()

      self.run_success_file_name = ".success"

      self.bp_candidates_bed_file = os.path.join(self.output_directory, 'bp_candidates.bed')

   ###########################################################################

   def prepare_jobs(self):

      self.jobs = list()

      if self.single_file_list:
          if not os.path.isfile(self.single_file_list):
              raise FileNotFoundError("Couldn't find the bed file list ", self.single_file_list)

          command = '{executable} -i {bed_file_list} -o {merged_bed_file}'.\
                    format(executable = self.executable ,
                           bed_file_list = self.single_file_list,
                           merged_bed_file = self.bp_candidates_bed_file
                           )
          job = Job.Job(command = command ,
                                       output_directory = self.cluster_out_directory ,
                                       job_name = "final_bed_merge",
                                       time_limit = 120 , cores = 1 , memory = 4096)
          self.jobs.append(job)
          return 0

      for directory in self.input_files:
         output_directory = os.path.join( self.output_directory ,os.path.basename(directory) )
         os.makedirs(output_directory, exist_ok = True)
         lib_bp_bed_files_list = os.path.join(directory,
                                           "candidate_bp_files_list.txt")
         merged_bed_file = os.path.join(output_directory, 'candidate_branchpoints.bed')
         self.list_of_merged_files.append(merged_bed_file)

         if not os.path.isfile(lib_bp_bed_files_list):
            raise FileNotFoundError("Could not find the bed file list", lib_bp_bed_files_list)

         command = '{executable} -i {bed_file_list} -o {merged_bed_file}'.\
                    format(executable = self.executable ,
                           bed_file_list = lib_bp_bed_files_list,
                           merged_bed_file = merged_bed_file
                           )

         job = Job.Job(command = command ,
                                       output_directory = self.cluster_out_directory ,
                                       job_name = os.path.basename(directory),
                                       time_limit = 120 , cores = 1 , memory = 4096)
         self.jobs.append(job)


      with open(self.merged_files, 'w') as list_stream:
        for bed_file in self.list_of_merged_files:
            print(bed_file, file = list_stream)

      

   ###########################################################################

   ##########################################################################
   def post_run(self):
       self.error_message = ""
       missing_files = list()

       if not self.single_file_list:
           if not os.path.isfile(self.merged_files):
               missing_files.append(self.merged_files)
           else:
               with open(self.merged_files, 'r') as merged_files_stream:
                   for line in merged_files_stream:
                       this_file = line.strip()
                       if not os.path.isfile(this_file):
                           missing_files.append( this_file )
       else:
           if not os.path.isfile(self.bp_candidates_bed_file):
               missing_files.append(self.bp_candidates_bed_file)

       if len(missing_files) > 0 :
           raise (FileNotFoundError("Error: The following files are missing\n" +
                 "\n".join(missing_files) ) )

       super().post_run()
      
   ###########################################################################

###############################################################################

###############################################################################

def main():
   parser = argparse.ArgumentParser(description=
   '''
   Bed File Merger
   ''')
   parser.add_argument("-i" ,
                       help = "Input directory.." ,
                       required = True ,
                       metavar = "input directory" ,
                       type = str)
   parser.add_argument("-o" ,
                       help = "Output directory" ,
                       required = True ,
                       metavar = "output directory",
                       type = str)
   args = parser.parse_args(JobGroup.correct_command_line_arguments("-a"))

   jobGroup = MergeBedJobGroup(
                 input_directory  = args.i ,
                 output_directory = args.o ,
                  )

   jobGroup.run_in_main()
   ####################################################################################

if __name__ == '__main__':
   main()
