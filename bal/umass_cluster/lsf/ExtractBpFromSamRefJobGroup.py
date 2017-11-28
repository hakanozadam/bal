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

class ExtractBpFromSamRefJobGroupException(Exception):
   pass

class ExtractBpFromSamRefJobGroup(JobGroup.JobGroup):

   def __init__(self , input_directory, output_directory ,
                  individual_lib_directory = "",
                  name = "Extract_Bp_From_Sam_Ref" , run_time = 2,
                  executable = "", memory = 16000   ,
                  **kwargs ):

      self.input_extensions = ['sam']

      super().__init__( name = name, input_directory = input_directory , output_directory = output_directory ,
                        executable = executable, input_extensions = self.input_extensions ,
                        **kwargs)

      self.run_time        = run_time * 60 # Convert hours to minutes

      if self.run_time < 60:
          self.run_time = 60

      if memory < 16000 :
          self.memory = 16000
      else:
          self.memory = memory

      self.exclude_system_files()
      self.jobs = list()

      self.run_success_file_name = ".success"

      self.bp_candidates_bed_file = os.path.join(self.output_directory, 'bp_candidates.bed')
      self.bp_candidates_lib_bed_files_directory = os.path.join(self.output_directory,
                                                                'bp_candidates_library_bed_files')
      os.makedirs(self.bp_candidates_lib_bed_files_directory, exist_ok=True)
      self.individual_lib_directory = individual_lib_directory

   ###########################################################################

   def prepare_jobs(self):
      self.jobs = list()

      if len( self.input_files ) > 2:
          raise( ExtractBpFromSamRefJobGroupException( " Found too many sam files:"\
              + "\n".join(self.input_files)) )

      arguments = " -1 " + self.input_files[0]

      if len(self.input_files) == 2:
          arguments += " -2 " + self.input_files[1]

      arguments += " -o " + self.output_directory

      command = self.executable + arguments

      job = Job.Job(command = command ,
                                       output_directory = self.cluster_out_directory ,
                                       job_name = "bp_from_bed_ref",
                                       time_limit = self.run_time ,
                                       cores = 1 , memory = self.memory)
      self.jobs.append(job)

      if self.individual_lib_directory != "":
         library_directories = glob.glob( self.individual_lib_directory + "/*")
         library_names = map(os.path.basename, library_directories)
         self.individual_output_directory = os.path.join( self.output_directory, "individual_libraries")
         os.makedirs(self.individual_output_directory, exist_ok=True)

         for library_directory in library_directories:
            library_name = os.path.basename(library_directory)
            this_output_directory = os.path.join(self.individual_output_directory, library_name)
            os.makedirs(this_output_directory, exist_ok=True)

            if len(self.input_files) == 1:
               sam_file = os.path.join(library_directory, 'alignment.sam' )
               if not os.path.isfile(sam_file):
                  raise(FileNotFoundError("Could not find the sam file " + sam_file ))
               arguments = " -1 " + sam_file + " -o " + this_output_directory

            elif len(self.input_files) == 2:
               sam_file_1 = os.path.join(library_directory, 'alignment_1.sam' )
               sam_file_2 = os.path.join(library_directory, 'alignment_2.sam' )
               if not os.path.isfile(sam_file_1):
                  raise(FileNotFoundError("Could not find the sam file " + sam_file_1 ))
               if not os.path.isfile(sam_file_2):
                  raise(FileNotFoundError("Could not find the sam file " + sam_file_2 ))

               arguments  = " -1 " + sam_file_1 + " -2 " + sam_file_2 + " -o " + this_output_directory

            command = self.executable + arguments
            job = Job.Job(          command          = command ,
                                       output_directory = self.cluster_out_directory ,
                                       job_name         = "bp_from_bed_ref",
                                       time_limit       = self.run_time ,
                                       cores            = 1 ,
                                       memory           = self.memory)
            self.jobs.append(job)

   ###########################################################################

   ##########################################################################
   def post_run(self):
       self.error_message = ''

       if not os.path.isfile(self.bp_candidates_bed_file):
           self.error_message = "Could not find the bed file: " + self.bp_candidates_bed_file

       if self.individual_lib_directory != "":
          library_directories = glob.glob( self.individual_lib_directory + "/*")
          for library_directory in library_directories:
             library_name          = os.path.basename(library_directory)
             this_output_directory = os.path.join(self.individual_output_directory, library_name)
             bed_file              = os.path.join(this_output_directory, "bp_candidates.bed" )

             if not os.path.isfile(bed_file):
                self.error_message += "\nCould not find the bed file " + bed_file
             else:
                new_bed_file = os.path.join( self.bp_candidates_lib_bed_files_directory, library_name + ".bed"  )
                os.rename(bed_file, new_bed_file)

       super().post_run()
      
   ###########################################################################

###############################################################################

###############################################################################

def main():
   pass
   ####################################################################################

if __name__ == '__main__':
   main()
