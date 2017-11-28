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

class BranchpointRefJobGroupException(Exception):
   pass

class BranchpointRefJobGroup(JobGroup.JobGroup):

   def __init__(self , input_directory, output_directory ,
                  bed_file , fasta_file,
                  executable,
                  number_of_nucleotides = 100,
                  name = "Ref_Branchpoint" , run_time = 3,
                  memory = 8000   ,
                  **kwargs ):

      self.input_extensions = []

      super().__init__( name = name, input_directory = input_directory , output_directory = output_directory ,
                        executable = executable, input_extensions = self.input_extensions ,
                        **kwargs)

      self.bed_file   = os.path.abspath(bed_file)
      self.fasta_file = os.path.abspath(fasta_file)

      self.run_time        = run_time * 60 # Convert hours to minutes

      if self.run_time < 60:
          self.run_time = 60
      self.exclude_system_files()
      self.jobs = list()
      self.merged_files = os.path.join(self.output_directory, "candidate_bp_files_list.txt")
      self.list_of_merged_files = list()

      self.number_of_nucleotides = number_of_nucleotides
      if self.number_of_nucleotides < 50 :
          self.number_of_nucleotides = 50

      self.run_success_file_name = ".success"
      self.memory = memory

      self.bp_candidates_bed_file = os.path.join(self.output_directory, 'bp_candidates.bed')

   ###########################################################################

   def prepare_jobs(self):
       command = "{executable} -f {fasta_file} -b {bed_file} " \
                 "-o {output_directory} -n {nucs}".format(
                     executable = self.executable,
                     fasta_file = self.fasta_file,
                     bed_file   = self.bed_file,
                     output_directory = self.output_directory,
                     nucs = self.number_of_nucleotides
                  )
       this_job = Job.Job(command = command ,
                                       output_directory = self.cluster_out_directory ,
                                       job_name = self.name,
                                       time_limit = self.run_time , cores = 1 ,
                                       memory = self.memory)
       self.jobs.append(this_job)

   ###########################################################################

   ##########################################################################
   def post_run(self):
       missing_files = list()
       suffixes = ('.1.bt2', '.2.bt2', '.3.bt2', '.4.bt2', '.rev.1.bt2',  '.rev.2.bt2')
       ref_base = os.path.join(self.output_directory,'bp')
       for s in suffixes:
           this_file = ref_base + s
           if not os.path.isfile(this_file):
               missing_files.append(s)

       if len(missing_files) > 0:
           self.error_message = "The following files could not be found:" +\
               "\n".join(missing_files)
           raise(FileNotFoundError(self.error_message))

       super().post_run()



      
   ###########################################################################

###############################################################################

###############################################################################

