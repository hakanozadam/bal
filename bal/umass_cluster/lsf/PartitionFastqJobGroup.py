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


###############################################################################
import sys
script_directory = os.path.dirname(os.path.realpath(__file__))
sys.path.append( os.path.split( script_directory)[0]  )
from lsf import JobGroup , Job

class PartitionFastqJobGroupException(Exception):
   pass

class PartitionFastqJobGroup(JobGroup.JobGroup):

   def __init__(self , input_directory, output_directory , reads_per_file,
                  alignment_mode, run_time = 1,
                  name = "Partition Fastq Files",
                  executable = "",
                  memory = 4096   , **kwargs ):

      self.input_extensions = ["fastq" , "fq"]

      super().__init__( name = name, input_directory = input_directory , output_directory = output_directory ,
                        executable = executable, input_extensions = self.input_extensions ,
                        **kwargs)


      self.alignment_mode  = alignment_mode
      self.run_threads     = 1
      self.run_time        = run_time * 60 # Convert hours to minutes

      if self.run_time < 230:
          self.run_time = 230
      if self.run_threads < 1:
          self.run_threads = 1

      self.reads_per_file = reads_per_file

      self.memory = memory if memory else 4096

      self.exclude_system_files()
      self.jobs = list()

      self._arrange_input_files()

      self.output_files = list()

      for file in self.arranged_input_files:
          lib_directory = os.path.join(self.output_directory , file.name)
          os.makedirs(lib_directory, exist_ok = True)

   ###########################################################################

   def prepare_jobs(self):

      self.jobs = list()

      if self.alignment_mode == "single":
          for file in self.arranged_input_files:
              input_fastq              = file.strand_1_file
              output_directory         = os.path.join(self.output_directory , file.name)
              command = "{executable} -i {input} -o {output_directory} -n {reads_per_file} "\
                                     .format(executable       = self.executable ,\
                                             output_directory = output_directory,\
                                             input            = input_fastq,
                                             reads_per_file   = self.reads_per_file)

              job = Job.Job(command = command ,
                                    output_directory = self.cluster_out_directory ,
                                    job_name = file.name,
                                    time_limit = self.run_time , cores = self.run_threads , memory = self.memory)
              self.jobs.append(job)
      else:
          for file in self.arranged_input_files:
              input_1_fastq              = file.strand_1_file
              input_2_fastq              = file.strand_2_file
              output_directory      = os.path.join(self.output_directory , file.name)

              command_1 = "{executable} -i {input} -o {output_directory} -n {reads_per_file} "\
                                     .format(executable       = self.executable ,\
                                             output_directory = output_directory,\
                                             input            = input_1_fastq,
                                             reads_per_file   = self.reads_per_file)

              command_2 = "{executable} -i {input} -o {output_directory} -n {reads_per_file} "\
                                     .format(executable       = self.executable ,\
                                             output_directory = output_directory,\
                                             input            = input_2_fastq,
                                             reads_per_file   = self.reads_per_file)

              job = Job.Job(command = command_1 ,
                                    output_directory = self.cluster_out_directory ,
                                    job_name = file.name + "_1",
                                    time_limit = self.run_time , cores = self.run_threads , memory = self.memory)
              self.jobs.append(job)
              job = Job.Job(command = command_2 ,
                                    output_directory = self.cluster_out_directory ,
                                    job_name = file.name + "_2",
                                    time_limit = self.run_time , cores = self.run_threads , memory = self.memory)
              self.jobs.append(job)

   ###########################################################################

   ##########################################################################
   def post_run(self):
       self.error_message = ""
       missing_files = list()

       super().post_run()
      
   ###########################################################################

   def _get_file_2(self, lib_name):
       ext2 = re.compile(r"(?P<ext2>.*)[._]2\.f(ast)?q$" , flags = re.IGNORECASE)
       for file in self.input_files:
           file_name = os.path.basename(file)
           ext2_search = ext2.search(file_name)
           if ext2_search and ext2_search.group("ext2") == lib_name:
               return file
       return ""

   ###########################################################################

   def _arrange_input_files(self):

       Library = namedtuple("Library" , ["name" , "strand_1_file" , "strand_2_file"])
       self.arranged_input_files = list()
       ext1 = re.compile(r"(?P<ext1>.*)[._]1\.f(ast)?q$" , flags = re.IGNORECASE)

       for file in sorted(self.input_files):
           file_name = os.path.basename(file)

           if self.alignment_mode == "single":
               self.arranged_input_files.append( Library(self.get_file_name(file), file, "") )
           else:
               ext1_search = ext1.search(file_name)
               if ext1_search:
                   lib_name = ext1_search.group("ext1")
                   file2 = self._get_file_2(lib_name)
                   if not file2:
                       raise(BalJobGroupException("Couldn't find the strand 2 file for %s"%file))
                   self.arranged_input_files.append(Library(lib_name, file, file2) )

###############################################################################
###############################################################################

def main():
   parser = argparse.ArgumentParser(description=
   '''
   BAL Pipeline For Umass Cluster
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
   parser.add_argument("-m" ,
                       help = "Alignment Mode" ,
                       required = True ,
                       metavar = "Mode",
                       choices = ['single' , 'paired'],
                       type = str)
   parser.add_argument("-n" ,
                       help = "Number of reads per file" ,
                       required = True ,
                       metavar = "Reads per file",
                       type = int)
   parser.add_argument("-t" ,
                       help = "Run time in HOURS" ,
                       required = False ,
                       metavar = "run_time",
                       default = 1,
                       type = int)
   args = parser.parse_args(JobGroup.correct_command_line_arguments("-a"))

   if not args.t:
       run_time = 1
   else:
       run_time = args.t

   jobGroup = PartitionFastqJobGroup(
                 input_directory  = args.i ,
                 output_directory = args.o ,
                 alignment_mode   = args.m ,
                 run_time         = args.t ,
                 reads_per_file   = args.n ,
                 memory           = 4096
                  )

   jobGroup.run_in_main()

######################################################################################
if __name__ == '__main__':
   main()
