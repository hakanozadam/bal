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
import datetime
import re
import glob
import subprocess
from functools import wraps
###############################################################################
import sys
script_directory = os.path.dirname(os.path.realpath(__file__))
sys.path.append( os.path.split( script_directory)[0]  )
from orca import Job
###############################################################################

def time_this(func):
    @wraps(func)
    def wrapper(obj, *args, **kwargs):
        setattr(obj,'start_time', datetime.datetime.now())
        result = func(obj, *args, **kwargs)
        setattr(obj,'finish_time', datetime.datetime.now())
        print('Running time is ' +\
            str( getattr(obj, 'finish_time') - getattr(obj, 'start_time') ) )
        return result
    return wrapper


###############################################################################
class JobGroupException(Exception):
   pass

class PostRunException(Exception):
   pass
#################################################################################
# Extend this class to particular classes such as Gzip, Bowtie and etc
# Implement the __init__() method annd don't forget to call the
# init method of the superclass but super()__init(****)
# Implement the prepare_jobs method and post_run mmethod and your script will be ready
# to be submitted

# Always run exclude_system_files before running the3 prepare jobs step!

class JobGroup:

   def __init__(self , name, input_directory, output_directory ,
                 executable , **kwargs ):

      ### Check the executable
      p = subprocess.Popen(["which " + executable], stdout=subprocess.PIPE, 
                                                  stderr=subprocess.PIPE,
                                                  shell = True)
      out_b , err_b = p.communicate()
      if p.returncode : 
         raise JobGroupException("Couldn't find the executable %s"%executable)
      else:
         self.executable = executable
         self.executable_arguments = kwargs.get("executable_arguments" , "")
      #### Check the input files
      self.input_extensions = kwargs.get("input_extensions" , list())
      self.file_selection_mode = kwargs.get("file_selection_mode" , "files")

      self.input_directory = os.path.abspath(input_directory)
      self.name = name

      if not os.path.isdir( input_directory):
         raise JobGroupException("Input directory %s doesn't exist!"%input_directory)

      input_files = list()
      if len(self.input_extensions) == 0:
         input_files = glob.glob(input_directory + "/*") 
      else:
         input_files = list()
         for extension in self.input_extensions:
            input_files += glob.glob(input_directory + "/*." + extension)

      self.input_files = list()
      if self.file_selection_mode == "directories" :
         input_files = glob.glob(self.input_directory + "/*")
         for file in input_files:
            if os.path.isdir( file  ):
               self.input_files.append(os.path.abspath(file))
      else: # File mode
         for file in input_files:
            if os.path.isfile( file  ):
               self.input_files.append(os.path.abspath(file))

      ## Done checking the input files
      
      self.output_directory = os.path.abspath(output_directory)
      os.makedirs(self.output_directory , exist_ok = True)
      self.cluster_out_directory = os.path.join(self.output_directory , ".cluster_out")
      os.makedirs(self.cluster_out_directory , exist_ok = True)
      self.lock_directory = os.path.join(self.cluster_out_directory , "lock")
      os.makedirs(self.lock_directory , exist_ok = True)
      ## Now we prepare the jobs to be submitted
      self.jobs = list()
  
      self.commands_log_file = os.path.join(self.cluster_out_directory,\
                                             "submitted_commands.log")
      self.run_error_file   = os.path.join(self.output_directory , 
                                        ".error.log")
      self.error_message    = ""
      self.run_success_file = os.path.join(self.output_directory ,\
                                        ".success.log")      
      self.do_cleanup = kwargs.get("do_cleanup" , False)
      self.report_run_error = kwargs.get("report_run_error" , True)

      self.make_report      = kwargs.get("make_report" , False)
      self.report_directory = os.path.join(self.output_directory , ".report")
      if not os.path.exists(self.report_directory) and self.make_report:
         os.makedirs(self.report_directory)
 
      self.make_latex      = kwargs.get("make_latex" , False)
      self.latex_directory = os.path.join(self.output_directory , ".latex")
      if not os.path.exists(self.latex_directory) and self.make_latex:
         os.makedirs(self.latex_directory)

      self.debug_mode = kwargs.get("debug_mode" , "")

      self.exclude_system_files()

      self.job_lock_files = dict()

   ############################################################################

   def get_file_name(self, file_path):
      if self.file_selection_mode == "files":
         return os.path.splitext(os.path.basename(file_path))[0]
      else :
         return os.path.basename(file_path)

   ############################################################################
   def prepare_jobs(self):
      self.jobs = list()
      for file in self.input_files:
         file_name = JobGroup.get_file_name(self, file)
         command   = self.executable + " " + self.executable_arguments +\
                     " -i " + self.input_directory + \
                     " -o " + self.output_directory
         this_job  = Job.Job(command = command , 
                             output_directory = self.cluster_out_directory ,
                             job_name = JobGroup.get_file_name(self , file) )
         self.jobs.append(this_job) 

   ############################################################################
   def submit_jobs(self):
      if os.path.exists(self.run_success_file):
         os.remove(self.run_success_file)
      if os.path.exists(self.run_error_file):
         os.remove(self.run_error_file)
      commands_log = open( self.commands_log_file , "w")
      lock_counter = 0

      for job in self.jobs:
         lock_counter += 1
         actual_command = job.command
         job.lock_file = os.path.join( self.lock_directory , 
                                  job.job_name + "_"  + str(lock_counter) + ".lock"  )
         lock_file_result = os.system("touch " + job.lock_file)
         if lock_file_result:
            raise JobGroupException("Couldn't create the lock file %s"%job.lock_file)
         job.command  += " ; sleep 5 ;  rm %s"%job.lock_file
         this_job_id = job.submit()
         self.job_lock_files[job.lock_file] = this_job_id
         ## We assume that there are less than 5 seconds between 
         ## the previous statement and this one
         if os.system("echo %s > %s"%(job.job_id , job.lock_file ) ):
            raise JobGroupException("Couldn't write to the lock file %s"%job.lock_file)
         job.command   = actual_command 
         commands_log.write( job.command + "\n" )
      commands_log.close()

   ############################################################################
   def get_active_job_ids(self):

       p = subprocess.Popen(['bjobs'], stdout=subprocess.PIPE,
                                                  stderr=subprocess.PIPE,
                                                  shell = True)
       out_b , err_b = p.communicate()
       out = out_b.decode("utf-8").rstrip("\n")
       err = err_b.decode("utf-8").rstrip("\n")

       if(err == "No unfinished job found"):
           job_ids = []
       else:
           job_lines = out.split("\n")[1:]
           job_ids = [line.strip().split(" ")[0] for line in job_lines ]
       return job_ids

   ###########################################################################
   # new version
   def wait_for_jobs_to_finish(self):
      os.system("sleep 10 ")
      print(self.name, "jobs:")
      while(True):
         exists_unfinished_jobs = False
         os.system("sleep 60 ")
         active_job_ids = list(set(self.get_active_job_ids() ) & set(self.job_lock_files.values() ) )
         print(" ", end= "\r")
         print( '\r' + str(len(active_job_ids)) + " jobs remaining in the LSF queue", end = "\r" )
         sys.stdout.flush()

         for job in self.jobs:

            if(os.path.exists(job.lock_file)):
               exists_unfinished_jobs = True
               this_job_id = self.job_lock_files[job.lock_file]
               if not (this_job_id in active_job_ids ):
                  print("WARNING: There may be a problem with the job having the id", this_job_id)
                  print("See the lsf out log file in the .cluster_out/out directory and .cluster_out/summary directory")
         if not exists_unfinished_jobs:
            break

   ############################################################################
   # old version
   # def wait_for_jobs_to_finish(self):
   #    while(True):
   #       exists_unfinished_jobs = False
   #       os.system("sleep 2 ")
   #       for job in self.jobs:
   #          if(os.path.exists(job.lock_file)):
   #             exists_unfinished_jobs = True
   #       if not exists_unfinished_jobs:
   #          break

   ############################################################################
   def verify_run(self):
      pass
 
   ###########################################################################

   def cleanup(self):
      pass
      #os.system("rm -r %s"%self.cluster_out_directory)

   ############################################################################
   def prepare_error_report(self):
      if self.error_message:
         error_file = open( self.run_error_file , "w" )
         error_file.write(self.error_message)
         error_file.close()
         raise(PostRunException(self.error_message))
      else: # SUCCESS
         if os.system("touch %s"%self.run_success_file):
            raise(JobGroupException("Couldn't create the file %s"%self.run_success_file))

   ############################################################################

   def prepare_report(self):
      pass

   ############################################################################

   def prepare_latex(self):
      pass

   ############################################################################
   def post_run(self):
      if self.make_report:
         self.prepare_report()
      if self.make_latex:
         self.prepare_latex()
      if self.do_cleanup:
         self.cleanup()
      if self.report_run_error:
         self.prepare_error_report()

   ############################################################################  
   # Excludes possible config and log files from input files
   # just in case
   def exclude_system_files(self):
      tmp_files = self.input_files
      self.input_files = list()
      restricted_list = [ os.path.abspath(self.run_error_file),
                          os.path.abspath(self.run_success_file),
                          os.path.abspath(self.cluster_out_directory),
                          os.path.abspath(self.report_directory),
                          os.path.abspath(self.latex_directory)
                        ]
      for file in tmp_files:
         if not os.path.abspath(file) in restricted_list:
            self.input_files.append(file)
      
      if len(self.input_files) == 0:
         raise(JobGroupException("No files found in the input directory %s"\
                                  %self.input_directory))
   
   ###########################################################################
   def remove_input_files(self):
      for file in self.input_files:
         os.remove(file)
   ###########################################################################
   def remove_output_files(self):
      pass
   ###########################################################################
   # This should be used to run and do post running operations
   def run(self):
      self.prepare_jobs()
      initial_directory = os.getcwd()
      os.chdir(self.output_directory)

      if self.debug_mode == "post_run":
         print("This module is running in post_run debug mode.\n" +
              "So no jobs are submitted. So just check the post_run files.")
      else:
         self.submit_jobs()
         self.wait_for_jobs_to_finish()

      os.chdir(initial_directory)
      self.post_run()

   ###########################################################################
   # Used when the class file is run as a stand alone script
   @time_this
   def run_in_main(self):
      if  os.path.isfile(self.run_success_file):
         print("This step has run before. Skipping...")
         return(0)
      self.prepare_jobs()
      initial_directory = os.getcwd()
      os.chdir(self.output_directory)

      if self.debug_mode == "post_run":
         print("This module is running in post_run debug mode.\n" +
              "So no jobs are submitted. So just check the post_run files.")
      else:    
         print("Running the jobs:")
         for job in self.jobs:
            print(job.command + "\n")
         JobGroup.submit_jobs(self)
         print("waiting for the jobs...")
         JobGroup.wait_for_jobs_to_finish(self)

      os.chdir(initial_directory)           
      try:
         self.post_run()
      except PostRunException as run_error:
         print(run_error.args[0])
      else:
         print("DONE: Output directory is %s"%self.output_directory)   

###############################################################################

   def merge_csv_files(self, input_list, output_file, delimeter = "," ):
       result_matrix = dict()
       header_list   = list()

       for file in input_list:
           with open(file , 'r') as csv_file:
               file_name = os.path.basename(file)
               header_list.append(file_name)
               for line in csv_file:
                   line_contents = line.strip().split(delimeter)
                   if len(line_contents) < 2:
                       continue
                   try:
                       result_matrix[line_contents[0]].append(line_contents[1])
                   except KeyError:
                       result_matrix[line_contents[0]] = [line_contents[1]]

       with open(output_file, "w") as out_csv:
            out_csv.write( delimeter +  delimeter.join(header_list) + "\n")
            for key in sorted( result_matrix.keys() ):
                out_csv.write(key + delimeter + delimeter.join(result_matrix[key]) + "\n" )



###############################################################################
def correct_command_line_arguments(parameter):
    import sys
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
###############################################################################
def main():
   print("This file contains the abstract class Job Group")

if __name__ == '__main__':
   main()
