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
import subprocess
import re
#################################################################

class JobException(Exception):
   pass

class Job:
   valid_cluster_types = ('lsf', )
   def __init__(self , command , output_directory , **kwargs):
      self.cluster_type       = kwargs.get('cluster_type' , 'lsf')
      if not self.cluster_type in Job.valid_cluster_types:
         raise JobException('%s is not a valid type. Valid types are %s'\
                %( self.cluster_type ,\
                    ','.join(Job.valid_cluster_types) ))

      self.command          = command

      self.output_directory    = os.path.abspath(output_directory)
      self.out_directory       = os.path.join(self.output_directory , "out")
      self.error_directory     = os.path.join(self.output_directory , "error")
      self.summary_directory   = os.path.join(self.output_directory , "summary")

      if not os.path.exists(self.output_directory):
         os.mkdir(self.output_directory)
      if not os.path.exists(self.out_directory):
         os.mkdir(self.out_directory)
      if not os.path.exists(self.error_directory):
         os.mkdir(self.error_directory)
      if not os.path.exists(self.summary_directory):
         os.mkdir(self.summary_directory)

      self.queue            = kwargs.get('queue' , 'short')
      self.time_limit       = kwargs.get('time_limit' , 230)
      self.memory           = kwargs.get('memory' , 4096)
      self.cores            = kwargs.get('cores' , 1)
      self.parameters       = kwargs.get('parameters' , '')
      self.job_name         = kwargs.get('job_name' ,\
                                os.path.split(self.output_directory.rstrip("/"))[1])
      if not self.job_name:
         self.job_name = os.path.split(self.output_directory.rstrip("/"))[1]
      
      self.submission_time = localtime()   
      time_stamp           = strftime("%m_%d_%Y_%H_%M_%S", self.submission_time)
      out_file_name        = r"%J" + "_%s_%s.out"%(self.job_name , time_stamp)
      self.out_file        = os.path.join(self.out_directory , out_file_name)
      error_file_name      = r"%J" + "_%s_%s.error"%(self.job_name , time_stamp)
      self.error_file      = os.path.join(self.error_directory , error_file_name)
      self.lock_file       = ""

      self.summary_file_name    = "_%s_%s.summary"%(self.job_name , time_stamp)

   ##########################################################

   def lsf_submit(self):
      if self.time_limit > 230:
         self.queue = 'long'

      r_string = "\"select[hname!=\'ghpcc-sgi\'] rusage[mem=%s] span[hosts=1]\""%self.memory
      lsf_command = "bsub -q {queue} -W {time_limit} -n {cores} -R \"select[idm]\" "\
                    "-R {R_string} -J {job_name} "\
                    "-e {error_file} -o {out_file} {parameters} \"{command}\"".format(queue=self.queue ,
                       time_limit = self.time_limit , cores = self.cores ,
                       command = self.command , R_string = r_string , 
                       job_name = self.job_name , error_file = self.error_file , 
                       out_file = self.out_file , parameters = self.parameters
                   )
      p = subprocess.Popen([lsf_command], stdout=subprocess.PIPE, 
                                                  stderr=subprocess.PIPE,
                                                  shell = True)
      out_b , err_b = p.communicate()
      out = out_b.decode("utf-8").rstrip("\n")
      err = err_b.decode("utf-8").rstrip("\n")
      if( err_b ):
         raise JobException(err)      
      submission_output = re.search( r"<([^>]*)>" , out)
      if not submission_output:
         raise JobException("There was a problem in getting the job number")
      self.job_id = submission_output.group(1)
      
      # Make Summmary
      submission_time = strftime("%m %d %Y , %H : %M : %S", self.submission_time)

      self.summary = "------      JOB SUBMISSION SUMMARY    -------\n" +\
                "Command            = %s\n\n"%self.command +\
                "Job Name           = %s\n"%self.job_name +\
                "Job Id             = %s\n"%self.job_id +\
                "Queue              = %s\n"%self.queue +\
                "Submission Time    = %s\n"%submission_time +\
                "Memory Requested   = %s\n"%self.memory +\
                "Time Limit         = %s\n"%self.time_limit +\
                "Cores              = %s\n"%self.cores +\
                "Additional Params  = %s\n\n"%self.parameters +\
                "bsub command       = %s\n"%lsf_command +\
                "--------------------------------------------------------\n"     
      return(self.job_id) 

   def submit(self):
      if(self.cluster_type == 'lsf'):
         this_job_id            = Job.lsf_submit(self)
         self.summary_file_name = str(this_job_id) + self.summary_file_name
         self.summary_file      = os.path.join(self.summary_directory , 
                                               self.summary_file_name)
         summary_file           = open(self.summary_file , "w")
         summary_file.write(self.summary)
         summary_file.close()
         return this_job_id

###################################################################

def main():
   parser = argparse.ArgumentParser(description=
   '''
HPC Cluster Submission classes and functions. This script submits the given
commannd to the given cluster system. Note that many of the arguments 
are optional their documentation can be found here.
Put the command argument (-c) inside quotes!
Put the -p argument inside quotes and if the first char is -
then put a space after the first quote. 
   ''')
   parser.add_argument("-c",
                       metavar = "command" ,
                       help = "Command to be submitted to the cluster" , 
                       required = True , 
                       type = str)
   parser.add_argument("-o",
                       metavar = "output directory" ,
                       help = "Output directory for cluster output files" , 
                       required = False ,
                       default = os.path.join(os.environ["HOME"] , "lsf_logs"),
                       type = str)
   parser.add_argument("-s",
                       metavar = "scheduler" ,
                       help = "Job scheduler system type. Default: lsf" , 
                       required = False ,
                       choices = ["lsf" , ] ,
                       default = "lsf",
                       type = str)
   parser.add_argument("-m",
                       metavar = "memory" ,
                       help = "Memory limit of the job in MB. Default: 4096" ,
                       required = False , 
                       default = 4096,
                       type = int)
   parser.add_argument("-t",
                       metavar = "time" ,
                       help = "Time limit of the job in minutes. Default: 230" ,
                       required = False ,
                       default = 230,
                       type = int)
   parser.add_argument("-r",
                       metavar = "cores" ,
                       help = "Number of cores to handle the job. Default: 1" ,
                       required = False ,
                       default = 1,
                       type = int)
   parser.add_argument("-n",
                       metavar = "name" ,
                       help = "Job name" ,
                       required = False ,
                       default = "" ,
                       type = str)
   parser.add_argument("-p",
                       metavar = "addnl. pars." ,
                       help = "Additional cluster parameters. Place this in quotes and put one space before the first -" ,
                       required = False ,
                       default = "" ,
                       type = str)
 
   arguments = parser.parse_args()
  
   this_job = Job(command          = arguments.c , 
                 output_directory  = arguments.o , 
                 cluster_type      = arguments.s ,
                 memory            = arguments.m ,
                 time_limit        = arguments.t ,
                 cores             = arguments.r , 
                 job_name          = arguments.n , 
                 parameters        = arguments.p ) 
   job_id = this_job.submit();  
   print(this_job.summary)


if __name__ == '__main__':
   main()


