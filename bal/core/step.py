
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

from abc import abstractmethod, ABCMeta
import subprocess
import os
import datetime

from .exceptions import *

#################################################################
#################################################################

class Step(metaclass = ABCMeta):

    #####################################################################
    def __init__(self, name, input_files, output_directory, executable = '' , executable_arguments = ''):
        self.name = name
        if not name:
            raise StepError("No name given!")

        self.input_files          = input_files
        self.output_directory     = os.path.abspath(output_directory)

        self.log_contents = list()
        self.log_file     = os.path.join(self.output_directory, self.name + ".bal.log")
        self.success_file = os.path.join(self.output_directory, "success.bal.log")
        self.failure_file = os.path.join(self.output_directory, "failure.bal.log")

        self.executable           = executable
        self.executable_arguments = executable_arguments

        self.command = ''
        self.module = ''

        self.error_messages = list()

    #######################################################################
    def __del__(self):
        pass

    #######################################################################

    @abstractmethod
    def prepare(self):
        pass

    ########################################################################

    # Since this function wil be overridden if there is a module run
    # we set it to fail to remind ours.elves it must be overridden in case of a module run.
    def _module_run(self):
        subprocess.call('touch ' + self.failure_file , shell=True )
        raise(StepError("You have to override the method _module_run."))

    #########################################################################

    def _executable_run(self):
        p = subprocess.Popen([self.command + self.executable_arguments], stdout=subprocess.PIPE,
                                                  stderr=subprocess.PIPE,
                                                  shell = True)
        std_out , std_err = p.communicate()
        self.std_out = std_out.decode("utf-8").rstrip()
        self.std_err = std_err.decode("utf-8").rstrip()
        self.returncode = p.returncode

    #########################################################################

    def _only_run(self):
        if self.command:
            return self._executable_run()
        else:
            return self._module_run()

    #########################################################################

    def run(self):
        if os.path.exists(self.success_file):
            os.remove(self.success_file)
        if os.path.exists(self.failure_file):
            os.remove(self.failure_file)
        os.makedirs(self.output_directory , exist_ok = True)
        self.start_time = datetime.datetime.now()
        self.prepare()
        self._only_run()
        self.post_run()
        self.end_time = datetime.datetime.now()
        self.write_log()

    ##########################################################################

    def post_run(self):
        if self.command:
            if self.returncode:
                subprocess.call('touch ' + self.failure_file , shell=True )
            else:
                subprocess.call('touch ' + self.success_file , shell=True )

    ##########################################################################

    @property
    def did_run(self):
        if self.did_fail or self.did_success:
            return True
        else:
            return False

    ##########################################################################

    @property
    def did_fail(self):
        if os.path.isfile(self.failure_file):
            return True
        else:
            return False

    ##########################################################################

    @property
    def did_success(self):
        if os.path.isfile(self.success_file):
            return True
        else:
            return False

    ###########################################################################

    @property
    def input_files(self):
        if hasattr(self, '_input_files'):
            return self._input_files
        else:
            return list()

    ###########################################################################

    @input_files.setter
    def input_files(self, files):
        self._input_files = list()
        missing_files = list()
        for file in files:
            file_path = os.path.abspath(file)
            if not os.path.isfile(file_path):
                missing_files.append(file_path)
            self._input_files.append(file_path)
        if len(missing_files) > 0:
            raise(InputError("Error: The following files couldn't be found\n{files}".\
                format(files = "\n".join(missing_files))))

    ###########################################################################

    def cleanup(self):
        pass

    ############################################################################

    @property
    def report(self):
        pass

    #############################################################################

    @property
    def running_time(self):
        return int( (self.end_time - self.start_time).seconds )

    #############################################################################

    @property
    def executable(self):
        if hasattr(self, '_executable'):
            return self._executable
        else:
            return ''

    #############################################################################

    @executable.setter
    def executable(self, executable):
        if executable == '':
            self._executable = ''
            return

        p = subprocess.Popen(["which " + executable], stdout=subprocess.PIPE,
                                                  stderr=subprocess.PIPE,
                                                  shell = True)
        p.communicate()

        if p.returncode:
            subprocess.call('touch ' + self.failure_file , shell=True )
            raise ExecutableError("Error: Couldn't find the executable {executable}.".format(executable = executable))

        else:
            self._executable = executable

    #############################################################################

    def write_log(self):
        time_format = "%m %d %Y , %H:%M:%S"
        run_result = 'Did not run'
        if self.did_success:
            run_result = 'Success'
        elif self.did_fail:
            run_result = 'Fail'

        with open(self.log_file, 'w') as log_stream:
            print("Step Name: " + self.name, file = log_stream)
            print("Result: ", run_result, file = log_stream )
            print("Start Time: " + str(self.start_time.strftime(time_format)) , file = log_stream)
            print("End Time: " + str(self.end_time.strftime(time_format)) , file = log_stream )
            print("Running Time: " + str( self.running_time ) + ' seconds' ,
                  file = log_stream)
            if(self.command):
                print("Command: " + self.command, file = log_stream)
                print("Return Code: " + str(self.returncode), file = log_stream)
                print("Std Out:\n" + self.std_out, file = log_stream)
                print("Std Error:\n" + self.std_err, file = log_stream)
            else:
                print("Module: " + self.module, file = log_stream )
            if( len(self.log_contents) > 0 ):
                print('-------------------', file = log_stream)
                for log_entry in self.log_contents:
                    print( log_entry, file = log_stream)
            if len(self.error_messages) > 0 :
                error_message = "\n".join(self.error_messages)
                print('-------------------', file = log_stream)
                print('Error Message:', file = log_stream)
                print(error_message, file = log_stream)


    #############################################################################

    def file_list_existence_check(file_list):
        missing_files = list()
        for f in file_list:
            if not os.path.exists(f):
                missing_files.append(f)
        return missing_files

    ##############################################################################

    def run_parallel_subprocess(self, n):
        active_processes = list()

        if len(self.commands) <= n:
            for c in self.commands:
                active_processes.append(subprocess.Popen([c + self.executable_arguments], stdout= subprocess.PIPE,
                                                      stderr=subprocess.PIPE, close_fds=True,
                                                      shell = True))
        else:
            for c in self.commands[0:n]:
                active_processes.append(subprocess.Popen([c + self.executable_arguments], stdout= subprocess.PIPE,
                                                      stderr=subprocess.PIPE, close_fds=True,
                                                      shell = True))

            for c in self.commands[n:]:
                process_removed = False
                while not process_removed:
                    for p in active_processes:
                        if p.poll() is not None:
                            p.stdout.close()
                            process_removed = True
                            active_processes.remove(p)
                            active_processes.append(subprocess.Popen([c + self.executable_arguments],
                                                                     stdout  = subprocess.PIPE,
                                                                     stderr  = subprocess.PIPE, close_fds=True,
                                                                     shell   = True))
                            break


        #### Wait for the remaining processes
        while active_processes:
                for p in active_processes:
                    if p.poll() is not None:
                        p.stdout.close()
                        active_processes.remove(p)

    ###################################################################################
