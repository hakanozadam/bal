
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


#################################################################
#################################################################

class StepError(Exception):
    pass

class InputError(StepError):
    pass

class ExecutableError(StepError):
    pass

class OutputError(StepError):
    pass
