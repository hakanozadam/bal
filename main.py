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


import os
import sys

##################################################################

def prep():
    try:
        import pysam
    except ImportError:
        print("Couldn't find the module pysam. Install the python3 package pysam before running bal.")
        exit(1)

##################################################################
def main():
    # Skip for now but later execute prep
    #prep()
    argument_1 = sys.argv[1]
    sys.argv = [ sys.argv[0] ] + sys.argv[2:]
    if argument_1 == 'align':
        print("Aligning... To be implemented")
    elif argument_1 == 'reference':
        print("Preparing refernce... To be implemented")
    else:
        print("Usage:\n./main.py align|reference arguments")


if __name__ == '__main__':
    main()