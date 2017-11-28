#!/bin/env python3

# AUTHORS:
#        Hakan Ozadam
#
#        Moore Laboratory
#        UMASS Medical School / HHMI
#        RNA Therapeutics Institute
#        Albert Sherman Center, ASC4-1009
#        368 Plantation Street
#        Worcester, MA 01605
#        USA
#
####################################################################

import argparse
import os
import operator
import sys

#####################################################################

def handle_arguments():
   parser = argparse.ArgumentParser(description=
   '''
   spombe gff file is not compatible with bal as it does not have gene_id/transcript_id fields.
   This scripts fixes that issue by outputting a compatible gtf file.
   ''')
   parser.add_argument("-i" ,
                       help = "Input gtf/gff file" ,
                       required = True ,
                       metavar = "input_gtf_file" ,
                       type = str)
   parser.add_argument("-o" ,
                       help = "Output gtf file" ,
                       required = True ,
                       metavar = "output_gtf_file" ,
                       type = str)

   arguments = parser.parse_args()
   return arguments

#####################################################################

def get_label_contents(label_string, imposed_identifier):
    contents = label_string.split(':')
    if len(contents) < 2 :
        raise(IOError("Error in splitting the label string :" + label_string ))
    if contents[0] != imposed_identifier:
        raise(IOError('Wrong identifier: ' + contents[0]))
    if contents[1] == '':
        raise(IOError('Empty contents in' + label_string))
    return(contents[1])

#####################################################################

def parse_attributes(attribute_string):
   result = dict()
   contents = attribute_string.split(";")
   if len(contents) < 1:
      print("Error in gtf file format:", attribute_string)
      exit(1)
   for c in contents:
      sides = c.split('=')
      if len(sides) < 2:
         print('Error in siedes')
         exit(2)
      result[sides[0]] = sides[1]
   return result

#####################################################################

def get_ID_type(parsed_att_dict):

    entry_type = ''

    for key, value in parsed_att_dict.items():
        if key == 'ID':
            contents = value.split(':')
            if len(contents) < 2:
                raise(IOError('Problem witd ID contents in the dictionary' + str(parsed_att_dict)))
            entry_type = contents[0]

    if entry_type == '':
        raise(IOError('No ID information found on ' + str(parsed_att_dict)))

    if entry_type not in ('gene', 'transcript', 'CDS'):
        raise(IOError('Invalid entry_type ' + entry_type + ' in ' + str(parsed_att_dict)))

    return entry_type

#####################################################################

def get_gene_id_from_gene(parsed_att_dict):
    gene_id = ''
    for key, value in parsed_att_dict.items():
        if key == 'ID':
            gene_id = get_label_contents(value, 'gene')
    if gene_id == '':
        raise(IOError('Error: Could not find gene_id in the att dict' + str()))

    return(gene_id)

#####################################################################

def get_transcript_labels_from_trans_attribute(parsed_att_dict):

    transcript_id, gene_id = '', ''

    for key, value in parsed_att_dict.items():
        if key == 'ID':
            transcript_id = get_label_contents(value, 'transcript')
        if key == 'Parent':
            gene_id = get_label_contents(value, 'gene')

    if transcript_id == '' or gene_id == '':
        raise (IOError('gene_id or transcript_id is empty in the dict:' + str(parsed_att_dict) ))

    result = dict()
    result['gene_id'] = gene_id
    result['transcript_id'] = transcript_id

    return result

#####################################################################

def get_exon_attributes(parsed_att_dict, genes_of_transcripts):
    result = dict()
    gene_id = ''
    transcript_id = ''

    for key, value in parsed_att_dict.items():
        if key == 'Parent':
            transcript_id = get_label_contents(value, 'transcript')

            gene_id = genes_of_transcripts.get(transcript_id, None)

            if gene_id == None:
                raise(IOError("Could not find the gene id for the entry "\
                              + str(parsed_att_dict)))

    if gene_id == '' or transcript_id == '':
        raise(IOError("gene_id or transcript_id is empty in the attributes"\
                      + str(parsed_att_dict)))
    result['gene_id'] = gene_id
    result['transcript_id'] = transcript_id

    return result

#####################################################################

def main():
   '''
   We assume that a transcript -> gene relation is given in the gtf file
   The given gtf file might not be sorted properly so in the first pass
   we find the transcript-> gene relation
   in the second pass we output the result using the transcript-> gene information
   '''
   arguments = handle_arguments()
   genes_of_transcripts = dict()

   # # First pass: extract the transcript -> gene relation
   # with GtfFile(arguments.i) as input_gtf_stream:
   #      for input_entry in input_gtf_stream:
   #        parsed_attributes = parse_attributes(input_entry.attribute)
   #        #print(parsed_attributes)
   #        if input_entry.feature == 'transcript':
   #           transcript_labels = get_transcript_labels_from_trans_attribute(parsed_attributes)
   #           genes_of_transcripts[ transcript_labels['transcript_id'] ] = transcript_labels['gene_id']
   #        else:
   #            continue


   # Second pass: determine the labels (gene_id, transcript_id) of the entries
   with GtfFile(arguments.i) as input_gtf_stream,\
        open(arguments.o, "w") as output_stream:

      for input_entry in input_gtf_stream:
          parsed_attributes = parse_attributes(input_entry.attribute)
          #print(parsed_attributes)
          if input_entry.feature == 'exon':
             exon_labels = get_exon_attributes(parsed_attributes, genes_of_transcripts)

             input_entry.attribute = "gene_id \"{gene_id}\"; transcript_id \"{transcript_id}\";"\
                                     .format(gene_id = exon_labels['gene_id'],
                                             transcript_id = exon_labels['transcript_id'])
             print(input_entry, file = output_stream)
          elif input_entry.feature == 'gene':
             gene_id = get_gene_id_from_gene(parsed_attributes)
             input_entry.attribute = "gene_id " + "\""+ gene_id +"\";"
             print(input_entry, file = output_stream)
          elif input_entry.feature == 'transcript':
             transcript_labels = get_transcript_labels_from_trans_attribute(parsed_attributes)
             genes_of_transcripts[ transcript_labels['transcript_id'] ] = transcript_labels['gene_id']
             input_entry.attribute = "gene_id \"{gene}\"; transcript_id \"{transcript}\";"\
                                       .format( gene = transcript_labels['gene_id'],
                                                transcript = transcript_labels['transcript_id'])
             print(input_entry, file = output_stream)
          elif input_entry.feature == 'CDS':
             exon_labels = get_exon_attributes(parsed_attributes, genes_of_transcripts)

             input_entry.attribute = "gene_id \"{gene_id}\"; transcript_id \"{transcript_id}\";"\
                                     .format(gene_id = exon_labels['gene_id'],
                                             transcript_id = exon_labels['transcript_id'])
             print(input_entry, file = output_stream)
          elif input_entry.feature  in ('repeat_region', 'biological_region') :
              continue
          else:
             entry_type = get_ID_type( parsed_attributes )
             extra_label = input_entry.feature


             if entry_type == 'transcript':
                 transcript_labels = get_transcript_labels_from_trans_attribute(parsed_attributes)
                 genes_of_transcripts[ transcript_labels['transcript_id'] ] = transcript_labels['gene_id']
                 input_entry.attribute = "gene_id \"{gene}\"; transcript_id \"{transcript}\"; " \
                                         "type_label \"{type_label}\";"\
                                         .format( gene = transcript_labels['gene_id'],
                                                  transcript = transcript_labels['transcript_id'],
                                                  type_label = extra_label )
             elif entry_type == 'gene':
                 gene_id = get_gene_id_from_gene(parsed_attributes)
                 input_entry.attribute = "gene_id " + "\""+ gene_id +"\"; " + "type_label \"" + extra_label + "\""
             else:
                 continue

             input_entry.feature = entry_type
             print(input_entry, file = output_stream)



#####################################################################

script_directory = os.path.abspath(os.path.dirname(os.path.realpath(__file__)) )
bal_dir = os.path.split( os.path.split( script_directory)[0] )[0]

if __name__ == '__main__':

    sys.path.append( bal_dir  )
    from genomic_io.gtf import GtfFile, GtfEntry, check_gtf_validity
    
    main()
else:
    exit(1)

####################################################################
