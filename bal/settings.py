###### BAL GENERAL SETTINGS ################
############################################
# Import this script in the other scripts
# to use these system-wide settings
############################################

settings = dict()

settings['field_separator']       = '__'
settings['fastq_field_separator'] = '___'

settings['max_N_in_intron'] = 5 # Maximum number of non nucleotide letters in the 5p intron sequence

### Paths are relative to ref_folder
### where ref_folder is given in the -x parameter
settings['five_prime_intron_reference_directory'] = 'five_prime_intron_reference'
settings['five_prime_intron_reference_name']      = 'intron'
settings['bowtie2_genome_directory']              = 'bowtie2_genome'
settings['bowtie2_gene_directory']                = 'bowtie2_gene_ref'
settings['bowtie2_gene_sub_directory']            = 'bt2_references'
settings['gene_sequence_directory']               = 'gene_sequence'
settings['hisat_genome_directory']                = 'hisat_genome'
settings['hisat_genome_name']                     = 'genome'
settings['bowtie2_genome_name']                   = 'genome'
settings['five_prime_intron_directory']           = 'five_prime_intron_sequence'
settings['hisat_known_splice_sites_file']         = 'known_splicesites.txt'
# This is the base name used in .fa and .bed files
settings['five_prime_intron_name']                = 'five_prime_introns'

# These are optimized for intron references of length 20
# We require at least 18 matches (with no mismatches or indels)
# Or more than 18 matches with some mismatches or indels
#settings['bowtie2_local_parameters'] = ""


## Parameters used to align trimmed reads
settings['bowtie2_trimmed_read_parameters'] = " --norc --no-unal -L 10 -N 1 -a "
settings['bp_candidate_file_name']          = "bp_candidates.bed"
settings['bp_reference_base']               = "bp"
settings['genome_fasta_file']               = "genome.fa"

settings["align_bp_bt2_arguments"] = " -L 10 -N 1 "

settings["trim_length_threshold"] = 8

settings['number_of_nucleotides'] = 100

# bp_side_in_coverage
settings['bp_side_coverage_threshold'] = 8
