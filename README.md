# Other_HTseq
A collection of scripts that don't belong anywhere specific

## genome_converter.sh ##
This script will take a bed file containing positions of regions of interest from a specific genome, and convert to a new bed file with coordinates updated for a different genome version or strain. This extracts sequence from orginal fasta using the positons, then blasts them against new genome to find new positions. 
