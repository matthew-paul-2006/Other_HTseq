#! /bin/bash
#SBATCH --verbose
#SBATCH --nodes=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=16GB
#SBATCH --time=4:00:00
#SBATCH --job-name=SNPfinder
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=mrp420@nyu.edu
#SBATCH --output=/scratch/%u/%x_%j.out
#SBATCH --error=/scratch/%u/%x_%j.err

#------------------------------------------------------------------------------#
#                                INSTRUCTIONS                                  #
#------------------------------------------------------------------------------#

# Adapted from Matts "SNPfinder.pbs" by Luis Silva.
# Use this script to map files then produce a sorted file with a separate index
# file for a region of interest.
# Can then use outputs in IGV to look for mismatches. For this later analysis
# you will need to load the file into IGV (ensure you also have the indexed file
# in the same directory). Once loaded in you need to zoom in on your gene of
# interest and change the settings to higlight mismatches (right-click). You
# will look for where you get sweeps (changes that occur in all of your reads).
# You can then turn on the translation setting to see what codons have changed
# and assess if it is synonymous. 

### Argument options:
# EXPID     Custom ID for output files
# RUNDIR    Path to directory to run script and save output in
# GENNAME   Basename of reference genome files preceded by absolute path to
#           directory containing the files. Must include a FASTA file and a
#           matching GFF file with the same basename.
#           An existing Bowtie index with a basename matching the files' is used
#           if found in the same directory; otherwise a new index is built
# FQ        Absolute path to input fastq file
# CHR       Name of chromosome in which the region of interest is located
# START     Start position of region of interest
# END       End position of region of interest

### EXAMPLE:
# sbatch --export \
# EXPID="ah9120a-0817",\
# RUNDIR="/scratch/lv38",\
# GENNAME="/home/lv38/Library/SK1Yue/Yue.SK1.genome.nuclear.mito.2micr",\
# FQ="/scratch/lv38/HN3C3AFXX_n01_ah9120a-0817.fastq.gz",\
# CHR="chrIX",START="211476",END="213293" \
# SNPfinder.slurm

#------------------------------------------------------------------------------#
#                                  Functions                                   #
#------------------------------------------------------------------------------#

function elapsed_time() {
    ENDTIME=$(date +%s)

    TIME=$(($ENDTIME - $1))
    if [ $TIME -lt 60 ]
    then
        echo "$TIME sec"
    elif [ $TIME -ge 60 ]  && [ $TIME -lt 3600 ]
    then
        echo "$(($TIME / 60)) min"
    else
        echo "$(($TIME / 60 / 60)) hr"
            fi
}

function check_arg() {
    if [ -z "$1" ]
    then
        echo ">>>>> Please provide values for all required arguments"
        exit 1
    fi
}

#------------------------------------------------------------------------------#
#                                  IO checks                                   #
#------------------------------------------------------------------------------#

# Check arguments
check_arg $EXPID
check_arg $RUNDIR 
check_arg $FQ
check_arg $CHR
check_arg $START
check_arg $END

# Check input files / dirs
[ -f $FQ ] || { echo "Could not find file: $FQ"; exit 1; }
[ -d $RUNDIR ] || { echo "Could not find directory: $RUNDIR"; exit 1; }
GENDIR=$(dirname "$GENNAME")
GENNAME=$(basename "$GENNAME")
[ -d $GENDIR ] || { echo "Could not find directory: $GENDIR"; exit 1; }

# Search for reference genome fasta file; exit if not found
FA=$(find $GENDIR -iname "${GENNAME}.fa*")

if [ -z "$FA" ]
then
    echo ">>>>> Could not find reference genome FASTA file"
    exit 1
fi

# Search for Bowtie index and build one if not found
# (a file named according to rule "fasta_base_name.1.ebwt")
# The following code will return the full basename
# (the provided $GENNAME might not include it in full)
IX=$(basename $FA)                              # Drop path to file
IX=${IX%.*}                                     # Drop extension
CHECKIX=$(find $GENDIR -iname "${IX}.1.ebwt")   # Search file
#IX=$(basename $IX | cut -d '.' -f 1)

if [ -z "$CHECKIX" ]
then
    echo ">>>>> Building Bowtie index..."
    module purge
    module load bowtie/gnu/1.2.1.1
    # Build index
    cd $GENDIR
    bowtie-build -f $FA $IX
fi

#------------------------------------------------------------------------------#
#                                                                              #
#                                Run pipeline                                  #
#                                                                              #
#------------------------------------------------------------------------------#

STARTTIME=$(date +%s)
echo \
"------------------------------------------------------------------------------"
echo ">>>>> Started SNPfinder: $EXPID"
echo \
"------------------------------------------------------------------------------"
date

#------------------------------------------------------------------------------#
#                  Align reads to reference genome with Bowtie                 #
#------------------------------------------------------------------------------#
cd $RUNDIR

# Abort if output directory already exists
if [ -d "$EXPID" ]
then
    echo ">>>>> Output directory already exists"
    exit 1
fi

mkdir ${EXPID}/
cd ${EXPID}/

echo ">>>>> Map reads with Bowtie:"
echo "      Strict mapping; discard perfectly matching reads..."
module purge
module load bowtie/gnu/1.2.1.1

bowtie -q -5 1 -v 0 -m 1 -p 8 \
    --un ${EXPID}_unmapped.fastq \
    --seed=123 \
    -S $GENDIR/$IX $FQ strictlymapped.sam

echo "      Loosely map reads containing mismatches (previously unaligned)..."

bowtie -q -5 1 -v 3 -p 8 \
    --seed=123 \
    -S $GENDIR/$IX ${EXPID}_unmapped.fastq \
    ${EXPID}_looselymapped.sam


#------------------------------------------------------------------------------#
#                           Get region of interest                             #
#------------------------------------------------------------------------------#
echo ">>>>> Convert SAM files to BAM format, index and sort..."
module purge
module load samtools/intel/1.3.1

samtools view -Sb ${EXPID}_looselymapped.sam > ${EXPID}_looselymapped.bam
samtools sort -o ${EXPID}_sorted.bam -O 'bam' \
    -T 'temp' ${EXPID}_looselymapped.bam
samtools index ${EXPID}_sorted.bam

echo ">>>>> Extract region of interest:"
echo "      $CHR:$START-$END"
samtools view ${EXPID}_sorted.bam -h $CHR:$START-$END \
    > ${EXPID}_region.bam
samtools sort -o ${EXPID}_final.bam -O 'bam' -T 'temp' ${EXPID}_region.bam
samtools index ${EXPID}_final.bam

#------------------------------------------------------------------------------#
#                                  Clean up                                    #
#------------------------------------------------------------------------------#
echo ">>>>> Delete unnecessary files..."
#rm ${EXPID}_sorted.bam
#rm ${EXPID}_sorted.bam.bai
#rm ${EXPID}_looselymapped.bam
#rm ${EXPID}_region.bam
#rm ${EXPID}_unmapped.fastq
#rm strictlymapped.sam
#rm ${EXPID}_looselymapped.sam

#------------------------------------------------------------------------------#
ELAPSEDTIME=$(elapsed_time $STARTTIME)
echo
echo "-----"
echo "-----"
echo "Completed pipeline in $ELAPSEDTIME"
echo \
"------------------------------------------------------------------------------"

