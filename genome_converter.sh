#!/bin/sh
#
#BATCH --verbose
#SBATCH --job-name=genomeconverter
#SBATCH --output=/scratch/mrp420/reports/slurm_%j.out
#SBATCH --error=/scratch/mrp420/reports/slurm_%j.err
#SBATCH --time=4:00:00
#SBATCH --nodes=1
#SBATCH --mem=10GB
#SBATCH --mail-type=ALL
#SBATCH --mail-user=mrp420@nyu.edu

module load bowtie/gnu/1.2.0 
#Start time
echo 'Start Time'
echo $(date +%x_%r)

#Current genome
TAG1=$arg1
#Target genome
TAG2=$arg2
#bed file of loci
TAG3=$arg3
#Identifier for sequence output
TAG4=$arg4

cd /scratch/mrp420/
mkdir /scratch/mrp420/genome_convert
cd /scratch/mrp420/genome_convert

#Convert my input bed file to a fasta
module load bedtools/intel/2.26.0
bedtools getfasta -name -fi ${TAG1} -bed ${TAG3} -fo ${TAG4}.fa.out

#Make blast library for comparison
module load blast+/intel/2.6.0
makeblastdb -in ${TAG2} -dbtype nucl
mv ${TAG2}.n* .

#Blast back to new genome
blastn -query test.fa.out -db ${TAG2} -outfmt 6 | sort -k1,1 -k12,12nr -k11,11n | sort -u -k1,1 --merge > ${TAG4}.blastn

#Convert output back into bed file
awk '{OFS="\t"} {print $2,$9,$10,$1,$3,$11}' ${TAG4}.blastn > ${TAG4}.bed
