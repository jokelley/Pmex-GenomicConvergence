################################################################################################
#
#           1. Read trimming and mapping                      
#
#
################################################################################################

# First part of the script locates all fastq files from current directory to all subfolders and then runs FastQC on them. Note, the output directory is hardcoded so you might want to change it
####################################################################################################################################

# Move into directory that contains raw files
cd raw_files

# Find all files ending in .fastq.gz from current directory and all subfolders. Output full path (`pwd` needed otherwise you'd only get the path from the current directory). -type specifies a normal file. -name allows us to look for things based on file name.
# find `pwd` -type f -name *.fastq.gz
# Saves names of directories to an array
fileList=($(find `pwd` -type f -name '*.fastq'))

# Loops through all elements in the array. Note that "i" holds the value for the element in the array during each iteration. In this case "i" holds the path to the fastq file. Path to output directory is fixed. Runs all FastQC on all fastq files located in array.
for i in ${fileList[@]}
do

fastqc --outdir fastqc_dir --noextract --nogroup $i

done

# Runs Trim Galore on the files specified in the "TrimGalore_fastq_Locations.txt"

# The for loop is hard coded to the number to lines in the "TrimGalore_fastq_Locations.txt" file
# "cat" opens the file
# "sed" is grabbing the ith line of the file
# "awk is grabbing the 1st or 2nd column from the line above
# 45 total pairs of reads (10 samples sequenced across 4 or 5 lanes each)

for i in {1..45}
 
do

# TrimGalore_fastq_Locations_psulph.txt has two columns
# Column one contains the paths to raw read one fastq files
# Column two contains the paths to raw read two fastq files

read1=$(cat TrimGalore_fastq_Locations.txt | sed -n ''$i'p' | awk '{print $1}')
read2=$(cat TrimGalore_fastq_Locations.txt | sed -n ''$i'p' | awk '{print $2}')
 
trim_galore --quality 20 --fastqc_args "--noextract --nogroup" --stringency 6 --gzip --length 50 --clip_R1 5 --clip_R2 5 --paired ${read1} ${read2}
 
done

# Index reference genome
bwa index <path to reference genome>
samtools faidx <path to reference genome>

# Map to reference genome using BWA-mem
# Use same loop structure seen above

for i in {1..45}

do

# trimmed_fastq_Locations.txt has three columns
# Column one contains the paths to trimmed read one fastq files
# Column two contains the paths to trimmed read two fastq files
# Column three contains the name of each sample

read1=$(cat trimmed_fastq_Locations.txt | sed -n ''$i'p' | awk '{print $1}')
read2=$(cat trimmed_fastq_Locations.txt | sed -n ''$i'p' | awk '{print $2}')
filename=$(cat trimmed_fastq_Locations.txt | sed -n ''$i'p' | awk '{print $3}')

bwa mem <path to reference genome> ${read1} ${read2} > ${filename}.sam

done


# Process the sam files after mapping to the reference genome

for i in {1..45}

do

# sam_locations.txt has five columns
# Column one contains the paths to sam files
# Column two contains the sample names
# Column three contains the lane number (samples were sequenced across multiple sequencing lanes)
# Column four contains the samplename and the lane number (to have unique file names if they are in the same directory)

samfile=$(cat sam_locations.txt | sed -n ''$i'p' | awk '{print $1}')
samplename=$(cat sam_locations.txt | sed -n ''$i'p' | awk '{print $2}')
lane=$(cat sam_locations.txt | sed -n ''$i'p' | awk '{print $3}')
samplewithlane$(cat sam_locations.txt | sed -n ''$i'p' | awk '{print $4}')

# First, convert to bam file to save space
# -u = Output uncompressed BAM, -b = Output in the BAM format, -h = Include the header in the output, -t = links to .fai index
samtools view -ubht <path to reference>.fai ${samfile} -o ${samplewithlane}.bam

# SortSam: Sort the bam files by coordinate
java -Xmx4G -jar picard.jar SortSam INPUT=${samplewithlane}.clean.fix.bam OUTPUT=${samplewithlane}.F.bam SORT_ORDER=coordinate VALIDATION_STRINGENCY=LENIENT

# Index the bam files
samtools index ${samplewithlane}.F.bam

# Check alignment statistics using bamtools
bamtools stats -in ${samplewithlane}.F.bam > ${samplewithlane}.STATS

# AddOrReplaceReadGroups: Add readgroups to the bam files for SNP calling
java -Xmx6G -jar /home/repository_software/picard-tools-1.133/picard.jar AddOrReplaceReadGroups INPUT=${samplewithlane}.F.bam OUTPUT=${samplewithlane}.F.readgroup.bam SORT_ORDER=coordinate RGPL=Illumina RGLB=${samplewithlane} RGSM=${samplename} RGPU=${lane}

done

# Use samtools merge to combine bam files from the same sample into one bam file (this is an example of merging 4 bams into one)
samtools merge ${samplename}_combined.bam ${samplename}1.F.bam ${samplename}2.F.bam ${samplename}3.F.bam ${samplename}4.F.bam