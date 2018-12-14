################################################################################################
#
#           3. Inferring phylogeny using SNPhylo                      
#
#
################################################################################################

# Use SNPhylo to infer maximum likelihood phylogeny using whole genome data
# First need to combine filtered vcf containing our samples with data from an outgroup
# Download outgroup data (guppy, Poecilia reticulata) using fastq-dump

fastq-dump --outdir sra_guppy_data --gzip --skip-technical --read-filter pass --dumpbase --split-3 --clip SRR1171024

# Go through the same steps as we did with the other samples (in scripts labelled 1 and 2)
# Merge vcf file from guppy with vcf file from all other samples (using vcf-merge from vcftools)

# Use vcf-merge (Perl module from vcftools) to merge the individual vcf files into one vcf
# Move into directory that contains vcftools Perl modules
# Use bgzip to zip the merged vcf file, use tabix to index the merged vcf file

./vcf-merge merged.vcf.gz guppy.vcf.gz | bgzip -c > merged_guppy.vcf.gz
tabix merged_guppy.vcf.gz

# Use bcftools to rename contig names to numbers because SNPhylo only recognizes numbers as valid contig names
# Chromosome_rename.txt contains two columns, one column contains the old names, the second contains the new names (numbered in order of appearance in the vcf)

bcftools annotate --rename-chrs chromosome_rename.txt merged_guppy.vcf.gz > \
merged_guppy_renamed_for_snphylo.vcf.gz

# Run Snphylo, -v is the vcf file, -d is the name of the gds file (snphylo converts into this format), -a should be the number of contigs present (23 is the default)
# -b indicates that the user wants to perform non-parametric bootstrapping, -B indicates the number of bootstraps

snphylo.sh -v merged_guppy_renamed_for_snphylo.vcf.gz -d snphylo.output.gds -o guppy -a <# of contigs> \
-b -B 1000