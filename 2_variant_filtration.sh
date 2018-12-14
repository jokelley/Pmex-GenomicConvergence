################################################################################################
#
#           2. Calling/Filtering SNPs                      
#
#
################################################################################################

# Use GATK to call SNPs on a per population basis
# Create sequence dictionary for reference prior to calling SNPs

java -jar picard.jar CreateSequenceDictionary \ 
      R=<path to reference> \
      O=<path to reference without .fa extension>.dict

# Use UnifiedGenotyper in GATK to call SNPs per individual
# Use EMIT_ALL_SITES to show all sites in the resulting vcf files (facilitates merging of vcf files later)
# Example: Call SNPs for the individual from Banos del azufre
	  
java -Xmx2048m -jar GenomeAnalysisTK.jar \
-T UnifiedGenotyper \
-R <path to reference> \
-o banos.vcf.gz \
--output_mode EMIT_ALL_SITES \
-I banos_combined.bam \

# Use vcf-merge (Perl module from vcftools) to merge the individual vcf files into one vcf
# Move into directory that contains vcftools Perl modules
# Use bgzip to zip the merged vcf file, use tabix to index the merged vcf file

./vcf-merge banos.vcf.gz bonita.vcf.gz gloria.vcf.gz rosita.vcf.gz vs.vcf.gz lluvia.vcf.gz vg.vcf.gz elazufre.vcf.gz esperanza.vcf.gz ixtapangajoya.vcf.gz | bgzip -c > merged.vcf.gz
tabix merged.vcf.gz

# Filter merged vcf file with GATK
# Use SelectVariants to pull out biallelic SNPs

java -jar GenomeAnalysisTK.jar \
    -T SelectVariants \
    -R <path to reference> \
    -V merged.vcf.gz \
    -selectType SNP \
    -restrictAllelesTo BIALLELIC \
    -o merged_just_biallelic_snps.vcf

# Use VariantFiltration to apply standard GATK hard filters

java -jar GenomeAnalysisTK.jar \
    -T VariantFiltration \
    -R <path to reference> \
    -V merged_just_biallelic_snps.vcf \
    --filterExpression "QD < 2.0 || FS > 60.0 || MQ < 40.0 || MQRankSum < -12.5 || ReadPosRankSum < -8.0" \
    --filterName "my_snp_filter" \
    -o filter_applied_biallelic_snps.vcf

# Use SelectVariants to exclude sites that didn't pass the hard filters

java -jar GenomeAnalysisTK.jar \
    -T SelectVariants \
    -R <path to reference> \
    -V filter_applied_biallelic_snps.vcf \
    --excludeFiltered \
    -o filter_excluded_biallelic_snps.vcf

# Use vcftools to further filter the SNPs. --minDP 10 only keeps genotypes supported by at least 8x coverage. --max-missing 1.0 only keeps sites where all individuals have an inferred genotype.
# The combination of min-alleles 2 and max-alleles 2 ensures that only biallelic sites are kept
# Vcftools adds a suffix to the end of the out file automatically (".recode.vcf")

vcftools --vcf filter_excluded_biallelic_snps.vcf --min-alleles 2 --max-alleles 2 --minQ 30 --minDP 10 --max-missing 1.0 --recode --recode-INFO-all --out merged_snps_minq30_minDP10_maxmiss1.0