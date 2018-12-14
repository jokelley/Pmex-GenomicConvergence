################################################################################################
#
#           4. Local ancestry analysis using Saguaro                      
#
#
################################################################################################

# Convert vcf file to HMM format (this is the required format to use Saguaro)
# Use the vcf file that does not contain the outgroup. Saguaro is unsupervised, meaning that there is no way to indicate an outgroup
# -i indicates the input vcf file, -o indicates the output HMM file

VCF2HMMFeature -i merged.vcf.gz \
-o merged.saguaro

# Run Saguaro, -f indicates the input HMM file, -o indicates the output directory, -iter indicates the number of iterations to go through (29 results in 30 cacti)

Saguaro -f merged.saguaro -o <path to desired output directory> -iter 29