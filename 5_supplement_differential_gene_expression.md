# Analysis pipeline used to generate Figure S3

## Trimming with TrimGalore

##### Script is /data/kelley/projects/kerry/mx_rna_seq/1_trimgalore/scripts/trim.sh
    #!/bin/bash
    #SBATCH --job-name=trim
    #SBATCH --partition=kamiak
    #SBATCH --output=/data/kelley/projects/kerry/mx_rna_seq/1_trimgalore/scripts/trim.out
    #SBATCH --error=/data/kelley/projects/kerry/mx_rna_seq/1_trimgalore/scripts/trim.err
    #SBATCH --workdir=/data/kelley/projects/kerry/mx_rna_seq/1_trimgalore
    #SBATCH --mail-user=kerry.mcgowan@wsu.edu
    #SBATCH --mail-type=BEGIN,END,FAIL
    
    # PURPOSE: Trim reads using TrimGalore default settings.
    
    while read line || [ -n "$line" ];
    do
    	read1=$(echo $line | awk '{print $1}')
    	read2=$(echo $line | awk '{print $2}')
    
    trim_galore \
    	--fastqc_args "--noextract --nogroup" \
    	--output_dir /data/kelley/projects/kerry/mx_rna_seq/1_trimgalore \
    	--paired /data/kelley/projects/pmex/field2010/$read1 /data/kelley/projects/pmex/field2010/$read2
    
    done < /data/kelley/projects/kerry/mx_rna_seq/1_trimgalore/scripts/MXsample_list.txt
    
    # Top rows of MXsample_list.txt
    
    MX02_1.q33.fastq.gz	MX02_2.q33.fastq.gz
    MX03_1.q33.fastq.gz	MX03_2.q33.fastq.gz
    MX04_1.q33.fastq.gz	MX04_2.q33.fastq.gz
    MX05_1.q33.fastq.gz	MX05_2.q33.fastq.gz
 
 ## Pre-processing with Cufflinks gffread
 
     # Need to change the GFF file into a GTF file otherwise it won't work with the HISAT2 python script in the next step
     # *Be careful not to overwrite the orignal GFF file (I made a copy)
     # Did this in an idev
     
     gffread GCF_001443325.1_P_mexicana-1.0_genomic_copy.gff -T -o GCF_001443325.1_P_mexicana-1.0_genomic_copy.gtf
     
     mv GCF_001443325.1_P_mexicana-1.0_genomic_copy.gtf GCF_001443325.1_P_mexicana-1.0_genomic.gtf
     
     chmod g+x GCF_001443325.1_P_mexicana-1.0_genomic.gtf

## Indexing the reference genome with HISAT2

    # Need to make fiiles using python scripts provided by HISAT2 for the --ss and --exon options

##### Script is /data/kelley/projects/kerry/mx_rna_seq/2_hisat2/scripts/hisat2_extract_splice_sites.sh 

    #!/bin/bash
    #SBATCH --job-name=extract
    #SBATCH --partition=kamiak
    #SBATCH --output=/data/kelley/projects/kerry/mx_rna_seq/2_hisat2/scripts/extract.out
    #SBATCH --error=/data/kelley/projects/kerry/mx_rna_seq/2_hisat2/scripts/extract.err
    #SBATCH --workdir=/data/kelley/projects/kerry/mx_rna_seq/2_hisat2
    #SBATCH --mail-user=kerry.mcgowan@wsu.edu
    #SBATCH --mail-type=BEGIN,END,FAIL
    
    # PURPOSE: Extract splice sites from pmex GFF file prior to building the index.
    
    module load python/2.7.10
    
    /data/kelley/projects/programs/hisat2-2.1.0/hisat2_extract_splice_sites.py \
    /data/kelley/projects/anthonys_projects/pmex_genome_resequencing/reference_genome_pmex/GCF_001443325.1_P_mexicana-1.0_genomic.gtf \
    > /data/kelley/projects/kerry/mx_rna_seq/2_hisat2/splice_sites.txt

##### Script is /data/kelley/projects/kerry/mx_rna_seq/2_hisat2/scripts/hisat2_extract_exons.sh 

    #!/bin/bash
    #SBATCH --job-name=extract2
    #SBATCH --partition=kamiak
    #SBATCH --output=/data/kelley/projects/kerry/mx_rna_seq/2_hisat2/scripts/extract2.out
    #SBATCH --error=/data/kelley/projects/kerry/mx_rna_seq/2_hisat2/scripts/extract2.err
    #SBATCH --workdir=/data/kelley/projects/kerry/mx_rna_seq/2_hisat2
    #SBATCH --mail-user=kerry.mcgowan@wsu.edu
    #SBATCH --mail-type=BEGIN,END,FAIL
    
    # PURPOSE: Extract exons from pmex GFF file prior to building the index.
    
    module load python/2.7.10
    
    /data/kelley/projects/programs/hisat2-2.1.0/hisat2_extract_exons.py \
    /data/kelley/projects/anthonys_projects/pmex_genome_resequencing/reference_genome_pmex/GCF_001443325.1_P_mexicana-1.0_genomic.gtf \
    > /data/kelley/projects/kerry/mx_rna_seq/2_hisat2/exons.txt

##### Script is /data/kelley/projects/kerry/mx_rna_seq/2_hisat2/scripts/index_the_reference_genome.sh

    #!/bin/bash
    #SBATCH --job-name=index
    #SBATCH --partition=kamiak
    #SBATCH --nodes=1
    #SBATCH --mem=150G    ### Amount of memory
    #SBATCH --output=/data/kelley/projects/kerry/mx_rna_seq/2_hisat2/scripts/index.out
    #SBATCH --error=/data/kelley/projects/kerry/mx_rna_seq/2_hisat2/scripts/index.err
    #SBATCH --workdir=/data/kelley/projects/kerry/mx_rna_seq/2_hisat2
    #SBATCH --mail-user=kerry.mcgowan@wsu.edu
    #SBATCH --mail-type=BEGIN,END,FAIL
    
    # PURPOSE: Builds an index for the reference genome.
    
    hisat2-build \
    	-f /data/kelley/projects/anthonys_projects/pmex_genome_resequencing/reference_genome_pmex/GCF_001443325.1_P_mexicana-1.0_genomic.fa \
    	--ss /data/kelley/projects/kerry/mx_rna_seq/2_hisat2/splice_sites.txt \
    	--exon /data/kelley/projects/kerry/mx_rna_seq/2_hisat2/exons.txt \
    	/data/kelley/projects/anthonys_projects/pmex_genome_resequencing/reference_genome_pmex/GCF_001443325.1_P_mexicana-1.0_genomic

## Mapping using HISAT2

##### Script is /data/kelley/projects/kerry/mx_rna_seq/2_hisat2/scripts/mapping.sh

    #!/bin/bash
    #SBATCH --job-name=mapping
    #SBATCH --partition=cas
    #SBATCH --output=/data/kelley/projects/kerry/mx_rna_seq/2_hisat2/scripts/mapping.out
    #SBATCH --error=/data/kelley/projects/kerry/mx_rna_seq/2_hisat2/scripts/mapping.err
    #SBATCH --workdir=/data/kelley/projects/kerry/mx_rna_seq/2_hisat2
    #SBATCH --mail-user=kerry.mcgowan@wsu.edu
    #SBATCH --mail-type=BEGIN,END,FAIL
    
    # PURPOSE: Map reads to the previously indexed reference genome.
    
    mkdir -p /data/kelley/projects/kerry/mx_rna_seq/2_hisat2/stats
    
    mkdir -p /data/kelley/projects/kerry/mx_rna_seq/2_hisat2/unmapped
    
    while read line || [ -n "$line" ];
    do
    	read1=$(echo $line | awk '{print $1}')
    	read2=$(echo $line | awk '{print $2}')
    	filename=$(echo $line | awk '{print $3}')
    
    hisat2 \
    	--phred33 \
    	-k 10 \
    	--met-file /data/kelley/projects/kerry/mx_rna_seq/2_hisat2/stats/$filename.stats \
    	--rg-id $filename \
    	--rg SM:$filename \
    	--rg PL:illumina \
    	-p 1 \
    	--rna-strandness RF \
    	--fr \
    	--dta \
    	--un-conc-gz /data/kelley/projects/kerry/mx_rna_seq/2_hisat2/unmapped/$filename.unmapped \
    	-x /data/kelley/projects/anthonys_projects/pmex_genome_resequencing/reference_genome_pmex/GCF_001443325.1_P_mexicana-1.0_genomic \
    	-1 /data/kelley/projects/kerry/mx_rna_seq/1_trimgalore/$read1 \
    	-2 /data/kelley/projects/kerry/mx_rna_seq/1_trimgalore/$read2 \
    	-S /data/kelley/projects/kerry/mx_rna_seq/2_hisat2/$filename.sam
    
    done < /data/kelley/projects/kerry/mx_rna_seq/2_hisat2/scripts/MXsamples.txt
    
    # Top rows of MXsample_list.txt
    MX02_1.q33_val_1.fq.gz	MX02_2.q33_val_2.fq.gz	MX02.q33
    MX03_1.q33_val_1.fq.gz	MX03_2.q33_val_2.fq.gz	MX03.q33
    MX04_1.q33_val_1.fq.gz	MX04_2.q33_val_2.fq.gz	MX04.q33

## Converting to BAM and sorting using SAMtools

##### Script is /data/kelley/projects/kerry/mx_rna_seq/3_samtools/scripts/samtools_view_sort_loop.sh

    #!/bin/bash
    #SBATCH --partition=popgenom
    #SBATCH --job-name=samtools
    #SBATCH --output=/data/kelley/projects/kerry/mx_rna_seq/3_samtools/scripts/samtools.out
    #SBATCH --error=/data/kelley/projects/kerry/mx_rna_seq/3_samtools/scripts/samtools.err
    #SBATCH --workdir=/data/kelley/projects/kerry/mx_rna_seq/3_samtools
    #SBATCH --mail-type=END,FAIL
    #SBATCH --mail-user=kerry.mcgowan@wsu.edu
    
    while read line || [ -n "$line" ];
    do
    	sampleName=$(echo $line | awk '{print $1}')
    
    samtools view -b -h -S /data/kelley/projects/kerry/mx_rna_seq/2_hisat2/$sampleName.sam > /data/kelley/projects/kerry/mx_rna_seq/3_samtools/$sampleName.bam
    
    samtools sort /data/kelley/projects/kerry/mx_rna_seq/3_samtools/$sampleName.bam -o /data/kelley/projects/kerry/mx_rna_seq/3_samtools/$sampleName.sorted.bam
    
    done < /data/kelley/projects/kerry/mx_rna_seq/3_samtools/scripts/sampleNames.txt
    
    # Top rows of sampleNames.txt
    MX02.q33
    MX03.q33
    MX04.q33

## Generating GTF files using StringTie

##### Script is /data/kelley/projects/kerry/mx_rna_seq/4_stringtie/scripts/stringtie.sh

    #!/bin/bash
    #SBATCH --partition=popgenom
    #SBATCH --job-name=stringtie
    #SBATCH --cpus-per-task=4
    #SBATCH --output=/data/kelley/projects/kerry/mx_rna_seq/4_stringtie/scripts/stringtie.out
    #SBATCH --error=/data/kelley/projects/kerry/mx_rna_seq/4_stringtie/scripts/stringtie.err
    #SBATCH --workdir=/data/kelley/projects/kerry/mx_rna_seq/4_stringtie
    #SBATCH --mail-type=END,FAIL
    #SBATCH --mail-user=kerry.mcgowan@wsu.edu
    
    # PURPOSE: Generate GTF files in a ballgown folder.
    
    while read line || [ -n "$line" ];
    do
    	sampleName=$(echo $line | awk '{print $1}')
    
    mkdir -p /data/kelley/projects/kerry/mx_rna_seq/4_stringtie/ballgown
    
    mkdir -p /data/kelley/projects/kerry/mx_rna_seq/4_stringtie/ballgown/$sampleName
    
    stringtie \
    	/data/kelley/projects/kerry/mx_rna_seq/3_samtools/$sampleName.sorted.bam \
    	-o /data/kelley/projects/kerry/mx_rna_seq/4_stringtie/ballgown/$sampleName/$sampleName.gtf \
    	-p 4 \
    	-G /data/kelley/projects/anthonys_projects/pmex_genome_resequencing/reference_genome_pmex/GCF_001443325.1_P_mexicana-1.0_genomic.gff \
    	-B \
    	-e
    
    done < /data/kelley/projects/kerry/mx_rna_seq/4_stringtie/scripts/sampleNames.txt

    # Top rows of sampleNames.txt
    MX02.q33
    MX03.q33
    MX04.q33

## Generating a gene counts matrix (gene_count_matrix.csv) using prepDE.py (from StringTie)

##### Re-downlowaded prepDE.py script from the link in the StringTie manual and called it prepDE_20190115.py

##### NOTE: Sample 78 was removed before stringtie.sh from /data/kelley/projects/kerry/mx_rna_seq/4_stringtie/ballgown into ballgown_removed_samples because it was from a cave population

##### Script is /data/kelley/projects/kerry/mx_rna_seq/5_prepDE/scripts/prepDE.sh

    #!/bin/bash
    #SBATCH --partition=popgenom
    #SBATCH --job-name=prepDE
    #SBATCH --output=/data/kelley/projects/kerry/mx_rna_seq/5_prepDE/scripts/prepDE.out
    #SBATCH --error=/data/kelley/projects/kerry/mx_rna_seq/5_prepDE/scripts/prepDE.err
    #SBATCH --workdir=/data/kelley/projects/kerry/mx_rna_seq/5_prepDE
    #SBATCH --mail-type=END,FAIL
    #SBATCH --mail-user=kerry.mcgowan@wsu.edu
    
    #Purpose: Runs python script provided by StringTie to generate 2 CSV files containing the count matrices for genes and transcripts
    
    # Note what the working directory is set to
    
    # -i <path to stringtie output>
    # -g <where to output gene matrix>
    # -t <where to output transcript matrix>
    
    
    
    module load python/2.7.10
    
    python /data/kelley/projects/kerry/mx_rna_seq/5_prepDE/scripts/prepDE_20190115.py \
    	-i /data/kelley/projects/kerry/mx_rna_seq/4_stringtie/ballgown \
    	-g /data/kelley/projects/kerry/mx_rna_seq/5_prepDE/gene_count_matrix.csv \
    	-t /data/kelley/projects/kerry/mx_rna_seq/5_prepDE/transcript_count_matrix.csv
        
## Differential gene expression analyses using EdgeR

    #############################
    #### Pmex EdgeR Analysis #### 
    #############################
    
    ####################################
    #### Read in CSV of Gene Counts ####
    ####################################
    
    setwd('~/Documents/WSU/RESEARCH/convergent_evo_MX/5_prepDE/')
    data <- read.csv("gene_count_matrix.csv",row.names=1) # this is w/o cave sample #78
    
    #############################
    #### Change Column Names ####
    #############################
    
    names(data)[1] <- "MX02"
    names(data)[2] <- "MX03"
    names(data)[3] <- "MX04"
    names(data)[4] <- "MX05"
    names(data)[5] <- "MX06"
    names(data)[6] <- "MX07"
    names(data)[7] <- "MX29"
    names(data)[8] <- "MX30"
    names(data)[9] <- "MX31"
    names(data)[10] <- "MX32"
    names(data)[11] <- "MX34"
    names(data)[12] <- "MX35"
    names(data)[13] <- "MX44"
    names(data)[14] <- "MX45"
    names(data)[15] <- "MX46"
    names(data)[16] <- "MX48"
    names(data)[17] <- "MX49"
    names(data)[18] <- "MX50"
    names(data)[19] <- "MX51"
    names(data)[20] <- "MX52"
    names(data)[21] <- "MX53"
    names(data)[22] <- "MX54"
    names(data)[23] <- "MX55"
    names(data)[24] <- "MX57"
    names(data)[25] <- "MX59"
    names(data)[26] <- "MX60"
    names(data)[27] <- "MX61"
    names(data)[28] <- "MX62"
    names(data)[29] <- "MX63"
    names(data)[30] <- "MX71"
    names(data)[31] <- "MX73"
    names(data)[32] <- "MX74"
    names(data)[33] <- "MX75"
    names(data)[34] <- "MX76"
    names(data)[35] <- "MX77"
    
    ###############################
    #### Set Working Directory ####
    ###############################
    
    setwd('~/Documents/WSU/RESEARCH/convergent_evo_MX/6_edgeR/')
    
    ##############################
    #### Package Installation ####
    ##############################
    
    # install EdgeR and limma, only have to do this the first time unless you update the version of R
    
    # source("https://bioconductor.org/biocLite.R")
    # biocLite("edgeR")
    # biocLite("limma")
    
    #######################
    #### Load Packages ####
    #######################
    
    # load EdgeR and limma
    
    require("limma")
    require("edgeR")
    
    ###################################
    ####  Create grouping factor  #####
    ###################################
    
    # to be incorporated later into the DGEList
    
    group = c(rep("Pich",12),rep("Puya",11),rep("Taco",12))
    
    ##########################################################
    ####  Remove rows where counts = 0 across all samples ####
    ##########################################################
    
    dim(data)
      # 31791 genes before removing the 0 counts
    
    total <- data[rowSums(data) > 0,] # trim out any rows with no reads counts
    
    dim(total)
      # 26817 genes after removing the 0 counts
    
    ###################################
    ####  Creating DGEList object  ####
    ###################################
    
    y <- DGEList(counts=total, group=group)
    
    #########################
    ####  Normalization  ####
    #########################
    
    # Normalization for RNA composition (i.e. prevents genes from appearing falsely downregulated)
    
    # Creates effective library sizes
    
    y <- calcNormFactors(y) # calculates normalization factors to scale raw library sizes, TMM is default method
    
    y$samples # shows library sizes and normalization factors
    
    #######################################################
    ######## Multidimensional Scaling (MDS) Plot  #########
    #######################################################
    
    # gene.selection = 'common' selects the same genes for all comparisons
    
    # method='logFC' required with gene.selection argument
    
    # top = 500 plots the top 500 genes (this is the default)
    
    pdf("MDSplot_500genes_colored_by_drainage.pdf")
    par(cex.axis = 1.1, cex = 1.25)
    mds500 <- plotMDS(y, gene.selection = 'common', top = 500, method='logFC', cex=1, 
                  col=c(rep("blue",12), rep("red",11), rep("black",12)))
                  #xlim=c(-2.5,2.5), xlab = "MDS axis 1",
                  #ylim=c(-2.5,2.5), ylab = "MDS axis 2")
    legend("topright", cex = 1, 
           legend=c("Pich", "Puya", "Taco"), 
           col=c("blue","red","black"), 
           pch = c(15),
           pt.cex=1.75)
    print(mds500)
    dev.off()
      
    pdf("MDSplot_10000genes_colored_by_drainage.pdf")
    par(cex.axis = 1.1, cex = 1.25)
    mds10000 <- plotMDS(y, gene.selection = 'common', top = 10000, method='logFC', cex=1, 
                   col=c(rep("blue",12), rep("red",11), rep("black",12)))
                   #pch = c(rep(15,12), rep(15,11), rep(15,13)),
                   #xlim=c(-2.5,2.5), xlab = "MDS axis 1",
                   #ylim=c(-2.5,2.5), ylab = "MDS axis 2")
    legend("topright", cex = 1, 
           legend=c("Pich", "Puya", "Taco"), 
           col=c("blue","red","black"), 
           pch = c(15),
           pt.cex=1.75)
    print(mds10000)
    dev.off()
    
    ######################################
    #### Now separated by S versus NS ####
    ######################################
    
    pdf("MDSplot_500genes_colored_by_habitat.pdf")
    par(cex.axis = 1.1, cex = 1.25)
    mds500 <- plotMDS(y, gene.selection = 'common', top = 500, method='logFC', cex=1, 
                      col=c(rep("blue",6), rep("yellow",11), rep("blue",12), rep("yellow",6)))
    legend("topright", cex = 1, 
           legend=c("fresh", "sulfur"), 
           col=c("blue","yellow"), 
           pch = c(15),
           pt.cex=1.75)
    print(mds500)
    dev.off()
    
    pdf("MDSplot_10000genes_colored_by_habitat.pdf")
    par(cex.axis = 1.1, cex = 1.25)
    mds10000 <- plotMDS(y, gene.selection = 'common', top = 10000, method='logFC', cex=1, 
                      col=c(rep("blue",6), rep("yellow",11), rep("blue",12), rep("yellow",6)))
    legend("topright", cex = 1, 
           legend=c("fresh", "sulfur"), 
           col=c("blue","yellow"), 
           pch = c(15),
           pt.cex=1.75)
    print(mds10000)
    dev.off()
    
    #########################################################################################################################################
    #### Now split matrices by drainage #####################################################################################################
    #########################################################################################################################################
    
    # Pichucalco drainage
    pichu <- data[,1:12]
    
    # Puyacatengo drainage
    puya <- data[,c(23,22,21,20,19,18,17,16,15,14,13)]
    
    # Tacotalpa drainage
    taco <- data[,24:35]
    
    ########################################################################################################################
    #### PICHUCALCO ########################################################################################################
    ########################################################################################################################
    
    ###################################
    ####  Create grouping factor  #####
    ###################################
    
    # to be incorporated later into the DGEList
    
    group_pichu = c(rep("fresh",6),rep("sulfur",6))
    
    ##########################################################
    ####  Remove rows where counts = 0 across all samples ####
    ##########################################################
    
    dim(pichu)
    # 31791 genes before removing the 0 counts
    
    total_pichu <- pichu[rowSums(pichu) > 0,] # trim out any rows with no reads counts
    
    dim(total_pichu)
    # 24430 genes after removing the 0 counts
    
    ###################################
    ####  Creating DGEList object  ####
    ###################################
    
    y_pichu <- DGEList(counts=total_pichu, group=group_pichu)
    
    #########################
    ####  Normalization  ####
    #########################
    
    # Normalization for RNA composition (i.e. prevents genes from appearing falsely downregulated)
    
    # Creates effective library sizes
    
    y_pichu <- calcNormFactors(y_pichu) # calculates normalization factors to scale raw library sizes, TMM is default method
    
    y_pichu$samples # shows library sizes and normalization factors
    
    #######################################################
    ######## Multidimensional Scaling (MDS) Plot  #########
    #######################################################
    
    # gene.selection = 'common' selects the same genes for all comparisons
    
    # method='logFC' required with gene.selection argument
    
    # top = 500 plots the top 500 genes (this is the default)
    
    pdf("MDSplot_500genes_PICHU_colored_by_habitat.pdf")
    par(cex.axis = 1.1, cex = 1.25)
    mds500_pichu <- plotMDS(y_pichu, gene.selection = 'common', top = 500, method='logFC', cex=1, 
                  col=c(rep("blue",6), rep("yellow",6)))
                  #xlim=c(-2.5,2.5), xlab = "MDS axis 1",
                  #ylim=c(-2.5,2.5), ylab = "MDS axis 2")
    legend("topright", cex = 1, 
           legend=c("fresh", "sulfur"), 
           col=c("blue","yellow"), 
           pch = c(15),
           pt.cex=1.75)
    print(mds500_pichu)
    dev.off()
    
    pdf("MDSplot_10000genes_PICHU_colored_by_habitat.pdf")
    par(cex.axis = 1.1, cex = 1.25)
    mds10000_pichu <- plotMDS(y_pichu, gene.selection = 'common', top = 10000, method='logFC', cex=1, 
                            col=c(rep("blue",6), rep("yellow",6)))
    #xlim=c(-2.5,2.5), xlab = "MDS axis 1",
    #ylim=c(-2.5,2.5), ylab = "MDS axis 2")
    legend("topright", cex = 1, 
           legend=c("fresh", "sulfur"), 
           col=c("blue","yellow"), 
           pch = c(15),
           pt.cex=1.75)
    print(mds10000_pichu)
    dev.off()
    
    #######################################
    ####### Graphing Library Size #########
    #######################################
    
    # create color vector for graph
    colors_pichu <- c(rep("blue",6),rep("yellow",6))
    
    # las=2 makes label text perpendicular to axis
    par(las=2)
    
    # generate plot
    pdf("library_size_PICHU.pdf")
    library_pichu <- barplot(y_pichu$samples$lib.size, col=colors_pichu, 
                       cex.axis=.8, ylab = "Library Size", xlab = "Sample")
    legend("topright", cex = 0.82, 
           legend=c("fresh", "sulfur"),
           bty = "n", pt.cex=1.5, 
           fill=c("blue", "yellow"))
    print(library_pichu)
    dev.off()
    
    ####################################
    ########## Design Matrix ###########
    ####################################
    
    fac_pichu <- c(rep("fresh",6), rep("sulfur",6)) # use this to get groups
    
    fac_pichu <- factor(fac_pichu) # renames fac names with numbers
    
    design_pichu <- model.matrix(~0+fac_pichu) # 0 is the intercept
    
    colnames(design_pichu) <- levels(fac_pichu)
    
    design_pichu
    
    ###################################################
    ####### Estimating Dispersions - CR method ########
    ###################################################
    
    # Uses GLM instead of qCML method because testing multiple factors here
    
    # Uses Cox-Reid profile-adjusted likelihood
    
    # estimates common dispersion and tagwise dispersions in one run
    y_pichu <- estimateDisp(y_pichu, design_pichu)
    
    # common dispersions
    y_pichu$common.dispersion
    
    # tagwise (gene-specific dispersions)
    summary(y_pichu$tagwise.dispersion)
    
    # estimated prior degrees of freedom
    y_pichu$prior.df
    
    ###########################
    ####### BCV Plot ##########
    ###########################
    
    # BCV is the coefficient of variation with which the (unknown) true abundance of the gene varies between replicate RNA samples
    
    pdf("BCV_PICHU.pdf")
    BCV_pichu <- plotBCV(y_pichu)
    print(BCV_pichu)
    dev.off()
    
    ###############################################
    ####### DGE - Quasi-Likelihood F-test #########
    ###############################################
    
    # quasi-likelihood F-test better for bulk RNA-seq data because stricter error rate control, accounts for
      # uncertainty in dispersion estimation
    
    group_pichu <- factor(c(1,1,1,1,1,1,2,2,2,2,2,2))
    
    design_pichu <- model.matrix(~0+group_pichu)
    
    # QL model representing the study design fitted to the data
    fit_pichu <- glmQLFit(y_pichu, design_pichu)
    
    # tests use FDR < 0.05
    
    # compare FRESH (1) vs SULFUR (2)
    qlf_pichu <- glmQLFTest(fit_pichu, contrast=c(-1,1))
    fresh.vs.sulf_pichu <- topTags(qlf_pichu, n = 5000, p.value=0.05)
    fresh.vs.sulf_pichu <- data.frame(fresh.vs.sulf_pichu)
    require(xlsx)
    write.csv(x=fresh.vs.sulf_pichu, file = "fresh.vs.sulf_PICHU.csv")
    
    #####################################################################################################################################
    #### PUYACATENGO ####################################################################################################################
    #####################################################################################################################################
    
    ###################################
    ####  Create grouping factor  #####
    ###################################
    
    # to be incorporated later into the DGEList
    
    group_puya = c(rep("fresh",6),rep("sulfur",5))
    
    ##########################################################
    ####  Remove rows where counts = 0 across all samples ####
    ##########################################################
    
    dim(puya)
    # 31791 genes before removing the 0 counts
    
    total_puya <- puya[rowSums(puya) > 0,] # trim out any rows with no reads counts
    
    dim(total_puya)
    # 25223 genes after removing the 0 counts
    
    ###################################
    ####  Creating DGEList object  ####
    ###################################
    
    y_puya <- DGEList(counts=total_puya, group=group_puya)
    
    #########################
    ####  Normalization  ####
    #########################
    
    # Normalization for RNA composition (i.e. prevents genes from appearing falsely downregulated)
    
    # Creates effective library sizes
    
    y_puya <- calcNormFactors(y_puya) # calculates normalization factors to scale raw library sizes, TMM is default method
    
    y_puya$samples # shows library sizes and normalization factors
    
    #######################################################
    ######## Multidimensional Scaling (MDS) Plot  #########
    #######################################################
    
    # gene.selection = 'common' selects the same genes for all comparisons
    
    # method='logFC' required with gene.selection argument
    
    # top = 500 plots the top 500 genes (this is the default)
    
    pdf("MDSplot_500genes_PUYA_colored_by_habitat.pdf")
    par(cex.axis = 1.1, cex = 1.25)
    mds500_puya <- plotMDS(y_puya, gene.selection = 'common', top = 500, method='logFC', cex=1, 
                            col=c(rep("blue",6), rep("yellow",5)))
    #xlim=c(-2.5,2.5), xlab = "MDS axis 1",
    #ylim=c(-2.5,2.5), ylab = "MDS axis 2")
    legend("topright", cex = 1, 
           legend=c("fresh","sulfur"), 
           col=c("blue","yellow"), 
           pch = c(15),
           pt.cex=1.75)
    print(mds500_puya)
    dev.off()
    
    pdf("MDSplot_10000genes_PUYA_colored_by_habitat.pdf")
    par(cex.axis = 1.1, cex = 1.25)
    mds10000_puya <- plotMDS(y_puya, gene.selection = 'common', top = 10000, method='logFC', cex=1, 
                            col=c(rep("blue",6), rep("yellow",5)))
    #xlim=c(-2.5,2.5), xlab = "MDS axis 1",
    #ylim=c(-2.5,2.5), ylab = "MDS axis 2")
    legend("topright", cex = 1, 
           legend=c("fresh","sulfur"), 
           col=c("blue","yellow"), 
           pch = c(15),
           pt.cex=1.75)
    print(mds10000_pichu)
    dev.off()
    
    ################################################################################################################################################
    #### Take out Puyacatengo sample 51 because looks like outlier #################################################################################
    ################################################################################################################################################
    
    puya_minus51 <- puya[,-5]
    
    ###################################
    ####  Create grouping factor  #####
    ###################################
    
    # to be incorporated later into the DGEList
    
    group_puya_minus51 = c(rep("fresh",5),rep("sulfur",5))
    
    ##########################################################
    ####  Remove rows where counts = 0 across all samples ####
    ##########################################################
    
    dim(puya_minus51)
    # 31791 genes before removing the 0 counts
    
    total_puya_minus51 <- puya_minus51[rowSums(puya_minus51) > 0,] # trim out any rows with no reads counts
    
    dim(total_puya_minus51)
    # 24905 genes after removing the 0 counts
    
    ###################################
    ####  Creating DGEList object  ####
    ###################################
    
    y_puya_minus51 <- DGEList(counts=total_puya_minus51, group=group_puya_minus51)
    
    #########################
    ####  Normalization  ####
    #########################
    
    # Normalization for RNA composition (i.e. prevents genes from appearing falsely downregulated)
    
    # Creates effective library sizes
    
    y_puya_minus51 <- calcNormFactors(y_puya_minus51) # calculates normalization factors to scale raw library sizes, TMM is default method
    
    y_puya_minus51$samples # shows library sizes and normalization factors
    
    #######################################################
    ######## Multidimensional Scaling (MDS) Plot  #########
    #######################################################
    
    # gene.selection = 'common' selects the same genes for all comparisons
    
    # method='logFC' required with gene.selection argument
    
    # top = 500 plots the top 500 genes (this is the default)
    
    pdf("MDSplot_500genes_PUYA_minus51_colored_by_habitat.pdf")
    par(cex.axis = 1.1, cex = 1.25)
    mds500_puya_minus51 <- plotMDS(y_puya_minus51, gene.selection = 'common', top = 500, method='logFC', cex=1, 
                           col=c(rep("blue",5), rep("yellow",5)))
    #xlim=c(-2.5,2.5), xlab = "MDS axis 1",
    #ylim=c(-2.5,2.5), ylab = "MDS axis 2")
    legend("topright", cex = 1, 
           legend=c("fresh","sulfur"), 
           col=c("blue","yellow"), 
           pch = c(15),
           pt.cex=1.75)
    print(mds500_puya_minus51)
    dev.off()
    
    pdf("MDSplot_10000genes_PUYA_minus51_colored_by_habitat.pdf")
    par(cex.axis = 1.1, cex = 1.25)
    mds10000_puya_minus51 <- plotMDS(y_puya_minus51, gene.selection = 'common', top = 10000, method='logFC', cex=1, 
                           col=c(rep("blue",5), rep("yellow",5)))
    #xlim=c(-2.5,2.5), xlab = "MDS axis 1",
    #ylim=c(-2.5,2.5), ylab = "MDS axis 2")
    legend("topright", cex = 1, 
           legend=c("fresh","sulfur"), 
           col=c("blue","yellow"), 
           pch = c(15),
           pt.cex=1.75)
    print(mds10000_pichu)
    dev.off()
    
    #######################################
    ####### Graphing Library Size #########
    #######################################
    
    # create color vector for graph
    colors_puya_minus51 <- c(rep("blue",5),rep("yellow",5))
    
    # las=2 makes label text perpendicular to axis
    par(las=2)
    
    # generate plot
    pdf("library_size_PUYA_minus51.pdf")
    library_puya_minus51 <- barplot(y_puya_minus51$samples$lib.size, col=colors_puya_minus51, 
                             cex.axis=.8, ylab = "Library Size", xlab = "Sample")
    legend("topright", cex = 0.82, 
           legend=c("fresh", "sulfur"),
           bty = "n", pt.cex=1.5, 
           fill=c("blue", "yellow"))
    print(library_puya_minus51)
    dev.off()
    
    ####################################
    ########## Design Matrix ###########
    ####################################
    
    fac_puya_minus51 <- c(rep("fresh",5), rep("sulfur",5)) # use this to get groups
    
    fac_puya_minus51 <- factor(fac_puya_minus51) # renames fac names with numbers
    
    design_puya_minus51 <- model.matrix(~0+fac_puya_minus51) # 0 is the intercept
    
    colnames(design_puya_minus51) <- levels(fac_puya_minus51)
    
    design_puya_minus51
    
    ###################################################
    ####### Estimating Dispersions - CR method ########
    ###################################################
    
    # Uses GLM instead of qCML method because testing multiple factors here
    
    # Uses Cox-Reid profile-adjusted likelihood
    
    # estimates common dispersion and tagwise dispersions in one run
    y_puya_minus51 <- estimateDisp(y_puya_minus51, design_puya_minus51)
    
    # common dispersions
    y_puya_minus51$common.dispersion
    
    # tagwise (gene-specific dispersions)
    summary(y_puya_minus51$tagwise.dispersion)
    
    # estimated prior degrees of freedom
    y_puya_minus51$prior.df
    
    ###########################
    ####### BCV Plot ##########
    ###########################
    
    # BCV is the coefficient of variation with which the (unknown) true abundance of the gene varies between replicate RNA samples
    
    pdf("BCV_PUYA_minus51.pdf")
    BCV_puya_minus51 <- plotBCV(y_puya_minus51)
    print(BCV_puya_minus51)
    dev.off()
    
    ###############################################
    ####### DGE - Quasi-Likelihood F-test #########
    ###############################################
    
    # quasi-likelihood F-test better for bulk RNA-seq data because stricter error rate control, accounts for
    # uncertainty in dispersion estimation
    
    group_puya_minus51 <- factor(c(1,1,1,1,1,2,2,2,2,2))
    
    design_puya_minus51 <- model.matrix(~0+group_puya_minus51)
    
    # QL model representing the study design fitted to the data
    fit_puya_minus51 <- glmQLFit(y_puya_minus51, design_puya_minus51)
    
    # tests use FDR < 0.05
    
    # compare FRESH (1) vs SULFUR (2)
    qlf_puya_minus51 <- glmQLFTest(fit_puya_minus51, contrast=c(-1,1))
    fresh.vs.sulf_puya_minus51 <- topTags(qlf_puya_minus51, n = 5000, p.value=0.05)
    fresh.vs.sulf_puya_minus51 <- data.frame(fresh.vs.sulf_puya_minus51)
    require(xlsx)
    write.csv(x=fresh.vs.sulf_puya_minus51, file = "fresh.vs.sulf_PUYA_minus51.csv")
    
    ########################################################################################################################
    #### TACOTALPA #########################################################################################################
    ########################################################################################################################
    
    ###################################
    ####  Create grouping factor  #####
    ###################################
    
    # to be incorporated later into the DGEList
    
    group_taco = c(rep("fresh",6),rep("sulfur",6))
    
    ##########################################################
    ####  Remove rows where counts = 0 across all samples ####
    ##########################################################
    
    dim(taco)
    # 31791 genes before removing the 0 counts
    
    total_taco <- taco[rowSums(taco) > 0,] # trim out any rows with no reads counts
    
    dim(total_taco)
    # 24988 genes after removing the 0 counts
    
    ###################################
    ####  Creating DGEList object  ####
    ###################################
    
    y_taco <- DGEList(counts=total_taco, group=group_taco)
    
    #########################
    ####  Normalization  ####
    #########################
    
    # Normalization for RNA composition (i.e. prevents genes from appearing falsely downregulated)
    
    # Creates effective library sizes
    
    y_taco <- calcNormFactors(y_taco) # calculates normalization factors to scale raw library sizes, TMM is default method
    
    y_taco$samples # shows library sizes and normalization factors
    
    #######################################################
    ######## Multidimensional Scaling (MDS) Plot  #########
    #######################################################
    
    # gene.selection = 'common' selects the same genes for all comparisons
    
    # method='logFC' required with gene.selection argument
    
    # top = 500 plots the top 500 genes (this is the default)
    
    pdf("MDSplot_500genes_TACO_colored_by_habitat.pdf")
    par(cex.axis = 1.1, cex = 1.25)
    mds500_taco <- plotMDS(y_taco, gene.selection = 'common', top = 500, method='logFC', cex=1, 
                            col=c(rep("blue",6), rep("yellow",6)))
    #xlim=c(-2.5,2.5), xlab = "MDS axis 1",
    #ylim=c(-2.5,2.5), ylab = "MDS axis 2")
    legend("topright", cex = 1, 
           legend=c("fresh", "sulfur"), 
           col=c("blue","yellow"), 
           pch = c(15),
           pt.cex=1.75)
    print(mds500_taco)
    dev.off()
    
    pdf("MDSplot_10000genes_TACO_colored_by_habitat.pdf")
    par(cex.axis = 1.1, cex = 1.25)
    mds10000_taco <- plotMDS(y_taco, gene.selection = 'common', top = 10000, method='logFC', cex=1, 
                              col=c(rep("blue",6), rep("yellow",6)))
    #xlim=c(-2.5,2.5), xlab = "MDS axis 1",
    #ylim=c(-2.5,2.5), ylab = "MDS axis 2")
    legend("topright", cex = 1, 
           legend=c("fresh", "sulfur"), 
           col=c("blue","yellow"), 
           pch = c(15),
           pt.cex=1.75)
    print(mds10000_taco)
    dev.off()
    
    #######################################
    ####### Graphing Library Size #########
    #######################################
    
    # create color vector for graph
    colors_taco <- c(rep("blue",6),rep("yellow",6))
    
    # las=2 makes label text perpendicular to axis
    par(las=2)
    
    # generate plot
    pdf("library_size_TACO.pdf")
    library_taco <- barplot(y_taco$samples$lib.size, col=colors_taco, 
                             cex.axis=.8, ylab = "Library Size", xlab = "Sample")
    legend("topright", cex = 0.82, 
           legend=c("fresh", "sulfur"),
           bty = "n", pt.cex=1.5, 
           fill=c("blue", "yellow"))
    print(library_taco)
    dev.off()
    
    ####################################
    ########## Design Matrix ###########
    ####################################
    
    fac_taco <- c(rep("fresh",6), rep("sulfur",6)) # use this to get groups
    
    fac_taco <- factor(fac_taco) # renames fac names with numbers
    
    design_taco <- model.matrix(~0+fac_taco) # 0 is the intercept
    
    colnames(design_taco) <- levels(fac_taco)
    
    design_taco
    
    ###################################################
    ####### Estimating Dispersions - CR method ########
    ###################################################
    
    # Uses GLM instead of qCML method because testing multiple factors here
    
    # Uses Cox-Reid profile-adjusted likelihood
    
    # estimates common dispersion and tagwise dispersions in one run
    y_taco <- estimateDisp(y_taco, design_taco)
    
    # common dispersions
    y_taco$common.dispersion
    
    # tagwise (gene-specific dispersions)
    summary(y_taco$tagwise.dispersion)
    
    # estimated prior degrees of freedom
    y_taco$prior.df
    
    ###########################
    ####### BCV Plot ##########
    ###########################
    
    # BCV is the coefficient of variation with which the (unknown) true abundance of the gene varies between replicate RNA samples
    
    pdf("BCV_TACO.pdf")
    BCV_taco <- plotBCV(y_taco)
    print(BCV_taco)
    dev.off()
    
    ###############################################
    ####### DGE - Quasi-Likelihood F-test #########
    ###############################################
    
    # quasi-likelihood F-test better for bulk RNA-seq data because stricter error rate control, accounts for
    # uncertainty in dispersion estimation
    
    group_taco <- factor(c(1,1,1,1,1,1,2,2,2,2,2,2))
    
    design_taco <- model.matrix(~0+group_taco)
    
    # QL model representing the study design fitted to the data
    fit_taco <- glmQLFit(y_taco, design_taco)
    
    # tests use FDR < 0.05
    
    # compare FRESH (1) vs SULFUR (2)
    qlf_taco <- glmQLFTest(fit_taco, contrast=c(-1,1))
    fresh.vs.sulf_taco <- topTags(qlf_taco, n = 5000, p.value=0.05)
    fresh.vs.sulf_taco <- data.frame(fresh.vs.sulf_taco)
    require(xlsx)
    write.csv(x=fresh.vs.sulf_taco, file = "fresh.vs.sulf_TACO.csv")
    
    ####################################################################################################################################################
    #### Separate significantly upregulated and downregulated genes ####################################################################################
    ####################################################################################################################################################
    
    # Pichucalco
    pichu_up = rownames(fresh.vs.sulf_pichu[fresh.vs.sulf_pichu$logFC > 0,])
    
    pichu_down = rownames(fresh.vs.sulf_pichu[fresh.vs.sulf_pichu$logFC < 0,])
    
    # Puyacatengo
    puya_up = rownames(fresh.vs.sulf_puya_minus51[fresh.vs.sulf_puya_minus51$logFC > 0,])
    
    puya_down = rownames(fresh.vs.sulf_puya_minus51[fresh.vs.sulf_puya_minus51$logFC < 0,])
    
    # Tacotalpa
    taco_up = rownames(fresh.vs.sulf_taco[fresh.vs.sulf_taco$logFC > 0,])
    
    taco_down = rownames(fresh.vs.sulf_taco[fresh.vs.sulf_taco$logFC < 0,])
    
    #######################
    #### Venn Diagrams ####
    #######################
    
    # Upregulated
    universe.up <- unique(c(pichu_up, puya_up, taco_up))
    
    GroupA <- universe.up %in% pichu_up
    GroupB <- universe.up %in% puya_up
    GroupC <- universe.up %in% taco_up
    
    input.df <- data.frame(Pichu=GroupA, Puya=GroupB, Taco=GroupC)
    head(input.df)
    
    pdf("upregulated.pdf")
    a <- vennCounts(input.df)
    vennDiagram(a)
    print(a)
    dev.off()
    
    Q1all.up <- universe.up[which(input.df["Pichu"] == T & input.df["Puya"] == T & input.df["Taco"] == T)]
    
    # Downregulated
    universe.down <- unique(c(pichu_down, puya_down, taco_down))
    
    GroupD <- universe.down %in% pichu_down
    GroupE <- universe.down %in% puya_down
    GroupF <- universe.down %in% taco_down
    
    input.df2 <- data.frame(Pichu=GroupD, Puya=GroupE, Taco=GroupF)
    head(input.df2)
    
    pdf("downregulated.pdf")
    b <- vennCounts(input.df2)
    vennDiagram(b)
    print(b)
    dev.off()
    
    Q1all.down <- universe.down[which(input.df2["Pichu"] == T & input.df2["Puya"] == T & input.df2["Taco"] == T)]
    
    ###################################
    #### Getting BLAST annotations ####
    ###################################
    
    BLAST <- read.csv("Table S2 - Poecilia mexicana Annotations.csv")
    
    BLAST_annot <- subset(BLAST, select = c(gene.ID, Subject.sequence.ID, Protein.annotations, gene.name, E.value))
    
    # Upregulated
    Q1all.UP <- data.frame(Q1all.up)
    
    # add column heading that matches the blast output (Table S2 - Poecilia mexicana Annotations.csv)
    col_headings <- c("gene.ID")
    names(Q1all.UP) <- col_headings
    
    # extract swissprot accessions (subject sequence IDs) and Protein annotations from Blast output based on transcript IDs (query sequence IDs)
    All_candidates.up <- merge(BLAST_annot[, c("gene.ID", "Subject.sequence.ID", "Protein.annotations", "gene.name", "E.value")], Q1all.UP, by = "gene.ID", all.y = TRUE)
    write.csv(x=All_candidates.up, file = "All_candidates.up.csv")
    
    # Downregulated
    Q1all.DOWN <- data.frame(Q1all.down)
    
    # add column heading that matches the blast output (Table S2 - Poecilia mexicana Annotations.csv)
    col_headings <- c("gene.ID")
    names(Q1all.DOWN) <- col_headings
    
    # extract swissprot accessions (subject sequence IDs) and Protein annotations from Blast output based on transcript IDs (query sequence IDs)
    All_candidates.down <- merge(BLAST_annot[, c("gene.ID", "Subject.sequence.ID", "Protein.annotations", "gene.name", "E.value")], Q1all.DOWN, by = "gene.ID", all.y = TRUE)
    write.csv(x=All_candidates.down, file = "All_candidates.down.csv")

