# Pipeline used to generate candidate gene set

## 1. Trim with TrimGalore
    
Trim reads using TrimGalore default settings.
    
    while read line || [ -n "$line" ];
    do
    	read1=$(echo $line | awk '{print $1}')
    	read2=$(echo $line | awk '{print $2}')
    
    trim_galore --fastqc_args "--noextract --nogroup" --output_dir <path to output directory> --paired $read1 $read2
    
    done < samples.txt
    
    # Top rows of samples.txt
    sample1_read1.fastq.gz   sample1_read2.fastq.gz
    sample2_read1.fastq.gz   sample2_read2.fastq.gz
    sample3_read1.fastq.gz   sample3_read2.fastq.gz

## 2. Pre-process with Cufflinks gffread

Change the reference genome GFF file into a GTF file otherwise it won't work with the HISAT2 python script in the next step.
     
     gffread <gff_in> -T -o <gtf_out>

## 3. Index the reference genome with HISAT2

Make files using python scripts provided by HISAT2 for the --ss and --exon options.

Extract splice sites from reference GTF file prior to building the index:
    
    module load python/2.7.10
    hisat2_extract_splice_sites.py <gtf_in> > splice_sites.txt

Also extract exons from pmex GTF file prior to building the index:
    
    module load python/2.7.10
    hisat2_extract_exons.py <gtf_in> > exons.txt

Build an index using the reference genome.
    
    hisat2-build -f <reference_in> --ss splice_sites.txt --exon exons.txt <ht2_base>

## 4. Map using HISAT2

Map reads to the previously indexed reference genome.
    
    mkdir -p stats
    
    mkdir -p unmapped
    
    while read line || [ -n "$line" ];
    do
    	read1=$(echo $line | awk '{print $1}')
    	read2=$(echo $line | awk '{print $2}')
    	filename=$(echo $line | awk '{print $3}')
    
    hisat2 --phred33 -k 10 --met-file $filename.stats --rg-id $filename --rg SM:$filename --rg PL:illumina -p 1 --rna-strandness RF --fr \
    	--dta --un-conc-gz $filename.unmapped -x <hisat2-idx> -1 $read1 -2 $read2 -S $filename.sam
    
    done < samples.txt
    
    # Top rows of samples.txt
    sample1_read1.fq.gz   sample1_read2.fq.gz   sample1
    sample2_read1.fq.gz   sample2_read2.fq.gz   sample2
    sample3_read1.fq.gz   sample3_read2.fq.gz   sample3

## 5. Convert to BAM and sort using SAMtools
    
    while read line || [ -n "$line" ];
    do
    	sampleName=$(echo $line | awk '{print $1}')
    
    samtools view -b -h -S $sampleName.sam > $sampleName.bam
    
    samtools sort $sampleName.bam -o $sampleName.sorted.bam
    
    done < samples.txt
    
    # Top rows of samples.txt
    sample1
    sample2
    sample3

## 6. Generate GTF files using StringTie

Generate GTF files in a ballgown folder.
    
    while read line || [ -n "$line" ];
    do
    	sampleName=$(echo $line | awk '{print $1}')
    
    mkdir -p ballgown
    
    mkdir -p ballgown/$sampleName
    
    stringtie $sampleName.sorted.bam -o $sampleName.gtf -p 4 -G reference.gff -B -e
    
    done < samples.txt

    # Top rows of samples.txt
    sample1
    sample2
    sample3

## 7. Generate a gene counts matrix (gene_count_matrix.csv) using prepDE.py (from StringTie)

Run python script provided by StringTie to generate 2 CSV files containing the count matrices for genes and transcripts.

    module load python/2.7.10
    
    python prepDE.py -i ballgown -g gene_count_matrix.csv -t transcript_count_matrix.csv
        
## 8. Differential gene expression analyses using EdgeR

    data <- read.csv("gene_count_matrix.csv", row.names=1)
    
Change column names

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
   
Package installation

    source("https://bioconductor.org/biocLite.R")
    biocLite("edgeR")
    biocLite("limma")
    
Load packages

    require("limma")
    require("edgeR")
    
Create grouping factor to be incorporated later into the DGEList, grouped by drainage

    group = c(rep("Pich",12),rep("Puya",11),rep("Taco",12))
    
Remove rows where counts = 0 across all samples

    total <- data[rowSums(data) > 0,]

Creating DGEList object 

    y <- DGEList(counts=total, group=group)
    
Normalization, creates effective library sizes

    y <- calcNormFactors(y)
    
Multidimensional Scaling (MDS) Plot

    # top 500 genes, colored by drainage
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
    
    # top 10,000 genes, colored by drainage
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

    # top 500 genes, colored by sulfidic vs. non-sulfidic
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
    
    # top 10,000 genes, colored by sulfidic vs. non-sulfidic
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
    
Split dataframes by drainage

    # Pichucalco drainage
    pichu <- data[,1:12]
    
    # Puyacatengo drainage
    puya <- data[,c(23,22,21,20,19,18,17,16,15,14,13)]
    
    # Tacotalpa drainage
    taco <- data[,24:35]
    
---- PICHUCALCO DRAINAGE
    
Create grouping factor to be incorporated later into the DGEList

    group_pichu = c(rep("fresh",6),rep("sulfur",6))
    
Remove rows where counts = 0 across all samples

    total_pichu <- pichu[rowSums(pichu) > 0,]
    
Creating DGEList object

    y_pichu <- DGEList(counts=total_pichu, group=group_pichu)
    
Normalization, creates effective library sizes

    y_pichu <- calcNormFactors(y_pichu)

Multidimensional Scaling (MDS) Plot
    
    # top 500 genes, colored by sulfidic vs. non-sulfidic   
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
    
    # top 10,000 genes, colored by sulfidic vs. non-sulfidic   
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
    
Graph library size, colored by sulfidic vs. non-sulfidic

    colors_pichu <- c(rep("blue",6),rep("yellow",6))
    par(las=2)
    pdf("library_size_PICHU.pdf")
    library_pichu <- barplot(y_pichu$samples$lib.size, col=colors_pichu, 
                       cex.axis=.8, ylab = "Library Size", xlab = "Sample")
    legend("topright", cex = 0.82, 
           legend=c("fresh", "sulfur"),
           bty = "n", pt.cex=1.5, 
           fill=c("blue", "yellow"))
    print(library_pichu)
    dev.off()
    
Design Matrix

    fac_pichu <- c(rep("fresh",6), rep("sulfur",6))
    fac_pichu <- factor(fac_pichu)
    design_pichu <- model.matrix(~0+fac_pichu)
    colnames(design_pichu) <- levels(fac_pichu)
    
Estimating Dispersions - CR method
Uses GLM instead of qCML method because testing multiple factors.
Uses Cox-Reid profile-adjusted likelihood.

    y_pichu <- estimateDisp(y_pichu, design_pichu)
    
BCV Plot

    pdf("BCV_PICHU.pdf")
    BCV_pichu <- plotBCV(y_pichu)
    print(BCV_pichu)
    dev.off()
    
Differential Gene Expression - Quasi-Likelihood F-test

    group_pichu <- factor(c(1,1,1,1,1,1,2,2,2,2,2,2))
    design_pichu <- model.matrix(~0+group_pichu)
    fit_pichu <- glmQLFit(y_pichu, design_pichu)
    
    # compare SULFIDIC (1) to NON-SULFIDIC (-1)
    qlf_pichu <- glmQLFTest(fit_pichu, contrast=c(-1,1))
    fresh.vs.sulf_pichu <- topTags(qlf_pichu, n = 5000, p.value=0.05)
    fresh.vs.sulf_pichu <- data.frame(fresh.vs.sulf_pichu)
    require(xlsx)
    write.csv(x=fresh.vs.sulf_pichu, file = "fresh.vs.sulf_PICHU.csv")
    
---- PUYACATENGO DRAINAGE
    
Create grouping factor to be incorporated later into the DGEList

    group_puya = c(rep("fresh",6),rep("sulfur",5))
    
Remove rows where counts = 0 across all samples

    total_puya <- puya[rowSums(puya) > 0,]
    
Creating DGEList object

    y_puya <- DGEList(counts=total_puya, group=group_puya)
    
Normalization, creates effective library sizes

    y_puya <- calcNormFactors(y_puya)
    
Multidimensional Scaling (MDS) Plot
    
    # top 500 genes, colored by sulfidic vs. non-sulfidic
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
    
    # top 10,000 genes, colored by sulfidic vs. non-sulfidic
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
    
Took out Puyacatengo sample 51 because looks like an outlier

    puya_minus51 <- puya[,-5]
    
Create grouping factor to be incorporated later into the DGEList

    group_puya_minus51 = c(rep("fresh",5),rep("sulfur",5))
    
Remove rows where counts = 0 across all samples ####

    total_puya_minus51 <- puya_minus51[rowSums(puya_minus51) > 0,]
    
Creating DGEList object

    y_puya_minus51 <- DGEList(counts=total_puya_minus51, group=group_puya_minus51)
    
Normalization, creates effective library sizes

    y_puya_minus51 <- calcNormFactors(y_puya_minus51)
  
Multidimensional Scaling (MDS) Plot
    
    # top 500 genes, colored by sulfidic vs. non-sulfidic
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
    
    # top 10,000 genes, colored by sulfidic vs. non-sulfidic
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
    
Graphing Library Size

    colors_puya_minus51 <- c(rep("blue",5),rep("yellow",5))
    par(las=2)
    pdf("library_size_PUYA_minus51.pdf")
    library_puya_minus51 <- barplot(y_puya_minus51$samples$lib.size, col=colors_puya_minus51, 
                             cex.axis=.8, ylab = "Library Size", xlab = "Sample")
    legend("topright", cex = 0.82, 
           legend=c("fresh", "sulfur"),
           bty = "n", pt.cex=1.5, 
           fill=c("blue", "yellow"))
    print(library_puya_minus51)
    dev.off()
    
Design Matrix

    fac_puya_minus51 <- c(rep("fresh",5), rep("sulfur",5))
    fac_puya_minus51 <- factor(fac_puya_minus51)
    design_puya_minus51 <- model.matrix(~0+fac_puya_minus51)
    colnames(design_puya_minus51) <- levels(fac_puya_minus51)
    
Estimating Dispersions
Uses GLM instead of qCML method because testing multiple factors.
Uses Cox-Reid profile-adjusted likelihood.
    
    y_puya_minus51 <- estimateDisp(y_puya_minus51, design_puya_minus51)
    
BCV Plot
    
    pdf("BCV_PUYA_minus51.pdf")
    BCV_puya_minus51 <- plotBCV(y_puya_minus51)
    print(BCV_puya_minus51)
    dev.off()
    
Differential Gene Expression - Quasi-Likelihood F-test

    group_puya_minus51 <- factor(c(1,1,1,1,1,2,2,2,2,2))
    
    design_puya_minus51 <- model.matrix(~0+group_puya_minus51)

    # compare SULFIDIC (1) to NON-SULFIDIC (-1)
    qlf_puya_minus51 <- glmQLFTest(fit_puya_minus51, contrast=c(-1,1))
    fresh.vs.sulf_puya_minus51 <- topTags(qlf_puya_minus51, n = 5000, p.value=0.05)
    fresh.vs.sulf_puya_minus51 <- data.frame(fresh.vs.sulf_puya_minus51)
    require(xlsx)
    write.csv(x=fresh.vs.sulf_puya_minus51, file = "fresh.vs.sulf_PUYA_minus51.csv")
    
---- TACOTALPA DRAINAGE
    
Create grouping factor to be incorporated later into the DGEList
    
    group_taco = c(rep("fresh",6),rep("sulfur",6))
    
Remove rows where counts = 0 across all samples
    
    total_taco <- taco[rowSums(taco) > 0,]

Creating DGEList object
    
    y_taco <- DGEList(counts=total_taco, group=group_taco)
    
Normalization, creates effective library sizes
    
    y_taco <- calcNormFactors(y_taco)
    
Multidimensional Scaling (MDS) Plot
    
    # top 500 genes, colored by sulfidic vs. non-sulfidic
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
    
    # top 10,000 genes, colored by sulfidic vs. non-sulfidic
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
    
Graphing Library Size

    colors_taco <- c(rep("blue",6),rep("yellow",6))
    par(las=2)
    pdf("library_size_TACO.pdf")
    library_taco <- barplot(y_taco$samples$lib.size, col=colors_taco, 
                             cex.axis=.8, ylab = "Library Size", xlab = "Sample")
    legend("topright", cex = 0.82, 
           legend=c("fresh", "sulfur"),
           bty = "n", pt.cex=1.5, 
           fill=c("blue", "yellow"))
    print(library_taco)
    dev.off()
    
Design Matrix

    fac_taco <- c(rep("fresh",6), rep("sulfur",6))
    fac_taco <- factor(fac_taco)
    design_taco <- model.matrix(~0+fac_taco)
    colnames(design_taco) <- levels(fac_taco)
    
Estimating Dispersions
Uses GLM instead of qCML method because testing multiple factors.
Uses Cox-Reid profile-adjusted likelihood.
    
    y_taco <- estimateDisp(y_taco, design_taco)
    
BCV Plot

    pdf("BCV_TACO.pdf")
    BCV_taco <- plotBCV(y_taco)
    print(BCV_taco)
    dev.off()
    
Differential Gene Expression - Quasi-Likelihood F-test
    
    group_taco <- factor(c(1,1,1,1,1,1,2,2,2,2,2,2))
    design_taco <- model.matrix(~0+group_taco)
    fit_taco <- glmQLFit(y_taco, design_taco)
    # compare SULFIDIC (1) to NON-SULFIDIC (-1)
    qlf_taco <- glmQLFTest(fit_taco, contrast=c(-1,1))
    fresh.vs.sulf_taco <- topTags(qlf_taco, n = 5000, p.value=0.05)
    fresh.vs.sulf_taco <- data.frame(fresh.vs.sulf_taco)
    require(xlsx)
    write.csv(x=fresh.vs.sulf_taco, file = "fresh.vs.sulf_TACO.csv")
    
Separate significantly upregulated and downregulated genes 
    
    # Pichucalco
    pichu_up = rownames(fresh.vs.sulf_pichu[fresh.vs.sulf_pichu$logFC > 0,])
    pichu_down = rownames(fresh.vs.sulf_pichu[fresh.vs.sulf_pichu$logFC < 0,])
    
    # Puyacatengo
    puya_up = rownames(fresh.vs.sulf_puya_minus51[fresh.vs.sulf_puya_minus51$logFC > 0,])
    puya_down = rownames(fresh.vs.sulf_puya_minus51[fresh.vs.sulf_puya_minus51$logFC < 0,])
    
    # Tacotalpa
    taco_up = rownames(fresh.vs.sulf_taco[fresh.vs.sulf_taco$logFC > 0,])
    taco_down = rownames(fresh.vs.sulf_taco[fresh.vs.sulf_taco$logFC < 0,])
    
Venn Diagram
    
Upregulated genes

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
    
Downregulated genes

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
    
Get BLAST annotations
    
    BLAST <- read.csv("Table S2 - Poecilia mexicana Annotations.csv")
    
    BLAST_annot <- subset(BLAST, select = c(gene.ID, Subject.sequence.ID, Protein.annotations, gene.name, E.value))
    
Upregulated genes

    Q1all.UP <- data.frame(Q1all.up)
    
    # add column heading that matches the blast output (Table S2 - Poecilia mexicana Annotations.csv)
    col_headings <- c("gene.ID")
    names(Q1all.UP) <- col_headings
    
    # extract swissprot accessions (subject sequence IDs) and Protein annotations from Blast output based on transcript IDs (query sequence IDs)
    All_candidates.up <- merge(BLAST_annot[, c("gene.ID", "Subject.sequence.ID", "Protein.annotations", "gene.name", "E.value")], Q1all.UP, by = "gene.ID", all.y = TRUE)
    write.csv(x=All_candidates.up, file = "All_candidates.up.csv")

Downregulated genes

    Q1all.DOWN <- data.frame(Q1all.down)
    
    # add column heading that matches the blast output (Table S2 - Poecilia mexicana Annotations.csv)
    col_headings <- c("gene.ID")
    names(Q1all.DOWN) <- col_headings
    
    # extract swissprot accessions (subject sequence IDs) and Protein annotations from Blast output based on transcript IDs (query sequence IDs)
    All_candidates.down <- merge(BLAST_annot[, c("gene.ID", "Subject.sequence.ID", "Protein.annotations", "gene.name", "E.value")], Q1all.DOWN, by = "gene.ID", all.y = TRUE)
    write.csv(x=All_candidates.down, file = "All_candidates.down.csv")

