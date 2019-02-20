# Pipeline used to identify structural variants

## 1. Set up SVMerge Configuration file

Here is an example SVMerge config file (configfile)
All cutoffs and scoring criteria were suggested by the SVMerge manual

	project=project # The name of your project directory
	name=sample_name # Sample name to be used in output files
	version=version # Subdirectory for SV calling analysis
	svdir=sv_calls # Subdirectory for filtering and merging
	chrOther= <list of all scaffolds here>
	projdir=/full/path/to/project/project # Full project path
	callerlist=pindel sec breakdancer # SV callers to integrate
	bam=/path/to/file.bam # Full path to BAM file
	bai=/path/to/file.bam.bai # Full path to BAM index file
	bamdir=/path/to/bamDirectory/ # Directory containing chromosome BAMs
	chrrefdir=/path/to/ref/chromDir/ # Dirctory containing chromosome FASTAs
	reffile=/path/to/ref/ref.fa # FASTA file with reference chromosomes
	BDscore=25 # BreakDancerMax score cutoff
	BDrs=2 # BreakDancerMax min. supporting RPs
	PDscore=30 # Pindel score cutoff
	PDsupports=10 # Pindel min. supporting reads
	bedexe=/path/to/BEDTools/bin/intersectBed # BEDTools intersect binary
	overlapexe=/path/to/SVMerge/overlapExons.pl # Included with SVMerge

## 2. Create project directory

Use the makeNewProject.sh script (included with SVMerge) to set up the project directory and subfolders

	makeNewProject.sh configfile

## 3. Run BreakDancer

	bam2cfg.pl bamfile > bd.config

For each scaffold, run BreakDancerMax using:

	BreakDancerMax.pl -f -o [scaff#] bd.config > name.scaff#.max

The output will be one file per scaffold that contains SV calls

## 4. Run Pindel

Set up a configuration file that contains the path to the bam, the insert size for the library (can be taken from the bd.config output from Breakdancer), and the sample name.
Here is an example file named pindel_conf:

	/data/bamfile.bam  400  sample_name

For each scaffold, run

	pindel -f ref.fa -i pindel_conf -c [scaff#] -o [sample_name]

## 5. Run SECluster

SECluster uses paired-end reads with one end mapped and one end unmapped to find candidate large insertions.
To filter out these reads from your BAM file(s), for each scaffold, use the samfilter.sh script included with SVMerge

	samfilter.sh bamfile [scaffold] 1

This produces sam files for each scaffold called [scaffold].se.sam
These files are used as input into SECluster (included with SVMerge)

	SECluster.pl -f [scaffold].se.sam -q 20 -m 5 -c 5 -r [scaffold] > name.[scaffold].clusters

q=read mapping quality [20 is the recommended value]
m=minimum reads to form a cluster [5 is the recommended value]
c=minimum reads in each forward and reverse cluster [5 is the recommended value]

## 6. Filter and merge calls

Filter and merge calls from the different SV callers
From your main working directory, use the filterAndMerge.pl script included with SVMerge

	filterAndMerge.pl -c configfile

This creates a merged, raw SV set: /project/version/svdir/merged/sample.merged.tab

## 7. De novo local assemblies and alignments

Create a configuration file using the coord2config.pl script included with SVMerge

	cd /project/version/svdir/localAssemblies/config/
	coord2config.pl -c configfile -b sample.merged.tab

The output from this are files names sv.config.1, sv.config2, etc.

Run the local assemblies and alignments using the svAssemble.pl script included with SVMerge
Depending on the number of sv.config files (which depends on the number of SVs in your raw set), it would be worthwhile to run this as an array or loop

	cd /project/version/svdir/localAssemblies/
	svAssemble.pl config/sv.config.1 <path_to_samtools>
	svAssemble.pl config/sv.config.2 <path_to_samtools>
	...

## 8. Parse alignments

The alignments are parsed for contigs which provided evidence for SV breakpoints
(e.g. a contig with a deletion will align with a large gap)

First, prepare the input files with the splitFile.sh script included with SVMerge
This splits the merged SV file into separate files (sample.merged.tab.1, sample.merged.tab.2, etc.) containing 500 SVs each

	cd /project/version/svdir/localAssemblies/
	splitFile.sh ../merged/sample.merged.bed 500

Next, parse the alignments with the runParser.pl script included with SVMerge

	cd /project/version/svdir/localAssemblies/
	runParser.pl -s sample.merged.tab.1 -c configfile > alignparse.1
	$EXEDIR/runParser.pl -s sample.merged.tab.2 -c configfile > alignparse.2
	...

Next, concatenate the files, then parse the output
The parseBoundary.pl script is included with SVMerge

	cd /project/version/svdir/localAssemblies/.
	cat alignparse.* > alignparse
	cd /project/version/svdir
	parseBoundary.pl -a localAssemblies/alignparse -o final/sv.final

The resulting files are as follows (found in /project/version/svdir/final/):

	sv.final.rank1.tab # Highest confidence set
	sv.final.rank2.tab # Local assembly results ambiguous
	sv.final.rank3.tab # Local assembly showed no breakpoints
	sv.final.small.tab # Local assembly SV found is <100bp
	sv.final.query.tab # SV is in a region with unusually high coverage

## 9. Merge final call set

This step is to identify overlapping SV calls (also called complex SVs)
The mergeFinalCalls.pl script is included with SVMerge

	cd /project/version/svdir/final/
	mergeFinalCalls.pl -f sv.final.rank1.tab -c configfile > sv.final.rank1.merged.tab
