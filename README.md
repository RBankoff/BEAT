BEAT v 1.0.0 Documentation

#################
# Introduction #
#################

BEAT is a freeware, copyleft-compliant (GNU GPL3, see LICENSE.txt) bioinformatics pipeline developed by R Bankoff, with significant inputs from L Kistler, E Lauterbur, M Jerjos, and B Hohman at the Pennsylvania State University. BEAT's scripts are written in perl, and may be freely modified to fit the user's needs. The purpose of the program is to reconstruct targeted areas of a genome by aligning short read files from an unassembled species to a relative with a well annotated genome OR orthologous sequence and extracting the reads which map to the coordinates of orthologous loci. It is structured as a pair of core programs (BEAT_assembly.pl and BEAT_consensus.pl) that can either be accessed individually by the user, or through a central portal program called "BEAT.pl”. 

BEAT_assembly.pl: This is the script which directs the file splitting, indexing, and mapping of short reads for both tracks of BEAT. First, BEAT subdivides input short read files into small (default 1GB) chunks in order to take advantage the time benefits of parallel processing. There are two tracks to BEAT_assembly.pl:

Track 1: fast - Retrieve only short reads which align to a predetermined query sequence from a raw read file using NCBI’s blast+; map matching reads to the query. 

Track 2: complete - Map all reads against an annotated reference genome from a related species; can be used for open-ended querying based on known coordinates of orthologs.

BEAT_consensus.pl: 
The resultant map (.bam) files from either track of BEAT_assembly can be fed into this script, which will generate a consensus sequence for the whole query sequence (fast track) or loci specified by coordinates on an annotated reference (complete track). This step is automatically completed for fast track queries.

What follows is the outline and general use guide to the command-line release version of the pipeline *without* the use of the BEAT.pl portal, whose operation should be self-explanatory because of its interactive nature.

################################
# Dependencies and Permissions #
################################

None of the scripts need root access to run, but most of them require the installation of three key dependencies to run properly (or at all), which may require root permissions to install depending on your hardware configuration. The dependencies are:

1) Burrows-Wheeler Aligner (BWA; Li & Durbin, 2009; http://bio-bwa.sourceforge.net/)

2) Samtools, v0.1.19 (Li et al., 2009; http://samtools.sourceforge.net/)

3) Entrez Direct (Kans, 2013; http://www.ncbi.nlm.nih.gov/books/NBK179288/)

BEAT checks that these programs are installed and in the user's path before executing; if they are not, BEAT will throw an error to the effect of "edirect not found" and exit. It is highly recommended that users add the directory containing BEAT to their path. This can be accomplished with the following commands in a bash terminal:

`PATH = $PATH:~/Downloads/BEAT_v1_0_0/; export PATH`

###################################
# Compatibility and Optimization #
###################################

The BEAT suite is optimized for use on a linux cluster running the TORQUE Resource Manager (http://www.adaptivecomputing.com/products/open-source/torque/), but it can be run on a multithreading-capable local linux or unix system as well. To determine if you are using a system with the TORQUE Resource Manager installed, enter the command "qstat" in the terminal; if a list of running and queued jobs is returned, you are. This package is accompanied by a template file ("template_qsub.qsub"), which you may use, modify or replace to alter the default settings of each subdivided job. If you would like to adapt BEAT for other resource management systems, there is a built-in switch called "$type" in BEAT_MAIN.pl which has only one option at present ("qsub"); adding another option that mimicks the use of "$type" along with a system-specific template file in place of "template_qsub.qsub" should enable relatively easy conversion.   

BEAT is designed to intake FASTQ-formatted paired-end reads of indeterminate length, but a key strength of the pipeline is its ability to paralellize jobs to maximize efficiency. We therefore recommend that you use BEAT_splitter.pl to break up large short read files prior to assembling your reads to a reference genome, provided you have a cluster capable of running multiple simultaneous instances of BWA mem/samtools.

##############################
# Basic Usage instructions ##
##############################

To use BEAT in the interactive mode, simply enter “perl BEAT.pl” into a bash terminal in the working directory. For command-line instructions, the arguments and flags are detailed below.

##BEAT_assembly.pl usage:

`perl BEAT_assembly.pl [arg1] [arg2] [mandatory flags] <optional flags>`

### Arguments must be given in order:

arg1: number of queries for BEAT_fast; if using BEAT_complete, enter 1 

arg2: number of paired-end short-read files (2 minimum)

###Mandatory flags; these must be given, but order doesn’t matter

FLAG | EXAMPLE | DESCRIPTION
:-----: | :---: | :--------- 
-sr| SR_1_1.fastq SR_1_2.fastq |filenames of FASTQ-formatted paired-end short-read files 
 | OR |
-list |	List_of_SR_files.txt |filename of list of short-reads in local directory 
-query|	Query1.fa Query2.fa |filename(s) of 
 | OR |
-querylist | List_of_query_files.fa |filename of list of queries in local directory 
-usr | abc123 |cluster-recognized username for qsubbing (TORQUE) clusters 
-chrom | chr3 |name of the chromosome for each query in order entered 
-ref | hg19.fa | name of reference genome in fasta format
-config	| config.qsub |name of configuration file with cluster header
-type | qsub |set up to receive “qsub”, but can be modified (see below) 
-mode |	fast OR complete |defines which track the user will use, targeted or complete (see below)
-species | Homo sapiens |the name of the species used as the reference file

###Semi-optional flags; these must be included if certain preconditions are not fulfilled

FLAG | DESCRIPTION
:-----: | :--------- 
-first| indicates to BEAT_fast to make BLAST databases of short read files if none are present


###Optional flags; these may be specified or left at default values

FLAG|EXAMPLE|DEFAULT|DESCRIPTION
:-------:|:----:|:------:|:------------
-samtools_qflag|20|1|defines custom quality flag value
-e_value|100|50|minimum e value to keep in blast alignment
-perc_identity|85|undefined|minimum percentage sequence identity to in blast alignment

####IMPORTANT NOTE: To ensure proper results if using the TORQUE Resource Manager, you MUST set up your BEAT_assembly.pl job within a qsub script.

###Example run walkthorugh

Included in the BEAT github repository are a set of example files to demonstrate a properly set up run on the BEAT_fast track. Everything needed to run the pipeline is in there except for a reference genome, which for these data can be retrived from [here](https://www.ncbi.nlm.nih.gov/nuccore/568815589?report=fasta). The demo short read files are a greatly subsetted group of reads from an aye-aye sequence repository (SRA066444; Perry *et al.*, 2012) that been spiked with reads that match to the human TMC1 orthologue in a previous run of BLAST+. To modify or inspect the parameters of the job, simply opne the Demo_Master.qsub file in a text editor. Two parameters need to be changed in the example job file: your username, after the `-usr` flag, and the name of the reference you download after the `-ref` flag. To run the job, simply enter `qsub Demo_Master.qsub` on a TORQUE queue manager compliant system.

##BEAT_consensus.pl usage:

`perl BEAT_consensus.pl [1] [2] [3] [4]`

###Arguments must be given in order

ORDER|NAME|EXAMPLE|DESCRIPTION
:---:|:--------:|:------:|:--------------:
1|Listfile|List.txt|list file in format specified in format section below	
2|BAMfile|Master.bam|mapped .bam file from BEAT_assembly 
3|Outfile|Out_folder|suffix to be added to query name for results directory
4|Reference|hg38.fa|FASTA-formatted reference 

####NOTE: Unlike BEAT_assembly.pl, BEAT_consensus.pl can be run either locally or as a qsub job without issue.

###########################
# BEAT formatting scripts #
###########################

##BEAT_splitter.pl
###Purpose: 
A script to subdivide FASTQ files over 2GB into smaller chunks (default 1GB) that can be run in parallel. After subdivision, BEAT_splitter.pl outputs a properly formatted file listing the subdivisions to be input using the BEAT_MAIN.pl --list option. 
###Usage: 
`perl BEAT_splitter.pl [1] [2] [3] <4>`

Order|Name|Example|Description
:---:|:--------:|:------:|:--------------:
1|Short read list|List.txt|short reads to be split in a newline-delimited file 
2|Job name|Monday_map|name of prefix to be associated with the split files 
3|Outfile|Out_list.txt|name of file containing the split short read file names
4|Split value|15000000|optional: number of FASTQ-formatted lines in each daughter file

####NOTE: 1 GB is approximately 15,313,576 lines; if you do not enter a 4th argument, this is the value that will be used

##blast2fq.pl/fqselect.pl
###Purpose: 
A pair of scripts that reassociate sequence data from reads matched to the sequence of interest by BLAST+ with FASTQ quality data. fqselect4.pl is handed data by blast2fq.pl; running it alone is useless and thus arguments will not be specified.

###Usage: 
`perl blast2fq.pl [1] [2] [3] [4] [5]`

Order|Name|Example
:---:|:------:|:------------:
1|Blast results filename| file.blast
2|Minimum score| 50
3|Forward reads| s_1_1.fastq(.gz allowed)
4|Reverse reads| s_1_2.fastq(.gz allowed)
5|Out prefix| TMC1

##Consensus_mpileup.pl
###Purpose:
Generates the actual consensus sequences produced from samtools mpileup data, hard to use without going through BEAT_consensus.pl

###Usage:
`perl Consensus_mpileup.pl [1] [2] [3] [4]`

Order|Name|Example
:---:|:------:|:------------:
1|FASTA file name| TMC1.fasta
2|Chromosome name| 9
3|Mpileup file name| TMC1.mpileup
5|Outfolder prefix| TMC1 

##Streamline_formatDB.pl
###Purpose:
Subscript used to batch-process the creation of BLAST+ databases with the makeblastdb utility.
###Usage:
`perl Streamline_FormatDB.pl [name of file to be databased]`

##############################
# Entrez Direct Integration ##
##############################

BEAT uses the Entrez Direct suite of tools from the National Center for Biotechnology Information (NCBI) to fetch data and sequences relevant to user's tasks using the Entrez_fetch.pl script. The script is designed to run as a subroutine within BEAT_consensus entirely autonomously of the user, but the script can also be used to directly query and refine data retrieved from any of the NCBI's databases that are available to Entrez Direct users. The syntax for commands to Entrez_fetch.pl is as follows:

`perl Entrez_fetch.pl [1] [2] [3]`

Order|Name|Example
:---:|:------:|:------------:
1|Gene symbol (NCBI)| TMC1
2|Generic/species name| Homo
3|Download sequences?| (Y/n)


This will produce a listfile in the format needed by BEAT_consensus.pl to generate a consensus for the targeted query/reference, and can also be used to download sequences directly from NCBI.

#######################
# BEAT File Formats ##
#######################

##Standard formats:

###.bam: 
BAM files are the binary versions of SAM files, and cannot be read by humans without being passed through the samtools view interpreter.

###.fa/.fasta: 
a text-based format for representing nucleotide sequences.

Example:

`>GENE_NAME`

`ACGTCCTGACACTGATGTATGCCACAATGTAACGTCCTGACACTGATGTATGCCACAATGTA`

###.fq/.fastq: 
a text-based format for storing read sequences along with their PHRED-formatted quality scores.

Example:

`@SEQ_ID`

`ACGTCCTGACACTGATGTATGCCACAATGTAACGTCCTGACACTGATGTATGCCACAATGTA`

`+`

`!!''(%%%+*******!!$$CCGF928))))))'''!*.15HHHHHH(H!!!$%$$!'',.$`

###.qsub: 
A generic suffix to denote a script that is meant to be queued by the TORQUE Resource Management system.

Example:

`#PBS -l walltime=4:00:00`

`#PBS -l nodes=1:ppn=1`

`#PBS -l pmem=2gb`

`cd $PBS_O_WORKDIR`

###.sam: 
Samtools alignment specification, see (https://samtools.github.io/hts-specs/SAMv1.pdf) for full specifications.

##BEAT-specific formats:

###Listfile format:
The listfile for BEAT_consensus [arg1] is an internal file which can be generated automatically by using the Entrez_fetch.pl script (see above), or can be created manually with the following format, where each field is separated by a tab:

Field 1|Field 2|Field 3|Field 4
:---------:|:-------:|:---------:|:---------:
name of .gbl file|chromosome/scaffold name|start coordinates along chromosome/scaffold|end coordinates along chromosome/scaffold

####NOTE: for “targeted” track jobs, input 1 for field [3] and the length of the sequence for field [4].


###.mapping_script: 
An internally-generated queueable script that contains mapping instructions for a given pair of short reads. 

###.merging_script: 
An internally-generated queueable script that merges and indexes all short reads mapped by the merging script. 

###.consensus.fa: 
FASTA-formatted consensus for specified query/region

###.out: 
The position-wise mapping scores from samtools mpileup used by BEAT_consensus to extract targeted consensus sequences.

Example:

Ref name|Ref Coord|Ref ID|Read Depth|Read ID|Read Quality
:----------:|:------:|:---:|:---:|:---------:|:-----------:
scaffold_774	72731	G	9	..,,,,,,.	FFHHH@HH6
scaffold_774	72732	T	9	GGgggg,gG	HFHHHEHH6

