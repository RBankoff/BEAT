BEAT v 0.9.2 Documentation

##################
## Introduction ##
##################

BEAT is a freeware, copyleft-compliant (GNU GPL3, see LICENSE.txt) bioinformatics pipeline developed by R Bankoff, with significant inputs from L Kistler, E Lauterbur, M Jerjos, and B Hohman at the Pennsylvania State University. BEAT's scripts are written in perl, and may be freely modified to fit the user's needs. The purpose of the program is to reconstruct targeted areas of a genome by aligning short read files from an unassembled species to a relative with a well annotated genome OR orthologous sequence and extracting the reads which map to the coordinates of orthologous loci. It is structured as a pair of core programs (BEAT_assembly.pl and BEAT_consensus.pl) that can either be accessed individually by the user, or through a central portal program called "BEAT.pl”. 

BEAT_assembly.pl: This is the script which directs the file splitting, indexing, and mapping of short reads for both tracks of BEAT. First, BEAT subdivides input short read files into small (default 1GB) chunks in order to take advantage the time benefits of parallel processing. There are two tracks to BEAT_assembly.pl:

Track 1: targeted - retrieve only short reads which align to a predetermined query sequence from a raw read file using NCBI’s blast+; map matching reads to the query. 
Track 2: complete - map all reads against an annotated reference genome from a related species; can be used for open-ended querying based on known coordinates of orthologs.

BEAT_consensus.pl: The resultant map (.bam) files from either track of BEAT_assembly can be fed into this script, which will generate a consensus sequence for the whole query sequence (targeted track) or loci specified by coordinates on an annotated reference (complete track).

What follows is the outline and general use guide to the command-line release version of the pipeline *without* the use of the BEAT.pl portal, whose operation should be self-explanatory because of its interactive nature.

##################################
## Dependencies and Permissions ##
##################################

None of the scripts need root access to run, but most of them require the installation of three key dependencies to run properly (or at all), which may require root permissions to install depending on your hardware configuration. The dependencies are:

1) Burrows-Wheeler Aligner (BWA; Li & Durbin, 2009; http://bio-bwa.sourceforge.net/)
2) Samtools, v0.1.19 (Li et al., 2009; http://samtools.sourceforge.net/)
3) Entrez Direct (Kans, 2013; http://www.ncbi.nlm.nih.gov/books/NBK179288/)

BEAT checks that these programs are installed and in the user's path before executing; if they are not, BEAT will throw an error to the effect of "edirect not found" and exit. It is highly recommended that users add the directory containing BEAT to their path. This can be accomplished with the following commands in a bash terminal:

PATH = $PATH:~/Downloads/BEAT_0_9_2/; 
export PATH;

##############################
## Basic Usage instructions ##
##############################

To use BEAT in the interactive mode, simply enter “perl BEAT.pl” into a bash terminal in the working directory. For command-line instructions, the arguments and flags are detailed below.

BEAT_assembly.pl usage:
perl BEAT_assembly.pl [arg1] [arg2] [mandatory flags] (optional flags)

#Arguments must be given in order
arg1: number of queries for BEAT_targeted; if using BEAT_complete, enter 1 
arg2: number of paired-end short-read files (2 minimum)

#Mandatory flags; these must be given, but order doesn’t matter 

FLAG		EXAMPLE				DESCRIPTION
—sr		SR_1_1.fastq SR_1_2.fastq	filenames of FASTQ-formatted paired-end short-read files 
					OR 
—list 		List_of_SR_files.txt		filename of list of short-reads in local directory 

—query		Query1.fa Query2.fa 		filename(s) of 
					OR 
-querylist	List_of_query_files.fa		filename of list of queries in local directory 

-usr		abc123				cluster-recognized username for qsubbing (TORQUE) clusters 
-chrom		chr3				name of the chromosome for each query in order entered 
-ref		hg19.fa				name of reference genome in fasta format
-config		config.qsub			name of configuration file with cluster header
-type 		qsub				set up to receive “qsub”, but can be modified (see below) 
-mode		target OR complete		defines which track the user will use, targeted or complete (see below)

#Optional flags; these may be specified or left at default values

FLAG			EXAMPLE		DEFAULT		DESCRIPTION
-samtools_qflag		20		1		defines custom quality flag value
-e_value		100		50		minimum e value to keep in blast alignment
-perc_identity		85		undefined	minimum percentage sequence identity to in blast alignment

IMPORTANT NOTE: To ensure proper results if using the TORQUE Resource Manager, you MUST set up your BEAT_assembly.pl job within a qsub script.

BEAT_consensus.pl usage:
perl BEAT_consensus.pl [arg1] [arg2] [arg3] [arg4]

#Arguments must be given in order
	NAME		EXAMPLE		DESCRIPTION
arg1:	Listfile	List.txt	list file in format specified in format section below	
arg2:	BAMfile		Master.bam	mapped .bam file from BEAT_assembly 
arg3:	Outfile		Out_folder	suffix to be added to query name for results directory
arg4:	Reference	hg19.fa		FASTA-formatted reference 

NOTE: Unlike BEAT_assembly.pl, BEAT_consensus.pl can be run either locally or as a qsub job without issue.

#############################
## BEAT formatting scripts ##
#############################

BEAT_splitter.pl
Purpose: A script to subdivide FASTQ files over 2GB into smaller chunks (default 1GB) that can be run in parallel. After subdivision, BEAT_splitter.pl outputs a properly formatted file listing the subdivisions to be input under the BEAT_MAIN.pl --list option. 
Usage: 


BEAT_batch.pl
Purpose:
Usage:


blast2fq.pl
Purpose: A script to 
Usage:


Consensus_mpileup.pl
Purpose:
Usage:


fq2faReads.pl
Purpose:
Usage:


fqselect.pl
Purpose:
Usage:


####################################
## Compatibility and Optimization ##
####################################

The BEAT suite is optimized for use on a linux cluster running the TORQUE Resource Manager (http://www.adaptivecomputing.com/products/open-source/torque/), but it can be run on a multithreading-capable local linux or unix system as well. To determine if you are using a system with the TORQUE Resource Manager installed, enter the command "qstat" in the terminal; if a list of running and queued jobs is returned, you are. This package is accompanied by a template file ("template_qsub.qsub"), which you may use, modify or replace to alter the default settings of each subdivided job. If you would like to adapt BEAT for other resource management systems, there is a built-in switch called "$type" in BEAT_MAIN.pl which has only one option at present ("qsub"); adding another option that mimicks the use of "$type" along with a system-specific template file in place of "template_qsub.qsub" should enable relatively easy conversion.   

BEAT is designed to intake FASTQ-formatted paired-end reads of indeterminate length, but a key strength of the pipeline is its ability to paralellize jobs to maximize efficiency. We therefore recommend that you use 

###############################
## Entrez Direct Integration ##
###############################

BEAT uses the Entrez Direct suite of tools from the National Center for Biotechnology Information (NCBI) to fetch data and sequences relevant to user's tasks using the Entrez_fetch.pl script. The script is designed to run as a subroutine within BEAT_consensus entirely autonomously of the user, but the script can also be used to directly query and refine data retrieved from any of the NCBI's databases that are available to Entrez Direct users. The syntax for commands to Entrez_fetch.pl is as follows:

perl Entrez_fetch.pl [1st argument: Gene symbol or name] [second argument: Taxonomic unit]

This will produce a listfile in the format needed by BEAT_consensus.pl to generate a consensus for the targeted query/reference.

#######################
## BEAT File Formats ##
#######################

Standard formats:

.bam: BAM files are the binary versions of SAM files, and cannot be read by humans without being passed through the samtools view interpreter.

.fa[sta]: a text-based format for representing nucleotide sequences.
Example:
>GENE_NAME
ACGTCCTGACACTGATGTATGCCACAATGTAACGTCCTGACACTGATGTATGCCACAATGTA

.fastq: a text-based format for storing read sequences along with their PHRED-formatted quality scores.
Example:
@SEQ_ID
ACGTCCTGACACTGATGTATGCCACAATGTAACGTCCTGACACTGATGTATGCCACAATGTA
+
!!''(%%%+*******!!$$CCGF928))))))'''!*.15HHHHHH(H!!!$%$$!'',.$

.qsub: A generic suffix to denote a script that is meant to be queued by the TORQUE Resource Management system.
Example:
#PBS -l walltime=4:00:00
#PBS -l nodes=1:ppn=1
#PBS -l pmem=2gb
cd $PBS_O_WORKDIR

.sam: Samtools alignment specification, see (https://samtools.github.io/hts-specs/SAMv1.pdf) for full specifications.

BEAT-specific formats:

Listfile format:
The listfile for BEAT_consensus [arg1] is an internal file which can be generated automatically by using the Entrez_fetch.pl script (see above), or can be created manually with the following format:

FIELD:	[1]			[2]					[3]						[4]
[name of .gbl file]\t[chromosome/scaffold name]\t[start coordinates along chromosome/scaffold]\t[end coordinates along chromosome/scaffold]\n

where “\t” indicates a tab and “\n” indicates a newline. NOTE: for “targeted” track jobs, input 1 for field [3] and the length of the sequence for field [4].

Scriptfile formats:
.mapping_script: An internally-generated queueable script that contains mapping instructions for a given pair of short reads. 

.merging_script: An internally-generated queueable script that merges and indexes all short reads mapped by the merging script. 

Outfile formats:

.consensus.fa: FASTA-formatted consensus for specified query/region

.out: The position-wise mapping scores from samtools mpileup used by BEAT_consensus to extract targeted consensus sequences.
Example:
scaffold_774	72731	G	9	..,,,,,,.	FFHHH@HH6
scaffold_774	72732	T	9	GGgggg,gG	HFHHHEHH6
