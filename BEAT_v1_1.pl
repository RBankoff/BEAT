#!/usr/bin/perl -w
use FileHandle;
use Getopt::Long;


my $program = $ARGV[0];
my $interactive = 1;

unless (scalar(@ARGV) > 0){
        die "Basic Exon Assembly Tool, v. 0.9.0\nNo arguments detected. For a list of available arguments, type \"perl BEAT.pl help\". Exiting...\n\n";
}

sub HELP{
	print "\nProgram: Basic Exon Assembly Tool, v. 1.0.0\n";
	print "Usage: perl BEAT.pl <command> [options]\n\n";
	print "Commands:\n";
	print "\t--Raw read mapping--\t\t\t\t\tSynopsis:\n";
	print "\t\traw_input_fast\t\trapid extraction of targeted sequences from raw read files\n";
	print "\t\traw_input_complete\tfully map read files to reference for random-access consensus generation\n\n";
	print "\t--Query and Consensus generation--\n";
	print "\t\tentrez_fetch\t\tretrieve full sequences from NCBI and generate coordinate files\n";
	print "\t\tconsensus\t\tgenerate consensus sequences and associated data for user-specified genes\n";
#	print "\t\tbatch\t\t\tgenerate consensus sequences for large numbers of sequences in parallel\n\n";
	print "\t\tsplit\t\t\tbreak large short read files into properly associated subfiles for more efficient parallelization\n\n";
	print "\t\thelp\t\t\tthis menu\n\n";
	print "NOTES:\n1) BEAT was developed for use on computing clusters with sufficient computational resources to handle sustained heavy loading.\nIf you're unsure if your system is equipped to handle it, consider not using the raw_input_complete command.\n\n2)All main pipeline scripts must be in your path in order for BEAT to function. Furthermore, you must have the following programs installed:\n\nsamtools\nbwa\nentrez direct\n\nFor more information on how to install these tools, see the README file.\n\n";  
	exit;
}

if ($program =~ m/help$|^h$|^-h$|^--h$/i){
	&HELP();
}

elsif ($program =~ m/^raw/i){
	if ($interactive == 1){
		print "This is the BEAT interactive job setup wizard. If you would like to simply enter your parameters into BEAT manually, use \nthe BEAT_assembly.pl script.\n\nContinue? (Y/n)\n";
		my $CONTINUE = <STDIN>;
		chomp $CONTINUE;
		if ($CONTINUE =~ m/n/i){
			die "Ok, exiting now...\n";
		}

		#global values
		my @config;
		my $QUEUE;
		my $MODE;
		my @Qs;
		my @SRs;
		my @CHRs;
		my $Q_counter_variable = 0;
		my $Sr_counter_variable = 0;		
		
		print "BEAT is optimized to handle jobs submitted through the TORQUE package manager at the moment. Is this the system you are working on?(Y/n)\n";
		my $TORQUE = <STDIN>;
		chomp $TORQUE;

		if ($TORQUE =~ m/y/i){
			printf "OK great! ";
			$QUEUE = "qsub";
		}

		elsif ($TORQUE =~ m/n/i){
			print "Are you using GNU Parallel? (Y/n)\n";
			my $GNU = <STDIN>;
			chomp $GNU;
			if ($GNU =~ m/y/i){
				$QUEUE = "GNU";
			}
			else{
			#for future use in developing alternate system calling strategies, modify this variable
				print "Ok, we can try to work with that. Enter the command that you would use to submit a job to your cluster (e.g. qsub):\n"; 
		 		$QUEUE = <STDIN>;
				chomp $QUEUE;
			}
		}

		print "Enter the names of the short read files you will be mapping today in the format SR1_1.fastq SR1_2.fastq..SRN_1.fastq SRN_2.fastq, or type \"list\" to enter the path to a file containing a list of your short reads. Reads must be paired-end; list forward and reverse reads back-to-back:\n";
		my $SR_NO = <STDIN>;
		chomp $SR_NO;

		#if file ends in .txt, check that it's not a list
		if ($SR_NO =~ m/.txt$/){
			print "Filename ends in .txt. Is this a list? (Y/n)\n";
			my $LIST_CHECK = <STDIN>;
			chomp $LIST_CHECK;
			if ($LIST_CHECK =~ m/y/i){
				my $SR_LIST_IN = $SR_NO;
				open (SR_LIST_IN, "<", $SR_LIST_IN) or die "Couldn't open List file $SR_LIST_IN, $!\n";			
				while (my $line = <SR_LIST_IN>){
					chomp $line;
					push @SRs, $line;
					$Sr_counter_variable++;
				}
				close Q_LIST_IN;
			}
			else{
				my @intermediate = split(/\s/, $SR_NO);
				if ((scalar(@intermediate))%2!=0){
					die "Odd number of short read inputs. BEAT is only set up to deal with paired reads at the moment, sorry!\n";
				}
				else{
					foreach $SRN(@intermediate){
						push @SRs, $SRN;
						$Sr_counter_variable++;
					}
				}
			}
		}
		
		#checks to see if the user has invoked the list function
		elsif ($SR_NO =~ m/list/i){
			print "Enter the name of the file containing the list of short read file names in a newline-delimited format.\n";
			my $SR_LIST_IN = <STDIN>;
			chomp $SR_LIST_IN;
			open (SR_LIST_IN, "<", $SR_LIST_IN) or die "Couldn't open List file $SR_LIST_IN, $!\n";			
			while (my $line = <SR_LIST_IN>){
				chomp $line;
				push @SRs, $line;
				$Sr_counter_variable++;
			}
			close SR_LIST_IN;
		}

		#assume that they've just entered the names
		else{
			my @intermediate = split(/\s/, $SR_NO);
				#check that values entered are A) paired end and B) non-zero
				if (scalar(@intermediate)%2!=0){
					die "Odd number of short read inputs. BEAT is only set up to deal with paired reads at the moment, sorry!\n";
				}
				elsif (scalar(@intermediate) == 0){
					die "No short read file names entered, exiting...\n";
				}
				else{
					foreach $SRN(@intermediate){
						push @SRs, $SRN;
						$Sr_counter_variable++;
					}
				}
		}

		
		print "Ok, you have loaded " . scalar(@SRs) . " short read files:\n";
		my $S_COUNTER = 1;
		foreach $S(@SRs){
			unless(-e $S){
				die "No such query file $S, $!\n";
			}	
			print "\#$S_COUNTER\t$S\n";
			$S_COUNTER++;
		}
		
		#check if makeblastdb has already been used on these short read files
		my $FIRST_SWITCH = 0;		
		print "Have these short read files already been indexed with the blast+ makeblastdb tool? (Y/n)\n";
		my $MAKEBLASTDB_CHECK = <STDIN>;
		chomp $MAKEBLASTDB_CHECK;
		if ($MAKEBLASTDB_CHECK =~ m/n/i){
			$FIRST_SWITCH++;
		}

		my @intermediate_Chrs;

		print "Ok. Enter the binomial or generic name of the organism you're using as a reference:\n";
		my $REF_ORG = <STDIN>;
		chomp $REF_ORG;
		foreach $Q (@Qs){
			my $BARE_Q;
			if ($Q =~ m/.fa$|.fasta$/i){
				$BARE_Q = (split(/\./, $Q))[0];
			}
			else{
				$BARE_Q = $Q;
			}
			my $response = qx(perl Entrez_fetch.pl $BARE_Q $REF_ORG);
			my $LOCAL_PRECHROM = (split(/\t/, $response))[1];
			$LOCAL_PRECHROM =~ s/chr//gi;
			push @CHRs, $LOCAL_PRECHROM;
		}


		#FAST TRACK STARTS HERE
		if ($program =~ m/fast$/i){
			print "Please enter the filenames for the query files, separated by spaces, or type \"list\" to enter the path to a file containing a newline-delimited list of your query file names. Queries must be in Fasta format.\n";
			my $Q_NO = <STDIN>;
			chomp $Q_NO;

			#if file ends in .txt, check that it's not a list
			if ($Q_NO =~ m/.txt$/){
				print "Filename ends in .txt. Is this a list? (Y/n)\n";
				my $LIST_CHECK = <STDIN>;
				chomp $LIST_CHECK;
				if ($LIST_CHECK =~ m/y/i){
					my $Q_LIST_IN = $Q_NO;
					open (Q_LIST_IN, "<", $Q_LIST_IN) or die "Couldn't open List file $Q_LIST_IN, $!\n";			
					while (my $line = <Q_LIST_IN>){
						chomp $line;
						push @Qs, $line;
						$Q_counter_variable++;
					}
					close Q_LIST_IN;
				}
				else{
					my @intermediate = split(/\s/, $Q_NO);
					foreach $QMN(@intermediate){
						push @Qs, $QMN;
						$Q_counter_variable++;
					}
				}
			}
		
			#checks to see if the user has invoked the list function
			elsif ($Q_NO =~ m/list/i){
				print "Enter the name of the file containing the list of query file names in a newline-delimited format.\n";
				my $Q_LIST_IN = <STDIN>;
				chomp $Q_LIST_IN;
				open (Q_LIST_IN, "<", $Q_LIST_IN) or die "Couldn't open List file $Q_LIST_IN, $!\n";			
				while (my $line = <Q_LIST_IN>){
					chomp $line;
					push @Qs, $line;
					$Q_counter_variable++;
				}
				close Q_LIST_IN;
			}
			
			#assume that they've just entered the names 
			else{
				my @intermediate = split(/\s/, $Q_NO);
				foreach $QMN(@intermediate){
					push @Qs, $QMN;
					$Q_counter_variable++;
				}
			}
			
			#Review statement
			print "Ok, you have loaded " . scalar(@Qs) . " query files:\n";
			my $Q_COUNTER = 1;
			foreach $Q(@Qs){
				#check that files are where they're supposed to be 
				unless(-e $Q){
					die "No such query file $Q, $!\n";
				}
				print "\#$Q_COUNTER\t$CHRs[($Q_COUNTER-1)]\t$Q\n";
				$Q_COUNTER++;
			}
		}

		#Collects reference genome name/directory, checks to make sure it exists
		print "Enter the path to the fasta-formatted reference genome:\n";
		my $REF = <STDIN>;
		chomp $REF;
		unless(-e $REF){
			die "No such file $REF, cannot proceed without a valid reference!\n";
		}
		
		print "Enter your system username (e.g. abc123):\n";
		my $USR = <STDIN>;
		chomp $USR;
		
		#import config file
		print "Enter the path to a file containing the header for your cluster-submitted script:\n";
		my $CONF = <STDIN>;
		chomp $CONF;
		unless(-e $CONF){
			die "No such file $CONF, cannot proceed without a valid header file!\n";
		}
		open (CONFIG, "<", $CONF) or die "Couldn't open header file, $!\n";
        	while (my $line = <CONFIG>){
                	chomp $line;
                	push @config, $line;
        	}
		close CONFIG;
        	my $config_header = join("\n", @config);
		
		#check if user wants to fine tune parameters
		print "Adjust mapping values for samtools view, blast+ e-values, or blast+ \%identity scores?(Y/n)\n";
		my $TINY_VALUES = <STDIN>;
		chomp $TINY_VALUES;
		my $SAM_Q;
		my $E_VALUE;
		my $ID_SCORE;

		#presence/absence switches for fine-tuning parameters
		my $TINY_SWITCH = 0;
		my $SAM_SWITCH = 0;
		my $E_SWITCH = 0;
		my $ID_SWITCH = 0;

		#accept fine-tune parameters from user
		if ($TINY_VALUES =~ m/y/i){
			$TINY_SWITCH++;
			print "Enter value for samtools view qvalue; if default, enter \"default\"\n";
			my $SAMTOOLS_QVALUE = <STDIN>;
			chomp $SAMTOOLS_QVALUE;
			if($SAMTOOLS_QVALUE !~ m/default/i){
				$SAM_Q = $SAMTOOLS_QVALUE;
				$SAM_SWITCH++;
			}
			print "Enter value for e-value cutoff; if default, enter \"default\"\n";
			my $EVAL = <STDIN>;
			chomp $EVAL;
			if($EVAL !~ m/default/i){
				$E_VALUE = $EVAL;
				$E_SWITCH++;
			}
			print "Enter value for sequence identity cutoff; if default, enter \"default\"\n";
			my $ID_VAL = <STDIN>;
			chomp $ID_VAL;
			if($ID_VAL !~ m/default/i){
				$ID_SCORE = $ID_VAL;
				$ID_SWITCH++;
			} 
		}

		#set mode variable
		if ($program =~ m/fast$/i){
			$MODE = "fast";
		}
		elsif ($program =~ m/complete$/i){
			$MODE = "complete";
		}
		
		#creat unique job name
		print "Enter the name of this job:\n";
		my $JOB_NAME = <STDIN>;
		chomp $JOB_NAME;

		my $MASTER_NAME = $JOB_NAME . "_Master.qsub";

		#assemble various variables for printing to the master script
		my $assembled_queries = join (' ', @Qs);
		my $assembled_shortrd = join (' ', @SRs);
		my $assembled_chromsm = join (' ', @CHRs);

		#master command
		my $COMMAND_STRING;
		if ($program =~ m/fast$/i){
			$COMMAND_STRING = "perl BEAT_assembly.pl $Q_counter_variable $Sr_counter_variable -query $assembled_queries -sr $assembled_shortrd -usr $USR -chrom $assembled_chromsm -ref $REF -config $CONF -mode $MODE -type $QUEUE -job_name $JOB_NAME -species $REF_ORG";
		}
		elsif ($program =~ m/complete$/i){
			$COMMAND_STRING = "perl BEAT_assembly.pl 1 $Sr_counter_variable -sr $assembled_shortrd -usr $USR -ref $REF -config $CONF -mode $MODE -type $QUEUE -job_name $JOB_NAME -species $REF_ORG";
		}
		#check if makeblastdb step is necessary
		if ($FIRST_SWITCH > 0){
			$COMMAND_STRING .= " -first";
		}
		
		#add in user-specified fine tuning, if applicable
		if ($TINY_SWITCH != 0){
			if ($SAM_SWITCH > 0 ){
				$COMMAND_STRING .= " -samtools_qflag $SAM_Q";
			}
			if ($E_SWITCH > 0){
				$COMMAND_STRING .= " -e_value $E_VALUE";
			}
			if ($ID_SWITCH > 0){
				$COMMAND_STRING .= " -perc_identity $ID_SCORE";
			}
		}
		
		#create master file
		open (MASTER, ">", $MASTER_NAME) or die "Couldnt make master file $MASTER_NAME, $!\n";
		print MASTER $config_header . "\n";
		print MASTER $COMMAND_STRING . "\n";
		close MASTER;
		
		#Autorun?
		print "Submit job to the queue manager? (Y/n)\n";
		my $SUBMIT = <STDIN>;
		chomp $SUBMIT;
		if ($SUBMIT =~ m/y/i){
			my $JOB_INFO = qx($QUEUE $MASTER_NAME);
			print "Job was submitted to the queue manager at $JOB_INFO; data processing will begin shortly. Goodbye!\n\n";
			exit;
		}
		else{
			print "OK! A queue-able job file named $MASTER_NAME was created and can be submitted to your queue manager at any time. Goodbye!\n";
			exit;
		}  
	}
}
elsif($program =~ m/entrez_fetch/i){
	my @geneids;	
	print "This program is used by the main BEAT assembly tracks to retrieve chromosomal coordinate information from NCBI. It can also be used to download sequence data for queries\n";  
	print "Enter the binomial or generic name of the organism you're using as a reference:\n";
	my $REF_ORG = <STDIN>;
	chomp $REF_ORG;
	print "Enter the NCBI GENEID names of the genes you want information on, or enter \"list\" to provide the name of a newline-delimited listfile:\n";
	my $GENES = <STDIN>;
	if ($GENES =~ m/list/i){
		print "Ok, enter the path to the listfile now:\n";
		my $LIST = <STDIN>;
		chomp $LIST;
		open (LIST, "<", $LIST) or die "Couldn't find listfile $LIST, $!\n";
		while (my $line = <LIST>){
			chomp $line;
			push @geneids, $line;
		}
		close LIST;
	}
	else{
		@geneids = split(/\s/, $GENES);
	}
	print "Would you like to download the sequence data for these genes as well?(Y/n)\n";
	my $DOWNLOAD = <STDIN>;
	chomp $DOWNLOAD;
		
	
	print "GENEID\tCHR\tSTART\tEND\n";
	foreach $GID (@geneids){
		my $response;
		if ($DOWNLOAD =~ m/y/i){
			$response = qx(perl Entrez_fetch.pl $GID $REF_ORG Y);
		}
		else{
			$response = qx(perl Entrez_fetch.pl $GID $REF_ORG);
		}
		print $response;
	}
	if ($DOWNLOAD =~ m/y/i){
		print "Exon-by-Exon information is in files named [GENEID].long.gb in this directory; Downloaded sequences can be found in files named $REF_ORG" . "_" . "$GID" . ".fasta in this folder. Goodbye!\n";
		}
	else{
		print "Exon-by-Exon information is in files named [GENEID].long.gb in this directory. Goodbye!\n";
	}
	
	exit;

}

elsif($program =~ m/consensus/i){

=pod
#BEAT_consensus.pl argumentative structure:
	my $in = $ARGV[0]; = list of 
	my $master_bam = $ARGV[1];
	my $basedir_name = $ARGV[2];
	my $ref = $ARGV[3];
=cut
	my @geneids;	
	print "This program produces consensus sequences from short reads that have been mapped to a reference with BEAT_complete.\n";  
	print "Enter the binomial or generic name of the organism you're using as a reference:\n";
	my $REF_ORG = <STDIN>;
	chomp $REF_ORG;
	print "Enter the NCBI GENEID names of the genes you want information on, or enter \"list\" to provide the name of a newline-delimited listfile:\n";
	my $GENES = <STDIN>;
	if ($GENES =~ m/list/i){
		print "Ok, enter the path to the listfile now:\n";
		my $LIST = <STDIN>;
		chomp $LIST;
		open (LIST, "<", $LIST) or die "Couldn't find listfile $LIST, $!\n";
		while (my $line = <LIST>){
			chomp $line;
			push @geneids, $line;
		}
		close LIST;
	}

	else{
		@geneids = split(/\s/, $GENES);
	}
	print "Enter the path to the master bam file (file produced by BEAT_complete) you'd like to extract a consensus from:\n";
	my $M_BAM = <STDIN>;
	chomp $M_BAM;

	print "Enter the path to the reference genome used in the production of the master bam file (file produced by BEAT_complete) you'd like to extract a consensus from:\n";
	my $REF = <STDIN>;
	chomp $REF;
	
	foreach $GID (@geneids){
		my $response = qx(perl Entrez_fetch.pl $GID $REF_ORG);
		chomp $response;
		my @fixer = split(/\t/, $response);
		$fixer[0] .= ".fasta";
		my $ready_for_consensus = join("\t", @fixer);
		my $OUT_LIST_FOR_CONSENSUS = $GID . ".consensus_list";
		open (OUT_LIST_FOR_CONSENSUS, ">", $OUT_LIST_FOR_CONSENSUS) or die "Couldn't make out list for consensus $OUT_LIST_FOR_CONSENSUS, $!\n";
		print OUT_LIST_FOR_CONSENSUS "$ready_for_consensus\n";
		close OUT_LIST_FOR_CONSENSUS;
		`perl BEAT_consensus_v1_1.pl $OUT_LIST_FOR_CONSENSUS $M_BAM $GID $REF`;
		print "Consensus sequences have been generated and can be found in the folder $GID.\n";
	}
	print "Goodbye!\n";
	exit;
}
=pod
elsif($program =~ m/batch/i){
	
}
=cut
elsif($program =~ m/split/i){
	my @SRs;
	
	print "This program enables you to break up large fastq-formatted short read files into smaller chunks that can be processed in parallel. To begin, enter the names of the short read files you will be splitting today in the format SR1_1.fastq SR1_2.fastq..SRN_1.fastq SRN_2.fastq, or type \"list\" to enter the path to a file containing a list of your short reads. Reads must be paired-end; list forward and reverse reads back-to-back:\n";
	my $SR_NO = <STDIN>;
	chomp $SR_NO;
	my $Sr_counter_variable = 0;
	my $list_switch = 0;
	my $SR_LIST_HOLDER;

	#if file ends in .txt, check that it's not a list
	if ($SR_NO =~ m/.txt$/){
		print "Filename ends in .txt. Is this a list? (Y/n)\n";
		my $LIST_CHECK = <STDIN>;
		chomp $LIST_CHECK;
		if ($LIST_CHECK =~ m/y/i){
			$SR_LIST_HOLDER = $SR_NO;
			$list_switch++; 
		}
		else{
			my @intermediate = split(/\s/, $SR_NO);
			if ((scalar(@intermediate))%2!=0){
				die "Odd number of short read inputs. BEAT is only set up to deal with paired reads at the moment, sorry!\n";
			}
			else{
				foreach $SRN(@intermediate){
					push @SRs, $SRN;
					$Sr_counter_variable++;
				}
			}
		}
	}
		
	#checks to see if the user has invoked the list function
	elsif ($SR_NO =~ m/list/i){
		print "Enter the name of the file containing the list of short read file names in a newline-delimited format.\n";
		my $SR_LIST_IN = <STDIN>;
		chomp $SR_LIST_IN;
		$SR_LIST_HOLDER = $SR_LIST_IN;
		$list_switch++; 
	}

	#assume that they've just entered the names
	else{
		my @intermediate = split(/\s/, $SR_NO);
			#check that values entered are A) paired end and B) non-zero
			if (scalar(@intermediate)%2!=0){
				die "Odd number of short read inputs. BEAT is only set up to deal with paired reads at the moment, sorry!\n";
			}
			elsif (scalar(@intermediate) == 0){
				die "No short read file names entered, exiting...\n";
			}
			else{
				foreach $SRN(@intermediate){
					push @SRs, $SRN;
					$Sr_counter_variable++;
				}
			}
	}
	
	
	print "Enter a job name to identify the subdivided short-read files:\n";
	my $SHORT_READ_NAME = <STDIN>;
	chomp $SHORT_READ_NAME;

	my $GB_SIZE;	
	print "Enter a size to restrict each sub-file to in Gigabytes, or type \"default\"(1GB):\n";
	my $GB_SIZE_IN = <STDIN>;
	chomp $GB_SIZE_IN;
	
	if ($GB_SIZE_IN =~ m/default/i){
		$GB_SIZE = 15313576;
	}
	else{
		$GB_SIZE = ($GB_SIZE_IN * 15313576);
	}
	my $list_to_split = $SHORT_READ_NAME . ".presplit_list";
	my $list_of_split = $SHORT_READ_NAME . ".split_list"; 	
	
	if ($list_switch == 0){	
		print "Ok, you have loaded " . scalar(@SRs) . " short read files:\n";
		my $S_COUNTER = 1;
		open (LIST_TO_SPLIT, ">", $list_to_split) or die "Couldn't open $list_to_split, $!\n";
		foreach $S(@SRs){
	#		unless(-e $S){
	#			die "No such query file $S, $!\n";
	#		}	
			print "\#$S_COUNTER\t$S\n";
			print LIST_TO_SPLIT "$S\n";
			$S_COUNTER++;
		}
		print "Splitting now. This may take some time...\n";
		`perl BEAT_splitter.pl $list_to_split $SHORT_READ_NAME $list_of_split $GB_SIZE`;
	}
	else{
		print "Splitting now. This may take some time...\n";
		`perl BEAT_splitter.pl $SR_LIST_HOLDER $SHORT_READ_NAME $list_of_split $GB_SIZE`;
	}
exit;
}














