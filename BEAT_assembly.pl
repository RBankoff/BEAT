#!/usr/bin/perl -w
use FileHandle;
use Getopt::Long;
use List::MoreUtils qw(uniq);

#####################
##  BEAT_Complete  ##
## Mod: 02/25/2016 ##
#####################

############################
##  Hard-coded arguments  ##
############################

#Number of Query files to be BLASTED
my $numQ = $ARGV[0];

#Number of Short-Read Files to be read in
my $numSR = $ARGV[1];

#########################
## flag initialization ##
#########################

#declare mandatory data-related flags
my $query = '';
my $sr = '';
my $chrom = '';
my $ref = '';

#declare flags specific to configuration mode
my $config = '';
my $type = '';
my $usr = '';
my $job_name = '';
my $run_local = '';
my $mode = '';
my $list = '';
my $species = '';
my $querylist = '';

#declare optional flags (defaults are preset if no input is registered)
my $samtools_qflag = '';
my $e_value = '';
my $perc_identity = '';
#my @SRs;

#load flags
GetOptions(
	"query=s{$numQ}" => \@query, 
	"usr=s" => \$usr, 
	"first" => \$first, 
	"chrom=s{$numQ}" => \@chrom, 
	"ref=s" => \$ref, 
	"config=s" => \$config, 
        "run_local=s" => \$run_local, 
        "type=s" => \$type, 
        "mode=s" => \$mode,
	"samtools_qflag=s" => \$samtools_qflag,
	"e_value=s" => \$e_value,
	"perc_identity=s" => \$perc_identity,
	"job_name=s" => \$job_name,
	"list=s" => \$list,
	"sr=s{$numSR}" => \@SRs,
	"species=s" => \$species,
	"querylist=s" => \$querylist
);


if($list){
	open (LIST_FILE, "<", $list) or die "Couldn't open list file $list, $!\n";
	while (my $line = <LIST_FILE>){
		chomp $line;
		#print "LINE: $line\n";
		push @SRs, $line;
	}
	close LIST_FILE;
}
if($querylist){
	open (QLIST_FILE, "<", $querylist) or die "Couldn't open query list file $querylist, $!\n";
	while (my $line = <QLIST_FILE>){
		chomp $line;
		#print "LINE: $line\n";
		push @query, $line;
	}
	close QLIST_FILE;
}


################################
## Global variable initiation ##
################################

#Short_read file arrays
my @odds;
my @evens;

#Output file arrays
my @blastfull;
my @forward_outs;
my @reverse_outs;

#Config file details
my @config;
my $config_header;
my @discard;

#client program settings
my $real_samtools_q;
my $real_evalue;

#populates the two arrays of even- and odd-input Short Read Filenames
push(@evens, $SRs[$_*2+1])
for 0..int(@SRs/2)-1;
push(@odds, $SRs[$_*2])
for 0..int(@SRs/2)-1;

#Hash table for odds, causes each @odds value to be associated with a irrelevent key
my %params = map { $_ => 1 }@odds;

#####################
##		   ##
##	MAIN	   ##
##                 ##
#####################

#check for programs necessary to run the analysis
&DEPENDENCIES();

sub DEPENDENCIES{
=pod
	my @dependencies = ('samtools', 'bwa', 'qsub', 'blastn', 'makeblastdb');
	foreach $d(@dependencies){
		my $check = qx(which $d);
		if ($check =~ /^\/usr\/bin\/which:\sno/){
			if ($d =~ m/blastn|makeblastdb/){
				$check = qx(which $d);
				if ($check =~ /^\/usr\/bin\/which:\sno/){
					die "Couldn't load module blast+. Check that it is installed, and try manually loading the module\n\n";
				}
			}
			else{
				$check = qx(which $d);
                                if ($check =~ /^\/usr\/bin\/which:\sno/){
                                        die "Couldn't load module $d. Check that it is installed, and try manually loading the module\n\n";
                                }
			}
		}
	}
=cut
	
	unless($config){
                unless ($run_local){
                        die "No config flag given. To run locally, use the --run_local flag[Ed note: not implemented in 1.0.0, you must use a queuing manager]\n\n";
                }
        }
	open (CONFIG, "<", $config) or die "Couldn't open config file, $!\n";
        while (my $line = <CONFIG>){
                chomp $line;
                push @config, $line;
        }
	close CONFIG;
        $config_header = join("\n", @config);
  	
	#check if user has entered valid ID if using a qsub system (pipeline will break otherwise)
	if ($type =~ m/qsub/i){
	        unless ($usr =~	/[A-Z0-9]/i){
                        die "No	valid user given for \"-usr\" flag, necessary to run as	a queue managed pipeline\n\n";
                }
        }
	
	if ($mode =~ m/complete/){
		unless ($job_name =~ m/[A-Z0-9]/){
			die "Invalid or missing job name $job_name, $!\n";
		}
	}

        #check reference presence/index status
        unless ($ref =~ m/[A-Z].*/i){
                die "No reference provided!\n\n";
        }
	unless (-e $ref){
		die "Reference file $ref not found at path given.\n";
	}
	my $ref_samtools_index = $ref . ".fai";
        my $ref_bwa_index = $ref . ".bwt";
        unless (-e $ref_samtools_index){
                `samtools faidx $ref`;
        }
	unless (-e $ref_bwa_index){
                my $indexing_script = $ref . ".bwa_index_script";
                open (BWA_IDX, ">", $indexing_script) or die "Couldn't make indexing script $indexing_script, $!\n";
                print BWA_IDX "$config_header\n";
                print BWA_IDX "module load bwa\n";
                print BWA_IDX "bwa index $ref\n";
                if ($type =~ m/qsub/){
                        `qsub $indexing_script`;
                }

#Example of how other queue-managing scripts should be inserted into the source code
=pod
		elsif ($type =~ m/other_format/i){
			`other_format $indexing_script`;
		}
=cut
                push @discard, $indexing_script;
        }

	#samtools filtering options
	if ($samtools_qflag =~ /\d/){
                $real_samtools_q = $samtools_qflag;
        }
	else{
                $real_samtools_q = 1;
        }
	
	#blast evalue settings
	if ($e_value =~ /\d/){
		$real_evalue = $e_value;
        }
        else{
		$real_evalue = 50;
        }

	#Mode_triggering
	if ($mode =~ m/fast/i){
		&FAST_TIDY();
	}
	elsif ($mode =~ m/complete/i){
		&COMPLETE();
	}
}

sub COMPLETE{
	my @pre_merged;
	my $SR_counter = 0;

	#############################
	###    MAPPING/SORTING    ###
	#############################

	foreach $O(@odds){
		my $SRQ;
		$SRQ = $O;
		$SRQ .= ".map_script";
		print "SRQ: $SRQ\n";
		my $local_bam_base = (split(/\./, $O))[0];
		my $unsorted = $local_bam_base;
		$unsorted .= ".bam";
                my $local_bam_sorted_base = $local_bam_base . ".sorted";
		my $local_bam_sorted = $local_bam_sorted_base . ".bam";
		open (MAP_MAKER, ">", $SRQ) or die "Couldn't make mapping script $SRQ, $!\n";
		print MAP_MAKER "$config_header\n";
		print MAP_MAKER "module load samtools\nmodule load bwa\n";
		print MAP_MAKER "bwa mem $ref $O $evens[$SR_counter] | samtools view -Sb -F4 -q $real_samtools_q -o $unsorted -\n"; 
		print MAP_MAKER "samtools sort $unsorted $local_bam_sorted_base\n";
		print MAP_MAKER "samtools index $local_bam_sorted\n";
		push @discard, $SRQ;
		push @pre_merged, $local_bam_sorted;
		if ($type =~ m/qsub/){
			`qsub $SRQ`;
		}
		$SR_counter++;
	}
	if ($type =~ m/qsub/i){
		my $qstatchecker = qx(qstat -u $usr | wc -l);
                if ($qstatchecker > 6){
			while($qstatchecker > 6){
                      		$qstatchecker = qx(qstat -u $usr | wc -l);
                        	sleep(1);
                	}
       		}
	}

	############################
	###   MERGING/INDEXING   ###
	############################
	unless ($numSR == 2){
		my $all_bams = join(' ', @pre_merged);
		my $merging_script = $job_name . ".merge_script";
		my $merged_name = $job_name . ".merged.bam";
		open (MERGING_SCRIPT, ">", $merging_script) or die "Couldn't make merging script $merging_script, $!\n";
		print MERGING_SCRIPT "$config_header\n";
		print MERGING_SCRIPT "module load samtools\n";
		print MERGING_SCRIPT "samtools merge $merged_name $all_bams\n";
		print MERGING_SCRIPT "samtools index $merged_name\n";
		push @discard, $merging_script;
	
		if ($type =~ m/qsub/){
			`qsub $merging_script`;
		        my $qstatchecker = qx(qstat -u $usr | wc -l);
		        if ($qstatchecker > 6){
		                while($qstatchecker > 6){
		                        $qstatchecker = qx(qstat -u $usr | wc -l);
		                        sleep(1);
		                }
		        }
		}
	}
	
	print "Assembly completed. To query your data further, use \"BEAT consensus\". Goodbye!\n";
=pod
        foreach $disc(@discard){
                `rm $disc`;
        }
=cut
	exit;
}

		
sub FAST_TIDY{
	#if there is a first flag given in the argument line
	if ($first){
		foreach $SR (@SRs){
			my $SRQ;
			 if ($SR =~ m/.gz$/){
                                $SRQ = $SR;
                                $SRQ =~ s/.gz/.dbscript/;
                        }
                        else{
                                $SRQ = $SR;
                                $SRQ =~ s/.fastq/.dbscript/;
                        }
			open (DB_MAKER, ">", $SRQ) or die "Couldn't make database script, $!\n"; 
			print DB_MAKER "$config_header\n";
			print DB_MAKER "module load blast+\n";
			print DB_MAKER "perl Streamline_FormatDB.pl $SR\n";
			if ($type =~ m/qsub/i){
				`qsub $SRQ`;
			}
			push @discard, $SRQ;
		}

		if ($type =~ m/qsub/i){
			my $qstatchecker = qx(qstat -u $usr | wc -l);
			if ($qstatchecker > 6){
				while($qstatchecker > 6){
					$qstatchecker = qx(qstat -u $usr | wc -l);
					sleep(1);
				}
			}
		}
	}
	&FAST_BLAST();
}


sub FAST_BLAST{
	#counter variable setup
	my $SRcounter = 0;

	foreach $GzedSRGs (@SRs){
		my $s = $GzedSRGs;
		$s =~ s/.gz//;
		my @localquery = @query;
		my $name;
		foreach $lq (@localquery){
			
			#strips extensions off of files (e.g. CDH23.fasta = CDH23[0] fasta[1])
			my @splitter = split (/\./, $lq);
			my $local = $s . "." . $splitter[0];
            		$name = "$local" . ".blast_script";
            		
			# make a blast .qsub file for each query within a directory for each short read file
            		my $fh = FileHandle->new();
            		open($fh, ">", $name) or die "Not able to create first blast script, $! \n";
            		
			#User-specified system commands to be written to each blast script
            		print $fh "$config_header\n";
			print $fh "module load blast+\n";
			
			if ($perc_identity =~ /\d/){
				#print specified percentage identity query-specific blast commands
				print $fh "blastn -query $lq -db $s" . " -perc_identity $perc_identity" . " > " . "$local" . ".blast\n";
			}
			else{
				#print default identity query-specific blast commands
                                print $fh "blastn -query $lq -db $s" . " > " . "$local" . ".blast\n";
			}
			#outfile name array
			my $blasted = "blasted." . "$local";
        	    	push @blastfull, $blasted;

			if(exists($params{$s})){
             			print $fh "perl ./blast2fq.pl " . "$local" . ".blast" . " $real_evalue " . $odds[$SRcounter] . " " . $evens[$SRcounter] . " $blasted" . "\n";
			}
	
        	    	else{
        	     		print $fh "perl ./blast2fq.pl " . "$local" . ".blast" . " $real_evalue " . $odds[($SRcounter-1)] . " $evens[($SRcounter-1)]" . " $blasted" . "\n";
        	    	}

			if ($type =~ m/qsub/){
				`qsub $name `;
			}
      		}

		if(exists($params{$s})){
      			$SRcounter++;
		}
		push @discard, $name;

	}

	if($type =~ m/qsub/){
		my $qstatchecker = qx(qstat -u $usr | wc -l);
		if ($qstatchecker > 6){
			while($qstatchecker > 6){
				$qstatchecker = qx(qstat -u $usr | wc -l);
				sleep(1);
			}
		}
	}
	
	my @revblast = @blastfull;

	foreach $revb(@revblast){
		$revb = "$revb" . ".flipReads.fastq";
	}

	foreach $blf(@blastfull){
		$blf = "$blf" . ".originalReads.fastq";
	}

	foreach $revb1(@revblast){
		my @localnames = @query;
		foreach $l(@localnames){
			my @splt = split (/\./, $l);
			$l = $splt[0];
		}
		my @splitter = split (/\./, $revb1);
        	my %split_hash = map { $_ => 1 }@splitter;
                foreach $ln(@localnames){
                        if (exists($split_hash{$ln})){
                                my $reverse_output = $ln . ".REVERSE.fastq";
                               	`cat $revb1 >> $reverse_output `;
                               	push @reverse_outs, $reverse_output;
                               	last;
                        }

=pod
			else{
				shift @localnames;
			}
=cut
		}

	}

	foreach $blf1(@blastfull){
		my @localnames = @query;
		foreach $l(@localnames){
			my @splt = split (/\./, $l);
			$l = $splt[0];
		}
		my @splitter = split (/\./, $blf1);
                my %split_hash = map { $_ => 1 }@splitter;
                foreach $ln(@localnames){
                        if (exists($split_hash{$ln})){
                                my $forward_output = $ln . ".FORWARD.fastq";
                                `cat $blf1 >> $forward_output `;
                                push @forward_outs, $forward_output;
                                last;
                        }
=pod
			else{
				shift @localnames;
			}
=cut
		}
	}

&FAST_MAP();
}

sub FAST_MAP{

	my @localnames = uniq @forward_outs;
	my @local_reverse = uniq @reverse_outs;
	my $finalcounter = 0;
	my @local_list;
	my @master_bams;
	foreach $l(@localnames){
		my $local_bam_base = (split(/\./, $l))[0];
		my $local_bam = $local_bam_base . ".bam";
		my $local_bam_presorted = $local_bam_base . ".sorted";
		my $local_bam_sorted = $local_bam_presorted . ".bam";
#		my $local_bam_rmdupd = $local_bam_presorted . ".rmdup.bam";
		push @master_bams, $local_bam_sorted;

		#my $local_chr = $chrom[$finalcounter];
		my $local_script_name = $l . ".bwa_script";
		open (BWA, ">", $local_script_name) or die "Couldn't open bwa script $local_script_name, $!\n";
		print BWA "$config_header\n";
		print BWA "module load bwa\nmodule load samtools\n";
		print BWA "bwa mem $ref $l $local_reverse[$finalcounter] | samtools view -Sb -F4 -q $real_samtools_q - >$local_bam \n"; 
		print BWA "samtools sort $local_bam $local_bam_presorted \n";
#		print BWA "samtools rmdup $local_bam_sorted $local_bam_rmdupd \n";
		print BWA "samtools index $local_bam_sorted \n";
		close BWA;
		push @discard, $local_script_name;
		$finalcounter++;

		if ($type =~ m/qsub/){
			`qsub $local_script_name`;
		}
		
		my $local_list_holder = qx(perl Entrez_fetch.pl $local_bam_base $species);
		chomp $local_list_holder;
		my @fixer = split(/\t/, $local_list_holder);
		$fixer[0] .= ".fasta";
		my $ready_for_consensus = join("\t", @fixer);
		my $OUT_LIST_FOR_CONSENSUS = $local_bam_base . ".consensus_list";
		open (OUT_LIST_FOR_CONSENSUS, ">", $OUT_LIST_FOR_CONSENSUS) or die "Couldn't make out list for consensus $OUT_LIST_FOR_CONSENSUS, $!\n";
		print OUT_LIST_FOR_CONSENSUS "$ready_for_consensus\n";
		close OUT_LIST_FOR_CONSENSUS;
		push @local_list, $OUT_LIST_FOR_CONSENSUS;
	}

	if($type =~ m/qsub/){
		my $qstatchecker = qx(qstat -u $usr | wc -l);
		if ($qstatchecker > 6){
			while($qstatchecker > 6){
				$qstatchecker = qx(qstat -u $usr | wc -l);
				sleep(1);
			}
		}
	}

	#print "Assembly completed.\n";
	my $consensus_counter = 0;
	foreach $master(@master_bams){
		my $base = (split(/\./, $master))[0];
		my $consensus_script = $base . ".consensus_script";
		open (CONSENSUS_SCRIPT, ">", $consensus_script) or die "Couldn't make consensus script $consensus_script, $!\n";
		print CONSENSUS_SCRIPT "$config_header\n";
		print CONSENSUS_SCRIPT "perl BEAT_consensus.pl $local_list[$consensus_counter] $master $base $ref \n";
		$consensus_counter++;		
		if ($type =~ m/qsub/){
			`qsub $consensus_script`;
		}
	}
	
	if($type =~ m/qsub/){
		my $qstatchecker = qx(qstat -u $usr | wc -l);
		if ($qstatchecker > 6){
			while($qstatchecker > 6){
				$qstatchecker = qx(qstat -u $usr | wc -l);
				sleep(1);
			}
		}
	}
	
=pod
	foreach $disc(@discard){
		`rm $disc`;
	}
=cut
	exit;
}
