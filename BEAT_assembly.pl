#!/usr/bin/perl -w
use FileHandle;
use Getopt::Long;

###################
## BEAT_Complete ##
## Mod: 8/4/2015 ##
###################

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

#declare optional flags (defaults are preset if no input is registered)
my $samtools_qflag = '';
my $e_value = '';
my $perc_identity = '';


#load flags
GetOptions(
	"query=s{$numQ}" => \@query, 
	"usr=s" => \$usr, 
	"first" => \$first, 
	"sr=s{$numSR}" => \@SRs,  
	"chrom=s{$numQ}" => \@chrom, 
	"ref=s" => \$ref, 
	"config=s" => \$config, 
        "run_local=s" => \$run_local, 
        "type=s" => \$type, 
        "mode=s" => \$mode,
	"samtools_qflag=s" => \$samtools_qflag,
	"e_value=s" => \$e_value,
	"perc_identity=s" => \$perc_identity,
	"job_name=s" => \$job_name
);

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
	my @dependencies = ('samtools', 'bwa', 'qsub', 'blastn', 'makeblastdb');
	foreach $d(@dependencies){
		my $check = qx(which $d);
		if ($check =~ /^\/usr\/bin\/which:\sno/){
			if ($d =~ m/blastn|makeblastdb/){
				`module load blast+`;
				$check = qx(which $d);
				if ($check =~ /^\/usr\/bin\/which:\sno/){
					die "Couldn't load module blast+. Check that it is installed, and try manually loading the module\n\n";
				}
			}
			else{
				`module load $d`;
				$check = qx(which $d);
                                if ($check =~ /^\/usr\/bin\/which:\sno/){
                                        die "Couldn't load module $d. Check that it is installed, and try manually loading the module\n\n";
                                }
			}
		}
	}
	
	unless($config){
                unless ($run_local){
                        die "No config flag given. To run locally, use the --run_local flag\n\n";
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
                        die "No	valid user given for \"-usr\" flag, necessary to run as	a qsubbing pipeline\n\n";
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
		if ($O =~ m/.gz$/){
			$SRQ = $O =~ s/.gz/.map_script/;
		}
		else{
			$SRQ = $O =~ s/.fastq/.map_script/;
		} 
		my $local_bam_base = (split(/\./, $O))[0];
                my $local_bam_sorted_base = $local_bam_base . ".sorted";
		my $local_bam_sorted = $local_bam_sorted_base . ".bam";
		open (MAP_MAKER, ">", $SRQ) or die "Couldn't make mapping script $SRQ, $!\n";
		print MAP_MAKER "$config_header\n";
		print MAP_MAKER "module load samtools\nmodule load bwa\n";
		print MAP_MAKER "bwa mem $ref $O $evens[$SR_counter] | samtools view -Sb -F4 -q $real_samtools_q - | samtools sort -o $local_bam_sorted_base -";
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
	
	###########################
	###  DUPLICATE REMOVAL  ###
	###########################
		
	#gotta figure out cross platform (non-aliased) calling of the old samtools...

	############################
	###   MERGING/INDEXING   ###
	############################

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
	
	print "Assembly completed. To query your data further, use \"BEAT_consensus.pl\" or \"BEAT_batch.pl\" \n";
        foreach $disc(@discard){
                `rm $disc`;
        }
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
		foreach $ln(@localnames){
			if ($splitter[3] =~ /$ln/){
				my $reverse_output = $ln . ".REVERSE.fastq";
				`cat $revb1 >> $reverse_output `;
				push @reverse_outs, $reverse_output;
				last;
			}
			else{
				shift @localnames;
			}
		}

	}

	foreach $blf1(@blastfull){
		my @localnames = @query;
		foreach $l(@localnames){
			my @splt = split (/\./, $l);
			$l = $splt[0];
		}
		my @splitter = split (/\./, $blf1);
		foreach $ln(@localnames){
			if ($splitter[3] =~ /$ln/){
				my $forward_output = $ln . ".FORWARD.fastq";
				`cat $blf1 >> $forward_output `;
				push @forward_outs, $forward_output;
				last;
			}
			else{
				shift @localnames;
			}
		}
	}

&FAST_MAP();
}

sub FAST_MAP{

	my @localnames = @forward_outs;
	my $finalcounter = 0;

	foreach $l(@localnames){
		my $local_bam_base = (split(/\./, $l))[0];
		my $local_bam_sorted = $local_bam_base . ".sorted";
		my $local_chr = substr($chrom[$finalcounter], 3);
		my $local_script_name = $l . ".bwa_script";
		open (BWA, ">", $local_script_name) or die "Couldn't open bwa script $local_script_name, $!\n";
		print BWA "$config_header";
		print BWA "module load bwa\nmodule load samtools\n";
		print BWA "bwa mem $ref $l $reverse_outs[$finalcounter] | samtools view -Sb -F4 -q $real_samtools_q - | samtools sort -o $local_bam_sorted -";
		push @discard, $local_script_name;
		$finalcounter++;
	}
	print "Assembly completed.\n";
	foreach $disc(@discard){
		`rm $disc`;
	}
	exit;
}
