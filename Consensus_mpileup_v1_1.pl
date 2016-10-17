#!/usr/bin/perl -w 
use FileHandle;
use List::Util qw(max);
use lib '/gpfs/cyberstar/ghp3/Richard/lib/';
use Statistics::Frequency;

my $in = $ARGV[0];
my $chr = $ARGV[1];
my $vcf = $ARGV[2];
my $dir_base = $ARGV[3];

$in =~ s/.fasta$//;
#$in =~ s/.fa$//;
$in .= ".long.gb";
my $base_n = (split(/\./, $in))[0];
my $dir = $base_n . $dir_base;
unless (-e $dir){
	`mkdir $dir`;
}
open (LOG, ">", "log.txt");
print LOG "IN: $in\n";
my @ranges;
my @final;
my @outfiles;
my @Cons_files;
open (IN, "<", $in) or die "no such file $in, $!\n";
while (my $line = <IN>){
	chomp $line;
	my $range = (split(/\t/, $line))[1];
	push @ranges, $range;
}
close IN;
my $name_count = 1;
my $local_option_switch = "y";
my $local_option = 0;

foreach $r (@ranges){
	my $local_option_2 = 0;
	my @calls;
	$r =~ s/\.\./-/;
	my @pure_range = split ('-', $r);
	foreach $pr(@pure_range){
		$pr += 1;
	}
	my $direction_check = $pure_range[1] - $pure_range[0];
	my $direction;
	my @position_wise;

	#direction checking (F/R)
	if ($direction_check < 0){
        	$direction = "R";
		$r = $pure_range[1] . "-" . $pure_range[0];
	}
	else{
        	$direction = "F";
    	}

	my $local_interval;
	if ($local_option_switch =~ m/y/i){
		if ($local_option == 0){
			$local_option += $pure_range[0];
			print LOG "LOCAL OPTION : $local_option\n";
			$local_interval = $local_option;
		}
		else{
			$local_option_2 = $pure_range[0];
			$local_interval = $local_option_2;
		}
	}

	print LOG "local interval: $local_interval\n";
	$r = $chr . ":" . $r;
	my $r_out = $base_n . "_exon_" . $name_count . ".out";
	push @outfiles, $r_out;
	open (MPI, "<", $vcf) or die "couldn't open $vcf, $!\n";
	open (WRITE_ROUT, ">", $r_out) or die "Couldn't open $r_out for output, $!\n";
	while (my $line = <MPI>){
		chomp $line;
#		print LOG "MPI $counter : $line\n";
		my @mpileup_split = split(/\t/, $line);
		my $local_mpileup_coord = $mpileup_split[1];
		my $prelocal_mpileup_coord_base = (split(/:/, $mpileup_split[0]))[1];
		my $local_mpileup_coord_base = (split(/-/, $prelocal_mpileup_coord_base))[0];
		$local_mpileup_coord += $local_mpileup_coord_base;
		print LOG "LOCAL_MPILEUP: $local_mpileup_coord\n";
		if ($direction =~ /F/){
			print LOG "DIR: $direction\n";
			if ($local_mpileup_coord >= $pure_range[0]){
				if($local_mpileup_coord <= $pure_range[1]){
					print LOG "MPI: $line\n";
					print WRITE_ROUT "$line\n";
				}
			}
		}
		elsif ($direction =~ /R/){
			if ($local_mpileup_coord >= $pure_range[1]){
                		if($local_mpileup_coord	<= $pure_range[0]){
					print LOG "MPI: $line\n";
					print WRITE_ROUT "$line\n";
                		}
            		}
		}
	}
	close MPI;
	close WRITE_ROUT;

	open (ROUT, "<", $r_out);
	my @values;
	while (my $line = <ROUT>){
		chomp $line;
		push @values, $line;
	}
	close ROUT;	
	
	$direction_check = abs($direction_check);

	for (my $i = 0; $i < $direction_check; $i++){
		my $x = "N";
		push @position_wise, $x;
	}
	my $counter = 0;
	my $initiation;
	my $increment;
	my $deletion_offset = 0;
	my $deviation = 0;
	my $depth_total = 0;
	foreach $v (@values){
		if ($deletion_offset != 0){
			until ($deletion_offset == 0){
				$counter--;
				$deletion_offset--;	
			}
		}
		elsif ($deletion_offset == 0){
			my @tv = split (/\t/, $v);
			my $pos = $tv[1];
			#print "pos: $pos at $counter\n";
			my $difference;
			if ($counter == 0){
				$initiation = $pos;
				$difference = $initiation;
			}
			else{
				$increment = $pos;
				$difference = ($increment - $counter);
			}
			my $ref = $tv[2];
			my $DP = $tv[3];
			my $alt = $tv[4];
			my $quality = $tv[5];
			if ($DP > 2){
	    	    print LOG "VALUE AT $counter: $v\n";
        	    print LOG "ACTUAL_QUALITY at $counter: $quality\n";

				$depth_total += $DP;		
				if ($difference == $initiation){
				}
				elsif($difference > $initiation){
					my $local_diff = ($difference - $initiation);
        	        $counter += $local_diff;
				}
			
				my @insertions;
				while ($alt =~ m/\+\d((?:[ACGT]+)).*?/gi){
					push @insertions, $1 ;
					print LOG "I1 at counter $counter: $1\n";
				}
		    	my @deletions;
				while ($alt =~ m/\-\d((?:[ACGT\-]+)).*?/gi){
					push @deletions, $1 ;
					print LOG "D1 at counter $counter: $1\n";
				}
		    	$alt =~ s/\-\d((?:[ACGT-]+))//gi;
		    	$alt =~ s/\+\d((?:[ACGT]+))//gi;
				my @split_alts = split('', $alt);
				my $non_inserts = ' ';
				foreach $split_alt (@split_alts){
					$non_inserts .= $split_alt;
				}
				my $non_insert_length = length($non_inserts);
				my @holder = split('', $quality);
				my $alt_positionless = $alt;
                $alt_positionless =~ s/[\$]//g;
				$alt_positionless =~ s/[\^].{1}//g;
				#print LOG "A_P: $alt_positionless\n";
				my @split_alt_2 = split('', $alt_positionless);
			
				print LOG "Num qualities at counter $counter: " . scalar(@holder) . "\n";
            	print LOG "Num alts at counter $counter: " . scalar(@split_alt_2) . "\n";
            	print LOG "ALT_po: $alt_positionless\n";
            	print LOG "QUAL: $quality\n";

				my $qual_counter = 0;
            	my @quality_scores;
            	foreach $spl_alt(@split_alt_2){
                	print LOG "HOLDER[QUAL_COUNTER] at $counter: $holder[$qual_counter] spl_alt: $spl_alt\n";
                	my $ordinal_quality = ord($holder[$qual_counter]);
					print LOG "ORD: $ordinal_quality\n";
                	push @quality_scores, $ordinal_quality;
                	$qual_counter++;
            	}
				#$alt =~ s/[\$]//g;
                #$alt =~ s/[\^].{1}//g;

				if (scalar(@insertions) > 0){
					print "INSERTION at $counter\n";
					my $insertion;
					my $insertion_size = 0;
					my $r = join(' ', @insertions);
                	my @y = split(/\s/, $r);
                	$freq_insertions = scalar(@y);
                	my $actual_insertion;
                	my $actual_insertion_length;
					print LOG "FREQUENCY OF INSERTIONS at counter $counter: $freq_insertions\nFREQUENCY OF NON-INSERTIONS at counter $counter: $non_insert_length\n";
					my @array_of_lengths;
					my @array_of_chars;
		            foreach $l(@y){
                		my $l1 = $l;
						#print LOG "insertion at $counter: $l1\n";
                		$actual_insertion = $l; # =~ s/\d//g;
                		$actual_insertion_length = $l =~ s/\D//g;
						print LOG "insertion at $counter: $actual_insertion with length $actual_insertion_length\n";
						push @array_of_chars, $actual_insertion;
						push @array_of_lengths, $actual_insertion_length;
               		}
					my $ins_freq = Statistics::Frequency->new ( @array_of_chars );
					my %ins_freq = reverse $ins_freq->frequencies;
					my $ins_freq_size = Statistics::Frequency->new ( @array_of_lengths );
		        	my %ins_freq_size = reverse $ins_freq_size->frequencies;
					$insertion = $ins_freq{$ins_freq->frequencies_max};
					$insertion_size = $ins_freq_size{$ins_freq_size->frequencies_max};
		        	print LOG "insertion at $counter is length $insertion_size with characters $insertion\n";
		
					if (($non_insert_length/$freq_insertions) > 1.1){
						if ($non_inserts !~ /[ACGT]/gi){
	               			my $call = $ref;
	               			$position_wise[$counter] = $call;
	          			}
               			else{
            				my @alternate = split('', $alt);
               				my $commas = 0;
               				my $periods = 0;
               				my $nucleotides = 0;
               				my @nucleotides;
               				foreach $alt_char (@alternate){
							#print LOG "ALT_CHAR_NOT_CONS at COUNTER $counter : $alt_char\n";
               					if ($alt_char =~ m/\./){
               						$periods++;
               					}
               					elsif($alt_char =~ m/,/){
               						$commas++;
               					}
               					elsif($alt_char =~ m/[ACGT]/i){
               						$nucleotides++;
               						push @nucleotides, uc($alt_char);
               					}
               				}
							my $nuc_freq = Statistics::Frequency->new ( @nucleotides );
		           			my %nuc_freq = reverse $nuc_freq->frequencies;
							my $highest_nuc;
							my @nuc_qualities;
							my $local_counter_nuc = 0;
							#my %counter_hash;

							foreach $n(@nucleotides){
								print LOG "NUC at counter $counter and localcounter $local_counter_nuc : $n\n";
								my $values = $quality_scores[$local_counter_nuc];
								print LOG "VALUES: $values\n"; 
								push @nuc_qualities, $values;
								#$counter_hash{$n} = $values;#$local_counter_nuc;
								$local_counter_nuc++;
							}

							my @max_quality = @nuc_qualities;
							my $local_max;
							my $local_max_pos = 0;
							my $local_max_count = 0;
							foreach $max(@max_quality){
								if ($local_max_count == 0){
									$local_max = $max;
									$local_max_count++;
								}
								else{
									if ($max > $local_max){
										$local_max = $max;
										$local_max_pos = $local_max_count;
									}
									$local_max_count++;
								}
							}
							my $max_qual = $local_max;
							print LOG "MAX_QUAL: $max_qual\n";
       		        		$highest_nuc = $nuc_freq{$nuc_freq->frequencies_max};
							print LOG "HIGHEST NUC: $highest_nuc\n";
							#my $placeholder = $counter_hash{$max_qual};
							#print LOG "PLACEHOLDER at counter $counter: $placeholder\n";
							my $actual_alt = $nucleotides[$local_max_pos];

							if ($highest_nuc =~ m/$actual_alt/){ 
           						if($nucleotides > 1){
          		  					my $non_nucleotides = ($periods + $commas);
           		  					if (($non_nucleotides/$nucleotides) > 1){
           		      					my $call = $ref;
           								$position_wise[$counter] = $call;
       		     					}
               						else{
               							my $actual_alt = $highest_nuc;
										$position_wise[$counter] = uc($actual_alt);
   		            					$deviation++;
                       				}
                    			}
                    			else{
                       				my $call = $ref;
                       				$position_wise[$counter] = $call;
                       			}
							}
                    		else{
                       			my $call = $ref;
                       			$position_wise[$counter] = $call;
                    		}		
						}
					}
					else{ 	
						print LOG "BOOYEAH, $insertion at $counter!\n";			
						my $call = $ref . $insertion;
						$position_wise[$counter] = uc($call);
						$deviation += $insertion_size ;
					}	
				}
				elsif(scalar(@deletions) > 0){
					print LOG "DELETION at $counter!\n";		
					my $deletion;
                	my $deletion_size = 0;
                	my $r = join(' ', @deletions);
                	my @y = split(/\s/, $r);
					foreach $y1(@y){
						print LOG "FINDME $y1\n";
					}
                	$freq_deletions = scalar(@y);
					my $actual_deletion;
                	my $actual_deletion_length;
					my @actual_deletion;
					my @dels;
                	foreach $l(@y){
                		$actual_deletion = $l; #=~ s/\d//g;
						print LOG "ACTUAL DELETION: $actual_deletion\n";
						$actual_deletion_length = length($actual_deletion);# =~ s/\D//g;
						push @dels, $actual_deletion;
					}
					print LOG "FREQ_DEL at counter $counter: $freq_deletions\n";
                	$deletion = $dels[0];
					print LOG "ACTUAL DELETION at dels[0]: " . $dels[0] . "\n";
                	$deletion_size = length($actual_deletion);
                	if (($non_insert_length/$freq_deletions) > 1.1){
						if ($non_inserts !~ /[ACGT]/gi){
	        		       	my $call = $ref;
	                	   	$position_wise[$counter] = $call;
	             		}
                		else{
             				my @alternate = split('', $alt);
                			my $commas = 0;
                			my $periods = 0;
                			my $nucleotides = 0;
                			my @nucleotides;
                			foreach $alt_char (@alternate){
								#print LOG "ALT_CHAR_NOT_CONS at COUNTER $counter : $alt_char\n";
                				if ($alt_char =~ m/\./){
                					$periods++;
                				}
                				elsif($alt_char =~ m/,/){
                					$commas++;
                				}
                    			elsif($alt_char =~ m/[ACGT]/i){
                    				$nucleotides++;
                       				push @nucleotides, uc($alt_char);
                       			}
                    		}
							my $nuc_freq = Statistics::Frequency->new ( @nucleotides );
		            		my %nuc_freq = reverse $nuc_freq->frequencies;
							my $highest_nuc;
							my @nuc_qualities;
							my $local_counter_nuc = 0;
							#my %counter_hash;

							foreach $n(@nucleotides){
								print LOG "NUC at counter $counter and localcounter $local_counter_nuc : $n\n";
								my $values = $quality_scores[$local_counter_nuc];
								print LOG "VALUES: $values\n"; 
								push @nuc_qualities, $values;
								#$counter_hash{$n} = $values;#$local_counter_nuc;
								$local_counter_nuc++;
							}

							my @max_quality = @nuc_qualities;
							my $local_max;
							my $local_max_pos = 0;
							my $local_max_count = 0;
							foreach $max(@max_quality){
								if ($local_max_count == 0){
									$local_max = $max;
									$local_max_count++;
								}
								else{
									if ($max > $local_max){
										$local_max = $max;
										$local_max_pos = $local_max_count;
									}
									$local_max_count++;
								}
							}
							my $max_qual = $local_max;
							print LOG "MAX_QUAL: $max_qual\n";
       		        		$highest_nuc = $nuc_freq{$nuc_freq->frequencies_max};
							print LOG "HIGHEST NUC: $highest_nuc\n";
							#my $placeholder = $counter_hash{$max_qual};
							#print LOG "PLACEHOLDER at counter $counter: $placeholder\n";
							my $actual_alt = $nucleotides[$local_max_pos];

                    		if ($highest_nuc =~ m/$actual_alt/){ 
           						if($nucleotides > 1){
          		  					my $non_nucleotides = ($periods + $commas);
           		  					if (($non_nucleotides/$nucleotides) > 1){
           		      					my $call = $ref;
           								$position_wise[$counter] = $call;
       		     					}
               						else{
               							my $actual_alt = $highest_nuc;
										$position_wise[$counter] = uc($actual_alt);
   		            					$deviation++;
                       				}
                    			}
                    			else{
                       				my $call = $ref;
                       				$position_wise[$counter] = $call;
                       			}
							}
                    		else{
                       			my $call = $ref;
                       			$position_wise[$counter] = $call;
                    		}		
						}
					}				
					else{
						print "Deletion size: $deletion_size\n";
						$deletion_offset += $deletion_size;
						#$local_counter = $counter;
						my $call = $ref;
	                    $position_wise[$counter] = $call;
						#for (my $q = 0 ; $q < $deletion_size ; $q++){
							#my $D = "-";
							#$counter++;
							#$position_wise[$counter] = $D;
						#}
						#$counter++;
						$deviation += $deletion_size;	
					}
				}
				else{
					if ($non_inserts !~ /[ACGT]/gi){
	        	       	my $call = $ref;
	                   	$position_wise[$counter] = $call;
	             	}
                	else{
             			my @alternate = split('', $alt);
                		my $commas = 0;
                		my $periods = 0;
                		my $nucleotides = 0;
                		my @nucleotides;
                		foreach $alt_char (@alternate){
							print LOG "ALT_CHAR_NOT_CONS at COUNTER $counter : $alt_char\n";
                			if ($alt_char =~ m/\./){
                				$periods++;
                			}
                			elsif($alt_char =~ m/,/){
                				$commas++;
                			}
                    		elsif($alt_char =~ m/[ACGT]/i){
                    			$nucleotides++;
                    	   		push @nucleotides, uc($alt_char);
                    	   	}
                    	}
						my $nuc_freq = Statistics::Frequency->new ( @nucleotides );
		            	my %nuc_freq = reverse $nuc_freq->frequencies;
						my $highest_nuc;
						my @nuc_qualities;
						my $local_counter_nuc = 0;
						#my %counter_hash;

						foreach $n(@nucleotides){
							print LOG "NUC at counter $counter and localcounter $local_counter_nuc : $n\n";
							my $values = $quality_scores[$local_counter_nuc];
							print LOG "VALUES: $values\n"; 
							push @nuc_qualities, $values;
							#$counter_hash{$n} = $values;#$local_counter_nuc;
							$local_counter_nuc++;
						}

						my @max_quality = @nuc_qualities;
						my $local_max;
						my $local_max_pos = 0;
						my $local_max_count = 0;
						foreach $max(@max_quality){
							if ($local_max_count == 0){
								$local_max = $max;
								$local_max_count++;
							}
							else{
								if ($max > $local_max){
									$local_max = $max;
									$local_max_pos = $local_max_count;
								}
								$local_max_count++;
							}
						}
						my $max_qual = $local_max;
						print LOG "MAX_QUAL: $max_qual\n";
       		        	$highest_nuc = $nuc_freq{$nuc_freq->frequencies_max};
						print LOG "HIGHEST NUC: $highest_nuc\n";
						#my $placeholder = $counter_hash{$max_qual};
						#print LOG "PLACEHOLDER at counter $counter: $placeholder\n";
						my $actual_alt = $nucleotides[$local_max_pos];

						if ($highest_nuc =~ m/$actual_alt/){ 
           					if($nucleotides > 1){
          		  				my $non_nucleotides = ($periods + $commas);
           		  				if (($non_nucleotides/$nucleotides) > 1){
           		      				my $call = $ref;
           							$position_wise[$counter] = $call;
       		     				}
               					else{
               						my $actual_alt = $highest_nuc;
									$position_wise[$counter] = uc($actual_alt);
   		            				$deviation++;
                       			}
                    		}
                    		else{
                       			my $call = $ref;
                       			$position_wise[$counter] = $call;
                       		}
						}
                    	else{
                       		my $call = $ref;
                       		$position_wise[$counter] = $call;
                    	}		
					}
				}
				$counter++;
			}
			else{		
				my $call = "N";
	            $position_wise[$counter] = $call;
				$counter++;
			}
		}		
	}

	my $loc_counter = 0;
	my $N_counter = 0;
	foreach $pw (@position_wise){
		print "POSITION_WISE $loc_counter : $pw \n";
		$loc_counter++;
		if ($pw =~ /N/){
			$N_counter++;
		}
	}
	my $final_seq = join ('', @position_wise);
	my $r_cons = $base_n . "_exon_" . $name_count . ".consensus.fa";
	push @Cons_files, $r_cons;
	open (RCONS, ">", $r_cons) or die "no such file as $r_cons, $!\n";
	print RCONS ">$r\n$final_seq\n";
	close RCONS;
	my $corrected_counter = ($counter - $N_counter); 
	my $percent_sim = (1 - ($deviation/$corrected_counter));
	my $final_depth = abs($depth_total/$corrected_counter);
	if ($final_depth != 0){
		print "\tAverage depth of coverage: $final_depth\n";
		print "\tOverall position-wise similarity of query to reference: $percent_sim\n\n";
	}
	else{
		print "\tNo mapping sequence for $r, sorry!\n\n";
	}
	$name_count++;
}
my $outcount = 0;
foreach $c(@Cons_files){
	`mv $c $dir`;
	`mv $outfiles[$outcount] $dir`;
	$outcount++;
}
close LOG;
exit;

