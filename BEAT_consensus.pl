#!/usr/bin/perl -w 
use FileHandle;

my $in = $ARGV[0];
my $chr = $ARGV[1];
my $vcf = $ARGV[2];
my $dir_base = $ARGV[3];

$in .= ".gb";
my $base_n = (split(/\./, $in))[0];
my $dir = $base_n . $dir_base;
unless (-e $dir){
	`mkdir $dir`;
}
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
foreach $r (@ranges){
	my @calls;
	$r =~ s/\.\./-/;
	my @pure_range = split ('-', $r);
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
	$r = $chr . ":" . $r;
        my $r_out = $base_n . "_exon_" . $name_count . ".out";
	push @outfiles, $r_out;
	open (MPI, "<", $vcf) or die "couldn't open $vcf, $!\n";
	open (WRITE_ROUT, ">", $r_out) or die "Couldn't open $r_out for output, $!\n";
	while (my $line = <MPI>){
		chomp $line;
		my @mpileup_split = split(/\t/, $line);
		my $local_mpileup_coord = $mpileup_split[1];
		if ($direction =~ /F/){
			if ($local_mpileup_coord >= $pure_range[0]){
				if($local_mpileup_coord <= $pure_range[1]){
					print WRITE_ROUT "$line\n";
				}
			}
		}
		elsif ($direction =~ /R/){
			if ($local_mpileup_coord >= $pure_range[1]){
                                if($local_mpileup_coord	<= $pure_range[0]){
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
			$counter++;
			$deletion_offset--;	
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
		if ($DP > 0){
			$depth_total += $DP;		
			if ($difference == $initiation){
			}
			elsif($difference > $initiation){
				my $local_diff = ($difference - $initiation);
                                $counter += $local_diff;
			}
				
			my @insertions = $alt =~ m/\+((?:[\dACGT])*)/gi;
			my @deletions = $alt =~ m/-((?:[\dACGT])*)/gi;
			my $non_inserts = join('', $alt =~ m/(\.|,|\^|\$)/g);
			my $non_insert_length = length($non_inserts);
			if (scalar(@insertions) > 0){
				my $insertion;
				my $insertion_size = 0;
				my $r = join(' ', @insertions);
                               	my @y = split(/\s/, $r);
                               	$freq_insertions = scalar(@y);
                               	my $actual_insertion;
                               	my $actual_insertion_length;
                               	foreach $l(@y){
                              		my $l1 = $l;
                              		if ($l1 =~ m/^\d(([A-Z]).*)/ig){
                                       		$actual_insertion = $1;
                                       		$actual_insertion_length = $l =~ s/\D//g;
                              		}
                        	}
				$insertion = $actual_insertion;
				$insertion_size = $actual_insertion_length;
				print "FREQ_INSERTION: $freq_insertions\n\n";
				#Modified 11/5/15
				if ($freq_insertions != 0){
					if (($non_insert_length/$freq_insertions) > 1.1){
						my $call = $ref;
						$position_wise[$counter] = $call;
					}
				
					else{ 				
						my $call = $ref . $insertion;
						$position_wise[$counter] = uc($call);
						$deviation += $insertion_size ;
					}
				}			
			}
			elsif(scalar(@deletions) > 0){
				my $deletion;

				foreach $d(@deletions){
					print "Deletion: $d\n";
				}

                                my $deletion_size = 0;
                                my $r = join(' ', @deletions);
                                my @y = split(/\s/, $r);
                                $freq_deletions = scalar(@y);
				if (($freq_deletions) > 0){
                                	my $actual_deletion;
                                	my $actual_deletion_length;
                                	foreach $l(@y){
                      	        		my $l1 = $l;
                                		if ($l1 =~ m/\d(([A-Z]).*)/ig){
                        	        		$actual_deletion = $1;
						}
					}
                                	$deletion = $actual_deletion;
                                	$deletion_size = length($actual_deletion);
					print "DELETION_SIZE: $deletion_size\nDELETION:$deletion\n\n"; 
					print "FREQ_DELETIONS: $freq_deletions\n\n";
                                	if (($non_insert_length/$freq_deletions) < 2){
						my $call = $ref;
                                        	$position_wise[$counter] = $call;
						$deletion_offset += $deletion_size;
						$local_counter = $counter;
						for (my $q = 0 ; $q < $deletion_size ; $q++){
						 	my $D = "N";
							$local_counter++;
							$position_wise[$local_counter] = $D;
						}
                                        	$deviation += $deletion_size;
					}
				}
				else{
					my $call = $ref;
					$position_wise[$counter] = $call;
				}
			}
			else{
				if ($alt !~ /[ACGT]/gi){
			                my $call = $ref;
                                        $position_wise[$counter] = $call;
				}
                		else{
					my @alternate = split(',', $alt);
					my $commas = 0;
					my $periods = 0;
					my $nucleotides = 0;
					my @nucleotides;
                     			foreach $alt_char (@alternate){
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
					my %string = map { $_, 1 } @nucleotides;
					if (keys %string == 1) {
 						my $non_nucleotides = ($periods + $commas);
						if (($non_nucleotides/$nucleotides) > 1){
							my $call = $ref;
			                                $position_wise[$counter] = $call;
						}
						else{
							my $actual_alt = (split('', $alternate[0]))[0];
							$position_wise[$counter] = uc($actual_alt);
							$deviation++;
						}
					}
					else{
						my $call = $ref;
						$position_wise[$counter] = $call;
					}
					#WRITE A HASH FOR ELSE STATMENT
                		}
			}
			$counter++;
		}
		else{
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
exit;

