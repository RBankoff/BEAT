#!/usr/bin/perl -w
my $in = $ARGV[0];
my $job_name = $ARGV[1];
my $list_name = $ARGV[2];
my $split_value = $ARGV[3];

if (scalar(@ARGV) < 3){
	if ($in =~ m/help/){
		die "Enter the following: perl BEAT_splitter.pl [name of infile containing newline-separated paired-end short reads] [arbitrary name of your specific job] [name of the list of outfiles to be generated]\n[custom size to split files to (optional)]\n";
	}
	else{
		die "Incorrect number of arguments, try typing \"perl BEAT_splitter help\" for more information\n\n";
	}
}

my $default_split_value = 15313576; #Number of lines to reach roughly one gigabyte;
my $real_split_value;
if (scalar(@ARGV) == 3){
	$real_split_value = $default_split_value;
}
elsif(scalar(@ARGV) == 4){
	print "You have entered a custom value $split_value\n";
	$real_split_value = $split_value;
}
my @files;
my @odds;
my @evens;
my @final_odd_names;
my @final_even_names;

open (IN, "<", $in) or die "No such file $in, $!\n";
while (my $line = <IN>){
	chomp $line;
	push @files, $line;
}

push(@evens, $files[$_*2+1])
for 0..int(@files/2)-1;
push(@odds, $files[$_*2])
for 0..int(@files/2)-1;

my $name_counter = 1;
my $counter = 0;
foreach $o(@odds){	
	open (ODD, "<", $o) or die "No such file $o, $!\n";
	my $local_counter = 0;
	while (my $line = <ODD>){	
		if ($local_counter != $real_split_value){
			my $name = $job_name . $name_counter . "_1.fastq";
			open (OUT, ">>", $name) or die "Couldn't make file $name, $!\n";
			print OUT "$line";
			close OUT;
			$local_counter++;
		}
		else{
			$local_counter = 0;
			my $name = $job_name . $name_counter . "_1.fastq";
			push @final_odd_names, $name;
			$name_counter++;
		}
	}

	open (EVEN, "<", $evens[$counter]) or die "No such file $evens[$counter], $!\n";
        $local_counter = 0;
	$name_counter = 1;
        while (my $line = <EVEN>){
                if ($local_counter != $real_split_value){
                        my $name = $job_name . $name_counter . "_2.fastq";
                        open (OUT2, ">>", $name) or die "Couldn't make file $name, $!\n";
                        print OUT2 "$line";
                        close OUT2;
                        $local_counter++;
                }
                else{
                        $local_counter = 0;
                        my $name = $job_name . $name_counter . "_2.fastq";
                        push @final_even_names, $name;
                        $name_counter++;
                }
        }
	$counter++;
}
$name_counter = 0;
open (LIST, ">", $list_name) or die "Couldn't open file $list_name, $!\n";
foreach $l(@final_odd_names){
	print LIST "$l\n";
	print LIST "$final_even_names[$name_counter]\n";
	$name_counter++;
}
close LIST;
exit;
