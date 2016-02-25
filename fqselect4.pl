#!/usr/bin/perl

use strict;
use warnings;

my $in = $ARGV[0];	#whole fastq file
my $keep = $ARGV[1];	#list of IDs to keep (EXACT MATCHES!!)
my $out = $ARGV[2];	#outfile, new fastq

unless (@ARGV==3) {
	die "\nRequired: file.fastq (no weird wrapping please), list of IDs to keep, outfile\n\n";
}

unless (-f $in) {
	die "\nCouldn't find $in\n\n";
}

unless (-f $keep) {
	die "\nCouldn't find $keep\n\n";
}

my %keep;
open (KEEP, "<", "$keep") or die "couldn't open file $keep, $!\n";
while (my $id = <KEEP>) {
	chomp $id;
	print "FQS4ID: $id\n";
	#my @ids = split (/:/, $id);
	my $c = $id;
	$c =~ s/\s+//g;
	my $nums = "$c";
	$keep{$nums} = 0;
}
close KEEP;

my $gz = 0;
if ($in =~ /.gz$/) {
	$gz = 1;
}

if ($gz==1) {
        open(IN, sprintf("zcat %s |", $in));
}

else {
	open IN, $in;
}

my $NR;
open OUT, ">$out";
while (<IN>) {
	$NR++;
	my $line = $_;
	chomp $line;
	if ($NR%4==1) {
		if ($line !~ /^\@/) {
			die "\nERROR!! Your fastq file is truncated, squelched, or wrapped at line $NR, please fix it and try again...\n\n";
		}
		my $header = substr $line, 1;
	        my @ids = split(/\#/, $header);
		my $c = $ids[0];
		#print "C: $c\n";
	        #$c =~ s/\s//g;
	        my $nums = "$c";
		if (exists ($keep{$nums})) {
			print OUT "\@$header\n";
			my $i = 0;
			until ($i==3) {
				my $rec = <IN>;
				chomp $rec;
				print OUT "$rec\n";
				$i++;
				$NR++;
			}
		}
	}
}
