#!/usr/bin/perl

use strict;
use warnings;

#takes a blast file as input, and creates a fq read file with the IDs having >= given score.
#
my $blast = $ARGV[0];
my $min = $ARGV[1];
my $fastq = $ARGV[2];
my $flip = $ARGV[3];
my $out = $ARGV[4];

unless (@ARGV==5) {
	die "\nRequired: file.blast, minimum score, reads.fastq (.gz allowed), otherReads.fastq (.gz allowed), out prefix\n\n";
}

my $check = 0;
my @hits;
open IN, "$blast";
my $switch = 0;
while (my $line = <IN>) {
	chomp $line;
	if ($line =~ "Sequences producing significant alignments") {
		$switch = 1;
	}
	if ($switch == 1){
		if($line !~ m/^>/){
			if ($line !~ "Sequences producing significant alignments"){
				print "LINE: $line\n";
				push @hits, $line;
			}
		}
		else{
			$switch = 0;
		}
	}		
}


open IDS, ">$out.IDs.txt";
for my $hit (@hits) {
	chomp $hit;
	unless ($hit =~ /^$/) {
		print "HIT: $hit\n";
		my $id = (split(/\s{1,}/, $hit))[1];
		my $id2 = (split(/\#/, $id))[0] ;
		print "ID: $id2\n"; 
		my $score = (split(/\s{1,}/, $hit))[2];
		print "SCORE: $score\n";
		if ($score >= $min) {
			print IDS "$id2\n";
		}
	}
}
close IDS;

`perl ./fqselect4.pl $fastq $out.IDs.txt $out.originalReads.fastq `;
`perl ./fqselect4.pl $flip $out.IDs.txt $out.flipReads.fastq`;


