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
while (<IN>) {
	my $line = $_;
	chomp $line;
	if ($line =~ "Sequences producing significant alignments") {
		my $nex = <IN>;
		$nex = <IN>;
		until ($nex =~ /^$/) {
			$nex = <IN>;
			chomp $nex;
			$nex =~ s/^\s+//;
			push @hits, $nex;
		}
	}
}


open IDS, ">$out.IDs.txt";
for my $hit (@hits) {
	chomp $hit;
	unless ($hit =~ /^$/) {
		my $id = (split(/\s{1,}/, $hit))[0];
		$id = $id . " " . (split(/\s{1,}/, $hit))[1];
		my $score = (split(/\s{1,}/, $hit))[3];
		if ($score >= $min) {
			print IDS "$id\n";
		}
	}
}
close IDS;

`perl ./fqselect3.pl $fastq $out.IDs.txt $out.originalReads.fastq `;
`perl ./fqselect3.pl $flip $out.IDs.txt $out.flipReads.fastq`;


