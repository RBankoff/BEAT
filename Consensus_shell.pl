#!/usr/bin/perl -w

#
my $in = $ARGV[0];
my $master_bam = $ARGV[1];
my $basedir_name = $ARGV[2];

open (IN, "<", $in) or die "No such file $in, $!\n";
my @files;
while (my $line = <IN>){
	chomp $line;
	push @files, $line;
}
close IN;
my @refs;
my @chroms;
my @starts;
#my @ends;
my @bases;


foreach $f(@files){
	my @values = split(/\t/, $f);
	push @refs, $values[0];
	push @chroms, $values[1];
	push @starts, $values[2];
	push @ends, $values[3];
	my $base = (split(/\./, $values[0]))[0];
	push @bases, $base;
}

my @mpileups;

my $count = 0;
foreach $c(@chroms){
	my $loc = $c . ":" . $starts[$count] . "-" . $ends[$count];
	push @mpileups, $loc;
	print "MPileup region: $loc\n";
	$count++;
}
foreach $r(@refs){
	print "REFS: $r\n";
}

$count = 0;
foreach $m(@mpileups){
	my $vcf_name = $bases[$count] . "_ALL.mpileup";
	`samtools mpileup -pAf hg38.fa -o $vcf_name -r $m $master_bam`;
	`perl Consensus_mpileup.pl $refs[$count] $chroms[$count] $vcf_name $basedir_name`;
	$count++;

} 

exit;
