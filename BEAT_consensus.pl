#!/usr/bin/perl -w

my $in = $ARGV[0];
my $master_bam = $ARGV[1];
my $basedir_name = $ARGV[2];
my $ref = $ARGV[3];

if (scalar(@ARGV) < 4){
	die "incorrect input! Program requires the following: [infile (see documentation for formatting)\n";
}
=pod
elsif(scalar(@ARGV) == 4){
	if (($ref !~ m/.fa$/)|($ref !~ m/.fasta$)){
		die "reference not in fasta format, please check file and try again!\n";
	}
	if ($master_bam !~ m/.bam$/){
                die "reference not in fasta format, please check file and try again!\n";
        }
	my $master_index = $master_bam;
	$master_index .= ".bai";
	unless (-e $master_index){
		die ".bam file is not indexed, cannot proceed!\n";
	}
}
=cut
my @files;
my @refs;
my @chroms;
my @starts;
my @ends;
my @bases;
my @mpileups;

open (IN, "<", $in) or die "No such file $in, $!\n";
while (my $line = <IN>){
	chomp $line;
	push @files, $line;
}
close IN;

foreach $f(@files){
	my @values = split(/\t/, $f);
	push @refs, $values[0];
	push @chroms, $values[1];
	push @starts, $values[2];
	push @ends, $values[3];
	my $base = (split(/\./, $values[0]))[0];
	push @bases, $base;
}

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
	`samtools mpileup -pAf $ref -o $vcf_name -r $m $master_bam`;
	`perl Consensus_mpileup.pl $refs[$count] $chroms[$count] $vcf_name $basedir_name`;
	$count++;

} 

exit;
