#!/usr/bin/perl
use warnings;

my $name = $ARGV[0];
my $taxon = $ARGV[1];

my @coords;
my @starts;
my @ends;
my @precoords = qx(esearch -db gene -query "$name [GENE] $taxon [ORGN] Exon [WORD]" |  efetch -format gene_table);

my $next = "mRNA";
my $single = 0;
my $Chrom;
my $num_exons = 0;
foreach $f (@precoords){
	if ($single < 2){
		chomp $f;
		if ($f =~ /\d{1,}-\d{1,}\t/){
			push @coords, $f;
			$num_exons++;
		}
		elsif ($f =~ /NC_\d{0,}\.\d{0,}\s/){
			my $preChrom = $f;
			my @splitter = split(/\s/, $preChrom);
			foreach $s(@splitter){
				if ($s =~ m/^NC/){
					my @splitter_2 = split(/\./, $s);
					my $local_Chrom = substr($splitter_2[0], (length($splitter_2[0])-2));
					#print "LOCAL_CHROM: $local_Chrom\n";
					if((int($local_Chrom)) < 23){
						my @chars = split('', $local_Chrom);
						if ($chars[0] =~ /0/){
							#print "CHARS_0: $chars[0]\nCHARS_1: $chars[1]\n";
							$Chrom = $chars[1];
						}
						else{
							#print "CHARS_0: $chars[0]\nCHARS_1: $chars[1]\n";
                                                        $Chrom = $chars[0] . $chars[1];
						}
					}
				}
			} 
		}

		elsif($f =~ /^$next/){
			my $field = 'FF';
           		push @coords, $field;
			$num_exons = 0;
        	}
	}
	else{
		last;
	}

}
my @pre_concat_coords;
foreach $coo (@coords){
	my $pre_concat_coord = (split (/\t/, $coo))[0];
	push @pre_concat_coords, $pre_concat_coord;
}

my $concat_coords = join(',', @pre_concat_coords);
my @sizes = split (/FF/, $concat_coords);
my @sorted = sort { length $a <=> length $b } @sizes;
my @true_coords = split(/,/, $sorted[(scalar(@sorted)-1)]);
shift @true_coords;

foreach $cr (@true_coords){
	my @c_plus = split(/-/, $cr);
	my $start = $c_plus[0];
	my $end = $c_plus[1];
	push @starts, $start;
	push @ends, $end;
}

my $direction_check = ($starts[0]) - ($ends[0]);
my $direction;
if ($direction_check < 0){
	$direction = "F";
}
if ($direction_check > 0){
	$direction = "R";
	@starts = reverse @starts;
	@ends = reverse @ends;
}

my $genbank = $name . ".long.gb";
open (GENBANK, ">", $genbank);
my $last_counter = 0;
foreach $x(@starts){
	print GENBANK " exon\t$x" . ".." . "$ends[$last_counter]\n";
	$last_counter++;
}
close GENBANK;

print "$name\tchr$Chrom\t$starts[0]\t$ends[(scalar(@ends)-1)]\n";

exit;

