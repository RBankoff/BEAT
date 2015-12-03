#!/usr/bin/perl -w
my $SRG = $ARGV[0];
my $SR = $SRG;
$SR =~ s/.gz//g;

#`module load blast+`;
if ($SRG =~ m/.gz$/){
	open($fh, sprintf("gunzip -dc %s |", $SRG)) or die "Broken gunzip $!\n";
}
else{
	open($fh, sprintf("cat $SRG|")) or die "NO, $!\n";
}
open ($fh2, "| makeblastdb -in - -title $SR -dbtype nucl -parse_seqids -out $SR") or die "no piping formatdb!, $!\n";

			
#Fastq => Fasta sub
my $localcounter = 0;
while (my $line = <$fh>){
	if ($. % 4==1){
		print $fh2 "\>" . substr($line, 1);
		$localcounter++;
	}
	elsif ($localcounter == 1){
		print $fh2 "$line";
		$localcounter = 0;
	}
	else{
	}
}
close $fh;
close $fh2;
exit;
