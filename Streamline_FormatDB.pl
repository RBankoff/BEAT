#!/usr/bin/perl -w
my $SRG = $ARGV[0];
my $SR = $SRG;

if ($SRG =~ m/.gz/){
	$SR =~ s/.gz//g;
	open($fh, sprintf("gunzip -dc %s |", $SRG)) or die "Broken gunzip $!\n";
	open ($fh2, "| makeblastdb -in - -title $SR -dbtype nucl -out $SR") or die "no piping formatdb!, $!\n";
				
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

}
else{
	$SR =~ s/.fastq//g;
	open($fh, sprintf("cat |", $SRG)) or die "Broken cat $!\n";
	open ($fh2, "| makeblastdb -in - -title $SR -dbtype nucl -out $SR") or die "no piping formatdb!, $!\n";
			
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
}
