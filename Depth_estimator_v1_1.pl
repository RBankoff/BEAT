#!/usr/bin/perl -w
use Math::Round;
use Term::ANSIColor;

#input file in .out format
my $file = $ARGV[0];
#size of window to scan through
my $window = $ARGV[1];
#import stats data from file
my $baseline_file = $ARGV[2];
#toggle baseline functionality
my $base = $ARGV[3];
#number of standard deviations acceptable 
my $max = $ARGV[4];


#depth array
my @depth;
my $counter = 0;
unless (defined($base)){
	$base = "n";
}

my @baseline;
my $baseline_a;
my $baseline_n;
my $baseline_s;
if ($base =~ m/y/i){
	open (BASELINE, "<", $baseline_file) or die "Couldn't open baseline file $baseline_file, $!\n";
	while (my $line = <BASELINE>){
		chomp $line;
		push @baseline, $line;
	}
	close BASELINE;

	$baseline_a = $baseline[0];
	$baseline_n = $baseline[1];
	$baseline_s = $baseline[2];
}

my %depth_at_position;
open (FILE, "<", $file) or die "Couldn't open file $file, $!\n";
while (my $line = <FILE>){
	chomp $line;
	my $local_depth = (split(/\t/, $line))[3];
	push @depth, $local_depth;
	$depth_at_position{$counter} = $local_depth;
	$counter++;
}
close FILE;
my $precounter = $counter;

my $running_average = 0;
my $sq_total = 0;
if ($base =~ m/y/i){
	for (my $i = 0; $i < $baseline_n; $i++){
		push @depth, $baseline_a;
		$counter++;
	}
}


my $average;
my $std_dev;

foreach $d(@depth){
	$running_average += $d ; 
}

$average = $running_average / ( $counter + 1 );

foreach $d(@depth){
	$sq_total += (($d - $average) ** 2);
}
$std_dev = ($sq_total / ($counter)) ** 0.5;
if ($base ~= m/n/){
	open (BASE, ">>", $baseline_file) or die "Couldn't open baseline file $baseline_file for output, $!\n";

	print BASE "$file\t$average\t$counter\n";

	close BASE;
}

elsif ($base ~= m/y/i){
	my $local_counter = 0;
	my $block_counter = 1;

	my $condition = int( $precounter / $window );
	#print "Standard_deviation of sample: $std_dev\n";
	#print "Average depth of sample: $average\n";


	$std_dev = $baseline_s;
	$average = $baseline_a;
	print "Standard_deviation of population: $std_dev\n";
	print "Average depth of population: $average\n";

 
	my $max_st_dev = ( $max * $std_dev ) + $average ; 
	my $min_st_dev = $average - ( $max * $std_dev ) ;

	my $prereport = (split(/\//, $file))[0];
	my $report = $prereport;
	$report .= ".report.txt";
	open (OUT_REPORT, ">", $report);
	for (my $i = $condition ; $i > 1 ; $i--){
		my $local_running_average = 0;

		for (my $j = 0 ; $j < $window; $j++){
        		$local_running_average += $depth[$local_counter];
			$local_counter++;
		}
		my $local_average = $local_running_average / $window ;
		my $local_to_five = nearest(5, $local_average);
		my $num_asterisks = $local_to_five / 5 ;
	
		if ($local_average > $max_st_dev){
			print color('bold red on_black');
			print "Block $block_counter (" . ($local_counter - $window ) . "-" . $local_counter . "):" . "*"x$num_asterisks . "\n";
			printf "Average: $local_average";
			print OUT_REPORT "OVER: Block $block_counter (" . ($local_counter - $window ) . "-" . $local_counter . "):" . "*"x$num_asterisks . "\n";
			print color('reset');
			print "\n";
		}
		elsif ($local_average < $min_st_dev){
			print color('bold blue on_white');
              		print "Block $block_counter (" . ($local_counter - $window ) . "-" . $local_counter . "):" . "*"x$num_asterisks . "\n";
                	printf "Average: $local_average";
                	print OUT_REPORT "UNDER: Block $block_counter (" . ($local_counter - $window ) . "-" . $local_counter . "):" . "*"x$num_asterisks . "\n";
                	print color('reset');
			print "\n";
		}

		else{
			print "Block $block_counter (" . ($local_counter - $window ) . "-" . $local_counter . "):" . "*"x$num_asterisks . "\n";
			print "Average: $local_average\n";
		}
	$block_counter++;
	}
}

exit;
