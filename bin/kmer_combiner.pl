#!/usr/bin/perl

# kmer_combiner.pl
# Plug-in script to read the profiles of kmer-counted genomes and combine them in a nice way into a single giant matrix.

use v5.24.1;
use strict;
use warnings;
use Getopt::Long qw(GetOptions);
use Cwd qw(cwd);
use List::Util qw(first shuffle uniq any);
use Scalar::Util qw(looks_like_number);
use File::Basename;

# Required input parameters
my $input_list;
my $genome_list;
my $sep;
my $header_flag;
my $rownames_flag;
my $column;
my $output = "--";
my $mink;
my $maxk;
my $min_freq;
my $max_freq;
my $min_tm;
my $max_tm;
my $min_gc;
my $max_gc;
my $max_gini;
my $min_bind_dist;
my $max_bind_dist;
my $allow_ambig;
my $v;

sub usage {
	my $usage = "kmer_combiner.pl\n
	PURPOSE:            Plug-in script to read the profiles of kmer-counted genomes and combine the frequency counts in a nice way into a single giant matrix.
                        \n
	USAGE: kmer_combiner.pl [options]
          
	=== INPUT PARAMETERS ===
	-i                       input CSV list of FASTA formatted genome files
	-il                      input file containing a list of genome files, one per line
	-sep                     separator for the input data matrices.  One of 'space', 'csv', 'tab'. 
							 **** NOTE: ALL input matrices MUST adhere to the same format, including same column order. ****
	-header					 INT flag, input data matrices include header lines (Default: ON).
	-rownames 				 INT flag, input data matrices include row names as the first element on each row (Default: ON).
	-combine				 Column name or index that we want to combine from each profile. 
							(Uses 0-indexing, respectful of row names --> ie. the first column of data is index 0 if using row names)

	=== OUTPUT PARAMETERS ===
	-out                     output file base name, not including extensions

	=== FILTERING PARAMETERS ===
	-mink                    INT; the minimum value of k [Default: 8]
	-maxk                    INT; the maximum value of k [Default: 12]
	-min_freq			int; discard a primer if it does not bind to the foreground genome with at least this frequency.  Default: 5
	-max_freq			int; discard a kmer if it appears more often than this.  Default: n/a
	-min_tm				int; discard a primer if the estimated melting temperature is below this threshold value.  Default: 15
	-max_tm				int; discard a primer if the estimated melting temperature is above this threshold value.  Default: 45
	-min_gc				float; discard a primer if the GC% is below this threshold value.  Default: 15
	-max_gc				float; discard a primer if the GC% is above this threshold value.  Default: 45
	-max_gini				float; discard a primer if the individual-based Gini index is above this threshold value.  Default: 0.75
	-allow_ambig			int; allow primer sequences to contain up to X ambiguous bases. Default: 0
	-min_bind_dist			int; discard a kmer if the mean binding distance is less than this threshold value.  Default: 5000
	-max_bind_dist			int; discard a kmer if the mean binding distance is greater than this threshold value.  Default: 1500000

	=== OPTIONAL EXTRAS ===
	-v 					verbose
	
	\n";
	print $usage;
}
if ( scalar @ARGV == 0 ) { die usage(); }

GetOptions( 'in|i=s' => \$genome_list,
			'inlist|il=s' => \$input_list,
            'mink=i' => \$mink,
            'maxk=i' => \$maxk,
			'combine=s' => \$column,
			'header|h=i' => \$header_flag,
			'rownames|r=i' => \$rownames_flag,
			'min_freq=i' => \$min_freq,
			'max_freq=i' => \$max_freq,
            'verbose|v' => \$v,
			'output|o=s' => \$output,
			'min_gc=f' => \$min_gc,
			'max_gc=f' => \$max_gc,
            'min_tm=i' => \$min_tm,
            'max_tm=i' => \$max_tm,
            'max_gini=f' => \$max_gini,
		  'min_bind_dist=i' => \$min_bind_dist,
		  'max_bind_dist=i' => \$max_bind_dist,
		  'allow_ambig=i' => \$allow_ambig,	
) or die usage(); 

# Parameter setups
$column = ( $column )? $column : "Count";
$mink = ( $mink && int($mink) >= 5 )? $mink : 8;
$maxk = ( $maxk && int($maxk) > int($mink) )? $maxk : 'inf';
$min_freq = ( $min_freq && $min_freq >= 1 )? $min_freq : 5;
$max_freq = ( $max_freq && $max_freq > $min_freq )? $max_freq : 'inf';
$min_tm = ( $min_tm && $min_tm >= 15 )? $min_tm : 15;
$max_tm = ( $max_tm && $max_tm <= 80 )? $max_tm : 45;		
$min_gc = ( $min_gc && $min_gc >= 0.0 )? $min_gc : 0.0;
$max_gc = ( $max_gc && $max_gc <= 1.0 )? $max_gc : 1.0;
$min_bind_dist = ( $min_bind_dist && $min_bind_dist >= 5000 )? $min_bind_dist : 5000;
$max_bind_dist = ( $max_bind_dist && $max_bind_dist <= 1500000 )? $max_bind_dist : 1500000;		# 1,500,000
$allow_ambig = ( $allow_ambig && $allow_ambig <= $maxk )? $allow_ambig : 0;
$max_gini = ( $max_gini && $max_gini <= 1.0 && $max_gini >= 0.0 )? $max_gini : 1.0;
$header_flag = ( $header_flag && $header_flag == 0 )? 0 : 1;
$rownames_flag = ( $rownames_flag && $rownames_flag == 0 )? 0 : 1;
$v = defined($v)? 1 : 0;

# Set the default separator
if    ( $sep && ($sep =~ /ssv/ || $sep =~ /space/) )	{	$sep = " ";		}
elsif ( $sep && ( $sep =~ /csv/ || $sep =~ /comma/ ) )	{	$sep = ",";		}
else 													{	$sep = "\t";	}	

# We are keeping all the parameters from the original kmer counting script so that we can be as strict or exhaustive as we want at that phase,
#  and use this script as an additional way of filtering the data.  
# NOTE: we are not recomputing any calculations here (eg. Gini, Free Energy, etc.), only reading the data file and acceptign or rejecting a line.

###########################################################
#       PHASE 0: Preprocessing
###########################################################

# If it appears that we just gave a file containing a list of genomes, one per line
my @proc_order;
if ( $input_list )	{	
	open INPUT, $input_list or die "ERROR -- Where tf are your genomes?!\n$!\n";
		@proc_order = <INPUT>;
	close INPUT;
}
# If we just passed a CSV list of genomes directly
if ( $genome_list )	{
	my @tmp = split(",", $genome_list);
	push @proc_order, $_ foreach (@tmp);		# Add these to the proc_order list
}
die "ERROR -- Where tf are your genomes?!\n$!\n" if ( scalar @proc_order == 0 );		# We dont have any genomes!

# Adding this block just to make the output prettier, but it's not really necessary
my @names;
foreach my $genome_file ( @proc_order )	{
	my $genome = fileparse($genome_file, qr/\.[^.]*/);
	$genome =~ s/\.kmer.*//g;
	push @names, $genome;
}
#@names = grep { $_ =~ /[\S]/ } @names;	

# Remove any random empty elements and hash the headers into an enumeration
my $h = 0;
my %Index = ();
my %ReverseIndex = ();
my @headers = qw(Kmer Count Length GC% A T C G N Tm MinDist MaxDist MeanDist StdevDist Gini GibbsFreeEnergy);
_columnToIndex( $_ ) foreach ( @headers );

# Verify that the requested column exists in the data provided
my $index = ( looks_like_number($column) )? $column : _columnToIndex($column);
die " === KMER COMBINER Error :: Invalid data request! I can't combine a column that isnt in the input data! === \n" if ( not $headers[ $index ] );

###########################################################
#       PHASE 1: Operate on the real data
###########################################################

# Read in all the data --> WARNING: This might be very resource intensive depending on how many genomes you have and how loose your filtering parameters are
my %Kmers = ();
my $x = 1;		# Just a counter for $v
# Operate on all the genomes for each value of k
foreach my $genome_file ( @proc_order )	{
	
	my $genome = fileparse($genome_file, qr/\.[^.]*/);
	$genome =~ s/\.kmer.*//g;
	print STDERR " ---  Processing $genome ( genome $x / ", scalar @proc_order, " ) --- \n" if ( $v == 1 );
	$x++;
	
	# Open and read the file
	open IN, "$genome_file" or warn "I cant open $genome_file -- $!\n";
		my @genome = <IN>;
	close IN;
	next if ( scalar @genome == 1 );
	my $trash = shift @genome;  # Throw away the header line
	
	# Parse the input data
	foreach my $line ( @genome )	{
		chomp $line;
		my @kmer_profile = split("$sep", $line);
		my ( $kmer, $count, $length, $gc, $ambigs, $tm, $mean_dist, $gini, $gibbs );
		
		# We probably just have kmers and a count, use a minimal filter here.
		if ( scalar @kmer_profile == 2 )	{
			( $kmer, $count ) = ( $kmer_profile[0], $kmer_profile[1] );
			next if ( $count < $min_freq || $count > $max_freq );
		}
		# This is where things are expected to strictly adhere to HS kmers.pl profiling format for filtering, but it's entirely possible to just skip these filters entirely I suppose.
		# The fields we actually care about for filtering
		elsif ( scalar @kmer_profile == 16 )	{
			( $kmer, $count, $length, $gc, $ambigs, $tm, $mean_dist, $gini, $gibbs ) = ( $kmer_profile[0], $kmer_profile[1], $kmer_profile[2], $kmer_profile[3], $kmer_profile[8], $kmer_profile[9], $kmer_profile[12], $kmer_profile[14], $kmer_profile[15] );
			
			# Apply our basic filters: the kmer is skipped if it fails ANY of the specified filters
			# --> if you want to apply a ranked filtering process, then run this script on a loop with the given parameters
			next if ( ($count < $min_freq || $count > $max_freq) || ($length < $mink || $length > $maxk) || ($gc < $min_gc || $gc > $max_gc) || ($ambigs > $allow_ambig) || ($tm < $min_tm || $tm > $max_tm) || ($mean_dist < $min_bind_dist || $mean_dist > $max_bind_dist) || ($gini > $max_gini) );
		}
		
		# If we make it this far, lets add the kmer to the big hash if we've never seen it before
		if ( exists $Kmers{$kmer} )		{
			$Kmers{$kmer}->{$genome} = $kmer_profile[$index];
		}
		else	{
			$Kmers{$kmer}->{$_} = 0 foreach ( @names );	# Set the count to 0 for all the genomes, then update the one we've found this kmer in.
			$Kmers{$kmer}->{$genome} = $kmer_profile[$index];
		}
	}
}

################################################
# Finally, let's print our output matrix

# Print the headers 
my $succout = open(OUT, ">", "$output.kmers.combined_profiles.tab") if ( $output ne "--" );
if ( $succout ) 	{ 	print OUT    join("$sep", @names), "\n";	}
else				{ 	print STDOUT join("$sep", @names), "\n";	}
foreach my $kmer ( keys %Kmers )	{
	my @data;
	foreach my $genome ( @names )	{
		push @data, $Kmers{$kmer}->{$genome};
	}
	if ( $succout ) 	{	print OUT    "$kmer$sep", join("$sep", @data), "\n";		}
	else				{	print STDOUT "$kmer$sep", join("$sep", @data), "\n";		}
}
close OUT if ( $succout );

exit;


################### SUBROUTINES ##############################

sub _columnToIndex		{
	my ( $header )	= @_;
	return if ( not $header );		# Sanity check
	my $index;
	
	if ( exists $Index{$header}  ) 	{
		$index = $Index{$header};
	}
	# If the index doesnt exist in the Index hash, then add it at the end and increment the value
	else		{
		$Index{$header} = $h;				
		$ReverseIndex{$h} = $header;
		$h++;
		$index = $Index{$header};
	}
	return $index;
}

# Converts key current index into character
sub _indexToColumn	{
	my ( $index )	= @_;
	my $header;
	if ( exists( $ReverseIndex{$index} ) ) 	{
		$header = $ReverseIndex{$index};
	}
	# If the index doesnt exist in the ReverseIndex hash, then there is nothing we can do...
	return $header;
}