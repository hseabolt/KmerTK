#!/usr/bin/perl

# kmer_filter.pl
# Plug-in script to read the profiles of kmer-counted genomes and filter out data as needed.
# NOTE: this script does NOT do extensive recalculation of basic properties of kmers, we are only filtering them here

use strict;
use warnings;
use Getopt::Long qw(GetOptions);
use Cwd qw(cwd);
use List::Util qw(first shuffle uniq any);
use Scalar::Util qw(looks_like_number);
use File::Basename;

# Required input parameters
my $input = "--";
my $input_sites;
my $sep;
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
my $exclude;
my $advanced_filters;

sub usage {
	my $usage = "kmer_filter.pl\n
	PURPOSE:            Plug-in script to read the profiles of kmer-counted genomes and combine the frequency counts in a nice way into a single giant matrix.
                        \n
	USAGE: kmer_filter.pl [options]
          
	=== INPUT PARAMETERS ===
	-i                       input kmer profile from kemrs.pl script.
	-b 						 OPTIONAL; input binding sites profile from kmers.pl.  
							 This script will automatically filter that file to contain the same kmers that pass filtering in the main profile.
	-e 						 OPTIONAL; input list of kmers to filter out/exclude explicitly.
	-sep                     separator for the input file(s).  One of 'space', 'csv', 'tab'. ** All input files must use the same separator! **

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
	-advanced 				int; use advanced filtering to catch kmers with runs of length 4 or greater, dimer repeats (4 or more), and triplet repeats (3 or more)? Default: OFF

	=== OPTIONAL EXTRAS ===
	-v 					verbose
	
	\n";
	print $usage;
}
if ( scalar @ARGV == 0 ) { die usage(); }

GetOptions( 'in|i=s' => \$input,
			'sites|b=s' => \$input_sites,
			'exclude|e=s' => \$exclude,
            'mink=i' => \$mink,
            'maxk=i' => \$maxk,
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
		  'advanced=i' => \$advanced_filters,
) or die usage(); 

# Parameter setups
$mink = ( $mink && int($mink) >= 5 )? $mink : 8;
$maxk = ( $maxk && int($maxk) > int($mink) )? $maxk : 'inf';
$min_freq = ( $min_freq && $min_freq >= 0 )? $min_freq : 5;
$max_freq = ( $max_freq && $max_freq >= $min_freq )? $max_freq : 'inf';
$min_tm = ( $min_tm && $min_tm >= 15 )? $min_tm : 15;
$max_tm = ( $max_tm && $max_tm <= 80 )? $max_tm : 45;		
$min_gc = ( $min_gc && $min_gc >= 0.0 )? $min_gc : 0.0;
$max_gc = ( $max_gc && $max_gc <= 1.0 )? $max_gc : 1.0;
$min_bind_dist = ( $min_bind_dist && $min_bind_dist >= 5000 )? $min_bind_dist : 5000;
$max_bind_dist = ( $max_bind_dist && $max_bind_dist <= 1500000 )? $max_bind_dist : 1500000;		# 1,500,000
$allow_ambig = ( $allow_ambig && $allow_ambig <= $maxk )? $allow_ambig : 0;
$max_gini = ( $max_gini && $max_gini <= 1.0 && $max_gini >= 0.0 )? $max_gini : 1.0;
$v = defined($v)? 1 : 0;
$advanced_filters = ( $advanced_filters && $advanced_filters != 0 )? 1 : 0;

# Set the default separator
if    ( $sep && ($sep =~ /ssv/ || $sep =~ /space/) )	{	$sep = " ";		}
elsif ( $sep && ( $sep =~ /csv/ || $sep =~ /comma/ ) )	{	$sep = ",";		}
else 													{	$sep = "\t";	}	

# We are keeping all the parameters from the original kmer counting script so that we can be as strict or exhaustive as we want at that phase,
#  and use this script as an additional way of filtering the data.  
# NOTE: we are not recomputing any calculations here (eg. Gini, Free Energy, etc.), only reading the data file and acceptign or rejecting a line.

my %Exclude = ();
if ( $exclude && -e $exclude )	{
	open EXCLUDE, "$exclude" or warn "ERROR: Cannot find the kmers to exclude -- $!\n";
	while ( <EXCLUDE> )	{
		chomp $_;
		my @line = split("$sep", $_);
		next if ( $line[0] eq 'Kmer' || not $line[0] );		
		$Exclude{ $line[0] } = 1;
	}
	close EXCLUDE;
}

# Read the main data
my $fh = *STDIN;
my $succin = open(IN, "<", "$input") if ( $input ne "--" && -e $input );
$fh = *IN if ( $succin );

# Open the output file and print the headers 
my $succout = open(OUT, ">", "$output.kmers.profiles.filtered.tab") if ( $output ne "--" );

my %PassingKmers = ();
while ( <$fh> )	{
	chomp $_;
	my @kmer_profile = split("$sep", $_);
	my ( $kmer, $count, $length, $gc, $ambigs, $tm, $mean_dist, $gini, $gibbs );
	$kmer = $kmer_profile[0];
	
	# This block just prints the header line
	if ( $. == 1 )	{
		if ( scalar @kmer_profile == 16 )	{
			my @names = qw( Kmer Count Length GC% A T C G N Tm MinDist MaxDist MeanDist StdevDist Gini GibbsFreeEnergy );
			if ( $succout ) 	{ 	print OUT    join("$sep", @names), "\n";	}
			else				{ 	print STDOUT join("$sep", @names), "\n";	}
		}
		elsif ( scalar @kmer_profile == 2 )	{
			my @names = qw( Kmer Count );
			if ( $succout ) 	{ 	print OUT    join("$sep", @names), "\n";	}
			else				{ 	print STDOUT join("$sep", @names), "\n";	}
		}
		elsif ( $input =~ /combined_profiles/ )		{
			if ( $succout ) 	{ 	print OUT    "$_\n";	}
			else				{ 	print STDOUT "$_\n";	}
		}
		next;
	}

	# We probably just have kmers and a count, use a minimal filter here.
	if ( scalar @kmer_profile == 2 )	{
		$count = $kmer_profile[1];
		next if ( exists $Exclude{$kmer} || ($count < $min_freq || $count > $max_freq) );
	}
	# This is where things are expected to strictly adhere to HS kmers.pl profiling format for filtering, but it's entirely possible to just skip these filters entirely I suppose.
	# The fields we actually care about for filtering
	elsif ( scalar @kmer_profile == 16 )	{
		( $count, $length, $gc, $ambigs, $tm, $mean_dist, $gini, $gibbs ) = ( $kmer_profile[1], $kmer_profile[2], $kmer_profile[3], $kmer_profile[8], $kmer_profile[9], $kmer_profile[12], $kmer_profile[14], $kmer_profile[15] );
		
		# Apply our basic filters: the kmer is skipped if it fails ANY of the specified filters
		# --> if you want to apply a ranked filtering process, then run this script on a loop with the given parameters
		next if ( exists $Exclude{$kmer} || ($count < $min_freq || $count > $max_freq) || ($length < $mink || $length > $maxk) || ($gc < $min_gc || $gc > $max_gc) || ($ambigs > $allow_ambig) || ($tm < $min_tm || $tm > $max_tm) || ($mean_dist < $min_bind_dist || $mean_dist > $max_bind_dist) || ($gini > $max_gini) );
	}
	# If we dont recognize the format, then the best we can do is just remove the excluded sequences we were give
	else	{
		next if ( exists $Exclude{ $kmer_profile[0] } );
	}
	
	# If we have some advanced filters toggled on, check for those here.
	# We are looking for runs of the same residue,
	next if ( $advanced_filters == 1 && ( hasRun($kmer) == 1 || hasDimerRepeat($kmer) == 1 || hasTripletRepeat($kmer) == 1 ) );	
	# Not using the subsequence filter currently for algorithm speed, but leaving the code here to build this into the program in the future
	
	# If we make it this far, lets add the kmer to the hash of kmers to filter out of the binding sites data
	$PassingKmers{$kmer} = 1;
	if ( $succout ) 	{ 	print OUT    "$_\n";	}
	else				{ 	print STDOUT "$_\n";	}
}
close OUT if ( $succout );
close IN if ( $succin );

################################################

# If we were given the data for the binding sites, filter that now with the kmers from %PassingKmers
if ( $input_sites && -e $input_sites )	{
	open SITES, "$input_sites" or warn "Cannot open input binding sites data! $!\n";
	open OUT2, ">", "$output.kmers.binding_sites.filtered.tab";
	while ( <SITES> )	{
		chomp $_;
		my @data = split("\t", $_);
		print OUT2 "$_\n" if ( exists $PassingKmers{ $data[0] } );
	}
	close SITES;
	close OUT2;
}

exit;

######################### SUBROUTINES FOR MORE ADVANCED FILTERING #############################

# Does the primer contain runs of any symbol of length 4 (maybe we will want to change the length of the allowable runs later...)
sub hasRun	{
	my ( $kmer ) = @_;
	$kmer =~ m/([ACGTYRWSKMDVHBNX])\1{3,}/;
	return ( $1 )? 1 : 0;
}

# Does the primer contain consecutive dinucleotide repeats?
sub hasDimerRepeat	{
	my ( $kmer ) = @_;
	if ( length $kmer > 8 )	{
		$kmer =~ m/([ACGTYRWSKMDVHBNX][ACGTYRWSKMDVHBNX])\1{2,}/;		# If the primer is shorter than 8, then use 3 consecutive repeats
		return ( $1 )? 1 : 0;
	}
	else	{
		$kmer =~ m/([ACGTYRWSKMDVHBNX][ACGTYRWSKMDVHBNX])\1{3,}/;
		return ( $1 )? 1 : 0;
	}
}

# Does the primer contain consecutive triplet repeats?  Using a threshold of 3 consecutive, repeated triplets
sub hasTripletRepeat	{
	my ( $kmer ) = @_;
	$kmer =~ m/([ACGTYRWSKMDVHBNX][ACGTYRWSKMDVHBNX][ACGTYRWSKMDVHBNX])\1{2,}/;
	return ( $1 )? 1 : 0;
}

# Returns true if one string is a subsequence of the other, in a variety of comparisons
sub isSubsequence	{
	my ( $string1, $string2 ) = @_;
	return 1 if ( $string1 =~ m/$string2/ || $string2 =~ m/$string1/ );						# The original strings
	return 1 if ( revcom($string1) =~ m/$string2/ || revcom($string2) =~ m/$string1/ );		# The reverse complements
	return 0;
}


