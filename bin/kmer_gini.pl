#!/usr/bin/perl

# kmer_gini.pl
# A plugin script for kmers.pl which calculates the minimum set of kmers required to achieve a threshold value of the Gini evenness index and cover the entire genome.

# Author: MH Seabolt
# Last updated: 5-10-2020

use v5.24.1;
use strict;
use warnings;
use Getopt::Long qw(GetOptions);

use lib '/media/hunter/Data/scripts';
use lib 'D:\scripts';
use lib '/scicomp/home/ngr8/Biolinux/scripts';

# Required input parameters
my $input_sites = "--";
my $output = "--";
my $convergence_threshold;
my $n_sets;
my $k;
my $p;
my $q;
my $v;
my $max_mean_dist;
my $max_gini;
my $evaluate;
my $include;

sub usage {
	my $usage = "kmer_gini.pl\n
	PURPOSE:          	A plugin script for kmers.pl which calculates the minimum set of kmers required to achieve a threshold value of the Gini evenness index and cover the entire genome.

                        \n
	USAGE: kmer_gini.pl [options]
          
	=== INPUT PARAMETERS ===
	-i                     	Input binding sites profile, generated by kmers.pl script.
	-evaluate 		A list of kmers that we explicitly want to evaluate as a set, given either as a CSV list, a file with one kmer per line, or the format from this script.
	-include		A list of mandatory kmers that must be included in set construction, given either as a CSV list or a file with one kmer per line.  The algorithm will build sets around these.

	
	=== SET CONSTRUCTION OPTIONS ===
	-convergence 		FLOAT; Threshold value of the Gini Index to use to estimate when a given set has reached an acceptable convergence.  Must be between 0 and 1.  [ Default: 0.10 ]
	-n_sets			INT; Number of kmer sets to generate (all passing convergence checks).  [Default: 1]
	-buffer_size		INT; Number of previous Gini indices to save for convergence checks.  [ Default: 3 ]						
	-min_set_size		INT; Minimum kmer set size to return [ Default: 2 ]
	-max_set_size 		INT; Maximum kmer set size to return [ Default: no maximum ]
	-max_mean_dist		INT; Maximum mean distance between all kmers in a set [Default: 1,500,000]
	-max_gini		INT; Maximum allowed Gini index for a set of kmers.  Default: 1.00
	
	=== OUTPUT PARAMETERS ===
	-out                    Output file base name, not including extensions.  

							
	=== OPTIONAL EXTRAS ===
	-v 			verbose
	
	\n";
	print $usage;
}
if ( scalar @ARGV == 0 ) { die usage(); }

GetOptions( 'in|i=s' => \$input_sites,
			'convergence=f' => \$convergence_threshold,
			'output|o=s' => \$output,
			'buffer_size|k=i' => \$k,
			'n_sets|n=i' => \$n_sets,
			'min_set_size|p=i' => \$p,
			'max_set_size|q=i' => \$q,
			'max_mean_dist|d=i' => \$max_mean_dist,
			'max_gini=f' => \$max_gini,
			'verbose|v' => \$v,
			'evaluate=s' => \$evaluate,
			'include=s' => \$include,
) or die usage(); 

# Parameter setups
$convergence_threshold = ( $convergence_threshold && $convergence_threshold > 0.0 && $convergence_threshold <= 1.0 )? $convergence_threshold : 0.10;
$n_sets = ( $n_sets && $n_sets >= 1 )? $n_sets : 1;
$max_mean_dist = ( $max_mean_dist && $max_mean_dist > 1 )? $max_mean_dist : 1500000;
$max_gini = ( $max_gini && $max_gini <= 1.0 && $max_gini >= 0.0 )? $max_gini : 1.0;
$k = ( $k && $k >= 1 )? $k : 3;
$p = ( $p && $p >= 1 )? $p : 2;		# The default parameters for $q are below	
$v = defined($v)? 1 : 0;									
my $eval_flag = ( $evaluate )? 1 : 0;	

###########################################################
#       PHASE 0: Preprocess input data
###########################################################v

# Read the binding sites and populate a hash
my $fh = *STDIN;
my $succin = open(SITES, "<", "$input_sites") if ( $input_sites ne "--" && -e $input_sites );
$fh = *SITES if ( $succin );

my %BindingSites = ();
while ( <$fh> )	{
	chomp $_;
	my @line = split("\t", $_);
	my @sites = split(",", $line[1]);
	$BindingSites{ $line[0] } = \@sites;	
}
close SITES if ( $succin );
$q = ( $q && $q <= scalar keys %BindingSites )? $q : scalar keys %BindingSites;		

##############################################
#       PHASE 1: The Kmer-Gini algorithm
##############################################

my $succout = open(OUT, ">", "$output.kmer_sets.gini.tab") if ( $output ne "--" );
if ( $succout ) 	{ 	print OUT    "Set#\tNKmers\tNSites\tMeanDist\tGini\tKmers\tBindingSites\n";		}
else				{	print STDOUT "Set#\tNKmers\tNSites\tMeanDist\tGini\tKmers\tBindingSites\n";		}
open(LOG, ">", "$output.kmer_gini.log") if ( $v == 1 );

# If we are just wanting to evaluate a given set of kmers as a set
if ( $eval_flag == 1 )	{
	# Read the given list
	my @input_data;
	my @records;
	my @names;
	if ( -e $evaluate )		{
		open EVAL, "$evaluate" or die " === KMER GINI Error :: I cannot evaluate anything! URK!\n$!\n";
			@input_data = <EVAL>;
		close EVAL;
	}
	else	{
		@input_data = split(";", $evaluate);
	}
	# Parse the inputs here
	foreach my $line ( @input_data )	{
		chomp $line;
		next if ( $line =~ /^#/ || $line !~ /[\S]/ );
		my @line = split("\t", $line);
		next if ( $line[-2] && $line[-2] =~ /Kmers/ );

		# If we appear to have input from kmer_gini.pl
		if ( scalar @line == 7 ) {
			push @names, $line[0];		# The first field should be the name of the set.
			my @set = split(",", $line[-2]);
			push @records, \@set;
		}
		# If it looks like we are coming from kmer_graph.pl
		elsif ( scalar @line == 3 )	{
			push @names, $line[0];
			my @set = split(",", $line[-1]);
			push @records, \@set;
		}
		# If we just have a simple csv list of kmers with no name
		elsif ( scalar @line == 1 )		{
			my @set = split(",", $line[0]);
			die if ( not exists $BindingSites{$set[0]} );
			push @records, \@set;
			push @names, "Set1";
		}
	}

	# Build out each set we want to evaluate
	for ( my $i=0; $i < scalar @records; $i++ )		{
		my $sublist = $records[$i];
		my $name = $names[$i];
		my %Set = ();
		foreach ( @{$sublist} )	{
			if ( not exists $BindingSites{$_} )	{
				warn " === KMER GINI Warning :: A kmer ( $_ ) in the process of being evaluated isn't present in the binding sites list given! === \n";
				$name = "$name*" if ( $name !~ /\*$/ );
				next;
			}
			$Set{$_} = $BindingSites{$_};
		}
		my @set = keys %Set;
	
		# Run the computation and print a report
		my @sites; 
		@sites = ( @sites, @{$Set{$_}} ) foreach ( keys %Set );
		@sites = sort { $a <=> $b } @sites;
		my @distances = sites2distances( @sites );
		
		my $mean_dist = mean( @distances );
		my $set_gini = gini_index( @distances );
		if ( $succout )		{ 	print OUT    "$name\t", scalar @set, "\t", scalar @sites, "\t$mean_dist\t$set_gini\t", join(",", @set), "\t", join(",", @sites), "\n";		}
		else				{ 	print STDOUT "$name\t", scalar @set, "\t", scalar @sites, "\t$mean_dist\t$set_gini\t", join(",", @set), "\t", join(",", @sites), "\n";		}
	}
}
# Otherwise, we are doing set construction, the typical intended use for this script
else	{
	my %Include = ();
	if ( $include )	{
		if ( -e $include )	{
			open INCLUDE, "$include" or warn " --- KMER GINI Error :: I cant include sites that you didn't give me (as a file)\n$!\n";
				my @must_include = <INCLUDE>;
			close INCLUDE;
			foreach ( @must_include )	{
				chomp $_;
				my @line = split("\t", $_);
				$Include{ $line[0] } = $BindingSites{ $line[0] };
			}
		}
		else	{
			my @must_include = split(",", $include);
			$Include{$_} = $BindingSites{$_} foreach ( @must_include );	
		}
	}
	
	# Generate a new set of suitable primers and add it to the list of found sets
	my @primer_sets;
	while ( scalar @primer_sets < $n_sets )	{
		my ( $mean_dist, $set_gini, $set ) = kmer_gini( \%Include );
		next if ( not $set );
		
		# Separate out the kmer set and all the binding sites
		my @set = keys %{$set};
		push @primer_sets, \@set if ( @set );
		
		# The final Gini index and the binding sites
		my @sites; 
		@sites = ( @sites, @{$_} ) foreach ( values %{$set} );
		
		# Print the output as we go
		if ( $succout )		{ 	print OUT    "$output", scalar @primer_sets, "\t", scalar @set, "\t", scalar @sites, "\t$mean_dist\t$set_gini\t", join(",", @set), "\t", join(",", @sites), "\n";		}
		else				{ 	print STDOUT "$output", scalar @primer_sets, "\t", scalar @set, "\t", scalar @sites, "\t$mean_dist\t$set_gini\t", join(",", @set), "\t", join(",", @sites), "\n";		}
	}
	close OUT if  ( $succout );
}
close LOG if ( $v == 1 );

exit;

######################################## SUBROUTINES ######################################################

# Arithmetic Mean
sub mean	{ 
	return int sum(@_)/@_;
}

sub sum	{
	my ( @numbers ) = @_;
	my $sum = 0;
	foreach ( @numbers )	{	$sum += $_;		}
	return $sum;
}

sub sites2distances 	{
	my ( @sites ) = @_;
	@sites = sort { $a <=> $b } @sites;
	my @distances;
	for ( my $d=1; $d < scalar @sites; $d++ )	{
		my $f = $d - 1;
		push @distances, $sites[$d] - $sites[$f]; 	# The current element minus the previous element
	}
	return @distances;
}


# Calculates the Gini coefficient from a 1D list of non-negative, ascending sorted numbers
# The index ranges from 0 to 1, with values approaching 1 indicating increasing unevenness (we want a value close to 0)
sub gini_index		{
	my ( @distances ) = @_;
	@distances = sort { $a <=> $b } @distances;
	my $n = scalar @distances;
	my $sum = sum( @distances );
	my $gini_index;
	for ( my $i=0; $i < $n; $i++ )	{
		$gini_index += ((2*($i  +1) - $n - 1) * $distances[$i]) / ($n * $sum);		# Function described: https://github.com/oliviaguest/gini
	}
	return sprintf("%1.3f", $gini_index);
}

# Returns TRUE (1) if all the scores are at or below the convergence threshold value 
sub check_convergence	{
	my ( @scores ) = @_;
	return 0 if ( scalar @scores < $k );
	foreach ( @scores )	{
		return 0 if ( $scores[$_] >= $convergence_threshold );
	}
	return 1;
}

# The Kmer-Gini algorithm
sub kmer_gini	{
	my ( $included_sites ) = @_;
	
	# Initialize an empty set S, which will hold the kmers in the set, and an empty array of distances D, which will store all of the combined distances
	my %S = ();  
	my @D = ();
	my @convergence = ();
	my @distances = ();
	my $gini = 1.0;
	my $mean; 
	
	# If we have pre-set sites that must be included:
	if ( $included_sites && scalar keys %{$included_sites} > 0 )	{
		%S = %{$included_sites};
		@D = ( @D, @{$_} ) foreach ( values %S );
		@distances = sites2distances( @D );
		$gini = gini_index( @distances );
		$mean = mean( @distances );
	}

	# Some counters
	my $x = 0;
	my $max = scalar keys %BindingSites;
	print STDERR "$x\t$gini\n" if ( $v == 1 );
	
	# Add the first Gini score to the list of gini scores.
	# This list will act as a buffer for the last k updated 
	# Gini scores and will be used to estimate convergence of the algorithm.
	push @convergence, $gini;

	# Recursion:
	# While we have not achieved convergence of Gini indices, continue to add kmers to the set S until we do.
	while ( check_convergence(@convergence) != 1 )	{
	
		#Choose next kmer to add to %S, check that it isnt already in S
		my $kmer = (keys %BindingSites)[rand keys %BindingSites];
		next if ( exists $S{$kmer} );
		my @E = ( @D, @{$BindingSites{$kmer}} );
		
		# Recalculate the Gini index for the updated @D
		@distances = sites2distances(@E);
		$mean = mean( @distances );
		$gini = gini_index( @distances );
		
		# Skip if we calculate a worse Gini index than previously
		$x++;
		last if ( $x == $max );	
		next if ( $gini > $convergence[-1] );
		
		# Add the chosen kmer to %S, add the distances to @D
		$S{$kmer} = $BindingSites{$kmer};
		delete $BindingSites{$kmer};
		@D = @E;
		print LOG    scalar keys %S, "\t$gini\t$mean\n" if ( $v == 1 );
		print STDERR scalar keys %S, "\t$gini\t$mean\n" if ( $v == 1 );
		
		# Save the updated Gini index to our buffer of Gini indices
		# If saving this score overflows the buffer, then get rid of the oldest score in it.
		push @convergence, $gini;
		shift @convergence if ( scalar @convergence > $k );	
		
		# At this point, the loop will automatically check for convergence before the next iteration.
		last if ( scalar keys %S == $q );
		$x = 0;
	}
	
	# Restore the keys in %S to %BindingSites
	$BindingSites{$_} = $S{$_} foreach ( keys %S );
	
	# Final checks and return
	( $mean <= $max_mean_dist && $gini <= $max_gini )? return ( $mean, $gini, \%S ) : return 0;
}