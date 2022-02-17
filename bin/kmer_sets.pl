#!/usr/bin/perl

# kmer_sets.pl
# A plugin script for kmers.pl which enables set operations to be calculated from an input data matrix.
# Compute a requested Set operation between all pairs of columns in the input matrix.
# Requires CPAN model Set::Scalar for this, we are using these implementations for out set operations API, this script is just a wrapper for that.

# NOTES on input data matrix formatting:
# -- We expect the format to have the kmers (independent variable) as the ROWS of the matrix (the x-axis) and the COLUMNS should be the different source sequences (the dependent, y-axis variable).
# TAB, SSV, CSV input matrices are allowed, specify with -sep option.  Defaults to TAB.

# Author: MH Seabolt
# Last updated: 4-28-2020

use strict;
use warnings;
use Getopt::Long qw(GetOptions);
use Cwd qw(cwd);
use List::Util qw(first shuffle uniq any);
use Scalar::Util qw(looks_like_number);
use File::Basename; 

# Essential CPAN modules
use Set::Scalar; 	

use lib '/media/hunter/Data/scripts';
use lib 'D:\scripts';
use lib '/scicomp/home/ngr8/Biolinux/scripts';

# Required input parameters
my $input = "--";
my $output = "--";
my $output2;
my $sep;
my $header_flag;
my $rownames_flag;
my $operation;
my $NArm;
my $min;
my $v;

sub usage {
	my $usage = "kmer_sets.pl\n
	PURPOSE:           A plugin script for kmers.pl which enables Set operations to be calculated from an input data matrix.
					   Compute a requested Set operation between all pairs of columns in the input matrix.
                        \n
	USAGE: kmer_sets.pl [options]
          
	=== INPUT PARAMETERS ===
	-i                       input data matrix, assumed to be in TAB format unless otherwise specified with -sep.
	-sep                     separator for the input data matrix.  One of 'space, 'csv', 'tab'.
	-header					 INT flag, input data matrix includes a header line (Default: ON).
	-rownames 				 INT flag, input data matrix includes row names as the first element on each row (Default: ON).

	=== SET OPERATION OPTIONS ===
	-operation 				the mathematical Set operation we want to compute. Note that the output of this code will vary slightly depending on the operation requested. 
							
							One of the following options: 
							DERIVATION OPTIONS --> union, intersection, difference, symmetric_difference, unique
							COMPARISON OPTIONS --> equal, disjoint, subset, superset, proper_intersection, proper_subset, proper_superset, compare
	-NArm					 OPTIONAL INT flag; Do you want to remove zeros and NA values in the data? [Default: OFF]	
	
	=== OUTPUT PARAMETERS ===
	-out                     Output file base name, not including extensions.  
							 The main output file is a pairwise symmetrical matrix keyed on the column names.
							 The values are the number of elements in the derived set after comparing two sets, if using a derivation operation.
							 Alternatively, if using a comparison operation, the values are the True/False result of the comparison.
	-output2 				 OPTIONAL; 2nd output file base name, not including extensions. 
							 This file contains the individual data points calculated by the set mathematics, if using a derivation operation
							 (ie. if you use the intersection method, this file contains the specific elements in the derived intersecting set)
							
	=== OPTIONAL EXTRAS ===
	-min 					INT, the minimum value (frequency) for a kmer to be considered part of a set (Default: 1)
	-v 					verbose
	
	\n";
	print $usage;
}
if ( scalar @ARGV == 0 ) { die usage(); }

GetOptions( 'in|i=s' => \$input,
			'sep|s=s' => \$sep,
            'header|h=i' => \$header_flag,
			'rownames|r=i' => \$rownames_flag,
			'operation|x=s' => \$operation,
			'NArm|z=i' => \$NArm,
			'output|o=s' => \$output,
			'output2|out2=s' => \$output2,
			'min=i' => \$min,
			'verbose|v' => \$v,
) or die usage(); 

# Parameter setups
$header_flag = ( $header_flag && $header_flag == 0 )? 0 : 1;
$rownames_flag = ( $rownames_flag && $rownames_flag == 0 )? 0 : 1;
$NArm = ( $NArm && $NArm == 1 )? 1 : 0;
$min = ( $min && looks_like_number($min) )? $min : 1;					# Allowing this value to take many forms, as long as it is a number
$v = defined($v)? 1 : 0;

# Set the expected input data entry separator
if    ( $sep && ($sep =~ /ssv/  || $sep =~ /space/ ) )	{	$sep = " ";		}
elsif ( $sep && ( $sep =~ /csv/ || $sep =~ /comma/ ) )	{	$sep = ",";		}
else 													{	$sep = "\t";	}	

# Sanity check to validate that we have a valid operation to compute
my %Methods = ( 
	'union' => 1, 'difference' => 1, 'intersection' => 1, 'unique' => 1, 'symmetric_difference' => 1, 'complement' => 1, 										# Derivation operations == method 1
	'subset' => 2, 'superset' => 2, 'equal' => 2, 'disjoint' => 2, 'compare' => 3, 'proper_intersection' => 2, 'proper_subset' => 2, 'proper_superset' => 2		# Comparison operations == method 2
);
my $method = ( exists $Methods{$operation} )? $Methods{$operation} : 0;
die " xxx KMER SETS :: Invalid set operation, choose one from the the given list! xxx \n" if ( $method == 0 );

##############################################
#       PHASE 0: Construct the Sets
##############################################

print STDERR " === KMER SETS :: Reading input data === \n" if ( $v == 1 );

# Read the data
my $fh = *STDIN;
my $succin = open(IN, "<", "$input") if ( $input ne "--" && -e $input );
$fh = *IN if ( $succin ); 

print STDERR " === KMER SETS :: Building the Sets === \n" if ( $v == 1 );

# Populate the Sets
my %Sets = ();
my @headers;
while ( <$fh> )	{
	chomp $_;
	my @row = split("$sep", $_);
	
	# Process if we have a HEADER line, tagged with -h option
	if ( $. == 1 )		{
		if ( $header_flag == 1 )	{	 
			$headers[$_]   = $row[$_] foreach ( 0 .. scalar @row - 1 );			
		}
		# If we dont have headers, the first line is data: Initialize the headers to integer values, then operate on the data
		elsif ( $header_flag == 0 )		{	
			if ( $rownames_flag == 1 )	{	$headers[$_-1] = $_ foreach ( 1 .. scalar @row );								}
			else						{	$headers[$_]   = $_ foreach ( 0 .. scalar @row );	my $trash = shift @headers;		}			
		}
		
		@headers = grep { /\S/ } @headers;
		$Sets{$_} = Set::Scalar->new() foreach ( @headers );
		next if ( $header_flag == 1 );
	}
	elsif ( $. == 2 )	{
		# Sanity check -- the header line should have either equal to OR one less element than the subsequence matrix rows
		# Be a little careful here with this, this is a likely place to introduce unwanted bugs --> also, this does not check that every line has the correct number of elements
		die "Kmer_sets.pl ERROR --> Header line contains an incorrect number of elements!\n" unless ( (scalar keys %Sets == (scalar @row - 1)) || (scalar keys %Sets == scalar @row) );
	}	
	
	# Set the row name, either the first element in the row, or just the line number
	my $rowname = ( $rownames_flag == 1 )? shift @row : ( $header_flag == 1 )? $. - 2 : $. - 1;	
	
	# Add this kmer to the correct Sets in the %Sets hash 
	for ( my $i=0; $i < scalar @row; $i++ )		{
		if ( $row[$i] =~ /NA/ )	{
			next if ( $NArm == 1 );
			$Sets{ $headers[$i] }->insert($rowname);
		}
		elsif ( looks_like_number($row[$i]) )	{
			$Sets{ $headers[$i] }->insert($rowname) if ( $row[$i] >= $min );
		}
		else	{
			warn " --- KMER SETS warning :: A value in line $. of your input doesn't look like a number...?!\n";
			next;
		}
	}
}
close IN if ( $input ne "--" );

##############################################
#       PHASE 1: Set Operations
##############################################

print STDERR " === KMER SETS :: Computing the Set mathematics === \n" if ( $v == 1 );

# Computes pairwise set mathematics on each non-redundant pair of columns
# Stores the SIZE of the derived set in a symmetrical NxN matrix 
# AND ALSO stores the actual elements to write out to a file
my $Matrix = [ ];
open OUT2, ">", "$output2.kmer_sets.data" if ( $output2 && $method == 1 );
for ( my $i=0; $i < scalar @headers; $i++ )	{
	for ( my $j=0; $j < scalar @headers; $j++ )	{

		# Initialize an empty set to hold the results
		my $a = $Sets{ $headers[$i] };
		my $b = $Sets{ $headers[$j] };
		my $derived_set = Set::Scalar->new();
		
		# This is just an if-else ladder for the set operations
		# Derivation operations:
		if ( $method == 1 )	{
			if    ( $operation eq "union" )					{	$derived_set = $a->union($b);						}
			elsif ( $operation eq "intersection" )			{	$derived_set = $a->intersection($b);				}
			elsif ( $operation eq "difference" )			{	$derived_set = $a->difference($b);					}
			elsif ( $operation eq "symmetric_difference" )	{	$derived_set = $a->symmetric_difference($b);		}
			elsif ( $operation eq "unique" ) 				{	$derived_set = $a->unique($b);						}
			# We are not including the complement function here since that only operates on a single set, not two as a comparison
			
			# Store the SIZE of the derived set, remember that we are only doing the math for a lower triagle matrix,
			# but we will flip the coordinates to fill the complete square matrix.
			$Matrix->[$i]->[$j] = $derived_set->size;
			
			#Build the list of pairwise differences if the user requested them
			print OUT2 "$headers[$i]\t$headers[$j]\t", join("$sep", $derived_set->members), "\n" if ( $derived_set && $output2 && $method == 1 && $headers[$i] ne $headers[$j] );
		}
		# Comparison operations  --> these methods return true/false, 
		# or compare() returns a string from the following list: "equal", "disjoint", "proper subset", "proper superset", "proper intersect"
		elsif ( $method == 2 )	{
			if    ( $operation eq "equal" ) 				{	$derived_set = $a->is_equal($b);					}
			elsif ( $operation eq "disjoint" ) 				{	$derived_set = $a->is_disjoint($b);					}
			elsif ( $operation eq "proper_intersection" ) 	{	$derived_set = $a->is_properly_intersecting($b);	}
			elsif ( $operation eq "proper_subset" ) 		{	$derived_set = $a->is_proper_subset($b);			}
			elsif ( $operation eq "proper_superset" ) 		{	$derived_set = $a->is_proper_superset($b);			}
			elsif ( $operation eq "subset" ) 				{	$derived_set = $a->is_subset($b);					}
			elsif ( $operation eq "superset" ) 				{	$derived_set = $a->is_superset($b);					}

			# Store the TRUE/FALSE result, remember that we are only doing the math for a lower triagle matrix,
			# but we will flip the coordinates to fill the complete square matrix.
			$Matrix->[$i]->[$j] = ( $derived_set == 1 )? 1 : 0;
		}
		elsif ( $method == 3 )	{
			# This is just the compare method, but it returns string values that are easier to handle in a 3rd case statement
			$derived_set = $a->compare($b);			# $operation eq "compare"
			
			# Store the TRUE/FALSE result, remember that we are only doing the math for a lower triagle matrix,
			# but we will flip the coordinates to fill the complete square matrix.
			$Matrix->[$i]->[$j] = $derived_set;
		}
	}
}
close OUT2 if ( $output2 && $method == 1 );

##############################################
#       PHASE 2: Print the outputs
##############################################

# Print $Matrix, using the same $sep as the inputs and including the headers and rownames
my $succout = open(MATRIX, ">", "$output.kmer_sets.symMatrix") if ( $output ne "--" );
if ( $succout ) 	{ print MATRIX join("$sep", @headers), "\n";	}
else				{ print STDOUT join("$sep", @headers), "\n";	}
for ( my $i=0; $i < scalar @{$Matrix}; $i++ )	{
	if ( $succout )		{ print MATRIX "$headers[$i]$sep", join("$sep", @{$Matrix->[$i]}), "\n";		}
	else				{ print STDOUT "$headers[$i]$sep", join("$sep", @{$Matrix->[$i]}), "\n";		}
}
close MATRIX if ( $succout );

exit;
