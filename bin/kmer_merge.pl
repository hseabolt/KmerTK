#!/usr/bin/perl

# kmer_merge.pl 
# A plugin script for kmers.pl which merges multiple frequency matrices (headers and row names are required for this) into a new, resized larger matrix.
# The intent here is to function a bit like R's cbind or rbind.

# Author: MH Seabolt
# Last updated: 6-3-2020

use strict;
use warnings;
use Getopt::Long qw(GetOptions);
use Cwd qw(cwd);
use List::Util qw(first shuffle uniq any);
use File::Basename; 
use Parallel::ForkManager;

use lib '/media/hunter/Data/scripts';
use lib 'D:\scripts';
use lib '/scicomp/home/ngr8/Biolinux/scripts';

# Required input parameters
my $input_list;
my $input_file;
my $output = "--";
my $v;
my $sep;

sub usage {
	my $usage = "kmer_merge.pl\n
	PURPOSE:           A plugin script for kmers.pl which merges multiple frequency matrices (headers and row names are required for this) into a new, resized larger matrix, in the vein of R's cbind or rbind.
                        \n
	USAGE: kmer_merge.pl [options]
          
	=== INPUT PARAMETERS ===
	-i                       input CSV list of frequency matrices
	-il                      input file containing a list of matrix files, one per line   
	-sep                     separator for the input data matrix.  One of 'space', 'csv', 'tab'.
	
	=== OUTPUT PARAMETERS ===
	-out                     output file base name, not including extensions
	
	=== OPTIONAL EXTRAS ===
	-v 					verbose
	
	\n";
	print $usage;
}
if ( scalar @ARGV == 0 ) { die usage(); }

GetOptions( 'in|i=s' => \$input_list,
			'infile|il=s' => \$input_file,
			'output|o=s' => \$output,	
			'verbose|v' => \$v,		
			'sep|s=s' => \$sep,				
) or die usage(); 

# Parameter setups
$v = defined($v)? 1 : 0;	

# Set the default separator
if    ( $sep && ($sep =~ /ssv/ || $sep =~ /space/) )	{	$sep = " ";		}
elsif ( $sep && ( $sep =~ /csv/ || $sep =~ /comma/ ) )	{	$sep = ",";		}
else 													{	$sep = "\t";	}	

###########################################################################################
# If it appears that we just gave a file containing a list of genomes, one per line
my @proc_order;
if ( $input_file && -e $input_file )	{	
	open INPUT, $input_file or die "ERROR -- Where tf are your genomes?!\n$!\n";
		@proc_order = <INPUT>;
	close INPUT;
}
# If we just passed a CSV list of genomes directly
else		{
	my @tmp = split(",", $input_list);
	push @proc_order, $_ foreach (@tmp);		# Add these to the proc_order list
}
die "ERROR -- Where tf are your genomes?!\n$!\n" if ( scalar @proc_order == 0 );		# We dont have any genomes!

# Initialize the multidimensional hashes and some simple lists
my @names;
foreach my $matrix_file ( @proc_order )	{
	# Get the genome's basename
	chomp $matrix_file;
	my $matrix = fileparse($matrix_file);
	$matrix =~ s/\.kmer.*//g;
	push @names, $matrix;
}
###########################################################################################

# This hash is associated with the _stateToIndex() subroutine, but it must be outside the subroutine scope to be maintainable and accessible.
my $h = 0;				# The initial index in the prebuilt hash of states below, we will increment it if needed
my %Index = ();
my %ReverseIndex = ();

###########################################################################################
# Foreach file in the list of matrices, open it, and add all of it's data
# Note: we are storing ALL kmers, regardless of lexocographic size, etc.  
my %Kmers = ();
my $x = 1;		# Just a counter for $v
# Operate on all the matrix files
foreach my $matrix_file ( @proc_order )	{
	my $matrix = fileparse($matrix_file, qr/\.[^.]*/);
	$matrix =~ s/\.kmer.*//g;
	print STDERR " ---  Processing $matrix ( genome $x / ", scalar @proc_order, " ) --- \n" if ( $v == 1 );
	$x++;
	
	# Open and read the file
	open MATRIX, "$matrix_file" if ( -e $matrix_file );		# Might not have any excluded sequences;
		my @cells = <MATRIX>;
		chomp $cells[0];
	close MATRIX;

	# Remove any random empty elements from the header row and hash the headers into an enumeration
	my @headers = split("$sep", shift @cells);
	@headers = grep { $_ =~ /[\S]/ } @headers;		
	_columnToIndex( $_ ) foreach ( @headers );
	
	foreach my $row ( @cells )	{
		chomp $row;
		my @row = split("$sep", $row);
		my $kmer = shift @row;
		
		#Quick sanity check: Is the length of the row now the same as the number of headers?
		if ( scalar @row != scalar @headers )		{
			warn " === KMER MERGE Warning :: Incorrect number of elements on row $kmer (line ", $. + 1, ") === \n";
			next;
		}
		
		# Resize the row to the correct length
		for ( my $j=0; $j < scalar keys %Index; $j++ )	{
			$Kmers{$kmer}->[$j] = ( exists $Kmers{$kmer}->[$j] )? $Kmers{$kmer}->[$j] : 0;
		}
		
		# Set the data in the appropriate columns
		for ( my $i=0; $i < scalar @row; $i++ )	{	
			my $index = _columnToIndex( $headers[$i] );
			$Kmers{$kmer}->[$index] = ( exists $Kmers{$kmer}->[$index] )? ($Kmers{$kmer}->[$index] + $row[$i]) : $row[$i];
		}
		# There will be alot of missing values in rows that we only see once and various other scenarios. 
		# We will address these when we go to print out the final product
	}
}
###########################################################################################

###########################################################################################
# Print out the newly resized/merged matrix
# Dimensions = ( rows=number of keys in %Kmers, cols=number of keys in %Index )
my @colnames;
for ( my $i=0; $i < scalar keys %Index; $i++ )	{
	push @colnames, _indexToColumn($i);
}

# Print the header line (will be number of Axis 1 elements minus 1)
my $succout = open(OUT, ">", "$output") if ( $output ne '--' );
if ( $succout ) 	{ 	print OUT 	 join("$sep", @colnames), "\n";		}
else 				{ 	print STDOUT join("$sep", @colnames), "\n";		}

# Print each row, adjusting for potential missing values (init anything missing to 0)
foreach my $n ( keys %Kmers )	{
	for ( my $m=0; $m < scalar @colnames; $m++ )	{
		my $index = _columnToIndex( $colnames[$m] );
		$Kmers{$n}->[$m] = ( exists $Kmers{$n}->[$m] )? $Kmers{$n}->[$m] : 0;
	}
	if ( $succout ) 	{	print OUT    "$n$sep", join("$sep", @{$Kmers{$n}}), "\n";	}
	else 				{	print STDOUT "$n$sep", join("$sep", @{$Kmers{$n}}), "\n";	}
}
close OUT if ( $succout );
###########################################################################################

exit;



##########################   SUBROUTINES    ###################################

# Converts key current character into index
sub _columnToIndex	{
	my ( $column )	= @_;
	return if ( not $column );		# Sanity check
	chomp $column;
	my $index;
	
	if ( exists $Index{$column}  ) 	{
		$index = $Index{$column};
	}
	# If the index doesnt exist in the Index hash, then add it at the end and increment the value
	else		{
		$Index{$column} = $h;				
		$ReverseIndex{$h} = $column;
		$h++;
		$index = $Index{$column};
	}
	return $index;
}

# Converts key current index into character
sub _indexToColumn	{
	my ( $index )	= @_;
	my $column;
	if ( exists( $ReverseIndex{$index} ) ) 	{
		$column = $ReverseIndex{$index};
	}
	# If the index doesnt exist in the ReverseIndex hash, then there is nothing we can do...
	return $column;
}