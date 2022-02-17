#!/usr/bin/perl

# kmer_unique.pl
# A plugin script for kmers.pl to read a matrix of the frequency counts and determine if a kmer is unique to a single genome.


# Author: MH Seabolt
# Last updated: 4-21-2020

use strict;
use warnings;
use Getopt::Long qw(GetOptions);
use Cwd qw(cwd);
use List::Util qw(first shuffle uniq any);
use Scalar::Util;
use File::Basename; 
use Parallel::ForkManager;

use lib '/media/hunter/Data/scripts';
use lib 'D:\scripts';
use lib '/scicomp/home/ngr8/Biolinux/scripts';

# Required input parameters
my $input = "--";
my $output = "--";
my $sep;
my $header_flag;
my $rownames_flag;
my $maxg;
my $v;
my $threads;

sub usage {
	my $usage = "kmer_unique.pl\n
	PURPOSE:           A plugin script for kmers.pl to read a matrix of the frequency counts and determine if a kmer is unique to a single genome.
                        \n
	USAGE: kmer_unique.pl [options]
          
	=== INPUT PARAMETERS ===
	-i                       input data matrix, assumed to be in TAB format unless otherwise specified with -sep.
	-sep                     separator for the input data matrix.  One of 'space;, 'csv', 'tab'.
	-header					 INT flag, input data matrix includes a header line (Default: ON).
	-rownames 				 INT flag, input data matrix includes row names as the first element on each row (Default: ON).

	=== OUTPUT PARAMETERS ===
	-out                     output file base name, not including extensions
	
	=== UNIQUENESS PARAMETERS ===
	-maxg                    INT; the maximum number of genomes that a given kmer can appear in to be considered 'unique'. [Default: 1]
	
	=== OPTIONAL EXTRAS ===
	-v 					verbose
	-t 					INT; number of threads to use for multi-threading
	
	\n";
	print $usage;
}
if ( scalar @ARGV == 0 ) { die usage(); }

GetOptions( 'in|i=s' => \$input,
			'sep|s=s' => \$sep,
            'header|h=i' => \$header_flag,
			'rownames|r=i' => \$rownames_flag,
            'maxg=i' => \$maxg,
			'output|o=s' => \$output,
			'verbose|v' => \$v,
		  'threads|t=i' => \$threads,
) or die usage(); 

# Parameter setups
$maxg = ( $maxg && int($maxg) >= 1 )? $maxg : 1;
$header_flag = ( $header_flag && $header_flag == 0 )? 0 : 1;
$rownames_flag = ( $rownames_flag && $rownames_flag == 0 )? 0 : 1;
$v = defined($v)? 1 : 0;
$threads = ( $threads && $threads >= 1 )? $threads : 1;

# Set the default separator
if    ( $sep && ($sep =~ /ssv/ || $sep =~ /space/) )	{	$sep = " ";		}
elsif ( $sep && ( $sep =~ /csv/ || $sep =~ /comma/ ) )	{	$sep = ",";		}
else 													{	$sep = "\t";	}	

# Preprocess the first line of data
my @headers;
my $line1;
if ( $input ne "--" && -e $input )	{
	open IN, "$input" or die "Cannot open the input file given -- $!\n";
		$line1 = <IN>;	
	close IN;
}
else	{
	$line1 = <STDIN>;
}	
chomp $line1;
my @line1 = split("$sep", $line1);
if ( $header_flag == 1 )	{	 
	$headers[$_]   = $line1[$_] foreach ( 0 .. scalar @line1 - 1 );			
}
# If we dont have headers, the first line is data: Initialize the headers to integer values, then operate on the data
elsif ( $header_flag == 0 )		{	
	if ( $rownames_flag == 1 )	{	$headers[$_-1] = $_ foreach ( 1 .. scalar @line1 );								}
	else						{	$headers[$_]   = $_ foreach ( 0 .. scalar @line1 );	my $trash = shift @headers;		}			
}
# Remove any random empty elements and hash the headers into an enumeration
my $h = 0;
my %Index = ();
my %ReverseIndex = ();
@headers = grep { $_ =~ /\S/ } @headers;
_columnToIndex( $_ ) foreach ( @headers );


###########################################################
#       PHASE 1: Read and operate on the real data
###########################################################


# Read the data
my $fh = *STDIN;
my $succin = open(IN, "<", "$input") if ( $input ne "--" && -e $input );
$fh = *IN if ( $succin ); 

# Instantiate out Fork Manager and open our output file 
my $manager = Parallel::ForkManager->new($threads);
my $succout = open(OUT, ">", "$output.kmers.unique.tab") if ( $output ne "--" );

if ( $header_flag == 1 )	{
	if ( $succout )		{	print OUT    "Count$sep","Appears_in\n";		}	
	else				{ 	print STDOUT "Count$sep","Appears_in\n";		}
}

# Using a while loop since I assume that the input file will be too large to read into RAM to work on
my $x = 0;
while ( <$fh> )	{
	$manager->start and next;
	chomp $_;
	my @row = split("$sep", $_);

	# Set the row name, either the first element in the row, or just the line number
	my $rowname = ( $rownames_flag == 1 )? shift @row : ( $header_flag == 1 )? $. - 2 : $. - 1;		
	
	# Count the non-zero elements in the row
	my ( $count, $appears_in ) = count(@row);
	
	# Print the line if the kmer is considered DISTINCT ( appears in less than the maximum number of GENOMES)
	if ( $succout ) 	{	print OUT    "$rowname$sep$count$sep", join(",", @{$appears_in}), "\n" if ( $count <= $maxg );		}
	else				{ 	print STDOUT "$rowname$sep$count$sep", join(",", @{$appears_in}), "\n" if ( $count <= $maxg );		}

	# Terminate the child process and increment x
	$x++;
	$manager->finish;
}
close IN if ( $input ne "--" );
close OUT if ( $succout );

exit;


################## SUBROUTINES ##################################

# Returns the number of non-empty elements in a given array
sub count	{
	my ( @data ) = ( scalar @_ == 0 )? return 0 : @_;
	my $appears_in = [ ];
	
	my @indices;
	for ( my $i=0; $i < scalar @data; $i++ )	{
		push @indices, $i if ( $data[$i] && ($data[$i] != 0 || $data[$i] !~ /(na|NA)/) );
	}
	
	my $count = scalar @indices;
	foreach my $index ( @indices )	{
		push @{$appears_in}, _indexToColumn($index);
	}
	return ( $count, $appears_in );
}

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





