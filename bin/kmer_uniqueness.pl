#!/usr/bin/perl

# kmer_uniqueness.pl
# A plugin script for kmers.pl to read a matrix of the frequency counts and rank the uniqueness of a kmer relative to a set of genomes.


# Author: MH Seabolt
# Last updated: 4-21-2020

use v5.24.1;
use strict;
use warnings;
use Getopt::Long qw(GetOptions);

use lib '/media/hunter/Data/scripts';
use lib 'D:\scripts';
use lib '/scicomp/home/ngr8/Biolinux/scripts';

# Required input parameters
my $input = "--";
my $output = "--";
my $sep;
my $header_flag;
my $rownames_flag;
my $max_genomes;
my $must_include;
my $dont_include;
my $v;
my $threads;
my $ranking_range;
my $tally;
my $ingroup;

sub usage {
	my $usage = "kmer_uniqueness.pl\n
	PURPOSE:           A plugin script for kmers.pl to read a matrix of the frequency counts and rank the uniqueness and specificity of a kmer relative to a set of genomes.
                        \n
	USAGE: kmer_uniqueness.pl [options]
          
	=== INPUT PARAMETERS ===
	-i                       input data matrix, assumed to be in TAB format unless otherwise specified with -sep.
	-sep                     separator for the input data matrix.  One of 'space', 'csv', 'tab'.
	-header					 INT flag, input data matrix includes a header line (Default: ON).
	-rownames 				 INT flag, input data matrix includes row names as the first element on each row (Default: ON).

	=== OUTPUT PARAMETERS ===
	-out                     output file base name, not including extensions
	
	=== UNIQUENESS PARAMETERS ===
	-max_genomes             INT; the maximum number of genomes that a given kmer can appear in to be considered 'unique'. [Default: all genomes] 
	-ranking_range			 The range of rankings that we want to keep.  Eg: -ranking_range 1,2,3
	-must_include			 A CSV list of genomes that a kmer must appear in at least one (and meet the other conditionals) to be returned as output.  Think of this like grep :)
	-dont_include			 A CSV list of genomes that a kmer must not appear in.  Think of this like grep -v.
	-ingroup 				 A CSV list of genomes that are in-group genomes for the purpose of calculating specifity ratios.
	
	=== OPTIONAL EXTRAS ===
	-v 					verbose
	-tally 				INT; Output a tallied list of the frequency that a genome appears in the rankings [Default: OFF]
	
	\n";
	print $usage;
}
if ( scalar @ARGV == 0 ) { die usage(); }

GetOptions( 'in|i=s' => \$input,
			'sep|s=s' => \$sep,
            'header|h=i' => \$header_flag,
			'rownames|r=i' => \$rownames_flag,
            'max_genomes=i' => \$max_genomes,
			'ranking_range=s' => \$ranking_range,
			'must_include=s' => \$must_include,
			'dont_include=s' => \$dont_include,
			'output|o=s' => \$output,
			'verbose|v' => \$v,
		  'tally=i' => \$tally,
		  'ingroup=s' => \$ingroup,
) or die usage(); 

# Parameter setups
$max_genomes = ( $max_genomes && int($max_genomes) >= 1 )? $max_genomes : 'inf';
$header_flag = ( $header_flag && $header_flag == 0 )? 0 : 1;
$rownames_flag = ( $rownames_flag && $rownames_flag == 0 )? 0 : 1;
$v = defined($v)? 1 : 0;
$tally = ( $tally && $tally == 1 )? 1 : 0;

# Set the default separator
if    ( $sep && ($sep =~ /ssv/ || $sep =~ /space/) )	{	$sep = " ";		}
elsif ( $sep && ( $sep =~ /csv/ || $sep =~ /comma/ ) )	{	$sep = ",";		}
else 													{	$sep = "\t";	}	

###################################################################
# The parser for the ranking range -- this is a terrible parser, but its very small, so I dont care :)
my %Ranks = ( '1-TotallyUnique' => 0, '2-MostlyUnique', => 0, '3-SomewhatUnique' => 0, '4-Variable' => 0, '5-SomewhatUbiquituous' => 0, '6-MostlyUbiquitous' => 0, '7-TotallyUbiquitous' => 0 );
if ( $ranking_range )	{
	my @range = split(",", $ranking_range);
	foreach my $config ( @range )	{
		# If we've got a range
		if ( $config =~ /:|-/ )		{
			my @ranks = split(/:|-/, $config);
			my ( $start, $stop ) = ( $ranks[0], $ranks[1] );
			for ( my $i=$start; $i <= $stop; $i++ )	{
				my @key = grep { $_ =~ /^$i/ } keys %Ranks;
				$Ranks{$key[0]} = 1;
			}
		}
		else	{
			my @key = grep { $_ =~ /^$config/ } keys %Ranks;
			$Ranks{$key[0]} = 1;
		}
	}
}
# Otherwise, assume we want them all
else	{
	$Ranks{$_} = 1 foreach ( keys %Ranks );
}
###################################################################

###################################################################
# The parser for the -targets input, if we have it.
# This builds a lookup table of ingroup genomes (if it ain't in this hash, it's not an ingroup genome)
my %Ingroup = ();
if ( $ingroup )		{
	my @ingroups;
	if ( -e $ingroup )	{
		open INGROUP, $ingroup or die " === KMER UNIQUENESS Error :: Cannot open -ingroups file!\n$!\n";
			@ingroups = <INGROUP>;
		close INGROUP;
	}
	else	{
		chomp $ingroup;
		@ingroups = split(",", $ingroup);
	}
	foreach ( @ingroups )	{
		chomp $_;
		$Ingroup{$_} = 1;
	}
}
###################################################################

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

###################################################################
# The parser for $must_include, which operates on the headers
my %Grep = ();
if ( $must_include )	{
	if ( -e $must_include )	{
		open INCL, "$must_include" or warn "$!\n";
			my @includes = <INCL>;
		close INCL;
		foreach ( @includes )	{
			chomp $_;
			my @line = split("$sep", $_);
			$Grep{$line[0]} = 1;
		}
	}
	else	{
		my @includes = split(",", $must_include);
		$Grep{$_} = 1 foreach (@includes);
	}
}
else	{
	$Grep{$_} = 1 foreach (@headers);		# Just get them all.
}
if ( $dont_include )	{
	if ( -e $dont_include )	{
		open NOPE, "$dont_include" or warn "$!\n";
			my @includes = <NOPE>;
		close NOPE;
		foreach ( @includes )	{
			chomp $_;
			my @line = split("$sep", $_);
			delete $Grep{$line[0]} if ( exists $Grep{$line[0]} );
		}
	}
	else	{
		my @includes = split(",", $dont_include);
		delete $Grep{$_} foreach (@includes);
	}
}
###################################################################

# Calculate the rankings thresholds:
my $percentile20 = int( scalar @headers * 0.20 );
my $percentile40 = int( scalar @headers * 0.40 );
my $percentile60 = int( scalar @headers * 0.60 );
my $percentile80 = int( scalar @headers * 0.80 );

###########################################################
#       PHASE 1: Read and operate on the real data
###########################################################


# Read the data
my $fh = *STDIN;
my $succin = open(IN, "<", "$input") if ( $input ne "--" && -e $input );
$fh = *IN if ( $succin ); 

# Instantiate out Fork Manager and open our output file 
my $succout = open(OUT, ">", "$output.kmers.uniqueness_rankings.tab") if ( $output ne "--" );

if ( $header_flag == 1 )	{
	if ( $succout )		{	print OUT    "TotalCount$sep","IngroupCount$sep","SpecificityRatio$sep","Ranking$sep","Appears_in\n";		}	
	else				{ 	print STDOUT "TotalCount$sep","IngroupCount$sep","SpecificityRatio$sep","Ranking$sep","Appears_in\n";		}	
}

# If we want to tally the genomes that we find
my $succout2 = open(TALLY, ">", "$output.tally") if ( $tally == 1 );
my %Tally = ();

# Open an output file to write a histogram of the 7 buckets we have
open HIST, ">", "$output.kmer_uniqueness.histogram";
my %Histogram = ( '1-TotallyUnique' => 0, '2-MostlyUnique', => 0, '3-SomewhatUnique' => 0, '4-Variable' => 0, '5-SomewhatUbiquituous' => 0, '6-MostlyUbiquitous' => 0, '7-TotallyUbiquitous' => 0 );

# Using a while loop since I assume that the input file will be too large to read into RAM to work on
my $x = 0;
while ( <$fh> )	{
	chomp $_;
	my @row = split("$sep", $_);
	next if ( $header_flag == 1 && $. == 1 );
	
	# Set the row name, either the first element in the row, or just the line number
	my $rowname = ( $rownames_flag == 1 )? shift @row : ( $header_flag == 1 )? $. - 2 : $. - 1;		
	
	# Count the non-zero elements in the row, skipping if by chance we encounter a row with a count of zero (indicates some filtering questions)
	my ( $count, $appears_in ) = count(@row);
	next if ( $count == 0 );
	
	# Rank the uniqueness of the kmer
	my $rank;
	if    ( $count == 1 )											{	$rank = "1-TotallyUnique";			}
	elsif ( $count >= 2 && $count < $percentile20 )					{	$rank = "2-MostlyUnique";			}
	elsif ( $count >= $percentile20 && $count < $percentile40 )		{	$rank = "3-SomewhatUnique";			}
	elsif ( $count >= $percentile40 && $count < $percentile60 )		{	$rank = "4-Variable";				}
	elsif ( $count >= $percentile60 && $count < $percentile80 )		{	$rank = "5-SomewhatUbiquituous";	}
	elsif ( $count >= $percentile80 && $count < scalar @headers )	{	$rank = "6-MostlyUbiquitous";		}
	elsif ( $count == scalar @headers )								{ 	$rank = "7-TotallyUbiquitous";		}
	else															{	$rank = "??";						}
	
	# Check if an 'appears in' list contains the required genomes if we are using the $must_include options
	my $not_included = 0;
	foreach ( @{$appears_in} )	{
		last if ( $not_included > 0 );
		$not_included++ if ( exists $Grep{$_} );
	}
	next if ( $not_included == 0 );
	
	# If we are tallying the genomes
	if ( $succout2 && $count <= $max_genomes && $Ranks{$rank} == 1 )	{
		foreach ( @{$appears_in} )	{
			$Tally{$_}++;
		}
	}
	
	# Update the frequency counts in the histogram
	$Histogram{$rank}++;
	
	# Calculate the specificity ratio --> use NA if we dont have any ingroup genomes
	# Specificity ratio is calculated as ( n Ingroups / m Outgroups )
	my $specificity_ratio = "NA";
	my $ingroup_count = "NA";
	if ( scalar keys %Ingroup > 0 )	{
		$ingroup_count = 0;
		foreach ( @{$appears_in} )	{
			$ingroup_count++ if ( exists $Ingroup{ $_ } );
		}
		$specificity_ratio = ($count == $ingroup_count)? sprintf("%3.3f", ($ingroup_count / ($count))) : sprintf("%3.3f", ($ingroup_count / ($count - $ingroup_count)));
	}
	
	# Print the line if the kmer is considered DISTINCT ( appears in less than the maximum number of GENOMES)
	if ( $succout ) 	{	print OUT    "$rowname$sep$count$sep$ingroup_count$sep$specificity_ratio$sep$rank$sep", join(",", @{$appears_in}), "\n" if ( $count <= $max_genomes && $Ranks{$rank} == 1 );		}
	else				{ 	print STDOUT "$rowname$sep$count$sep$ingroup_count$sep$specificity_ratio$sep$rank$sep", join(",", @{$appears_in}), "\n" if ( $count <= $max_genomes && $Ranks{$rank} == 1 );		}

	# Increment x
	$x++;
}
close IN if ( $input ne "--" );
close OUT if ( $succout );

# Print out the tally file
if ( $succout2 ) 	{
	my $valuesum = sum( values %Tally );
	print TALLY "Genome\tTally\tPercent\n";
	foreach ( sort { $Tally{$b} <=> $Tally{$a} } keys %Tally )	{
		my $percent = sprintf("%3.2f", (100 * ( $Tally{$_} / $valuesum )) );
		print TALLY "$_\t$Tally{$_}\t$percent\n" ;		# Sorts the %Tally hash in descending order based on the values.
	}
}
close TALLY if ( $succout2 );

# Print out the histogram file
print HIST "Rank\tCount\n";
foreach my $key ( sort keys %Histogram )	{
	print HIST "$key\t$Histogram{$key}\n";
}
close HIST;

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

# Simple arithmetic sum of a list of numbers
sub sum 	{
	my ( @list ) = @_;
	my $sum = 0;
	return "NA" if ( scalar @list == 0 );
	foreach ( @list )	{
		$sum += $_;
	}
	return $sum;
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





