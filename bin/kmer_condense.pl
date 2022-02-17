#!/usr/bin/perl

# kmer_condense.pl
# A plugin script for kmers.pl to read a data matrix and collapse multiple columns into a single column in order to reduce computational complexity downstream.
# Requires that the input data matrix has a header line!

# Author: MH Seabolt
# Last updated: 5-6-2020

use strict;
use warnings;
use Getopt::Long qw(GetOptions);
use Cwd qw(cwd);
use List::Util qw(first shuffle uniq any);
use Scalar::Util qw(looks_like_number);
use File::Basename; 
use Parallel::ForkManager;

use lib '/media/hunter/Data/scripts';
use lib 'D:\scripts';
use lib '/scicomp/home/ngr8/Biolinux/scripts';

# Required input parameters
my $input = "--";
my $output = "--";
my $sep;
my $cfg;
my $compute;
my $zeroAsNA;
my $header_flag;
my $rownames_flag;
my $v;
my $threads;

sub usage {
	my $usage = "kmer_condense.pl\n
	PURPOSE:           A plugin script for kmers.pl to read a data matrix and collapse multiple columns into a single column in order to reduce computational complexity downstream.
					   Requires that the input data matrix has a header line!
                        \n
	USAGE: kmer_condense.pl [options]
          
	=== INPUT PARAMETERS ===
	-i                       Input data matrix, assumed to be in TAB format unless otherwise specified with -sep.
	-sep                     Separator for the input data matrix.  One of 'space', 'csv', 'tab'.
	-header					 INT flag, input data matrix includes a header line (Default: ON).
	-rownames 				 INT flag, input data matrix includes row names as the first element on each row (Default: ON).
	-cfg					 Input configuration file (expecting TAB format) in the form: 
									COL_HEADER		NEW_COL_NAME
							 Command line syntax is the range of columns to collapse, given in brackets with a colon, as [a:b],
									followed by two colons to denote the new column name.  Multiple instructions can be passed as a semicolon-delimited list.
									Eg. [FirstCol:LastCol]::NewCol;...  ==> -cfg [fA:fC,fe]::Giardia;[B1,B2:B4,B5]::Microbiome	
							 Can be given as indices or as column names as they appear in the input data.

	=== OUTPUT PARAMETERS ===
	-out                     Output file base name, not including extensions
	
           
	=== CONDENSING PARAMETERS ===
	-compute                 Calculation to compute when condensing columns -- one of 'sum', 'mean', median', 'count', 'max', 'min', 'random', 'choose'.  [Default: sum]
	-zeroAsNA				 OPTIONAL INT flag; Do you want to treat zeros in the data as undefined/NA values? [Default: OFF]	

	=== OPTIONAL EXTRAS ===
	-v 					Verbose
	-t 					INT; number of threads to use for multi-threading
	
	\n";
	print $usage;
}
if ( scalar @ARGV == 0 ) { die usage(); }

GetOptions( 'in|i=s' => \$input,
			'sep|s=s' => \$sep,
			 'cfg=s' => \$cfg,
			 'header|h=i' => \$header_flag,
			'rownames|r=i' => \$rownames_flag,
            'compute=s' => \$compute,
			'zeroAsNA|z=i' => \$zeroAsNA,
			'output|o=s' => \$output,
			'verbose|v' => \$v,
		  'threads|t=i' => \$threads,
) or die usage(); 

# Parameter setups
$v = defined($v)? 1 : 0;
$header_flag = ( $header_flag && $header_flag == 0 )? 0 : 1;
$rownames_flag = ( $rownames_flag && $rownames_flag == 0 )? 0 : 1;
$zeroAsNA = ( $zeroAsNA && $zeroAsNA == 1 )? 1 : 0;
$threads = ( $threads && $threads >= 1 )? $threads : 1;

# Set the default separator
if    ( $sep && ($sep =~ /ssv/ || $sep =~ /space/) )	{	$sep = " ";		}
elsif ( $sep && ( $sep =~ /csv/ || $sep =~ /comma/ ) )	{	$sep = ",";		}
else 													{	$sep = "\t";	}	

###########################################################
#       PHASE 0: Read and validate config instructions
###########################################################

# Set the default collapsing computational
my %Allowable = ( 'sum' => 1, 'mean' => 1, 'median' => 1, 'min' => 1, 'max' => 1, 'random' => 1, 'choose' => 1, 'count' => 1, 'product' => 1 );
if    ( $compute && exists $Allowable{$compute} )	{ ;	}		# Great, we are happy as is, move on.
else 	{	$compute = "sum";	}


###############################################
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

###################################################

# The parser for the configuration options -->
# Command line syntax is the range of columns to collapse, given in brackets with a colon, as [a:b], followed by two colons to denote the new column name.
# Multiple instructions can be passed as a semicolon-delimited list
# Eg. [FirstCol:LastCol]::NewCol;...  ==> -cfg [fA:fC,fe]::Giardia;[B1,B2:B4,B5]::Microbiome
my %Cfg = ();
my %SeenCols = ();
if ( not -e $cfg )	{
	my @cfg_raw = split(";", $cfg);
	foreach my $config ( @cfg_raw )	{
		
		# Determine what sort of instructions we have
		my ( $old_names, $new_name ) = ($config =~ /\[(.*)\]::(.*)/ );
		warn " xxx KMER CONDENSE Error :: Malformed config instructions given: $config xxx \n" if ( not $old_names || not $new_name );
		my @old_names = split(",", $old_names);
		
		foreach my $old_name ( @old_names )	{
			# If we've got a range of columns
			if ( $old_name =~ /:|-/ )		{
				my @cols = split(/:|-/, $old_name);
				my ( $start, $stop ) = ( $cols[0], $cols[1] );
				# If they look like numbers, then this is relatively easy
				if ( exists $ReverseIndex{$start} && exists $ReverseIndex{$stop} )	{
					( $start, $stop ) = ( $stop, $start ) if ( $stop < $start );
					foreach ( $start .. $stop )	{
						$Cfg{$new_name}{$headers[$_]} = 1;
						$SeenCols{$headers[$_]} = 1;
					}
				}
				# If they dont look like numbers, then add the index of the column name we passed --> be careful here, the order matters!
				elsif ( exists $Index{$start} && exists $Index{$stop} )	{
					my ( $start_index, $stop_index ) = ( _columnToIndex($start), _columnToIndex($stop) );
					foreach ( $start_index .. $stop_index )	{
						$Cfg{$new_name}{$headers[$_]} = 1;
						$SeenCols{$headers[$_]} = 1;
					}
				}
				# Otherwise... yikes, don't know what to do here
				else	{
					die " xxx KMER CONDENSE Error :: Malformed config instructions given: $config xxx \n";
				}
			}
			# Else if we find a column that we don't recognize because it doesn't exist in the original data
			elsif ( not exists $Index{$old_name} && not exists $ReverseIndex{$old_name} )	{
				warn " xxx KMER CONDENSE Error :: A requested column is not present in the original data! xxx \n" if ( $v == 1 && $output ne "--" );
				next;
			}
			# Else, if we just have a list of discrete csv columns, add them directly to the list
			else 	{
				$Cfg{$new_name}{$old_name} = 1;
				$SeenCols{$old_name} = 1;
			}
		}
	}
}
# Alternatively, it's probably easier to just use a config file if we have very complex instructions to pass.
# Open the configuration file and read it
elsif ( -e $cfg )	{
	open CFG, "$cfg" or die "Cannot read the input configuration file!\n$!\n";
		my @cfg_raw = <CFG>;
	close CFG;
	foreach my $line ( @cfg_raw )	{
		chomp $line;
		my @line = split("\t", $line);
		$Cfg{$line[1]}{$line[0]} = 1;			# Adds the column header to the factor list in the hash, where the 'factor' is the new column header that we are collapsing into.
		$SeenCols{$line[0]} = 1;
	}
}
# Otherwise, URK something is wrong! DIE!
else	{
	die " xxx KMER CONDENSE Error :: No valid configuration instructions given!  RIP. xxx \n";
}
# Now we need to account for any columns that were not covered in the config instructions 
foreach my $colname ( @headers )	{
	# If the element in header, the existing colname, exists within the %Cfg hash, then do nothing.
	next if ( exists $SeenCols{$colname} );	
	
	# Otherwise, add this column to the %Cfg hash so that it will be processed downstream
	$Cfg{$colname}{$colname} = 1;
}
my @factors = sort keys %Cfg;

###########################################################
#       PHASE 1: Read and operate on the real data
###########################################################

# Read the data
my $fh = *STDIN;
my $succin = open(IN, "<", "$input") if ( $input ne "--" && -e $input );
$fh = *IN if ( $succin ); 

# Instantiate out Fork Manager and open our output file 
my $manager = Parallel::ForkManager->new($threads);
my $succout = open(OUT, ">", "$output.kmers.condensed.tab") if ( $output ne "--" );

if ( $header_flag == 1 )	{
	if ( $succout )		{	print OUT    join("$sep", @factors), "\n";		}	
	else				{ 	print STDOUT join("$sep", @factors), "\n";		}
}

# Condense the columns
my $x = 0;
while ( <$fh> )	{
	$manager->start and next if ( $threads > 1 );
	chomp $_;
	my @row = split("$sep", $_);
	next if ( $header_flag == 1 && $. == 1 );
		
	# Set the row name, either the first element in the row, or just the line number
	my $rowname = ( $rownames_flag == 1 )? shift @row : ( $header_flag == 1 )? $. - 2 : $. - 1;	

	# @new_line here is an array containing the collapsed values that we want, we will join it when we print out the output line
	my @new_line;
	
	# Collapse the columns that we want and perform the requested computation
	for ( my $i=0; $i < scalar @factors; $i++ )	{
		my @cols_to_collapse = sort keys %{ $Cfg{$factors[$i]} };
		
		# Get all the data that we need into a new list
		my @collapsable_data;
		foreach my $col ( @cols_to_collapse )	{
			my $index = _columnToIndex( $col );
			if ( $zeroAsNA == 1 && ($row[$index] == 0 || $row[$index] =~ /(na|NA)/) ) 		{ 					;						}	# Ignore this value if it is zero or NA
			else																			{	push @collapsable_data, $row[$index];	}
		}
		
		# Perform the requested computation
		my $calculation = "NA";
		if 	  ( $compute eq "mean"    )	{	$calculation = sprintf("%3.3f", mean(@collapsable_data));		}
		elsif ( $compute eq "median"  )	{	$calculation = median(@collapsable_data);	}
		elsif ( $compute eq "random"  )	{	$calculation = choose(@collapsable_data);	}
		elsif ( $compute eq "choose"  )	{	$calculation = choose(@collapsable_data);	}
		elsif ( $compute eq "count"   ) {	$calculation = count(@collapsable_data);	}
		elsif ( $compute eq "min" 	  ) {	$calculation = min(@collapsable_data);		}
		elsif ( $compute eq "max" 	  )	{	$calculation = max(@collapsable_data);		}
		elsif ( $compute eq "product" )	{	$calculation = product(@collapsable_data);	}
		else							{	$calculation = sum(@collapsable_data);		}
		
		# Add the new calculation to @new_line
		$new_line[$i] = $calculation;
	}
	
	# Print the line if the kmer is considered unique.
	# Note that if we didn't have row names with our input, then we are not adding them here, since I assume that omitting them was intentional in the first place.
	if ( $rownames_flag == 1 )	{
		if ( $succout ) 	{ 	print OUT    "$rowname$sep", join("$sep", @new_line), "\n";		}
		else 				{ 	print STDOUT "$rowname$sep", join("$sep", @new_line), "\n";		}
	}
	else	{
		if ( $succout ) 	{ 	print OUT    join("$sep", @new_line), "\n";		}
		else 				{ 	print STDOUT join("$sep", @new_line), "\n";		}
	}
	
	# Increment $x and terminate the child process
	$x++;
	$manager->finish if ( $threads > 1 );
}
close IN if ( $input ne "--" );
close OUT if ( $succout );

exit;



##################### SUBROUTINES #########################################

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

########  Descriptive Statistics Subroutines  ############

# Computes the mean value of a list
sub mean	{ 
	( scalar @_ == 0 )? return 0 : return sum(@_)/@_;
}

# Returns the sum of a numerical array
sub sum	{
	my ( @numbers ) = ( scalar @_ == 0 )? return 0 : @_;
	my $sum = 0;
	foreach ( @numbers )	{	$sum += $_;		}
	return $sum;
}

# Returns the sum of a numerical array
sub product	{
	my ( @numbers ) = ( scalar @_ == 0 )? return 0 : @_;
	my $prod = 0;
	foreach ( @numbers )	{	$prod *= $_;		}
	return $prod;
}

# Returns the max value from a set of numerical args
sub max 	{
	my (@l) = ( scalar @_ == 0 )? return "NA" : @_;
	my $max = $l[0];
	foreach my $x ( @l )	{	$max = $x if ( $x > $max );		}
	return $max;
}

# Returns the min value from a set of numerical args
sub min	{
	my (@l) = ( scalar @_ == 0 )? return "NA" : @_;
	my $min = $l[0];
	foreach my $x ( @l )	{
		$min = $x if ( $x < $min );
	}
	return $min;
}

# Returns the median value of an array list
sub median {
    my (@data) = ( scalar @_ == 0 )? return 0 : sort { $a <=> $b } @_;
    if ( scalar(@data) % 2 ) {
        return ( $data[ @data / 2 ] );
    } 
    else {
        my ( $upper, $lower );
        $lower = $data[ @data / 2 ];
        $upper = $data[ @data / 2 - 1 ];
        return ( mean( $lower, $upper ) );
    }
}

# Randomly selects an element from a given array
sub choose	{
	my (@data) = ( scalar @_ == 0 )? return "NA" : @_;
	my $index = int(rand( scalar @data ));
	return $data[$index];
}

# Returns the number of non-empty elements in a given array
sub count	{
	( scalar @_ == 0 )? return 0 : return scalar( grep { defined($_) && $_ > 0 } @_ );
}
