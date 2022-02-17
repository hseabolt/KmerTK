#!/usr/bin/perl

# kmer_fetch.pl
# A plugin script for kmers.pl to read a data matrix and slice-n-dice based on matching rownames or colnames

# Author: MH Seabolt
# Last updated: 5-6-2020

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
my $input_sites;
my $sep;
my $header_flag;
my $rownames_flag;
my $col_cfg;
my $row_cfg;
my $cdrop;
my $rdrop;
my $where;
my $v;

sub usage {
	my $usage = "kmer_fetch.pl\n
	PURPOSE:           A plugin script for kmers.pl to read an input data matrix and slice-n-dice based on matching rownames or colnames.
                        \n
	USAGE: kmer_fetch.pl [options]
          
	=== INPUT PARAMETERS ===
	-i                       		Input data matrix, assumed to be in TAB format unless otherwise specified with -sep.
	-sep                     		Separator for the input data matrix.  One of 'space;, 'csv', 'tab'.
	-header					INT flag, input data matrix includes a header line (Default: ON).
	-rownames 				INT flag, input data matrix includes row names as the first element on each row (Default: ON).
	-b 					OPTIONAL; input binding sites profile from kmers.pl.  
						This script will automatically filter that file to contain the same kmers that pass filtering in the main profile.

	=== OUTPUT PARAMETERS ===
	-out                     		Output file base name, not including extensions
           
	=== FETCH PARAMETERS ===
	-rows           			The range of rows to retrieve, given as a:b. Multiple instructions can be passed as a csv list.
						Can be given as indices or as column names as they appear in the input data.
								Eg. FirstCol:LastCol,ThirdCol;...  ==> -rows fA:fC,B5
	-cols					Use the same syntax as for -rows.
	-drop_rows				The range of rows to drop, using the same syntax as -rows.
	-drop_cols				The range of columns to drop, using the same syntax as -rows.
	-where					An optional mathematical/logical function as an additional filter --> UNDER CONSTRUCTION!
	

	=== OPTIONAL EXTRAS ===
	-v 					verbose
	
	\n";
	print $usage;
}
if ( scalar @ARGV == 0 ) { die usage(); }

GetOptions( 'in|i=s' => \$input,
			'sites|b=s' => \$input_sites,
			'sep|s=s' => \$sep,
            'header|h=i' => \$header_flag,
			'rownames|r=i' => \$rownames_flag,
            'cols=s' => \$col_cfg,
			'rows=s' => \$row_cfg,
			'output|o=s' => \$output,
			'drop_cols|cdrop=s' => \$cdrop,
			'drop_rows|rdrop=s' => \$rdrop,
			'where=s' => \$where,
			'verbose|v' => \$v,
) or die usage(); 

# Parameter setups
$header_flag = ( $header_flag && $header_flag == 0 )? 0 : 1;
$rownames_flag = ( $rownames_flag && $rownames_flag == 0 )? 0 : 1;
$v = defined($v)? 1 : 0;

# Set the default separator
if    ( $sep && ($sep =~ /ssv/ || $sep =~ /space/) )	{	$sep = " ";		}
elsif ( $sep && ( $sep =~ /csv/ || $sep =~ /comma/ ) )	{	$sep = ",";		}
else 													{	$sep = "\t";	}	

###########################################################
#       PHASE 0: Read and validate config instructions
###########################################################

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
@headers = grep { $_ =~ /[\S]/ && $_ !~ /Kmer/ } @headers;		
_columnToIndex( $_ ) foreach ( @headers );


###################################################

# The parser for the configuration options -->
# Command line syntax is the range of columns to collapse, given in brackets with a colon, as [a:b], followed by two colons to denote the new column name.
# Multiple instructions can be passed as a semicolon-delimited list
# Eg. $ ... -cfg fA:fC,fE,B1,B2:B4,B5
die " xxx KMER FETCH Error :: No valid retrieval instructions given!  RIP. xxx \n" unless ( $col_cfg || $row_cfg || $cdrop || $rdrop || $where );

#########################
#####   COLUMNS   #######
#########################
my %Cfg = ();
if ( $col_cfg )	{
	if ( -e $col_cfg )	{
		open CGET, $col_cfg;
			my @cgets = <CGET>;
		close CGET;
		foreach my $record ( @cgets )	{
			chomp $record;
			my @record = split("$sep", $record);
			$Cfg{ $record[0] } = 1;
		}
	}
	else	{
		my @cfg_raw = split(",", $col_cfg);
		foreach my $config ( @cfg_raw )	{
			# If we've got a range of columns
			if ( $config =~ /:|-/ )		{
				my @cols = split(/:|-/, $config);
				my ( $start, $stop ) = ( $cols[0], $cols[1] );
				# If they look like numbers, then this is relatively easy
				if ( exists $ReverseIndex{$start} && exists $ReverseIndex{$stop} )	{
					( $start, $stop ) = ( $stop, $start ) if ( $stop < $start );
					foreach ( $start .. $stop )	{
						$Cfg{$headers[$_]} = 1;
					}
				}
				# If they dont look like numbers, then add the index of the column name we passed --> be careful here, the order matters!
				elsif ( exists $Index{$start} && exists $Index{$stop} )	{
					my ( $start_index, $stop_index ) = ( _columnToIndex($start), _columnToIndex($stop) );
					foreach ( $start_index .. $stop_index )	{
						$Cfg{$headers[$_]} = 1;
					}
				}
				# Otherwise... yikes, don't know what to do here
				else	{
					warn " xxx KMER FETCH Error :: Malformed config instructions given: $config xxx \n";
				}
			}
			# Else if we find a column that we don't recognize because it doesn't exist in the original data
			elsif ( not exists $Index{$config} && not exists $ReverseIndex{$config} )	{
				warn " xxx KMER CONDENSE Error :: A requested column is not present in the original data! xxx \n" if ( $v == 1 && $output ne "--" );
				next;
			}
			# Else, if we just have a list of discrete csv columns, add them directly to the list
			else 	{
				$Cfg{$config} = 1;
			}
		}
	}
}
# If we don't have any instructions to retrieve columns, then assume we want to get them all
else	{	$Cfg{$_} = 1 foreach ( @headers );		}

# If we are want to drop any columns, do it here
# The code here reads exactly the same syntax as the above retrieval code, just in reverse (ie. we delete keys from %Cfg instead of adding them)
if ( $cdrop )	{
	if ( -e $cdrop )	{
		open CDROP, $cdrop;
			my @cdrops = <CDROP>;
		close CDROP;
		foreach my $record ( @cdrops )	{
			chomp $record;
			my @record = split("$sep", $record);
			delete $Cfg{ $record[0] };
		}
	}
	else	{
		my @cfg_raw = split(",", $cdrop);
		foreach my $config ( @cfg_raw )	{
			# If we've got a range of columns
			if ( $config =~ /:|-/ )		{
				my @cols = split(/:|-/, $config);
				my ( $start, $stop ) = ( $cols[0], $cols[1] );
				# If they look like numbers, then this is relatively easy
				if ( exists $ReverseIndex{$start} && exists $ReverseIndex{$stop} )	{
					( $start, $stop ) = ( $stop, $start ) if ( $stop < $start );
					foreach ( $start .. $stop )	{
						delete $Cfg{$headers[$_]};
					}
				}
				# If they dont look like numbers, then add the index of the column name we passed --> be careful here, the order matters!
				elsif ( exists $Index{$start} && exists $Index{$stop} )	{
					my ( $start_index, $stop_index ) = ( _columnToIndex($start), _columnToIndex($stop) );
					foreach ( $start_index .. $stop_index )	{
						delete $Cfg{$headers[$_]};
					}
				}
			}
			# Else, if we just have a list of discrete csv columns, add them directly to the list
			else 	{	delete $Cfg{$config};	}
		}
	}
}
my @columns = sort keys %Cfg;

##########################
#####     ROWS     #######
##########################

###################################################
# The parsers here are more or less the same for the columns, but tweaked a bit
my %Rfg = ();
if ( $row_cfg )		{
	if ( -e $row_cfg )	{
		open RGET, $row_cfg;
			my @rgets = <RGET>;
		close RGET;
		foreach my $record ( @rgets )	{
			chomp $record;
			my @record = split("$sep", $record);
			$Rfg{ $record[0] } = 1;
		}
	}
	else	{
		my @rfg_raw = split(",", $row_cfg);
		foreach my $config ( @rfg_raw )	{
			# If we've got a range of row indices or names
			if ( $config =~ /:|-/ )		{
				my @rows = split(/:|-/, $config);
				my ( $start, $stop ) = ( $rows[0], $rows[1] );
				
				# If they look like numbers/indices, then this is the easy part
				if ( $start =~ /^[0-9]*$/ && $stop =~ /^[0-9]*$/ )	{
					( $start, $stop ) = ( $stop, $start ) if ( $stop < $start );
					$Rfg{$_} = 1 foreach ( $start .. $stop );
				}
				
				# If they dont look like numbers, then add the two flanking names and set the hash value as a start flag (2) and a stop flag (3).
				# When we are in the main loop of the program, we will have controllers to grab all rows between the start and stop flags when they are active
				elsif ( $start =~ /[A-Za-z]/ && $stop =~ /[A-Za-z]/ )	{
					$Rfg{$start} = 2;
					$Rfg{$stop}  = 3;
				}
				# Otherwise... yikes, don't know what to do here
				else	{
					warn " xxx KMER FETCH Error :: Malformed config instructions given: $config xxx \n";
				}
			}
			# Else, if we just have a list of discrete csv row indices or names, add them directly to the list
			else 	{
				$Rfg{$config} = 1;
			}
		}
	}
}
# If we are want to drop any rows, delete the keys here
if ( $rdrop )	{
	if ( -e $rdrop )	{
		open RDROP, $rdrop;
			my @rdrops = <RDROP>;
		close RDROP;
		foreach my $record ( @rdrops )	{
			chomp $record;
			my @record = split("$sep", $record);
			$Rfg{ $record[0] } = 6;
		}
	}
	else	{
		my @rfg_raw = split(",", $rdrop);
		foreach my $config ( @rfg_raw )	{
			# If we've got a range of rows
			if ( $config =~ /:|-/ )		{
				my @rows = split(/:|-/, $config);
				my ( $start, $stop ) = ( $rows[0], $rows[1] );
				
				# If they look like numbers/indices, then this is the easy part
				if ( $start =~ /^[0-9]*$/ && $stop =~ /^[0-9]*$/ )	{
					( $start, $stop ) = ( $stop, $start ) if ( $stop < $start );
					delete $Rfg{$_} foreach ( $start .. $stop );
				}
				
				# If they dont look like numbers, then add the two flanking names and set the hash value as a delete-start flag (4) and a delete-stop flag (5).
				# When we are in the main loop of the program, we will have controllers to grab all rows between the start and stop flags when they are active
				elsif ( $start =~ /[A-Za-z]/ && $stop =~ /[A-Za-z]/ )	{
					$Rfg{$start} = 4;
					$Rfg{$stop}  = 5;
				}
			}
			# Else, if we just have a discrete csv list, set another flag (6), which says to only delete this entry during the deletion logic block.
			# If this element happens to be the start of a retrieval range, then set it as (7)
			else 	{	$Rfg{$config} = ( $Rfg{$config} == 2 )? 7 : 6;	}
		}
	}
}

##########################
#####    WHERE     #######
##########################

# This parsers reads the where statement and prepares the code for computing the requested logic
# It should be given in the form -where row_name=Giardia --> this would extract all rows where the row name matches Giardia 
my $axis;
my $function;
my $operator;
my $operand;
my %Allowable = ( 'sum' => 1, 'mean' => 1, 'median' => 1, 'min' => 1, 'max' => 1, 'count' => 1, 'product' => 1, 'name' => 1 );
if ( $where )	{
	my @wheres = split(",", $where);
	foreach my $statement ( @wheres )	{
		my ( $axis, $function, $operator, $operand ) = ( $where =~ /(row|col)_(\w*)(>|<|=|!=|>=|<=)(\w*)/ );
		print "$axis $function $operator $operand\n\n\n";
		
	
		# Is method name legal?
		warn "Conditional method $where is not in the recognized form!  Will skip this criteria for fetching data.\n" unless ( $axis && $function && $operator && $operand );
	}
}

###########################################################
#       PHASE 1: Read and operate on the real data
###########################################################

# Read the data
my $fh = *STDIN;
my $succin = open(IN, "<", "$input") if ( $input ne "--" && -e $input );
$fh = *IN if ( $succin ); 

# Open our output file 
my $succout = open(OUT, ">", "$output.kmers.fetched.tab") if ( $output ne "--" );
if ( $header_flag == 1 )	{
	if ( $succout )		{	print OUT    join("$sep", @columns), "\n";		}	
	else				{ 	print STDOUT join("$sep", @columns), "\n";		}
}

# Retrieve the requested data from the input matrix
my $loop_flag = 0;
my $delete_flag = 0;
my %PassingKmers = ();
while ( <$fh> )	{
	chomp $_;
	my @row = split("$sep", $_);
	next if ( $header_flag == 1 && $. == 1 );
	
	# Set the row name and index, either the first element in the row, or the line number corrected for Perl's 0-indexing and for the pres/absence of the header line
	my $rowname = ( $rownames_flag == 1 )? shift @row : ( $header_flag == 1 )? $. - 2 : $. - 1;	
	my $rowindex = ( $header_flag == 1 )? $. - 2 : $. - 1;	
	
	# Check to see if this row is one that we think we want to retrieve based on the retrieval instructions for the rows.
	# If no rows are specified, then we assume that we want them all (ie. this step is skipped).
	if ( $row_cfg )	{
		no warnings 'uninitialized';
		#### The retrieval flag logic ladder ####
		if    ( $Rfg{$rowname} == 1 || $Rfg{$rowindex} == 1 )	{					 	$PassingKmers{$rowname} = 1; 			}		# Do nothing, we want this row. 
		elsif ( $Rfg{$rowname} == 7 || $Rfg{$rowindex} == 7 )	{	$loop_flag = 1;		$PassingKmers{$rowname} = 1;	next;	}
		elsif ( $Rfg{$rowname} == 2 || $Rfg{$rowindex} == 2 )	{	$loop_flag = 1;		$PassingKmers{$rowname} = 1;			}		# If we find a row that we want which starts a range of rows, then activate the flag and keep the row
		elsif ( $Rfg{$rowname} == 3 || $Rfg{$rowindex} == 3 )	{	$loop_flag = 0;		$PassingKmers{$rowname} = 1;			}		# De-activate the flag, but keep this row
		elsif ( $loop_flag == 1 )								{						$PassingKmers{$rowname} = 1;	 		}		# The row name won't exist explicitly in the %Rfg hash, but we want to print it as part of a range of implicitly requested rows.
		else													{														next; 	}		# We do not want this row.
		use warnings;
	}
	if ( $rdrop )	{
		#### The delete flag logic ladder ####
		no warnings 'uninitialized';
		if    ( $Rfg{$rowname} == 4 || $Rfg{$rowindex} == 4 )	{	$delete_flag = 1;	next;		}
		elsif ( $Rfg{$rowname} == 5 || $Rfg{$rowindex} == 5 )	{	$delete_flag = 0;	next;		}
		elsif ( $Rfg{$rowname} == 6 || $Rfg{$rowindex} == 6 )	{						next;		}
		elsif ( $delete_flag == 1 )								{						next;		}
		else													{	$PassingKmers{$rowname} = 1;	}		# Do nothing here, this row is one that we want to keep.
		use warnings;
	}
	
	# @new_line here is an array containing the selected column values that we want, we will join it when we print out the output line
	my @new_line;
	
	# If we match a row, then get the columns that we need and print it out 
	foreach my $col ( @columns )	{
		my $index = _columnToIndex( $col );
		push @new_line, $row[$index];
	}
	
	# If we have a 'where' conditional, compute those options here
	# UNDER CONSTRUCTION!
	
	

	# Note that if we didn't have row names with our input, then we are not adding them here, since I assume that omitting them was intentional in the first place.
	if ( $rownames_flag == 1 )	{
		if ( $succout ) 	{ 	print OUT    "$rowname$sep", join("$sep", @new_line), "\n";		}
		else 				{ 	print STDOUT "$rowname$sep", join("$sep", @new_line), "\n";		}
	}
	else	{
		if ( $succout ) 	{ 	print OUT    join("$sep", @new_line), "\n";		}
		else 				{ 	print STDOUT join("$sep", @new_line), "\n";		}
	}
}
close IN if ( $input ne "--" );
close OUT if ( $succout );

# Filter the binding sites as well if we have them.
if ( $input_sites && -e $input_sites )		{
	open SITES, "$input_sites" or warn "Cannot open input binding sites data! $!\n";
	open OUT2, ">", "$output.kmers.binding_sites.fetched.tab";
	while ( <SITES> )	{
		chomp $_;
		my @data = split("\t", $_);
		print OUT2 "$_\n" if ( exists $PassingKmers{$data[0]} );
	}
	close SITES;
	close OUT2;
}

exit;


##################### SUBROUTINES ###################################

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


