#!/usr/bin/perl

# kmer_graph.pl
# Plug-in script to read the profiles of kmer-counted genomes and generate a graph network of compatible kmers linked by edges.

use strict;
use warnings;
use Getopt::Long qw(GetOptions);
use Cwd qw(cwd);
use List::Util qw(first shuffle uniq any);
use File::Basename;

# Required input parameters
my $input = "--";
my $output = "--";
my $outfmt;
my $sep;
my $max_gibbs;
my $header_flag;
my $v;
my $oligo_conc;
my $mg_conc;
my $monovalent_cation_conc;
my $dntp_conc;
my $temperature;
my $outsep;
my $evaluate;
my $eval_true_flag;

sub usage {
	my $usage = "kmer_graph.pl\n
	PURPOSE:            	Plug-in script to read the profiles of kmer-counted genomes and generate a graph network of compatible kmers linked by edges.
                        \n
	USAGE: kmer_graph.pl [options]
          
	=== INPUT PARAMETERS ===
	-i                      Input data matrix, assumed to be in TAB format unless otherwise specified with -sep.
	-sep                    Separator for the input data matrix.  One of 'space, 'csv', 'tab'.
	-header			INT flag, input data matrix includes a header line (Default: ON).
	-evaluate 		A list of kmers that we explicitly want to evaluate as a set, given either as a CSV list or a file with one kmer per line.
	

	=== OUTPUT PARAMETERS ===
	-out                    Output file base name, not including extensions
	-outfmt 		Output file format (either 'sif' or 'sites' [Default]).
	-outsep 		Separator for the output data.  One of 'space', 'csv', 'tab'.

	=== ADJACENCY PARAMETERS ===
	-gibbs			FLOAT; the minimum allowable Gibbs free energy calculation between a pair of aligned kmers to consider them compatable.
				Typically, the lower numbers indicate greater dimer stability.  [ Default: -6.0 ]

	=== OPTIONAL EXTRAS ===
	-v 			Verbose
	-eval_true		INT flag, used with -evaluate; only output sets that evaluate to true (don't print incompatible sets)? [ Default: OFF ]
	
	\n";
	print $usage;
}
#if ( scalar @ARGV == 0 ) { die usage(); }

GetOptions( 'in|i=s' => \$input,
			'sep|s=s' => \$sep,
            'header|h=i' => \$header_flag,
            'verbose|v' => \$v,
			'output|o=s' => \$output,
			'outfmt=s' => \$outfmt,
			'outsep=s' => \$outsep,
			'gibbs=f' => \$max_gibbs,
			'evaluate=s' => \$evaluate,
			'temperature=f' => \$temperature,
			'mgcl=f' => \$mg_conc,
			'oligo_conc=f' => \$oligo_conc,
			'dntps=f' => \$dntp_conc,
			'monovalent=f' => \$monovalent_cation_conc,
			'eval_true=i' => \$eval_true_flag,
) or die usage(); 

# Parameter setups
$header_flag = ( $header_flag && $header_flag == 0 )? 0 : 1;
$temperature = ( $temperature )? $temperature : 30;
$outfmt = ( $outfmt && $outfmt =~ /sif/ )? 'sif' : 'sites';
$max_gibbs = ( $max_gibbs && $max_gibbs != -6.0 )? $max_gibbs : -6.0;
$v = defined($v)? 1 : 0;
my $eval_flag = ( $evaluate )? 1 : 0;	
$eval_true_flag = ( $eval_flag && $eval_true_flag && $eval_true_flag != 0 )? 1 : 0;	

# Set the default separator
if    ( $sep && ($sep =~ /ssv/ || $sep =~ /space/) )	{	$sep = " ";		}
elsif ( $sep && ( $sep =~ /csv/ || $sep =~ /comma/ ) )	{	$sep = ",";		}
else 													{	$sep = "\t";	}

# Set the default output separator
if    ( $outsep && ($outsep =~ /ssv/ || $outsep =~ /space/) )	{	$outsep = " ";		}
elsif ( $outsep && ( $outsep =~ /csv/ || $outsep =~ /comma/ ) )	{	$outsep = ",";		}
else 															{	$outsep = "\t";		}

# Some variables for primer dimer thermodynamics calculations --> the user can edit these, but they are currently directly set for MEPI lab@CDC standard protocols
# Starting ionic concentration variables
$oligo_conc = ( $oligo_conc )? $oligo_conc : 250; #in nM
$mg_conc = ( $mg_conc )? $mg_conc : 2.0; #in mM
$monovalent_cation_conc = ( $monovalent_cation_conc )? $monovalent_cation_conc : 50; #in mM
$dntp_conc = ( $dntp_conc )? $dntp_conc : 0.1; #in mM

# deltaG (kcal/mol)
my %oligo_dG =( 'initC' => 0.98, 'initG' => 0.98, 'initA' => 1.03, 'initT' => 1.03 );

# deltaH (kcal/mol)
my %oligo_dH = ( 
	AA => -7.9, TT => -7.9, AT => -7.2, TA => -7.2, CA => -8.5, TG => -8.5, GT => -8.4, AC => -8.4, 
	CT => -7.8, AG => -7.8, GA => -8.2, TC => -8.2, CG => -10.6, GC => -9.8, GG => -8.0, CC => -8.0, 
	initC => 0.1, initG => 0.1, initA => 2.3, initT => 2.3,
);

my %oligo_dH_full = ( 
	AATT => -7.9, TTAA => -7.9, ATTA => -7.2, TAAT => -7.2, CAGT => -8.5, TGAC => -8.5, GTCA => -8.4, ACTG => -8.4, 
	CTGA => -7.8, AGTC => -7.8, GACT => -8.2, TCAG => -8.2, CGGC => -10.6, GCCG => -9.8, GGCC => -8.0, CCGG => -8.0,
	initC => 0.1, initG => 0.1, initA => 2.3, initT => 2.3,
		
	# Like pair mismatches 
	AATA => 1.2, ATAA => 1.2, CAGA => -0.9, AGAC => -0.9, GACA => -2.9, ACAG => -2.9, TAAA => 4.7,AAAT => 4.7, 
	ACTC => 0.0, CTCA => 0.0, CCGC => -1.5, CGCC => -1.5, GCCC => 3.6, CCCG => 3.6, TCAC => 6.1, CACT => 6.1 ,
	AGTG => -3.1, GTGA => -3.1,	CGGG => -4.9, GGGC => -4.9,	GGCG => -6.0, GCGG => -6.0,	TGAG => 1.6, GAGT => 1.6, 
	ATTT => -2.7, TTTA => -2.7,	CTGT => -5.0, TGTC => -5.0,	GTCT => -2.2, TCTG => -2.2,	TTAT => 0.2, TATT => 0.2,
		
	# G.T mismatches 
	AGTT => 1.0, TTGA => 1.0, ATTG  => -2.5, GTTA  => -2.5,	CGGT  => -4.1, TGGC  => -4.1, CTGG  => -2.8, GGTC  => -2.8,
	GGCT  => 3.3, TCGG  => 3.3,	GGTT  => 5.8, TTGG  => 5.8,	GTCG  => -4.4, GCTG  => -4.4, GTTG  => 4.1,	GTTG  => 4.1,
	TGAT  => -0.1, TAGT  => -0.1, TGGT  => -1.4, TGGT  => -1.4,	TTAG  => -1.3, 	GATT  => -1.3, 
		
	# G.A mismatches 
	AATG  => -0.6, GTAA  => -0.6, AGTA => -0.7, ATGA => -0.7, CAGG  => -0.7, GGAC  => -0.7, CGGA  => -4.0, AGGC  => -4.0,
	GACG  => -0.6, GCAG  => -0.6, GGCA  => 0.5, ACGG  => 0.5, TAAG  => 0.7, GAAT  => 0.7, TGAA  => 3.0, AAGT  => 3.0, 
		
	# C.T mismatches 
	ACTT => 0.7, TTCA => 0.7, ATTC => -1.2, CTTA => -1.2, CCGT => -0.8, TGCC => -0.8, CTGC => -1.5, CGTC => -1.5,
	GCCT => 2.3, TCCG => 2.3, GTCC => 5.2, CCTG => 5.2, TCAT => 1.2, TACT => 1.2, TTAC => 1.0, CATT => 1.0, 
		
	# A.C mismatches 
	AATC => 2.3, CTAA => 2.3, ACTA => 5.3, ATCA => 5.3, CAGC => 1.9, CGAC => 1.9, CCGA => 0.6, AGCC => 0.6, 
	GACC => 5.2, CCAG => 5.2, GCCA => -0.7, ACCG => -0.7, TAAC => 3.4, CAAT => 3.4, TCAA => 7.6, AACT => 7.6,
);
	
# deltaS (cal/K.mol)
my %oligo_dS = (
	AA => -22.2, TT => -22.2, AT => -20.4, TA => -21.3, CA => -22.7, TG => -22.7, GT => -22.4, AC => -22.4, 
	CT => -21.0, AG => -21.0, GA => -22.2, TC => -22.2, CG => -27.2, GC => -24.4, GG => -19.9, CC => -19.9, 
	initC => -2.8, initG => -2.8, initA => 4.1, initT => 4.1, sym => -1.4,
);
	
my %oligo_dS_full=(
	AATT => -22.2, 	TTAA => -22.2, ATTA => -20.4, TAAT => -21.3, CAGT => -22.7, TGAC => -22.7, 
	GTCA => -22.4, 	ACTG => -22.4, CTGA => -21.0, AGTC => -21.0, GACT => -22.2, TCAG => -22.2, 
	CGGC => -27.2, 	GCCG => -24.4, GGCC => -19.9, CCGG => -19.9, initC => -2.8, initG => -2.8,
	initA => 4.1, initT => 4.1, sym => -1.4,
		
	# Like pair mismatches	
	AATA => 1.7, ATAA => 1.7, CAGA => -4.2, AGAC => -4.2, GACA => -9.8, ACAG => -9.8, 
	TAAA => 12.9, AAAT => 12.9,	ACTC => -4.4, CTCA => -4.4, CCGC => -7.2, CGCC => -7.2, 
	GCCC => 8.9, CCCG => 8.9, TCAC => 16.4, CACT => 16.4, AGTG => -9.5, GTGA => -9.5, 
	CGGG => -15.3, GGGC => -15.3, GGCG => -15.8, GCGG => -15.8,	TGAG => 3.6, GAGT => 3.6, 
	ATTT => -10.8, TTTA => -10.8, CTGT => -15.8, TGTC => -15.8,	GTCT => -8.4, TCTG => -8.4, 
	TTAT => -1.5, TATT => -1.5,
		
	# G.T mismatches
	AGTT => 0.9, TTGA => 0.9,ATTG  => -8.3, GTTA  => -8.3,CGGT  => -11.7, TGGC  => -11.7,
	CTGG  => -8.0, GGTC  => -8.0,GGCT  => 10.4, TCGG  => 10.4, GGTT  => 16.3, TTGG  => 16.3,
	GTCG  => -12.3, GCTG  => -12.3,	GTTG  => 9.5, GTTG => 9.5, TGAT  => -1.7, TAGT => -1.7,
	TGGT  => -6.2, TGGT =>  -6.2, TTAG  => -5.3, GATT =>  -5.3, 
	
	# G.A mismatches
	AATG => -2.3, GTAA => -2.3,	AGTA => -2.3, ATGA => -2.3,	CAGG => -2.3, GGAC => -2.3,
	CGGA => -13.2, AGGC => -13.2, GACG => -1.0, GCAG => -1.0, GGCA => 3.2, ACGG => 3.2,
	TAAG => 0.7, GAAT => 0.7, TGAA => 7.4, AAGT => 7.4, 
	
	# C.T mismatches
	ACTT => 0.2, TTCA => 0.2, ATTC => -6.2, CTTA => -6.2, CCGT => -4.5, TGCC => -4.5,
	CTGC => -6.1, CGTC => -6.1,	GCCT => 5.4, TCCG => 5.4, GTCC => 13.5, CCTG => 13.5,
	TCAT => 0.7, TACT => 0.7, TTAC => 0.7, CATT => 0.7, 
	
	# A.C mismatches
	AATC => 4.6, CTAA => 4.6, ACTA => 14.6, ATCA => 14.6, CAGC => 3.7, CGAC => 3.7,
	CCGA => -0.6, AGCC => -0.6,	GACC => 14.2, CCAG => 14.2,	GCCA => -3.8, ACCG => -3.8,
	TAAC => 8.0,  CAAT => 8.0, 	TCAA => 20.2, AACT => 20.2,
);

# Recalculate oligo_dG (kcal/mol) for PCR salt conditions (the initiation values are salt independent)
# This script automatically reculates the init values for 30 degrees for WGA/MDA reactions
recalculate_dG($temperature);

print STDERR " === KMER GRAPH :: Processing inputs === \n" if ( $v == 1 );

###################################################################
# Evaluate if a given set of kmers is compatible as a clique 
# --> this is done without graph construction to make it fast, since this is easy to compute and building a very large graph is time consuming
if ( $eval_flag == 1 )	{
	# Read the given list
	my @records;
	if ( -e $evaluate )		{
		if ( $evaluate =~ /\.kmer_sets\.gini\.tab/ )	{
			open EVAL, "$evaluate" or die " === KMER GRAPH Error :: I cannot evaluate anything! URK!\n$!\n";
				@records = <EVAL>;
			close EVAL;
		}
		else	{
			$/ = "\/\/\n";
			open EVAL, "$evaluate" or die " === KMER GRAPH Error :: I cannot evaluate anything! URK!\n$!\n";
				@records = <EVAL>;
			close EVAL;
			$/ = "\n";
		}
		my $trash = shift @records if ( $header_flag == 1 );  # Throw away the header line
	}
	else	{
		@records = split(";", $evaluate);
	}
	# Parse, including removal of the annoying // record separator if needed
	for ( my $i=0; $i < scalar @records; $i++ )		{
		chomp $records[$i];
		my @sublist;
		if ( $evaluate =~ /\.kmer_sets\.gini\.tab/ )	{
			my @line = split("\t", $records[$i]);
			@sublist = split(",", $line[5]);
		}
		else	{
			@sublist = split(/,|\n/, $records[$i]);
			pop @sublist if ( $sublist[-1] eq "\/\/"  );
			for ( my $j=0; $j < scalar @sublist; $j++ )		{
				my @line = split("\t", $sublist[$j]);
				$sublist[$j] = $line[0];
			}	
		}
		$records[$i] = \@sublist;
	}
	
	my $succout = open(OUT, ">", "$output.kmers.graph.evaluated.tab") if ( $output ne "--" );
	if ( $succout ) 	{	print OUT    "Set#\tIsClique?\tVertices\n";		}
	else				{	print STDOUT "Set#\tIsClique?\tVertices\n";		}
	
	# For each list of kmers in @records, evaluate whether that set of kmers would form a clique
	# This is a binary yes or no answer, since that is the fastest evaluation we can do (still n^2 though....)
	for ( my $i=0; $i < scalar @records; $i++ )		{
		my $result = evaluate_clique( $records[$i] );
		next if ( $eval_true_flag == 1 && $result == 0 );
		if ( $succout ) 	{	print OUT    "Eval", $i+1, "\t$result\t", join(",", @{$records[$i]}), "\n";		}
		else				{	print STDOUT "Eval", $i+1, "\t$result\t", join(",", @{$records[$i]}), "\n";		}
	}

	exit;
}
###################################################################

##################################################################
# Open and read the file and populate a large hash of data
my $fh = *STDIN;
my $succin = open(IN, "<", "$input") if ( $input ne "--" && -e $input );
$fh = *IN if ( $succin ); 
my @profile = <$fh>;
my $trash = shift @profile if ( $header_flag == 1 );  # Throw away the header line
close IN if ( $succin );

# Parse the input data
my %Kmers = ();
foreach my $line ( @profile )	{
	chomp $line;
	my @kmer_profile = split("\t", $line);
	my @binding_sites;
	
	# If the format looks like it came from kmer_gini.pl script --> note that this parser ignores the pooled binding sites since we cant really tease them apart
	if ( scalar @kmer_profile == 7 )	{
		my @kmers = split(",", $kmer_profile[-2]);
		$Kmers{$_} = 1 foreach ( @kmers );
	}
	# Otherwise, assume we are given a standard binding sites profile
	else	{
		@binding_sites = split(",", $kmer_profile[1]);
		$Kmers{ $kmer_profile[0] } = \@binding_sites;
	}
	
}
###################################################################

print STDERR " === KMER GRAPH :: Constructing graph === \n" if ( $v == 1 );

###################################################################
# Generate a graph network where edges are unwieghted (by default w=1) and undirected (self-loops are not allowed).
# A pair of kmers are compatable if they do not form heterodimers, nor is one primer a subsequence of the other.
my %G = ();
foreach my $kmer ( keys %Kmers )	{
	foreach my $adjkmer ( keys %Kmers )	{
		# Skip if we've already handled this particular pair
		next if ( exists $G{$kmer}{$adjkmer} || exists $G{$adjkmer}{$kmer} );
	
		# Initialize the kmer pair as 0
		$G{$kmer}{$adjkmer} = 0; 
		
		# Skip self-loops if the two kmers are the same, or if one kmer is a subsequence of the other 
		# OR if we expect the two kmers to form an overlapping heterodimer (ie. a reasonably strong semi-global alignment)
		# This syntax should capture both outcomes.
		next if ( isSubsequence($kmer, $adjkmer) == 1 || primer_dimer($kmer, $adjkmer) < $max_gibbs );	
		
		# Otherwise, if we passed all of our checkpoints to this point, we assume this is a compatable pair of kmers
		$G{$kmer}{$adjkmer} = 1;
	}
}
###################################################################

print STDERR " === KMER GRAPH :: Analyzing the graph === \n" if ( $v == 1 );

###################################################################

# Want to do any other work/algorithms with the graph?
# Do it here

###################################################################




###################################################################

###################################################################
# Print out the SIF file --> we are NOT producing an adjacency matrix output for this because it would likely be far too enormous to be useful
# For the SIF format, we are using a SPACE SEPARATOR. Can adjust to a tab if needed, but usually Cytoscape prefers the space.
# Additionally, we are NOT reporting any lone kmers, since they are not going to be useful to us if they arent connected/compatable with anything else.
if ( $outfmt eq 'sif' )	{
	my $succout = open(OUT, ">", "$output.kmers.graph.sif") if ( $output ne "--" );
	foreach my $kmer ( keys %G )	{
		foreach my $adjkmer ( keys %G )	{
			if ( $succout )		{
				print OUT "$kmer$sep", "1$sep$adjkmer\n" if ( exists $G{$kmer}{$adjkmer} && $G{$kmer}{$adjkmer} == 1 );
			}
			else	{
				print STDOUT "$kmer$sep", "1$sep$adjkmer\n" if ( exists $G{$kmer}{$adjkmer} && $G{$kmer}{$adjkmer} == 1 );
			}
			delete $G{$adjkmer}{$kmer};
		}
	}
	close OUT if ( $succout );
}
else	{
	my %Seen = ();
	my $succout = open(OUT, ">", "$output.kmers.compatible_binding_sites.tab") if ( $output ne "--" );
	foreach my $kmer ( keys %G )	{
		next if ( exists $Seen{$kmer} );
		if ( $succout ) 	{	print OUT    "$kmer\t", join(",", @{$Kmers{$kmer}}), "\n";		}
		else				{	print STDOUT "$kmer\t", join(",", @{$Kmers{$kmer}}), "\n";		}			
		$Seen{$kmer} = 1;
		foreach my $adjkmer ( keys %G )	{
			next if ( exists $Seen{$adjkmer} );
			if ( $succout )		{ 	print OUT    "$adjkmer\t", join(",", @{$Kmers{$adjkmer}}), "\n";	}
			else				{	print STDOUT "$adjkmer\t", join(",", @{$Kmers{$adjkmer}}), "\n";	}
			$Seen{$adjkmer} = 1;
		}
	}
	close OUT if ( $succout );
}
###################################################################

exit;


#################### SUBROUTINES ####################################

# Reverse complements a nucleotide sequence, including ambiguous base calls
sub revcom	{
	my ( $seq ) = @_;
	my $revcom = reverse $seq;
	$revcom =~ tr/ACGTYRWSKMDVHBNX-/TGCARYWSMKHBDVNX-/ ;
	return ( uc $revcom );
}

# Calculates the Gini coefficient from a 1D list of non-negative, ascending sorted numbers
# The index ranges from 0 to 1, with values approaching 1 indicating increasing unevenness (we want a value close to 0)
sub gini_index		{
	my ( @distances ) = @_;
	my $n = scalar @distances;
	my $sum = sum( @distances );
	my $gini;
	for ( my $i=0; $i < $n; $i++ )	{
		$gini += ((2*($i  +1) - $n - 1) * $distances[$i]) / ($n * $sum);		# Function described: https://github.com/oliviaguest/gini
	}
	return sprintf("%1.4f", $gini);
}

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

# Computes the GC% of a putative primer and returns 1 if the GC content exceeds 80% of the primer sequence
sub gc_content		{
	my ( $primer ) = @_;
	my $gc = $primer =~ tr/GCS//;
	return ($gc / length $primer);
}

# Evaluate a clique and return TRUE (1) if a set of vertices (kmers) is indeed a clique in the graph network
sub evaluate_clique	{
	my ( $vertices ) = @_;
	my @vertices = @{$vertices};
	for ( my $i=0; $i < scalar @vertices; $i++ )	{
		for ( my $j=$i; $j < scalar @vertices; $j++ )	{
			next if ( $vertices[$i] eq $vertices[$j] );
			return 0 if ( isSubsequence($vertices[$i], $vertices[$j]) == 1 || primer_dimer($vertices[$i], $vertices[$j]) < $max_gibbs );
		}
	}
	return 1;
}


################################################################################################################################################################################# 
#
# The following code and subroutines are modified from "PerlPrimer", a Perl app written by Owen Marshall and licensed under the GNU Public License.
# Citation: Marshall OJ. PerlPrimer: cross-platform, graphical primer design for standard, bisulphite and real-time PCR. Bioinformatics 2004 20(15):2471-2472
#
################################################################################################################################################################################# 

# Returns only the complement of a nucleotide sequence, not the reverse
sub complement {
	$_ = shift;
	tr/AGCTagct/TCGAtcga/;
	return $_;
}

sub recalculate_dG {
	# because dG = dH - TdS, and dS is dependent on the salt concentration ...
	my ( $temperature ) = @_;
	
	# Big problems if [dNTPs] > [Mg++] !!  This is a quick fix ...
	my $salt_correction;
	if ($mg_conc > $dntp_conc) 		{	$salt_correction = sqrt($mg_conc - $dntp_conc);		} 
	else							{	$salt_correction = 0;								}
	
	my $na_eq=($monovalent_cation_conc + 120 * $salt_correction)/1000;
	
	# the length of each NN dimer is 2, therefore the modifier is 1
	my $entropy_adjust = (0.368 * log($na_eq));
		
	foreach my $key ( keys %oligo_dH_full ) {
		next if $key =~ /init/; # the length of each monomer is 1, thus the modifier of dS is 0 and the values are precalulated
		
		my $dS = $oligo_dS_full{$key} + $entropy_adjust;
		my $dG = $oligo_dH_full{$key}-((273.15+$temperature)*($dS/1000));
		$oligo_dG{$key} = $dG;
	}
}

sub primer_dimer {
	# This is my second attempt at implementing a primer-dimer system:
	# The first attempt aligned each combination together explicitly; this attempt
	# creates a binding matrix and then reads each primer combination from the
	# matrix.  It's not significantly faster, but it does have the advantage of
	# extending this algorithm to allow for unequal loops (although this has not
	# been implemented as yet) and providing the possiblity of reading all
	# primer-dimers (not just 3') by adjusting a single variable (which is used when
	# displaying primer information.
		
	my ($primer_f, $primer_r, $pd_full) = @_;
	return unless ($primer_f) && ($primer_r);
			
	my ($k, $l);
	my @score=();
	my %primer_hash=();
	my @score_sort=();
	my @bind_string=();
	my %rating_hash=();
	my $pd_extensible;
	
	# $pl = greatest length
	my $pfl=length($primer_f);
	my $prl=length($primer_r);
	my $pl = ($pfl>$prl ? $pfl : $prl);
	
	my $rcompr = reverse(complement($primer_r));
	my $rcomprlc = lc($rcompr);
	my $fprimer_r=lc(reverse($primer_f));
	my $rprimer_r=reverse($primer_r);
	
	# Scan the primers against each other:
	# The default is to only consider 3' primer-dimers, for speed concerns (5'
	# pd's will reduce the primer population, but won't cause extendible dimers) -
	# however, setting $pd_full to 1 will calculate all primer-dimers.  This is
	# currently used only when viewing individual primers, for both speed concerns
	# and because it's 3' primer-dimers that are the real problem in PCR.
	
	# create a binding array for each of the four bases
	for $l (0 .. $pfl-1) {
		my $mbase = substr($fprimer_r, $l, 1);
		$primer_hash{$mbase}[$l]=1;
		for $k ( 'a', 'g', 'c', 't' ) {
			$primer_hash{$k}[$l] ||=0;
		}
	}
		
	# create the primer matrix
	my @primer_comp;
	for $k (0 .. $prl-1) {
		$primer_comp[$k]=$primer_hash{substr($rcomprlc, $k, 1)};
	}
		
	# read each combination from the matrix, calculate dG for each dimer
	my $pd_len = ($pd_full ? $pfl+$prl-1 : $pl-2);
	for $k (0 .. $pd_len) {
		$score[$k]=0;
		my $bind;
		my $score_p=0;
		
		# extensible primer short-circuit - ignore all primers that will
		# not create extensible (i.e. amplifiable) dimers
		my $start = $k>$pfl-1 ? $pfl-1 : $k;
		my $end = $k>$prl-1 ? $prl-1 : $k;
		if ($pd_extensible && !$pd_full) {
			next unless $primer_comp[0][$start] == 1;
			next unless $primer_comp[$end][$start-$k] == 1;
		}
		
		# elsif ($pd_extensible) {
			# # no point reconsidering them!
			# next if $primer_comp[0][$start] == 1 && $primer_comp[$end][$start-$k] == 1;
		# }
		
		# read the binding data
		for $l (0 .. $prl-1) {
			if (($k-$l)<$pfl) {
				$bind .= $primer_comp[$l][$k-$l] if ($k-$l)>=0;
			} else {
				# spacer
				$bind .= "2";
			}
		}
		
		# Single matched bases surrounded by mismatches are unstable,
		# so we remove them with the regexp (look ahead is needed otherwise
		# strings of consecutive match/mismatches are not caught)
		$bind =~ s/01(?=[^1])/00/gx;
		
		# Short circuit if there's nothing to bind
		next unless $bind =~ /[1]/;
		
		# Find start and end of similarity
		my ($pb_init,$pb_end);
		for $l (0 .. length($bind)-1) {
			# at first I tried finding the initiating terminal bases with
			# regexps, but that was much slower ...
			if (substr($bind, $l, 1) eq "1") {
				defined($pb_init) || ($pb_init = $l);
				$pb_end=$l;
			}
		}
				
		if (defined($pb_init)) {
			# deltaG calculation
			for $l ($pb_init .. $pb_end-1) {
				next if substr($bind, $l, 2) eq "00";
				next if substr($bind, $l, 1) eq "2";
				
				#print"1: ", substr($primer_f, $pfl-$k+$l-1, 2).substr($rprimer_r, $l, 2), "  oligo_dG -> ", $oligo_dG{ substr($primer_f, $pfl-$k+$l-1, 2) . substr($rprimer_r, $l, 2) }, "  3: $score_p\n";
				$score_p += $oligo_dG{substr($primer_f, $pfl-$k+$l-1, 2).substr($rprimer_r, $l, 2)};
			}
			
			# init term corrections
			my $initterm = "init" . substr($rprimer_r, $pb_init, 1);
			$score_p += $oligo_dG{$initterm};
			
			my $endterm = "init" . substr($rprimer_r, $pb_end, 1);
			$score_p += $oligo_dG{$endterm};
			
			# add to the hash ...
			$score[$k] = sprintf("%.2f",$score_p);
			$bind_string[$k] = $bind;
			$rating_hash{$score[$k]} = $k;
		}
	}
	
	# sort the dimers to give the most stable:	
	@score_sort = sort { $a <=> $b } @score;
		
	# Returns the most stable dimer
	return $score_sort[0];
}
