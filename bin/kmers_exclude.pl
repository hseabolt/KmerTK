#!/usr/bin/perl

# kmers.pl
# A naive kmer-counting wrapper that constructs a SQLite database of kmers and their binding sites within genomes.
# Be warned that this no effort has been made here to make this space or memory efficient, 
# so expect this to eat up alot of system resources that scales with the size of the genomes

# Author: MH Seabolt
# Last updated: 4-10-2020

use strict;
use warnings;
use Getopt::Long qw(GetOptions);
use Cwd qw(cwd);
use List::Util qw(first shuffle uniq any);
use File::Basename; 
use Parallel::ForkManager;

# Required input parameters
my $input_list;
my $genome_list;
my $scoring;
my $output;
my $mink;
my $maxk;
my $step;
my $min_freq;
my $min_tm;
my $max_tm;
my $max_gini;
my $min_bind_dist;
my $max_bind_dist;
my $allow_ambig;
my $organelles;
my $v;
my $threads;

sub usage {
	my $usage = "kmers.pl\n
	PURPOSE:            A naive kmer-counting wrapper that outputs TAB files of kmers and their binding sites in a series of genomes.
                        \n
	USAGE: kmers.pl [options]
          
	=== INPUT PARAMETERS ===
	-i                       input CSV list of FASTA formatted genome files
	-il                      input file containing a list of genome files, one per line   

	=== OUTPUT PARAMETERS ===
	-out                     output file base name, not including extensions
	-score			INT flag; compute detailed kmer scoring? [Default: OFF]
           
	=== KMER COUNTING PARAMETERS ===
	-mink                    INT; the minimum value of k [Default: 8]
	-maxk                    INT; the maximum value of k [Default: 12]
	-step                    INT; the step size to use for comparing values of k. Default = 1
	
	=== FILTERING PARAMETERS ===
	-min_freq			int; discard a primer if it does not bind to the foreground genome with at least this frequency.  Default: 5
	-min_tm				int; discard a primer if the estimated melting temperature is below this threshold value.  Default: 15
	-max_tm				int; discard a primer if the estimated melting temperature is below this threshold value.  Default: 45
	-max_gini				int; discard a primer if the individual-based Gini index is above this threshold value.  Default: 1.00
	-allow_ambig			int; allow primer sequences to contain up to X ambiguous bases. Default: 0
	-min_bind_dist			int; discard a kmer if the mean binding distance is less than this threshold value.  Default: 5000
	-max_bind_dist			int; discard a kmer if the mean binding distance is greater than this threshold value.  Default: 1500000
	-organelles				INT flag; How to handle contigs in the genomes that are identified as organelles (mitochondria, plasmids, apicoplasts, chloroplasts )?
								0 == SKIP organelles,  1 == count kmers on ALL contigs [Default],  2 == count kmers ONLY on organelles

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
            'step=i' => \$step,
			'min_freq=i' => \$min_freq,
            'verbose|v' => \$v,
			'output|o=s' => \$output,
			'score=i' => \$scoring,
            'min_tm=i' => \$min_tm,
            'max_tm=i' => \$max_tm,
            'max_gini=f' => \$max_gini,
		  'min_bind_dist=i' => \$min_bind_dist,
		  'max_bind_dist=i' => \$max_bind_dist,
		  'allow_ambig=i' => \$allow_ambig,
		  'organelles=i' => \$organelles,
		'threads|t=i' => \$threads,
			
) or die usage(); 

# Parameter setups
$mink = ( $mink && int($mink) >= 5 )? $mink : 8;
$maxk = ( $maxk && int($maxk) > int($mink) && int($maxk) <= 30 )? $maxk : $mink;
$step = ( $step && int($step) > 1 )? $step : 1;
$scoring = ( $scoring && $scoring == 1 )? 1 : 0;
$organelles = ( $organelles && $organelles =~ /[0-2]/ )? $organelles : 1;
$min_freq = ( $min_freq >= 1 )? $min_freq : 5;
$min_tm = ( $min_tm && $min_tm >= 15 )? $min_tm : 15;
$max_tm = ( $max_tm && $max_tm <= 80 )? $max_tm : 45;		
$min_bind_dist = ( $min_bind_dist && $min_bind_dist >= 5000 )? $min_bind_dist : 5000;
$max_bind_dist = ( $max_bind_dist && $max_bind_dist <= 1500000 )? $max_bind_dist : 1500000;		# 1,500,000
$allow_ambig = ( $allow_ambig && $allow_ambig <= $maxk )? $allow_ambig : 0;
$max_gini = ( $max_gini && $max_gini <= 1.0 && $max_gini >= 0.0 )? $max_gini : 1.0;
$v = defined($v)? 1 : 0;
$threads = ( $threads && $threads >= 1 )? $threads : 1;

# Notes on the "organelles" handling flag:
# 	0 = skip non-chromosomal contigs entirely, 
#	1 = include ALL contigs, disregarding labelling, 
# 	2 = count kmers ONLY on non-chromosomal contigs

##############################################
#       PHASE 0: INITIAL SETUP
##############################################

# Set up the possible values of k if we are using multiple ones
my @ks;
my @proc_order;
for ( my $k=$mink; $k <= $maxk; $k+=$step )                      {
	push @ks, $k;
}

# If it appears that we just gave a file containing a list of genomes, one per line
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

##############################################
#       PHASE 1: KMER COUNTS
##############################################

# The approach to kmer counting currently in use is a naive one, which may be able to be sped up with 
# a better algorithm (eg. Jellyfish, BFcounter, DSK), however, regardless of the algorithm used, 
# kmer counting is not a space or time-efficient process.

# This phase of the pipeline is computationally expensive to compute so be prepared to be a little patient here.
print STDERR " ===        KMER COUNTS        === \n" if ( $v == 1 );

# Instantiate the Fork Manager
my $manager = Parallel::ForkManager->new($threads);

# Initialize the multidimensional hashes and some simple lists
my %Genomes = ();
my %BindingSites = ();
my @names;
my @sites;
foreach my $genome_file ( @proc_order )	{
	# Get the genome's basename
	chomp $genome_file;
	my $genome = fileparse($genome_file);
	$genome =~ s/\.fasta//g;
	push @names, $genome;
	push @sites, "$genome\_binding_sites";
	
	$Genomes{$genome} = {};
	$BindingSites{$genome} = {};
}

my %Kmers = ();
my $x = 1;		# Just a counter for $v
# Operate on all the genomes for each value of k
foreach my $genome_file ( @proc_order )	{
	$manager->start and next;	

	my $genome = fileparse($genome_file, qr/\.[^.]*/);
	$genome =~ s/(\.fasta|\.[0-9].*)//g;
	print STDERR " ---  Processing $genome ( genome $x / ", scalar @proc_order, " ) --- \n" if ( $v == 1 );
	$x++;
	
	# Open and read the file
	$/ = ">";
	open FASTA, "$genome_file" if ( -e $genome_file );		# Might not have any excluded sequences;
		my @fastas = <FASTA>;
		my $trash = shift @fastas;	# Get rid of the first element, which will be a lone ">" symbol
	close FASTA;
	$/ = "\n";

	# Open the initial output files for this genome
	open BINDSITES, ">", "$genome.EXCLUDE.kmers.binding_sites.tab";
	
	# This one holds all the relevant metadata about the kmer
	open PROFILE, ">", "$genome.EXCLUDE.kmers.profiles.tab";
	print PROFILE "Kmer\tCount\tLength\tGC%\tA\tT\tC\tG\tN\tTm\tMinDist\tMaxDist\tMeanDist\tStdevDist\tGini\tGibbsFreeEnergy\n";
	
	$Genomes{$genome}->{total_len} = 0;		# Reset this for each genome
	foreach my $record ( @fastas )	{
		my ($header, @seq) = split "\n", $record;
		$header=~ s/ /_/g;				# Convert any spaces to underscores
		$header =~ s/\s+//g;			# Remove any other odd whitespace
		$header =~ s/\|//g;				# Remove pipe chars (often seen in Genbank headers)
		$header =~ s/>//g;
		
		# Skip this if we think we have an organelle, and the flag is activated
		if ( $organelles == 0 && ( $header =~ m/(plasmid|oplast|mitochond)/ ) )		{
			next;		# If we are skipping organelles and we come across one, skip.
		}
		elsif ( $organelles == 2 && ( $header !~ m/(plasmid|oplast|mitochond)/ ) )		{
			next;		# Skip if we are only interested in organelles and we dont think this contig is one
		}
		
		
		### NOTE: we are not storing the FASTA sequences in this program to save some memory space
		my $seq = join '', @seq;
		$seq =~ s/>//g;						# Remove any lingering ">" symbols
		$seq = uc $seq;						# Convert everything to uppercase 
	
		# Count the kmers of size(s) $k in the genomes
		foreach my $k ( @ks )	{
			print STDERR " ---  Checkpoint 1: Identifying $k-mers for $header  --- \n" if ( $v == 1 );
			for ( my $i=0; $i < length $seq; $i++ )	{
				last if ( length($seq) - $k < $i );
				
				# Get the kmer itself and some basic attributes 
				my $kmer = substr($seq, $i, $k );	
				my $revcom = revcom($kmer);

				# Increment the kmer if we've seen it before --> 
				# note that we are including the reverse complement of the same as 
				# the same one and only storing one of them (whichever is encountered first)
				if ( exists $Kmers{$kmer} && not exists $Kmers{$revcom}	) {
					$Kmers{$kmer}++;
					push @{ $BindingSites{$kmer} }, $i + $Genomes{$genome}->{total_len};
				}
				elsif ( not exists $Kmers{$kmer} && exists $Kmers{$revcom}	) {
					$Kmers{$revcom}++;
					push @{ $BindingSites{$revcom} }, $i + $Genomes{$genome}->{total_len};
				}
				# Otherwise, we assume this is the first time we've seen this kmer, so initialize it's count to 1.
				# Additionally, initialize the binding sites if needed and update the first instance
				else	{
					# Increment the current genome's kmer count and binding sites, where the current mapping site is the global length plus the current position
					$Kmers{$kmer} = 1;
					push @{ $BindingSites{$kmer} }, $i + $Genomes{$genome}->{total_len};							
				}
			}
		}
		
		# Update the global total length for the genome and commit changes to the database
		$Genomes{$genome}->{total_len} += length $seq;
	}

	print STDERR " ---  Checkpoint 2: Computing stats for $genome  --- \n" if ( $v == 1 );
	foreach my $kmer ( keys %Kmers )		{				

		# Grab some simple attributes
		my ( $a, $t, $c, $g, $n ) = base_composition($kmer);
		my $length = length $kmer;
		my $gc_content = gc_content($kmer);
		my $tm = Tm($kmer);
		
		# Some decisions to make for filtering at this point:
		#  -- less than minimum frequency
		#  -- Tm out of specified ranges
		#  -- Ambiguity count too high
		next if ( $Kmers{$kmer} < $min_freq || $tm > $max_tm || $tm < $min_tm || $n > $allow_ambig );
	
		# Do we expect this primer to dimerize with itself (a homodimer)?
		# This is estimated with Gibbs free energy, allowing up to -6 kcal/mol
		my $gibbs = gibbs_free_energy( $kmer );

		# Transform the input list of binding sites into a list of distances between each binding site and calculate the mean distance between binding sites
		#( we will use these distances for the final Gini calculation)
		my @distances;
		my $max_dist = "NA";
		my $min_dist = "NA";
		my $mean_dist = "NA";
		my $stdev_dist = "NA";
		for ( my $p=1; $p < scalar @{ $BindingSites{$kmer} }; $p++ )	{
			my $q = $p - 1;
			push @distances, @{$BindingSites{$kmer}}[$p] - @{$BindingSites{$kmer}}[$q]; 	# The current element minus the previous element
		}
		$max_dist = max( @distances );
		$min_dist = min( @distances );
		$mean_dist = int(mean( @distances ));	
		$stdev_dist = int(stdev( @distances ));

		# Calculate the individual Gini index of the binding sites in the genome (index of the even-ness of breadth coverage)
		# Using a minimum value of 3 sites here, since only 2 sites == 1 spacer, which is useless to calculate
		my $gini = ( scalar @distances > 2 )? gini_index( @distances ) : "NA";

		# Further filtering decisions to make:
		# -- mean binding distance is out of specified ranges
		# -- Gini index is too high
		next if ( ($gini ne "NA" && $gini > $max_gini) || ($mean_dist ne "NA" && ($mean_dist < $min_bind_dist || $mean_dist > $max_bind_dist)) );
		
		# If we made it this far with this kmer, then 
		# Write the data out to the two output files per genome
		print PROFILE "$kmer\t$Kmers{$kmer}\t$length\t$gc_content\t$a\t$t\t$c\t$g\t$n\t$tm\t$min_dist\t$max_dist\t$mean_dist\t$stdev_dist\t$gini\t$gibbs\n";
		print BINDSITES "$kmer\t", join(",", @{$BindingSites{$kmer}}), "\n";
	}
	
	# Reset %BindingSites and %Kmers for the next genome
	%BindingSites = ();
	%Kmers = ();
	
	close BINDSITES;
	close PROFILE;

	# Exit the child process
	$manager->finish;
}



exit;




#################### SUBROUTINES ###################################3

# Descriptive Statistics 
sub mean	{ 
	return sum(@_)/@_;
}

sub stdev	{
	my @data = @_;
	if ( @data == 1 ) {  return 0;    }
	my $mean = mean(@data);
	my $std = 0;
	foreach ( @data ) {	$std += ( $_ - $mean ) ** 2;  }
	$std = ( $std / scalar @data ) ** 0.5;
	return $std;
}

sub sum	{
	my ( @numbers ) = @_;
	my $sum = 0;
	foreach ( @numbers )	{	$sum += $_;		}
	return $sum;
}

# Returns the max value from a set of numerical args
sub max 	{
	my (@l) = @_;
	my $max = $l[0];
	foreach my $x ( @l )	{	$max = $x if ( $x > $max );	}
	return $max;
}

# Returns the min value from a set of numerical args
sub min	{
	my (@l) = @_;
	my $min = $l[0];
	foreach my $x ( @l )	{	$min = $x if ( $x < $min );	}
	return $min;
}

# Calculates the Gini coefficient from a 1D list of non-negative, ascending sorted numbers
# The index ranges from 0 to 1, with values approaching 1 indicating increasing unevenness (we want a value close to 0)
sub gini_index		{
	my ( @distances ) = @_;
	@distances = sort { $a <=> $b } @distances;
	my $n = scalar @distances;
	my $sum = sum( @distances );
	my $gini;
	for ( my $i=0; $i < $n; $i++ )	{
		$gini += ((2*($i  +1) - $n - 1) * $distances[$i]) / ($n * $sum);		# Function described: https://github.com/oliviaguest/gini
	}
	return sprintf("%1.4f", $gini);
}

# Count the number of ambiguous bases in a sequence
# Returns the ratio of the number of ambiguous bases divided by the length of the sequence
sub base_composition {
	my ( $primer ) = @_;
	my $Acount = ( $primer =~ tr/Aa// );
	my $Tcount = ( $primer =~ tr/Tt// );
	my $Gcount = ( $primer =~ tr/Gg// );
	my $Ccount = ( $primer =~ tr/Cc// );
	my $Ncount = ( $primer =~ tr/YyRrWwSsKkMmDdVvHhBbNnXx// );
	return ( $Acount, $Tcount, $Ccount, $Gcount, $Ncount );
}

# Reverse complements a nucleotide sequence, including ambiguous base calls
sub revcom	{
	my ( $seq ) = @_;
	my $revcom = reverse $seq;
	$revcom =~ tr/ACGTYRWSKMDVHBNX-/TGCARYWSMKHBDVNX-/ ;
	return ( uc $revcom );
}

# Returns only the complement of a nucleotide sequence, not the reverse
sub complement	{
	my ( $seq ) = @_;
	my $com = uc $seq;
	$com =~ tr/ACGTYRWSKMDVHBNX-/TGCARYWSMKHBDVNX-/;
	return $com;
}

# Simple method used by MHS to estimate the annealing temperature of a primer sequence
# Not going to fancy with this calculation, even though much more elaborate models exist for estimating Tm.
sub Tm		{
	my ( $primer ) = @_;
	my $at = $primer =~ tr/ATYRWKMDVHBNX//;
	my $gc = $primer =~ tr/GCS//;
	return ( (2*$at) + (4*$gc) );		# Tm = 2(at) + 4(gc), here allowing for ambiguities
}

# Computes the GC% of a putative primer and returns 1 if the GC content exceeds 80% of the primer sequence
sub gc_content		{
	my ( $primer ) = @_;
	my $gc = ( $primer =~ tr/GCS// );
	return sprintf("%3.2f", $gc / length $primer );
}

# Calculates the Gibbs free energy (G) for a primer, related to secondary structure
# G is essentially a measure of the spontaneity of a reaction, in this case
# the energy required to break the secondary structure of a hairpin loop in a primer.
# Higher negative values of G indicate stable, undesirable hairpins.
# Formula:
# 	deltaG = deltaH - T*deltaS, where H is the Enthalpy and S is the Entropy
# The deltaG here is calculated by taking into account the longest stretch of complementary bases
sub gibbs_free_energy	{
	my ( $primer ) = @_;
	my $G = 0;
	my $T = 303.15; 				# T = the assumed Ta of the rxn in Celsuis, converted to Kelvin by adding 273.15 (thus, our T here is 30C )
	
	# Hash of nearest-neighbor dinucleotide thermodynamics
	# The units here are in kcal/mol (S is cal/mol*K), with the T = 298 Kelvin (just below our ideal temp) and using 1M MgCl as the salt
	# The entries "G25mM" are the free energy after salt correction for 25 mM MgCl2 that we use in Molecular Epi lab ( G25mM = G - ( mi * log(0.025)) )
	# This data is empirically derived and described in Table 2 of this publication:
	# 		Huguet, J. M., Ribezzi-Crivellari, M., Bizarro, C. V., & Ritort, F. (2017). Derivation of nearest-neighbor DNA parameters in 
	#		magnesium from single molecule experiments. Nucleic Acids Research, 45(22), 12921.
	my %NNBP = (
		"AATT"   => { "H" => -8.46, "S" => -22.7, "G" => -1.69, "mi" => 0.086, "G25mM" => -1.5522  },
		"ACTG"   => { "H" => -8.67, "S" => -22.7, "G" => -1.91, "mi" => 0.074, "G25mM" => -1.7914  },
		"AGTC"   => { "H" => -7.65, "S" => -19.6, "G" => -1.81, "mi" => 0.070, "G25mM" => -1.1214  },
		"ATTA"   => { "H" => -8.62, "S" => -23.7, "G" => -1.55, "mi" => 0.093, "G25mM" => -1.4010  },
		"CAGT"   => { "H" => -9.99, "S" => -26.2, "G" => -2.17, "mi" => 0.079, "G25mM" => -2.0434  },
		"CCTT"   => { "H" => -9.90, "S" => -25.9, "G" => -2.18, "mi" => 0.032, "G25mM" => -2.1287  },
		"CGGC"   => { "H" => -11.03,"S" => -28.1, "G" => -2.65, "mi" => 0.058, "G25mM" => -2.5571  },
		"CCGG"   => { "H" => -11.03,"S" => -28.1, "G" => -2.65, "mi" => 0.058, "G25mM" => -2.5571  },			# Note: HS added this one, it is the same as the CGGC line
		"GACT"   => { "H" => -8.59, "S" => -22.5, "G" => -1.88, "mi" => 0.075, "G25mM" => -1.7598  },
		"GCCG"   => { "H" => -11.22,"S" => -28.4, "G" => -2.74, "mi" => 0.058, "G25mM" => -2.6471  },
		"TAAT"   => { "H" => -7.49, "S" => -20.4, "G" => -1.38, "mi" => 0.088, "G25mM" => -1.2390  },
		"ATinit" => { "H" => 2.77, "S" => 6.2, "G" => 0.91, "mi" => 0 },
		"GCinit" => { "H" => 2.65, "S" => 5.6, "G" => 0.97, "mi" => 0 },
	);
	
	# Gibbs free energy is the sum of the individual (salt corrected) NNBP free energies ( here, just the sum of the G25mM terms )
	for ( my $i=0; $i < length $primer; $i++ )	{
		my $nearneighbor = substr($primer, $i, 2);
		last if ( length $nearneighbor == 1 );
		$nearneighbor .= complement($nearneighbor);
		# Add the deltaG if we have a match.
		if ( exists $NNBP{$nearneighbor}{G25mM} )	{
			$G += $NNBP{$nearneighbor}{G25mM};		
		}
		else	{
			$G += -1.8241;		# If we dont have any matches, the add the average G25mM as our best guess estimate.
		}
	}
	return $G;
}

# Dynamically computes a semi-global alignment between two primers using Watson-Crick complementary base pairing.
# Both primers are assumed to be in 5' -> 3' orientation
# Returns the predicted Gibbs free energy of the longest stretch of complementary base pairings
sub dimer		{
	my ( $query, $subject, $semiglobal )	= @_;
	$semiglobal = ( $semiglobal )? $semiglobal: 0;
	my @matrix;
	
	# Our scoring function
	my $match = 1;
	my $mismatch = -1;
	my $gapopen = -99;
	my $gapextend = -5;
	my $qry_end_gap = ( $semiglobal == 0 )? -99 : 0;
	my $subj_end_gap = ( $semiglobal == 0 )? -99 : 0;
	my $qry_initial_gap = ( $semiglobal == 0 )? -99 : 0;
	my $subj_initial_gap = ( $semiglobal == 0 )? -99 : 0;
	
	# Set some flags and initialize the list of possible longest subsequences
	my $gap_flag = 0;
	my $longest_flag = 1;
	my $longest_substring = "";
	my $longest_gc = 0.0;		# We will keep track of the GC content of the longest subsequence, since we want to prioritize CG pairings over others.
	
	# Watson-Crick canonical base pairings, including ambiguity codes
	my %Complement = (
		"A" => "T", "C" => "G",	"G" => "C",	"T" => "A",
		"Y" => "R", "R" => "Y", "W" => "W", "S" => "S",
		"K" => "M", "M" => "K", "D" => "H", "V" => "B",
		"H" => "D", "B" => "V", "X" => "N", "N" => "N",
	);

	# Initialize the matrix
	$matrix[0][0]{score} = 0;
	$matrix[0][0]{pointer} = "none";

	# Initializing the first column and row to 0 enables us to compute semi-global alignments if we want
	# (there is no gap penalty at the beginning of the sequences)
	for ( my $j = 1; $j <= length($query); $j++ )	{
		$matrix[0][$j]{score} = $qry_initial_gap * $j;
		$matrix[0][$j]{pointer} = "left";
	}
	for ( my $i = 1; $i <= length($subject); $i++ )	{
		$matrix[$i][0]{score} = $subj_initial_gap * $i;			
		$matrix[$i][0]{pointer} = "up";
	}

	# Fill the matrix
	my $tmp_longest; my $tmp_gc;
	for ( my $i = 1; $i <= length($subject) ; $i++ ) 	{
		for ( my $j = 1; $j <= length($query) ; $j++ ) 	{
			my ( $diagonal_score, $left_score, $up_score );
			
			# Calculate the match score
			my $letter1 = substr($query, $j-1, 1);
			my $letter2 = substr($subject, $i-1, 1);
			
			# If we have a Watson-Crick match, count the match and update the longest substring and GC content
			if ( $letter1 eq $Complement{$letter2} )	{		
				$diagonal_score = $matrix[$i-1][$j-1]{score} + $match;
				
				# Update the longest substring data and the flags
				$longest_flag = 1;
				$tmp_longest = ( $longest_flag == 0 )? "" : "$tmp_longest" . "$letter1";
				$tmp_gc = ( length $tmp_longest == 0 )? 0.0 : gc_content($tmp_longest);
				if ( length $tmp_longest >= length $longest_substring && $tmp_gc >= $longest_gc )	{
					$longest_substring = $tmp_longest;
					$longest_gc = $tmp_gc;
				}
				
			}
			else	{		
				$diagonal_score = $matrix[$i-1][$j-1]{score} + $mismatch;	
				$longest_flag = 0;			# Turn the flag off for the next comparison
			}
	
			# Calculate the gap scores, including the terminal gaps
			if ( $i == length($subject) )	{
				$up_score = $matrix[$i-1][$j]{score};
				$left_score = $matrix[$i][$j-1]{score} + $qry_end_gap;
			}
			elsif ( $j == length($query) )	{
				$up_score = $matrix[$i-1][$j]{score} + $subj_end_gap;
				$left_score = $matrix[$i][$j-1]{score};
			}
			else	{
				$up_score   = ( $gap_flag == 1 )? $matrix[$i-1][$j]{score} + $gapextend : $matrix[$i-1][$j]{score} + $gapopen;
				$left_score = ( $gap_flag == 1 )? $matrix[$i][$j-1]{score} + $gapextend : $matrix[$i][$j-1]{score} + $gapopen;
			}

			# Choose the best score
			if ( $diagonal_score >= $up_score )	{
				if ( $diagonal_score >= $left_score ) 	{
					$matrix[$i][$j]{score} = $diagonal_score;
					$matrix[$i][$j]{pointer} = "diagonal";
					$gap_flag = 0;
				}
				else {
					$matrix[$i][$j]{score} = $left_score;
					$matrix[$i][$j]{pointer} = "left";
					$gap_flag = 1;
				}
			}
			else 	{
				if ( $up_score >= $left_score ) 	{
					$matrix[$i][$j]{score} = $up_score;
					$matrix[$i][$j]{pointer} = "up";
					$gap_flag = 1;
				}
				else	{
					$matrix[$i][$j]{score} = $left_score;
					$matrix[$i][$j]{pointer} = "left";
					$gap_flag = 1;
				}
			}	
		}
	}

	# Calculate and return the expected Gibbs free energy of this subsequence, if the length is 2 or more.  ( Just return 0 if there is only 1 base pairing).
	( length $longest_substring  >= 2 )? return sprintf("%3.2f", gibbs_free_energy($longest_substring )) : return "0.000";
}
