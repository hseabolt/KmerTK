#!/usr/bin/perl

# kmers.pl
# A naive kmer-counting wrapper that counts and computes some basic metrics of kmers and their binding sites within genomes.
# Be warned that this no effort has been made here to make this space or memory efficient, 
# so expect this to eat up alot of system resources that scales with the size of the genomes.
# This script has options to filter kmers based on several parameters, although it is intended to be used as a more 
# exhaustive search space, followed by more rigorous filtering using kmer_combiner's options for filtering.

# Several plug-in scripts have been written to help facilitate further operation on kmers:
#  -- kmer_count.pl 	 ==>  a very simplified version of this script that only counts the kmer frequencies, does not compute profiles.
#  -- kmer_combiner.pl   ==>  combines the frequency counts of a given list of genomes into a nice (but very large...) TAB matrix. Includes options for filtering kmers.
#  -- kmer_unique.pl     ==>  determines if a kmer is unique to a subset of genomes
#  -- kmer_condense.pl   ==>  collapses columns in a data matrix from kmer_combiner into single columns, given a config file
#  -- kmer_sets.pl       ==>  a wrapper script for doing Set operations given a data matrix
#  -- kmer_fetch.pl      ==>  fetches columns and/or row slices from an input data matrix
#  -- kmer_graph.pl      ==>  construct a graph network of thermodynamically-compatable primers, given a kmer profile computed by this script. Outputs either SIF or BindingSites formats.
#  -- kmer_gini.pl 		 ==>  Constructs sets of kmers to minimize the Gini idnex given a set of parameters.
# PLUGINS IN PROGRESS:
# *-- kmer_statistics.pl ==> computes scoring and statistics based on the detailed profile generated by kmers.pl.
# *-- kmer_merge.pl		 ==> merges two matrices into one bigger one, can be done with either a combined matrix, binding sites, or profiles
# *-- kmer_search.pl	 ==> uses a Trie data structure to search large numbers of kmers against multiple targets in a fast, parallelized manner.


# Author: MH Seabolt
# Last updated: 4-10-2020

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
my $genome_list;
my $outdir;
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
my $oligo_conc;
my $mg_conc;
my $monovalent_cation_conc;
my $dntp_conc;
my $temperature;
my $advanced_filters;

sub usage {
	my $usage = "kmers.pl\n
	PURPOSE:            A naive kmer-counting wrapper that outputs TAB files of kmers and their binding sites in a series of genomes.
                        \n
	USAGE: kmers.pl [options]
          
	=== INPUT PARAMETERS ===
	-i                       input CSV list of FASTA formatted genome files
	-il                      input file containing a list of genome files, one per line   

	=== OUTPUT PARAMETERS ===
	-outdir                  directory name for writing output files
           
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
	-filter_repeats			int; use advanced filtering to catch kmers with runs of length 4 or greater, dimer repeats (4 or more), and triplet repeats (3 or more)? Default: OFF
	-organelles				INT flag; How to handle contigs in the genomes that are identified as organelles (mitochondria, plasmids, apicoplasts, chloroplasts )?
								0 == SKIP organelles,  1 == count kmers on ALL contigs [Default],  2 == count kmers ONLY on organelles

	=== OPTIONAL EXTRAS ===
	-v 					verbose
	-t 					INT; number of threads to use for multi-threading
	
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
			'outdir|o=s' => \$outdir,
            'min_tm=i' => \$min_tm,
            'max_tm=i' => \$max_tm,
            'max_gini=f' => \$max_gini,
		  'min_bind_dist=i' => \$min_bind_dist,
		  'max_bind_dist=i' => \$max_bind_dist,
		  'allow_ambig=i' => \$allow_ambig,
		  'organelles=i' => \$organelles,
		  'threads|t=i' => \$threads,
		  	'temperature=f' => \$temperature,
			'mgcl=f' => \$mg_conc,
			'oligo_conc=f' => \$oligo_conc,
			'dntps=f' => \$dntp_conc,
			'monovalent=f' => \$monovalent_cation_conc,
			'filter_repeats=i' => \$advanced_filters,
			
) or die usage(); 

# Parameter setups
$mink = ( $mink && int($mink) >= 5 )? $mink : 8;
$maxk = ( $maxk && int($maxk) > int($mink) )? $maxk : $mink;
$step = ( $step && int($step) > 1 )? $step : 1;
$organelles = ( $organelles && $organelles =~ /[0-2]/ )? $organelles : 1;
$min_freq = ( $min_freq && $min_freq > -1 )? $min_freq : 5;
$min_tm = ( $min_tm && $min_tm >= 15 )? $min_tm : 15;
$max_tm = ( $max_tm && $max_tm <= 80 )? $max_tm : 45;		
$min_bind_dist = ( $min_bind_dist && $min_bind_dist >= 5000 )? $min_bind_dist : 5000;
$max_bind_dist = ( $max_bind_dist && $max_bind_dist <= 1500000 )? $max_bind_dist : 1500000;		# 1,500,000
$allow_ambig = ( $allow_ambig && $allow_ambig <= $maxk )? $allow_ambig : 0;
$max_gini = ( $max_gini && $max_gini <= 1.0 && $max_gini >= 0.0 )? $max_gini : 1.0;
$v = defined($v)? 1 : 0;
$threads = ( $threads && $threads >= 1 )? $threads : 1;
$advanced_filters = ( $advanced_filters && $advanced_filters != 0 )? 1 : 0;
$outdir = ( $outdir )? $outdir : cwd();

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

mkdir $outdir if ( not -d $outdir );

############ THERMODYNAMICS SETUP ##################

# Some variables for primer dimer thermodynamics calculations --> the user can edit these, but they are currently directly set for MEPI lab@CDC standard protocols
# Starting ionic concentration variables
$oligo_conc = ( $oligo_conc )? $oligo_conc : 250; #in nM
$mg_conc = ( $mg_conc )? $mg_conc : 2.0; #in mM
$monovalent_cation_conc = ( $monovalent_cation_conc )? $monovalent_cation_conc : 50; #in mM
$dntp_conc = ( $dntp_conc )? $dntp_conc : 0.1; #in mM
$temperature = ( $temperature )? $temperature : 30;

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
	$manager->start and next if ( $threads > 1 );

	# The actual work we want to do
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
	open BINDSITES, ">", "$outdir\/$genome.kmers.binding_sites.tab";
	
	# This one holds all the relevant metadata about the kmer
	open PROFILE, ">", "$outdir\/$genome.kmers.profiles.tab";
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
				$kmer = canonical( $kmer, revcom($kmer) );
				
				# If we are using the repeat filters:
				next if ( $advanced_filters == 1 && ( hasRun($kmer) == 1 || hasDimerRepeat($kmer) == 1 || hasTripletRepeat($kmer) == 1 ) );	
				
				# Increment the kmer if we've seen it before --> 
				# note that we are including the reverse complement of the same as 
				# the same one and only storing one of them (whichever is encountered first)
				if ( exists $Kmers{$kmer} ) {
					$Kmers{$kmer}++;
					push @{ $BindingSites{$kmer} }, $i + $Genomes{$genome}->{total_len};
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
		my $gibbs = primer_dimer( $kmer, $kmer );

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
	
	# End the child process
	$manager->finish if ( $threads > 1 );
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

# Returns the lexicographically smaller (ie. first in alphabetical order) of a kmer and it's reverse complement.
sub canonical	{
	return ( sort { $a cmp $b } @_ )[0];
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

################################################################################################################################################################################# 
#
# The following code and subroutines are modified from "PerlPrimer", a Perl app written by Owen Marshall and licensed under the GNU Public License.
# Citation: Marshall OJ. PerlPrimer: cross-platform, graphical primer design for standard, bisulphite and real-time PCR. Bioinformatics 2004 20(15):2471-2472
#
#################################################################################################################################################################################

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