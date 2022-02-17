#!/usr/bin/perl

# kmers.pl
# A naive kmer-counting wrapper that constructs a SQLite database of kmers and their binding sites within genomes.
# Be warned that this no effort has been made here to make this space or memory efficient, 
# so expect this to eat up alot of system resources that scales with the size of the genomes

# Author: MH Seabolt
# Last updated: 3-6-2020

use strict;
use warnings;
use Getopt::Long qw(GetOptions);
use Cwd qw(cwd);
use List::Util qw(first shuffle uniq);
use File::Basename; 

# Required CPAN modules
use DBI;

# Required input parameters
my $input_list;
my $genome_list;
my $scoring;
my $output;
my $mink;
my $maxk;
my $step;
my $v;
my $dbname;
my $journalmode;

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
	-db 			name of the database file you want to access or create, requires DBI module from CPAN and SQLite to be installed
           
	=== KMER PARAMETERS ===
	-mink                    INT; the minimum value of k [Default: 8]
	-maxk                    INT; the maximum value of k [Default: 12]
	-step                    INT; the step size to use for comparing values of k. Default = 1
	
	\n";
	print $usage;
}
if ( scalar @ARGV == 0 ) { die usage(); }

GetOptions( 'in|i=s' => \$genome_list,
			'inlist|il=s' => \$input_list,
            'mink=i' => \$mink,
            'maxk=i' => \$maxk,
            'step=i' => \$step,
            'verbose|v' => \$v,
			'output|o=s' => \$output,
			'score=i' => \$scoring,
			'dbname=s' => \$dbname,
			'journal=s' => \$journalmode,
			
) or die usage(); 

# Parameter setups
$mink = ( $mink && int($mink) >= 5 )? $mink : 8;
$maxk = ( $maxk && int($maxk) > int($mink) && int($maxk) <= 30 )? $maxk : $mink;
$step = ( $step && int($step) > 1 )? $step : 1;
$scoring = ( $scoring && $scoring == 1 )? 1 : 0;
$v = defined($v)? 1 : 0;
$dbname = ( $dbname )? $dbname : $output;
$journalmode = ( $journalmode )? $journalmode : "WAL";

##############################################
#       PHASE 0: INITIAL SETUP
##############################################

# Hook up our SQlite database handle
my $dbh = DBI->connect("dbi:SQLite:dbname=$dbname", "", "", { 
	PrintError => 0,
	RaiseError => 1,
	AutoCommit => 0,
} ) or die "ERROR: Cannot connect to the requested database!\n $DBI::errstr\n";
$dbh->do( "PRAGMA journal_mode = $journalmode" );

# Generate the database structure if one does not already exist
$dbh->do( "CREATE TABLE IF NOT EXISTS kmers ( 
	Sequence CHAR(25) NOT NULL,
	Length INT NOT NULL,
	GC_percent decimal(3,2) NOT NULL,
	A INT NOT NULL,
	T INT NOT NULL,
	C INT NOT NULL,
	G INT NOT NULL,
	N INT NOT NULL,
	Tm decimal(3,1) NOT NULL,
	PRIMARY KEY (Sequence)
)" );

$dbh->do( "CREATE TABLE IF NOT EXISTS genomes (
	Name CHAR(25) NOT NULL,
	File CHAR(50) NOT NULL,
	PRIMARY KEY ( Name )
)" );
	
$dbh->do( "CREATE TABLE IF NOT EXISTS contigs (
	Name CHAR(25) NOT NULL,
	Genome CHAR(25) NOT NULL,
	PRIMARY KEY ( Name, Genome ),
	CONSTRAINT contigs_ibfk_1 FOREIGN KEY ( Genome ) REFERENCES genomes ( Name )
)" );	
	
$dbh->do( "CREATE TABLE IF NOT EXISTS counts (
	Kmer CHAR(25) NOT NULL,
	Genome CHAR(25) NOT NULL,
	Count INT NOT NULL,
	PRIMARY KEY ( Kmer, Genome ),
	CONSTRAINT counts_ibfk_1 FOREIGN KEY ( Kmer ) REFERENCES kmers (Sequence),
	CONSTRAINT counts_ibfk_2 FOREIGN KEY ( Genome ) REFERENCES genomes ( Name )
)" );

$dbh->do( "CREATE TABLE IF NOT EXISTS locations (
	Kmer CHAR(25) NOT NULL,
	Genome CHAR(25) NOT NULL,
	Contig CHAR(25) NOT NULL,
	Position INT NOT NULL,
	GlobalIndex INT NOT NULL,
	PRIMARY KEY ( Kmer, Genome, Contig, Position ),
	CONSTRAINT locations_ibfk_1 FOREIGN KEY ( Genome ) REFERENCES genomes ( Name ),
	CONSTRAINT locations_ibfk_2 FOREIGN KEY ( Contig ) REFERENCES contigs ( Name ),
	CONSTRAINT locations_ibfk_3 FOREIGN KEY ( Kmer ) REFERENCES kmers ( Sequence )
)" );

$dbh->do( "CREATE TABLE IF NOT EXISTS stats (
	Kmer CHAR(25) NOT NULL,
	Genome CHAR(25) NOT NULL,
	MinDist INT,
	MaxDist INT,
	MeanDist INT,
	StdevDist INT,
	Gini DECIMAL(1,4),
	Gibbs DECIMAL(9,3),
	PRIMARY KEY(Kmer, Genome),
	CONSTRAINT stats_ibfk_1 FOREIGN KEY ( Kmer ) REFERENCES kmers ( Sequence ),
	CONSTRAINT stats_ibfk_2 FOREIGN KEY ( Genome ) REFERENCES genomes ( Name )
)" );

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
print STDERR " ===        PHASE 1: KMER COUNTS        === \n" if ( $v == 1 );

# Initialize the multidimensional hashes and some simple lists
my %Genomes = ();
my %Primers = ();
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
	$dbh->do( "INSERT INTO genomes ( Name, File ) VALUES ( '$genome', '$genome_file' )" );
	
	$Genomes{$genome} = {};
	$BindingSites{$genome} = {};
	foreach my $k ( @ks )	{
		$Primers{$genome}->{$k} = {};
	}
}

my %Kmers = ();
# Operate on all the genomes for each value of k
foreach my $genome_file ( @proc_order )	{

	my $genome = fileparse($genome_file);
	$genome =~ s/\.fasta//g;
	print STDERR " ---  Processing $genome genome  --- \n" if ( $v == 1 );
	
	# Open and read the file
	$/ = ">";
	open FASTA, "$genome_file" if ( -e $genome_file );		# Might not have any excluded sequences;
		my @fastas = <FASTA>;
		my $trash = shift @fastas;	# Get rid of the first element, which will be a lone ">" symbol
	close FASTA;
	$/ = "\n";

	$Genomes{$genome}->{total_len} = 0;		# Reset this for each genome
	foreach my $record ( @fastas )	{
		my ($header, @seq) = split "\n", $record;
		$header=~ s/ /_/g;				# Convert any spaces to underscores
		$header =~ s/\s+//g;			# Remove any other odd whitespace
		$header =~ s/\|//g;				# Remove pipe chars (often seen in Genbank headers)
		$header =~ s/>//g;
		#$header =~ s/\..*//g;
		$dbh->do( "INSERT INTO contigs ( Name, Genome ) VALUES ( '$header', '$genome' )" );
		
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

				# Increment the kmer if we've seen it before.
				if ( exists $Kmers{$kmer} )		{
					push @{ $BindingSites{$genome}->{$kmer} }, $i + $Genomes{$genome}->{total_len};				
					$dbh->do( "UPDATE counts SET Count = (Count + 1) WHERE ( Kmer = '$kmer' AND Genome = '$genome' ) " );
					$dbh->do( "INSERT INTO locations ( Kmer, Genome, Contig, Position, GlobalIndex ) VALUES ( '$kmer', '$genome', '$header', '$i', $i+$Genomes{$genome}->{total_len})" );
				}
				# Otherwise, we assume this is the first time we've seen this kmer, so initialize it's count to 1.
				# Additionally, initialize the binding sites if needed and update the first instance
				else	{
					# Initialize the kmer in all genomes to 0 if it's never been seen before
					foreach ( @names )	{
						$dbh->do( "INSERT INTO counts ( Kmer, Genome, Count ) VALUES ( '$kmer' , '$_', '0' )" ) if ( not exists $Kmers{$kmer}{$_} );		#$Primers{$_}->{$k}->{$kmer} );
						$Kmers{$kmer}{$_} = 1 if ( not exists $Kmers{$kmer}{$_} );
					}
					
					# Increment the current genome's kmer count and binding sites
					$dbh->do( "UPDATE counts SET Count = '1' WHERE Kmer = '$kmer' AND Genome = '$genome' " );
					push @{ $BindingSites{$genome}->{$kmer} }, $i + $Genomes{$genome}->{total_len};			# The current mapping site is the global length plus the current position
					
					# Insert the kmer and some simple attributes into the SQLite database
					my ( $a, $t, $c, $g, $n ) = base_composition($kmer);
					my $gc_content = gc_content($kmer);
					my $tm = Tm($kmer);
					$dbh->do( "INSERT INTO kmers ( Sequence, Length, GC_percent, A, T, C, G, N, Tm ) VALUES ( '$kmer', '$k', '$gc_content', '$a', '$t', '$c', '$g', '$n', '$tm' )" );	
					$dbh->do( "INSERT INTO locations ( Kmer, Genome, Contig, Position, GlobalIndex ) VALUES ( '$kmer', '$genome', '$header', '$i', $i+$Genomes{$genome}->{total_len})" );
				}
			}
		}
		
		# Update the global total length for the genome and commit changes to the database
		$Genomes{$genome}->{total_len} += length $seq;
		$dbh->commit;
	}

	print STDERR " ---  Checkpoint 2: Computing stats for $genome  --- \n" if ( $v == 1 );
	foreach my $kmer ( keys %Kmers )		{				
		next if ( not exists $BindingSites{$genome}->{$kmer} );
		
		# Do we expect this primer to dimerize with itself (a homodimer)?
		# This is estimated with Gibbs free energy, allowing up to -6 kcal/mol
		
		
		# Do we expect internal or 3' hairpin secondary structures to form?
		# Estimate this with Gibbs free energy computation -- we will tolerate 3' end deltaG of -2 kcal/mol or an internal hairpin of -3 kcal/mol
		
		
		# Transform the input list of binding sites into a list of distances between each binding site and calculate the mean distance between binding sites
		#( we will use these distances for the final Gini calculation)
		my @distances;
		my $max_dist = "NA";
		my $min_dist = "NA";
		my $mean_dist = "NA";
		my $stdev_dist = "NA";
		@{$BindingSites{$genome}->{$kmer}} = sort { $a <=> $b } @{$BindingSites{$genome}->{$kmer}};
		if ( scalar @{$BindingSites{$genome}->{$kmer}} > 2 ) 	{
			for ( my $p=1; $p < scalar @{ $BindingSites{$genome}->{$kmer} }; $p++ )	{
				my $q = $p - 1;
				push @distances, @{$BindingSites{$genome}->{$kmer}}[$p] - @{$BindingSites{$genome}->{$kmer}}[$q]; 	# The current element minus the previous element
			}
			@distances = sort { $a <=> $b } @distances;
			$max_dist = max( @distances );
			$min_dist = min( @distances );
			$mean_dist = int(mean( @distances ));	
			$stdev_dist = int(stdev( @distances ));
		}

		# Calculate the individual Gini index of the binding sites in the foreground genome (index of the even-ness of breadth coverage)
		my $gini = ( scalar @distances > 2 )? gini_index( @distances ) : "NA";
			
		# Calculate the Gibbs free energy for the kmer
		# --- UNDER CONSTRUCTION
		my $gibbs = gibbs_free_energy( $kmer );
			
		# Write the data out to the database
		$dbh->do( "INSERT INTO stats ( Kmer, Genome, MinDist, MaxDist, MeanDist, StdevDist, Gini, Gibbs ) VALUES ( '$kmer', '$genome', '$min_dist', '$max_dist', '$mean_dist', '$stdev_dist', '$gini', '$gibbs' )" );
	
		# Delete the list in the %BindingSites hash to help clean up RAM
		delete $BindingSites{$genome}->{$kmer};
	}
	$dbh->commit;	
}

# Commit changes and disconnect from the database
$dbh->disconnect or warn "ERROR: Database disconnection failed!\n $DBI::errstr\n";


exit;




#################### SUBROUTINES ###################################3

# Descriptive Statistics 
sub mean	{ 
	return sum(@_)/@_;
}

sub stdev	{
	my @data = @_;
        if ( @data == 1 ) {
                return 0;
        }
        my $mean = mean(@data);
        my $std = 0;
        foreach ( @data ) {
                $std += ( $_ - $mean ) ** 2;
        }
        $std = ( $std / scalar @data ) ** 0.5;
        return $std;
}

sub sum	{
	my ( @numbers ) = @_;
	my $sum = 0;
	
	foreach ( @numbers )	{
		$sum += $_;
	}
	
	return $sum;
}

# Returns the max value from a set of numerical args
sub max 	{
	my (@l) = @_;
	my $max = $l[0];
	foreach my $x ( @l )	{
		$max = $x if ( $x > $max );
	}
	return $max;
}

# Returns the min value from a set of numerical args
sub min	{
	my (@l) = @_;
	my $min = $l[0];
	foreach my $x ( @l )	{
		$min = $x if ( $x < $min );
	}
	return $min;
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
		$G += $NNBP{$nearneighbor}{G25mM};
	}
	return $G;
}
