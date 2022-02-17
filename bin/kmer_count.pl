#!/usr/bin/perl

# kmer_count.pl
# A naive kmer-counting wrapper that counts kmers awithin genomes.
# Be warned that this no effort has been made here to make this space or memory efficient, 
# so expect this to eat up alot of system resources that scales with the size of the genomes.
# This script only COUNTS frequencies and has minimal filtering options, use kmer_profile.pl to compute a more detailed profile of a genome, including binding sites and extensive filtering options.

# Author: MH Seabolt
# Last updated: 5-15-2020

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
my $max_freq;
my $organelles;
my $v;
my $threads;

sub usage {
	my $usage = "kmer_count.pl\n
	PURPOSE:            A naive kmer-counting wrapper that outputs TAB files of kmer frequencies in a genome.
                        \n
	USAGE: kmer_count.pl [options]
          
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
	-min_freq			INT; discard a kmer if it occurs with less than this frequency.  
	-max_freq			INT; discard a kmer if it occurs with more than this frequency.  
	-organelles			INT flag; How to handle contigs in the genomes that are identified as organelles (mitochondria, plasmids, apicoplasts, chloroplasts )?
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
			'max_freq=i' => \$max_freq,
            'verbose|v' => \$v,
			'outdir|o=s' => \$outdir,
		  'organelles=i' => \$organelles,
		  'threads|t=i' => \$threads,
			
) or die usage(); 

# Parameter setups
$mink = ( $mink && int($mink) >= 5 )? $mink : 8;
$maxk = ( $maxk && int($maxk) > int($mink) && int($maxk) <= 30 )? $maxk : $mink;
$step = ( $step && int($step) > 1 )? $step : 1;
$organelles = ( $organelles && $organelles =~ /[0-2]/ )? $organelles : 1;
$min_freq = ( $min_freq && $min_freq > -1 )? $min_freq : 5;
$max_freq = ( $max_freq && $max_freq > $min_freq )? $max_freq : 'inf';
$v = defined($v)? 1 : 0;
$threads = ( $threads && $threads >= 1 )? $threads : 1;
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
my @names;
foreach my $genome_file ( @proc_order )	{
	# Get the genome's basename
	chomp $genome_file;
	my $genome = fileparse($genome_file);
	$genome =~ s/\.fasta//g;
	push @names, $genome;
	$Genomes{$genome} = {};
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
	
	# This one holds all the relevant metadata about the kmer
	open KMERS, ">", "$outdir\/$genome.kmers.count.tab";
	print KMERS "Kmer\tCount\n";
	
	
	# Count the kmers of size(s) $k in the genomes
	foreach my $k ( @ks )	{
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
	
			print STDERR " ---  Checkpoint 1: Identifying $k-mers for $header  --- \n" if ( $v == 1 );
	
			for ( my $i=0; $i < length $seq; $i++ )	{
				last if ( length($seq) - $k < $i );
				
				# Get the kmer itself and some basic attributes 
				my $kmer = substr($seq, $i, $k );	
				$kmer = canonical( $kmer, revcom($kmer) );
				
				# Increment the kmer if we've seen it before --> 
				# note that we are including the reverse complement of the same as 
				# the same one and only storing one of them (whichever is encountered first)
				if ( exists $Kmers{$kmer} ) {
					$Kmers{$kmer}++;
				}
				# Otherwise, we assume this is the first time we've seen this kmer, so initialize it's count to 1.
				# Additionally, initialize the binding sites if needed and update the first instance
				else	{
					# Increment the current genome's kmer count and binding sites, where the current mapping site is the global length plus the current position
					$Kmers{$kmer} = 1;					
				}
			}
		}
		# Check that the kmer falls within our min and max frequency parameters, then print
		foreach my $kmer ( keys %Kmers )		{				
			next if ( $Kmers{$kmer} < $min_freq || $Kmers{$kmer} > $max_freq );
			print KMERS "$kmer\t$Kmers{$kmer}\n";
		}
		
		# Reset %Kmers for the next genome
		%Kmers = ();
	}
	close KMERS;
	
	# End the child process
	$manager->finish if ( $threads > 1 );
}

exit;




#################### SUBROUTINES ###################################3

# Reverse complements a nucleotide sequence, including ambiguous base calls
sub revcom	{
	my ( $seq ) = @_;
	my $revcom = reverse $seq;
	$revcom =~ tr/ACGTYRWSKMDVHBNX-/TGCARYWSMKHBDVNX-/ ;
	return ( uc $revcom );
}

# Returns the lexicographically smaller (ie. first in alphabetical order) of a kmer and it's reverse complement.
sub canonical	{
	return ( sort { $a cmp $b } @_ )[0];
}