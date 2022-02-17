#!/usr/bin/perl

# hash_kmers.pl
# Accepts a genome sequence in FASTA format and generates a hash of kmer sequences of user-defined size.

use strict;
use warnings;
use Getopt::Long qw(GetOptions);
use Cwd qw(cwd);
use List::Util qw(first);
use File::Basename;
use POSIX qw(floor);

# Required input parameters
my $input = "--";
my $output = "--";
my $k;
my $mink;
my $maxk;
my $step;

sub usage {
	my $usage = "hash_kmers.pl\n
	PURPOSE: 	Accepts a genome sequence in FASTA format and generates a hash of kmer sequences of user-defined size.				
				\n
	USAGE:	hash_kmers.pl -i genome.fasta -d dir -o output -r reads.txt -v
	-i		input genome sequence file in FASTA format
	-o		output file name (not including extensions!)
	-k		(INT) kmer size to report (Default: 30)
	-mink	(INT) starting kmer size, if using steps.  Dont use -k if using steps.
	-maxk	(INT) maximum kmer size, if using steps.  Dont use -k if using steps.
	-step	(INT) step size, if using.
	\n";
	print $usage;
}
if ( scalar @ARGV == 0 )	{ die usage(); }

GetOptions(	'input|i=s' => \$input,
			'out|o=s' => \$output,
			'k=i' => \$k,
			'mink=i' => \$mink,
			'maxk=i' => \$maxk,
			'step|s=i' => \$step,
) or die usage();

# Parameter setups
die "INPUT FILE ERROR: $!\n" if ( not -e $input );

### Set the kmer values
if ( $k && $k > 0 )		{
	$mink = $k;
	$maxk = $k;
	$step = 0;
}
elsif ( not $k && ( ($mink && $mink > 0) || ($maxk && $maxk > 0) || ($step && $step > 0) ) )	{
	# Set minimum kmer values
	if ( $mink && $mink > 0 )	{
		if ( $maxk && $maxk < $mink )	{
			( $mink, $maxk ) = ( $maxk, $mink );		# Reverse the min and max if max is less than min
		}
	}	
	else		{
		$mink = 20;
	}
	
	# Set maximum kmer values
	if ( $maxk && $maxk > 0 )	{
		if ( $mink && $mink > $maxk )	{
			( $mink, $maxk ) = ( $maxk, $mink );		
		}
	}	
	else		{
		$maxk = 40;
	}

	# Set the step size
	if ( $step && $step > 0 )	{
		# Do nothing
	}
	else		{
		$step = 10;
	}	
	$k = 30;		# Just set this to a default value -- we wont be using it anyways, but we can avoid unitialized warnings :)
}
else		{
	$k = 30;
	$mink = $k;
	$maxk = $k;
	$step = 1;
}

#print "K:		$k\n";
#print "MINK:	$mink\n";
#print "MAXK:	$maxk\n";
#print "STEP:	$step\n"; 
#print "\n\n\n";


# Store the FASTA sequences in a hash
$/ = ">";
my @fastas;
if ( -e $input && $input ne "--" )	{
	open DATA, $input or die "FASTA FILE ERROR: $!\n";
		@fastas = <DATA>;
		my $trash = shift @fastas;	# Get rid of the first element, which will be a lone ">" symbol
	close DATA;
}
elsif ( $input eq "--" )	{
	@fastas = <STDIN>;
	my $trash = shift @fastas;
}
$/ = "\n";

my %Alignment = ();
foreach my $record ( @fastas )	{
	my ($header, @seq) = split "\n", $record;
	my $seq = join '', @seq;
	my @headers = split " ", $header;
	
	$headers[0] =~ s/ /_/g;				# Convert any spaces to underscores
	$headers[0] =~ s/-//g;				# Convert any hyphens in SEQUENCE NAMES to underscores (Paup really doesnt like these)
	$headers[0] =~ s/\s+//g;			# Remove any other odd whitespace
	$headers[0] =~ s/\|//g;				# Remove pipe chars (often seen in Genbank headers)
	$seq =~ s/>//g;						# Remove any lingering ">" symbols
	$seq = uc $seq;						# Convert everything to uppercase 
	
	# Store the sequences as a hash
	$Alignment{$headers[0]} = $seq;
}

my %Kmers = ();
foreach my $sequence ( values %Alignment )	{
	chomp $sequence;
	my $rcseq = revcom( $sequence );
	
	# For each kmer size:
	for ( my $i = $mink; $i <= $maxk; $i += $step )	{
#		print "I:		$i\n";
#		print "SeqLen:	", length($sequence), "\n";
		
		for ( my $j = 0; $j < length($sequence); $j += $i )	{
#			print "\t\t\tJ:		$j\n";
			
			# Get the substr of the sequence and its revcomp
			my $kmer = substr($sequence, $j, $i);
			my $rckmer = substr($rcseq, $j, $i);
			next if ( length($kmer) < $i || length($rckmer) < $i );
			
			# Increment the counters
			if ( exists $Kmers{$i}{$kmer} )	{
				$Kmers{$i}{$kmer} += 1;
			}
			else	{
				$Kmers{$i}{$kmer} = 1;
			}
			if ( exists $Kmers{$i}{$rckmer} )	{
				$Kmers{$i}{$rckmer} += 1;
			}
			else	{
				$Kmers{$i}{$rckmer} = 1;
			}
			
		}
	}
}

# Print out the output
if ( $output ne "--" )	{
	foreach my $kmersize ( sort { $a <=> $b } keys %Kmers )	{
		open OUT, ">", "$output.kmers.$kmersize.tab" or die "OUTPUT FILE ERROR: $!\n";
		my %Size = %{ $Kmers{$kmersize} };
		foreach my $subseq ( sort { $a cmp $b } keys %Size )	{
			print OUT "$subseq\t$Size{$subseq}\n";
		}
	}
}
else		{
	foreach my $kmersize ( sort { $a <=> $b } keys %Kmers )	{
		my %Size = %{ $Kmers{$kmersize} };
		foreach my $subseq ( sort { $a cmp $b } keys %Size )	{
			print STDOUT "$subseq\t$Size{$subseq}\n";
		}
	}
}


exit;

################# SUBROUTINES #########################

# Reverse complements a sequence
sub revcom	{
	my ( $seq ) = @_;
	my $revcom = $seq;
	
	$revcom =~ tr/ACGTacgt/TGCAtgca/;
	$revcom = reverse $revcom;
	
	return ( uc $revcom );
}

