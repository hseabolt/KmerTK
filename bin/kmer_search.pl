#!/usr/bin/perl

# kmer_search_sites.pl
# A plugin script...
# This one searches for kmers that have binding sites within a given range of sites (ie, we want to know what kmers are between sites 1000-2000)







# Need a set of binding sites, nothing else

###################################################################
# Read the binding sites and populate a hash
# Be mindful that depending on how the original set of kmers was constructed/filtered, not all sites may be covered.
my $fh = *STDIN;
my $succin = open(SITES, "<", "$input_sites") if ( $input_sites ne "--" && -e $input_sites );
$fh = *SITES if ( $succin );

my %BindingSites = ();
while ( <$fh> )	{
	chomp $_;
	my @line = split("\t", $_);
	my @sites = split(",", $line[1]);
	$BindingSites{ $line[0] } = \@sites;	
}
close SITES if ( $succin );
###################################################################


###################################################################
# The parser for the binding sites requested
my %Cfg = ();
if ( $cfg )	{
	my @raw;
	if ( -e $cfg )	{
		open CGET, $cfg;
			@raw = <CGET>;
		close CGET;
		chomp $_ foreach ( @raw );
	}
	else	{
		@raw = split(",", $cfg);
	}
	die " === KMER SEARCH Error :: Cannot search for anything if you don't tell me where to look! Bad search instructions given! ===\n" if ( scalar @raw == 0 );
	
	# Foreach interval, get the start position and the stop position and store them as a key-value pair in %Cfg.
	# If we only have a single position, start and stop are the same.
	foreach my $config ( @raw )	{
		# If we've got a range of columns
		if ( $config =~ /:|-/ )		{
			my @cols = split(/:|-/, $config);
			my ( $start, $stop ) = ( $cols[0], $cols[1] );
			$stop = $start if ( not $stop );
			( $start, $stop ) = ( $stop, $start ) if ( $stop < $start );
			$Cfg{$start} = $stop;
		}
		else	{
			$Cfg{$config} = $config;
		}
	}
}
###################################################################


###################################################################
# For each site in the range we want to search, grep that value from the hash values
# Return the list of kmers from the range
my @kmers;
foreach my $interval_start ( keys %Cfg )	{
	
}
###################################################################

###################################################################
# Output the same binding site format that we usually want.
# Note that each kmer is expected to have binding sites outside the only ones that we want, but we have to acknowledge those.
my $succout = open(OUT, ">", "$output.kmers.binding_sites.searched.tab") if ( $output ne "--" );
foreach my $kmer ( @kmers )	{
	if ( $succout ) 	{	print OUT    "$kmer\t", join(",", @{$BindingSites{$kmer}}), "\n";		}
	else 				{	print STDOUT "$kmer\t", join(",", @{$BindingSites{$kmer}}), "\n";		}
}	
close OUT if ( $succout );
###################################################################

exit;