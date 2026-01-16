#!/usr/bin/perl
use strict;
use warnings;
use File::Basename;
use List::Util qw(first);
use 5.010;
local $SIG{__WARN__} = sub {
	my $message =shift;
	die $message;
};
#
# Author: Adapted for OpenCRAVAT annotation format
# This script takes annotated files from OpenCRAVAT and creates a table of all variants
# in the format to be loaded to the database
# Output: Chr\tStart\tEnd\tRef\tAlt\tSample\tCaller\tQUAL\tFS\tTotalReads\tAltReads
#
# OpenCRAVAT format has columns:
# - QUAL at column 60
# - Sample columns starting with .GT suffix (e.g., trim_sample.GT)
# - Pattern: .GT, TotalCoverage, RefCoverage, VarCoverage, Variant Allele Freq
#

my (%CALLER, %QUAL, %FS, %READS);

foreach my $file(@ARGV){
	unless (open(FH, $file)){
		print STDERR "Can not open the file $file\n";
		exit;
	}

	# Extract caller from filename (e.g., HC_tumor_RNA from filename.HC_tumor_RNA.annotated.txt)
	my @a = split(/[.]/, basename($file));
	my $caller = $a[1];

	my %samples;
	my $index_qual;
	my $index_fs;

	while(<FH>){
		chomp;
		my @line = split("\t", $_);

		# Process header line
		if($. ==1 and $_ =~ /^Chr/){
			# Find QUAL column
			$index_qual = first { $line[$_] eq 'QUAL' } 0..$#line;
			$index_fs   = first { $line[$_] eq 'QUAL' } 0..$#line;

			if(!defined $index_qual){
				print STDERR "ERROR: QUAL column not found in file $file\n";
				exit;
			}

			# Find all .GT columns (sample columns)
			# OpenCRAVAT pattern: .GT (col i), TotalCoverage (col i+1), RefCoverage (col i+2), VarCoverage (col i+3)
			for (my $i=0; $i<=$#line; $i++){
				if($line[$i] =~ /\.GT$/){
					my $tmp = $line[$i];
					$tmp =~ s/\.GT$//;
					$tmp =~ s/^trim_//;  # Remove trim_ prefix if present
					$samples{$i} = $tmp;
					print STDERR "Found sample: $tmp at column $i (TotalCoverage at ".($i+1).", VarCoverage at ".($i+3).")\n";
				}
			}

			if(scalar keys %samples == 0){
				print STDERR "ERROR: No sample columns (.GT) found in file $file\n";
				exit;
			}

			next;
		}

		# Skip lines with underscores in chromosome name
		if($line[0] =~ /_/){
			next;
		}

		# Create variant key: Chr, Start, End, Ref, Alt
		my $key  = "$line[0]\t$line[1]\t$line[2]\t$line[3]\t$line[4]";

		# Process each sample
		foreach my $sample_col (keys %samples){
			my $sample_name = $samples{$sample_col};
			my $total_reads = $line[$sample_col+1];  # TotalCoverage is at GT column + 1
			my $alt_reads   = $line[$sample_col+3];  # VarCoverage is at GT column + 3

			# Skip if coverage data is not numeric
			next unless ($total_reads =~ /^\d+$/ and $alt_reads =~ /^\d+$/);

			my $variant_key = "$key\t$sample_name";

			if ( not exists $READS{$variant_key}){
				$READS{$variant_key} = "$total_reads\t$alt_reads";
				$QUAL{$variant_key}  = "$line[$index_qual]";
				$FS{$variant_key}    = "$line[$index_fs]";
				$CALLER{$variant_key} = "$caller";
			}
			else{
				# If variant already seen (from another caller), append with semicolon
				$QUAL{$variant_key}   = $QUAL{$variant_key}.";$line[$index_qual]";
				$FS{$variant_key}     = $FS{$variant_key}.";$line[$index_fs]";
				$CALLER{$variant_key} = $CALLER{$variant_key}.";$caller";
			}
		}
	}
	close FH;
}

# Print the data
# Format: Chr\tStart\tEnd\tRef\tAlt\tSample\tCaller\tQUAL\tFS\tTotalReads\tAltReads
foreach (sort keys %READS) {
	my ($total, $alt) = split("\t", $READS{$_});
	if ($total =~ /^\d+$/ and $alt =~ /^\d+$/){
		print "$_\t$CALLER{$_}\t$QUAL{$_}\t$FS{$_}\t$READS{$_}\n";
	}
}
