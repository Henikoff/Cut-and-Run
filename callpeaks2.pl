#!/usr/bin/perl -w

# Siva Kasinathan, Henikoff Lab/FHCRC (skasin@uw.edu)
# Improved peakcaller: requires peaks to have a minimum width (i.e.,
# a minimum number of continuous bases above threshold). This should
# allow the usage of a lower threshold and calling broader peaks

# Input BED file should be sorted: sort -k1,1 -k2,2n BED_file > BED_file_sorted

use strict;
use autodie;

die "Usage: callpeaks2.pl <sorted bedgraph> [threshold] [interpeak distance] [min width] [max width]\n"
  if (! defined $ARGV[4]);

my $bed_file = $ARGV[0];
my $threshold = $ARGV[1];
my $inter_peak = $ARGV[2];
my $minw = $ARGV[3];
my $maxw = $ARGV[4];

print STDERR "BED File:            ", $bed_file, "\n";
print STDERR "Threshold:           ", $threshold, "\n";
print STDERR "Inter-peak distance: ", $inter_peak, "\n";
print STDERR "Minimum peak width:  ", $minw, "\n";
print STDERR "Maximum peak width:  ", $maxw, "\n";

# Track number of peaks and peaks
# that are thrown out
my $npeaks = 0;
my $nfail = 0;

my @peak;
my @current;
reset_peak(\@peak); # Initialize peak

open BED, '<', $bed_file
  or die "Could not open input BED file\n";

while(<BED>){
  chomp;
  my @line = split/\s+/;

  # Assume bedgraph: chr start end value
  # (can change these depending on input BED format)
  $current[0] = $line[0];
  $current[1] = $line[1];
  $current[2] = $line[2];
  $current[3] = $line[3];

  if ($current[3] >= $threshold){ # process if above thresh
    my $add_success = update_peak(\@peak,\@current, \$inter_peak);

    if (! $add_success  ){

      print_peak(\@peak, \$minw, \$maxw, \$npeaks, \$nfail);
      reset_peak(\@peak);
      update_peak(\@peak,\@current);

    }
  }

}

close BED;

# Print last peak
print_peak(\@peak, \$minw, \$maxw, \$npeaks, \$nfail);

print STDERR "Number of called peaks:     ", $npeaks, "\n";
print STDERR "Number of peaks thrown out: ", $nfail, "\n";
print STDERR "DONE.\n";

exit;

sub update_peak{
  my $peak = $_[0];
  my $current = $_[1];
  my $inter_peak = $_[2];

  my $result = 0;

  # Return failure if peak and the
  # current entry are not on the same
  # chromosome or if the current
  # entry is too far from peak
  if ($$peak[0] ne "*"){
    return 0 if ($$peak[0] ne $$current[0]);
    return 0 if ($$current[1]- $$peak[2] >= $$inter_peak);
  } else {
    $$peak[0] = $$current[0];
    $$peak[1] = $$current[1];
  }

  # Add current entry to peak
  $$peak[2] = $$current[2];   # Extend peak end:
  $$peak[3] += $$current[3];  # Increment occupancy

  return 1;
}

sub print_peak{
  my $peak = $_[0];
  my $minw = $_[1];
  my $maxw = $_[2];
  my $npeaks = $_[3];
  my $nfail = $_[4];

  # Only print if peak width is within
  # bound specified by minw and maxw
  if ($$peak[2] - $$peak[1] >= $$minw && $$peak[2] - $$peak[1] <= $$maxw){
    $$npeaks++;
    print $$peak[0], "\t",
          $$peak[1], "\t",
          $$peak[2], "\t",
          $$peak[3], "\n";
  } else {
    $$nfail++;
  }

  return;
}

sub reset_peak{
  my $peak = $_[0];
  $$peak[0] = "*";
  $$peak[1] = 0;
  $$peak[2] = 0;
  $$peak[3] = 0;
  return;
}
