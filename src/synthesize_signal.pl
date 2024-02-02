#!/usr/bin/perl -w

use strict;

die "usage: $0 OUTPUT PROMENHfile ENHSIGNALfile (opt)ENHINPUTfile"
  unless $#ARGV>=2;

my @lineArr;
my %signalHash;
my $regionID;
my $prom;
my @enhArr;
my %signalHashPerGene;
my %inputHash;

open(SIGNAL, "<$ARGV[2]");
while(<SIGNAL>) {
  chomp($_);
  @lineArr = split("\t", $_);
  $regionID=$lineArr[0];
  $signalHash{$regionID} = $lineArr[1];
}
close(SIGNAL);

if($ARGV[3]) {
  open(INPUT, "<$ARGV[3]");
  while(<INPUT>) {
    chomp($_);
    @lineArr = split("\t", $_);
    $regionID=$lineArr[0];
#     print $regionID . "\n";
    $signalHash{$regionID} = $signalHash{$regionID} - $lineArr[1];
  }
  close(INPUT);
}


open(PROMENH, "<$ARGV[1]");
while(<PROMENH>) {
  chomp($_);
  @lineArr= split("\t", $_);
  $prom = $lineArr[0];
  $signalHashPerGene{$prom} = 0;

  if(exists($lineArr[1])) {
    @enhArr= split(",", $lineArr[1]);
    
    foreach my $enh (@enhArr) {
      $signalHashPerGene{$prom} += $signalHash{$enh};
    }
  }
}

open(OUT, ">$ARGV[0]");
foreach my $prom (keys %signalHashPerGene) {
  print OUT $prom . "\t" . $signalHashPerGene{$prom} . "\n";
}
close(OUT);


exit(0);