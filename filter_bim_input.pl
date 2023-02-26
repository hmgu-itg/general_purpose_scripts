#!/usr/bin/perl -w

use strict;

$\="\n";

my $fh;
my @N;
my %D;

open($fh,"<", $ARGV[0]);
while(<$fh>){
    chomp;
    my @a=split(/\t/);
    push @N,$a[1];
}
close($fh);

open($fh,"<", $ARGV[1]);
while(<$fh>){
    chomp;
    my @a=split(/ /);
    $D{$a[0]}=join(" ",$a[0],$a[3],$a[4],$a[5],$a[6],$a[7],$a[8],$a[9]);
}
close($fh);

print "SNP A1 A2 freq beta se p N";
foreach my $x (@N){
    print $D{$x};
}
