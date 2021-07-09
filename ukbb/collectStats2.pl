#!/usr/bin/env perl

use strict;

my %i2miss;
while(<STDIN>){
    chomp;
    my @a=split(/\t/);
    for (my $i=0;$i<scalar(@a);$i++){
	$i2miss{$i}++ if ($a[$i] eq "NA"); 
    }
}

$,="\t";
$\="\n";
foreach my $i (sort keys %i2miss){
    print $i,$i2miss{$i};
}
