#!/usr/bin/env perl

# two header lines in each chunk

use strict;

my $nr=0;
my %name2miss;
my %i2name;
while(<STDIN>){
    chomp;
    my @a=split(/\t/,$_,-1);
    if ($nr eq 0){
	for (my $i=0;$i<scalar(@a);$i++){
	    $i2name{$i}=$a[$i]; 
	}
    }elsif($nr gt 1){
	for (my $i=0;$i<scalar(@a);$i++){
	    $name2miss{$i2name{$i}}++ if ($a[$i] eq "NA"); 
	}
    }
    $nr++;
}

$,="\t";
$\="\n";
foreach my $n (keys %name2miss){
    print $n,$name2miss{$n};
}
