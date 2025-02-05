#!/usr/bin/perl -w

use strict;

$\="\n";
while(<STDIN>){
    chomp;
    my @a=split(/\t/,$_,-1);
    for (my $i=0;$i<scalar(@a);$i++){$a[$i]="NA" if length($a[$i])==0;}
    print join("\t",@a);
}
