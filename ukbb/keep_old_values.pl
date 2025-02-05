#!/usr/bin/perl -w

use strict;

$\="\n";
my $colstr=$ARGV[0]; # 1-based, comma separated column numbers, alternating between file1 and file2
my %H=split(/,/,$colstr,-1);
my $nr=0;
while(<STDIN>){
    chomp;
    if ($nr == 0){
	print STDOUT $_;
    }else{
	my @a=split(/\t/,$_,-1);
	if ($a[1] ne "NA"){
	    print STDOUT $_;
	    next;
	}
	foreach my $k (keys %H){
	    $a[$H{$k}-1]=$a[$k-1];
	}
	print STDOUT join("\t",@a);
    }
    $nr++;
}
