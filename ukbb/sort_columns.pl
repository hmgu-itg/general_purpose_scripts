#!/usr/bin/perl -w

use strict;

$\="\n";
my $nr=0;
my %name2i;
my @S;
while(<STDIN>){
    chomp;
    my @a=split(/\t/,$_,-1);
    if ($nr == 0){
	push @S,$a[0];
	for (my $i=0;$i<scalar(@a);$i++){
	    $name2i{$a[$i]}=$i;
	}
	foreach my $c (sort keys %name2i){
	    if ($c ne $S[0]){
		push @S,$c;
	    }
	}
	print STDOUT join("\t",@S);
    }else{
	my @L=();
	foreach my $c (@S){
	    push @L,$a[$name2i{$c}];
	}
	print STDOUT join("\t",@L);
    }
    $nr++;
}
