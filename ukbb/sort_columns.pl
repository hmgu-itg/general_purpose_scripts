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
	for (my $i=1;$i<scalar(@a);$i++){
	    push @{$name2i{$a[$i]}},$i;
	}
	foreach my $c (sort keys %name2i){
	    push @S,$c;
	}
	print STDOUT join("\t",@S);
    }else{
	my @L=();
	foreach my $c (@S){
	    if ($c eq "f.eid"){
		push @L,$a[0];
	    }
	    else{
		push @L,$a[${$name2i{$c}}[0]];
		push @L,$a[${$name2i{$c}}[1]];
	    }
	}
	print STDOUT join("\t",@L);
    }
    $nr++;
}
