#!/usr/bin/env perl

use strict;

my %f2c;
my %c2miss;
my $i2c;
my $iid;
my $nr=0;
while(<STDIN>){
    chomp;
    my @a=split(/\t/);
    if ($nr==0){
	for (my $i=0;$i<scalar(@a);$i++){
	    my $c=$a[$i];
	    $i2c{$i}=$c;
	    if ($c eq "f.eid"){
		$iid=$i;
		next;
	    }
	    $c2miss{$c}=0;
	    my $f=undef;
	    if ($c=~/^f\.(\d+)\.\d+\.\d+$/){
		$f=$1;
	    }
	    else{
		print STDERR "ERROR: wrong column format: $c";
		exit(1);
	    }
	    push @{$f2c{$f}},$c;
	}
    }
    elsif($nr>1){
	for (my $i=0;$i<scalar(@a);$i++){
	    next if ($i==$iid);
	    my $c=$a[$i];
	    $c2miss{$i2c{$i}}++ if ($c eq "NA"); 
	}
    }
    $nr++;
}

$,="\t";
foreach my $f (sort keys %f2c){
    foreach my $c (@{$f2c{$f}}){
	print $f,$c,$c2miss{$c};
    }
}
