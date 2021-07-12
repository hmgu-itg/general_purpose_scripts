#!/usr/bin/env perl

use strict;

my %i2f;
my %i2c;
my %skip;
my $nr=0;
$\="\n";
$,="\t";
while(<STDIN>){
    chomp;
    my @a=split(/\t/,$_,-1);
    if ($nr==0){
	for (my $i=0;$i<scalar(@a);$i++){
	    my $c=$a[$i];
	    $i2c{$i}=$c;
	    if ($c eq "f.eid" || $c eq "RELEASE" || $c eq "CREATED"){
		$skip{$i}=1;
		next;
	    }
	    my $f=undef;
	    if ($c=~/^f\.(\d+)\.\d+\.\d+$/){
		$f=$1;
	    }
	    else{
		print STDERR "ERROR: wrong column format: $c";
		exit(1);
	    }
	    $i2f{$i}=$f;
	}
    }
    else{
	next if (defined($skip{$a[0]}));
	print $i2f{$a[0]},$i2c{$a[0]},$a[1];
    }
    $nr++;
}
exit(0);
