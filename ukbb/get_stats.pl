#!/usr/bin/perl -w

use strict;

$\="\n";
my $fh;
my $nr=0;
my %name2i;
my %tmp;
my %S;
my @L=("NA,NA","x,NA","NA,x","x,x","x,y");
open($fh,"<",$ARGV[0]);
while(<$fh>){
    chomp;
    my @a=split(/\t/,$_,-1);
    if ($nr == 0){
	for (my $i=0;$i<scalar(@a);$i++){
	    $tmp{$a[$i]}++;
	}
	for (my $i=0;$i<scalar(@a);$i++){
	    if ($tmp{$a[$i]}==2){
		push @{$name2i{$a[$i]}},$i;
		foreach my $l (@L){
		    $S{$a[$i]}{$l}=0;
		}
	    }
	}
    }else{
	foreach my $k (keys %name2i){
	    my $i=${$name2i{$k}}[0];
	    my $j=${$name2i{$k}}[1];
	    my $x1=$a[$i];
	    my $x2=$a[$j];
	    my $label;
	    if ($x1 ne "NA" && $x2 eq "NA"){
		$label="x,NA";
	    }elsif($x1 eq "NA" && $x2 ne "NA"){
		$label="NA,x";
	    }elsif($x1 ne "NA" && $x2 ne "NA" && $x1 eq $x2){
		$label="x,x";
	    }elsif($x1 ne "NA" && $x2 ne "NA" && $x1 ne $x2){
		$label="x,y";
	    }else{
		$label="NA,NA";
	    }
	    $S{$k}{$label}++;
	}
    }
    $nr++;
}
close($fh);

my @tmp=("STAT","C");
foreach my $l (@L){
    push @tmp,$l;
}
print join("\t",@tmp);
foreach my $k (keys %S){
    @tmp=("STAT",$k);
    foreach my $l (@L){
	push @tmp,$S{$k}{$l};
    }
    print join("\t",@tmp);
}

