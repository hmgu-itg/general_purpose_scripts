#!/usr/bin/perl -w

use strict;

$\="\n";
$,="\t";

if ($ARGV[0] ne "mean" && $ARGV[0] ne "majority" && $ARGV[0] ne "cc"){
    print STDERR "ERROR: unrecognized mode: $ARGV[0]";
    exit(1);
}

# first field is ID, skipping when collapsing
while(<STDIN>){
    chomp;
    my @a=split(/\t/,$_,-1);
    if ($ARGV[0] eq "mean"){
	my $sum=0;
	my $count=0;
	for (my $i=1;$i<scalar(@a);$i++){if ($a[$i] ne "NA"){$count++;$sum+=$a[$i];}}
	if ($count==0){
	    print $a[0],"NA";
	}
	else{
	    print $a[0],$sum/$count;
	}
    }elsif($ARGV[0] eq "majority"){
	my %H=();
	for (my $i=1;$i<scalar(@a);$i++){if ($a[$i] ne "NA"){$H{$a[$i]}++;}}
	if (scalar(keys %H)==0){
	    print $a[0],"NA";
	}
	else{
	    my @b=(reverse sort {$H{$a} <=> $H{$b}} keys %H);
	    print $a[0],$b[0];
	}
    }elsif($ARGV[0] eq "cc"){
	my $val=$ARGV[1];
	my $flag=0;
	foreach my $x (@a){if ($x eq $val){$flag=1;last;}}
	print $a[0],$flag;
    }
}
