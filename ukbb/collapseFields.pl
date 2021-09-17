#!/usr/bin/perl -w

# ARGV[0]: "mean" OR "majority" OR "cc"
# if ARGV[0]=="cc", then ARGV[1] is the "case" value; ARGV[1] must not be "NA"
# if ARGV[1] occurs among the input fileds, then output is "1" ("case")
# otherwise if all input fields are "NA", the output is "NA"
# otherwise output is "0"

use strict;

$\="\n";
$,="\t";

if ($ARGV[0] ne "mean" && $ARGV[0] ne "majority" && $ARGV[0] ne "cc"){
    print STDERR "ERROR: unrecognized mode: $ARGV[0]";
    exit(1);
}

# first field is ID, skipping when collapsing
# fields are tab-separated
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
	my $val=$ARGV[1]; # value corresponding to "case"
	if ($val eq "NA"){
	    print STDERR "ERROR: collapseFields.pl: for case/control the case value must not be \"NA\"";
	    exit 1;
	}
	my $flag=0;
	my $flagNA=1; # if all fields are NA
	for (my $i=1;$i<scalar(@a);$i++){if ($a[$i] eq $val){$flag=1;$flagNA=0;last;} if ($a[$i] ne "NA"){$flagNA=0;}}
	if ($flag==1){
	    print $a[0],"1";
	}
	elsif($flagNA==1){
	    print $a[0],"NA";
	}
	else{
	    print $a[0],"0";
	}
    }
}
