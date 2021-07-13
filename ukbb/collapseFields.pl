#!/usr/bin/perl

$\="\n";
$,="\t";

if ($ARGV[0] ne "mean" && $ARGV[0] ne "majority" && $ARGV[0] ne "cc"){
    print STDERR "ERROR: unrecognized mode: $ARGV[0]";
    exit(1);
}



# first field is ID, skipping when collapsing
while(<STDIN>){
    chomp;
    @a=split(/\t/,$_,-1);
    if ($ARGV[0] eq "mean"){
	$sum=0;
	$count=0;
	for ($i=1;$i<scalar(@a);$i++){if ($a[$i] ne "NA"){$count++;$sum+=$a[$i];}}
	if ($count==0){
	    print $a[0],"NA";
	}
	else{
	    print $a[0],$sum/$count;
	}
    }elsif($ARGV[0] eq "majority"){
	%H=();
	for ($i=1;$i<scalar(@a);$i++){if ($a[$i] ne "NA"){$H{$a[$i]}++;}}
	if (scalar(keys %H)==0){
	    print $a[0],"NA";
	}
	else{
	    @b=(reverse sort {$H{$a} <=> $H{$b}} keys %H);
	    print $a[0],$b[0];
	}
    }elsif($ARGV[0] eq "cc"){
	$val=$ARGV[1];
	$flag=0;
	foreach $x (@a){if ($x eq $val){$flag=1;break;}}
	print $a[0],$flag;
    }
}
