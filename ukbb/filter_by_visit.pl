#!/usr/bin/perl -w

use strict;

# -1 : d1<d2
sub compare_dates{
    my $x1=shift(@_);
    my $x2=shift(@_);

    my @m1=$x1=~m/(\d{2})\/(\d{2})\/(\d{4})/;
    my @m2=$x2=~m/(\d{2})\/(\d{2})\/(\d{4})/;

    return -1 if ($m1[2]<$m2[2]);
    return 1 if ($m1[2]>$m2[2]);
    return -1 if ($m1[1]<$m2[1]);
    return 1 if ($m1[1]>$m2[1]);
    return -1 if ($m1[0]<$m2[0]);
    return 1 if ($m1[0]>$m2[0]);
    return 0;
}

my $infile=$ARGV[0];
my $visits=$ARGV[1];
my $instance=$ARGV[2];

$\="\n";
$,="\t";
my $fh;
my %VC; # instance --> column
my %VCD; # ID --> [ visit_0, visit_1, ... ]
my $min_inst=undef;
my $max_inst=undef;

# visits
open($fh,"<",$visits);
while(<$fh>){
    chomp;
    my @a=split(/\t/);
    if ($.==1){
	for (my $i=1;$i<scalar(@a);$i++){
	    $a[$i]=~m/^f\.\d+\.(\d)\.\d+/;
	    if (!defined($min_inst)){
		$min_inst=$1;
	    }
	    else{
		$min_inst=$1 if ($min_inst>$1);
	    }
	    if (!defined($max_inst)){
		$max_inst=$1;
	    }
	    else{
		$max_inst=$1 if ($max_inst<$1);
	    }	    
	    $VC{$1}=$i;
	}
	next;
    }

    # $VCD{$a[0]}=();
    my @b=();
    # recode date
    for (my $i=1;$i<scalar(@a);$i++){
	if ($a[$i]=~m/(\d{4})-(\d{2})-(\d{2})/){
	    push @b,$3."/".$2."/".$1;
	}
	else{
	    push @b,$a[$i];
	}
    }
    # print "array b=".join(",",@b);
    for (my $i=$min_inst;$i<=$max_inst;$i++){
	# print "instance=".$i." column=".$VC{$i}." b=".$b[$VC{$i}-1];
	push @{$VCD{$a[0]}},$b[$VC{$i}-1];
    }
    # print $a[0],join(",",@{$VCD{$a[0]}});
    
}
close($fh);

open($fh,"<",$infile);
while(<$fh>){
    chomp;
    my @a=split(/\t/);
    my $id=$a[0];
    my $d1=$a[1];
    if (defined($VCD{$id})){
	my $r=$VCD{$id};
	for (my $i=$min_inst;$i<=$instance;$i++){
	    my $d2=$$r[$i-$min_inst];
	    next if ($d2 eq "NA");
	    if (compare_dates($d1,$d2)!=1){
		# print $id,$d1."<=".$d2;
		print $id;
		last;
	    }
	}
    }
}
close($fh);
