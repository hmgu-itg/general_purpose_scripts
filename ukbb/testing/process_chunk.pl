#!/usr/bin/perl -w

use strict;

$\="\n";
$,="\t";

my $max_visits=10;
my $max_arrays=10;

my $case1_str=$ARGV[0];
my $case2_str=$ARGV[1];
my $icd9_str=$ARGV[2];
my $icd10_str=$ARGV[3];
my $opcodes_str=$ARGV[4];

my (%case1,%case2,%sicd9,%sicd10,%ICD9,%ICD10,%OP,%sop);

$case1{$_}=1 for (split(/,/,$case1_str,-1));
$case2{$_}=1 for (split(/,/,$case2_str,-1));
$sicd9{$_}=1 for (split(/,/,$icd9_str,-1));
$sicd10{$_}=1 for (split(/,/,$icd10_str,-1));
$sop{$_}=1 for (split(/,/,$opcodes_str,-1));

for my $x ("A" .. "B"){
    for my $i (0 .. 999){
	my $z=$x.sprintf("%03d",$i);
	$OP{$z}=1 if (!defined($sop{$z}));
    }
}

for my $x ("K" .. "N"){
    for my $i (0 .. 99){
	my $z=$x.sprintf("%02d",$i);
	$ICD10{$z}=1 if (!defined($sicd10{$z}));
    }
}

for my $x ("Q" .. "T"){
    for my $i (0 .. 99){
	my $z=$x.sprintf("%02d",$i);
	$ICD9{$z}=1 if (!defined($sicd9{$z}));
    }
}

# --------------------------------------------------------------------------------

while(<STDIN>){
    chomp;
    my $ID=$_;
    if (defined($case1{$ID})){

    }
    elsif(defined($case2{$ID})){


    }
    else{


    }
    print $_,$case1_str;
}

print $_ for (sort keys %sop);

