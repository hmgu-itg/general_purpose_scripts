#!/usr/bin/perl -w

use strict;

$\="\n";
$,="\t";

sub checkF{
    my $expr_str=shift;
    my $r=shift;
    eval $expr_str;
}

my $max_visits=10;
my $max_arrays=10;

my $case1_str=$ARGV[0];
my $case2_str=$ARGV[1];
my $icd9_str=$ARGV[2];
my $icd10_str=$ARGV[3];
my $opcodes_str=$ARGV[4];

my (%case1,%case2,%sicd9,%sicd10,%ICD9,%ICD10,%OP,%fullOP,%sop);

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
	$fullOP{$z}=1;
    }
}

for my $x ("Q" .. "T"){
    for my $i (0 .. 99){
	my $z=$x.sprintf("%02d",$i);
	$ICD9{$z}=1 if (!defined($sicd9{$z}));
    }
}

my @code=();
foreach my $x (keys %sop){
    $x=~s/or/||/ig;
    $x=~s/and/&&/ig;
    $x=~s/([a-zA-Z\d]+)/defined\(\$r->\{\1\}\)/g;
    push @code, "(".$x.")";
}
my $exs=join(" || ",@code);

# --------------------------------------------------------------------------------

while(<STDIN>){
    chomp;
    my $ID=$_;
    my $inst=int(rand($max_visits));
    my $arrays=int(rand($max_arrays));
    my $hr=\%ICD10;
    my $icd9=0;
    if (rand()<0.1){
	$hr=\%ICD9;
	$icd9=1;
    }
    
    if (defined($case1{$ID})){

    }
    elsif(defined($case2{$ID})){


    }
    else{
	for (my $i=0;$i<=$inst;$i++){
	    for (my $j=0;$j<=$arrays;$j++){
		if ($icd9){
		    print $ID,$i,$j,(keys %$hr)[rand keys %$hr],"NA",(keys %fullOP)[rand keys %fullOP];
		}
		else{
		    print $ID,$i,$j,"NA",(keys %$hr)[rand keys %$hr],(keys %fullOP)[rand keys %fullOP];
		}
	    }
	}

    }
}

