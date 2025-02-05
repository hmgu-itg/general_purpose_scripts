#!/usr/bin/perl -w

use strict;
use List::Util qw(shuffle);

$\="\n";
$,="\t";

sub checkF{
    my $expr_str=shift;
    my $r=shift;
    eval $expr_str;
}

sub getSingleKey{
    my $r=shift;
    my $x=(keys %$r)[rand keys %$r];
    while($x=~/\s/){
	$x=(keys %$r)[rand keys %$r];
    }

    return $x;
}

my $max_visits=10;
my $max_arrays=10;

my $case1_str=$ARGV[0];
my $case2_str=$ARGV[1];
my $icd9_str=$ARGV[2];
my $icd10_str=$ARGV[3];
my $opcodes_str=$ARGV[4];

my (%case1,%case2,%sicd9,%sicd10,%ICD9,%ICD10,%fullICD9,%fullICD10,%OP,%fullOP,%sop);

$case1{$_}=1 for (split(/,/,$case1_str,-1));
$case2{$_}=1 for (split(/,/,$case2_str,-1));
$sicd9{$_}=1 for (split(/,/,$icd9_str,-1));
$sicd10{$_}=1 for (split(/,/,$icd10_str,-1));
for my $x (split(/,/,$opcodes_str,-1)){
    $x=~tr/_:?/ ()/;
    $sop{$x}=1;
}

for my $x ("A" .. "B"){
    for my $i (0 .. 999){
	my $z=$x.sprintf("%03d",$i);
	$OP{$z}=1 if (!defined($sop{$z}));
	$fullOP{$z}=1;
    }
}

for my $x ("K" .. "N"){
    for my $i (0 .. 99){
	my $z=$x.sprintf("%02d",$i);
	$ICD10{$z}=1 if (!defined($sicd10{$z}));
	$fullICD10{$z}=1;
    }
}

for my $x ("Q" .. "T"){
    for my $i (0 .. 99){
	my $z=$x.sprintf("%02d",$i);
	$ICD9{$z}=1 if (!defined($sicd9{$z}));
	$fullICD9{$z}=1;
    }
}

my @code=();
foreach my $x (keys %sop){
    $x=~s/or/||/ig;
    $x=~s/and/&&/ig;
    $x=~s/([a-zA-Z\d]+)/defined\(\$r->\{$1\}\)/g;
    push @code, "(".$x.")";
}
my $exs=join(" || ",@code);

sub selectOpTuple{
    my $n=shift;
    my $href=shift;
    my $expr_str=shift;

    my %D=();
    $D{$_}=1 for ((shuffle keys %$href)[0..$n-1]);
    
    while(checkF($expr_str,\%D)){
	%D=();
	$D{$_}=1 for ((shuffle keys %$href)[0..$n-1]);
    }

    return join(",",keys %D);
}

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
    my $selected_icd9=(keys %sicd9)[rand keys %sicd9];
    my $selected_icd10=(keys %sicd10)[rand keys %sicd10];
    
    if (defined($case2{$ID})){
	# case2: subset of case1
	for (my $i=0;$i<=$inst;$i++){
	    for (my $j=0;$j<=$arrays;$j++){
		if ($i==0 && $j==0){
		    my $z=getSingleKey(\%sop);
		    if ($icd9){
			print $ID,$i,$j,$selected_icd9,"NA",$z,"2";
		    }
		    else{
			print $ID,$i,$j,"NA",$selected_icd10,$z,"2";
		    }
		    next;
		}

		if ($icd9){
		    print $ID,$i,$j,(keys %fullICD9)[rand keys %fullICD9],"NA",(keys %fullOP)[rand keys %fullOP],"2";
		}
		else{
		    print $ID,$i,$j,"NA",(keys %fullICD10)[rand keys %fullICD10],(keys %fullOP)[rand keys %fullOP],"2";
		}
	    }
	}
    }
    elsif(defined($case1{$ID})){
	for (my $i=0;$i<=$inst;$i++){
	    my @icd9values=();
	    my @icd10values=();
	    if ($i==0){
		push @icd9values,$selected_icd9;
		for (my $j=1;$j<=$arrays;$j++){
		    push @icd9values,(keys %fullICD9)[rand keys %fullICD9];
		}
		push @icd10values,$selected_icd10;
		for (my $j=1;$j<=$arrays;$j++){
		    push @icd10values,(keys %fullICD10)[rand keys %fullICD10];
		}
	    }
	    else{
		for (my $j=0;$j<=$arrays;$j++){
		    push @icd9values,(keys %fullICD9)[rand keys %fullICD9];
		}
		for (my $j=0;$j<=$arrays;$j++){
		    push @icd10values,(keys %fullICD10)[rand keys %fullICD10];
		}
	    }
	    my $str=selectOpTuple($arrays+1,\%OP,$exs);
	    my @opcvalues=split(/,/,$str,-1);
	    for (my $j=0;$j<=$arrays;$j++){
		if ($icd9){
		    print $ID,$i,$j,$icd9values[$j],"NA",$opcvalues[$j],"1";
		}
		else{
		    print $ID,$i,$j,"NA",$icd10values[$j],$opcvalues[$j],"1";
		}
	    }	    
	}
    }
    else{
	for (my $i=0;$i<=$inst;$i++){
	    for (my $j=0;$j<=$arrays;$j++){
		if ($icd9){
		    print $ID,$i,$j,(keys %$hr)[rand keys %$hr],"NA",(keys %fullOP)[rand keys %fullOP],"0";
		}
		else{
		    print $ID,$i,$j,"NA",(keys %$hr)[rand keys %$hr],(keys %fullOP)[rand keys %fullOP],"0";
		}
	    }
	}
    }
}

