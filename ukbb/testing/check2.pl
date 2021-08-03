#!/usr/bin/perl -w

use strict;
use Algorithm::Combinatorics qw(combinations);

$\="\n";

sub check{
    my $str=shift;
    my $f=shift;
    
    my %h=();
    foreach my $x (split(/,/,$str,-1)){
	$h{$x}=1;
    }

    $f=~s/or/||/ig;
    $f=~s/and/&&/ig;
    $f=~s/([a-zA-Z\d]+)/defined\(\$h\{$1\}\)/g;

    eval $f;
}

sub minimal{
    my $ar=shift;
    my $to_add=shift;

    my %new_element=map { $_ => 1 } @$to_add;
    
    my $flag=1;
    foreach my $x (@$ar){
	my %h=map { $_ => 1 } @$x;
	my $f=1;
	foreach my $z (keys %h){
	    if (!defined($new_element{$z})){
		$f=0;
		last;
	    }
	}
	if ($f==1){
	    $flag=0;
	    last;
	}
    }

    return $flag;
}

my $expression=$ARGV[0];
my @terms=($expression=~/([A-Z]\d+)/g);
exit(0) if (scalar(@terms)==0);
my %T;
foreach my $x (@terms){
    $T{$x}=1;
}
my @output;
my @A=keys %T;
for (my $i=1;$i<=scalar(@A);$i++){
    foreach my $c (combinations(\@A,$i)){
	if (check(join(",",@$c),$expression) && minimal(\@output,$c)){
	    push @output,$c;
	}
    }
}

foreach my $c (@output){
    print join(",",@$c);
}
