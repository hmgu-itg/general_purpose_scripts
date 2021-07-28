#!/usr/bin/perl

use strict;

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
    $f=~s/([a-zA-Z\d]+)/defined\(\$h\{\1\}\)/g;

    eval $f;
}

my $fh;
my @expressions;

# expressions
open($fh,"<",$ARGV[0]);
while(<$fh>){chomp;next if /^#/;next if /^\s*$/;push @expressions, $_;}
close($fh);

# data
open($fh,"<",$ARGV[1]);
while(<$fh>){
    chomp;
    my @a=split(/\t/,$_,-1);
    my $flag=0;
    foreach my $expr (@expressions){
	if (check($a[$ARGV[2]],$expr)){
	    $flag=1;
	    last;
	}
    }

    print $a[0] if ($flag==1);
}
close($fh);
