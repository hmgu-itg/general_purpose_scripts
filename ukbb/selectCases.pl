#!/usr/bin/perl -w

use strict;

$\="\n";

my $fh;
my %h=();
my $valcol=$ARGV[1];
open($fh,"<",$ARGV[0]);
while(<$fh>){chomp;next if /^#/;next if /^\s*$/;$h{$_}=1;}
close($fh);
my @expressions=();
foreach my $x (keys %h){
    $x=~s/or/||/ig;
    $x=~s/and/&&/ig;
    $x=~s/([a-zA-Z\d]+)/defined\(\$r->\{$1\}\)/g;
    push @expressions, "(".$x.")";
}
my $code=join(" || ",@expressions);
sub F{
    my $r=shift;
    my $c=shift;
    eval $c;
}

while(<STDIN>){
    chomp;
    my %D=();
    my @a=split(/\t/,$_,-1);
    foreach my $x (split(/,/,$a[$valcol],-1)){$D{$x}=1;}
    print $a[0] if F(\%D,$code);
}

