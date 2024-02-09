#!/usr/bin/perl -w

use strict;
use List::Util qw( min );

$,="\t";
$\="\n";

my %D;

sub get_children{
    my $r=shift;
    my $parent=shift;
    my @a;
    foreach my $k (keys %{$r}){
	if ($k=~/^$parent\/[^\/]+$/){
	    push @a,{"name"=>$k,"size"=>$r->{$k},"children"=>get_children($r,$k)};
	}
    }
    return \@a;
}

sub hash2str{
    my $r=shift;
    my @a;
    for my $c (@{$r->{children}}){
	push @a,hash2str($c);
    }
    if (scalar(@a)==0){
	return "\"".$r->{name}." (".$r->{size}.")\":{}";
    }
    else{
	return "\"".$r->{name}." (".$r->{size}.")\":{".join(",",@a)."}";
#	return "\"".$r->{name}."\":{\"size\":\"".$r->{size}."\",".join(",",@a)."}";
    }
}

while(<STDIN>){
    chomp;
    my @a=split(/\s+/,$_,2);
    if ($a[1]=~/(.*)\/$/){
	$a[1]=$1;
    }
    $D{$a[1]}=$a[0];
}

my %H;
foreach my $k (keys %D){
    $H{$k}=$k=~tr/\///;
    #print $k,$D{$k},$H{$k};
}

my $minv=min values(%H);
#print $minv;
my @mink;

foreach my $k (keys %H){
    if ($H{$k}==$minv){
	push @mink,$k;
    }
}

#print join(",",@mink);
my $root=$mink[0];
#print $root;

my $r=get_children(\%D,$root);
# foreach my $c (@$r){
#     print $c->{name},$c->{size};
#     print hash2str($c);
# }

print hash2str({name=>$root,size=>$D{$root},children=>$r});

exit(0);

