#!/usr/bin/perl

$\="\n";

BEGIN{
    %h=();
    open(fh,"<",$ARGV[0]);
    while(<fh>){chomp;$h{$_}=1;}
    close(fh);
    @code=();
    foreach $x (keys %h){
	$x=~s/or/||/ig;
	$x=~s/and/&&/ig;
	$x=~s/([a-zA-Z]\d{3})/defined\(\$r->\{\1\}\)/g;
	push @code, "(".$x.")";
    }
    $c=join(" || ",@code);
    #print "\n".$c."\n\n";
    sub F{
	$r=shift;
	eval $c;
    }
}
{
    while(<STDIN>){
	chomp;
	%D=();
	@a=split(/\t/,$_,-1);
	foreach $x (split(/,/,$a[2],-1)){$D{$x}=1;}
	print $a[0] if F(\%D);
    }
}
