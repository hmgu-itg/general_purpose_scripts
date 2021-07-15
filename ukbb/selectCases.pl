#!/usr/bin/perl

$\="\n";

BEGIN{
    %h=();
    $valcol=$ARGV[1];
    # print STDERR "$ARGV[0]\n";
    # print STDERR "$ARGV[1]\n";
    open(fh,"<",$ARGV[0]);
    while(<fh>){chomp;$h{$_}=1;}
    close(fh);
    # print STDERR scalar(keys %h)."\n";
    @code=();
    foreach $x (keys %h){
	$x=~s/or/||/ig;
	$x=~s/and/&&/ig;
	$x=~s/([a-zA-Z\d]+)/defined\(\$r->\{\1\}\)/g;
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
	foreach $x (split(/,/,$a[$valcol],-1)){$D{$x}=1;}
	print $a[0] if F(\%D);
    }
}
