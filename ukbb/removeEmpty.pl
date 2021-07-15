#!/usr/bin/perl

$\="\n";
while(<STDIN>){
    chomp;
    @a=split(/\t/,$_,-1);
    for ($i=0;$i<scalar(@a);$i++){$a[$i]="NA" if length($a[$i])==0;}
    print join("\t",@a);
}
