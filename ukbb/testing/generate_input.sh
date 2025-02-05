#!/usr/bin/env bash

##################################################################
#           CREATE TEST INPUT FILES WITH RANDOM ENTRIES
#
# $1: average number of rows in output files
# $2: number of columns in output files
# $3: number of output files
# $4: output files prefix
#
##################################################################

rows=$1
cols=$2
N=$3
prefix=$4

m=$((cols*N))
n=${#m}
declare -a rows
mr=0
for i in $(seq -w 1 $N);do
    del=$((RANDOM%2))
    # increase/decrease by up to 5%
    max_mod=$(echo $rows/20|bc)
    to_mod=$((RANDOM%max_mod))
    if [[ $del -eq 1 ]];then
	to_mod="-"$to_mod
    fi
    r=$((rows+to_mod))
    rows+=($r)
    if [[ $r -gt $mr ]];then
	mr=$r
    fi
done

   
for i in $(seq -w 1 $N);do
    echo -n "Generating file $i/$N" 1>&2
    fname=${prefix}_${i}".txt.gz"
    r=${rows[$i]}
    echo " ($r rows)"
    
    echo -ne "f.eid\t" | gzip -c > $fname
    for j in $(seq $(((i-1)*cols+1)) $((i*cols)));do echo $j;done|perl -slne 'BEGIN{@a=();}{push @a,sprintf("COL%0".$w."d",$_);}END{print join("\t",@a);}' -- -w=$n | gzip -c >> $fname
    for j in $(seq -w 1 $mr);do
	echo -ne "$j\t"
	tr -dc 0-9 </dev/urandom | head -c $(echo "6*$cols" | bc -l)|perl -lne 'print join("\t",/(.{1,6})/g);'|gawk -v seed=$RANDOM 'BEGIN{FS=OFS="\t";srand(seed);}{x=rand();n=1+int(NF*rand());if (x>=0 && x<0.05){$n="";} if (x>=0.05 && x<0.1){$n="\""$n"\"";} if(x>=0.1 && x<0.15){$n=rand();} if(x>=0.15 && x<0.2){$n=rand()", "rand();} if(x>=0.2 && x<0.25){$n="\""rand()", "rand()"\"";} if(x>=0.25 && x<0.3){$n=" ";} print $0;}'
    done | shuf -n $r | gzip - -c >> $fname
done
