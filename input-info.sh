#!/bin/bash
p=$1
path=$( echo ${p%/*} )
file=$( echo ${p##*/} )

echo $path, $file

wc -l $file | awk '{print "number of genes: " $1; exit}'

awk -F' ' '{print "number of experiments: " NF-1; exit}' $file