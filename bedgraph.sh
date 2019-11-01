#!/bin/sh -e
mkdir -p igv
cd cnt
for i in *.bed.gz; do
	cat <(echo 'track type=bedGraph visibility=full color=200,100,0 altColor=0,100,200 priority=20 height=200 name='$i) <(gunzip -c $i) |\
		gzip -c >../igv/`echo $i | sed 's/\.bed\./.bedgraph./'`
done
