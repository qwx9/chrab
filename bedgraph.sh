#!/bin/sh -e
f=/tmp/bedgraph.tmp
mkdir -p igv
cd cnt
for i in *.bed.gz; do
	echo 'track type=bedGraph visibility=full color=200,100,0 altColor=0,100,200 priority=20 height=200 name='$i >$f
	gunzip -c $i >>$f
	gzip -c <$f >../igv/`echo $i | sed 's/\.bed\./.bedgraph./'`
done
rm -f $f
