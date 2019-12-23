#!/bin/sh
# export R sessionInfo with all used packages for README
grep -h 'library(' *.R |\
	sort |\
	uniq >/tmp/session.R
echo 'sessionInfo()' >>/tmp/session.R
Rscript /tmp/session.R 2>/dev/null
rm /tmp/session.R
