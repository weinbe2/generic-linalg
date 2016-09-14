#! /bin/bash

# This program outputs the time for a solve on its own line.
# This little perl command places the time back on the same
# line as the other inversion information.

infile=$1
outfile=$2

perl -00pe 's/\n(?=Time)/ /g' $infile > $outfile
