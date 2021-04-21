#! /bin/bash

# run a succession of scripts on the raw mdat data files

filename=$1
basefilename=${filename%.mdat}

mdatfile=$filename
rootfile=${basefilename}".root"
sortedfile=${basefilename}"_sorted.root"

root -q -b 'mdat_conv.C("'$mdatfile'")'

root -q -b 'Sorter.C("'$rootfile'")'

root -q -b 'Correlator.C("'$sortedfile'")'
