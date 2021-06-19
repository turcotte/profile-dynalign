#! /bin/sh

SCRIPT=`pwd`/../../../scripts/randprogtree.prl
INPDIR=`pwd`/../../data/tRNA
OUTDIR=`pwd`/06

exec $SCRIPT $INPDIR $OUTDIR -m 25 -g 4
