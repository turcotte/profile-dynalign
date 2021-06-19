#! /bin/sh

SCRIPT=`pwd`/../../../scripts/randprogtree.prl
INPDIR=`pwd`/../../data/5S
OUTDIR=`pwd`/06

exec $SCRIPT $INPDIR $OUTDIR -m 15 -g 4
