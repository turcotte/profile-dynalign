#! /bin/sh

SCRIPT=`pwd`/../../../scripts/allpairs.prl
INPDIR=`pwd`/../../data/tRNA
OUTDIR=`pwd`/01

exec $SCRIPT $INPDIR $OUTDIR -m 25 -g 4
