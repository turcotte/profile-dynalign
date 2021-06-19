#! /bin/sh

SCRIPT=`pwd`/../../../scripts/allpairs.prl
INPDIR=`pwd`/../../data/5S
OUTDIR=`pwd`/01

exec $SCRIPT $INPDIR $OUTDIR -m 15 -g 4
