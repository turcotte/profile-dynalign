#! /bin/sh

SCRIPT=`pwd`/../../../scripts/besttree.prl
INPDIR=`pwd`/../../data/tRNA
OUTDIR=`pwd`/02
ALIDIR=`pwd`/01

exec $SCRIPT $INPDIR $OUTDIR $ALIDIR -m 25 -g 4
