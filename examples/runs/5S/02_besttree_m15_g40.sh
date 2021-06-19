#! /bin/sh

SCRIPT=`pwd`/../../../scripts/besttree.prl
INPDIR=`pwd`/../../data/5S
OUTDIR=`pwd`/02
ALIDIR=`pwd`/01

exec $SCRIPT $INPDIR $OUTDIR $ALIDIR -m 15 -g 4
