#! /bin/sh

SCRIPT=`pwd`/../../../scripts/bestprog.prl
INPDIR=`pwd`/../../data/tRNA
OUTDIR=`pwd`/03
ALIDIR=`pwd`/01

exec $SCRIPT $INPDIR $OUTDIR $ALIDIR -m 25 -g 4
