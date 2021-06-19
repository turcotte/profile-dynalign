#! /bin/sh

SCRIPT=`pwd`/../../../scripts/randbaltree.prl
INPDIR=`pwd`/../../data/tRNA
OUTDIR=`pwd`/05

exec $SCRIPT $INPDIR $OUTDIR -m 25 -g 4
