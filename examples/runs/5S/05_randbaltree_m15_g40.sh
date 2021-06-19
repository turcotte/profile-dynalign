#! /bin/sh

SCRIPT=`pwd`/../../../scripts/randbaltree.prl
INPDIR=`pwd`/../../data/5S
OUTDIR=`pwd`/05

exec $SCRIPT $INPDIR $OUTDIR -m 15 -g 4
