#! /bin/sh

SCRIPT=`pwd`/../../../scripts/clustalwtree.prl
INPDIR=`pwd`/../../data/5S
OUTDIR=`pwd`/04

exec $SCRIPT $INPDIR $OUTDIR -m 15 -g 4
