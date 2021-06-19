#! /bin/sh

SCRIPT=`pwd`/../../../scripts/clustalwtree.prl
INPDIR=`pwd`/../../data/tRNA
OUTDIR=`pwd`/04

exec $SCRIPT $INPDIR $OUTDIR -m 25 -g 4
