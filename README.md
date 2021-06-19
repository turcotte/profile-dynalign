# Profile-Dynalign 1.0

## INTRODUCTION

This is the software  system Profile-Dynalign. We extended the program
Dynalign (Mathews & Turner  (2002) JMB 317(2):191-203 --- an algorithm
for finding  the secondary structure  common to two RNA  sequences) to
align profiles  of multiple sequences.  We also  explored the relative
successes  of different  approaches to  designing the  tree  that will
guide progressive alignments of sequence profiles to create a multiple
sequence alignment and prediction of conserved structure.

## HISTORY

Profile-Dynalign  is  an  extension  of  Dynalign,  written  by  David
Mathews. An  initial prototype for  the alignment of two  profiles was
developed by Luke  Cen (Winter 2006). It was  then further extended by
Amelia  Bellamy-Royds (Summer  2006)  who wrote  the  source code  and
analysed  the  data  for  the  various approches  to  designing  guide
trees. The  work was  done under the  supervision of  Marcel Turcotte.

## PUBLICATIONS

Amelia   B.  Bellamy-Royds   and  Marcel   Turcotte   (Submitted)  Can
Clustal-Style Progressive Pairwise  Alignment of Multiple Sequences Be
Used in RNA Secondary Structure Prediction?

## INSTALLATION

PREREQUISITE:  emma from  EMBOSS is  required for  building  tree from
sequence information.

In order  to build the main  program, go to the  source code directory
(src), edit the make file to suite your needs and type make.  Copy the
executable  (pdynalign) to  a directory  specified by  the environment
PATH.

## ENVIRONMENT  VARIABLES:

DYNALIGN_DIR  should be pointing  at the  directory that  contains the
energy tables (datafiles).

PERLLIB should be pointing at the scripts directory.

## SCRIPTS

[ Main scripts ]

allpairs.prl: 

     Creates a set of all pairwise Profile-Dynalign alignments.

bestprog.prl: 

     Creates a Profile-dynalign alignment  in a progressive fashion by
     joining sequences  with the best energy scores  from all pairwise
     alignments (requires running allpairs.prl first).

besttree.prl:

     Creates  a Profile-dynalign  alignment using  a nearest-neighbour
     tree  based on  the energy  scores from  all  pairwise alignments
     (requires running allpairs.prl first).

clustalwtree.prl: 

    Creates  a  Profle-dynalign alignment  using  a phylogenetic  tree
    structure determined from a Clustal W primary sequence alignment.

randbaltree.prl:

    Creates a  Profile-Dynalign alignment using  a randomized balanced
    tree structure.

randprogtree.prl:

    Creates a  Profile-Dynalign alignment by  progressively adding one
    sequence in random order.

[ Utilities ]

seq2prf.prl:

    Converts a  Dynalign-style .seq  file into Stockholm  profile .prf
    format.

ct2mfoldaux.prl:

    Extracts a set of constraints for mfold from a connect (ct) file.

bracket2ct:

    Converts a file in bracket notation to ct format.

ct_compare:

    Compares  two ct  files (a  reference and  a predction)  so  as to
    produce statistics: ppv, coverage and MCC.

## USAGE

The directory  examples contains sample  scripts for each of  the five
approaches for building the guide  trees, as well as computing all the
pairwise alignments.

## CONTACT

Marcel Turcotte
marcel.turcotte@uottawa.ca
