#! /usr/bin/perl
#                              -*- Mode: Perl -*- 
# TreeBuilder.pm --- 
# Author          : Amelia Bellamy-Rods
# Created On      : Summer 2006
# Last Modified By: Marcel Turcotte
# Last Modified On: Tue Dec 12 12:33:17 2006
# 
# Amelia B. Bellamy-Royds and Marcel Turcotte (2006) Can Clustal-Style 
#   Progressive Pairwise Alignment of Multiple Sequences Be Used in RNA 
#   Secondary Structure Prediction?
#
# Copyright (C) 2006 University of Ottawa
#
# This copyrighted source code is freely distributed under the terms
# of the GNU General Public License. 
#
# See LICENSE for details.

package TreeBuilder;

use warnings;
use File::Spec::Functions;
use File::Basename;
use File::Copy;
# use POSIX;
use Cwd;
use Sys::Hostname;

##### Change the following constants as necessary #####

# File pattern to use to identify input files
$ProfilePattern = "*.prf";
$ProfileSuffix = ".prf";

# Path of hdynalign program for various hosts
%HDPath = ( "default" => "pdynalign" );

# Directory with hdynalign data files:

$HDDir  = $ENV{ "DYNALIGN_DIR" };

# "Approximately Infinity" (as far as any energy values are concerned):

$BigNumber = 30000;

##### These default HDynalign parameters can be changed 

# by setting the specified command-line option when
# running one of the tree-building scripts #####

$defaultMaxSep = 15;
$maxSepOption = "-m";

$defaultGapPenalty = 4.0; 
$gapPenaltyOption = "-g";

$defaultBPInsertionState = 1;
$BPInsertionOption = "-i";

################################################################

# SUBROUTINE &checkInput:
# Confirms that the first parameter is a valid existing directory
# and that the parameter value is a newly-made directory,
# making it if necessary.
# Returns the input array of parameters, if there was no error

sub checkInput (@) {

    # Check at least two parameters given
    unless (defined($_[0]) && defined($_[1]) ) { 
	die "Need input and output directories given as parameters\n--";
    }

    # Check input directory:
    unless (-d $_[0] ) { 
	die "$_[0] not a valid directory\n--";
    }

    # Check output directory:
    if ( -d $_[1] ) {
	# if given directory exists, create a unique sub-directory
	# to use for the results 
	# that way, the file tree will not get mixed up with anything else
	my $temp;
	while ( -e ($temp = catfile($_[1], "hd-".(int(rand(100000)) ) ) )) {}
	$_[1] = $temp; 
    }

    # ? restrict permissions on directory??
    mkdir $_[1] or die "Cannot make directory $_[1] :\n $!\n--";
    
    return @_;
}

# SUBROUTINE &getFiles:
# get all the matching input file names
# Requires a string naming a file directory as input
# Returns an array of filenames for all files in that
# directory which match $ProfilePattern

sub getFiles ($) {
	return glob (catfile($_[0], $ProfilePattern));
}

# SUBROUTINE &randomize:
# Randomize the order of the input parameters
# Returns an array containing the values
# in random order.

sub randomize (@) {

	# rename @_ so that changes don't affect original data
	my @Input = @_; 

	# create and set length of randomized array of file names
	my @Random;
	$#Random = $#Input; 

	# assign the values of @Input to @Random
	# in a random order:
	
	# while (@Input){
	# 	my $r = int (rand(@Input) );
	# 	$Random[@Input] = $Input[$r];
	# 	$Input[$r] = $Input[0];
	# 	shift @Input; 
	# }
	
	for (my $i = @Input - 1; $i >= 0; $i--) {
		my $r = int (rand ($i + 1) );
		$Random[$i] = $Input[$r];
		$Input[$r] = $Input[$i];

		# i.e.  Until there are no unassigned input files
		# left: Randomly pick an index of an unassigned
		# filename (initially all values in the Input array)
		# Assign the file name at this random index to a value
		# in the randomized array; Replace the used-up value
		# with an unassigned value so that all the unassigned
		# values will be at indices <= $i of param array at
		# the next loop (Note: still true if $r = $i, since
		# value at index $i won't be used in next loop)

	}
	return @Random;
}

# SUBROUTINE &copyFile
# Copies a file into a particular directory
# Paramaters (in order) are the directory and the filename w/ path
# The basename of the copied file is the same as the original.

sub copyFile ($$) {
	my ($dir, $file) = @_;
	copy($file, catfile($dir, basename($file) ))
		or die "Could not copy $file into $dir:\n $!\n--";
		
}

##The following routines are different ways to build
# a directory tree for a list of input files.
# In all cases, the result is a binary tree,
# with the files stored in the "leaf" directories.
# All require that the first parameter is the name
# of the top-level directory within which to build the
# tree, and the remaining parameters are the files
# to be copied into the leaves. ##

# SUBROUTINE &buildProgTree:
# Creates a tree such that each level consists of a single file
# added to a progressively larger tree, based on input order.

sub buildProgTree ($@){
	my $dir = shift @_;
	
	while (@_) {
	
		# create a directory for the first file, and save it
		mkdir ($NewDir = catfile($dir, "L1" )) 
			or die "Cannot make required output directory $NewDir:\n $!\n--";
		copyFile($NewDir, (shift @_) );
		
		# create a directory for remaining files
		mkdir ($NewDir = catfile($dir, "R".scalar(@_) ) )
			or die "Cannot make required output directory $NewDir:\n $!\n--";
			
		if (@_ == 1) {
			# then save single file in new directory
			copyFile($NewDir, (shift @_) );
		}
		$dir = $NewDir;
	}		

}

# SUBROUTINE &buildBalTree:
# Creates a balanced tree, based on the input order.

sub buildBalTree ($@); # pre-declare for recursion

sub buildBalTree ($@) {
	my $dir = shift @_;
	my $half = int (@_ / 2);
	if ($half) {
		# then at least two elements in @_, 
		# create directories and recurse
		# using complementary slices of the input array
		mkdir ($NewDir = catfile($dir, "L".($half) )) 
			or die "Cannot make required output directory $NewDir:\n $!";
		buildBalTree($NewDir, @_[0..($half-1)]);
		
		mkdir ($NewDir = catfile($dir, "R".(@_ - $half))) 
			or die "Cannot make required output directory $NewDir:\n $!";
		buildBalTree($NewDir, @_[ $half .. $#_ ]);
	}
	else {
		# $half was 0, therefore only one filename in @_
		# save the file in the $dir directory
		copyFile($dir, $_[0]);
	}
	
}

# SUBROUTINE &buildClustalWTree:
# Builds a file tree based on a Clustal W phylogeny
# (i.e. a phylogeny derived from a primary-sequence-based
# multiple alignment).  
# Uses these two sub methods (defined below):

sub getEmmaDendogram($@);

sub buildTreeFromDendogram($$@);

sub buildClustalWTree($@) {
    my $dir = shift @_;
	# change working dir so that extra emma output files
	# are stored in the output directory
	chdir $dir;
	Cwd::chdir $dir; #just in case
	buildTreeFromDendogram(getEmmaDendogram($dir,@_),$dir, @_);
}

#### End of treebuilding routines ####

# SUBROUTINE &getEmmaDendogram:
# Calls EMBOSS's "emma" program, which is an interface 
# to Clustal W. 
# Assumes that emma is in the system path.
# Requires the name of an output directory in which to
# save files, and the array of filenames to use
# Creates (1) a fasta format multiple sequence file
# (using EMBOSS's seqret program to read and convert
# the sequences)
# (2) the ClustalW alignment file
# (3) the ClustalW dendogram file
# Returns the name/path to the dendogram file

# The following paramters are used in the emma program
# Comment out all but one option for each variable

$slow = "-slow"; #Use dynamic programming for initial pair-wise alignments
# $slow = "-noslow"; #Use fast approximations for pair-wise alignments

# default, create dendogram and alignment
# $onlyd = "-onlydend"; #only create dendogram, not alignment file

$alignformat = "clustal";
# other format options are fasta, phylip, etc.

sub getEmmaDendogram($@){
	my $dir = shift @_;
	
	# Create a file to contain all sequences
	$multiseqfile = catfile($dir, "allSequences.fa");
	open MULTISEQ, ">", $multiseqfile;
	
	# Read each sequence using seqret and append
	# the results to the multi-sequence file
	foreach $seqfile (@_) {
		print MULTISEQ ">".basename($seqfile)."\n";
		# Use seqret to parse the sequence and output it in raw format
		print MULTISEQ `seqret stockholm::$seqfile stdout -osformat2 raw`;
		print MULTISEQ ("\n\n");
	
	}
	close MULTISEQ;
	
	$dendfile = catfile($dir, "clustalw.dnd");
	$alignfile = catfile($dir, "clustalw.aln");
	
	# test:
	# print $multiseqfile, $alignfile, $dendfile;

	# Call emma
	system "emma", ( "-sequence", "$multiseqfile",
			 "-outseq", "$alignfile",
			 "-dendoutfile", "$dendfile",
			 "-osformat2", $alignformat, 
			 $slow,
			 #$onlyd,
			 "-auto"
			 );
			# add any other parameters, as necessary
	unless (-e $dendfile) { die "Error running emma:\n $?\n--"; }
	
	# Return dendogram filename:
	return $dendfile;
}

# SUBROUTINE &buildTreeFromDendogram
# Takes a dendogram file created by an EMBOSS program
# and uses it to build the corresponding file 
# directory tree.
# Requires (in order) the path of the dendogram file,
# the output directory in which to root the tree,
# and then the list of file names
# The names of the sequences in the dendogram should
# match the names of the files.

sub buildTreeFromDendogram($$@){
	my $dendfile = shift @_;
	my $dir = shift @_;
	
	my $dendstring;
	# open dendogram file and store contents as string
	{
	    open DENDFILE, "<", $dendfile;
	    local $/;  #undefined, i.e. slurp whole file
	    $dendstring = <DENDFILE>;
	    close DENDFILE;
	} # braces to contain "local"
	
	# remove the outside-most brackets and semicolon
	# by replacing them with empty string
	$dendstring =~ s{^\(\n?}{}; # ( then maybe newline at beginning of string
	$dendstring =~ s{\);?\n?$}{}; # ), then maybe ; and newline at end of string
	
	# split the dendogram at ( ) , or :
	# possibly followed by a  newline
	# include the separators but not the  newline in the list
	my @pieces = split /([\)\(,:])\n?/, $dendstring;
	# print @pieces, "\n"; #test
	
	my @dirs; # a stack containing a hierachy of directories
		  # ending in the "current" directory within the tree structure
	my @dircount; # the no. of subdirectories for each directory in @dirs

	my @filecount; # the total no. of profile files contained within each directory in @dirs 

	push @dirs, $dir;
	push @dircount, 0; 
	push @filecount, 0;
	
	PIECE: while (@pieces) {
		my $piece = shift @pieces;
		next PIECE unless (defined($piece) && $piece =~ /\S/ );
		# print "Piece is ", $piece, "\n";
		# print "dirs are ", @dirs, "\n";
		$currentdir = $dirs[$#dirs];
		my $newdir;
		if ($piece eq "(") {
					
			# create a new directory within the current directory
			# directory name is "R" if a sub-directory already exists
			# in current directory, or "L" otherwise
			$newdir = catfile($currentdir, ($dircount[$#dirs]? "R":"L") );
			mkdir ($newdir) 
				or die "Cannot make required output directory $newdir:\n $!\n--";
			
			# update dircount
			$dircount[$#dirs]++;
			
			# push the new directory onto the stacks
			push @dirs, $newdir;
			push @dircount, 0;
			push @filecount, 0;
			
			# print "Created directory $newdir", "\n"; #test
			
			next PIECE;
		}
		elsif ($piece eq ")") {
		    # rename the current directory to indicate the number of files it contains
		    move ($currentdir, $currentdir.$filecount[$#dircount]);
		    # close off the current directory and move up
		    # in the hierarchy by popping the stacks
		    pop @dirs;
		    pop @dircount;
		    pop @filecount;
		    # print "popped\n";
		    next PIECE;
		}
		elsif ($piece eq ",") {
			# Check to see that the current directory
			# does not already have two sub-directories
			# (this happens at the "root" of the tree")

			if ($dircount[$#dirs] == 2) {
				# find the subdirectories within the current directory
				my ($sub1, $sub2) = glob catfile($currentdir, '[R,L]*');

				unless (defined($sub1) && defined($sub2) ) {
					die "Unexpected problems within $currentdir\n--";
				}

				# create a new sub directory, temporarily named "newL"
				my $newSub;
				mkdir ($newSub = catfile($currentdir, "newL")) 
					or die "Cannot make required output directory $newSub:\n $!\n--";

				# move the old sub directories into the new sub
				move($sub1, catfile($newSub, basename($sub1)) )
					or die "Cannot move directory $sub1 to $newSub:\n $!\n--";
				move($sub2, catfile($newSub, basename($sub2)) )
					or die "Cannot move directory $sub1 to $newSub:\n $!\n--";

				# rename the new subdirectory 
				# the number of files in thenew sub is the same as the total in the current dir
				move ($newSub, catfile($currentdir, "L".$filecount[$#dirs] ))
					or die "Cannot rename directory $newSub:\n $!\n--";
				
				# adjust the dircount for the currentdir
				$dircount[$#dirs] = 1; 
					#newSub should be the only direct sub-directory
				# print "rearranged directories\n"
			}
			# otherwise, nothing needs to be done
			next PIECE;
		}
		elsif ($piece eq ":") {
			# discard the following token
			# which represents a branch length
			shift @pieces;
		}
		else { 
			# $piece is a sequence name (base file name)
			# look for a match in the set of filenames in @_
			# create a new directory and copy the file into it
			my $filepath;
			foreach $path (@_) {
				if ($path =~ /$piece/) {
					$filepath = $path;
					last;
				}
			}
			unless (defined($filepath)) {die "Cannot find a file for sequence $piece .\n--";}
			
			mkdir ($newdir = catfile($currentdir, ($dircount[$#dirs] ? "R1":"L1") )) 
				or die "Cannot make required output directory $NewDir:\n $!\n--";
			$dircount[$#dirs]++; # mark that a sub-directory has been added to the current directory
			map {$_++} (@filecount); #increase the filecount for all directories in current hierarchy
			copyFile($newdir, $filepath);
			# print "copied $filepath into $newdir\n";
		}
		#print "\n";
	}
}

# SUBROUTINE &getEnergyTableFromPairs
# Uses a directory of pairwise dynalign results,
# as created by the allpairs.prl script,
# to generate a table of energy values for all the pairwise
# alignments.
#
# Requires as input
# (1) the top directory containing the output from the allpairs.prl script
# (2) the list of full file paths for all the sequences (profile files)
# which were used to generate the pairwise alignments
#
# Returns a two-dimensional table, in the form of a 
# hash in which the values are references to further hashes.
# The keys to both the outer and inner hashes are the full sequence paths
# and the values of the inner hashes are the energy scores from the pairwise
# alignment files.

sub getEnergyTableFromPairs($@){
	my $PairsDir = shift @_;
	
	# create the hash table, with one sub-hash for each
	# entry in the input array of filenames
	my %ETable;
	foreach $seq (@_) {
		$ETable{$seq} = {}; # set the value of each "row" in the table 
				    # to be a reference to an (initially empty) hash
	}
	
	
	opendir TABLE, $PairsDir;
	foreach $level1 (readdir TABLE)  {
		next if ($level1 eq "." || $level1 eq "..");
		$level1 = catfile($PairsDir, $level1);
		next if (!(-d $level1));
		# print "row: ", $level1, "\t";
		
		
		# find the profile in the input files which corresponds 
		# to this directory:
		my $prof1;
		my $seqname = basename($level1);
		foreach $path (@_) {
			if ($path =~ m/$seqname/) {
				$prof1 = $path;
				last;
			}
		}
		unless (defined($prof1)) {die "Cannot find a file for sequence $seqname.\n--";}
		
		opendir ROW, $level1;
		foreach $level2 (readdir ROW) {
			next if ($level2 eq "." || $level2 eq "..");
			$level2 = catfile($level1, $level2);
			next if (!(-d $level2));
			#print $level2, "\t";
		
			
			my $prof2;
			$seqname = basename($level2);
			foreach $path (@_){
				if ($path =~ m/$seqname/) {
					$prof2 = $path;
					last;
				}
			}
			unless (defined($prof2)) {die "Cannot find a file for sequence $seqname.\n--";}
			
			# read the corresponding pairwise alignment profile
			my ($alignment) = glob (catfile($level2, $ProfilePattern));
			open ALIGN, $alignment;
			my ($line, $en);
			do {
				$line = readline(ALIGN);
			} until ((!defined($line))||( ($en) = ($line =~ m/^#=GF\sEN\s([0-9\.\-]+)/)) );
			# i.e. looking for a line of form:
			#	#=(whitespace)EN(whitespace)(some decimal number)(possibly other stuff, like units, here)
			close ALIGN;
			
			unless (defined($line)) {die "Cannot find energy value in profile $alignment.\n--";}
			
			# store the energy value from the pattern-match in the table
			${$ETable{$prof1}}{$prof2} = $en + 0; 
			# also fill the other half of the table:
			${$ETable{$prof2}}{$prof1} = $en + 0;
			
		}
		closedir ROW;
		#print "\n";
	}
	closedir TABLE;
	
	# check: print table
	# print "Energy Table:\n";
	# foreach $first (sort(keys(%ETable))) {
	# 	print $first, "\t";
	# 	my %row = %{$ETable{$first}};
	# 	foreach $second (sort(keys(%row))) {
	# 		print $second, ":\t", $row{$second}, "\t";
	# 	}
	# 	print "\n";
	# }
	
	return %ETable;
}

# SUBROUTINE &NearestNeighbourTree
# Takes as input a hash table created by the 
# &getEnergyTableFromPairs subroutine,
# and uses it to build a tree, such that the pairwise
# scores are used to build the tree in nearest-neigbour fashion
# (i.e. for any sub-tree of the final tree, the average pairwise score 
# for any of the sequences when aligned with each of the others  
# will be less than the score for any of those sequences
# when aligned with a sequence not included in the sub-tree).
#
# returns the root node of a tree, as a reference to a two-element array
# containing the depth 1 nodes/leaves of the tree, with the nodes
# of the same structure as the root, and the leaves being strings containing
# the paths to the original profile files.

sub NearestNeighbourTree(@) {
    my %ETable = @_;


    my ($minE, $p1, $p2);
    my %NodeValues;

# Add all the keys (sequence paths) from the Pairs table
# to the NodeValues hash as both key and value
    map {$NodeValues{$_} = $_; } (keys(%ETable));

    while (scalar(keys(%ETable)) > 2) {
		# find the current minimum energy
		$minE = $BigNumber;
		foreach $prof1 (keys(%ETable)) {
			my %row = %{$ETable{$prof1}};
			foreach $prof2 (keys(%row)) {
			if ($minE > $row{$prof2}){
				$minE = $row{$prof2};
				$p1 = $prof1;
				$p2 = $prof2;
			}
			}
		}

		# create a two-element array of the profiles to align
		# the reference to which is a node in the tree
		# replace string versions of nodes with actual values
		# (References and other scalars are converted to strings when
		# used as keys in the hash)
		my $node = [$NodeValues{$p1}, $NodeValues{$p2}];

		# for every other "row" in the table,
		# replace the "columns" for $p1, $p2
		# with a column for the node, with a value equal to
		# the average of the values for the individual columns
		# also create a new "row" for the node
		# and then delete the rows for the individual elements
		my %newrow;
		foreach $r (keys(%ETable)) {
			if ($r eq $p1 || $r eq $p2) {
				delete $ETable{$r};
				next;
				}
			my $row = $ETable{$r};

			$newrow{$r} = (${$row}{$p1} + ${$row}{$p2})/2;
			delete ${$row}{$p1};
			delete ${$row}{$p2};
			${$row}{$node} = $newrow{$r};
		}
		$ETable{$node} = \%newrow;
		$NodeValues{$node} = $node; #store the node reference without converting it to a string

    }	

# NOTE: at this point, ETable should consist of two "rows"
# (sub-hashes), the keys to which are the final nodes/profiles
# to join to form the tree; form a new node to group together these,
# and replace string versions of nodes with true references:
    
    my $root = [keys(%ETable)];
    map {$_ = $NodeValues{$_};} (@$root);

    return $root;
}

# SUBROUTINE &ProgressiveAlignment 
# Takes as input a hash table created by the 
# &getEnergyTableFromPairs subroutine,
# and uses it to build a tree, such that the best scoring
# sequences are aligned first, and then the remaining sequences
# added one at a time, based on their average score from
# pairwise alignments with each sequence currently in the tree.
#
# returns the root node of a tree, as a reference to a two-element array
# containing the depth 1 nodes/leaves of the tree, with the nodes
# of the same structure as the root, and the leaves being strings containing
# the paths to the original profile files.

sub ProgressiveAlignment(@) {
    my %ETable = @_;


    my ($minE, $p1, $p2);
    my %NodeValues;

    # Add all the keys (sequence paths) from the Pairs table
    # to the NodeValues hash as both key and value
    map {$NodeValues{$_} = $_; } (keys(%ETable));

    # find the current minimum energy
    $minE = $BigNumber;
    foreach $prof1 (keys(%ETable)) {
	my %row = %{$ETable{$prof1}};
	foreach $prof2 (keys(%row)) {
	    if ($minE > $row{$prof2}){
		$minE = $row{$prof2};
		$p1 = $prof1;
		$p2 = $prof2;
	    }
	}
    }

    # create a node (reference to a two-element array) joining
    # these profiles:
    my $node = [$p1, $p2];

    # determine the average energy for any of the other profiles
    # aligned with the two profiles in the node:
    my %avgE;
    {
	my %row1 = %{delete $ETable{$p1}};
	my %row2 = %{delete $ETable{$p2}};
	map { $avgE{$_}=($row1{$_}+$row2{$_})/2; } (keys(%ETable)) ;
    }

    my $newnode;
    my $seqcount = 2;

    while (scalar(keys(%ETable)) > 0) {
	# find the current average minimum energy of another profile 
	# when aligned with the sequences joined by the current node
	$minE = $BigNumber;
	foreach $prof (keys(%ETable)) {
	    if ($minE > $avgE{$prof}){
		$minE = $avgE{$prof};
		$p1 = $prof;
	    }
	}
	
	
	# create a node joining the old node to the new profile
	$newnode = [$p1, $node];
	
	# remove the scores for the added profile from the table:
	my %row = %{delete $ETable{$p1}};
	# calculate a new set of averages:
	my %newavg;
	foreach $prof (keys(%ETable)) {	
	    $newavg{$prof} = ($avgE{$prof}*$seqcount + $row{$prof})/($seqcount + 1);
	}
	%avgE = %newavg;
	++$seqcount;
	$node = $newnode;
	
    }	
    # return the final node
    return $newnode;
}

# SUBROUTINE &printTree
# recursive method to create string descriptions of a tree
# structure based on a root node, as created by the 
# &ProgressiveAlignment or &NearestNeigbour subroutines.
#
# Initial arguments should be the root node,
# followed by the value 0, then two empty strings
# (representing zero indent in the tree,
# and the initial values for the two output strings)
#
# Returns (in order) a string containing a visual depiction of the tree
# and a string containing a bracket-notation dendogram of the tree.
#
sub printTree($$$$);
sub printTree($$$$){
    my ($node, $indent, $pretty, $dendogram) = @_;
   
    my $headers = "";
	if ($indent > 0) {
		$headers = "|".(" "x7);
		$headers x= ($indent-1);
		$headers .= "+".("-"x7);
	}
	if (ref($node)) {
	    $dendogram .= "(\n";
	    ($pretty, $dendogram) = printTree($$node[0], $indent+1, $pretty, $dendogram);
	    $pretty .= ( ("|".(" "x7))x($indent+1)."\n" ) ;
	    $dendogram .= ":1"; #so dendogram can be used by programs which expect a branch-length (not calculated)
	    
	    $dendogram .= ",\n";
	    $pretty .= $headers."*"."\n";

	    $pretty .= ( ("|".(" "x7))x($indent+1)."\n" ) ;
	    ($pretty, $dendogram) = printTree($$node[1], $indent+1, $pretty, $dendogram);
	    $dendogram .= ":1)\n";
	}
	else {
	    $pretty .= $headers."Sequence:".basename($node)."\n";
	    $dendogram .= basename($node);
	}
	return ($pretty, $dendogram);
}

# SUBROUTINE &traverseTree:
# Call the HDynalign routine on a file tree
# Requires a directory name as input
# plus the list of arguments to send to the HDynalign program
# If the directory contains a profile file,
# returns the name of that file
# Else, if the directory has subdirectories, recursively
# calls this method on two of them,
# then executes the hdynalign program on the two
# returned filenames, and itself returns the 
# name of the profile file created by the program.
#
# If there are more than one profile files, or
# more than two subdirectories, just takes the first it gets
# If there is no profile and less than two sub-directories,
# throws an error.

sub traverseTree ($@); #pre-declare for recursion

sub callHD($$$$@);

sub traverseTree ($@) {

	my $dir = shift @_;
	
	# look for a profile-format file
	if (my ($profile) = glob (catfile($dir, $ProfilePattern)) ){
		return $profile;
	}
	
	# (else) read entries from directory until two sub-directories found
	my @subdir;
	opendir DIR, $dir;
	while (my $x = (readdir DIR) ) {
		next if ($x eq "." || $x eq "..");
		$x = catfile($dir, $x);
		$subdir[@subdir] = $x if (-d $x );
		
		last if (@subdir == 2);
	}
	closedir DIR;
	
	# exit if didn't find sub-directories
	if (@subdir != 2) { die "$dir does not have expected tree structure\n--";}
	
	# get profiles from sub-trees:
	my $a = traverseTree($subdir[0], @_);
	my $b = traverseTree($subdir[1], @_);

	
	# create a name for an output profile file
	# based on the name of the directory
	my $outfile;
	$outfile = basename($dir).$ProfileSuffix;
	$outfile = catfile($dir, $outfile);
		
	return callHD($a, $b, $dir, $outfile, @_);
}

# SUBROUTINE &callHD:
# Runs the hdynalign program
# Requires as input (in order):
# the full paths of the two input profiles,
# the path of an output directory in which to save ct files,
# the path and name of an output file for the alignment profile
# the switches and values to set hdynalign parameters
#
# determines the path of the program to call based on the %HDPath hash
# and the executing computer's hostname (from the Sys::Hostname::hostname() function)

sub callHD($$$$@)  {
	my ($infile1, $infile2, $outdir, $outfile, %hdopts) = @_;
	
	# change working directory so that hdynalign can find its data files
	# chdir $HDDir;
	# Cwd::chdir $HDDir; #just in case
	# print "Current directory: ", getcwd();
	
	my $HDBinary = $HDPath{ hostname() } || $HDPath{ "default" };
	
	print "Aligning $infile1 \n with $infile2 \n using $HDBinary\n\n";
	
	# open hdynalign as a process to pipe to
	open (HD, "| ".$HDBinary)
			or die "Cannot pipe to hdynalign:\n $! $?\n--";
			
	# Pipe hdynalign parameters in
	# and echo them to standard out
	# NOTE: hdynalign asks for parameters interactively in this order:
	
	# Enter the name of the first profile:
	print HD $infile1, "\n"; 
	# print "First profile: ", $infile1, "\n"; 
	
	# Enter the name of the second profile:
	print HD $infile2, "\n";
	# print "Second profile: ", $infile2, "\n"; 
	
	# Enter the max separation:
	my $maxsep = ( defined($hdopts{$maxSepOption}) ? $hdopts{$maxSepOption} : $defaultMaxSep );
	print HD $maxsep, "\n";
	# print "Maximum separation: ", $maxsep, "\n"; 
	
	# Enter the gap penalty:
	my $gap = (defined($hdopts{$gapPenaltyOption}) ? $hdopts{$gapPenaltyOption} : $defaultGapPenalty);
	print HD $gap, "\n";
	# print "Gap penalty", $gap, "\n"; 
	
	# Allow single BP inserts into one sequence? (1/0)
	my $insert = ( defined($hdopts{$BPInsertionOption}) ? $hdopts{$BPInsertionOption} : $defaultBPInsertionState );
	print HD $insert, "\n";
	# print "BP inserts? ", ($insert ? "true" : "false"), "\n"; 
	
	# Enter the directory for output ct files of the first profile:
	print HD $outdir, "\n";
	# print "First output directory: ", $outdir, "\n"; 
	
	# Enter the directory for output ct files of the second profile:
	print HD $outdir, "\n";
	# print "Second output directory: ", $outdir, "\n"; 
	
	# Enter the file name for the output of the alignment:
	print HD $outfile, "\n";
	# print "Output profile: ",$outfile, "\n"; 
	
	# Close (i.e. wait for) hdynalign
	close HD; # or die "Error with HDynalign:\n Pipe Error: $! \n Child Error: $? \n--";
	
	print "**********************\n\n";
	
	return $outfile;
}

#Succesful completion when importing module
1;

