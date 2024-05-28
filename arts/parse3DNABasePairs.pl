#!/usr/bin/perl

##########################################################################
#
# GOAL: Convert the output of the 3DNA find_apir utility to the
# following simple format that describe base pairs:
# <chain id> <reside id> <chain id> <residue id>
#
# INPUT:  An output file of the 3DNA find_pair utlity
#
# OUTPUT: The following simple format that describe base pairs:
# <chain id> <reside id> <chain id> <residue id>
# e.g.: R 1 72 R
#       R 2 71 R
#       R 3 70 R
#       R 4 69 R
#       R 5 68 R
#
# AUTHOR:
# Oranit Dror (contact: oranit@post.tau.ac.il)
#
##########################################################################

# Find the directory where this perl script is located
use FindBin;

use strict;


sub printUsage {
  print "Convert the output of the 3DNA find_apir utility to the following simple format that describe base pairs: \n";
  print "<chain id> <reside id> <chain id> <residue id> \n \n";
  print "Contact: oranit\@post.tau.ac.il \n";
  print "Usage: $FindBin::Script <3DNA output file> \n";
  print "Usage: $FindBin::Script --help \n";
}


# Check the arguments
if ($ARGV[0] =~ "--help") {
    printUsage();
    exit(0);
}

if( $#ARGV < 0) {
    printUsage();
    exit(0);
}

open (DATA, $ARGV[0]) or die "Cannot open input file $ARGV[0]";

print "# base-pairs for $ARGV[0]:\n";
print "# FORMAT: <chain>:<residue>-<chain>:<residue>\n";

while(<DATA>){
  if( /#[ ]*[0-9]* [|x+] ([-A-Z0-9]):[.]*([0-9]+)_:[^:]*:[.]*([0-9]+)_:([-A-Z0-9])/) {
    my $chain1 = $1;
    my $residue1 = $2;
    my $residue2 = $3;
    my $chain2 = $4;

    if ($chain1 eq '-') {
      $chain1 = ' ';
    }

    if ($chain2 eq '-') {
      $chain2 = ' ';
    }

    print "$chain1:$residue1-$chain2:$residue2\n";
  }

  next;
}
