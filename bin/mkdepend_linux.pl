#!/usr/bin/perl -w

#Extract module dependencies
foreach $file (<*.F *.F90>) {

  $objfile=$file;
  $objfile=~s/\..*$/\.o/;

  # Module dependencies
  open PIPE,"cat $file |" or die "Cannot open pipe \n";
  while (<PIPE>) {
     chop;
     if (/^[ ]*use /i) {
        s/^[ ]*use[ ]*//i;
        s/[^a-z0-9_].*$//i;
        {
           print "$objfile:      $_.o\n";
        }
     }
  }
  close PIPE;

  #makefile dependencies (all files)
  print "$objfile:      makefile\n";
}
