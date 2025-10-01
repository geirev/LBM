#!/usr/bin/perl

use strict;
use warnings;

# Get all Fortran source files
my @files = glob("*.F *.F90");

foreach my $file (@files) {
    # Create object file name
    (my $objfile = $file) =~ s/\..*$/\.o/;

    # Store found module dependencies
    my %dependencies;

    open my $fh, '<', $file or die "Cannot open file $file: $!\n";

    while (my $line = <$fh>) {
        chomp $line;

        # Look for 'use' statements (case-insensitive)
        if ($line =~ /^\s*use\s+([a-z0-9_]+)/i) {
            my $module = $1;
            $dependencies{"$module.o"} = 1;
        }
    }

    close $fh;

    # Print module dependencies
    foreach my $dep (sort keys %dependencies) {
        print "$objfile:      $dep\n";
    }

    # Always depend on the makefile
    print "$objfile:      makefile\n";
}
