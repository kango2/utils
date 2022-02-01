#!/usr/bin/perl
use strict;
use warnings;

my ($infasta, $outfasta, $maskids) = @ARGV;

open (F, "<$infasta") or die $!;
open (O, ">$outfasta") or die $!;

my %masks = ();

open (M, "<$maskids") or die $!;
while (<M>){
	chomp $_;
	$masks{$_}="";
}
close M;

my $seqid = "";

while (my $line = <F>) {
	if ($line =~ />(\S+)/){
		$seqid = $1;
		print O $line;
	}
	else {
		if (! exists $masks{$seqid}) {
			print O $line;
		}
		else {
			$line =~ tr/[A-Z]/N/;
			$line =~ tr/[a-z]/n/;
			print O $line;
		}
	}
}
close F;
close O;
