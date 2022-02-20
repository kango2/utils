#!/usr/bin/perl
use strict;
use warnings;

my ($infasta, $outfasta, $inputann, $outputann, $annformat, $maxchrsize) = @ARGV;

open (F, "<$infasta") or die $!;
open (O, ">$outfasta") or die $!;

my $seqid = "";
my $description = "";
my $sequence = "";
my %brokenchr = ();
my @characters = ("A".."ZZZ"); ##limitations for more than 18278 fragments

if ($inputann =~ /*.gz$/) {
	open (IA, "zcat $inputann |") or die $!;
}
else {
	open (IA, "<$inputann") or die $!;
}
if ($outputann =~ /*.gz$/) {
	open (OA, " | gzip > $outputann") or die $!;
}
else {
	open (OA, ">$outputann") or die $!;
}

while (my $ann = <IA>){
	print OA $ann if (substr($ann,0,1) eq "#"); ## print comment lines as is
	my @a = split ("\t", $ann); ## annotations should be tab separated
	my $end = -1;
	my $start = -1;
	if ($annformat =~ /bed/i) {
		$start = $a[1];
		$end = $a[2];
	}
	elsif ($annformat =~ /gtf/i){
		$start = $a[3];
		$end = $a[4];
	}
	elsif ($annformat =~ /gff/i){
		$start = $a[3];
		$end = $a[4];
	}
	
	if ($end < $maxchrsize) {
		print OA $ann;
	}
	else {
		$brokenchr{$a[0]} = $start;
		$
	}
	
}


while (my $line = <F>) {
	chomp $line;
	if ($line =~ />(\S+)(.*)/){
		if ($sequence && length($sequence) > 0) {
			if (length($sequence) > $maxchrsize) {
				$brokenchr{$seqid} = "";
				##split into two with A, B, C and so on qualifiers
				for (my $j = 0 ; $j < $maxchrsize)
				print O ">$seqid".($description ? $description : "")."\n";
				for (my $i=0; $i<length($sequence); $i+=60) { 
					print O substr($sequence, $i, 60) ."\n"
				}
			
			}
			else{
				print O ">$seqid".($description ? $description : "")."\n";
				for (my $i=0; $i<length($sequence); $i+=60) { 
					print O substr($sequence, $i, 60) ."\n"
				}
			}

		}
		$seqid = $1;
		$description = $2;
		$sequence = "";
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
