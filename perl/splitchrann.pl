#!/usr/bin/perl
use strict;
use warnings;

my ($infasta, $outfasta, $inputann, $outputann, $annformat, $breaks) = @ARGV;

my %breaks = ();
foreach my $b (split(",",$breaks)) {
	if ($b =~ /(\S+):(\d+)/){
		$breaks{$1} = $2;
	}
	else {
		print "Supply breaks in format: [chr:coordintate,chr:coordinate]\n";
		exit(1);
	}
}


#===|===
#123|456 #old, gff/gtf
#123|123 #new, gff/gtf
#012|345 #old, bed
#012|012 #new, bed

my $seqid = "";
my $description = "";
my $sequence = "";
my %brokenchr = ();
my @characters = ("A".."ZZZ"); ##limitations for more than 18278 fragments

if ($inputann =~ /.*gz$/) {
	open (IA, "zcat $inputann |") or die $!;
}
else {
	open (IA, "<$inputann") or die $!;
}
if ($outputann =~ /.*gz$/) {
	open (OA, " | gzip > $outputann") or die $!;
}
else {
	open (OA, ">$outputann") or die $!;
}

while (my $ann = <IA>){
	print OA $ann if (substr($ann,0,1) eq "#"); ## print comment lines as is
	my @a = split ("\t", $ann); ## annotations should be tab separated	
	if (exists $breaks{$a[0]}) {
		my $end = -1;
		my $start = -1;
		if ($annformat =~ /bed/i) {
			$start = $a[1];
			$end = $a[2];
		}
		elsif ($annformat =~ /g[tf]f/i){
			$start = $a[3];
			$end = $a[4];
		}
		if ($end <= $breaks{$a[0]}) {
			print OA $ann;
		}
		else {
			$start = $start - $breaks{$a[0]};
			$end = $end - $breaks{$a[0]};
			if ($annformat =~ /bed/i) {
				$a[0] = "$a[0].B"; ## this need to change if multiple splits on one chromosome
				$a[1] = $start;
				$a[2] = $end;	
			}
			else {
				$a[0] = "$a[0].B";
				$a[3] = $start;
				$a[4] = $end;
			}
			print OA join("\t", @a);
		}
	}	
	else {
		print OA $ann;
	}
}
close IA;
close OA;

if ($infasta =~ /.*gz$/){
	open (F, "zcat $infasta |") or die $!;
}
else{
	open (F, "<$infasta") or die $!;
}
if ($outfasta =~ /.*gz$/){
	open (O, "| gzip >$outfasta") or die $!;
}
else {
	open (O, ">$outfasta") or die $!;
}

while (my $line = <F>) {
	chomp $line;
	if ($line =~ />(\S+)(.*)/){
		if ($sequence && length($sequence) > 0) {
			if (exists $breaks{$seqid}) {
				print O ">$seqid".($description ? $description : "")."\n";
				my $tmpseq = substr($sequence, 0, $breaks{$seqid});
				for (my $i=0; $i<length($tmpseq); $i+=60) { 
					print O substr($tmpseq, $i, 60) ."\n";
				}
				print O ">$seqid.B".($description ? $description : "")."\n";
				$tmpseq = substr($sequence, $breaks{$seqid});
				for (my $i=0; $i<length($tmpseq); $i+=60) { 
					print O substr($tmpseq, $i, 60) ."\n";
				}
			}
			else {
				print O ">$seqid".($description ? $description : "")."\n";
				for (my $i=0; $i<length($sequence); $i+=60) { 
					print O substr($sequence, $i, 60) ."\n";
				}
			}
		}
		$seqid = $1;
		$description = $2;
		$sequence = "";
	}
	else {
		$sequence = $line;
	}
}
if ($sequence && length($sequence) > 0) {
	if (exists $breaks{$seqid}) {
		print O ">$seqid".($description ? $description : "")."\n";
		my $tmpseq = substr($sequence, 0, $breaks{$seqid});
		for (my $i=0; $i<length($tmpseq); $i+=60) { 
			print O substr($tmpseq, $i, 60) ."\n";
		}
		print O ">$seqid.B".($description ? $description : "")."\n";
		$tmpseq = substr($sequence, $breaks{$seqid});
		for (my $i=0; $i<length($tmpseq); $i+=60) { 
			print O substr($tmpseq, $i, 60) ."\n";
		}
	}
	else {
		print O ">$seqid".($description ? $description : "")."\n";
		for (my $i=0; $i<length($sequence); $i+=60) { 
			print O substr($sequence, $i, 60) ."\n";
		}
	}
}
close F;
close O;

exit;

