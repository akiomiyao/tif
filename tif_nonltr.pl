#!/usr/bin/perl
#
# tif_nonltr.pl
#
# Copyright(C) 2020 MIYAO Akio. All Rights Reserved.
#
# Use of non-profit research purpose is permitted.

$target = $ARGV[0]; 
$blastdb = $ARGV[1];
$head = $ARGV[2];

if ($head eq ""){
    print "Usage:
  perl tif_nonltr.pl target blastdb head_sequence

  perl tif_nonltr.pl SRR835080 dmel626.fasta GAACCCTCTGTCGTAGAACACTACTA

  Before run this script, installing of NCBI BLAST+ and sratoolkit (fastq-dump) is required.
  URLs for download files is
  https://blast.ncbi.nlm.nih.gov/Blast.cgi?CMD=Web&PAGE_TYPE=BlastDocs&DOC_TYPE=Download
  https://trace.ncbi.nlm.nih.gov/Traces/sra/sra.cgi?view=software

  Save fasta file of reference genome into the tif directory. 

  If target is assesson number, directory will be created and then
  fastq files will be downloaded into target/read directory.

  Result will be saved to 'target/result' file.  
";
    exit;
}

if (! -e $blastdb){
    print "
  $blastdb NOT found.
  Save $blastdb in the tif directory and run again.

";
    exit;
}

if (! -e "$blastdb.nsq"){
    system("makeblastdb -in $blastdb -dbtype nucl");
}

if (! -e $target){
    system("mkdir $target $target/read") if ! -e "$target/read";
    open(REPORT, "> $target/log");
    print REPORT "Target: $target
Reference: $blastdb
Head sequence: $head

";
    close(REPORT);
    report("Job begin.");
    report("Downloading fastq files.");
    system("cd $target/read && fastq-dump --split-files $target");
    report("Downloading fastq files. Done.");
}else{
    open(REPORT, "> $target/log");
    print REPORT "Target: $target
Reference: $blastdb
Head sequence: $head

";
    close(REPORT);
    report("Job begin.");
}

$file_list = "$target/read/*";

report("Searching flanking sequences of $head.");
open(IN, "cat $file_list |grep $head|");
open(OUT, "> $target/upstream");
while(<IN>){
    $pos = index($_, $head);
    $upstream = substr($_, 0, $pos);
    $upstream =~ s/A*$//;
    $upstream =~ s/T*$//;
    
    if (length($upstream) > 20){
	$i++;
	print OUT ">$i
$upstream\n"
    }
}
close(IN);
close(OUT);

report("Searching flanking sequences of polyA.");
$i = 0;
open(IN, "cat $file_list |grep AAAAAAAAAAAAAAA|");
open(OUT, "> $target/downstream");
while(<IN>){
    next if $_ !~/^[ACGT]/;
    chomp;
    $pos = index($_, "AAAAAAAAAAAAAAA");
    $downstream = substr($_, $pos + 15);
    $downstream =~ s/^A*//;
    $tmp = $downstream;
    $tmp =~ y/[ACGT]//d;
    if (length($downstream) > 20 and $tmp eq ""){
	$i++;
	print OUT">$i
$downstream\n"
    }
}
close(IN);
close(OUT);


open(OUT, "|sort -k 1 -k 2 -n |uniq > $target/blast");
open(IN, "blastn -db $blastdb -outfmt 6 -num_alignments 1 < $target/upstream|");
while(<IN>){
    @row = split;
    print OUT "$row[1]\t$row[9]\thead\n";
}
open(IN, "blastn -db $blastdb -outfmt 6 -num_alignments 1 < $target/downstream|");
while(<IN>){
    @row = split;
    print OUT "$row[1]\t$row[8]\ttail\n";
}
close(IN);
close(OUT);

open(IN, "$target/blast");
open(OUT, "> $target/result");
while(<IN>){
    chomp;
    @row = split;
    if ($pchr eq $row[0] and $row[2] ne $prev and abs($ppos - $row[1]) < 100){
	print OUT "$pline
$_\n";
	print "$pline
$_\n";
    }
    $prev = $row[2];
    $ppos = $row[1];
    $pchr = $row[0];
    $pline = $_;
}

if (-s "$target/result" == 0){
    report("No transposition found.");
}

sub report{
    my $message = shift;
    my $now = `date`;
    chomp($now);
    $message = "$now : $message\n";
    open(REPORT, ">> $target/log");
    flock(REPORT,2);
    print $message;
    print REPORT $message;
    close(REPORT);
}
