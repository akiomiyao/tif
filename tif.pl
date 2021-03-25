#!/usr/bin/perl
#
#  Copyright (c) 2012 - 2015
#       National Institute of Agrobiological Sciences.  All rights reserved.
#  Copyright (c) 2016
#       National Agriculture and Food Research Organization.  All rights reserved.
#
#       
# This code is derived from software contributed to NIAS and NARO by
# Akio Miyao.
#
# Redistribution and use in source forms, with or without modification,
# are permitted provided that the following conditions
# are met:
# 1. Redistributions of source code must retain the above copyright
#    notice, this list of conditions and the following disclaimer.
#
# THIS SOFTWARE IS PROVIDED BY THE REGENTS AND CONTRIBUTORS ``AS IS'' AND
# ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
# IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
# ARE DISCLAIMED.  IN NO EVENT SHALL THE REGENTS OR CONTRIBUTORS BE LIABLE
# FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
# DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS
# OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION)
# HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
# LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY
# OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF
# SUCH DAMAGE.
#
#
#

$usage = "    $0 - Transposon Insertion Finder

    perl $0 reference.fasta target head_sequence tail_sequence (max_process_number)

    For endogenous retrotransposon Tos17 in rice,
    perl $0 IRGSP-1.0_genome.fasta.gz ttm5 TGTTAAATATATATACA TTGCAAGTTAGTTAAGA

    For P-element of Drosophila malanogaster
    perl $0 dmel-all-chromosome-r6.26.fasta target CATGATGAAATAACAT ATGTTATTTCATCATG

    Example
    In the directory of $0 script,
    mkdir target
    mkdir target/read
    cp somewhere/fastq_files target/read
    cp somewhere/reference_fasta .
    perl $0 reference_fasta target head_sequence tail_sequence

    Result and log are saved in target directory.

    Compressed read and genome files with gz, bz2 and xz extentions
    will be analyzed without pre-decompression. 

    reference.fasta is the reference sequence with multi-fasta format.
    head_sequence is a short sequence at 5'-end of target transposon.
    tail_sequence is a short sequence at 3'-end of target transposon.
    Length of head_seq and short_seq shold be from 17 to 21 bp.
    max_process_number is optional. Default number of maximun process
    is number of CPUs.
    

Author: Akio Miyao <miyao\@affrc.go.jp>

";

$ref        = $ARGV[0];
$target     = $ARGV[1];
$head       = $ARGV[2];
$tail       = $ARGV[3];
$maxprocess = $ARGV[4];

if ($tail eq ""){
    print $usage;
    exit;
}

if ($ref =~ /:/){
    ($sub, $name) = split(':', $ref);
    system("touch $target/child/$sub.$head.$tail.$name");
    &$sub();
    system("rm $target/child/$sub.$head.$tail.$name");
    exit;
}

$uname = `uname`;
chomp($uname);
if ($uname eq "Linux"){
    open(IN, "/proc/cpuinfo");
    while(<IN>){
        $maxprocess ++  if /processor/;
    }
    close(IN);
}elsif ($uname eq "FreeBSD"){
    open(IN, "sysctl kern.smp.cpus |");
    while(<IN>){
	chomp;
	$maxprocess = (split(': ', $_))[1];
    }
    close(IN);
}elsif($uname eq "Darwin"){
    $zcat = "gzcat";
    open(IN, "sysctl hw.logicalcpu |");
    while(<IN>){
	chomp;
	$maxprocess = (split(': ', $_))[1];
    }
    close(IN);
}
$maxprocess = 4 if $maxprocess eq "";

$start = time();
open(LOG, "> $target/log.$head.$tail");
print LOG "Program: $0
Reference: $ref
Target: $target
Head: $head
Tail: $tail
Max process: $maxprocess

";
&log("TIF start.");

if (-d "$target/child"){
    system("rm -r $target/child");
}

system("mkdir $target/child");

$rhead = complement($head);
$rtail = complement($tail);

if (! -e "$target/selected.$head.$tail"){
    opendir(DIR, "$target/read");
    foreach $file (sort readdir(DIR)){
	next if $file =~ /^\./;
	if ($file =~ /.gz$/){
	    $catcmd = "zcat";
	}elsif($file =~ /.bz2$/){
	    $catcmd = "bzcat";
	}elsif($file =~ /.xz$/){
	    $catcmd = "xzcat";
	}else{
	    $catcmd = "cat";
	}
	&waitFork;
	&log("searching $head in $file");
	system("touch $target/child/$head.$file && $catcmd $target/read/$file | grep $head > $target/tmp.$file.$head && rm $target/child/$head.$file &");
	&waitFork;
	&log("searching $tail in $file");
	system("touch $target/child/$tail.$file && $catcmd $target/read/$file | grep $tail > $target/tmp.$file.$tail && rm $target/child/$tail.$file &");
	&waitFork;
	&log("searching $rhead in $file");
	system("touch $target/child/$rhead.$file && $catcmd $target/read/$file | grep $rhead > $target/tmp.$file.$rhead && rm $target/child/$rhead.$file &");
	&waitFork;
	&log("searching $rtail in $file");
	system("touch $target/child/$rtail.$file && $catcmd $target/read/$file | grep $rtail > $target/tmp.$file.$rtail && rm $target/child/$rtail.$file &");
    }
    &join;
    system("cat $target/tmp.* > $target/grep.$head.$tail && rm $target/tmp.*");
}

&log("Reading chromosome sequences.");

if (-e "chr"){
    system("rm -r chr");
}
system("mkdir chr");

if ($ref =~ /.gz$/){
    open(IN, "zcat $ref |");
}elsif($ref =~ /.bz2$/){
    open(IN, "bzcat $ref |");
}elsif($ref =~ /.xz$/){
    open(IN, "xzcat $ref |");
}else{
    open(IN, $ref);
}
while(<IN>){
    chomp;
    if (/^>/){
	$pos += 0;
	if ($seq eq ""){
	    $pos = 0;
	}else{
	    while(1){
		$filename = "000000000" . $pos;
		$filename = substr($filename, length($filename) - 10, 10);
		$filename = $name . "-" . $filename;
		if ($pos == 0){
		    push(@chr, $filename);
		    for($i = 0; $i < 50; $i++){
			$fragment .= "NNNNNNNNNN";
		    }
		    $fragment .= substr($seq, $pos, 1000500);
		    open(OUT, "> chr/$filename");
		    print OUT $fragment;
		    close(OUT);
		}else{
		    $fragment = substr($seq, $pos - 500, 1001000);
		    if ($fragment eq ""){
			$seq = "";
			$pos = 0;
			last;
		    }
		    push(@chr, $filename);
		    open(OUT, "> chr/$filename");
		    print OUT $fragment;
		    close(OUT);
		}
		$pos += 1000000;
	    }
	}
	s/^>//;
        $name = (split)[0];
    }else{
	y/a-z/A-Z/;
	y/ACGT/N/c;
	$seq .= $_;
    }
}
close(IN);
while(1){
    $filename = "000000000" . $pos;
    $filename = substr($filename, length($filename) - 10, 10);
    $filename = $name . "-" . $filename;
    if ($pos == 0){
	push(@chr, $filename);
	for($i = 0; $i < 50; $i++){
	    $fragment .= "NNNNNNNNNN";
	}
	$fragment .= substr($seq, $pos, 1000500);
	open(OUT, "> chr/$filename");
	print OUT $fragment;
	close(OUT);
    }else{
	$fragment = substr($seq, $pos - 500, 1001000);
	if ($fragment eq ""){
	    $seq = "";
	    $pos = 0;
	    last;
	}
	push(@chr, $filename);
	open(OUT, "> chr/$filename");
	print OUT $fragment;
	close(OUT);
    }
    $pos += 1000000;
}

foreach $name (@chr){
    &waitFork;
    &log("Selecting reads in $name.");
    system("perl tif.pl select:$name $target $head $tail &");
}

&join;

open(IN, "cat $target/tmp.$head.$tail.* |");
while(<IN>){
    chomp;
    @row = split;
    $count{$row[0]} += $row[1];
}
close(IN);

open(OUT, "> $target/selected.$head.$tail");
foreach (sort keys %count){
    if ($count{$_} == 0){
	print OUT "$_\n";
    }
}
close(OUT);
system("rm $target/tmp.$head.$tail.*");

foreach $name (@chr){
    &waitFork;
    &log("Mapping of junction in $name.");
    system("perl tif.pl map:$name $target $head $tail &");
}

&join;

open(IN, "cat $target/map.$head.$tail.* |");
while(<IN>){
    chomp;
    @row = split('\t', $_);
    if ($row[0] eq "head"){
	if (length($row[4]) > length($maphead{$row[1]}{$row[2]}{$row[3]})){
	    $maphead{$row[1]}{$row[2]}{$row[3]} = $row[4];
	    $mhc{$row[1]}{$row[2]} ++;
	}
    }else{
	if (length($row[4]) > length($maptail{$row[1]}{$row[2]}{$row[3]})){
	    $maptail{$row[1]}{$row[2]}{$row[3]} = $row[4];
	    $mtc{$row[1]}{$row[2]} ++;
	}
    }
}
close(IN);
system("rm $target/map.$head.$tail.*");

open(OUT, "> $target/result.$head.$tail");
foreach $chr (sort keys %maphead){
    foreach $pos (sort bynumber keys %{$maphead{$chr}}){
        foreach $direction (sort keys %{$maphead{$chr}{$pos}}){
            $upstream = $maphead{$chr}{$pos}{$direction};
            for($i = $pos - 20 ; $i <= $pos + 20; $i++){
                if ($maptail{$chr}{$i}{$direction} ne ""){
                    $downstream = $maptail{$chr}{$i}{$direction};
                    $tsd_size = abs($pos - $i) + 1;
                    $tsd_head =  substr($upstream, length($upstream) - $tad_size, $tsd_size);
                    $tsd_tail =  substr($downstream, 0, $tsd_size);
                    print STDERR "$chr\t$pos\t$i\t$tsd_size\t$tsd_head\t$tsd_tail\t$direction\t$upstream\t$downstream\t$mhc{$chr}{$pos}\t$mtc{$chr}{$i}\n";
		    print OUT "$chr\t$pos\t$i\t$tsd_size\t$tsd_head\t$tsd_tail\t$direction\t$upstream\t$downstream\t$mhc{$chr}{$pos}\t$mtc{$chr}{$i}\n";
               }
            }
        }
    }
}
close(OUT);

$end = time();

&log("TIF end.");

$elapsed = $end - $start;
$hour = int($elapsed / 3600);
$min = $elapsed % 3600;
$sec = $min % 60;
$min = int($min / 60);
if ($hour >= 24){
    $day = int($hour / 24);
    $hour = $hour % 24;
}
if ($day > 1){
    $etime .= "$day days ";
}elsif($day == 1){
    $etime .= "$day day ";
}
if ($hour > 1){
    $etime .= "$hour hours ";
}elsif($hour == 1){
    $etime .= "$hour hour ";
}
if ($min > 1){
    $etime .= "$min minutes ";
}elsif($min == 1){
    $etime .= "$min minute ";
}
if ($sec > 1){
    $etime .= "$sec seconds ";
}elsif($sec == 1){
        $etime .= "$sec second ";
}
print STDERR "$etime ($elapsed seconds) elapsed.\n";
print LOG    "$etime ($elapsed seconds) elapsed.\n";
close(LOG);

sub select{
    $rhead = complement($head);
    $rtail = complement($tail);
    open(IN, "$target/grep.$head.$tail");
    while(<IN>){
	chomp;
	$read{$_} = 0;
    }
    close(IN);
    open(CHR, "chr/$name");
    my $seq = <CHR>;
    close(CHR);
    foreach $read (sort keys %read){
	foreach ($head, $tail, $rhead, $rtail){
	    $pos = index($read, $_);
	    if ($pos > -1){
		$tmp = substr($read, $pos - 20, 40);
		last;
	    }
	}
	$pos = index($seq, $tmp, 0);
	if ($pos != -1){
	    $read{$read} = 1;
	}
	$comp = complement($tmp);
	$pos = index($seq, $comp, 0);
	if ($pos != -1){
	    $read{$read} = 1;
	}
    }
    open(OUT, "> $target/tmp.$head.$tail.$name");
    foreach (sort keys %read){
	print OUT "$_\t$read{$_}\n";
    }
    close(OUT);
}

sub map{
    open(IN, "$target/selected.$head.$tail");
    while(<IN>){
	chomp;
	if (/$head/){
	    &addHead($_);
	}elsif(/$tail/){
	    &addTail($_);
	}
	$comp = &complement($_);
	if ($comp =~ /$head/){
	    &addHead($comp);
	}elsif($comp =~ /$tail/){
	    &addTail($comp);
	}
    }
    close(IN);
    
    open(OUT, "> $target/map.$head.$tail.$name");
    open(CHR, "chr/$name");
    my $seq = <CHR>;
    close(CHR);
    my ($chr, $bp) = split('-', $name);
    $bp += 0;
    foreach $upstream (sort keys %head){
	$junction = substr($upstream, length($upstream) - 20, 20);
	$rjunction = complement($junction);
	while (1) {
	    $pos = index($seq, $junction, $pos + 1);
	    if ($pos > 500){
		$length = length($upstream);
		$ref = substr($seq, $pos - $length + 20, $length);
		if ($ref eq $upstream){
		    $tpos = $bp + $pos + 20 - 500;
		    if ($length > length($maphead{$chr}{$tpos}{forward})){
			print OUT "head\t$chr\t$tpos\tforward\t$upstream\n";
		    }
		}
	    }elsif($pos == -1){
		last;
	    }
	}
	while(1){
	    $rpos = index($seq, $rjunction, $rpos + 1);
	    if ($rpos > 500){
		$cupstream = complement($upstream);
		$length = length($cupstream);
		$ref = substr($seq, $rpos, $length);
		if ($ref eq $cupstream){
		    $tpos = $bp + $rpos + 1 - 500;
		    if ($length > length($maphead{$chr}{$tpos}{reverse})){
			print OUT "head\t$chr\t$tpos\treverse\t$upstream\n";
		    }
		}
	    }elsif($rpos == -1){
		last;
	    }
	}
    }
    
    foreach $downstream (sort keys %tail){
	$junction = substr($downstream, 0, 20);
	$rjunction = complement($junction);
	while (1) {
	    $pos = index($seq, $junction, $pos + 1);
	    if ($pos > 500){
		$length = length($downstream);
		$ref = substr($seq, $pos, $length);
		if ($ref eq $downstream){
		    $tpos = $bp + $pos + 1 - 500;
		    if ($length > length($maptail{$chr}{$tpos}{forward})){
			print OUT "tail\t$chr\t$tpos\tforward\t$downstream\n";
		    }
		}
	    }elsif($pos == -1){
		last;
	    }
	}
	while(1){	
	    $rpos = index($seq, $rjunction, $rpos + 1);
	    if ($rpos > 500){
		$cdownstream = complement($downstream);
		$length = length($cdownstream);
		$ref = substr($seq, $rpos - $length + 20, $length);
		if ($ref eq $cdownstream){
		    $tpos = $bp + $rpos + 20 - 500;
		    if ($length > length($maptail{$chr}{$tpos}{reverse})){
			print OUT "tail\t$chr\t$tpos\treverse\t$downstream\n";
		    }
		}
	    }elsif($rpos == -1){
		last;
	    }
	}
    }
    close(OUT);
}

sub addHead{
    my $seq = shift;
    $pos = index($seq, $head);
    $upstream = substr($seq, 0, $pos);
    if (length($upstream) > 20){
	$head{$upstream} = 1;
    }
}

sub addTail{
    my $seq = shift;
    $pos = index($seq, $tail);
    $downstream = substr($seq, $pos + length($tail));
    if (length($downstream) > 20){
	$tail{$downstream} = 1;
    }
}

sub complement{
    my $seq = shift;
    my @seq = split('', $seq);
    my $i = 0;
    my $out = "";
    for ($i = length($seq) - 1 ; $i >= 0; $i--){
        if ($seq[$i] eq "A"){
            $out .= "T";
        }elsif($seq[$i] eq "C"){
            $out .= "G";
        }elsif($seq[$i] eq "G"){
            $out .= "C";
        }elsif($seq[$i] eq "T"){
            $out .= "A";
        }
    }
    return $out;
}

sub waitFork{
    my $count;
    while(1){
	select undef, undef, undef, 0.1;
	$count = 0;
	opendir(DIR, "$target/child");
	foreach(sort readdir(DIR)){
	    if ($_ !~ /^\./){
		$count++;
	    }
	}
	return if $count < $maxprocess;
    }
}

sub join{
    my $count;
    while(1){
	sleep 1;
	$count = 0;
	opendir(DIR, "$target/child");
	foreach(sort readdir(DIR)){
	    if ($_ !~ /^\./){
		$count++;
	    }
	}
	return if $count == 0;
    }
}

sub getTimestamp{
    my $time = shift;
    $time = time() if $time eq "";
    ($sec, $min, $hour, $mday, $mon, $year, $wday) = localtime($time);
    my $timestamp = sprintf ("%04d/%02d/%02d %02d:%02d:%02d", $year + 1900, $mon + 1, $mday, $hour, $min, $sec);
    return $timestamp;
}

sub log{
    my $comment = shift;
    print STDERR &getTimestamp() . " $comment\n";
    print LOG    &getTimestamp() . " $comment\n";
}

sub bynumber{
    $a <=> $b;
}
