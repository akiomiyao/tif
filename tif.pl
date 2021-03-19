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

    For example,
    perl $0 reference.fasta target head_sequence tail_sequence

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

    Compressed file with gz, bz2 and xz extentions will be analyzed
    without pre-decompression. 

    reference.fasta is the reference sequence with multi-fasta format.
    head_sequence is a short sequence at 5'-end of target transposon.
    tail_sequence is a short sequence at 3'-end of target transposon.
    Length of head_seq and short_seq shold be from 17 to 21 bp.

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

$uname = `uname`;
chomp($uname);
if ($uname eq "Linux"){
    open(IN, "/proc/cpuinfo");
    while(<IN>){
        $maxprocess ++  if /processor/;
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

$hsize = length($head);

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
    &waitAll;
    system("cat $target/tmp.* > $target/selected.$head.$tail && rm $target/tmp.*");
}

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

&log("Reading chromosome sequences.");

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
        s/^>//;
        $name = (split)[0];
    }else{
	y/acgtn/ACGTN/;
	$chr{$name} .= $_;
    }
}
close(IN);

&log("Mapping of junction of head.");

foreach $upstream (sort keys %head){
    $junction = substr($upstream, length($upstream) - 20, 20);
    $rjunction = complement($junction);
    foreach $name (sort keys %chr){
	while (1) {
	    $pos = index($chr{$name}, $junction, $pos + 1);
	    if ($pos > -1){
		$length = length($upstream);
		$ref = substr($chr{$name}, $pos - $length + 20, $length);
		if ($ref eq $upstream){
		    $tpos = $pos + 20;
		    $mhc{$name}{$tpos} ++;
		    if ($length > length($maphead{$name}{$tpos}{forward})){
			$maphead{$name}{$tpos}{forward} = $upstream;
		    }
		}
	    }
	    $rpos = index($chr{$name}, $rjunction, $rpos + 1);
	    if ($rpos > -1){
		$cupstream = complement($upstream);
		$length = length($cupstream);
		$ref = substr($chr{$name}, $rpos, $length);
		if ($ref eq $cupstream){
		    $tpos = $rpos + 1;
		    $mhc{$name}{$tpos} ++;
		    if ($length > length($maphead{$name}{$tpos}{reverse})){
			$maphead{$name}{$tpos}{reverse} = $upstream;
		    }
		}
	    }
	    last if $pos == -1 and $rpos == -1;
	}
    }
}

&log("Mapping of junction of tail.");

foreach $downstream (sort keys %tail){
    $junction = substr($downstream, 0, 20);
    $rjunction = complement($junction);
    foreach $name (sort keys %chr){
	while (1) {
	    $pos = index($chr{$name}, $junction, $pos + 1);
	    if ($pos > -1){
		$downstream = $downstream;
		$length = length($downstream);
		$ref = substr($chr{$name}, $pos, $length);
		if ($ref eq $downstream){
		    $tpos = $pos + 1;
		    $mtc{$name}{$tpos} ++;
		    if ($length > length($maptail{$name}{$tpos}{forward})){
			$maptail{$name}{$tpos}{forward} = $downstream;
		    }
		}
	    }
	    $rpos = index($chr{$name}, $rjunction, $rpos + 1);
	    if ($rpos > -1){
		$cdownstream = complement($downstream);
		$length = length($cdownstream);
		$ref = substr($chr{$name}, $rpos - $length + 20, $length);
		if ($ref eq $cdownstream){
		    $tpos = $rpos + 20;
		    $mtc{$name}{$tpos} ++;
		    if ($length > length($maptail{$name}{$tpos}{reverse})){
			$maptail{$name}{$tpos}{reverse} = $downstream;
		    }
		}
	    }
	    last if $pos == -1 and $rpos == -1;
	}
    }
}

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
    $read_length = length($seq);
    $pos = index($seq, $tail);
    $downstream = substr($seq, $pos + $hsize, $read_length);
    if (length($downstream) > 20){
        $junction = substr($downstream, 0, 20);
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
	sleep 1;
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

sub waitAll{
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
