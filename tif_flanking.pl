#!/usr/bin/perl
#
#  Copyright (c) 2019
#       National Agriculture and Food Research Organization.  All rights reserved.
#
#       
# This code is derived from software contributed to NARO by
# Akio Miyao.
#
# Redistribution and use in source and binary forms, with or without
# modification, are permitted provided that the following conditions
# are met:
# 1. Redistributions of source code must retain the above copyright
#    notice, this list of conditions and the following disclaimer.
# 2. Redistributions in binary form must reproduce the above copyright
#    notice, this list of conditions and the following disclaimer in the
#    documentation and/or other materials provided with the distribution.
# 3. Neither the name of the Institute nor the names of its contributors
#    may be used to endorse or promote products derived from this software
#    without specific prior written permission.
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
# Version 1.0
#

$usage = "    tif_flanking.pl - Select flanking sequences of transposon insertion
    
    This sciript has same function of tif_basic.pl based on algorism 1.

    For example,
    perl tif_flanking.pl head_sequence tail_sequence TSD_size
    
    For Tos17,
    perl tif_flanking.pl TGTTAAATATATATACA TTGCAAGTTAGTTAAGA 5

    For P-element,
    perl tif_flanking.pl CATGATGAAATAACATA TATGTTATTTCATCATG 8

    Short reads in read directory will be analyzed.
    Compressed fastq file with gz, bz2 and xz extentions will be analyzed
    without pre-decompression. All files in read directory should be
    same format.

Author: Akio Miyao <miyao\@affrc.go.jp>

";

$s = {};

$head     = $ARGV[0];
$tail     = $ARGV[1];
$tsd_size = $ARGV[2];

if ($tsd_size eq ""){
    print $usage;
    exit;
}

$tail_size = length($tail);

$start = time();
($sec, $min, $hour, $mday, $mon, $year, $wday) = localtime($start);
printf STDERR ("TIF Start: %04d/%02d/%02d %02d:%02d:%02d\n", $year + 1900, $mon + 1, $mday, $hour, $min, $sec);

opendir(DIR, "read");
foreach (readdir(DIR)){
    if (/gz$/){
	$command = "zcat";
    }elsif(/bz2$/){
	$command = "bzcat";
    }elsif(/xz$/){
	$command = "xzcat";
    }
}

$command = "cat" if $command eq "";

open(IN, "$command read/* |");
while(<IN>){
    $line = 0 if $line ++ == 3;
    $total++;
    if ($total % 1000000 == 0){
	print STDERR "$total reads analyzed.\n";
    }
    if ($line == 2){
	chomp;
	$comp = &complement($_);
	if (/$head/){
	    &addHead($_);
	    next;
	}elsif(/$tail/){
	    &addTail($_);
	    next;
	}
	if ($comp =~ /$head/){
	    &addHead($comp);
	    next;
	}elsif($comp =~ /$tail/){
	    &addTail($comp);
	    next;
	}
    }
}
close(IN);

open(OUT, ">tif.fasta");
foreach $tsd (sort keys %{$s->{head}}){
    if ($s->{tail}{$tsd} ne ""){
	$hj = substr($s->{head}{$tsd}, length($s->{head}{$tsd}) - 20, 20) . "_head";
	$tj = substr($s->{tail}{$tsd}, 0, 20) . "_tail";

        $result = ">$hj
$s->{head}{$tsd}
>$tj
$s->{tail}{$tsd}\n";
    }
    print $result;
    print OUT $result;
}
close(OUT);

print STDERR "\n";
($sec, $min, $hour, $mday, $mon, $year, $wday) = localtime($start);
printf STDERR ("TIF Start: %04d/%02d/%02d %02d:%02d:%02d\n", $year + 1900, $mon + 1, $mday, $hour, $min, $sec);
$end = time();
($sec, $min, $hour, $mday, $mon, $year, $wday) = localtime($end);
printf STDERR ("TIF End: %04d/%02d/%02d %02d:%02d:%02d\n", $year + 1900, $mon + 1, $mday, $hour, $min, $sec);

$elapsed = $end - $start;

print STDERR "$elapsed seconds elapsed.

Result was saved to tif.fasta\n";

sub addHead{
    my $seq = shift;
    my $pos = index($seq, $head);
    my $upstream = substr($seq, 0, $pos);
    my $tsd;
    if (length($upstream) > 20){
        $tsd = substr($upstream, length($upstream) - $tsd_size, $tsd_size);
        if (length($s->{head}{$tsd}) < length($upstream)){
            $s->{head}{$tsd} = $upstream;
        }
    }
}

sub addTail{
    my $seq = shift;
    my $read_length = length($seq);
    my $pos = index($seq, $tail);
    my $downstream = substr($seq, $pos + $tail_size, $read_length);
    my $tsd;
    if (length($downstream) > 20){
        $tsd =  substr($downstream, 0, $tsd_size);
        if (length($s->{tail}{$tsd}) < length($downstream)){
            $s->{tail}{$tsd} = $downstream;
        }
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
