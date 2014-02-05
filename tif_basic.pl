#!/usr/bin/perl
#
#  Copyright (c) 2012 - 2014
#       National Institute of Agrobiological Sciences.  All rights reserved.
#
# This code is derived from software contributed to NIAS by
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

$s = {};

# Tos17
#$head = "TGTTAAATATATATACAAGCT"; #21
#$tail = "TAGGTTGCAAGTTAGTTAAGA";
$head = "TGTTAAATATATATACA"; #17
$tail = "TTGCAAGTTAGTTAAGA";

$tsd_size = 5;
$file_list = "read/SRR556173_?.fastq"; # for ttm2
#$file_list = "read/SRR556174_?.fastq read/SRR556175_?.fastq"; # for ttm5

#P-element
#$head = "CATGATGAAATAACATAAGG"; #21
#$tail = "CCTTATGTTATTTCATCATG"; 
#$head = "CATGATGAAATAACAT"; #17
#$tail = "ATGTTATTTCATCATG";
#$tsd_size = 8;
#$file_list = "read/SRR823377_?.fastq"; # for D. melanogaster
#$file_list = "read/SRR823382_?.fastq"; # for D. melanogaster

$tail_size = length($tail);

open(IN, "cat $file_list |grep $head|");
while(<IN>){
    $pos = index($_, $head);
    $upstream = substr($_, 0, $pos);
    if (length($upstream) > 20){
        $tsd = substr($upstream, length($upstream) - $tsd_size, $tsd_size);
        if (length($s->{head}{$tsd}) < length($upstream)){
            $s->{head}{$tsd} = $upstream;
        }
    }
}
close(IN);

open(IN, "cat $file_list |grep $tail |");
while(<IN>){
    chomp;
    $pos = index($_, $tail);
    $downstream = substr($_, $pos + $tail_size, 100);
    if (length($downstream) > 20){
        $tsd =  substr($downstream, 0, $tsd_size);
        if (length($s->{tail}{$tsd}) < length($downstream)){
            $s->{tail}{$tsd} = $downstream;
        }
    }
}
close(IN);

open(OUT, ">tif.fasta");
foreach $tsd (sort keys %{$s->{head}}){
    if ($s->{tail}{$tsd} ne ""){
	$hj = substr($s->{head}{$tsd}, length($s->{head}{$tsd}) - 20, 20) . "_head";
	$tj = substr($s->{tail}{$tsd}, 0, 20) . "_tail";

        print OUT ">$hj
$s->{head}{$tsd}
>$tj
$s->{tail}{$tsd}\n";
    }
}
close(OUT);
