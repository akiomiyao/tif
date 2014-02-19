#!/usr/bin/perl
#
# TIF Program (Original version)
# Script name: tif_original.pl
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

$head = "TGTTAAATATATATACA";  # 5'-end sequence of  Tos17
$tail = "TTGCAAGTTAGTTAAGA"; # 3'-end sequence of Tos17
$tsd_size = 5; # size of target site duplication

$tail_size = length($tail);

$file_list = "/directory_of_fastq_file/*.fastq";

open(IN, "cat $file_list |grep $head|");
while(<IN>){
    $pos = index($_, $head);
    $upstream = substr($_, 0, $pos);
    if (length($upstream) > 20){
        $tsd = substr($upstream, length($upstream) - $tsd_size, $tsd_size);
        if (length($head{$tsd}) < length($upstream)){
            $head{$tsd} = $upstream;
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
        if (length($tail{$tsd}) < length($downstream)){
            $tail{$tsd} = $downstream;
        }
    }
}
close(IN);

foreach $tsd (sort keys %head){
    if ($tail{$tsd} ne ""){
        print ">TSD $tsd head
$head{$tsd}
>TSD $tsd tail
$tail{$tsd}\n";
    }
}
