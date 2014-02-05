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

#tos17
#$head = "TGTTAAATATATATACAAGCT"; # 21 base
#$tail = "TAGGTTGCAAGTTAGTTAAGA";
$head = "TGTTAAATATATATACA"; # 17 base
$tail = "TTGCAAGTTAGTTAAGA";

$file_list = "./read/SRR556173_?.fastq"; #for ttm2;
#$file_list = "./read/SRR556174_?.fastq ./read/SRR556175_?.fastq"; #for ttm5;

$blast_command = "blastall -p blastn -d ./blast/IRGSP1.0 -m 8 -b 1"; # legacy blast
#$blast_command = "~/ncbi-blast-2.2.29+/bin/blastn -db ./blast/IRGSP1.0 -outfmt 6 -num_alignments 1"; # blast plus

#P-element
#$head = "CATGATGAAATAACATAAGG"; # 21 base
#$tail = "CCTTATGTTATTTCATCATG"; 
#$head = "CATGATGAAATAACATA"; # 17 base
#$tail = "TATGTTATTTCATCATG";
#$tsd_size = 8;
#$file_list = "read/SRR823377_?.fastq"; # for D. melanogaster
#$file_list = "read/SRR823382_?.fastq"; # for D. melanogaster

#$blast_command = "blastall -p blastn -d ./blast/dmel-r5.55 -m 8 -b 1"; # legacy blast
#$blast_command = "~/ncbi-blast-2.2.29+/bin/blastn -db ./blast/dmel-r5.55 -outfmt 6 -num_alignments 1"; # blast plus

$hsize = length($head);

$s = {};

open(IN, "cat $file_list |grep $head|");
while(<IN>){
    $pos = index($_, $head);
    $upstream = substr($_, 0, $pos);
    if (length($upstream) > 20){
        $junction = substr($upstream, length($upstream) - 20, 20);
        if (length($upstream) > length($head{$junction})){
            $head{$junction} = $upstream;
        }
    }
}
close(IN);

open(IN, "cat $file_list |grep $tail |");
while(<IN>){
    chomp;
    $pos = index($_, $tail);
    $downstream = substr($_, $pos + $hsize, 100);
    if (length($downstream) > 20){
        $junction = substr($downstream, 0, 20);
        if (length($downstream) > length($tail{$junction})){
            $tail{$junction} = $downstream;
        }
    }
}
close(IN);

open(OUT, ">tif.fasta");
foreach $junction (sort keys %head){
    $hj = $junction . "_head";
    print OUT ">$hj
$head{$junction}
";
}

foreach $junction (sort keys %tail){
    $tj = $junction . "_tail";
    print OUT ">$tj
$tail{$junction}
";
}
close(OUT);

open(IN, "cat tif.fasta | $blast_command |");
while(<IN>){
    @row = split;
    if ($row[2] > 96){
        if (/head/){
            $s->{head}{$row[1]}{$row[9]} = $row[0];
        }elsif (/tail/){
            $s->{tail}{$row[1]}{$row[8]} = $row[0];
        }
    }
}

open(OUT, "> tif.position");
foreach $chr (sort keys %{$s->{head}}){
    foreach $head (sort bynumber keys %{$s->{head}{$chr}}){
        foreach $tail (sort bynumber keys %{$s->{tail}{$chr}}){
            $tsd_size = abs($tail - $head) +1;
            if ($tsd_size < 12 and $tsd_size > 1){
                $seq = (split('_', $s->{head}{$chr}{$head}))[0];
                $tsd_head = substr($seq, length($seq) - $tsd_size, $tsd_size);
                $seq = (split('_', $s->{tail}{$chr}{$tail}))[0];
                $tsd_tail = substr($seq, 0, $tsd_size);
                if ($tail < $head){
                    $direction = "forward";
                }else{
                    $direction = "reverse";
                }
                print OUT "$chr\t$tail\t$head\t$tsd_size\t$tsd_tail\t$direction\n" if $tsd_tail eq $tsd_head;
            }
        }
    }
}
close(OUT);

sub bynumber{
    $a <=> $b;
}
