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
# Version 1.1
#

$file_list = "./read/*.fastq"; 

#tos17
$head = "TGTTAAATATATATACA"; # 17 base
$tail = "TTGCAAGTTAGTTAAGA";

$blast_command = "blastn -db IRGSP-1.0_genome.fasta -outfmt 6 -num_alignments 1"; # blast plus
#$blast_command = "blastall -p blastn -d ./blast/IRGSP1.0 -m 8 -b 1"; # legacy blast

#P-element
#$head = "CATGATGAAATAACATA"; # 17 base
#$tail = "TATGTTATTTCATCATG";

#$blast_command = "blastall -p blastn -d ./dmel-r5.55 -m 8 -b 1"; # legacy blast
#$blast_command = "~/blastn -db dmel-r5.55.fasta -outfmt 6 -num_alignments 1"; # blast plus

$start = time();
($sec, $min, $hour, $mday, $mon, $year, $wday) = localtime($start);
printf STDERR ("TIF Start: %04d/%02d/%02d %02d:%02d:%02d\n", $year + 1900, $mon + 1, $mday, $hour, $min, $sec);

$tsize = length($tail);

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
    $read_length = length($_);
    $pos = index($_, $tail);
    $downstream = substr($_, $pos + $tsize, $read_length);
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
            if ($tsd_size < 10 and $tsd_size > 1){
                $seq = (split('_', $s->{head}{$chr}{$head}))[0];
                $tsd_head = substr($seq, length($seq) - $tsd_size, $tsd_size);
                $seq = (split('_', $s->{tail}{$chr}{$tail}))[0];
                $tsd_tail = substr($seq, 0, $tsd_size);
                if ($tail < $head){
                    $direction = "forward";
                }else{
                    $direction = "reverse";
                }
                print     "$chr\t$tail\t$head\t$tsd_size\t$tsd_tail\t$direction\n" if $tsd_tail eq $tsd_head;
                print OUT "$chr\t$tail\t$head\t$tsd_size\t$tsd_tail\t$direction\n" if $tsd_tail eq $tsd_head;
            }
        }
    }
}
close(OUT);

($sec, $min, $hour, $mday, $mon, $year, $wday) = localtime($start);
printf STDERR ("TIF Start: %04d/%02d/%02d %02d:%02d:%02d\n", $year + 1900, $mon + 1, $mday, $hour, $min, $sec);
$end = time();
($sec, $min, $hour, $mday, $mon, $year, $wday) = localtime($end);
printf STDERR ("TIF End: %04d/%02d/%02d %02d:%02d:%02d\n", $year + 1900, $mon + 1, $mday, $hour, $min, $sec);

$elapsed = $end - $start;

print STDERR "$elapsed seconds elapsed.\n";

sub bynumber{
    $a <=> $b;
}
