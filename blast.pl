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

#
# This script is used for result of tif_basic.pl
#

$s = {};

#for rice
#$blast_command = "blastall -p blastn -d ./blast/IRGSP1.0 -m 8 -b 1"; # legacy blast
$blast_command = "ncbi-blast-2.9.0+/bin/blastn -db blast/IRGSP-1.0_genome.fasta -outfmt 6 -num_alignments 1"; # blast plus

#for Drosophila melanogaster
#$blast_command = "blastall -p blastn -d ./blast/dmel-r5.55 -m 8 -b 1"; # legacy blast
#$blast_command = "~/ncbi-blast-2.2.29+/bin/blastn -db ./blast/dmel-r5.55 -outfmt 6 -num_alignments 1"; # blast plus

open(IN, "cat tif.fasta | $blast_command |"); # tif.fasta is a result of tif_basic.pl
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
                print     "$chr\t$tail\t$head\t$tsd_size\t$tsd_tail\t$direction\n" if $tsd_tail eq $tsd_head;
                print OUT "$chr\t$tail\t$head\t$tsd_size\t$tsd_tail\t$direction\n" if $tsd_tail eq $tsd_head;
            }
        }
    }
}
close(OUT);

sub bynumber{
    $a <=> $b;
}
