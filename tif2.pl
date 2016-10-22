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
#

# Usage: perl tif2.pl ref.fasta head_seq tail_seq
#
# Example: perl tif2.pl IRGSP1.0.fasta TGTTAAATATATATACA TTGCAAGTTAGTTAAGA
#
# ref.fasta is the path of reference in multi-fasta format.
# head_seq is a short sequence at 5'-end of target transposon.
# tail_seq is a short sequence at 3'-end of target transposon.
# Length of head_seq and short_seq is from 17 to 21 bp.
#
# All Short-read sequences with fastq format in ./read directory are analyzed.
# 

$ref = $ARGV[0];
$head = $ARGV[1];
$tail = $ARGV[2];


#tos17
#$head = "TGTTAAATATATATACAAGCT"; # 21 base
#$tail = "TAGGTTGCAAGTTAGTTAAGA";

# default
$ref  = "ref.fasta" if $ref eq "";
$head = "TGTTAAATATATATACA" if $head eq ""; # 17 base
$tail = "TTGCAAGTTAGTTAAGA" if $tail eq "";


$hsize = length($head);

open(IN, $ref);
while(<IN>){
    chomp;
    if (/^>/){
	$name = join('', (split)[0, 1]);
	$name =~ s/^>//;
    }else{
	y/acgtn/ACGTN/;
	$chr{$name} .= $_;
    }
}
close(IN);

open(IN, "cat ./read/*.fastq |grep $head|");
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

open(IN, "cat ./read/*.fastq |grep $tail |");
while(<IN>){
    chomp;
    $read_length = length($_);
    $pos = index($_, $tail);
    $downstream = substr($_, $pos + $hsize, $read_length);
    if (length($downstream) > 20){
        $junction = substr($downstream, 0, 20);
        if (length($downstream) > length($tail{$junction})){
	    $tail{$junction} = $downstream;
        }
    }
}
close(IN);

foreach $junction (sort keys %head){
    $rjunction = complement($junction);
    foreach $name (sort keys %chr){
	while (1) {
	    $pos = index($chr{$name}, $junction, $pos + 1);
	    $rpos = index($chr{$name}, $rjunction, $rpos + 1);
	    if ($pos > -1){
		$tpos = $pos + 20;
		$maphead{$name}{$tpos}{forward} = $head{$junction};
	    }
	    if ($rpos > -1){
		$tpos = $rpos + 1;
		$maphead{$name}{$tpos}{reverse} = $head{$junction};
	    }
	    last if $pos == -1 and $rpos == -1;
	}
    }
}

foreach $junction (sort keys %tail){
    $rjunction = complement($junction);
    foreach $name (sort keys %chr){
	while (1) {
	    $pos = index($chr{$name}, $junction, $pos + 1);
	    $rpos = index($chr{$name}, $rjunction, $rpos + 1);
	    if ($pos > -1){
		$tpos = $pos + 1;
		$maptail{$name}{$tpos}{forward} = $tail{$junction};
	    }
	    if ($rpos > -1){
		$tpos = $rpos + 20;
		$maptail{$name}{$tpos}{reverse} = $tail{$junction};
	    }
	    last if $pos == -1 and $rpos == -1;
	}
    }
}

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
                    print "$chr\t$pos\t$i\t$tsd_size\t$tsd_head\t$tsd_tail\t$direction\t$upstream\t$downstream\n";
                }
            }
        }
    }
}

sub bynumber{
    $a <=> $b;
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