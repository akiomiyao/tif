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

$usage = "    $0 - Transposon Insertion Finder (Algorithm 2)

    perl $0 reference.fasta target head_sequence tail_sequence

    Another format
    perl $0 ref=reference_fasta,target=target,head=head_sequence,tail=tail_sequence,max_process=8,nogenotype
    max_process and nogenotype are optional.

    For endogenous retrotransposon Tos17 in rice,
    perl $0 IRGSP-1.0_genome.fasta.gz ttm5 TGTTAAATATATATACA TTGCAAGTTAGTTAAGA
    perl $0 ref=IRGSP-1.0_genome.fasta.gz,target=ttm5,head=TGTTAAATATATATACA,tail=TTGCAAGTTAGTTAAGA
    perl $0 ref=IRGSP-1.0_genome.fasta.gz,target=ttm5,head=TGTTAAATATATATACA,tail=TTGCAAGTTAGTTAAGA,max_process=8,nogenotype

    For P-element of Drosophila malanogaster
    perl $0 dmel-all-chromosome-r6.26.fasta target CATGATGAAATAACAT ATGTTATTTCATCATG

    Example
    In the directory of $0 script,
    mkdir target
    mkdir target/read
    cp somewhere/fastq_files target/read
    cp somewhere/reference_fasta .
    perl $0 reference_fasta target head_sequence tail_sequence

    Result and log will be saved in target directory.

    Compressed read and genome files with gz, bz2 and xz extentions
    will be analyzed without pre-decompression. 

    reference.fasta is the reference sequence with multi-fasta format.
    head_sequence is a short sequence at 5'-end of target transposon.
    tail_sequence is a short sequence at 3'-end of target transposon.
    Length of head_seq and short_seq shold be from 17 to 21 bp.
    max_process_number is optional. Default number of maximun process
    is number of CPUs. Calculation of genotype requires long time.
    If you do not want to genotype data, add 'nogenotype' in the option.

    Result is tif_result.head_sequence.tail_sequence file in the target directory.

    Result format: chromosome, junction of head_seq, junction of tail_seq,
    TSD size, TSD, direction of insertion, flanking sequence of head,
    flanking sequence of tail, number of flanking sequence of head,
    number of flanking sequence of tail, number of wild type reads, genotype.
    Genotypes M and H are homozygous and heterozygous mutation, respectively.

    Reference:
    Nakagome, M., Solovieva, E., Takahashi, A., Yasue, H., Hirochika. H., Miyao, A.
    Transposon Insertion Finder (TIF): a novel program for detection of de novo
    transpositions of transposable elements.
    BMC Bioinformatics 15, 71 (2014).
    https://doi.org/10.1186/1471-2105-15-71
    GitHub: https://github.com/akiomiyao/tif

    Author: Akio Miyao <miyao\@affrc.go.jp>

";

$ref        = $ARGV[0];
$target     = $ARGV[1];
$head       = $ARGV[2];
$tail       = $ARGV[3];

foreach (@ARGV){
    if (/,/){
	foreach (split(',', $_)){
	    next if $_ eq "";
	    if (/=/){
		($name, $value) = split('=', $_);
		$$name = $value;
	    }else{
		$$_ = 1;
	    }
	}
    }else{
	if (/=/){
	    ($name, $value) = split('=', $_);
	    $$name = $value;
	}else{
	    $$_ = 1;
	}
    }
}

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
}else{
    if (-d "$target/child"){
	system("rm -r $target/child");
    }
    system("mkdir $target/child");
     if (! -d "$target/tmp"){
	system("mkdir $target/tmp");
    }
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

$maxprocess = $max_process if $max_process ne "";

$start = time();
open(LOG, "> $target/tif_log.$head.$tail");
print LOG "Program: $0
Reference: $ref
Target: $target
Head: $head
Tail: $tail
Max process: $maxprocess

";
&log("TIF start.");

$rhead = complement($head);
$rtail = complement($tail);

if (! -e "$target/tif_grep.$head.$tail"){
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
	if ($catcmd eq "cat"){
	    system("touch $target/child/$head.$file.h && grep $head $target/read/$file > $target/tmp/tmp.$file.$head && rm $target/child/$head.$file.h || rm $target/child/$head.$file.h &");
	}else{
	    system("touch $target/child/$head.$file.h && $catcmd $target/read/$file | grep $head > $target/tmp/tmp.$file.$head && rm $target/child/$head.$file.h || rm $target/child/$head.$file.h &");
	}
	&waitFork;
	&log("searching $tail in $file");
	if ($catcmd eq "cat"){
	    system("touch $target/child/$tail.$file.t && grep $tail $target/read/$file > $target/tmp/tmp.$file.$tail && rm $target/child/$tail.$file.t || rm $target/child/$tail.$file.t&");
	}else{
	    system("touch $target/child/$tail.$file.t && $catcmd $target/read/$file | grep $tail > $target/tmp/tmp.$file.$tail && rm $target/child/$tail.$file.t || rm $target/child/$tail.$file.t &");
	}
	&waitFork;
	&log("searching $rhead in $file");
	if ($catcmd eq "cat"){
	    system("touch $target/child/$rhead.$file.rh && grep $rhead $target/read/$file > $target/tmp/tmp.$file.$rhead && rm $target/child/$rhead.$file.rh || rm $target/child/$rhead.$file.rh &");
	}else{
	    system("touch $target/child/$rhead.$file.rh && $catcmd $target/read/$file | grep $rhead > $target/tmp/tmp.$file.$rhead && rm $target/child/$rhead.$file.rh || rm $target/child/$rhead.$file.rh &");
	}
	&waitFork;
	&log("searching $rtail in $file");
	if ($catcmd eq "cat"){
	    system("touch $target/child/$rtail.$file.rt && grep $rtail $target/read/$file > $target/tmp/tmp.$file.$rtail && rm $target/child/$rtail.$file.rt ||rm  $target/child/$rtail.$file.rt &");
	}else{
	    system("touch $target/child/$rtail.$file.rt && $catcmd $target/read/$file | grep $rtail > $target/tmp/tmp.$file.$rtail && rm $target/child/$rtail.$file.rt || rm $target/child/$rtail.$file.rt &");
	}
    }
    &joinAfterGrep;
    system("cat $target/tmp/tmp.* > $target/tif_grep.$head.$tail && rm $target/tmp/tmp.*");
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
    &log("Mapping of junction in $name.");
    system("perl $0 map:$name $target $head $tail &");
}

&join;

open(IN, "cat $target/tmp/map.$head.$tail.* |");
opendir (DIR, "chr");
foreach $chr (sort readdir(DIR)){
    if (-s "$target/tmp/map.$head.$tail.$chr" > 0){
	open(IN, "chr/$chr");
	$seq = <IN>;
	close(IN);
        $pos = (split("-", $chr))[1];
	open(IN, "$target/tmp/map.$head.$tail.$chr");
	while(<IN>){
	    chomp;
	    $wt = "";
	    @row = split('\t', $_);
	    if ($row[0] eq "head"){
		if (length($row[4]) > length($maphead{$row[1]}{$row[2]}{$row[3]})){
		    $maphead{$row[1]}{$row[2]}{$row[3]} = $row[4];
		    $mhc{$row[1]}{$row[2]} ++;
		}
		if ($wt{"$row[1] $row[2] f"} eq ""){
		    $wt = substr($seq, $row[2] - $pos + 500 - 21, 40);
		    if (length($wt) == 40 or $wt !~ /NN|\ /){
			$cwt = complement($wt);
			$wt{"$row[1] $row[2] f"} = $wt;
			$wt{"$row[1] $row[2] r"} = $cwt;
			$wtcount{$wt} = 0;
			$wtcount{$cwt} = 0;
		    }
		}
	    }else{
		if (length($row[4]) > length($maptail{$row[1]}{$row[2]}{$row[3]})){
		    $maptail{$row[1]}{$row[2]}{$row[3]} = $row[4];
		    $mtc{$row[1]}{$row[2]} ++;
		}
		if ($wt{"$row[1] $row[2] f"} eq ""){
		    $wt = substr($seq, $row[2] -$pos + 500 - 21, 40);
		    if (length($wt) == 40 or  $wt !~ /NN|\ /){
			$cwt = complement($wt);
			$wt{"$row[1] $row[2] f"} = $wt;
			$wt{"$row[1] $row[2] r"} = $cwt;
			$wtcount{$wt} = 0;
			$wtcount{$cwt} = 0;
		    }
		}
	    }
	}
	close(IN);
    }
}

if (! $nogenotype){
    &log("Counting wild type reads.");
    opendir(DIR, "$target/read");
    foreach $file (sort readdir(DIR)){
	next if $file =~ /^\./;
	&log("Processing $file");
	if ($file =~ /.gz$/){
	    open(IN, "zcat $target/read/$file |");
	    $fastq = 1;
	}elsif($file =~ /.bz2$/){
	    open(IN, "bzcat $target/read/$file |");
	    $fastq = 1;
	}elsif($file =~ /.xz$/){
	    open(IN, "xzcat $target/read/$file |");
	    $fastq = 1;
	}elsif($file =~ /.fastq$|.fq$/){
	    open(IN, "$target/read/$file");
	    $fastq = 1;
	}else{
	    open(IN, "$target/read/$file");
	}
	while(<IN>){
	    chomp;
	    if ($fastq){
		$count++;
		if ($count == 2){
		    $complement = &complement($_);
		    for($i = 0; $i <= length($_) - 40; $i++){
			$fragment = substr($_, $i, 40);
			if (defined $wtcount{$fragment}){
			    $wtcount{$fragment}++;
			}
			if ($file !~/sort_uniq/){
			    $fragment = substr($complement, $i, 40);
			    if (defined $wtcount{$fragment}){
				$wtcount{$fragment}++;
			    }
			}
		    }
		    $reads ++;
		    if ($reads % 1000000 == 0){
			&log("Counting wild type. $reads reads processed.");
		    }
		}elsif($count == 4){
		    $count = 0;
		}
	    }else{
		for($i = 0; $i <= length($_) - 40; $i++){
		    $fragment = substr($_, $i, 40);
		    if (defined $wtcount{$fragment}){
			$wtcount{$fragment}++;
		    }
		    if ($file !~/sort_uniq/){
			$fragment = substr($complement, $i, 40);
			if (defined $wtcount{$fragment}){
			    $wtcount{$fragment}++;
			}
		    }
		}
		$reads ++;
		if ($reads % 1000000 == 0){
		    &log("Counting wild type. $reads reads processed.");
		}
	    }
	}
    }
}

open(OUT, "> $target/tif_result.$head.$tail");
foreach $chr (sort keys %maphead){
    foreach $pos (sort bynumber keys %{$maphead{$chr}}){
        foreach $direction (sort keys %{$maphead{$chr}{$pos}}){
            $upstream = $maphead{$chr}{$pos}{$direction};
            for($i = $pos - 20 ; $i <= $pos + 20; $i++){
                if ($maptail{$chr}{$i}{$direction} ne ""){
                    $downstream = $maptail{$chr}{$i}{$direction};
                    $tsd_size = abs($pos - $i) + 1;
                    $tsd_head =  substr($upstream, length($upstream) - $tsd_size, $tsd_size);
                    $tsd_tail =  substr($downstream, 0, $tsd_size);
		    if ($tsd_head eq $tsd_tail){
			if ($nogenotype){
			    my $output = "$chr\t$pos\t$i\t$tsd_size\t$tsd_head\t$direction\t$upstream\t$downstream\t$mhc{$chr}{$pos}\t$mtc{$chr}{$i}\n";
			    print STDERR $output;
			    print OUT $output;
			}else{
			    $wt = 0;
			    $fragment = $wt{"$chr $pos f"};
			    $wt += $wtcount{$fragment} if $fragment ne "";
			    $fragment = $wt{"$chr $pos r"};
			    $wt += $wtcount{$fragment} if $fragment ne "";
			    if ($wt > 0){
				$genotype = "H";
			    }else{
				$genotype = "M";
			    }
			    my $output = "$chr\t$pos\t$i\t$tsd_size\t$tsd_head\t$direction\t$upstream\t$downstream\t$mhc{$chr}{$pos}\t$mtc{$chr}{$i}\t$wt\t$genotype\n";
			    print STDERR $output;
			    print OUT $output;
			}
		    }
		}
            }
        }
    }
}
close(OUT);

$timestamp = `date '+%Y-%m-%d %H:%M:%S %z'`;
chomp($timestamp);
$filedate = (split('\ ', $timestamp))[0];
$filedate =~ y/-//d;
open(IN,"$target/tif_result.$head.$tail");
open(OUT,"> $target/tif_result.$head.$tail.vcf");
if ($nogenotype){
    print OUT "##fileformat=VCFv4.3
##fileDate=$filedate
##source=<PROGRAM=tif.pl,target=$target,reference=$ref>
##INFO=<ID=SVTYPE,Number=1,Type=String,Description=\"Type of structural variant\">
##INFO=<ID=MEINFO,Number=9,Type=String,Description=\"Movile element info of the form ME_HEAD_SEQ,ME_TAIL_SEQ,JUNCTION_POS_OF_HEAD,JUNCTION_POS_OF_TAIL,TSD_SIZE,TSD_SEQUENCE,DIRECTION,COUNT_OF_READS_WITH_JUNCTION_OF_HEAD,COUNT_OF_READS_WITH_JUNCTION_OF_TAIL\">
##ALT=<ID=INS,Description=\"Insertion of a mobile element\">
##created=<TIMESTAMP=\"$timestamp\">
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n";
}else{
    print OUT "##fileformat=VCFv4.3
##fileDate=$filedate
##source=<PROGRAM=tif.pl,target=$target,reference=$ref>
##INFO=<ID=GT,Number=1,Type=String,Description=\"Genotype\">
##INFO=<ID=AF,Number=.,Type=Float,Description=\"Allele Frequency\">
##INFO=<ID=DP,Number=1,Type=Integer,Description=\"Approximate read depth\">
##INFO=<ID=SVTYPE,Number=1,Type=String,Description=\"Type of structural variant\">
##INFO=<ID=MEINFO,Number=9,Type=String,Description=\"Movile element info of the form ME_HEAD_SEQ,ME_TAIL_SEQ,JUNCTION_POS_OF_HEAD,JUNCTION_POS_OF_TAIL,TSD_SIZE,TSD_SEQUENCE,DIRECTION,COUNT_OF_READS_WITH_JUNCTION_OF_HEAD,COUNT_OF_READS_WITH_JUNCTION_OF_TAIL\">
##ALT=<ID=INS,Description=\"Insertion of a mobile element\">
##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">
##FORMAT=<ID=AD,Number=.,Type=String,Description=\"Allelic depths for the reference and alternate alleles in the order listed\">
##FORMAT=<ID=DP,Number=1,Type=String,Description=\"Read Depth\">
##created=<TIMESTAMP=\"$timestamp\">
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t$target\n";
}
while(<IN>){
    chomp;
    @row = split;
    $rnuc = substr($row[6], length($row[6]) - 1, 1);
    my $chr;
    ($chr = $row[0]) =~ s/chr//;
    $chr += 0 if $chr =~/^[0-9]+$/;
    if ($nogenotype){
	print OUT "$chr\t$row[1]\t.\t$rnuc\t$rnuc<INS>\t.\t.\tMEINFO=$head,$tail,$row[1],$row[2],$row[3],$row[4],$row[5],$row[8],$row[9];SVTYPE=INS\n";
    }else{
	if ($row[10] > 0){
	    $gt = "1/0";
	}else{
	    $gt = "1/1";
	}
	my $total = $row[8] + $row[9] + $row[10];
	my $alt = $row[8] + $row[9];
	my $wt = $row[10];
	my $ad = "$wt,$alt";
	my $af = int($alt * 1000 / $total) / 1000;
	print OUT "$chr\t$row[1]\t.\t$rnuc\t$rnuc<INS>\t.\t.\tGT=$gt;AF=$af;DP=$total;MEINFO=$head,$tail,$row[1],$row[2],$row[3],$row[4],$row[5],$row[8],$row[9];SVTYPE=INS\tGT:AD:DP\t$gt:$ad:$total\n";
    }
}
close(IN);
close(OUT);

$end = time();

&log("TIF end. $target");

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

sub map{
    open(IN, "$target/tif_grep.$head.$tail");
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
    
    open(OUT, "> $target/tmp/map.$head.$tail.$name");
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

sub joinAfterGrep{
    my $count;
    while(1){
	sleep 1;
	$count = 0;
	open(IN, "ps xaww |");
	while(<IN>){
	    $count ++ if /grep/;
	}
	close(IN);
	if ($count == 0){
	    system("rm $target/child/* > /dev/null 2>&1");
	    return;
	}
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
    $comment = &getTimestamp() . " $comment\n";
    print STDERR $comment;
    print LOG    $comment;
}

sub bynumber{
    $a <=> $b;
}
