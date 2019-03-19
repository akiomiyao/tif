# TIF

### Transposon Insertion Finder

Transposon Insertion Finder (TIF) is a search program to detect insertions of transposable element from short reads of next generation sequencer. The program is written in a perl script. The program runs on Unix (Linux) platform. The TIF requires short sequences of both ends of target transposable element and length of target site duplication (TSD). Basic TIF (tif_basic.pl) does not require reference genome sequence to select short reads containing target sites, whereas extended TIF (tif_extended require the reference sequence and BLAST program.

TIF is one of the fastest and the smallest program among analysis programs of next generation sequencing (NGS). The distinctive feature of TIF is direct selection containing end sequences of the target transposable element from short reads of NGS.

### Update
Now, new script tif.pl is implemented. (2019-03-19)

     e.g. perl tif.pl ref.fasta TGTTAAATATATATACA TTGCAAGTTAGTTAAGA
     First argument is the path of reference sequence with multi-fasta format.
     Second argument is the head sequence of transposon.
     Third argument is the tail sequence of transposon.
     All short reads (e.g name_r1.fastq, name_r2.fastq) in './read' directory will be analyzed.
     Run without argument, help will be shown.

This version does not depend on BLAST search. Search script was included in tif.pl.

### Download TIF
Download zip file of PED from https://github.com/akiomiyao/tif and extract.  
or  
% git clone https://github.com/akiomiyao/tif.git  
If you got scripts from github, update to newest version is very easy using pull command of git.  
% git pull  

### Static data required by TIF (for demonstration)

You can exchange any head and tail sequences of transposon which you want to detect transpositions.

For Tos17 retrotransposon of rice

      Head of Tos17: TGTTAAATATATATACA
      Tail of Tos17: TTGCAAGTTAGTTAAGA
      Size of TSD: 5

For mPing transposon of rice (DNA type transposon)

      Head of mPing: GGCCAGTCACAATGGGG
      Tail of mPing: AGCCATTGTGACTGGCC
      Size of TSD: 3
      Because TSD of mPing is too short, extended TIF is recommended. 

For nDart transposon of rice (DNA type transposon)

      Head of nDart: TAGAGGTGGCCAAACGGGC
      Tail of nDart: GCCCGTTTGGCCACCTCTA
      Size of TSD: 8


For P-element of Drosophila melanogaster

      Head of P-element: CATGATGAAATAACAT
      Tail of P-element: ATGTTATTTCATCATG
      Size of TSD: 8


### To obtain short read data

Download sra tool kit from

      https://trace.ncbi.nlm.nih.gov/Traces/sra/sra.cgi?view=software

In your home directory,

      tar xvfz sratoolkit.2.9.6-centos_linux64.tar.gz
      copy fastq-dump in bin directory to executable directory.

For ttm2 (Rice mutant)

        cd tif/read
        fastq-dump --split-files -A SRR556173
    
For ttm5 (Rice mutant)

        cd tif/read
        fastq-dump --split-files -A SRR556174
        fastq-dump --split-files -A SRR556175
    
For D. melanogaster

        cd tif/read
        fastq-dump --split-files -A SRR823377
        fastq-dump --split-files -A SRR823382



### BLAST programs

Download BLAST programs 

      ftp://ftp.ncbi.nlm.nih.gov/blast/executables/release/LATEST/blast-2.2.26-x64-linux.tar.gz
      ftp://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/ncbi-blast-2.2.29+-x64-linux.tar.gz

New version of BLAST can be downloaded from 
      https://ftp.ncbi.nlm.nih.gov/blast/executables/

To make blast data base

      mkdir tif/blast

For Rice

      cd tif/blast
      wget http://rapdb.dna.affrc.go.jp/download/archive/irgsp1/IRGSP-1.0_genome.fasta.gz
      gzip -d IRGSP-1.0_genome.fasta.gz
      formatdb -t IRGSP1.0 -p F -n IRGSP1.0 -i IRGSP-1.0_genome.fasta

For Drosophira melanogaster

      cd tif/blast
      wget ftp://ftp.flybase.net/genomes/Drosophila_melanogaster/dmel_r5.55_FB2014_01/fasta/dmel-all-chromosome-r5.55.fasta.gz
      gzip -d dmel-all-chromosome-r5.55.fasta.gz
      formatdb -t dmel-r5.55 -p F -n dmel-r5.55 -i dmel-all-chromosome-r5.55.fasta


### Search targets of transposon

Edit and select head, tail, tsd_size and file list in tif_basic.pl and blast.pl, or tif_extended.pl.

To test TIF algorithm 1

      cd tif
      perl tif_basic.pl
      perl blast.pl
      
Output of tif_basic.pl is tif.fasta, a multiple FASTA file.
  
The blast.pl reads tif.fasta and returns tif.position containing location and direction of TE insertion sites.

To test TIF algorithm 2

      cd tif
      perl tif_extended.pl
      
The tif_extended.pl returns both tif.fasta and tif.position files.

For new extended tif

      cd tif
      perl tif2.pl

The tif2.pl reads nucleotide sequence of rice genome saved in chr directory, position of junction will be detected by text search against the genome sequence. This enable to detect insertions even on repetitive loci.

### Citing TIF

- Cite this article as: Nakagome M, Solovieva E, Takahashi A, Yasue H, Hirochika H, Miyao A (2014) Transposon Insertion Finder (TIF): a novel program for detection of de novo transpositions of transposable elements. BMC Bioinformatics 5:71.
- doi:10.1186/1471-2105-15-71

### Update

- 1.5 tif.pl is implemented. 2019-03-19
- 1.4 tif2.pl is improved. 2016-10-22
- 1.3 Add new extended version tif2.pl 2015-03-02
- 1.2 Update README.md 2014-10-09
- 1.1 Update link of SRA-toolkit in README.md 2014-08-01
- 1.0 Inital version 2014-02-05
