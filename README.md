# TIF

### Transposon Insertion Finder

Transposon Insertion Finder (TIF) is a search program to detect insertions of transposable element from short reads of next generation sequencer. The program is written in a perl script. The program runs on Unix (Linux) platform. The TIF requires short sequences of both ends of target transposable element. Basic TIF (tif_basic.pl) does not require reference genome sequence to select short reads containing target sites, whereas extended TIF (tif_extended require the reference sequence and BLAST program). Length of target site duplication (TSD) is required by tif_basic.pl and tif_original.pl.  

TIF (tif_original.pl) is one of the fastest and the smallest program among analysis programs of next generation sequencing (NGS). The distinctive feature of TIF is direct selection containing end sequences of the target transposable element from short reads of NGS.

Web page: https://akiomiyao.github.io/tif/

### Update
- tif.pl genotyping function has been added. (2022-01-04)  
     Run tif.pl without an argument, help will be shown.

- tif.pl has been improved. (2021-03-23)  

     *e.g.* perl tif.pl reference_genome.fasta target_directory TE_head_sequence TE_tail_sequence  
     First argument is the path of reference sequence with multi-fasta format.  
     Second argument is the target directory containing read directory.  
     Third argument is the head sequence of transposon.  
     Fourth argument is the tail sequence of transposon.  
     Fifth argument is maximun number of process (optional).  
     All short reads (*e.g* name_r1.fastq, name_r2.fastq) in 'target/read' directory will be analyzed.  
     TE_head and TE_tail sequences are 17-23 bases of 5'-end, and 3'-end of the transposable element.  
     Run without argument, help will be shown. 
     
     *e.g.*
     perl tif.pl IRGSP-1.0_genome.fasta ttm5 TGTTAAATATATATACA TTGCAAGTTAGTTAAGA

     This version does not depend on BLAST search. Search script was included in tif.pl.  
     This update is for the multi-core environment.  
- New script tif_nonltr.pl is implemented. (2020-08-19)  
  tif_nonltr.pl detects insertions of non-LTR retrotransposons.  
  perl tif_nonltr.pl  
  will show the usage.  
- New script tif_flanking.pl is implemented. (2019-03-21)  
     tif_flanking is update of tif_basic.pl.  
     If you do not have reference genome sequnce, try tif_flanking.pl.  
     tif_flanking outputs flanking sequence of transposon insertion in fasta format.  
     Run without argument, help will be shown.  
- TIF is a powerful tool.  
     Genomic impact of stress-induced transposable element mobility in *Arabidopsis*  
     David Roquis, Marta Robertson, Liang Yu, Michael Thieme, Magdalena Julkowska, Etienne Bucher  
     Nucleic Acids Res. 49(18):10431-10447 (2021)  
     https://doi.org/10.1093/nar/gkab828  
     
     Sensitive detection of pre-integration intermediates of long terminal repeat retrotransposons in crop plants  
     Jungnam Cho, Matthias Benoit, Marco Catoni, Hajk-Georg Drost, Anna Brestovitsky, Matthijs Oosterbeek, Jerzy Paszkowski  
     Nature Plants, 5:26–33 (2019)  
     https://doi.org/10.1038/s41477-018-0320-9  

     Mobilization of Pack-CACTA transposons in *Arabidopsis* suggests the mechanism of gene shuffling  
     Marco Catoni, Thomas Jonesman, Elisa Cerruti, Jerzy Paszkowski  
     Nucleic Acids Research, 47(3):1311–1320 (2019)  
     https://doi.org/10.1093/nar/gky1196  

### Download TIF
Download zip file of PED from https://github.com/akiomiyao/tif and extract.  

or  

      git clone https://github.com/akiomiyao/tif.git  

If you got scripts from github, update to newest version is very easy using git pull command.  

      git pull  

### Demonstration of tif.pl

For example,  

      perl tif.pl IRGSP-1.0_genome.fasta target TGTTAAATATATATACA TTGCAAGTTAGTTAAGA

or  

      perl tif.pl TAIR10_chr_all.fas target GAGGGATCATCTCTTGTGTC GACTGGCCAGACGATTATTC

or

      perl tif.pl dmel-all-chromosome-r6.29.fasta target CATGATGAAATAACAT ATGTTATTTCATCATG

Before run tif.pl, download fastq file in target/read directory.  
For example, in the case of target name is ttm2  

      cd tif  
      mkdir ttm2  
      mkdir ttm2/read  
      cp somewhere/ttm2.fastq ttm2/read  
      perl tif.pl IRGSP-1.0_genome.fasta ttm2 TGTTAAATATATATACA TTGCAAGTTAGTTAAGA  

Result will be saved to tif_result.head_sequence.tail_sequence and vcf files in the target directory.  

The tif.pl is easy to use and has high sensitivity rather than old programs.  

For *Tos17* retrotransposon of rice

      Head of Tos17: TGTTAAATATATATACA
      Tail of Tos17: TTGCAAGTTAGTTAAGA
      Size of TSD: 5
      fastq: SRR556173 SRR556174 SRR556175
      reference: https://rapdb.dna.affrc.go.jp/download/archive/irgsp1/IRGSP-1.0_genome.fasta.gz

For *mPing* transposon of rice (DNA type transposon)

      Head of mPing: GGCCAGTCACAATGGGG
      Tail of mPing: AGCCATTGTGACTGGCC
      Size of TSD: 3

For *nDart* transposon of rice (DNA type transposon)

      Head of nDart: TAGAGGTGGCCAAACGGGC
      Tail of nDart: GCCCGTTTGGCCACCTCTA
      Size of TSD: 8

For *P*-element of *Drosophila melanogaster*

      Head of P-element: CATGATGAAATAACAT
      Tail of P-element: ATGTTATTTCATCATG
      Size of TSD: 8
      fastq: SRR823377 SRR823382
      reference: ftp://ftp.flybase.net/genomes/Drosophila_melanogaster/dmel_r6.29_FB2019_04/fasta/dmel-all-chromosome-r6.29.fasta.gz

For *Hi* of *Arabidopsis thaliana*

      Head of Hi: GAGGGATCATCTCTTGTGTC
      Tail of Hi: GACTGGCCAGACGATTATTC
      Size of TSD: 9
      doi: https://doi.org/10.1038/emboj.2013.169
      fastq: DRR001193 (ddm1 mutant)
      reference: https://www.arabidopsis.org/download_files/Genes/TAIR10_genome_release/TAIR10_chromosome_files/TAIR10_chr_all.fas

For *ONSEN* of *Arabidopsis thaliana*

      Head of ONSEN: TGTTGAAAGTTAAACTTGAT
      Tail of ONSEN: AAAAGAATTTTACTCTAACA
      Size of TSD: 5
      ONSEN accessions: AT1G11265, AT1G21945, AT1G48710, AT1G58140, AT3G32415, AT3G59720, AT3G61330, AT5G13205
      doi: https://doi.org/10.1093/nar/gkab828
      fastq: https://zenodo.org/record/5052019#.Yd-Fu_7P2Ul
      reference: https://www.arabidopsis.org/download_files/Genes/TAIR10_genome_release/TAIR10_chromosome_files/TAIR10_chr_all.fas
      
### To obtain short read data

Download sra tool kit from  

https://trace.ncbi.nlm.nih.gov/Traces/sra/sra.cgi?view=software

In your home directory,  

      tar xvfz sratoolkit.2.x.x-xxxxxx_linux64.tar.gz
      copy fastq-dump in bin directory to executable directory.

Detail of fastq-dump setup is described in  

https://akiomiyao.github.io/ped/sratoolkit/index.html (English)  

https://akiomiyao.github.io/ped/sratoolkit/index_ja.html (Japanese)

For ttm2 (Rice mutant)

        fastq-dump --split-files -A SRR556173
    
For ttm5 (Rice mutant)

        fastq-dump --split-files -A SRR556174
        fastq-dump --split-files -A SRR556175
    
For *D. melanogaster*

        fastq-dump --split-files -A SRR823377
        fastq-dump --split-files -A SRR823382

for tif.pl,  
fastq files saved in tif/target/read directory are analyzed.  

target name can be changed to your fevorite.  

for old programs,  
fastq files saved in tif/read directory are analyzed.  


### BLAST

BLAST is required by blast.pl and tif_extended.pl.  
New script, tif.pl, does not require BLAST.  

Download BLAST programs 

ftp://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/ncbi-blast-2.9.0+-x64-linux.tar.gz

New version of BLAST can be downloaded from 
https://ftp.ncbi.nlm.nih.gov/blast/executables/
      
Copy blastn and makeblastdb to executable directory.    
     cp ncbi-blast-2.9.0+/bin/blastn ~/bin    
     cp ncbi-blast-2.9.0+/bin/makeblastdb ~/bin  

To make blast data base

      makeblastdb -in reference_genome.fasta -dbtype nucl

For Rice

      cd tif
      wget http://rapdb.dna.affrc.go.jp/download/archive/irgsp1/IRGSP-1.0_genome.fasta.gz
      gzip -d IRGSP-1.0_genome.fasta.gz
      makeblastdb -in IRGSP-1.0_genome.fasta -dbtype nucl
      
For *Drosophira melanogaster*

      cd tif
      wget ftp://ftp.flybase.net/genomes/Drosophila_melanogaster/dmel_r6.29_FB2019_04/fasta/dmel-all-chromosome-r6.29.fasta.gz
      gzip -d dmel-all-chromosome-r6.29.fasta.gz
      makeblastdb -in dmel-all-chromosome-r6.29.fasta -dbtype nucl 


### Search insertions of transposon by tif_basic.pl, blast.pl and tif_extended.pl

Run without any arguments, help message will be shown.

Save fastq files in read directory.

     cd tif
     cp somewhere/foo.fastq read

To test TIF algorithm 1

      cd tif
      perl tif_basic.pl head_sequence tail_sequence TSD_size
      perl blast.pl blatdb_name

      For example,
      cd tif
      perl tif_basic.pl TGTTAAATATATATACA TTGCAAGTTAGTTAAGA 5
      perl blast.pl IRGSP-1.0_genome.fasta
           
Output of tif_basic.pl is tif.fasta, a multiple FASTA file.  
tif.fasta will be saved in tif directory.  
  
The blast.pl reads tif.fasta and returns tif.position containing location and direction of TE insertion sites.

To test TIF algorithm 2

      cd tif
      perl tif_extended.pl reference_fasta_file head_sequence tail_sequence

      For example,
      cd tif
      perl tif_extended.pl IRGSP-1.0_genome.fasta TGTTAAATATATATACA TTGCAAGTTAGTTAAGA
      
The tif_extended.pl returns both tif.fasta and tif.position files in the tif directory.

### Citing TIF

- Cite this article as: Nakagome M, Solovieva E, Takahashi A, Yasue H, Hirochika H, Miyao A (2014) Transposon Insertion Finder (TIF): a novel program for detection of de novo transpositions of transposable elements. BMC Bioinformatics 5:71.  
- doi:10.1186/1471-2105-15-71
- https://bmcbioinformatics.biomedcentral.com/articles/10.1186/1471-2105-15-71


### Update

- 1.7 tif.pl has been improved for multi-core environment. 2021-03-23  
- 1.6 tif_flanking.pl is implemented 2019-03-21.
- 1.5 tif.pl is implemented. 2019-03-19
- 1.4 tif2.pl is improved. 2016-10-22
- 1.3 Add new extended version tif2.pl 2015-03-02
- 1.2 Update README.md 2014-10-09
- 1.1 Update link of SRA-toolkit in README.md 2014-08-01
- 1.0 Inital version 2014-02-05
