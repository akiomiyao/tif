#!/bin/sh
#
# A shell script to download NGS data for test
#
#
# Download sratoolkit from
# http://www.ncbi.nlm.nih.gov/Traces/sra/?view=software
# e.g. tar xvfz sratoolkit.2.3.5-2-centos_linux64.tar.gz in your home directory.
# 
~/sratoolkit.2.3.5-2-centos_linux64/bin/fastq-dump --split-files -A SRR823377
~/sratoolkit.2.3.5-2-centos_linux64/bin/fastq-dump --split-files -A SRR823382
~/sratoolkit.2.3.5-2-centos_linux64/bin/fastq-dump --split-files -A SRR556173
~/sratoolkit.2.3.5-2-centos_linux64/bin/fastq-dump --split-files -A SRR556174
~/sratoolkit.2.3.5-2-centos_linux64/bin/fastq-dump --split-files -A SRR556175
