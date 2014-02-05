#!/bin/sh
#
# A shell script to download NGS data for test
#
#
# Download sratoolkit from
# http://trace.ncbi.nlm.nih.gov/Traces/sra/sra.cgi?cmd=show&f=software&m=software&s=software
# e.g. tar xvfz sratoolkit.2.3.4-2-centos_linux64.tar.gz in your home directory.
# 
~/sratoolkit.2.3.4-2-centos_linux64/bin/fastq-dump --split-files -A SRR823377
~/sratoolkit.2.3.4-2-centos_linux64/bin/fastq-dump --split-files -A SRR823382
~/sratoolkit.2.3.4-2-centos_linux64/bin/fastq-dump --split-files -A SRR556173
~/sratoolkit.2.3.4-2-centos_linux64/bin/fastq-dump --split-files -A SRR556174
~/sratoolkit.2.3.4-2-centos_linux64/bin/fastq-dump --split-files -A SRR556175
