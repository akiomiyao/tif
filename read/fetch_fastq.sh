#!/bin/sh
#
# A shell script to download NGS data for test
#
#
# Download sratoolkit from
# http://www.ncbi.nlm.nih.gov/Traces/sra/?view=software
# e.g. tar xvfz sratoolkit.2.3.5-2-centos_linux64.tar.gz in your home directory.
#      copy bin/fastq-dump to executable directory.
#
fastq-dump --split-files -A SRR823377
fastq-dump --split-files -A SRR823382
fastq-dump --split-files -A SRR556173
fastq-dump --split-files -A SRR556174
fastq-dump --split-files -A SRR556175
