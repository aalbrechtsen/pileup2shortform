
g++ -O3 -o pileup2shortform pileup2shortform.cpp


samtools mpileup NA06985.mapped.ILLUMINA.bwa.CEU.low_coverage.20111114.bam | ./pileup2shortform -outfile temp -fai index.fai


samtools mpileup ../test/smallBam/smallNA07000.mapped.ILLUMINA.bwa.CEU.low_coverage.20111114.bam | ./pileup2shortform -outfile temp -fai /home/rowdy/github/angsd/test/hg19.fa.fai

