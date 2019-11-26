D2Lbowtie2index=/home/liurui/databases/D2L_V4.0/D2L_V4_0.fa

date;echo -ne "bowtie2\n"

inputdatafilesrootpath=/home/liurui/originaldata
cmdline=bowtie2 -p 8 -x $D2Lbowtie2index --rg-id ID --rg-id PL --rg-id PU --rg-id LB --rg-id SM --rg 'ID:${tag}' --rg 'PL:illumina' --rg 'PU:indvd' --rg 'LB:ninglab' --rg 'SM:${tag}' -1 ${1.fq.gz} -2 ${2.fq.gz}|samtools view -@ 8 -bS - -o ${output=/home/liurui/data/bamfiles/bowtie2bam|suffix=bam}
