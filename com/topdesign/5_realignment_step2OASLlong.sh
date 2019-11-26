picard_tools_path=/home/liurui/software/picard-tools-1.119
GATK_path=/home/liurui/software
D2Lbowtie2index=/home/liurui/databases/OASLlong.fa
inputdatafilesrootpath=/home/liurui/data/bamfiles/OASLref

cmdline=java -Xmx40g -Djava.io.tmpdir=/home/liurui/tmp  -jar $GATK_path/GenomeAnalysisTK.jar -T IndelRealigner  -R $D2Lbowtie2index -I ${1.bam}  -o ${output=/home/liurui/data/bamfiles/bowtie2bam|suffix=sorted.dedup.realigned.bam}
