picard_tools_path=/home/liurui/software/picard-tools-1.119
GATK_path=/home/liurui/software
D2Lbowtie2index=/home/liurui/databases/D2L_V4.0/D2L.V4.0.fa
inputdatafilesrootpath=/home/liurui/data/bamfiles/bowtie2bam

cmdline=java -Xmx40g -Djava.io.tmpdir=/home/liurui/tmp  -jar $GATK_path/GenomeAnalysisTK.jar -T IndelRealigner  -R $D2Lbowtie2index -I ${sorted.dedup.uniqmapped.bam} -targetIntervals ${.realigner.intervals} -o ${output=/home/liurui/data/bamfiles/bowtie2bam|suffix=sorted.dedup.realigned.0.3uniqmapped.bam}
