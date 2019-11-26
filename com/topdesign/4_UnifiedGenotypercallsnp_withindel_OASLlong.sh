bowtie2_path=/pub/tool/bowtie2-2.1.0
samtools_path=/pub/tool/samtools-0.1.19/bin
fastqc_path=/home/bioinfo/liurui/software/FastQC
picard_tools_path=/home/bioinfo/liurui/software/picard-tools-1.119
GATK_path=/home/liurui/software
duckbowtie2index=/home/liurui/databases/OASLlong.fa


date;echo -ne "snp calling\n"
inputdatafilesrootpath=/home/liurui/data/bamfiles/OASLref
cmdline=java -Xmx40g -Djava.io.tmpdir=/home/liurui/tmp  -jar $GATK_path/GenomeAnalysisTK.jar -T UnifiedGenotyper -R $duckbowtie2index -I ${mergesorted.bam} -glm BOTH -o ${output=/home/liurui/data/vcffiles/OASLref|suffix=withindel.vcf}

