bowtie2_path=/pub/tool/bowtie2-2.1.0
samtools_path=/pub/tool/samtools-0.1.19/bin
fastqc_path=/home/bioinfo/liurui/software/FastQC
duckbowtie2index=/home/bioinfo/liurui/databases/bowtie2idx/duck_1_0_77_genome
picard_tools_path=/home/bioinfo/liurui/software/picard-tools-1.119
GATK_path=/home/liurui/software
duckbowtie2index=/home/liurui/databases/Anas_platyrhynchos.BGI_duck_1.0.dna_sm.toplevel.fa


date;echo -ne "snp calling\n"
inputdatafilesrootpath=/home/liurui/data/bamfiles/beijingduckref
cmdline=java -Xmx40g -Djava.io.tmpdir=/home/liurui/tmp  -jar $GATK_path/GenomeAnalysisTK.jar -T UnifiedGenotyper -R $duckbowtie2index -I ${sorted.group.dedup.realigned.bam} -glm BOTH -o ${output=/home/liurui/data/vcffiles/beijingref|suffix=withindel.vcf}

