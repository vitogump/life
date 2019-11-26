GATK_path=/home/liurui/software
duckbowtie2index=/home/liurui/databases/D2L_V4.0/D2L.V4.0.fa
inputdatafilesrootpath=/home/liurui/data
cmdline=java -Xmx40g -Djava.io.tmpdir=/home/liurui/tmp  -jar $GATK_path/GenomeAnalysisTK.jar -T ReadBackedPhasing -R $duckbowtie2index -I ${dedup.realigned.uniqmapped.bam} --variant ${withindel.vcf} -o ${output=/home/liurui/data/vcffiles/D2L|suffix=withindel.phased.vcf}
