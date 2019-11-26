picard_tools_path=/home/liurui/software/picard-tools-1.119
GATK_path=/home/liurui/software
duckbowtie2index=/home/liurui/databases/Anas_platyrhynchos.BGI_duck_1.0.dna_sm.toplevel.fa
inputdatafilesrootpath=/home/liurui/data/bamfiles/beijingduckref

cmdline=java -Xmx40g -Djava.io.tmpdir=/home/liurui/tmp  -jar $GATK_path/GenomeAnalysisTK.jar -T IndelRealigner -allowPotentiallyMisencodedQuals -R $duckbowtie2index -I ${dedup.uniqmapped.bam} -targetIntervals ${.realigner.intervals} -o ${output=/home/liurui/data/bamfiles/beijingduckref|suffix=sorted.sampled.dedup.realigned.uniqmapped.bam}

