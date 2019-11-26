picard_tools_path=/home/liurui/software/picard-tools-1.119
GATK_path=/home/liurui/software
inputdatafilesrootpath=/home/liurui/data/bamfiles/beijingduckref
cmdline=java -Xmx40g -Djava.io.tmpdir=/home/liurui/tmp -jar $picard_tools_path/BuildBamIndex.jar I=${uniqmapped.bam}
