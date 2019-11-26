picard_tools_path=/home/liurui/software/picard-tools-1.119

inputdatafilesrootpath=/home/liurui/data/bamfiles/bowtie2bam
cmdline=java -Xmx40g -Djava.io.tmpdir=/home/liurui/tmp  -jar  $picard_tools_path/BuildBamIndex.jar I=${uniqmapped.bam}
