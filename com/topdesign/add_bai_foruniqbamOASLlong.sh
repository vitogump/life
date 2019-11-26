picard_tools_path=/home/liurui/software/picard-tools-1.119

inputdatafilesrootpath=/home/liurui/data/bamfiles/OASLref
cmdline=java -Xmx40g -Djava.io.tmpdir=/home/liurui/tmp  -jar  $picard_tools_path/BuildBamIndex.jar I=${sorted.bam}
