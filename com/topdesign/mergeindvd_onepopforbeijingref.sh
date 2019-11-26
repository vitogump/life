inputdatafilesrootpath=/home/liurui/data/bamfiles/beijingduckref
cmdline=java -Xmx40g -Djava.io.tmpdir=/home/liurui/tmp -jar $picard_tools_path/MergeSamFiles.jar ASSUME_SORTED=true INPUT=${dedup.realigned.bam} OUTPUT=${output=/home/liurui/data/bamfiles/beijingduckref|suffix=dedup.realigned.mergesorted.bam} 
