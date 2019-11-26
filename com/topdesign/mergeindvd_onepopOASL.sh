picard_tools_path=/home/liurui/software/picard-tools-1.119
inputdatafilesrootpath=/home/liurui/data/bamfiles/OASLref
cmdline=java -Xmx40g -Djava.io.tmpdir=/home/liurui/tmp -jar $picard_tools_path/MergeSamFiles.jar ASSUME_SORTED=true INPUT=${dedup.bam} OUTPUT=${output=/home/liurui/data/bamfiles/OASLref|suffix=dedup.realigned.mergesorted.bam} 
