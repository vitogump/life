inputdatafilesrootpath=/home/liurui/data/bamfiles/beijingduckref
cmdline=samtools view -h  -@ 8 ${sorted.bam} |awk '$0!~/XS:i/{print $0}'|samtools view -@ 8 -bS - -o ${output=/home/liurui/data/bamfiles/beijingduckref|suffix=sorted.uniqmapped.bam}
