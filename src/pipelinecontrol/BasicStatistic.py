'''
Created on 2015-1-23

@author: liurui
'''
from optparse import OptionParser
import re, os
parser = OptionParser()

#"output data name is defined as 'inputdatapath folder name'+'is subfolder name'+'is subfolder name'+..."
parser.add_option("-o", "--outputname", dest="outputname",help="outputname")
parser.add_option("-s", "--scriptstorepath", dest="scriptstorepath",help="")
parser.add_option("-l", "--readlength", dest="readlength", help="readlength")
parser.add_option("-g", "--genomelength", dest="genomelength", help="bam bai sam sorted.bam vcf blast and so on. note this is just used in the cmdline output parameter")
parser.add_option("-q", "--quiet",
                  action="store_false", dest="verbose", default=True,
                  help="don't print status messages to stdout")
(options, args) = parser.parse_args()

scriptsstoredir=options.scriptstorepath
readlength=int(options.readlength.strip())
genomelength=int(options.genomelength.strip())
outputfile=open(options.outputname,'w')
totalreads="";averagedeth="";alignmentrate=""
sampleInfoMap={}
if __name__ == '__main__':
    scriptfiles = os.listdir(path=scriptsstoredir)
    for filename in scriptfiles:
        if re.search(r".*Script\.out$",filename)!=None:
            print(filename)
            samplename=re.search(r".*\.([^\.]*)Script\.out$",filename).group(1)
            reportfile=open(scriptsstoredir+"/"+filename,'r')
            for line in reportfile:
                if re.search(r"^\[samopen\] SAM header is present:",line)!=None:
                    try:
                        totalreads=reportfile.readline()
                        print(totalreads,readlength,genomelength)
                        totalreads=int(re.search(r"^(\d+) reads;",totalreads).group(1))
                    except:
                        print("except")
                        continue
                    
                    averagedeth=totalreads*readlength*2/genomelength
            try:
                alignmentrate=float(re.search(r"^([\d\.]*)%",line).group(1).strip())
            except:
                print("sssssssss",line)
                continue
            print(alignmentrate,"alignmentrate*averagedeth",averagedeth)
            sampleInfoMap[samplename]=[totalreads,averagedeth,alignmentrate,alignmentrate*averagedeth/100]
            reportfile.close()
    print("samplename","total readspair","average depth","alignment rate","alignmented depth",sep="\t",file=outputfile)
    for samplename in sorted(sampleInfoMap.keys()):
        print(samplename+"\t"+str(sampleInfoMap[samplename][0])+"\t"+'%.3f'%(sampleInfoMap[samplename][1])+"\t"+str(sampleInfoMap[samplename][2])+"\t"+'%.3f'%(sampleInfoMap[samplename][3]),file=outputfile)
    outputfile.close()       