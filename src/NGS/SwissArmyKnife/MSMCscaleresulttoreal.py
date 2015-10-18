'''
Created on 2015-8-22

@author: liurui
'''
from optparse import OptionParser
import re
parser = OptionParser()



#scriptDir,mode="series",logfile
parser.add_option("-g", "--generation_time", dest="generation_time",help="scriptDir")
parser.add_option("-u", "--mutationrate", dest="mutationrate",help="oneline scriptexamplefile")
parser.add_option("-i","--inputfilename",dest="inputfilename",help="a little note message")
# parser.add_option("-l", "--logfile", dest="logfile", help="bam bai sam sorted.bam vcf blast and so on. note this is just used in the cmdline output parameter")
parser.add_option("-o", "--outputfilename", dest="outputfilename",help="p:parallel s:series")

                                                                                                                                                          
(options, args) = parser.parse_args()
infile=open(options.inputfilename,'r')
outfile=open(options.outputfilename,'w')
mu=float(options.mutationrate)
print("mu",mu)
g=float(options.generation_time)
if __name__ == '__main__':
    title=infile.readline().strip()
    if len(re.split(r"\t",title))==4:
        print(title.strip()+"\trealstarttime\trealendtime\trealeffectivepopsize",file=outfile)
        for line in infile:
            linelist=re.split(r"\s+",line.strip())
            print(*linelist,sep="\t",end="\t",file=outfile)
            try:
                print(str((float(linelist[1])*g)/(mu/2)),str((float(linelist[2])*g)/(mu/2)),str((1/float(linelist[3]))/mu),sep="\t",file=outfile)
            except:
                print(linelist)
    elif len(re.split(r"\t",title))==6:
        print(title.strip()+"\trealstarttime\trealendtime\trelative_cross_coalescence_rate",file=outfile)
        for line in infile:
            linelist=re.split(r"\s+",line.strip())
            print(*linelist,sep="\t",end="\t",file=outfile)
            print(str((float(linelist[1])*g)/(mu/2)),str((float(linelist[2])*g)/(mu/2)),str(((2*float(linelist[4]))/(float(linelist[3])+float(linelist[5])))),sep="\t",file=outfile)
    outfile.close()
    infile.close()