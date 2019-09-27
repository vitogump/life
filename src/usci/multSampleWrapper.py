'''
Created on 2019年9月18日

@author: liurui
'''

from optparse import OptionParser
import os,re

parser = OptionParser()
parser.add_option("-1", "--p1fq", dest="p1fq", help="")
parser.add_option("-2", "--p2fq", dest="p2fq",default=None, help="this is an optional, if present, then use the seed as -1")
parser.add_option("-t", "--timestosample", dest="timestosample",default="25", help="")
parser.add_option("-s", "--startseed", dest="startseed")
parser.add_option("-c", "--sampletoreadscount", dest="sampletoreadscount")
parser.add_option("-o", "--sampleedfq1", dest="sampleedfq1")
parser.add_option("-O", "--sampleedfq2", dest="sampleedfq2")

parser.add_option("-i", "--insertindex", dest="insertindex",default=0,help="this parameter indicate the lenght of the suffix,to insert roundNumber before it ")

parser.add_option("-C", "--checkcounts", dest="checkcounts",action="store_true",default=False)
parser.add_option("-T", "--sampleTools", dest="sampleTools",default="seqtk sample")
(options, args) = parser.parse_args()
if __name__ == '__main__':
    seed=int(options.startseed);iji=int(options.insertindex)
    if options.checkcounts: c=os.popen("less "+options.p1fq+"|wc -l");print(options.p1fq,"reads counts:",int(c.readline().strip())/4)
    j=0
    for i in range(seed,seed+int(options.timestosample)):
        j+=1
        cmd=options.sampleTools+" -s"+str(i)+" "+options.p1fq+" "+str(options.sampletoreadscount)+" > "+options.sampleedfq1[:-iji]+"_"+str(j)+"_"+options.sampleedfq1[-iji:]
        print(cmd)
        if os.path.isfile(options.sampleedfq1+str(j)):
            print(options.sampleedfq1,"exist")
            continue
        a=os.system(cmd)
        if options.checkcounts:c=os.popen("wc -l "+options.sampleedfq2[:-iji]+"_"+str(j)+"_"+options.sampleedfq2[-iji:]);print(options.sampleedfq1+str(j),"reads counts:",int(re.split(r'\s+',c.readline().strip())[0])/4)
        
        if a==0:print("seed ",i,"for "+options.sampleedfq1+" sample successfully!")
        if options.p2fq:
            cmd=options.sampleTools+" -s"+str(i)+" "+options.p2fq+" "+str(options.sampletoreadscount)+" > "+options.sampleedfq2[:-iji]+"_"+str(j)+"_"+options.sampleedfq2[-iji:]
            a=os.system(cmd)
            if options.checkcounts:c=os.popen("wc -l "+options.sampleedfq2[:-iji]+"_"+str(j)+"_"+options.sampleedfq2[-iji:]);print(options.sampleedfq2+str(j),"reads counts:",int(re.split(r'\s+',c.readline().strip())[0])/4)
            if a==0:print("seed ",i,"for "+options.sampleedfq2+" sample successfully!")