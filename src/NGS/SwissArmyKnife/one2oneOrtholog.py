'''
Created on 2022年3月28日

@author: RuiLiu
'''
from optparse import OptionParser
import os, re
import pickle

from NGS.BasicUtil import Util


parser = OptionParser()
parser.add_option("-i", "--groupfile", dest="groupfile", help="input group file")
parser.add_option("-s", "--steps", dest="steps", help="1: selected assigend species 1-to-1 gene family into group file\n,2:extract seq.fa from previous out,3,4")
parser.add_option("-c", "--compliantFasta",dest="compliantFasta")
parser.add_option("-n", "--speciesname", dest="speciesname", help="")
(options, args) = parser.parse_args()
if __name__ == '__main__':
    groupfile=open(options.groupfile,"r")
    if options.steps.strip()=="1":
        
        of=open(options.groupfile+".out","w")
        for line in groupfile:
            #linelist=re.split(r"\s+",line.strip())
            if line.count("Lox|")==line.count("GRCh38|")==line.count("MMUL|")==line.count("PapA|")==line.count("GSM|")==line.count("Cjac|")==line.count("CanF|")==line.count("Pika|")==line.count("PIG|")==line.count("turT|")==line.count("GRCm|")==line.count("UMD|")==line.count("Equ|")==line.count("Oar|")==1:
                print(line,end="",file=of)
            else:
                print(line[:20])
        of.close()
    elif options.steps.strip()=="2":
        try:
            os.mkdir(options.groupfile.replace(".","_"))
        except:
            pass
        fileHmap={"Lox":0,"GRCh38":0,"MMUL":0,"PapA":0,"GSM":0,"Cjac":0,"CanF":0, "Pika":0,"PIG":0, "turT":0 , "GRCm":0, "UMD":0, "Equ":0, "Oar":0}
        for fn in os.listdir(options.compliantFasta):
            pathTofn=os.path.join(options.compliantFasta.strip(),fn.strip())
            print("build/read indexfile for:",fn)
            if fn.replace(".fasta","") in fileHmap and os.path.isfile(pathTofn):
                try:
                    duckrefindex = pickle.load(open(pathTofn + ".myfasteridx", 'rb'))
            #             originalspeciesindex = pickle.load(open(originalspeciesref + ".myindex", 'rb'))
                except IOError:
                    Util.generateFasterRefIndex(pathTofn, pathTofn + ".myfasteridx",chrsignal="|")
                    duckrefindex = pickle.load(open(pathTofn + ".myfasteridx", 'rb'))                 
                fileHmap[fn.replace(".fasta","")]=[open(pathTofn,'r'),duckrefindex]
                print()
            else:
                print("pasd",fn,pathTofn)
        for line in groupfile:
            linelist= re.split(r"\s+",line.strip())
            for nameOf_s_p in linelist[1:]:
                of_tmp=open(os.path.join(options.groupfile.replace(".","_"),linelist[0][1:]+".fa"),"w")
                s_p=re.search(r"(^[\d\w]+)\|(.*)$",nameOf_s_p)
                if s_p.group(1) not in fileHmap:
                    print("skip seq:",s_p.group(1),s_p.group(2))
                    continue
                RefSeqMap=Util.getRefSeqBypos_faster(fileHmap[s_p.group(1)][0], fileHmap[s_p.group(1)][1], s_p.group(2), 1, 9999999999)
                
                print(">"+s_p.group(0)+"\n"+"".join(RefSeqMap[s_p.group(2)][1:]),file=of_tmp)
                of_tmp.close()
        
 
        for fn in fileHmap.keys():
            fileHmap[fn][0].close()
    
    groupfile.close()