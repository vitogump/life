# -*- coding: UTF-8 -*-
'''
Created on 2018年6月8日

@author: Dr.liu
'''
from optparse import OptionParser
import os,re,sys

parser = OptionParser()
parser.add_option("-l", "--seqlib", dest="seqlib",default=None,# action="callback",type="string",callback=useoptionvalue_previous2,
                  help="first col corresponding vcf's, second corresponding new. all vcf's chrom is not necessary. transchr only occur in outfile")
parser.add_option("-i","--indgenotypefile",dest="indgenotype",help="default infile1_infile2")

parser.add_option("-o", "--output", dest="output", default=True,
                  help="don't print status messages to stdout")
(options, args) = parser.parse_args()

DIPLOTYPES = ['A', 'C', 'G', 'K',"k","m", 'M', 'N', 'S',"s","r", 'R', 'T', 'W',"w","y", 'Y']
PAIRS = ['AA', 'CC', 'GG', 'GT',"TG","CA", 'AC', 'NN', 'CG',"GC","GA",'AG', 'TT', 'AT',"TA",'TC','CT']
diploHaploDict = dict(zip(DIPLOTYPES,PAIRS))
haploDiploDict = dict(zip(PAIRS,DIPLOTYPES))
tft=["非滩羊","滩羊"]
fo=open(options.output,'w')
#lib
seqlibMapBy4Pos={}
NumlibMapBy4Pos={}
#blind
indnamelist=[]
positionlist=[]
genotypeOfeachInd=[]
#result
judgeIndlist=[]
if __name__ == '__main__':
    #readlib
    libf=open(options.seqlib,'r')
    print(libf.readline())
    for seqcline in libf:
        linelist=re.split(r"\s+",seqcline.strip())
        tetramerPos=tuple(re.split(r"\+",linelist[0]))
        tetramerSeq=""
        for b in re.split(r"\_",linelist[1]):
            tetramerSeq+=haploDiploDict[b.upper()].upper()
            
        if tetramerPos in seqlibMapBy4Pos.keys():
            seqlibMapBy4Pos[tetramerPos].append(tetramerSeq)
            NumlibMapBy4Pos[tetramerPos].append((int(linelist[2]),int(linelist[3])))
        else:
            seqlibMapBy4Pos[tetramerPos]=[tetramerSeq]
            NumlibMapBy4Pos[tetramerPos]=[(int(linelist[2]),int(linelist[3]))]
    print(len(seqlibMapBy4Pos),seqlibMapBy4Pos,sep="\n")
    libf.close()
    
    toBeTestindf=open(options.indgenotype,'r')
    titleline=toBeTestindf.readline()
    #read tobetest ind info
    for pos in re.split(r"\t",titleline.strip())[4:]: positionlist.append(pos.strip())
    for indGenoType in toBeTestindf:
        linelist=re.split(r"\t",indGenoType.strip())
        indnamelist.append(linelist[0].strip())
        seq=[]
        for g in linelist[4:]:
            if len(g.strip())==1:
                seq.append(g.upper())
            elif len(g.strip())==0:
                seq.append("N")
            else:
                seq.append(haploDiploDict[g.strip().upper()])
        seq+=["N"]*(18-len(seq))
        genotypeOfeachInd.append(seq)
    print(*genotypeOfeachInd,sep="\n")
    print(len(genotypeOfeachInd),len(positionlist))
    #start test
    for indidx in reversed(range(len(indnamelist))):
        print(indidx)
        for pos1,pos2,pos3,pos4 in seqlibMapBy4Pos.keys():#each recombination
            pos1idx,pos2idx,pos3idx,pos4idx=positionlist.index(pos1),positionlist.index(pos2),positionlist.index(pos3),positionlist.index(pos4)
            indseq=genotypeOfeachInd[indidx][pos1idx]+genotypeOfeachInd[indidx][pos2idx]+genotypeOfeachInd[indidx][pos3idx]+genotypeOfeachInd[indidx][pos4idx]
            if indseq.upper() in seqlibMapBy4Pos[(pos1,pos2,pos3,pos4)] and sum(NumlibMapBy4Pos[(pos1,pos2,pos3,pos4)][seqlibMapBy4Pos[(pos1,pos2,pos3,pos4)].index(indseq.upper())])>=10 and 0 in NumlibMapBy4Pos[(pos1,pos2,pos3,pos4)][seqlibMapBy4Pos[(pos1,pos2,pos3,pos4)].index(indseq.upper())]:
                print(NumlibMapBy4Pos[(pos1,pos2,pos3,pos4)][seqlibMapBy4Pos[(pos1,pos2,pos3,pos4)].index(indseq.upper())])
                print(indnamelist[indidx],tft[NumlibMapBy4Pos[(pos1,pos2,pos3,pos4)][seqlibMapBy4Pos[(pos1,pos2,pos3,pos4)].index(indseq.upper())].index(0)],"high confidence",file=fo)
                indnamelist.pop(indidx);genotypeOfeachInd.pop(indidx)
                break
            elif indseq.upper() in seqlibMapBy4Pos[(pos1,pos2,pos3,pos4)]:
                print("not",NumlibMapBy4Pos[(pos1,pos2,pos3,pos4)][seqlibMapBy4Pos[(pos1,pos2,pos3,pos4)].index(indseq.upper())])
            else:
                print(pos1,pos2,pos3,pos4,indseq)
    print(indnamelist,len(indnamelist))
    #loser test
    for indidx in reversed(range(len(indnamelist))):
        for pos1,pos2,pos3,pos4 in seqlibMapBy4Pos.keys():#each recombination
            pos1idx,pos2idx,pos3idx,pos4idx=positionlist.index(pos1),positionlist.index(pos2),positionlist.index(pos3),positionlist.index(pos4)
            indseq=genotypeOfeachInd[indidx][pos1idx]+genotypeOfeachInd[indidx][pos2idx]+genotypeOfeachInd[indidx][pos3idx]+genotypeOfeachInd[indidx][pos4idx]
            if indseq.upper() in seqlibMapBy4Pos[(pos1,pos2,pos3,pos4)] :
                print(NumlibMapBy4Pos[(pos1,pos2,pos3,pos4)][seqlibMapBy4Pos[(pos1,pos2,pos3,pos4)].index(indseq.upper())])
                if NumlibMapBy4Pos[(pos1,pos2,pos3,pos4)][seqlibMapBy4Pos[(pos1,pos2,pos3,pos4)].index(indseq.upper())][0]>NumlibMapBy4Pos[(pos1,pos2,pos3,pos4)][seqlibMapBy4Pos[(pos1,pos2,pos3,pos4)].index(indseq.upper())][1]:
                    print(indnamelist[indidx],"滩羊","low confidence",file=fo)
                else:
                    print(indnamelist[indidx],"非滩羊","low confidence",file=fo)
#                 print((pos1,pos2,pos3,pos4),NumlibMapBy4Pos[(pos1,pos2,pos3,pos4)][seqlibMapBy4Pos[(pos1,pos2,pos3,pos4)].index(indseq.upper())],indidx)
                indnamelist.pop(indidx);genotypeOfeachInd.pop(indidx)
                break
        else:
            print(indnamelist[indidx],"非滩羊\tverylowconfidence",file=fo)
    print(indnamelist,len(indnamelist))
    toBeTestindf.close();fo.close()
#             print(pos1,pos2,pos3,pos4)