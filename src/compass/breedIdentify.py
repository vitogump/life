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
parser.add_option("-i","--indgenotypefile",dest="indgenotype",nargs=2,help="default infile1_infile2")
parser.add_option("-u","--update",dest="identifiedresultoupdatelib",nargs=2,default=None,help="if this option exist, then updatelib, otherwise identify samples")
parser.add_option("-o", "--output", dest="output", default=True,
                  help="don't print status messages to stdout")
(options, args) = parser.parse_args()

No_ofsites=18
DIPLOTYPES = ['A', 'C', 'G', 'K',"k","m", 'M', 'N', 'S',"s","r", 'R', 'T', 'W',"w","y", 'Y']
PAIRS = ['AA', 'CC', 'GG', 'GT',"TG","CA", 'AC', 'NN', 'CG',"GC","GA",'AG', 'TT', 'AT',"TA",'TC','CT']
diploHaploDict = dict(zip(DIPLOTYPES,PAIRS))
haploDiploDict = dict(zip(PAIRS,DIPLOTYPES))
tft=["非滩羊","滩羊"]

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
    libtitle=libf.readline()
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
    print(len(seqlibMapBy4Pos),*seqlibMapBy4Pos,sep="\n")
    for pos1,pos2,pos3,pos4 in seqlibMapBy4Pos.keys():
        print(pos1,pos2,pos3,pos4)
        print(seqlibMapBy4Pos[(pos1,pos2,pos3,pos4)])
        print(NumlibMapBy4Pos[(pos1,pos2,pos3,pos4)])
    libf.close()
#     exit()
    toBeTestindf=open(options.indgenotype[0],'r')
    print(options.indgenotype)
    titleline=toBeTestindf.readline()
    #read tobetest ind info
    for pos in re.split(r"\t",titleline.strip())[int(options.indgenotype[1])-1:]: positionlist.append(pos.strip())
    for indGenoType in toBeTestindf:
        linelist=re.split(r"\t",indGenoType.strip())
        indnamelist.append(linelist[0].strip())
        seq=[]
        for g in linelist[int(options.indgenotype[1])-1:]:
            if len(g.strip())==1:
                seq.append(g.upper())
            elif len(g.strip())==0:
                seq.append("N")
            else:
                seq.append(haploDiploDict[g.strip().upper()])
        seq+=["N"]*(No_ofsites-len(seq))
        genotypeOfeachInd.append(seq)
    print(*genotypeOfeachInd,sep="\n")
    print(positionlist)
    print(len(genotypeOfeachInd),len(positionlist))
    if options.identifiedresultoupdatelib!=None:
        answerfile=open(options.identifiedresultoupdatelib[0],'r')
#         answertitle=answerfile.readline()
        for ans in answerfile:
            anslist=re.split(r"\t",ans.strip())
            if anslist[0].strip() in indnamelist:
                indidx=indnamelist.index(anslist[0].strip())#used to find genotype
                for pos1,pos2,pos3,pos4 in seqlibMapBy4Pos.keys():
                    pos1idx,pos2idx,pos3idx,pos4idx=positionlist.index(pos1),positionlist.index(pos2),positionlist.index(pos3),positionlist.index(pos4)
                    indseq=genotypeOfeachInd[indidx][pos1idx]+genotypeOfeachInd[indidx][pos2idx]+genotypeOfeachInd[indidx][pos3idx]+genotypeOfeachInd[indidx][pos4idx]
                    if indseq.upper() in seqlibMapBy4Pos[(pos1,pos2,pos3,pos4)]:
                        i=seqlibMapBy4Pos[(pos1,pos2,pos3,pos4)].index(indseq.upper())
                        NoOfT,NoOfNT=NumlibMapBy4Pos[(pos1,pos2,pos3,pos4)][i]
                        if len(anslist)>=int(options.identifiedresultoupdatelib[1]) and anslist[int(options.identifiedresultoupdatelib[1])-1].lower()=="y":
                            NoOfT+=1
                        elif len(anslist)>=int(options.identifiedresultoupdatelib[1]) and anslist[int(options.identifiedresultoupdatelib[1])-1].lower()=="n":
                            NoOfNT+=1
                        NumlibMapBy4Pos[(pos1,pos2,pos3,pos4)][i]=(NoOfT,NoOfNT)
                        print("update",anslist[0].strip())
                    elif "N" not in indseq.upper():
                        seqlibMapBy4Pos[(pos1,pos2,pos3,pos4)].append(indseq.upper())
                        if len(anslist)>=int(options.identifiedresultoupdatelib[1]) and anslist[int(options.identifiedresultoupdatelib[1])-1].lower()=="y":
                            NumlibMapBy4Pos[(pos1,pos2,pos3,pos4)].append((1,0))
                        elif len(anslist)>=int(options.identifiedresultoupdatelib[1]) and anslist[int(options.identifiedresultoupdatelib[1])-1].lower()=="n":
                            NumlibMapBy4Pos[(pos1,pos2,pos3,pos4)].append((0,1))
        print("writing new lib...")
        nlibf=open(options.output+".lib",'w')
        print(libtitle.strip(),file=nlibf)
        for pos1,pos2,pos3,pos4 in seqlibMapBy4Pos.keys():
            for fourPosBase in seqlibMapBy4Pos[(pos1,pos2,pos3,pos4)]:
                print(pos1,pos2,pos3,pos4,sep="+",end="\t",file=nlibf)
                print(diploHaploDict[fourPosBase[0]],diploHaploDict[fourPosBase[1]],diploHaploDict[fourPosBase[2]],diploHaploDict[fourPosBase[3]],sep="_",end="\t",file=nlibf)
                print(*NumlibMapBy4Pos[(pos1,pos2,pos3,pos4)][seqlibMapBy4Pos[(pos1,pos2,pos3,pos4)].index(fourPosBase)],sep="\t",file=nlibf)
        nlibf.close()
        print("writing new lib done",options.output+".lib")
        exit()
    #start test
    fo=open(options.output,'w')
    for indidx in reversed(range(len(indnamelist))):
        #exclude all combination has a N  
        for pos1,pos2,pos3,pos4 in seqlibMapBy4Pos.keys():
            pos1idx,pos2idx,pos3idx,pos4idx=positionlist.index(pos1),positionlist.index(pos2),positionlist.index(pos3),positionlist.index(pos4)
            indseq=genotypeOfeachInd[indidx][pos1idx]+genotypeOfeachInd[indidx][pos2idx]+genotypeOfeachInd[indidx][pos3idx]+genotypeOfeachInd[indidx][pos4idx]      
            if "N" in indseq.upper():
                continue
            else:
                break
        else:
            print(indnamelist.pop(indidx),"can't determine",sep="\t",file=fo)
            print(genotypeOfeachInd.pop(indidx))
        #test for each 4sites combination
        for pos1,pos2,pos3,pos4 in seqlibMapBy4Pos.keys():
            pos1idx,pos2idx,pos3idx,pos4idx=positionlist.index(pos1),positionlist.index(pos2),positionlist.index(pos3),positionlist.index(pos4)
            indseq=genotypeOfeachInd[indidx][pos1idx]+genotypeOfeachInd[indidx][pos2idx]+genotypeOfeachInd[indidx][pos3idx]+genotypeOfeachInd[indidx][pos4idx]
            if indseq.upper() in seqlibMapBy4Pos[(pos1,pos2,pos3,pos4)] and sum(NumlibMapBy4Pos[(pos1,pos2,pos3,pos4)][seqlibMapBy4Pos[(pos1,pos2,pos3,pos4)].index(indseq.upper())])>=10 and 0 in NumlibMapBy4Pos[(pos1,pos2,pos3,pos4)][seqlibMapBy4Pos[(pos1,pos2,pos3,pos4)].index(indseq.upper())]:
                print(NumlibMapBy4Pos[(pos1,pos2,pos3,pos4)][seqlibMapBy4Pos[(pos1,pos2,pos3,pos4)].index(indseq.upper())])
                print(indnamelist[indidx],tft[NumlibMapBy4Pos[(pos1,pos2,pos3,pos4)][seqlibMapBy4Pos[(pos1,pos2,pos3,pos4)].index(indseq.upper())].index(0)],"highconfidence",sep="\t",file=fo)
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
                    print(indnamelist[indidx],"滩羊","lowconfidence",sep="\t",file=fo)
                else:
                    print(indnamelist[indidx],"非滩羊","lowconfidence",sep="\t",file=fo)
#                 print((pos1,pos2,pos3,pos4),NumlibMapBy4Pos[(pos1,pos2,pos3,pos4)][seqlibMapBy4Pos[(pos1,pos2,pos3,pos4)].index(indseq.upper())],indidx)
                indnamelist.pop(indidx);genotypeOfeachInd.pop(indidx)
                break
        else:
            print(indnamelist[indidx],"非滩羊\tverylowconfidence",sep="\t",file=fo)
    print(indnamelist,len(indnamelist))
    toBeTestindf.close();fo.close()
#             print(pos1,pos2,pos3,pos4)