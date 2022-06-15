# -*- coding: UTF-8 -*-
'''
Created on 2013-7-2

@author: rui
'''

from itertools import combinations
import random
import re, copy, math,  time, pysam
import numpy as np
from NGS.BasicUtil import VCFutil

seqerrorrate=0.008
outgidx=7;outg2idx=9
class Caculator():
    def __init__(self):
        #every Caculator which need two or more vcf have the follow two variable
        self.vcfnamelist=[]
        self.vcfnameKEY_vcfobj_pyBAMfilesVALUE={}
    def process(self, T):
        pass
    def getResult(self):
        pass
class Caculator_selectsitesForillumina():
    def __init__(self,of,tilingorderidx,best_recommendation,T,VIPidx):
        #every Caculator which need two or more vcf have the follow two variable
        self.vcfnamelist=[]
        self.vcfnameKEY_vcfobj_pyBAMfilesVALUE={}
        self.contained=[]
        self.VIPidx=VIPidx+1
        self.m={"FORWARD":"Plus","REVERSE":"Minus"}
        self.greaterTsort={"1":[],"2":[],"3":[],"4":[],"5":[],"6":[],"7":[]}# insertSort
        self.source_version=tilingorderidx+1#first element is position info
        self.rcmidx=best_recommendation+1#score
        self.classedrecs=[]
        self.lowthanT=[]#lower than threshold
        self.Tvalue=T
        self.pf=of
        self.spf=None
        self.count=0
    def process(self, T):
        self.count+=1
        self.classedrecs.append(T)
        if T[self.VIPidx]=="VIP":
            print(T,file=open("vipppp",'a'))
            self.contained.append(len(self.classedrecs)-1)
        if float(T[self.rcmidx])<self.Tvalue:
            #only make sure first value is the bigest
            if self.lowthanT!=[] and float(T[self.rcmidx])>float(self.classedrecs[self.lowthanT[0]][self.rcmidx]):
                self.lowthanT.insert(0,len(self.classedrecs)-1)
            else:
                self.lowthanT.append(len(self.classedrecs)-1)
        else:
            if self.greaterTsort[T[self.source_version]]!=[] and float(T[self.rcmidx])>float(self.classedrecs[self.greaterTsort[T[self.source_version]][0]][self.rcmidx]):
                self.greaterTsort[T[self.source_version]].insert(0,len(self.classedrecs)-1)
            else:
                self.greaterTsort[T[self.source_version]].append(len(self.classedrecs)-1)
    def getResult(self):
#         print(self.count)
#         print(self.restATCT,self.restNonATCT,self.contained,len(self.classedrecs))
        
        if len(self.contained)!=0:
            for recidx in reversed(self.contained):
                e=self.classedrecs.pop(recidx)
                print(*e[1:],sep="\t",file=self.spf)
                print(e[1],"SNP",e[2],e[4],e[5],e[3],e[6],e[7],e[8],self.m[e[8].upper()],"Soybean","FALSE",sep=",",file=self.pf)
        else:
            for k in ["1","2","3"]:#,"4","5","6","7"]:#["2","3","4","5"]
                if self.greaterTsort[k]!=[] :#
                    e=self.classedrecs.pop(self.greaterTsort[k][0])
                    print(*e[1:],sep="\t",file=self.spf)
                    print(e[1],"SNP",e[2],e[4],e[5],e[3],e[6],e[7],e[8],self.m[e[8].upper()],"Soybean","FALSE",sep=",",file=self.pf)
                    break
            else:
                e_value=0#score of the site
                for k in ["4","5","6"]:
                    if self.greaterTsort[k]!=[] and float(self.classedrecs[self.greaterTsort[k][0]][self.rcmidx])>e_value:
                        eidx=self.greaterTsort[k][0];e_value=float(self.classedrecs[self.greaterTsort[k][0]][self.rcmidx])
                if e_value!=0:#found the max socre in 4 5 6
                    e=self.classedrecs.pop(eidx)
                    print(*e[1:],sep="\t",file=self.spf)
                    print(e[1],"SNP",e[2],e[4],e[5],e[3],e[6],e[7],e[8],self.m[e[8].upper()],"Soybean","FALSE",sep=",",file=self.pf)
                elif self.lowthanT!=[]:# all score are less than self.Tvalue
                    e=self.classedrecs.pop(self.lowthanT[0])
                    print(*e[1:],sep="\t",file=self.spf)
                    print(e[1],"SNP",e[2],e[4],e[5],e[3],e[6],e[7],e[8],self.m[e[8].upper()],"Soybean","FALSE",sep=",",file=self.pf)
        self.contained=[]
        self.greaterTsort={"1":[],"2":[],"3":[],"4":[],"5":[],"6":[],"7":[]}# insertSort
        self.pf.flush()
        self.classedrecs=[]
        self.lowthanT=[]#lower than threshold
        return 1,2
class Caculator_addpriority():
    def __init__(self,of,tilingorderidx=2,best_recommendation=32,rmdupmap={},mustin=set()):
        #every Caculator which need two or more vcf have the follow two variable
        self.vcfnamelist=[]
        self.vcfnameKEY_vcfobj_pyBAMfilesVALUE={}
        self.contained=[]
        self.mmm={"2":"2","3":"3","4":"4","5":"1"}
        self.restNonATCT={"1":[],"2":[],"3":[],"4":[],"5":[],"6":[],"7":[],"10":[]}#{"2":[],"3":[],"4":[],"5":[]}# "10":[] for records ignoring to add in chips
        self.restATCT={"1":[],"2":[],"3":[],"4":[],"5":[],"6":[],"7":[],"10":[]}#{"2":[],"3":[],"4":[],"5":[]}
        self.tidx=tilingorderidx+1#first element is position info
        self.rcmidx=best_recommendation+1
        self.rmdupmap=rmdupmap
        self.classedrecs=[]
        self.notrecommand=[]
        self.pf=of
        self.count=0   
        self.temp=open("dupidnotequal.txt",'w')
        self.musctincludes=mustin
    def process(self, T):
        self.count+=1
        seql=re.split(r"\[.+\]",T[2].strip())
        if seql[0]+seql[1]  in self.rmdupmap:
            if len(set(self.rmdupmap[seql[0]+seql[1]][1:]))>1:
                print(T,file=self.temp)
                return
            else:
                T[self.tidx]==str(self.rmdupmap[seql[0]+seql[1]][0])
        if T[1] in self.musctincludes:
            self.classedrecs.append(T)
            self.contained.append(len(self.classedrecs)-1)
            return

        if (T[self.rcmidx]=="not_recommended" and  T[self.tidx].strip()!="1") or T[self.rcmidx]=="not_possible":#or T[-1]=="not_recommended"
            self.notrecommand.append(T)
        else:
#             print("should not be here",T)
            self.classedrecs.append(T)
            if T[self.tidx].strip()=="1": #and T[self.rcmidx]=="recommended":#"1"
                self.contained.append(len(self.classedrecs)-1)
            elif "A/T" in T[2] or "C/G" in T[2]:
                if T[self.rcmidx]=="recommended":#keep the recommended first
                    self.restATCT[T[self.tidx].strip()].insert(0,len(self.classedrecs)-1)
                else:
                    self.restATCT[T[self.tidx].strip()].append(len(self.classedrecs)-1)
            else:
                if T[self.rcmidx]=="recommended":
                    self.restNonATCT[T[self.tidx].strip()].insert(0,len(self.classedrecs)-1)
                else:
                    self.restNonATCT[T[self.tidx].strip()].append(len(self.classedrecs)-1)
    def getResult(self):
#         print(self.count)
#         print(self.restATCT,self.restNonATCT,self.contained,len(self.classedrecs))

        if len(self.contained)!=0:
            for recidx in reversed(self.contained):
                e=self.classedrecs.pop(recidx)
                if "unaln" not in e[1]:e[self.addVIPidx+1]="VIP"
                if "unaln" in e[1] and e[self.rcmidx]=="not_recommended":
                    print(*e[1:],"9",sep="\t",file=self.pf);continue
                print(*e[1:],"1",sep="\t",file=self.pf)
#             for e in self.classedrecs:
#                 if e[self.rcmidx]=="recommended":
#                     print(*e[1:],"8",sep="\t",file=self.pf)
#                 else:
#                     print(*e[1:],"9",sep="\t",file=self.pf)
        else:
            for k in ["1","2","3","4","5","6","7"]:#["2","3","4","5"]
                if self.restNonATCT[k]!=[] and self.classedrecs[self.restNonATCT[k][0]][self.rcmidx]=="recommended":#
                    e=self.classedrecs.pop(self.restNonATCT[k][0])
                    print(*e[1:],k,sep="\t",file=self.pf)#self.mmm[k]
                    break
            else:
                for k in ["1","2","3","4","5","6","7"]:
                    if self.restNonATCT[k]!=[]:
                        e=self.classedrecs.pop(self.restNonATCT[k][0])
                        print(*e[1:],k,sep="\t",file=self.pf)#self.mmm[k]
                        break
                else:
                    for k in ["1","2","3","4","5","6","7"]:
                        if self.restATCT[k]!=[] :
                            e=self.classedrecs.pop(self.restATCT[k][0])
                            print(*e[1:],k,sep="\t",file=self.pf)#self.mmm[k]
                            break
        for e in self.classedrecs:
            if e[self.rcmidx]=="recommended":
                print(*e[1:],"8",sep="\t",file=self.pf)
            else:
                print(*e[1:],"9",sep="\t",file=self.pf)
        for ne in self.notrecommand:
            print(*ne[1:],"9",sep="\t",file=self.pf)
        self.contained=[]
        self.restNonATCT={"1":[],"2":[],"3":[],"4":[],"5":[],"6":[],"7":[],"10":[]}#{"2":[],"3":[],"4":[],"5":[]}
        self.restATCT={"1":[],"2":[],"3":[],"4":[],"5":[],"6":[],"7":[],"10":[]}#{"2":[],"3":[],"4":[],"5":[]}
        self.pf.flush()
        self.classedrecs=[]
        self.notrecommand=[]
        return 1,2
#             self.markedrecs
class Caculate_ddaf(Caculator):
    def __init__(self,delta_DerAftotal,absdelta_DerAftotal):
        super().__init__()
        self.count = 0
        self.ddaf = 0
        self.absddaf=0
        self.delta_DerAftotal=copy.deepcopy(delta_DerAftotal)
        self.absdelta_DerAftotal=copy.deepcopy(absdelta_DerAftotal)
#         self.T=""
    def process(self, T):
        self.count+=1
        self.ddaf += float(T[4])-float(T[5])
        self.absddaf += abs(float(T[4])-float(T[5]))
#         self.T=T
    def getResult(self):
        c=self.count
        absddaf=self.absddaf
        ddaf=self.ddaf
        self.count = 0;self.ddaf = 0;self.absddaf=0
#         print(c,[absddaf,ddaf],self.T)
        return c,[absddaf,ddaf]
#Bits for intyerpreting and manipulating sequence data

DIPLOTYPES = ['A', 'C', 'G', 'K', 'M', 'N', 'S', 'R', 'T', 'W', 'Y']
PAIRS = ['AA', 'CC', 'GG', 'GT', 'AC', 'NN', 'CG', 'AG', 'TT', 'AT', 'CT']
diploHaploDict = dict(zip(DIPLOTYPES,PAIRS))
haploDiploDict = dict(zip(PAIRS,DIPLOTYPES))

def haplo(diplo): return diploHaploDict[diplo]

def diplo(pair): return haploDiploDict[pair]


#convert one ambiguous sequence into two haploid pseudoPhased sequences

def pseudoPhase(sequence, seqType = "diplo"):
    if seqType == "diplo": pairs = [haplo(s) for s in sequence]
    else: pairs = sequence
    return [[p[0] for p in pairs], [p[1] for p in pairs]]


################################################################################################################

#modules for working with and analysing alignments

numSeqDict = {"A":0,"C":1,"G":2,"T":3,"N":np.NaN}
SeqnumDict = {0:"A",1:"C",2:"G",3:"T"}
def numHamming(numArrayA, numArrayB):
    dif = numArrayA - numArrayB
    return np.nanmean(dif[~np.isnan(dif)] != 0)
def distMatrix(numArray):
    N,l = numArray.shape
    distMat = np.zeros((N,N))
    for i in range(N - 1):
        for j in range(i + 1, N):
            distMat[i,j] = distMat[j,i] = numHamming(numArray[i,:], numArray[j,:])
    return distMat
class Caculate_popDiv(Caculator):
    def __init__(self,considerINDEL,tvcfconfig,rvcfconfig,outputname):
        super().__init__()
        self.outputname=outputname
        self.considerINDEL=considerINDEL
        self.indnamesOfEachPop=[]
        self.MethodToSeqpoplist=[]
        self.listOfpopvcfRecsmapByAChr=[]
        for vcfconfigf in [tvcfconfig[0],rvcfconfig[0]]:
            vcfconfig=open(vcfconfigf,"r")
            self.listOfpopvcfRecsmapByAChr.append({})
            for line in vcfconfig:
                vcffilename_obj=re.search(r"vcffilename=(.*)",line.strip())
                if vcffilename_obj!=None:
                    vcfname=vcffilename_obj.group(1).strip()
                    self.vcfnamelist.append(vcfname)
                    self.outputname+=("_"+re.split(r"\.",re.search(r"[^/]*$",vcfname).group(0))[0])[:3]
                    self.vcfnameKEY_vcfobj_pyBAMfilesVALUE[vcfname]=[]
                    self.vcfnameKEY_vcfobj_pyBAMfilesVALUE[vcfname].append(VCFutil.VCF_Data(vcfname))
                    self.indnamesOfEachPop.append(self.vcfnameKEY_vcfobj_pyBAMfilesVALUE[vcfname][0].VcfIndexMap["title"][9:])
                elif line.split():
                    self.vcfnameKEY_vcfobj_pyBAMfilesVALUE[vcfname].append(pysam.Samfile(line.strip(),'rb'))
            vcfconfig.close()
            if re.search(r"indvd[^/]+",vcfname)!=None:
                self.MethodToSeqpoplist.append("indvd")
    
            elif re.search(r"pool[^/]+",vcfname)!=None:
                self.MethodToSeqpoplist.append("pool")
    
            else:
                print("vcfname must with 'pool' or 'indvd'")
                exit(-1)    
            
        self.minAN=40;self.minAC=2
        self.GQthreshold=30
        self.DPindthreshold=6
        print(self.indnamesOfEachPop)
        self.N=len(self.indnamesOfEachPop[0])+len(self.indnamesOfEachPop[1])
        self.popIndices=[[],[]]
        for ind_idx in range(len(self.indnamesOfEachPop[0])):
            self.popIndices[0].append(2*ind_idx)
            self.popIndices[0].append(2*ind_idx+1)
        for ind_idx in range(len(self.indnamesOfEachPop[0]),len(self.indnamesOfEachPop[0])+len(self.indnamesOfEachPop[1])):
            self.popIndices[1].append(2*ind_idx)
            self.popIndices[1].append(2*ind_idx+1)
        print(self.popIndices)
        #below variable should be changed every win
        self.positions=[]
        self.seqs=[[] for e in range(self.N)]# pop1_inds pop2_inds
        self.numArray=[]
        self.array=[]
    def process(self,T):
        """T is like (pos,REF,ALT,(INFO,FORMAT,sampleslist),(INFO,FORMAT,sampleslist)) should have only two pop
        """
        if self.considerINDEL == "no" and (len(T[1]) != 1 or len(T[2]) != 1):
            return
        site=[]#should have 
        AN=AC=0
        for popidx in range(3,5):
            if T[popidx]==None:
                # check depth ,if passed treat as fix as ref
                if len(self.vcfnameKEY_vcfobj_pyBAMfilesVALUE[self.vcfnamelist[popidx-3]])==1:
#                     print("skip this pos",T)
                    return
                else:
#                     depth_linelist=self.depthobjlist[tpopidx-3].getdepthByPos_optimized(self.currentchrID,T[0])
                    sum_depth=0
                    for samfile in self.vcfnameKEY_vcfobj_pyBAMfilesVALUE[self.vcfnamelist[popidx-3]][1:]:
                        ACGTdep=samfile.count_coverage(self.currentchrID,T[0]-1,T[0])
                        for dep in ACGTdep:
                            sum_depth+=dep[0]

                    if sum_depth>self.DPindthreshold*len(self.indnamesOfEachPop[popidx-3])*1.2:#use 1.2 is ribitrarily, beacuse we don't have GQ to filter some ind so use more strength thread in depth
                        site+=[T[1].upper()*2 for x in range(len(self.indnamesOfEachPop[popidx-3]))]
                        AN+=(len(self.indnamesOfEachPop[popidx-3])*2)
                    else:
                        return
            else:
                if self.MethodToSeqpoplist[popidx-3]=="indvd":
#                     AN += int(re.search(r"AN=(\d+)[;,]", T[popidx][0]).group(1))
                    AC += int(re.search(r"AC=(\d+)[;,]", T[popidx][0]).group(1))
                    GT_idx = (re.split(":", T[popidx][1])).index("GT")
                    GQ_idx=(re.split(":", T[popidx][1])).index("GQ")
                    DP_idx=(re.split(":", T[popidx][1])).index("DP")
                    for ind in T[popidx][2]:
                        if (len(re.split(":",ind))==1 or "./." in ind) or int(re.split(":", ind)[GQ_idx])<self.GQthreshold or int(re.split(":", ind)[DP_idx])<self.DPindthreshold:# ./.
                            site.append("NN")
                            continue
                        GT01 = re.split("/", re.split(":", ind)[GT_idx])
                        GT_TGT=T[int(GT01[0])+1]+T[int(GT01[1])+1]
#                         AC +=1
                        AN +=2
                        site.append(GT_TGT.upper())
#                     else:
#                         print("pass",self.vcfnamelist[popidx-3])
                elif self.MethodToSeqpoplist[popidx-3]=="pool":
                    print("unfinished")
        if AN>=self.minAN and AC>=self.minAC and AC<=(self.N*2-self.minAC):
            for x in range(self.N):
                self.seqs[x].append(site[x])
            self.positions.append(T[0])
    def getResult(self):
        pseudoPhasedSeqs=[]
        Nsites=len(self.positions)
        if Nsites<20:
            return 0,np.NaN
        for x in range(self.N):
            pseudoPhasedSeqs+= pseudoPhase(self.seqs[x], "pairs")
        if pseudoPhasedSeqs is not None:
            self.array = np.array([list(seq) for seq in pseudoPhasedSeqs])
            self.numArray = np.array([[numSeqDict[b] for b in seq] for seq in pseudoPhasedSeqs])
        else:
            self.array = np.empty((0,self.N))
            self.numArray = np.empty((0,self.N))
        self.nanMask = ~np.isnan(self.numArray)
        #get distMatrix
        distMat=distMatrix(self.numArray)
        np.fill_diagonal(distMat, np.NaN)
        dxy=np.nanmean(distMat[np.ix_(self.popIndices[0],self.popIndices[1])])
        print(len(pseudoPhasedSeqs[random.randint(0,self.N)]),Nsites,dxy)
        self.positions=[]
        self.seqs=[[] for e in range(self.N)]# pop1_inds pop2_inds
        self.numArray=[]
        self.array=[]
        return Nsites, dxy
class Caculate_popPI(Caculator):
    def __init__(self,considerINDEL,tvcfconfig,outputname):
        super().__init__()
        self.outputname=outputname
        self.considerINDEL=considerINDEL
        self.indnamesOfEachPop=[]
        self.MethodToSeqpoplist=[]
        self.listOfpopvcfRecsmapByAChr=[]
        for vcfconfigf in [tvcfconfig[0]]:
            vcfconfig=open(vcfconfigf,"r")
            self.listOfpopvcfRecsmapByAChr.append({})
            for line in vcfconfig:
                vcffilename_obj=re.search(r"vcffilename=(.*)",line.strip())
                if vcffilename_obj!=None:
                    vcfname=vcffilename_obj.group(1).strip()
                    self.vcfnamelist.append(vcfname)
                    self.outputname+=("_"+re.split(r"\.",re.search(r"[^/]*$",vcfname).group(0))[0])[:3]
                    self.vcfnameKEY_vcfobj_pyBAMfilesVALUE[vcfname]=[]
                    self.vcfnameKEY_vcfobj_pyBAMfilesVALUE[vcfname].append(VCFutil.VCF_Data(vcfname))
                    self.indnamesOfEachPop.append(self.vcfnameKEY_vcfobj_pyBAMfilesVALUE[vcfname][0].VcfIndexMap["title"][9:])
                elif line.split():
                    self.vcfnameKEY_vcfobj_pyBAMfilesVALUE[vcfname].append(pysam.Samfile(line.strip(),'rb'))
            vcfconfig.close()
            if re.search(r"indvd[^/]+",vcfname)!=None:
                self.MethodToSeqpoplist.append("indvd")
    
            elif re.search(r"pool[^/]+",vcfname)!=None:
                self.MethodToSeqpoplist.append("pool")
    
            else:
                print("vcfname must with 'pool' or 'indvd'")
                exit(-1)    
            
        self.minAN=40;self.minAC=2
        self.GQthreshold=30
        self.DPindthreshold=6
        print(self.indnamesOfEachPop)
        self.N=len(self.indnamesOfEachPop[0])#+len(self.indnamesOfEachPop[1])
        self.popIndices=[[]]
        for ind_idx in range(len(self.indnamesOfEachPop[0])):
            self.popIndices[0].append(2*ind_idx)
            self.popIndices[0].append(2*ind_idx+1)
#         for ind_idx in range(len(self.indnamesOfEachPop[0]),len(self.indnamesOfEachPop[0])+len(self.indnamesOfEachPop[1])):
#             self.popIndices[1].append(2*ind_idx)
#             self.popIndices[1].append(2*ind_idx+1)
        print(self.popIndices)
        #below variable should be changed every win
        self.positions=[]
        self.seqs=[[] for e in range(self.N)]# pop1_inds pop2_inds
        self.numArray=[]
        self.array=[]
    def process(self,T):
        """T is like (pos,REF,ALT,(INFO,FORMAT,sampleslist),(INFO,FORMAT,sampleslist)) should have only NOE VCF POP
        """
        if self.considerINDEL == "no" and (len(T[1]) != 1 or len(T[2]) != 1):
            return
        site=[]#should have 
        AN=AC=0
        for popidx in range(3,4):
            if T[popidx]==None:
                # check depth ,if passed treat as fix as ref
                if len(self.vcfnameKEY_vcfobj_pyBAMfilesVALUE[self.vcfnamelist[popidx-3]])==1:
#                     print("skip this pos",T)
                    return
                else:
#                     depth_linelist=self.depthobjlist[tpopidx-3].getdepthByPos_optimized(self.currentchrID,T[0])
                    sum_depth=0
                    for samfile in self.vcfnameKEY_vcfobj_pyBAMfilesVALUE[self.vcfnamelist[popidx-3]][1:]:
                        ACGTdep=samfile.count_coverage(self.currentchrID,T[0]-1,T[0])
                        for dep in ACGTdep:
                            sum_depth+=dep[0]

                    if sum_depth>self.DPindthreshold*len(self.indnamesOfEachPop[popidx-3])*1.2:#use 1.2 is ribitrarily, beacuse we don't have GQ to filter some ind so use more strength thread in depth
                        site+=[T[1].upper()*2 for x in range(len(self.indnamesOfEachPop[popidx-3]))]
                        AN+=(len(self.indnamesOfEachPop[popidx-3])*2)
                    else:
                        return
            else:
                if self.MethodToSeqpoplist[popidx-3]=="indvd":
#                     AN += int(re.search(r"AN=(\d+)[;,]", T[popidx][0]).group(1))
                    AC += int(re.search(r"AC=(\d+)[;,]", T[popidx][0]).group(1))
                    GT_idx = (re.split(":", T[popidx][1])).index("GT")
                    GQ_idx=(re.split(":", T[popidx][1])).index("GQ")
                    DP_idx=(re.split(":", T[popidx][1])).index("DP")
                    for ind in T[popidx][2]:
                        if (len(re.split(":",ind))==1 or "./." in ind) or int(re.split(":", ind)[GQ_idx])<self.GQthreshold or int(re.split(":", ind)[DP_idx])<self.DPindthreshold:# ./.
                            site.append("NN")
                            continue
                        GT01 = re.split("/", re.split(":", ind)[GT_idx])
                        GT_TGT=T[int(GT01[0])+1]+T[int(GT01[1])+1]
#                         AC +=1
                        AN +=2
                        site.append(GT_TGT.upper())
#                     else:
#                         print("pass",self.vcfnamelist[popidx-3])
                elif self.MethodToSeqpoplist[popidx-3]=="pool":
                    print("unfinished")
        if AN>=self.minAN and AC>=self.minAC and AC<=(self.N*2-self.minAC):
            for x in range(self.N):
                self.seqs[x].append(site[x])
            self.positions.append(T[0])
    def getResult(self):
        pseudoPhasedSeqs=[]
        Nsites=len(self.positions)
        if Nsites<20:
            return 0,np.NaN
        for x in range(self.N):
            pseudoPhasedSeqs+= pseudoPhase(self.seqs[x], "pairs")
        if pseudoPhasedSeqs is not None:
            self.array = np.array([list(seq) for seq in pseudoPhasedSeqs])
            self.numArray = np.array([[numSeqDict[b] for b in seq] for seq in pseudoPhasedSeqs])
        else:
            self.array = np.empty((0,self.N))
            self.numArray = np.empty((0,self.N))
        self.nanMask = ~np.isnan(self.numArray)
        #get distMatrix
        distMat=distMatrix(self.numArray)
        np.fill_diagonal(distMat, np.NaN)
        pi=np.nanmean(distMat[np.ix_(self.popIndices[0],self.popIndices[0])])
        print(len(pseudoPhasedSeqs[random.randint(0,self.N)]),Nsites,pi)
        self.positions=[]
        self.seqs=[[] for e in range(self.N)]# pop1_inds pop2_inds
        self.numArray=[]
        self.array=[]
        return Nsites, pi
class Caculate_SNPsPerBIN(Caculator):
    def __init__(self, winwidth, considerINDEL="no", MethodToSeq="pool"):
        self.considerINDEL = considerINDEL.lower()
        self.winwidth = winwidth
        self.COUNTED = 0
        self.MethodToSeq = MethodToSeq
    def process(self, T, seqerrorrate=0.01):
        self.COUNTED += 1
        return
        if self.considerINDEL == "no" and (len(T[1]) != 1 or len(T[2]) != 1):
            return
        if self.considerINDEL == "just" and (len(T[1]) == 1 and len(T[2]) == 1):
            return
        dp4 = re.search(r"DP4=(\d*),(\d*),(\d*),(\d*)", T[3])
        refdep = 0;altalleledep = 0
        if dp4 != None:  # vcf from samtools 
            refdep = int(dp4.group(1)) + int(dp4.group(2))
            altalleledep = int(dp4.group(3)) + int(dp4.group(4))    
        else:
            if self.MethodToSeq == "pool":
                AD_idx = (re.split(":", T[4])).index("AD")  # gatk GT:AD:DP:GQ:PL
                for sample in T[5]:
                    if len(re.split(":", sample)) == 1:  # ./.
                        continue
                    AD_depth = re.split(",", re.split(":", sample)[AD_idx])
                    try :
                        refdep += int(AD_depth[0])
                        altalleledep += int(AD_depth[1])
                    except ValueError:
                        print(sample, end="|")

                if refdep <= seqerrorrate * (refdep + altalleledep):  # not fixed
                    return
            elif self.MethodToSeq == "indvd":
                AN = int(re.search(r"AN=(\d+);", T[3]).group(1))
                AF = int(float(re.search(r"AF=([\d\.e-]+);", T[3]).group(1)))
                if AF == 1:
                    return
#                 refdep=AN-AC
#                 altalleledep=AC                
#        if refdep+altalleledep<10:
#            return
        self.COUNTED += 1
    def getResult(self):
        snpsinthiswin = self.COUNTED
        snpsdensity = 1000 * snpsinthiswin / self.winwidth
        self.COUNTED = 0
        return snpsinthiswin, snpsdensity
class Caculate_phastConsValue(Caculator):
    def __init__(self):
        super().__init__()
        self.conservationvalue = 0
        self.totalPostionsAwin = 0
    def process(self, T, NumOfPositions):
        self.conservationvalue += T[2] * NumOfPositions
        self.totalPostionsAwin += NumOfPositions
    def getResult(self):
        winvalue = "NA"
        if self.totalPostionsAwin == 0:
            self.conservationvalue = 0
            print("getResult")
            return "NA"
        else:
            winvalue = self.conservationvalue
            print(self.conservationvalue, self.totalPostionsAwin)
            winvalue = self.conservationvalue / self.totalPostionsAwin
            self.conservationvalue = 0
            self.totalPostionsAwin = 0
            return winvalue
def siteFreqs(numArray,nanMask, sites):# i don't quite understand this especially when sites is 0
    N,l=numArray.shape
    if not sites and sites!=0: sites = range(l)
    if type(sites) is not list: sites = [sites]
    return [binFreqs(numArray[:,x][nanMask[:,x]].astype(int)) for x in sites]
def binFreqs(numArr):
    n = len(numArr)
    if n == 0: return np.array([np.NaN]*4)
    else: return 1.* np.bincount(numArr, minlength=4) / n
class Caculate_ABB_BAB_BBAA(Caculator):
    def __init__(self,considerINDEL,p1vcfconfig,p2vcfconfig,p3vcfconfig,Ovcfconfig,outputname,cpid,pseudoPoolAN):
        super().__init__()
        self.outputname=outputname
        self.considerINDEL=considerINDEL
        self.pseudoAN=pseudoPoolAN# must be even
        self.indnameOfEachPop=[[],[],[],[]]
        self.MethodToSeqpoplist=[]
        self.listOfpopvcfRecsmapByAChr=[]
        self.vcflistOf4Pop=[[],[],[],[]]
        i=0
        for vcfconfigflist in [p1vcfconfig,p2vcfconfig,p3vcfconfig,Ovcfconfig]:# notice the order
            self.MethodToSeqpoplist.append([])
            for vcfconfigf in vcfconfigflist:
                vcfconfig=open(vcfconfigf,"r")
                for line in vcfconfig:
                    self.listOfpopvcfRecsmapByAChr.append({})
                    vcffilename_obj=re.search(r"vcffilename=(.*)", line.strip())
                    if vcffilename_obj!=None:
                        vcfname=vcffilename_obj.group(1).strip()
                        self.vcfnamelist.append(vcfname)
                        self.vcflistOf4Pop[i].append(vcfname)
                        self.outputname+=("_"+re.split(r"\.",re.search(r"[^/]*$",vcfname).group(0))[0])[:3]
                        self.vcfnameKEY_vcfobj_pyBAMfilesVALUE[vcfname]=[]
                        self.vcfnameKEY_vcfobj_pyBAMfilesVALUE[vcfname].append(VCFutil.VCF_Data(vcfname))
                        self.indnameOfEachPop[i]+=list(self.vcfnameKEY_vcfobj_pyBAMfilesVALUE[vcfname][0].VcfIndexMap["title"][9:])
                    elif line.split():
                        self.vcfnameKEY_vcfobj_pyBAMfilesVALUE[vcfname].append(pysam.Samfile(line.strip(),"rb"))
                vcfconfig.close()
                if re.search(r'indvd[^/]+',vcfname)!=None:
                    self.MethodToSeqpoplist[i].append("indvd")
                elif re.search(r"pool[^/]+",vcfname)!=None:
                    self.MethodToSeqpoplist[i].append("pool")
                else:
                    print("vcfname must with 'pool' or 'indvd'")
                    exit(-1)
            if "pool" in self.MethodToSeqpoplist[i]:
                self.indnameOfEachPop[i]=["P"+str(i+1)+"_pseudoindvd"+str(j) for j in range(1,int(self.pseudoAN/2)+1)]
            i+=1
        
        self.minAN=12;self.minAC=1;self.DPpoolthreshold=20;self.DPindthreshold=6;self.GQthreshold=30
        print(self.indnameOfEachPop)
        self.tempfile=open(self.outputname+"mallard_spotbilled_domestic_fanya"+str(cpid)+".snp",'w')
        print("chrNo\tANC\tDER\tP1derFreq\tP2derFreq\tP3derFreq\tP4derFreq\tBBAA\tABBA\tBABA",file=self.tempfile)
        self.positions=[]#all four population are the same
        self.popseqs=[[],[],[],[]];self.popnumArray=[[],[],[],[]];self.poparray=[[],[],[],[]]
        self.popseqs[0]=[[] for e in range(len(self.indnameOfEachPop[0]))]; self.popnumArray[0]=[]; self.poparray[0]=[]#P1
        self.popseqs[1]=[[] for e in range(len(self.indnameOfEachPop[1]))]; self.popnumArray[1]=[]; self.poparray[1]=[]#P2
        self.popseqs[2]=[[] for e in range(len(self.indnameOfEachPop[2]))]; self.popnumArray[2]=[]; self.poparray[2]=[]#P3
        self.popseqs[3]=[[] for e in range(len(self.indnameOfEachPop[3]))]; self.popnumArray[3]=[]; self.poparray[3]=[]#P4
        self.all4popseqs=[[] for e in range(len(self.indnameOfEachPop[0]+self.indnameOfEachPop[1]+self.indnameOfEachPop[2]+self.indnameOfEachPop[3]))];self.all4popnumArray=[];self.all4poparray=[]
    def process(self,T):
        """T is like (pos,REF,ALT,(INFO,FORMAT,sampleslist),(INFO,FORMAT,sampleslist),(),(),) should have more than four populations
        """       
        if self.considerINDEL=="no" and (len(T[1])!=1 or len(T[2])!=1):
            return
        startvcfidx=3
        all4popidx=0
        siteP4=[]
        for Pidx in range(4):
            site=[]
            AN=AC=0#N will not count into AN,so the AN will be used as a threshold,like DPpool 
            AFpool=[];DPpool=[]#event indvd vcf compain with DPpool,because sometimes the indvd will be regard as pool when there is pool in the same group
            for popidx in range(startvcfidx,startvcfidx+len(self.vcflistOf4Pop[Pidx])):
                if T[popidx]==None:
                    if len(self.vcfnameKEY_vcfobj_pyBAMfilesVALUE[self.vcflistOf4Pop[Pidx][popidx-startvcfidx]])==1:
                        print("no sam is avalible,treat as fixed")
                        DPpool.append(self.DPpoolthreshold);AFpool.append(0)
                        continue
                    else:
                        sum_depth=0
                        for samfile in self.vcfnameKEY_vcfobj_pyBAMfilesVALUE[self.vcflistOf4Pop[Pidx][popidx-startvcfidx]][1:]:
                            ACGTdep=samfile.count_coverage(self.currentchrID,T[0]-1,T[0])
                            for dep in ACGTdep:
                                sum_depth+=dep[0]
                        if (self.MethodToSeqpoplist[Pidx][popidx-startvcfidx]=="pool" and sum_depth>self.DPpoolthreshold):
                            DPpool.append(self.DPpoolthreshold);AFpool.append(0)
                        elif (sum_depth>self.DPindthreshold*len(self.indnameOfEachPop[Pidx][popidx-startvcfidx])*1.2 and self.MethodToSeqpoplist[Pidx][popidx-startvcfidx]=="indvd"):
                            DPpool.append(self.DPpoolthreshold);AFpool.append(0)
#                             print("produce site 1",site,"\n",self.vcfnameKEY_vcfobj_pyBAMfilesVALUE[self.vcfnamelist[popidx-3]][0].VcfIndexMap["title"][9:])
                            site+=[T[1].upper()*2 for x in range(len(self.vcfnameKEY_vcfobj_pyBAMfilesVALUE[self.vcfnamelist[popidx-3]][0].VcfIndexMap["title"][9:]))]
                            AN+=(len(self.indnameOfEachPop[Pidx][popidx-startvcfidx])*2)
                        else:
                            continue#AN will be zero,this time site is 0,no AFpool no DPpool,
                else:
                    if self.MethodToSeqpoplist[Pidx][popidx-startvcfidx]=="indvd":
                        AFpool.append(int(float(re.search(r"AF=([\d\.e-]+);", T[popidx][0]).group(1))))
                        DPpool.append(int(re.search(r"AN=(\d+)[;,]", T[popidx][0]).group(1)))#just take a position
                        AC += int(re.search(r"AC=(\d+)[;,]", T[popidx][0]).group(1))
                        GT_idx = (re.split(":", T[popidx][1])).index("GT")
                        GQ_idx=(re.split(":", T[popidx][1])).index("GQ")
                        DP_idx=(re.split(":", T[popidx][1])).index("DP")
                        for ind in T[popidx][2]:
                            if (len(re.split(":",ind))==1 or "./." in ind) or int(re.split(":", ind)[GQ_idx])<self.GQthreshold or int(re.split(":", ind)[DP_idx])<self.DPindthreshold:# ./.
                                site.append("NN")
                                continue
                            GT01 = re.split("/", re.split(":", ind)[GT_idx])
                            GT_TGT=T[int(GT01[0])+1]+T[int(GT01[1])+1]
    #                         AC +=1
                            AN +=2
                            site.append(GT_TGT.upper())
#                         print("in dvd",self.vcflistOf4Pop[Pidx][popidx-startvcfidx],"should equal",self.vcfnamelist[popidx-3],site)
                    elif self.MethodToSeqpoplist[Pidx][popidx-startvcfidx]=="pool":
                        AD_idx = (re.split(":", T[popidx][1])).index("AD")
                        refdep = 0;altalleledep = 0
                        for indsample in T[popidx][2]:
                            if len(re.split(":",indsample))==1:#./.
                                continue
                            AD_depth=re.split(",", re.split(":", indsample)[AD_idx])
                            refdep += int(AD_depth[0])
                            altalleledep += int(AD_depth[1])
                        DP=altalleledep+refdep
                        DPpool.append(DP)
                        AF=altalleledep/DP
                        AFpool.append(AF)
                        
            if "pool" in self.MethodToSeqpoplist[Pidx]:#
                for afidx in reversed(range(len(AFpool))):
                    if DPpool[afidx]<self.DPpoolthreshold:
                        AFpool.pop(afidx)
                if len(AFpool)<len(self.vcflistOf4Pop[Pidx])/4:
#                     print("return")
                    return
                AC=pseudoAC=round(np.mean(AFpool)*self.pseudoAN)
                pseudoAC=int(pseudoAC)
                sitestr=T[2]*(pseudoAC)+T[1]*(self.pseudoAN-pseudoAC)
                AN=len(sitestr)
                site=["".join(sitestr[s:s+2]) for s in range(0,len(sitestr),2)]
            if AN<self.minAN:#site==[] only when AN==0
#                 print("return")
                return
            siteP4+=site
            startvcfidx+=len(self.vcflistOf4Pop[Pidx])
#             print(site,len(self.popseqs[Pidx]),Pidx)
#             if Pidx==3 or AN>=self.minAN and AC>self.minAC and AC<=(len(self.popseqs[Pidx])*2-self.minAC):
        bases=set("".join(siteP4))
        if "N" in bases: bases.remove("N")
        if len(siteP4)==len(self.all4popseqs) and len(bases)==2:
            for Pidx in range(4):
                for x in range(len(self.popseqs[Pidx])):
#                     print(all4popidx,len(self.all4popseqs),x,self.indnameOfEachPop[Pidx][x],end="|")
                    self.popseqs[Pidx][x].append(siteP4[all4popidx])#self.pop1seqs or self.pop2seqs ....
                    self.all4popseqs[all4popidx].append(siteP4[all4popidx])
                    all4popidx+=1
#                 print()
            self.positions.append(T[0])
    def getResult(self):
        Nsites=len(self.positions)
        if Nsites==0:
            print(Nsites)
            return 0,[np.NaN,np.NaN] #ABBAsum,BABAsum,BBAAsum,D,fd
        #initialize array from self.xxxseqs for popall ,pop1, pop2, pop3, pop4
        pseudoPhasedSeqs=[]
        for x in range(len(self.all4popseqs)):
#             print("self.all4popseqs[]",x,self.all4popseqs[x])
            pseudoPhasedSeqs+=pseudoPhase(self.all4popseqs[x],"pairs")
#         print(pseudoPhasedSeqs)
        if pseudoPhasedSeqs is not None:
            self.all4poparray=np.array([list(seq) for seq in pseudoPhasedSeqs])
            self.all4popnumArray=np.array([[numSeqDict[b] for b in seq] for seq in pseudoPhasedSeqs])
        else:
            self.all4poparray=np.empty((0,len(self.all4popseqs)))
            self.all4popnumArray=np.empty((0,len(self.all4popseqs)))
#         print(self.all4popnumArray)
        nanMask_all4pop = ~np.isnan(self.all4popnumArray)
        nanMasklist=[[],[],[],[]]
        for Pidx in range(4):
            pseudoPhasedSeqs=[]
            for x in range(len(self.popseqs[Pidx])):
                pseudoPhasedSeqs+=pseudoPhase(self.popseqs[Pidx][x], "pairs")
            if pseudoPhasedSeqs is not None:
                self.poparray[Pidx]=np.array([list(seq) for seq in pseudoPhasedSeqs])
                self.popnumArray[Pidx]=np.array([[numSeqDict[b] for b in seq] for seq in pseudoPhasedSeqs])
            else:
                self.poparray[Pidx]=np.empty((0,len(self.popseqs[Pidx])))
                self.popnumArray[Pidx]=np.empty((0,len(self.popseqs[Pidx])))
            nanMasklist[Pidx] = ~np.isnan(self.popnumArray[Pidx])
        #initialize array finished 
        BBAAsum=ABBAsum= BABAsum = maxABBAsum = maxBABAsum = 0.0
        
        for i in range(Nsites):
            allFreqs=siteFreqs(self.all4popnumArray, nanMask_all4pop, i)[0]
            P1Freqs,P2Freqs,P3Freqs,P4Freqs = [siteFreqs(refertonumArray,refertonanMask,i)[0] for refertonumArray,refertonanMask in ((self.popnumArray[0],nanMasklist[0]), (self.popnumArray[1],nanMasklist[1]), (self.popnumArray[2],nanMasklist[2]), (self.popnumArray[3],nanMasklist[3]))]    
            if np.any(np.isnan(P1Freqs)) or np.any(np.isnan(P2Freqs)) or np.any(np.isnan(P3Freqs)) or np.any(np.isnan(P4Freqs)): 
                print(self.currentchrID,self.positions[i],"N","N",P1Freqs[np.argsort(allFreqs)[-1]],P2Freqs[np.argsort(allFreqs)[-1]],P2Freqs[np.argsort(allFreqs)[-1]],P4Freqs[np.argsort(allFreqs)[-1]],0,0,0,sep="\t",file=self.tempfile)
                continue
                
            #if the outgroup is fixed, then that is the ancestral state - otherwise the ancestral state is the most common allele overall
#             print(P1Freqs,P2Freqs,P3Freqs,P4Freqs,allFreqs,file=self.tempfile)
#             print(allFreqs,np.where(allFreqs > 0))
            if np.max(P4Freqs) == 1.:
                anc = np.where(P4Freqs == 1)[0][0]#ancetral allele is which is fixed (get the index)
                der = [i for i in np.where(allFreqs > 0)[0] if i != anc][0]  # derived is the index that is > 0 but not anc
            else:der = np.argsort(allFreqs)[-2];anc= np.argsort(allFreqs)[-1]# the less common base overall
            #derived allele frequencies
            P1derFreq = P1Freqs[der]
            P2derFreq = P2Freqs[der]
            P3derFreq = P3Freqs[der]
            P4derFreq = P4Freqs[der]
            PDderFreq = max(P2derFreq,P3derFreq)
            BBAA=P1derFreq*P2derFreq*(1 - P3derFreq)*(1 - P4derFreq)
            ABBA=(1 - P1derFreq) * P2derFreq * P3derFreq * (1 - P4derFreq)
            BABA=P1derFreq * (1 - P2derFreq) * P3derFreq * (1 - P4derFreq)
            BBAAsum+=BBAA
            ABBAsum += ABBA
            BABAsum += BABA
            maxABBAsum += (1 - P1derFreq) * PDderFreq * PDderFreq * (1 - P4derFreq)
            maxBABAsum += P1derFreq * (1 - PDderFreq) * PDderFreq * (1 - P4derFreq)
            try:
                D=(ABBAsum - BABAsum) / (ABBAsum + BABAsum)
                print(self.currentchrID,self.positions[i],SeqnumDict[anc],SeqnumDict[der],P1derFreq,P2derFreq,P3derFreq,P4derFreq,BBAA,ABBA,BABA,sep="\t",file=self.tempfile)
            except:D=np.NaN
            try: 
                if D>=0 :fd=(ABBAsum - BABAsum) / (maxABBAsum - maxBABAsum)
                else: fd=np.NaN
            except: fd=np.NaN
            
        for Pidx in range(4):
            self.popseqs[Pidx]=[[] for e in range(len(self.indnameOfEachPop[Pidx]))]
            self.popnumArray[Pidx]=[]
            self.poparray[Pidx]=[]
        self.all4popseqs=[[] for e in range(len(self.indnameOfEachPop[0]+self.indnameOfEachPop[1]+self.indnameOfEachPop[2]+self.indnameOfEachPop[3]))];self.all4popnumArray=[];self.all4poparray=[]
        self.positions=[]
        if Nsites<8:
            print(Nsites)
            return 0,[np.NaN,np.NaN] #ABBAsum,BABAsum,BBAAsum,D,fd
        return Nsites,[D,fd]
         
         
class Caculate_Dstatistics(Caculator):
    def __init__(self, considerFixed=False):
        super().__init__()
        self.ABBA = 0
        self.BABA = 0
        self.numerator_fixed = 0
        self.denominator_fixed = 0
        self.numerator_snp = 0
        self.denominator_snp = 0
        self.considerFixed = considerFixed
        self.COUNTEDforSNP_notonlyfixed = 0
    def process(self, T, seqerrorrate=0.01):
        """T:(pos,"a,b","c,d","e,f",A_base_idx)     1 - A_base_idx= B_base_idx ie. T[4] is the A idx . 1-T[4] is the B idx
        """
        p1A = int(re.split(r",", T[1])[T[4]]);p1B = int(re.split(r",", T[1])[1 - T[4]])
        p2A = int(re.split(r",", T[2])[T[4]]);p2B = int(re.split(r",", T[2])[1 - T[4]])
        p3A = int(re.split(r",", T[3])[T[4]]);p3B = int(re.split(r",", T[3])[1 - T[4]])
        if p3A == 0 and p3B != 0:  # p3 fixed as B
            if p2A == 0 and p2B != 0:  # p2 fixed as B
                if p1B == 0 and p1A != 0:  # p1 fixed as A
                    self.ABBA += 1
                    print(T, "abba")
                    self.numerator_fixed += 1
                    self.denominator_fixed += 1
            elif p2B == 0 and p2A != 0:  # p2 fixed as A
                if p1A == 0 and p1B != 0:  # p1 fixed as B
                    self.BABA += 1
                    print(T, "baba")
                    self.numerator_fixed += -1
                    self.denominator_fixed += 1
        if (not self.considerFixed) and(p1A == 0 or p1B == 0 or p2A == 0 or p2B == 0 or p3A == 0 or p3B == 0):
            return
        try:
            self.numerator_snp += p3B / (p3B + p3A) * ((p1A / (p1A + p1B)) * (p2B / (p2A + p2B)) - (p1B / (p1A + p1B)) * (p2A / (p2A + p2B)))
            self.denominator_snp += p3B / (p3B + p3A) * ((p1A / (p1A + p1B)) * (p2B / (p2A + p2B)) + (p1B / (p1A + p1B)) * (p2A / (p2A + p2B)))
            self.COUNTEDforSNP_notonlyfixed += 1
        except ZeroDivisionError:
            print(self.denominator_snp, self.numerator_snp, T)
    def getResult(self):
        ABBAcount = self.ABBA
        BABAcount = self.BABA
        try:
            D_fixed = self.numerator_fixed / self.denominator_fixed
            D_snp = self.numerator_snp / self.denominator_snp
        except ZeroDivisionError:
            D_fixed = 'NA'
            D_snp = 'NA'
        self.numerator_fixed = 0;self.denominator_fixed = 0;self.ABBA = 0;self.BABA = 0
        self.numerator_snp = 0;self.denominator_snp = 0;noofsnps = copy.deepcopy(self.COUNTEDforSNP_notonlyfixed);self.COUNTEDforSNP_notonlyfixed = 0
        return ABBAcount, BABAcount, D_fixed, D_snp, noofsnps
        
class Caculate_df(Caculator):
    def __init__(self,pop1idxlist,pop2idxlist,filehanlder):
        super().__init__()
        self.pop1idxlist=pop1idxlist
        self.pop2idxlist=pop2idxlist
        self.COUNTED = 0# nooffixdiff
        self.COUNTEDadditional=[0,[0,0]]#noofheterozygosity (pop1recs,pop2recs)
        self.unsufficentfixediff=0
        self.pop1_indvds=None#=6#when pop1 is none at a pos,and no depth information
        self.pop2_indvds=None#=6
        self.filepoiner=filehanlder
        self.currentchrID=None
    def process(self,T):
        #no matter the present of the multiple allele or not,it's still work correct.
        pop1reffixed=0;pop1altfixed=0;pop2reffixed=0;pop2altfixed=0
        pop1het=0;pop2het=0
        pop1recs=0;pop2recs=0
        GT_idx = (re.split(":", T[3 + 0][1])).index("GT")  # gatk GT:AD:DP:GQ:PL
        for pop1idx in self.pop1idxlist:
            sample=T[3+0][2][pop1idx]
            if len(re.split(":", sample)) == 1:  # ./.
                continue
            pop1recs+=1
            if re.split(":", sample)[GT_idx]=="0/0":
                pop1reffixed+=1
            elif re.split(":", sample)[GT_idx]=="1/1":
                pop1altfixed+=1
            elif re.split(":", sample)[GT_idx]=="0/1":
                pop1het+=1
        for pop2idx in self.pop2idxlist:
            sample=T[3+0][2][pop2idx]
            if len(re.split(":", sample)) == 1:  # ./.
                continue
            pop2recs+=1
            if re.split(":", sample)[GT_idx]=="0/0":
                pop2reffixed+=1
            elif re.split(":", sample)[GT_idx]=="1/1":
                pop2altfixed+=1
            elif re.split(":", sample)[GT_idx]=="0/1":
                pop2het+=1
#                                  pop1 fixed as alt ,pop2 fixed as ref                                                      pop1 fixed as ref ,pop2 fixed as alt
#         print(self.pop2_indvds,pop2reffixed+pop2altfixed,pop2reffixed*pop2altfixed,self.pop1_indvds,pop1altfixed+pop1reffixed,pop1reffixed*pop1altfixed)
#         print(pop1reffixed,"==0",pop1altfixed,"==",self.pop1_indvds,pop2reffixed,"==",self.pop2_indvds,pop2altfixed,"==0",(pop1reffixed==0 and pop1altfixed==self.pop1_indvds and pop2reffixed==self.pop2_indvds and pop2altfixed==0 ),pop1reffixed,"==",self.pop1_indvds,pop1altfixed,"==0",pop2reffixed,"==0",pop2altfixed,"==",self.pop2_indvds,(pop1reffixed==self.pop1_indvds and pop1altfixed==0 and pop2reffixed==0 and pop2altfixed==self.pop2_indvds))
        if (pop1reffixed==0 and pop1altfixed==self.pop1_indvds and pop2reffixed==self.pop2_indvds and pop2altfixed==0 ) or (pop1reffixed==self.pop1_indvds and pop1altfixed==0 and pop2reffixed==0 and pop2altfixed==self.pop2_indvds):
            self.COUNTED+=1
            print(self.currentchrID,T,file=self.filepoiner)
#             print("bingle",file=open("tttttt.txt",'a'))
        elif (pop1reffixed*pop1altfixed!=0 or pop1het!=0) or (pop2reffixed*pop2altfixed!=0 or pop2het!=0 ) :
            self.COUNTEDadditional[0]+=1
        elif (pop1reffixed*pop1altfixed==0 and pop1het==0) and (pop2reffixed*pop2altfixed==0 and pop2het==0 ) :
            self.COUNTEDadditional[1][0]+=pop1recs/self.pop1_indvds;self.COUNTEDadditional[1][1]+=pop2recs/self.pop2_indvds;self.unsufficentfixediff+=1
    def getResult(self):
        if self.unsufficentfixediff!=0:
            pop1unsufficentfixed=self.COUNTEDadditional[1][0]/self.unsufficentfixediff
            pop2unsufficentfixed=self.COUNTEDadditional[1][1]/self.unsufficentfixediff
        else:
            pop1unsufficentfixed=0;pop2unsufficentfixed=0
        noofhet=self.COUNTEDadditional[0]
        nooffixediff=self.COUNTED
        self.COUNTED=0
        self.COUNTEDadditional=[0,[0,0]]
        self.unsufficentfixediff=0
        return [noofhet,(pop1unsufficentfixed,pop2unsufficentfixed)],nooffixediff #self.COUNTEDadditional,self.COUNTED

class CaculatorToFindTAGs(Caculator):
    def __init__(self,mod,Interferingf):
        self.TAGSforAwin=[]
        self.cwinNo=0
        self.rand2snpForEveryWin=[]
        self.SNPwithAFforAwin=[]#element (pos,af,upstreamdistance,downstreamdistance)
        self.mod=mod#randomvcf,selectTAG
        self.previousPos=1
        self.Interferingf=Interferingf
        self.curchrom=""
    def process(self,T):
        "T=(pos, REF, ALT, INFO,FORMAT,sampleslist)"
        if self.mod=="randomvcf" and re.search(r"AF=([\d\.e-]+)[;,]",T[3])!=None:
            af=float(re.search(r"AF=([\d\.e-]+)[;,]",T[3]).group(1))
            if len(self.SNPwithAFforAwin)>=1: self.SNPwithAFforAwin[-1][3]=T[0]-self.previousPos 
            self.SNPwithAFforAwin.append([T[0],af,T[0]-self.previousPos,35,T[1],T[2]])
            self.previousPos=T[0]
        elif self.mod=="selectTAG":
            self.TAGSforAwin.append(T[0])
    def getResult(self):
        if self.mod=="randomvcf":
            passedsites=[e6 for e6 in self.SNPwithAFforAwin if (e6[2]>=35 or e6[3]>=35)]
            passedsites.sort(key=lambda x:abs(x[2]-0.5))
            if len(passedsites)>=2:
                snp1=passedsites[0]
                if snp1[2]<35 :
                    idx1=self.SNPwithAFforAwin.index(snp1)
                    for e6 in reversed(self.SNPwithAFforAwin[:idx1]):
                        if e6[0]>snp1[0]-35 and e6[0]!=snp1[0]: print(self.curchrom,*e6,file=self.Interferingf)
                elif snp1[3]<35:
                    idx1=self.SNPwithAFforAwin.index(snp1)
                    for e6 in self.SNPwithAFforAwin[idx1:]:
                        if e6[0]<snp1[0]+35 and e6[0]!=snp1[0]: print(self.curchrom,*e6,file=self.Interferingf)
                snp2=passedsites[1]
                if snp2[2]<35 :
                    idx2=self.SNPwithAFforAwin.index(snp2)
                    for e6 in reversed(self.SNPwithAFforAwin[:idx2]):
                        if e6[0]>snp2[0]-35 and e6[0]!=snp2[0]: print(self.curchrom,*e6,file=self.Interferingf)
                elif snp2[3]<35:
                    idx2=self.SNPwithAFforAwin.index(snp2)
                    for e6 in self.SNPwithAFforAwin[idx2:]:
                        if e6[0]<snp2[0]+35 and e6[0]!=snp2[0]: print(self.curchrom,*e6,file=self.Interferingf)                
                self.rand2snpForEveryWin.append((snp1,snp2))
            elif len(passedsites)==1:
                self.rand2snpForEveryWin.append((passedsites[0],"NA"))
                if passedsites[0][2]<35 :
                    idx1=self.SNPwithAFforAwin.index(passedsites[0])
                    for e6 in reversed(self.SNPwithAFforAwin[:idx1]):
                        if e6[0]>passedsites[0][0]-35 and e6[0]!=passedsites[0][0]: print(self.curchrom,*e6,file=self.Interferingf)
                elif passedsites[0][3]<35:
                    idx1=self.SNPwithAFforAwin.index(passedsites[0])
                    for e6 in self.SNPwithAFforAwin[idx1:]:
                        if e6[0]<passedsites[0][0]+35 and e6[0]!=passedsites[0][0]: print(self.curchrom,*e6,file=self.Interferingf)                
            else:
                self.rand2snpForEveryWin.append(("NA","NA"))
            self.SNPwithAFforAwin=[]
            return -1,self.rand2snpForEveryWin[-1]#-1,(snp1_6info,snp2_6info)
        elif self.mod=="selectTAG":
            self.TAGSforAwin=list(set(self.TAGSforAwin))
            if len(self.TAGSforAwin)>=2:
                idxlist=random.sample([j for j in range(len(self.TAGSforAwin))],2)
                tag1=self.TAGSforAwin[idxlist[0]]
                tag2=self.TAGSforAwin[idxlist[1]]
                nooftags=2
                valuetoReturn=[tag1,tag2]
            elif len(self.TAGSforAwin)==1:
                nooftags=1
                valuetoReturn=[self.TAGSforAwin[0],self.rand2snpForEveryWin[self.cwinNo][0]]
            else:
                nooftags=0
                self.rand2snpForEveryWin[self.cwinNo]
                valuetoReturn=[self.rand2snpForEveryWin[self.cwinNo][0],self.rand2snpForEveryWin[self.cwinNo][1]]
            self.TAGSforAwin=[]
            self.cwinNo+=1
            return nooftags,valuetoReturn
class Caculate_Hp_master_slave(Caculator):
    def __init__(self, listOftargetpopvcfconfig,outfileprewithpath, minsnps=10,depth=10):
        super().__init__()
        self.minsnps = minsnps
        self.depth=depth
        self.SeqMethodlist=[]
        self.vcfnamelist=[]
        self.vcfnameKEY_vcfobj_pyBAMfilesVALUE={}
        self.outputname=outfileprewithpath
        self.listOfpopvcfRecsmapByAChr=[]
        for vcfconfigfilename in listOftargetpopvcfconfig[:]:
            self.listOfpopvcfRecsmapByAChr.append({})
            vcfconfig=open(vcfconfigfilename,"r")
            for line in vcfconfig:
                vcffilename_obj=re.search(r"vcffilename=(.*)",line.strip())
                if vcffilename_obj!=None:
                    vcfname=vcffilename_obj.group(1).strip()
                    self.vcfnamelist.append(vcfname)
                    self.outputname+=("_"+re.split(r"\.",re.search(r"[^/]*$",vcfname).group(0))[0])[:3]
                    self.vcfnameKEY_vcfobj_pyBAMfilesVALUE[vcfname]=[]
                    self.vcfnameKEY_vcfobj_pyBAMfilesVALUE[vcfname].append(VCFutil.VCF_Data(vcfname))
                elif line.split():
                    self.vcfnameKEY_vcfobj_pyBAMfilesVALUE[vcfname].append(pysam.Samfile(line.strip(),'rb'))
            vcfconfig.close()
            if re.search(r"indvd[^/]+",vcfname)!=None:
                self.SeqMethodlist.append("indvd")
    
            elif re.search(r"pool[^/]+",vcfname)!=None:
                self.SeqMethodlist.append("pool")
    
            else:
                print("vcfname must with 'pool' or 'indvd'")
                exit(-1) 
        print(self.SeqMethodlist)
        self.COUNTED = [0] * len(self.SeqMethodlist)
        self.CNMI = [0] * len(self.SeqMethodlist)
        self.CNMA = [0] * len(self.SeqMethodlist)
        self.sum_mean_2pq = 0  
    def process(self, T, seqerrorrate=0.008, mode=1):
        if len(T[1]) != len(T[2]) or len(T[2])!=1 or len(T[2])!=1:
            return
        for MethodToSeq_idx in range(len(self.SeqMethodlist)):
            MethofToSeq = self.SeqMethodlist[MethodToSeq_idx]
            if T[3 + MethodToSeq_idx] == None:
                continue
            if MethofToSeq == "pool":
                refdep = 0;altalleledep = 0
                AD_idx = (re.split(":", T[3 + MethodToSeq_idx][1])).index("AD")  # gatk GT:AD:DP:GQ:PL
                for sample in T[3 + MethodToSeq_idx][2]:
                    if len(re.split(":", sample)) == 1:  # ./.
                        continue
                    AD_depth = re.split(",", re.split(":", sample)[AD_idx])
                    try :
                        refdep += int(AD_depth[0])*0.7
                        altalleledep += int(AD_depth[1])*0.7
                    except ValueError:
                        print(sample, end="|")
            elif MethofToSeq == "indvd":
                AF = float(re.search(r"AF=([\d\.e-]+)[;,]", T[3 + MethodToSeq_idx][0]).group(1))
                AN = int(re.search(r"AN=(\d+)[;,]", T[3 + MethodToSeq_idx][0]).group(1))
                AC = int(re.search(r"AC=(\d+)[;,]", T[3 + MethodToSeq_idx][0]).group(1))
                refdep = AN - AC
                altalleledep = AC
            if refdep <= seqerrorrate * (refdep + altalleledep):  # skip fixed as altallele ,ie refdep == 0
                continue
            if refdep + altalleledep < self.depth:
                continue
            self.COUNTED[MethodToSeq_idx] += 1
            if refdep < altalleledep:
                self.CNMI[MethodToSeq_idx] += refdep
                self.CNMA[MethodToSeq_idx] += altalleledep
            else:
                self.CNMA[MethodToSeq_idx] += refdep
                self.CNMI[MethodToSeq_idx] += altalleledep
    def getResult(self):
        HETEROZY = ['NA'] * len(self.SeqMethodlist)
        for MethodToSeq_idx in range(len(self.SeqMethodlist)):
            try:
                HETEROZY[MethodToSeq_idx] = self.CNMA[MethodToSeq_idx] * self.CNMI[MethodToSeq_idx] * 2 / ((self.CNMA[MethodToSeq_idx] + self.CNMI[MethodToSeq_idx]) ** 2)
            except ZeroDivisionError:
                # print("the Heterozigosity value of currentwindow is dividsion by zero,so set it to be NA")
                HETEROZY[MethodToSeq_idx] = 'NA'
        het_count = 0;het_sum = 0;pop_idx=0
        for pop_idx in range(len(HETEROZY)-1,-1,-1) :
            if HETEROZY[pop_idx] != 'NA' and self.COUNTED[pop_idx]>=self.minsnps:
                het_count += 1
                het_sum += HETEROZY[pop_idx]
            else:
                self.COUNTED.pop(pop_idx)
#                 print(HETEROZY[pop_idx],end="\t")
#         print()
        if self.COUNTED==[]:
            noofsnpcount=0
        else:
            noofsnpcount= min(self.COUNTED)
        if het_count == 0 :
            HETEROZY_toreturn = 'NA'
        else:
            HETEROZY_toreturn = het_sum / het_count
        self.COUNTED = [0] * len(self.SeqMethodlist)
        self.CNMA = [0] * len(self.SeqMethodlist)
        self.CNMI = [0] * len(self.SeqMethodlist)
        return noofsnpcount, HETEROZY_toreturn
class Caculate_Hp(Caculator):
    def __init__(self, SeqMethodlist=["pool"], minsnps=10,depth=10):
        super().__init__()
        self.minsnps = minsnps
        self.depth=depth
        self.COUNTED = [0] * len(SeqMethodlist)
        self.CNMI = [0] * len(SeqMethodlist)
        self.CNMA = [0] * len(SeqMethodlist)
        self.sum_mean_2pq = 0
        self.SeqMethodlist = SeqMethodlist
    def process(self, T, seqerrorrate=0.008, mode=1):
        if len(T[1]) != len(T[2]) or len(T[2])!=1 or len(T[2])!=1:
            return
        for MethodToSeq_idx in range(len(self.SeqMethodlist)):
            MethofToSeq = self.SeqMethodlist[MethodToSeq_idx]
            if T[3 + MethodToSeq_idx] == None:
                continue
            if MethofToSeq == "pool":
                refdep = 0;altalleledep = 0
                AD_idx = (re.split(":", T[3 + MethodToSeq_idx][1])).index("AD")  # gatk GT:AD:DP:GQ:PL
                for sample in T[3 + MethodToSeq_idx][2]:
                    if len(re.split(":", sample)) == 1:  # ./.
                        continue
                    AD_depth = re.split(",", re.split(":", sample)[AD_idx])
                    try :
                        refdep += int(AD_depth[0])*0.7
                        altalleledep += int(AD_depth[1])*0.7
                    except ValueError:
                        print(sample, end="|")
            elif MethofToSeq == "indvd":
                AF = float(re.search(r"AF=([\d\.e-]+)[;,]", T[3 + MethodToSeq_idx][0]).group(1))
                AN = int(re.search(r"AN=(\d+)[;,]", T[3 + MethodToSeq_idx][0]).group(1))
                AC = int(re.search(r"AC=(\d+)[;,]", T[3 + MethodToSeq_idx][0]).group(1))
                refdep = AN - AC
                altalleledep = AC
            if refdep <= seqerrorrate * (refdep + altalleledep):  # skip fixed as altallele ,ie refdep == 0
                continue
            if refdep + altalleledep < self.depth:
                continue
            self.COUNTED[MethodToSeq_idx] += 1
            if refdep < altalleledep:
                self.CNMI[MethodToSeq_idx] += refdep
                self.CNMA[MethodToSeq_idx] += altalleledep
            else:
                self.CNMA[MethodToSeq_idx] += refdep
                self.CNMI[MethodToSeq_idx] += altalleledep
    def getResult(self):
        HETEROZY = ['NA'] * len(self.SeqMethodlist)
        for MethodToSeq_idx in range(len(self.SeqMethodlist)):
            try:
                HETEROZY[MethodToSeq_idx] = self.CNMA[MethodToSeq_idx] * self.CNMI[MethodToSeq_idx] * 2 / ((self.CNMA[MethodToSeq_idx] + self.CNMI[MethodToSeq_idx]) ** 2)
            except ZeroDivisionError:
                # print("the Heterozigosity value of currentwindow is dividsion by zero,so set it to be NA")
                HETEROZY[MethodToSeq_idx] = 'NA'
        het_count = 0;het_sum = 0;pop_idx=0
        for pop_idx in range(len(HETEROZY)-1,-1,-1) :
            if HETEROZY[pop_idx] != 'NA' and self.COUNTED[pop_idx]>=self.minsnps:
                het_count += 1
                het_sum += HETEROZY[pop_idx]
            else:
                self.COUNTED.pop(pop_idx)
#                 print(HETEROZY[pop_idx],end="\t")
#         print()
    
        noofsnpcount = min(self.COUNTED)
        if het_count == 0 :
            HETEROZY_toreturn = 'NA'
        else:
            HETEROZY_toreturn = het_sum / het_count
        self.COUNTED = [0] * len(self.SeqMethodlist)
        self.CNMA = [0] * len(self.SeqMethodlist)
        self.CNMI = [0] * len(self.SeqMethodlist)
        return noofsnpcount, HETEROZY_toreturn
class Caculate_depth_judge(Caculator):
    def __init__(self, total_samples, winsize, mindepth, speciesorder=[], sampleidxlisttocount={}):
        self.mindepth = int(mindepth)
        self.total_samples = total_samples
        self.winsize = winsize
        if speciesorder == [] and (not sampleidxlisttocount):
            self.COVERED_COUNT = [0] * total_samples
            self.AVERAGE_DEPTH = [0] * total_samples
            self.speciesorder = None;
            self.sampleidxlisttocount = None
        elif len(speciesorder) != 0 and len(sampleidxlisttocount.keys()) != 0:
            self.COVERED_COUNT = [0] * len(speciesorder)
            self.AVERAGE_DEPTH = [0] * len(speciesorder)
            self.speciesorder = speciesorder
            self.sampleidxlisttocount = sampleidxlisttocount
    def process(self, T, seqerrorrate=0.01):
        """
        T=(pos,sample1dp,sample2dp,,,,,,)
        """
#         print(T,"\n",self.AVERAGE_DEPTH)
        if self.speciesorder == [] and (not self.sampleidxlisttocount):
            for sampleidx in range(1, len(T)):
                self.AVERAGE_DEPTH[sampleidx - 1] += int(T[sampleidx])
                if int(T[sampleidx]) >= self.mindepth:
                    self.COVERED_COUNT[sampleidx - 1] += 1
        elif len(self.speciesorder) != 0 and len(self.sampleidxlisttocount.keys()) != 0:
            for species in self.speciesorder:
                totaldepth = 0
                for sampleidx in self.sampleidxlisttocount[species]:
                    totaldepth += int(T[sampleidx])
                self.AVERAGE_DEPTH[self.speciesorder.index(species)] += totaldepth
                if totaldepth >= self.mindepth:
                    self.COVERED_COUNT[self.speciesorder.index(species)] += 1
                
            
    def getResult(self):
        """pecentage of cover,average depth
        """
        countlist = copy.deepcopy(self.COVERED_COUNT);average = copy.deepcopy(self.AVERAGE_DEPTH)
        del self.AVERAGE_DEPTH[:]
        del self.COVERED_COUNT[:]
        if self.speciesorder == [] and (not self.sampleidxlisttocount):
            self.COVERED_COUNT = [0] * self.total_samples
            self.AVERAGE_DEPTH = [0] * self.total_samples
        elif len(self.speciesorder) != 0 and len(self.sampleidxlisttocount.keys()) != 0:
            self.COVERED_COUNT = [0] * len(self.speciesorder)
            self.AVERAGE_DEPTH = [0] * len(self.speciesorder)
        return "empty", ([a / self.winsize for a in countlist], [a / self.winsize for a in average])
class Caculate_S_ObsExp_difference(Caculator):
    def __init__(self,mindepthtojudefixed,listOftargetpopvcfconfig,listOfrefpopvcffileconfig,dbvariantstoolstojudgeancestral,toplevelTablejudgeancestralname,outfileprewithpath):
        super().__init__()
        self.dbvariantstoolstojudgeancestral=dbvariantstoolstojudgeancestral
        self.topleveltablejudgeancestralname=toplevelTablejudgeancestralname
        self.MethodToSeqpoplist=[]
        self.mindepthtojudefixed=20
        self.flankseqfafile=None
        self.N_of_targetpop=len(listOftargetpopvcfconfig)
        self.N_of_refpop=len(listOfrefpopvcffileconfig)
        self.listOfpopvcfRecsmapByAChr=[]
        self.vcfnamelist=[]
        self.vcfnameKEY_vcfobj_pyBAMfilesVALUE={}
        self.outputname=outfileprewithpath
        for vcfconfigfilename in listOftargetpopvcfconfig[:]+listOfrefpopvcffileconfig[:]:
            self.listOfpopvcfRecsmapByAChr.append({})
            vcfconfig=open(vcfconfigfilename,"r")
            for line in vcfconfig:
                vcffilename_obj=re.search(r"vcffilename=(.*)",line.strip())
                if vcffilename_obj!=None:
                    vcfname=vcffilename_obj.group(1).strip()
                    self.vcfnamelist.append(vcfname)
                    self.outputname+=("_"+re.split(r"\.",re.search(r"[^/]*$",vcfname).group(0))[0])[:3]
                    self.vcfnameKEY_vcfobj_pyBAMfilesVALUE[vcfname]=[]
                    self.vcfnameKEY_vcfobj_pyBAMfilesVALUE[vcfname].append(VCFutil.VCF_Data(vcfname))
                elif line.split():
                    self.vcfnameKEY_vcfobj_pyBAMfilesVALUE[vcfname].append(pysam.Samfile(line.strip(),'rb'))
            vcfconfig.close()
            if re.search(r"indvd[^/]+",vcfname)!=None:
                self.MethodToSeqpoplist.append("indvd")
    
            elif re.search(r"pool[^/]+",vcfname)!=None:
                self.MethodToSeqpoplist.append("pool")
    
            else:
                print("vcfname must with 'pool' or 'indvd'")
                exit(-1)   

        self.currentchrID=None
        self.COUNT=0
        self.obsseq=[]
        self.CEXP=[]
        self.CfixedDerived=0
        self.freq_xaxisKEY_yaxisVALUERelation=None
        self.minsnps=10
        self.alignedSNP_absentinfo={}#{chrNo:[[pos1,REF,ALT,(),(),(),,],[pos1,REF,ALT,(),(),(),,],[],,]}
        self.dynamicIU_toptable_obj=None
    def __del__(self):
        self.flankseqfafile.close()
    def process(self,T):
        """T=[pos,ref,alt,pop1,pop2,.....,popn]
        in ordered as the self.vcfnamelise
        """
        if len(T[1]) != len(T[2]) or len(T[2])!=1  or len(T[2])!=1:
            return
        snp=self.dbvariantstoolstojudgeancestral.operateDB("select","select * from "+self.topleveltablejudgeancestralname+" where chrID='"+self.currentchrID+"' and snp_pos='"+str(T[0])+"'")
        if not snp or  snp[0][outgidx]==None or snp[0][outg2idx]==None or snp[0][5]==None:#needed info of SNPs are absent, snp[0][9] and snp[0][5] are dependent on the fellow code segment
#             print(self.currentchrID,T,"snp not find,skip")
#             print("append in to alignedSNP_absentinfo",snp,T)
            self.alignedSNP_absentinfo[self.currentchrID].append(T)
            return
        else:
            A_base_idx=100
            fanyadepthlist=re.split(r",",snp[0][outgidx]);taihudepthlist=re.split(r",",snp[0][outg2idx])
            if fanyadepthlist and taihudepthlist and (fanyadepthlist[0].strip()=="0" and int(fanyadepthlist[1]) >=self.mindepthtojudefixed and int(taihudepthlist[0])<=int(taihudepthlist[1])*seqerrorrate or (taihudepthlist[0].strip()=="0" and int(taihudepthlist[1])>=self.mindepthtojudefixed and int(fanyadepthlist[0])<=int(fanyadepthlist[1])*seqerrorrate) ):
                A_base_idx=1
            elif fanyadepthlist and taihudepthlist and( fanyadepthlist[1].strip()=="0" and int(fanyadepthlist[0])>=self.mindepthtojudefixed and int(taihudepthlist[1])<=int(taihudepthlist[0])*seqerrorrate or (taihudepthlist[1].strip()=="0" and int(taihudepthlist[0])>=self.mindepthtojudefixed and int(fanyadepthlist[1])<=int(fanyadepthlist[0])*seqerrorrate)):
                A_base_idx=0
            else:
#                 print("skip snp",snp[0][1],snp[0][7],snp[0][9],snp[0][11],snp[0][13])
                return

        ancestrallcontext=snp[0][5].strip()[0].upper()+snp[0][3+A_base_idx].strip().upper()+snp[0][5].strip()[2].upper()
        if "CG" in ancestrallcontext or "GC" in ancestrallcontext:
#             print("skip CG site",ancestrallcontext)
            return
        ##########x-axis
        countedAF=0;target_DAF_sum=0
        for tpopidx in range(3,self.N_of_targetpop+3):
            if T[tpopidx]==None:
                if len(self.vcfnameKEY_vcfobj_pyBAMfilesVALUE[self.vcfnamelist[tpopidx-3]])==1:
#                     print("skip this pos",T)
                    continue
                else:
#                     depth_linelist=self.depthobjlist[tpopidx-3].getdepthByPos_optimized(self.currentchrID,T[0])
                    sum_depth=0
                    for samfile in self.vcfnameKEY_vcfobj_pyBAMfilesVALUE[self.vcfnamelist[tpopidx-3]][1:]:
                        ACGTdep=samfile.count_coverage(self.currentchrID,T[0]-1,T[0])
                        for dep in ACGTdep:
                            sum_depth+=dep[0]
#                     for idx in self.species_idx_list[tpopidx-3][:]:
#                         sum_depth+=int(depth_linelist[idx])
                    if sum_depth>self.mindepthtojudefixed:
                        AF=0
                    else:
#                         print(sum_depth,"low coverage skip")
                        continue
            else:
                if self.MethodToSeqpoplist[tpopidx-3]=="indvd":
                    AF=float(re.search(r"AF=([\d\.e-]+)[;,]", T[tpopidx][0]).group(1))
                    AN = float(re.search(r"AN=([\d]+)[;,]", T[tpopidx][0]).group(1))
                    if AN<5:
                        continue
                elif self.MethodToSeqpoplist[tpopidx-3]=="pool":
                    refdep = 0;altalleledep = 0
                    AD_idx = (re.split(":", T[tpopidx][1])).index("AD")
                    for sample in T[tpopidx][2]:
                        if len(re.split(":", sample)) == 1:  # ./.
                            continue
                        AD_depth = re.split(",", re.split(":", sample)[AD_idx])
                        try :
                            refdep += int(AD_depth[0])
                            altalleledep += int(AD_depth[1])
                        except ValueError:
                            print(sample, end="|")
                    if (refdep==altalleledep and altalleledep==0) or altalleledep+ refdep<10:
                        continue
                    AF=altalleledep/(altalleledep+refdep)
            if A_base_idx==0:
                DAF=1-AF
            elif A_base_idx==1:
                DAF=AF
            target_DAF_sum+=DAF;countedAF+=1
        if  countedAF==0 or target_DAF_sum/countedAF==0:#:
#             print("skip this snp,because it fiexd as ancestral or no covered in this pos in target pops",T,snp)
            return
        target_DAF=target_DAF_sum/countedAF
        #########y-axis
        countedAF=0;rer_DAF_sum=0
        for rpopidx in range(3+self.N_of_targetpop,self.N_of_refpop+self.N_of_targetpop+3):
            if T[rpopidx]==None:
                if len(self.vcfnameKEY_vcfobj_pyBAMfilesVALUE[self.vcfnamelist[rpopidx-3]])==1:
#                     print("skip this snp",T)
                    continue
                else:
#                     depth_linelist=self.depthobjlist[rpopidx-3-self.N_of_targetpop].getdepthByPos_optimized(self.currentchrID,T[0])
                    sum_depth=0
                    for samfile in self.vcfnameKEY_vcfobj_pyBAMfilesVALUE[self.vcfnamelist[rpopidx-3]][1:]:
                        ACGTdep=samfile.count_coverage(self.currentchrID,T[0]-1,T[0])
                        for dep in ACGTdep:
                            sum_depth+=dep[0]
#                     for idx in self.species_idx_list[rpopidx-3-self.N_of_targetpop][:]:
#                         sum_depth+=int(depth_linelist[idx])
                    if sum_depth>self.mindepthtojudefixed:
                        AF=0
                    else:
#                         print(sum_depth,"low depth skip")
                        continue
            else:
                if self.MethodToSeqpoplist[rpopidx-3]=="indvd":
                    AF=float(re.search(r"AF=([\d\.]+)[;,]", T[rpopidx][0]).group(1))
                elif self.MethodToSeqpoplist[rpopidx-3]=="pool":
                    refdep = 0;altalleledep = 0
                    AD_idx = (re.split(":", T[rpopidx][1])).index("AD")
                    for sample in T[rpopidx][2]:
                        if len(re.split(":",sample))==1:
                            continue
                        AD_depth = re.split(",", re.split(":", sample)[AD_idx])
                        try :
                            refdep += int(AD_depth[0])
                            altalleledep += int(AD_depth[1])
                        except ValueError:
                            print(sample, end="|")
                    if refdep==altalleledep and altalleledep==0:
                        continue
                    AF=altalleledep/(altalleledep+refdep)
                if A_base_idx==0:
                    DAF=1-AF
                elif A_base_idx==1:
                    DAF=AF
                rer_DAF_sum+=DAF;countedAF+=1
        if  countedAF==0 :#or rer_DAF_sum==0
#             print("skip this snp,because it  no covered in this pos in ref pops",T,snp)
            return
        for a,b in sorted(self.freq_xaxisKEY_yaxisVALUERelation.keys()):
            if target_DAF>a and target_DAF<=b:
                self.CEXP.append(self.freq_xaxisKEY_yaxisVALUERelation[(a,b)])
                self.obsseq.append(rer_DAF_sum/countedAF)
                break
        
        self.COUNT+=1
        if rer_DAF_sum/countedAF==1:
            self.CfixedDerived+=1
    def getResult(self):
#         for chrom in self.alignedSNP_absentinfo.keys():
        if self.alignedSNP_absentinfo[self.currentchrID]!=[]:
            self.dynamicIU_toptable_obj.insertorUpdatetopleveltable(self.alignedSNP_absentinfo,self.flankseqfafile,50)
            No_Of_snpT=len(self.alignedSNP_absentinfo[self.currentchrID])
            print("self.CEXP SNP sites before process",len(self.CEXP))
            for whatever in range(No_Of_snpT):
                snpT=self.alignedSNP_absentinfo[self.currentchrID].pop(0)
#                     print(snpT,len(self.alignedSNP_absentinfo[chrom]))
                print(snpT[0:3],"process",len(self.alignedSNP_absentinfo[self.currentchrID]),self.COUNT)
                self.process(snpT)
            if len(self.alignedSNP_absentinfo[self.currentchrID])!=0:print("length of snpT",len(self.alignedSNP_absentinfo[self.currentchrID]));self.alignedSNP_absentinfo[self.currentchrID]=[]
            print("self.CEXP SNP sites processed",len(self.CEXP))    
        S1=0
        S2="NA"
        noofsnp=self.COUNT
        for i in range(len(self.obsseq)):

            S1+=self.obsseq[i]/self.CEXP[i]

        else:
            try:
                S1=math.log(S1/noofsnp)
            except ValueError:
                S1="MAX"#may be use the smallest self.obsseq[x] to get S1
            except:
                S1="NA"
        try:
            S2=math.log(np.sum(self.obsseq)/np.sum(self.CEXP))#(np.sum(self.obsseq)-self.CEXP)/np.std(self.obsseq,ddof=1)
        except:
            S2="NA"
        
        self.COUNT=0
        self.CEXP=[]
        self.obsseq=[]
        self.CfixedDerived=0
        """S1 is the corrected value, S2 is the not very appropriate value I used before 
        """       
        if S1=="NA"  or noofsnp<self.minsnps:
            return noofsnp,["NA","NA"]
        return noofsnp,[S1,S2]
# class Caculate_pairFst(Caculator):
#     def __init__(self,mindepthtojudefixed,listOftargetpopvcfconfig,listOfrefpopvcffileconfig,outfileprewithpath,minsnps=10):
#         super().__init__()
#         self.minsnps=minsnps
#         self.considerfixdiffinfst=False
#         self.MethodToSeqpoplist=[]
#         self.mindepthtojudefixed=mindepthtojudefixed
#         self.N_of_targetpop=len(listOftargetpopvcfconfig)
#         self.N_of_refpop=len(listOfrefpopvcffileconfig)
#         self.outputname=outfileprewithpath
#         self.vcfnamelist=[]#target ref
#         self.listOfpopvcfRecsmapByAChr=[]
#         self.vcfnameKEY_vcfobj_pyBAMfilesVALUE={}
#         for vcfconfigfilename in listOftargetpopvcfconfig[:]+listOfrefpopvcffileconfig[:]:
#             self.listOfpopvcfRecsmapByAChr.append({})
#             vcfconfig=open(vcfconfigfilename,"r")
#             for line in vcfconfig:
#                 vcffilename_obj=re.search(r"vcffilename=(.*)",line.strip())
#                 if vcffilename_obj!=None:
#                     vcfname=vcffilename_obj.group(1).strip()
#                     self.vcfnamelist.append(vcfname)
#                     self.outputname+=("_"+re.split(r"\.",re.search(r"[^/]*$",vcfname).group(0))[0])[:3]
#                     self.vcfnameKEY_vcfobj_pyBAMfilesVALUE[vcfname]=[]
#                     self.vcfnameKEY_vcfobj_pyBAMfilesVALUE[vcfname].append(VCFutil.VCF_Data(vcfname))
#                 elif line.split():
#                     self.vcfnameKEY_vcfobj_pyBAMfilesVALUE[vcfname].append(pysam.Samfile(line.strip(),'rb'))
#             vcfconfig.close()
#             if re.search(r"indvd[^/]+",vcfname)!=None:
#                 self.MethodToSeqpoplist.append("indvd")
#     
#             elif re.search(r"pool[^/]+",vcfname)!=None:
#                 self.MethodToSeqpoplist.append("pool")
#     
#             else:
#                 print("vcfname must with 'pool' or 'indvd'")
#                 exit(-1)   
#         self.currentchrID=None
#         self.COUNT=[(0,0)]*(self.N_of_refpop*self.N_of_targetpop)
#         self.CNK=[0]*(self.N_of_refpop*self.N_of_targetpop)
#         self.CDK=[0]*(self.N_of_refpop*self.N_of_targetpop)
#         self.CfixedDerived=0
#         self.freq_xaxisKEY_yaxisVALUERelation=None
#         self.minsnps=10
#     def process(self,T):
#         """T=[pos,ref,alt,pop1,pop2,.....,popn]
#         T is in the order as self.vcfnamelist
#         """
#         if len(T[1]) != len(T[2]) or len(T[2])!=1  or len(T[2])!=1:
#             return
# #         snp=self.dbvariantstoolstojudgeancestral.operateDB("select","select * from "+self.topleveltablejudgeancestralname+" where chrID='"+self.currentchrID+"' and snp_pos='"+str(T[0])+"'")
#         for TT_idx in range(3,3+self.N_of_targetpop):
#             TT=T[TT_idx]
#             refdep_1=0;altalleledep_1=0
#             if self.MethodToSeqpoplist[TT_idx-3]=="pool":
#                 if TT==None:
#                     if len(self.vcfnameKEY_vcfobj_pyBAMfilesVALUE[self.vcfnamelist[TT_idx-3]])==1:
#                         print("skip this pos",T)
#                         continue
#                     else:
#                         sum_depth=0
#                         for samfile in self.vcfnameKEY_vcfobj_pyBAMfilesVALUE[self.vcfnamelist[TT_idx-3]][1:]:
#                             ACGTdep=samfile.count_coverage(self.currentchrID,T[0]-1,T[0])
#                             for dep in ACGTdep:
#                                 sum_depth+=dep[0]
#                         if sum_depth>self.mindepthtojudefixed:
#                             refdep_1=sum_depth;altalleledep_1=0
#                         else:
#                             continue
#                 else:
#                     AD_idx_1 = (re.split(":", TT[1])).index("AD")  # gatk GT:AD:DP:GQ:PL
#             for RT_idx in range(3+self.N_of_targetpop,len(T)):
#                 RT=T[RT_idx]
#                 refdep_2=0;altalleledep_2=0
#                 
#         ##########x-axis
#         countedAF=0;target_DAF_sum=0
#         for tpopidx in range(3,self.N_of_targetpop+3):
#             if T[tpopidx]==None:
#                 if self.depthobjlist==[]:
#                     print("skip this pos",T)
#                     continue
#                 else:
#                     depth_linelist=self.depthobjlist[tpopidx-3].getdepthByPos_optimized(self.currentchrID,T[0])
#                     sum_depth=0
#                     for idx in self.species_idx_list[tpopidx-3][:]:
#                         sum_depth+=int(depth_linelist[idx])
#                     if sum_depth>self.mindepthtojudefixed:
#                         AF=0
#                     else:
#                         continue
#             else:
#                 if self.MethodToSeqpoplist[tpopidx-3]=="indvd":
#                     AF=float(re.search(r"AF=([\d\.e-]+)[;,]", T[tpopidx][0]).group(1))
#                     AN = float(re.search(r"AN=([\d]+)[;,]", T[tpopidx][0]).group(1))
#                     if AN<5:
#                         continue
#                 elif self.MethodToSeqpoplist[tpopidx-3]=="pool":
#                     refdep = 0;altalleledep = 0
#                     AD_idx = (re.split(":", T[tpopidx][1])).index("AD")
#                     for sample in T[tpopidx][2]:
#                         if len(re.split(":", sample)) == 1:  # ./.
#                             continue
#                         AD_depth = re.split(",", re.split(":", sample)[AD_idx])
#                         try :
#                             refdep += int(AD_depth[0])
#                             altalleledep += int(AD_depth[1])
#                         except ValueError:
#                             print(sample, end="|")
#                     if (refdep==altalleledep and altalleledep==0) or altalleledep+ refdep<10:
#                         continue
#                     AF=altalleledep/(altalleledep+refdep)
#             if A_base_idx==0:
#                 DAF=1-AF
#             elif A_base_idx==1:
#                 DAF=AF
#             target_DAF_sum+=DAF;countedAF+=1
#         if  countedAF==0 :#or target_DAF_sum/countedAF==0:
# #             print("skip this snp,because it fiexd as ancestral or no covered in this pos in target pops",T,snp)
#             return
#         target_DAF=target_DAF_sum/countedAF
#         #########y-axis
#         countedAF=0;rer_DAF_sum=0
#         for rpopidx in range(3+self.N_of_targetpop,self.N_of_refpop+self.N_of_targetpop+3):
#             if T[rpopidx]==None:
#                 if self.depthobjlist==[]:
#                     print("skip this snp",T)
#                     continue
#                 else:
#                     depth_linelist=self.depthobjlist[rpopidx-3-self.N_of_targetpop].getdepthByPos_optimized(self.currentchrID,T[0])
#                     sum_depth=0
#                     for idx in self.species_idx_list[rpopidx-3-self.N_of_targetpop][:]:
#                         sum_depth+=int(depth_linelist[idx])
#                     if sum_depth>self.mindepthtojudefixed:
#                         AF=0
#                     else:
#                         continue
#             else:
#                 if self.MethodToSeqpoplist[rpopidx-3-self.N_of_targetpop]=="indvd":
#                     AF=float(re.search(r"AF=([\d\.]+)[;,]", T[rpopidx][0]).group(1))
#                 elif self.MethodToSeqpoplist[rpopidx-3-self.N_of_targetpop]=="pool":
#                     refdep = 0;altalleledep = 0
#                     AD_idx = (re.split(":", T[rpopidx][1])).index("AD")
#                     for sample in T[rpopidx][2]:
#                         if len(re.split(":",sample))==1:
#                             continue
#                         AD_depth = re.split(",", re.split(":", sample)[AD_idx])
#                         try :
#                             refdep += int(AD_depth[0])
#                             altalleledep += int(AD_depth[1])
#                         except ValueError:
#                             print(sample, end="|")
#                     if refdep==altalleledep and altalleledep==0:
#                         continue
#                     AF=altalleledep/(altalleledep+refdep)
#                 if A_base_idx==0:
#                     DAF=1-AF
#                 elif A_base_idx==1:
#                     DAF=AF
#                 rer_DAF_sum+=DAF;countedAF+=1
#         if  countedAF==0 or rer_DAF_sum==0:
# #             print("skip this snp,because it  no covered in this pos in ref pops",T,snp)
#             return
#         for a,b in sorted(self.freq_xaxisKEY_yaxisVALUERelation.keys()):
#             if target_DAF>a and target_DAF<=b:
#                 self.CEXP+=self.freq_xaxisKEY_yaxisVALUERelation[(a,b)]
#                 break
#         self.obsseq.append(rer_DAF_sum/countedAF)
#         self.COUNT+=1
#         if rer_DAF_sum/countedAF==1:
#             self.CfixedDerived+=1
#     def getResult(self):
#         S1="NA"
#         S2="NA"
#         try:
#             S1=math.log(np.sum(self.obsseq)/self.CEXP)
#             S2=0#(np.sum(self.obsseq)-self.CEXP)/np.std(self.obsseq,ddof=1)
#         except:
#             S1="NA"
#             S2="NA"
#         noofsnp=self.COUNT
#         self.COUNT=0
#         self.CEXP=0
#         self.obsseq=[]
#         self.CfixedDerived=0
#         if S1=="NA" and S2=="NA" or noofsnp<self.minsnps:
#             return noofsnp,"NA"
#         return noofsnp,[S1,S2]
class Caculate_IS(Caculator):
    def __init__(self,mindepthtojudefixed,listOftargetpopvcfconfig,listOfrefpopvcffileconfig,minsnps=10):
        super().__init__()
        self.minsnps=minsnps
        self.considerfixdiffinfst=False
        self.MethodToSeqpoplist=[]
        self.mindepthtojudefixed=mindepthtojudefixed
        self.N_of_targetpop=len(listOftargetpopvcfconfig)
        self.N_of_refpop=len(listOfrefpopvcffileconfig)
        self.vcfnamelist=[]#target ref order
        self.listOfpopvcfRecsmapByAChr=[]
        self.vcfnameKEY_vcfobj_pyBAMfilesVALUE={}
        for vcfconfigfilename in listOftargetpopvcfconfig[:]+listOfrefpopvcffileconfig[:]:
            self.listOfpopvcfRecsmapByAChr.append({})
            vcfconfig=open(vcfconfigfilename,"r")
            for line in vcfconfig:
                vcffilename_obj=re.search(r"vcffilename=(.*)",line.strip())
                if vcffilename_obj!=None:
                    vcfname=vcffilename_obj.group(1).strip()
                    self.vcfnamelist.append(vcfname)
                    self.vcfnameKEY_vcfobj_pyBAMfilesVALUE[vcfname]=[]
                    self.vcfnameKEY_vcfobj_pyBAMfilesVALUE[vcfname].append(VCFutil.VCF_Data(vcfname))
                elif line.split():
                    self.vcfnameKEY_vcfobj_pyBAMfilesVALUE[vcfname].append(pysam.Samfile(line.strip(),'rb'))
            vcfconfig.close()
            if re.search(r"indvd[^/]+",vcfname)!=None:
                self.MethodToSeqpoplist.append("indvd")
    
            elif re.search(r"pool[^/]+",vcfname)!=None:
                self.MethodToSeqpoplist.append("pool")
    
            else:
                print("vcfname must with 'pool' or 'indvd'")
                exit(-1)   
        self.currentchrID=None
#         self.COUNT=[0]*((self.N_of_refpop+self.N_of_targetpop)*(self.N_of_refpop+self.N_of_targetpop-1)/2)
#         re.search(r"[^/]*$",vcfname).group(0)
        self.IS_Tinner={}#[0]*((self.N_of_targetpop-1)*self.N_of_targetpop/2)
        self.IS_Rinner={}#[0]*(self.N_of_refpop*(self.N_of_refpop-1)/2)
        self.IS_TR={}#[0]*self.N_of_refpop*self.N_of_targetpop
        self.combination_idx_list=list(combinations([i for i in range(self.N_of_refpop+self.N_of_targetpop)],2))
        print(self.combination_idx_list,"self.combination_idx_list")
        self.vcfname_combination=[];vcfname_combination=[]
        for pop_1_idx,pop_2_idx in self.combination_idx_list:
            if pop_1_idx<self.N_of_targetpop and pop_2_idx<self.N_of_targetpop:#
                
                vcfname_combination.append(re.split(r"\.",re.search(r"[^/]*$",self.vcfnamelist[pop_1_idx]).group(0))[0]+"_"+re.split(r"\.",re.search(r"[^/]*$",self.vcfnamelist[pop_2_idx]).group(0))[0])
                self.IS_Tinner[(pop_1_idx,pop_2_idx)]=[vcfname_combination[-1]]
            elif (pop_1_idx<self.N_of_targetpop and pop_2_idx>=self.N_of_targetpop) or (pop_2_idx<self.N_of_targetpop and pop_1_idx>=self.N_of_targetpop):
                
                vcfname_combination.append(re.split(r"\.",re.search(r"[^/]*$",self.vcfnamelist[pop_1_idx]).group(0))[0]+"_"+re.split(r"\.",re.search(r"[^/]*$",self.vcfnamelist[pop_2_idx]).group(0))[0])
                self.IS_TR[(pop_1_idx,pop_2_idx)]=[vcfname_combination[-1]]
            elif pop_1_idx>=self.N_of_targetpop and pop_2_idx>=self.N_of_targetpop:
                
                vcfname_combination.append(re.split(r"\.",re.search(r"[^/]*$",self.vcfnamelist[pop_1_idx]).group(0))[0]+"_"+re.split(r"\.",re.search(r"[^/]*$",self.vcfnamelist[pop_2_idx]).group(0))[0])
                self.IS_Rinner[(pop_1_idx,pop_2_idx)]=[vcfname_combination[-1]]
            else:
                print("what's wrong with CaculatorIS")
        print("order the vcfname_combination")

        for pop_1_idx,pop_2_idx in sorted(self.IS_Tinner.keys()):
            idx=self.combination_idx_list.index((pop_1_idx,pop_2_idx))
            self.vcfname_combination.append(vcfname_combination[idx])

        for pop_1_idx,pop_2_idx in sorted(self.IS_TR.keys()):
            idx=self.combination_idx_list.index((pop_1_idx,pop_2_idx))
            self.vcfname_combination.append(vcfname_combination[idx])
        for pop_1_idx,pop_2_idx in sorted(self.IS_Rinner.keys()):
            idx=self.combination_idx_list.index((pop_1_idx,pop_2_idx))
            self.vcfname_combination.append(vcfname_combination[idx])
                #RIinner            
        self.minsnps=7
        self.considerdepth=True
        print(self.vcfname_combination)
        print(self.IS_Rinner)
        print(self.IS_Tinner)
        print(self.IS_TR)
    def process(self,T):
        """T=[pos,ref,alt,pop1,pop2,.....,popn]
        T is in the order as self.vcfnamelist
        """
        if len(T[1]) != len(T[2]) or len(T[2])!=1  or len(T[2])!=1:
            return
        for pop_1_idx,pop_2_idx in self.combination_idx_list:
            TT=T[pop_1_idx+3]
            if TT==None:
                if len(self.vcfnameKEY_vcfobj_pyBAMfilesVALUE[self.vcfnamelist[pop_1_idx]])==1:
                    print("skip this pos",T)
                    continue
                elif self.considerdepth:
                    sum_depth=0
                    for samfile in self.vcfnameKEY_vcfobj_pyBAMfilesVALUE[self.vcfnamelist[pop_1_idx]][1:]:
                        ACGTdep=samfile.count_coverage(self.currentchrID,T[0]-1,T[0])
                        for dep in ACGTdep:
                            sum_depth+=dep[0]
                    times=1
                    if self.MethodToSeqpoplist[pop_1_idx]=="indvd" and self.vcfnamelist[pop_1_idx]!="/home/liurui/data/vcffiles/beijingref/shaoxing/shaoxingqingkeegg27.indvd.withindel.vcf":
                        times=7#this threshold is just for mallard and spotbilled
                    if sum_depth>=self.mindepthtojudefixed*times:
                        AF_1=0
                    else:
#                         print(self.vcfnamelist[pop_1_idx],sum_depth,"low coverage skip ,samfile",T[0],end="\t")
                        AF_1="NA"
                else:
                    AF_1=0
            else:
                if self.MethodToSeqpoplist[pop_1_idx]=="pool":
                    refdep_1=0;altalleledep_1=0
                    AD_idx_1 = (re.split(":", TT[1])).index("AD")  # gatk GT:AD:DP:GQ:PL
                    for sample in TT[2][:]:
                        if len(re.split(":", sample)) == 1 or re.split(":", sample)[AD_idx_1]==".":  # ./.
                            continue
                        AD_depth = re.split(",", re.split(":", sample)[AD_idx_1])
                        try:
                            refdep_1 += int(AD_depth[0])
                            altalleledep_1 += int(AD_depth[1])
                        except ValueError:
                            print("an ValueError except in int(AD_depth[]) in sample ",sample, end="|")
                    if  altalleledep_1+refdep_1<self.mindepthtojudefixed*0.8:
                        AF_1="NA"
                    AF_1=altalleledep_1/(altalleledep_1+refdep_1)
                elif self.MethodToSeqpoplist[pop_1_idx]=="indvd":
                    AF_1=float(re.search(r"AF=([\d\.e-]+);", TT[0]).group(1))
                    AN = float(re.search(r"AN=([\d]+);", TT[0]).group(1))
                    if AN<5:
                        AF_1="NA"
            RT=T[pop_2_idx+3]
            if RT==None:
                if len(self.vcfnameKEY_vcfobj_pyBAMfilesVALUE[self.vcfnamelist[pop_2_idx]])==1:
                    print("skip this pos",T)
                    continue
                else:
                    sum_depth=0
                    for samfile in self.vcfnameKEY_vcfobj_pyBAMfilesVALUE[self.vcfnamelist[pop_2_idx]][1:]:
                        ACGTdep=samfile.count_coverage(self.currentchrID,T[0]-1,T[0])
                        for dep in ACGTdep:
                            sum_depth+=dep[0]
                    times=1
                    if self.MethodToSeqpoplist[pop_2_idx]=="indvd" and self.vcfnamelist[pop_2_idx] !="/home/liurui/data/vcffiles/beijingref/shaoxing/shaoxingqingkeegg27.indvd.withindel.vcf":
                        times=7
                    if sum_depth>=self.mindepthtojudefixed*times:
                        AF_2=0
                    else:
#                         print(self.vcfnamelist[pop_2_idx],sum_depth,"low coverage skip ,samfile pos",T[0],end="\t")
                        AF_2="NA"
            else:
                if self.MethodToSeqpoplist[pop_2_idx]=="pool":
                    refdep_2=0;altalleledep_2=0
                    AD_idx_1 = (re.split(":", RT[1])).index("AD")  # gatk GT:AD:DP:GQ:PL
                    for sample in RT[2][:]:
                        if len(re.split(":", sample)) == 1 or re.split(":", sample)[AD_idx_1]==".":  # ./.
                            continue
                        AD_depth = re.split(",", re.split(":", sample)[AD_idx_1])
                        try:
                            refdep_2 += int(AD_depth[0])
                            altalleledep_2 += int(AD_depth[1])
                        except ValueError:
                            print("an ValueError except in int(AD_depth[]) in sample ",sample, end="|")
                    if  altalleledep_2+refdep_2<self.mindepthtojudefixed*0.8:
                        AF_2="NA"
                    AF_2=altalleledep_2/(altalleledep_2+refdep_2)
                elif self.MethodToSeqpoplist[pop_2_idx]=="indvd":
                    AF_2=float(re.search(r"AF=([\d\.e-]+);", RT[0]).group(1))
                    AN = float(re.search(r"AN=([\d]+);", RT[0]).group(1))
                    if AN<5:
                        AF_2="NA"
            #this snp is useable
            if AF_1!="NA" and AF_2!="NA":
                IS=1-abs(AF_1-AF_2)
            else:
                IS="NA"
            if pop_1_idx<self.N_of_targetpop and pop_2_idx<self.N_of_targetpop:
                self.IS_Tinner[(pop_1_idx,pop_2_idx)].append(IS)
            elif (pop_1_idx<self.N_of_targetpop and pop_2_idx>=self.N_of_targetpop) or (pop_2_idx<self.N_of_targetpop and pop_1_idx>=self.N_of_targetpop):
                self.IS_TR[(pop_1_idx,pop_2_idx)].append(IS)
            elif pop_1_idx>=self.N_of_targetpop and pop_2_idx>=self.N_of_targetpop:
                #RIinner
                self.IS_Rinner[(pop_1_idx,pop_2_idx)].append(IS)


    def getResult(self):
        noofsnp=[0]*(int((self.N_of_refpop+self.N_of_targetpop)*(self.N_of_refpop+self.N_of_targetpop-1)/2))
        ISlist=["NA"]*(int((self.N_of_refpop+self.N_of_targetpop)*(self.N_of_refpop+self.N_of_targetpop-1)/2))        
#         for k in range(len(self.vcfname_combination)):
#             noofsnp[k]=0;ISlist[k]

#         try:
        for pop_1_idx,pop_2_idx in sorted(self.IS_Tinner.keys()):
            templist=[]
            idx_as_title=self.vcfname_combination.index(self.IS_Tinner[(pop_1_idx,pop_2_idx)][0])
            for e in self.IS_Tinner[(pop_1_idx,pop_2_idx)][1:]:
                if e!="NA":
                    noofsnp[idx_as_title]+=1
                    templist.append(e)
            if noofsnp[idx_as_title]>self.minsnps:
                ISlist[idx_as_title]=np.mean(templist)
            self.IS_Tinner[(pop_1_idx,pop_2_idx)]=[self.IS_Tinner[(pop_1_idx,pop_2_idx)][0]]
        for pop_1_idx,pop_2_idx in sorted(self.IS_TR.keys()):
            templist=[]
            idx_as_title=self.vcfname_combination.index(self.IS_TR[(pop_1_idx,pop_2_idx)][0])
            for e in self.IS_TR[(pop_1_idx,pop_2_idx)][1:]:
                if e!="NA":
                    noofsnp[idx_as_title]+=1
                    templist.append(e)
            if noofsnp[idx_as_title]>self.minsnps:
                ISlist[idx_as_title]=np.mean(templist)
            self.IS_TR[(pop_1_idx,pop_2_idx)]=[self.IS_TR[(pop_1_idx,pop_2_idx)][0]]
        for pop_1_idx,pop_2_idx in sorted(self.IS_Rinner.keys()):
            templist=[]
            idx_as_title=self.vcfname_combination.index(self.IS_Rinner[(pop_1_idx,pop_2_idx)][0])
            for e in self.IS_Rinner[(pop_1_idx,pop_2_idx)][1:]:
                if e!="NA":
                    noofsnp[idx_as_title]+=1
                    templist.append(e)
            if noofsnp[idx_as_title]>self.minsnps:
                ISlist[idx_as_title]=np.mean(templist)
            self.IS_Rinner[(pop_1_idx,pop_2_idx)]=[self.IS_Rinner[(pop_1_idx,pop_2_idx)][0]]                          
#         for pop_1_idx,pop_2_idx in self.combination_idx_list:
#             if pop_1_idx<self.N_of_targetpop and pop_2_idx<self.N_of_targetpop:#
# 
#             elif (pop_1_idx<self.N_of_targetpop and pop_2_idx>=self.N_of_targetpop) or (pop_2_idx<self.N_of_targetpop and pop_1_idx>=self.N_of_targetpop):
# 
#             elif pop_1_idx>=self.N_of_targetpop and pop_2_idx>=self.N_of_targetpop:
# 
#             i+=1
        print(noofsnp)
        print(ISlist)
        return noofsnp,ISlist                    
class Caculate_Fst(Caculator):
    def __init__(self, MethodToSeqpop1="pool", MethodToSeqpop2="indvd", minsnps=10):
        super().__init__()
        self.minsnps = minsnps
        self.considerfixdiffinfst=False
        self.CNk = 0
        self.CDk = 0
        self.COUNTED = [0,0]#fst used snp,fixed difference snp
#         self.considerFixed = considerFixed
        self.MethodToSeqpop1 = MethodToSeqpop1
        self.MethodToSeqpop2 = MethodToSeqpop2
#         self.depthforcurrentchrom=None
        self.depthobjmap=None
        self.species_idx_map=None
        self.currentchrID=None
        self.pop1_indvdsormediandepth=None#=6#when pop1 is none at a pos,and no depth information
        self.pop2_indvdsormediandepth=None#=6
    def process(self, T, seqerrorrate=0.01):
        """T=[pos,ref,alt,pop1,pop2]"""
        if len(T[1]) != len(T[2]) or len(T[2])!=1  or len(T[2])!=1:
            return
        refdep_1 = 0;refdep_2 = 0
        altalleledep_1 = 0;altalleledep_2 = 0
#        T=(pos,REF,ALT,(INFO,FORMAT,sampleslist),(INFO,FORMAT,sampleslist))
        pop1 = T[3]
        pop2 = T[4]
#         dp4_1 = re.search(r"DP4=(\d*),(\d*),(\d*),(\d*)", pop1[0])
#         if dp4_1 != None:  # vcf from samtools
#         if None !=None:  # vcf from samtools
#             pass
#             refdep_1 = int(dp4_1.group(1)) + int(dp4_1.group(2))
#             altalleledep_1 = int(dp4_1.group(3)) + int(dp4_1.group(4))
#             dp4_2 = re.search(r"DP4=(\d*),(\d*),(\d*),(\d*)", pop2[0])
#             refdep_2 = int(dp4_2.group(1)) + int(dp4_2.group(2))
#             altalleledep_2 = int(dp4_2.group(3)) + int(dp4_2.group(4))
#         else:  # vcf from gatk
        if self.MethodToSeqpop1 == "pool":              
            if pop1==None:
                if self.depthobjmap==None:
                    refdep_1=self.pop1_indvdsormediandepth
                    altalleledep_1=0
                else:
                    depth_linelist=self.depthobjmap["vcfpop1_ref"].getdepthByPos_optimized(self.currentchrID,T[0])#  re.split(r"\t",self.depthforcurrentchrom["vcfpop1_ref"][int(T[0])-1])
                    sum_depth=0
                    for idx in self.species_idx_map["vcfpop1_ref"]:
                        sum_depth+=int(depth_linelist[idx])
                    if sum_depth>=self.pop1_indvdsormediandepth:
                        refdep_1=self.pop1_indvdsormediandepth
                        altalleledep_1=0
                    else:
                        return
            else:
                AD_idx_1 = (re.split(":", pop1[1])).index("AD")  # gatk GT:AD:DP:GQ:PL
                for sample in pop1[2][:]:
                    if len(re.split(":", sample)) == 1:  # ./.
                        continue
                    AD_depth = re.split(",", re.split(":", sample)[AD_idx_1])
                    try:
                        refdep_1 += int(AD_depth[0])
                        altalleledep_1 += int(AD_depth[1])
                    except ValueError:
                        print(sample, end="|")
        elif self.MethodToSeqpop1 == "indvd":
            if pop1==None:
                if self.depthobjmap==None:
                    refdep_1=self.pop1_indvdsormediandepth
                    altalleledep_1=0
                else:
                    depth_linelist=self.depthobjmap["vcfpop1_ref"].getdepthByPos_optimized(self.currentchrID,int(T[0]))#re.split(r"\t",self.depthforcurrentchrom["vcfpop1_ref"][int(T[0])-1])
                    sum_depth=0
                    for idx in self.species_idx_map["vcfpop1_ref"]:
                        if int(depth_linelist[idx])>4:
                            sum_depth+=1
#                         sum_depth+=int(depth_linelist[idx])
#                     if sum_depth>=self.pop1_indvdsormediandepth:
                    refdep_1=sum_depth*2
                    altalleledep_1=0
#                     else:
#                         return
            else:
                AN = int(re.search(r"AN=(\d+);", pop1[0]).group(1))
                AC = int(re.search(r"AC=(\d+);", pop1[0]).group(1))
                refdep_1 = AN - AC
                altalleledep_1 = AC
        if self.MethodToSeqpop2 == "pool":
            
            if pop2==None:
                if self.depthobjmap==None:
                    refdep_2=self.pop2_indvdsormediandepth
                    altalleledep_2=0
                else:
                    depth_linelist=self.depthobjmap["vcfpop2"].getdepthByPos_optimized(self.currentchrID,int(T[0]))#re.split(r"\t",self.depthforcurrentchrom["vcfpop2"][int(T[0])-1])
                    sum_depth=0
                    for idx in self.species_idx_map["vcfpop2"]:
                        sum_depth+=int(depth_linelist[idx])
                    if sum_depth>=self.pop2_indvdsormediandepth:
                        refdep_2=self.pop2_indvdsormediandepth
                        altalleledep_2=0
                    else:
                        return
            else:
                AD_idx_2 = (re.split(":", pop2[1])).index("AD")
                for sample in pop2[2][:]:
                    if len(re.split(":", sample)) == 1:  # ./.
                        continue
                    AD_depth = re.split(",", re.split(":", sample)[AD_idx_2])
                    try:
                        refdep_2 += int(AD_depth[0])
                        altalleledep_2 += int(AD_depth[1])
                    except ValueError:
                        print(sample, end="|")
        elif self.MethodToSeqpop2 == "indvd":
            if pop2==None:
                if self.depthobjmap==None:
                    refdep_2=self.pop2_indvdsormediandepth
                    altalleledep_2=0
                else:
                    depth_linelist=self.depthobjmap["vcfpop2"].getdepthByPos_optimized(self.currentchrID,int(T[0]))#re.split(r"\t",self.depthforcurrentchrom["vcfpop2"][int(T[0])-1])

                    sum_depth=0
                    for idx in self.species_idx_map["vcfpop2"]:
                        if int(depth_linelist[idx])>1:
                            sum_depth+=1

                    refdep_2=sum_depth*2
                    altalleledep_2=0

            else:
                AN = int(re.search(r"AN=(\d+);", pop2[0]).group(1))
                AC = int(re.search(r"AC=(\d+);", pop2[0]).group(1))
                refdep_2 = AN - AC
                altalleledep_2 = AC

              
                
#         if  (refdep_1 <= seqerrorrate * (refdep_1 + altalleledep_1) or refdep_2 <= seqerrorrate * (refdep_2 + altalleledep_2)):
#             return  # NOTICT HERE
        if refdep_1==0 and refdep_2==0:#skip both fixed as alt
            return
        if ((refdep_1 + altalleledep_1 - 1) * (refdep_1 + altalleledep_1))==0 or  ((refdep_2 + altalleledep_2 - 1) * (refdep_2 + altalleledep_2))==0:
            return
        if (refdep_1==0 and altalleledep_2==0 and altalleledep_1>=self.pop1_indvdsormediandepth and refdep_2>=self.pop2_indvdsormediandepth) or (altalleledep_1==0 and refdep_2==0 and refdep_1>=self.pop1_indvdsormediandepth and altalleledep_2>=self.pop2_indvdsormediandepth):#fixed difference
            self.COUNTED[1]+=1
            if self.considerfixdiffinfst:
                print(T,"fixdifferent in Fst")
                pass
            else:
                print(T,"fixdiffernet not in Fst")
                return
        self.COUNTED[0] += 1
        h_1 = refdep_1 * altalleledep_1 / ((refdep_1 + altalleledep_1 - 1) * (refdep_1 + altalleledep_1))
        h_2 = refdep_2 * altalleledep_2 / ((refdep_2 + altalleledep_2 - 1) * (refdep_2 + altalleledep_2))
        Nk = ((refdep_1 / (refdep_1 + altalleledep_1) - refdep_2 / (refdep_2 + altalleledep_2)) ** 2 - h_1 / (refdep_1 + altalleledep_1) - h_2 / (refdep_2 + altalleledep_2))
        self.CNk += Nk
        self.CDk += (Nk + h_1 + h_2)
#         print("self.CNk",self.CNk,"NK",Nk,"self.CDk",self.CDk)
    def getResult(self):
        Fst = 'NA'
        try:
            Fst = self.CNk / self.CDk
        except ZeroDivisionError:
            # print("the Fst value of currentwindow is dividsion by zero,so set it to be NA")
            Fst = 'NA'
#         if self.COUNTED<=self.minsnps:
#             Fst='NA'
        self.CDk = 0
        self.CNk = 0
        noofsnp = self.COUNTED[0]
        nooffixdifference=self.COUNTED[1]
        self.COUNTED = [0,0]
        if noofsnp<self.minsnps:
            Fst="NA"
        return [noofsnp,nooffixdifference], Fst
