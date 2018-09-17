'''
Created on 2018年9月14日

@author: Dr.liu
'''
from optparse import OptionParser
import os,re,copy,random,itertools
from sklearn.metrics import accuracy_score
import pandas as pd
parser = OptionParser()
parser.add_option("-g", "--trainningSetGenotype_MapPedlikeFormat", dest="TsetGeno",action="append",nargs=2,help="default sep by '\t' ; warnning! this file may not sorted by pos,ped each col is a pos genotype")
parser.add_option("-l", "--samplelabelfile", dest="samplelabel",help="only use 'y' and 'n' samples")
(options, args) = parser.parse_args()

DIPLOTYPES = ['A', 'C', 'G', 'K',"k","m", 'M', 'N', 'S',"s","r", 'R', 'T', 'W',"w","y", 'Y']
PAIRS = ['AA', 'CC', 'GG', 'GT',"TG","CA", 'AC', 'NN', 'CG',"GC","GA",'AG', 'TT', 'AT',"TA",'TC','CT']
diploHaploDict = dict(zip(DIPLOTYPES,PAIRS))
haploDiploDict = dict(zip(PAIRS,DIPLOTYPES))

if __name__ == '__main__':
    totalPositionlist=[];intersectionPoslist=[];intersectionPed={};totalAnswerMap={}
    for mappedPair,answerfile in options.TsetGeno:
        print(mappedPair)
        curposlist=[];curpeddict={};curanswermap={}
        mapf=open(mappedPair+".map","r");pedf=open(mappedPair+".ped","r")
        answerF=open(answerfile,'r')
        for ansline in answerF:
            anslist=re.split(r"\t",ansline.strip())
            curanswermap[anslist[0]]=anslist[-1]
        
        for mline in mapf:
            curposlist.append(re.split(r"\t",mline.strip())[1].strip())
        for indgeno in pedf:
            pedlist=re.split(r"\t",indgeno.strip())

            indname=pedlist[1].strip()+mappedPair if pedlist[1].strip() in intersectionPed.keys() else pedlist[1].strip()
            if indname not in curanswermap.keys():
                if indname.find("_dup")!=-1:
                    totalAnswerMap[indname]=curanswermap[indname[:indname.find("_dup")]]
                elif mappedPair in indname:
                    totalAnswerMap[indname]=curanswermap[indname[:indname.find(mappedPair)]]
            else:
                totalAnswerMap[indname]=  curanswermap[indname]
            curpeddict[indname]=pedlist[6:]
            intersectionPed[indname]=[]
        
            

        #uniq samples
        pass
        #intersection
        if intersectionPoslist==[]:
            intersectionPoslist=copy.deepcopy(curposlist)
            intersectionPed=copy.deepcopy(curpeddict)
        else:
            #add new ind in intersectionPed, remove posD in each ind of intersectionPed
            setpos=set(intersectionPoslist).intersection(set(curposlist))
            print(setpos)
            for pos in reversed(intersectionPoslist):
                if pos in setpos:# insert or pop
                    for indI in curpeddict.keys():
                        intersectionPed[indI].insert(0,curpeddict[indI][curposlist.index(pos)])
                else:# insert or pop
                    for indD in intersectionPed.keys():
                        if indD in curpeddict.keys():continue
                        intersectionPed[indD].pop(intersectionPoslist.index(pos))
                    intersectionPoslist.remove(pos)
    fml=open("formachinelearnning",'w')
    print("name","\t".join(intersectionPoslist),"breed",sep="\t",file=fml)
    for k in intersectionPed.keys():
        if  totalAnswerMap[k]!="n" and totalAnswerMap[k]!="y":
            print(k,totalAnswerMap[k]);continue
        print(k,"\t".join(intersectionPed[k][:],),totalAnswerMap[k],sep="\t",file=fml)
    fml.close()
    #    print merged  pos
    tempmapfile=open("temp.map",'w');temppedfile=open("temp.ped",'w')
    GenoMapBypos={}
    for p in intersectionPoslist:
        GenoMapBypos[p]=[]
        print(p.split("_")[0],p,p.split("_")[1],p.split("_")[1],sep="\t",file=tempmapfile)
    for ind in intersectionPed.keys():
        for gidx in range(len(intersectionPed[ind])):
            GenoMapBypos[intersectionPoslist[gidx]].append(haploDiploDict[intersectionPed[ind][gidx]].upper())
        print(ind,ind,"0\t0\t1\t1","\t".join([e[0]+"\t"+e[1] for e in intersectionPed[ind]]),sep="\t",file=temppedfile)
    #  random select a elem in intersectionPoslist  
    #print(totalAnswerMap)
    def bestAccur(genolistOfApos,turelist):
#     for pos in GenoMapBypos.keys():
        glistOfApos=genolistOfApos.tolist()
        a=set(glistOfApos);
        if "NN" in a: a.remove("NN")
        triG=list(a)
        comb_triG=itertools.permutations(triG,2)
        t=0;Y=triG[0];N=triG[1]
        for a,b in comb_triG:
            glistOfApos=genolistOfApos.tolist()
            for d_idx in range(len(glistOfApos)):
                if glistOfApos[d_idx]==a:
                    glistOfApos[d_idx]="y"
                elif glistOfApos[d_idx]==b:
                    glistOfApos[d_idx]="n"
            ca=accuracy_score(glistOfApos,turelist)
            print(a,b,ca)
            if t<=ca:
                t=ca;Y=a;N=b
#                 print(t,a,b)
                
        return t,Y,N
    fff=pd.read_csv("formachinelearnning",sep="\t")
    AccurLaccordingPos=[]
    for pos in fff.columns:
        if pos !="name" and pos !="breed":
            AccurLaccordingPos.append((bestAccur(fff[pos],fff["breed"]),pos))
    AccurLaccordingPos.sort(key=lambda x:x[0])
    undistinguished=list(intersectionPed.keys())
    fourposcombin = itertools.combinations(AccurLaccordingPos[-8:],4)
    print(len(AccurLaccordingPos[-8:]),AccurLaccordingPos[-8:])
    print(len(list(fourposcombin)))
    for pos1,pos2,pos3,pos4 in fourposcombin:
        for i in range(len(fff["breed"])):
            print(fff.loc[i,[pos1[1],pos2[1],pos3[1],pos4[1]]],fff.loc[i,["breed"]])

#         print(AccurLaccordingPos.pop())