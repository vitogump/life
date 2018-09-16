'''
Created on 2018年9月14日

@author: Dr.liu
'''
from optparse import OptionParser
import os,re,copy,random

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
    print("name","\t".join(intersectionPoslist),file=fml)
    for k in intersectionPed.keys():
        if  totalAnswerMap[k]!="n" and totalAnswerMap[k]!="y":
            print(totalAnswerMap[k]);continue
        print("".join(intersectionPed[k][6:]),totalAnswerMap[k],file=fml)
    fml.close()
    #    print merged  pos
    tempmapfile=open("temp.map",'w');temppedfile=open("temp.ped",'w')
    for p in intersectionPoslist:
        print(p.split("_")[0],p,p.split("_")[1],p.split("_")[1],sep="\t",file=tempmapfile)
    for ind in intersectionPed.keys():
        print(ind,ind,"0\t0\t1\t1","\t".join([e[0]+"\t"+e[1] for e in intersectionPed[ind]]),sep="\t",file=temppedfile)
    #  random select a elem in intersectionPoslist  
    #print(totalAnswerMap)
    undistinguished=list(intersectionPed.keys())
    while undistinguished!=[] and intersectionPoslist!=[]:
        print(intersectionPoslist.pop())