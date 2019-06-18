'''
Created on 2019年6月4日

@author: RuiLiu
'''
import math,re,random

inputsnpfile=open("all.snp",'r')
genofile = open( "pseudoForAdmixedTools.geno", "w")
snpfile = open( "pseudoForAdmixedTools.snp", "w")
indfile = open("pseudoForAdmixedTools.ind", 'w')

mallard=["Mallard","Mallard","Mallard","Mallard","Mallard","Mallard","Mallard","Mallard","Mallard","Mallard","Mallard","Mallard","Mallard","Mallard"]
spotbilled=["Spot-billed","Spot-billed","Spot-billed","Spot-billed","Spot-billed","Spot-billed","Spot-billed","Spot-billed","Spot-billed","Spot-billed","Spot-billed","Spot-billed","Spot-billed"]
domestic=["Domestic","Domestic","Domestic","Domestic","Domestic","Domestic","Domestic","Domestic","Domestic","Domestic"]
fanya=["Muscovy","Muscovy","Muscovy","Muscovy","Muscovy","Muscovy","Muscovy","Muscovy","Muscovy","Muscovy"]

if __name__ == '__main__':
    inputsnpfile.readline()
    #print ind file
    for ind in mallard+spotbilled+domestic+fanya:
        print(ind,"U",ind,sep="\t",file=indfile)
    chrom=-1;c='KB742833'
    for line in inputsnpfile:
        genolistOrderbySamplelist = []
        linelist = re.split(r"\s+", line)
        #print snp file
        if c!=linelist[0]:
            chrom+=1;c=linelist[0]        
        print(linelist[0][:-2]+"_"+linelist[1],str(chrom),str(round(5.4696217209617786e-8 * int(linelist[1]),6)),linelist[1],linelist[2],linelist[3],sep="\t",file=snpfile)

        #print geno file mallard
        P1NofAlle1 = round(float(linelist[4])*len(mallard)*2);print(P1NofAlle1,end="mallard:\t")
        s=0
        for dind in range(1,P1NofAlle1,2):
            genolistOrderbySamplelist.append("0");s=dind+1
        if P1NofAlle1 != s:# True when P1NofAlle1 is an odd number
            genolistOrderbySamplelist.append("1");s+=1
        print(s,"==",P1NofAlle1)
        genolistOrderbySamplelist.extend(["2"]*(len(mallard)-math.ceil(s/2)))#;random.shuffle(genolistOrderbySamplelist) s/2 has no remainder or 0.5
        #P2derFreq    spotbilled
        P2NofAlle1 = round(float(linelist[5])*len(spotbilled)*2);print(P2NofAlle1,end="spotbilled:\t")
        s=0
        for dind in range(1,P2NofAlle1,2):
            genolistOrderbySamplelist.append("0");s=dind+1
        if P2NofAlle1 != s:# True when P2NofAlle1 is an odd number
            genolistOrderbySamplelist.append("1");s+=1
        print(s,"==",P2NofAlle1)
        genolistOrderbySamplelist.extend(["2"]*(len(spotbilled)-math.ceil(s/2)))#s/2 has no remainder or 0.5
        #P3derFreq    domestic
        s=0
        P3NofAlle1 = round(float(linelist[6])*len(domestic)*2);print(P3NofAlle1,end="domestic:\t")
        for dind in range(1,P3NofAlle1,2):
            genolistOrderbySamplelist.append("0");s=dind+1
        if P3NofAlle1 != s:# True when P3NofAlle1 is an odd number
            genolistOrderbySamplelist.append("1");s+=1
        print(s,"==",P3NofAlle1)
        genolistOrderbySamplelist.extend(["2"]*(len(domestic)-math.ceil(s/2)))#s/2 has no remainder or 0.5
        #P4derFreq    Muscovy
        s=0
        P4NofAlle1 = round(float(linelist[7])*len(fanya)*2);print(P4NofAlle1,end="Muscovy:\t")
        for dind in range(1,P4NofAlle1,2):
            genolistOrderbySamplelist.append("0");s=dind+1
        if P4NofAlle1 != s:# True when P4NofAlle1 is an odd number
            genolistOrderbySamplelist.append("1");s+=1
        print(s,"==",P4NofAlle1)
        genolistOrderbySamplelist.extend(["2"]*(len(fanya)-math.ceil(s/2)))#s/2 has no remainder or 0.5
        print("".join(genolistOrderbySamplelist),file=genofile)
        if len(genolistOrderbySamplelist)>47:
            break
    inputsnpfile.close()
    genofile.close()
    snpfile.close()
    indfile.close()
        
        
        
        