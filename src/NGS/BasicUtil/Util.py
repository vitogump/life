# -*- coding: UTF-8 -*-
import copy
import re,sys
import numpy,math
from multiprocessing.dummy import Pool
import pickle
import src.NGS.BasicUtil.DBManager as dbm
import config
from config import vcfdbname

'''
Created on 2013-6-30

@author: rui
'''
ip=config.ip
username=config.username
password=config.password
webdbname=config.webdbname
vcfdbname=config.vcfdbname
def alinmultPopSnpPos(vcfMaplist,jointmode="i"):
    """input:
    two or more map fomart like this [chrNo:[(pos,REF,ALT,INFO,FORMAT,sample,...),(pos,REF,ALT,INFO,FORMAT,sample,...),,,,,],{chrNo:[]},,,,,,]
    output:
    one map like this {chrNo:[(pos,REF,ALT,(INFO,FORMAT,sample,...),(INFO,FORMAT,sample,...)),(,,,(),()),,,,,],chrNo:[],,,}
                                            from pop1                        from pop2
    pop in the order of vcfMaplist
    """
    multipleVcfMap={}
    if len(vcfMaplist)==1 or jointmode=="o" or jointmode=="l":

        for currentChrom in vcfMaplist[0].keys():
            multipleVcfMap[currentChrom]=[]
            for SNPrec in vcfMaplist[0][currentChrom]:
                posInPop1 = SNPrec[0]#;print(posInPop1,file=open("testpos_8rep.txt0",'a'))
                RefInPop1 = SNPrec[1]
                AltInPop1 = SNPrec[2]
                multipleVcfMap[currentChrom].append([posInPop1,RefInPop1,AltInPop1,SNPrec[3:]])
        if len(vcfMaplist)==1:
            return copy.deepcopy(multipleVcfMap)
        vcfMap_obj_idx=0
        for vcfMap in vcfMaplist[1:]:
            vcfMap_obj_idx+=1
            for currentChrom in vcfMap.keys():
                for SNPrec in vcfMap[currentChrom]:
                    posInPop1 = SNPrec[0]#;print(posInPop1,file=open("testpos_8rep.txt"+str(vcfMap_obj_idx),'a'))
                    RefInPop1 = SNPrec[1]
                    AltInPop1 = SNPrec[2]
                    low=0;high=len(multipleVcfMap[currentChrom])-1
                    while low<=high:
                        mid = (low + high)>>1
                        if multipleVcfMap[currentChrom][mid][0]<posInPop1:
                            low=mid+1
                        elif multipleVcfMap[currentChrom][mid][0]>posInPop1:
                            high=mid-1
                        else:
                            #this is a tempera way, it would make a mislead in allele freq,as some pop's freq refer to alt1,other pop's freq refer to alt2,it cant tell, one can use alinmultPopSnpPos_diffrefalt() func as a alternate strategy
                            if AltInPop1 != multipleVcfMap[currentChrom][mid][2]:#differ alt alle
                                multipleVcfMap[currentChrom][mid][2]=(multipleVcfMap[currentChrom][mid][2]+","+RefInPop1+AltInPop1)
                                
                            fillNoneNum=vcfMap_obj_idx-(len(multipleVcfMap[currentChrom][mid])-3)
                            for i in range(fillNoneNum):
                                multipleVcfMap[currentChrom][mid].append(None)
                            multipleVcfMap[currentChrom][mid].append(SNPrec[3:])
                            break
                    else:
                        if jointmode=="o":
                            insertelem=[posInPop1,RefInPop1,AltInPop1]
                            for i in range(0,vcfMap_obj_idx):
                                insertelem.append(None)
                            insertelem.append(SNPrec[3:])
                            multipleVcfMap[currentChrom].insert(low,insertelem)
    #list(multipleVcfMap.keys())[0]==currentChrom
    #when a pos only exist in the former several pops,but not exist in the rear several pops (ie the length of each element of the  multipleVcfMap[currentChrom] are not same in some case),the loop block under are neccessary
        for REC_idx in range(0,len(multipleVcfMap[list(multipleVcfMap.keys())[0]])):
            #
            for i in range(len(vcfMaplist)+3-len(multipleVcfMap[list(multipleVcfMap.keys())[0]][REC_idx])):
                multipleVcfMap[list(multipleVcfMap.keys())[0]][REC_idx].append(None)

        return copy.deepcopy(multipleVcfMap)
    
    for currentChrom in vcfMaplist[0].keys():
#             self.FstMapByChrom[currentChrom] = []
        multipleVcfMap[currentChrom]=[]
        for SNPrec in vcfMaplist[0][currentChrom]:
            posInPop1 = SNPrec[0]
            RefInPop1 = SNPrec[1]
            AltInPop1 = SNPrec[2]
            elementToAppend=[posInPop1,RefInPop1,AltInPop1,SNPrec[3:]]
            if len(vcfMaplist)==1:
                multipleVcfMap[currentChrom].append(elementToAppend)
                continue
            for vcfMap_obj_idx in range(1,len(vcfMaplist[:])):
                vcfMap_obj=vcfMaplist[vcfMap_obj_idx]
                if currentChrom not in vcfMap_obj or len(vcfMap_obj[currentChrom])==0:
                    print("alinmultPopSnpPos",currentChrom,"didn't find in vcfMap2")
#                     if jointmode=="i":
                    break
#                     elif jointmode=="o":
#                         if vcfMap_obj_idx!=len(vcfMaplist)-1:
#                             elementToAppend.append(None)
#                         else:
#                             elementToAppend.append(None)
#                             multipleVcfMap[currentChrom].append(elementToAppend)
                low = 0
                high = len(vcfMap_obj[currentChrom]) - 1
                
                if re.search(r"[A-Za-z]+,[A-Za-z]+", AltInPop1) != None:  # multiple allels
                    continue
#                dp4 = re.search(r"DP4=(\d*),(\d*),(\d*),(\d*)", SNPrec[3])
#                 print(dp4.group(0))
                
                while low <= high:
                    mid = (low + high)>>1
                    if vcfMap_obj[currentChrom][mid][0]<posInPop1:
                        low=mid+1
                    elif vcfMap_obj[currentChrom][mid][0]>posInPop1:
                        high=mid-1
                    else:
                        if AltInPop1 == vcfMap_obj[currentChrom][mid][2]:#same alt alle
                            if vcfMap_obj_idx!=len(vcfMaplist)-1:
                                elementToAppend.append(vcfMap_obj[currentChrom][mid][3:])
                            elif vcfMap_obj_idx==len(vcfMaplist)-1:
                                elementToAppend.append(vcfMap_obj[currentChrom][mid][3:])
                                multipleVcfMap[currentChrom].append(elementToAppend)
#                         elif jointmode=="i":
#                         print("skip the different allele rec",currentChrom,posInPop1,AltInPop1,vcfMap_obj[currentChrom][mid][2])
#                             print(currentChrom,posInPop1,AltInPop1,vcfMap_obj[currentChrom][mid][2],"different alt allele,should skip this rec,but i have no time to improve this now")
#                         elif innerjoin_outjoin=="o":
#                             if vcfMap_obj_idx!=len(vcfMaplist)-1:
#                                 elementToAppend.append(None)
#                             elif vcfMap_obj_idx==len(vcfMaplist)-1:
#                                 elementToAppend.append(None)
#                                 multipleVcfMap[currentChrom].append(elementToAppend)
                        break
                else:
                    if jointmode=="i":
                        #ignore the rec
                        break
#                     elif jointmode=="o":
#                         if vcfMap_obj_idx!=len(vcfMaplist)-1:
#                             elementToAppend.append(None)
#                         elif vcfMap_obj_idx==len(vcfMaplist)-1:
#                             elementToAppend.append(None)
#                             multipleVcfMap[currentChrom].append(elementToAppend)                              
#                     print("snp not found in vcfMap2",SNPrec)
#                     self.doubleVcfMap[currentChrom].append(SNPrec+)
    return copy.deepcopy(multipleVcfMap)
def alinmultPopSnpPos_diffrefalt(vcfMaplist,jointmode="i"):
    """input:
    two or more map fomart like this [chrNo:[(pos,REF,ALT,INFO,FORMAT,sample,...),(pos,REF,ALT,INFO,FORMAT,sample,...),,,,,],chrNo:[],,,,,,]
    output:
    one map like this {chrNo:[(pos,(REF,ALT,INFO,FORMAT,sample,...),(REF,ALT,INFO,FORMAT,sample,...)),(,,,(),()),,,,,],chrNo:[],,,}
                                            from pop1                        from pop2
    """
    multipleVcfMap={}
    if len(vcfMaplist)==1 or jointmode=="o" or jointmode=="l":

        for currentChrom in vcfMaplist[0].keys():
            multipleVcfMap[currentChrom]=[]
            for SNPrec in vcfMaplist[0][currentChrom]:
                posInPop1 = SNPrec[0]#;print(posInPop1,file=open("testpos_8rep.txt0",'a'))

                multipleVcfMap[currentChrom].append([posInPop1,SNPrec[1:]])
        if len(vcfMaplist)==1:
            return copy.deepcopy(multipleVcfMap)
        vcfMap_obj_idx=0
        for vcfMap in vcfMaplist[1:]:
            vcfMap_obj_idx+=1
            for currentChrom in vcfMap:
                for SNPrec in vcfMap[currentChrom]:
                    posInPop1 = SNPrec[0]#;print(posInPop1,file=open("testpos_8rep.txt"+str(vcfMap_obj_idx),'a'))

                    low=0;high=len(multipleVcfMap[currentChrom])-1
                    while low<=high:
                        mid = (low + high)>>1
                        if multipleVcfMap[currentChrom][mid][0]<posInPop1:
                            low=mid+1
                        elif multipleVcfMap[currentChrom][mid][0]>posInPop1:
                            high=mid-1
                        else:
#                             if AltInPop1 == multipleVcfMap[currentChrom][mid][2]:#same alt alle
                            fillNoneNum=vcfMap_obj_idx-(len(multipleVcfMap[currentChrom][mid])-3)
                            for i in range(fillNoneNum):
                                multipleVcfMap[currentChrom][mid].append(None)
                            multipleVcfMap[currentChrom][mid].append(SNPrec[1:])
                            break
                    else:
                        if jointmode=="o":
                            insertelem=[posInPop1]
                            for i in range(0,vcfMap_obj_idx):
                                insertelem.append(None)
                            insertelem.append(SNPrec[1:])
                            multipleVcfMap[currentChrom].insert(low,insertelem)

    #list(multipleVcfMap.keys())[0]==currentChrom
    #when a pos only exist in the former several pops,but not exist in the rear several pops,the loop block under are neccessary
#         print(multipleVcfMap[list(multipleVcfMap.keys())[0]],len(multipleVcfMap[list(multipleVcfMap.keys())[0]]))
        for REC_idx in range(0,len(multipleVcfMap[list(multipleVcfMap.keys())[0]])):
            #
            for i in range(len(vcfMaplist)+1-len(multipleVcfMap[list(multipleVcfMap.keys())[0]][REC_idx])):
                multipleVcfMap[list(multipleVcfMap.keys())[0]][REC_idx].append(None)

        return copy.deepcopy(multipleVcfMap)
    
    for currentChrom in vcfMaplist[0].keys():
#             self.FstMapByChrom[currentChrom] = []
        multipleVcfMap[currentChrom]=[]
        for SNPrec in vcfMaplist[0][currentChrom]:
            posInPop1 = SNPrec[0]

            elementToAppend=[posInPop1,SNPrec[1:]]
            if len(vcfMaplist)==1:
                multipleVcfMap[currentChrom].append(elementToAppend)
                continue
            for vcfMap_obj_idx in range(1,len(vcfMaplist[:])):
                vcfMap_obj=vcfMaplist[vcfMap_obj_idx]
                if currentChrom not in vcfMap_obj or len(vcfMap_obj[currentChrom])==0:
                    print("alinmultPopSnpPos",currentChrom,"didn't find in vcfMap2")

                    break

                low = 0
                high = len(vcfMap_obj[currentChrom]) - 1
                
                if re.search(r"[A-Za-z]+,[A-Za-z]+", SNPrec[2]) != None:  # multiple allels
                    continue
         
                while low <= high:
                    mid = (low + high)>>1
                    if vcfMap_obj[currentChrom][mid][0]<posInPop1:
                        low=mid+1
                    elif vcfMap_obj[currentChrom][mid][0]>posInPop1:
                        high=mid-1
                    else:
#                         if AltInPop1 == vcfMap_obj[currentChrom][mid][2]:#same alt alle
                        if vcfMap_obj_idx!=len(vcfMaplist)-1:
                            elementToAppend.append(vcfMap_obj[currentChrom][mid][1:])
                        elif vcfMap_obj_idx==len(vcfMaplist)-1:
                            elementToAppend.append(vcfMap_obj[currentChrom][mid][1:])
                            multipleVcfMap[currentChrom].append(elementToAppend)

                        break
                else:
                    if jointmode=="i":
                        #ignore the rec
                        break

    return copy.deepcopy(multipleVcfMap)
def bedfiletools(bedfilename, withtitle=True):
    """
        return m={chr1:[(startpos,endpos,[optional_fields]),(),,,],chr2:[],,,,,}
    """
    m = {}
    f = open(bedfilename, 'r')
    if withtitle:
        f.readline()
    for line in f:
        linelist = re.split(r"\s+", line.strip())
        if len(linelist) < 3:
            continue
        if linelist[0] in m:
            m[linelist[0].strip()].append((int(linelist[1]), int(linelist[2]), linelist[3:]))
        else:
            m[linelist[0].strip()] = [(int(linelist[1]), int(linelist[2]), linelist[3:])]
    f.close()
    for chrom in m.keys():
        m[chrom].sort(key=lambda listRec:listRec[0])
    return m
def interval_setOperation(bedlikefile1, bedlikefile2):
    """
    note : no overlap region within bedlikefile1 ,bedlikefile2 each are required
    return 
    """
    aaa = open("aaalldone.txt", 'w')
    intersectionRegions = {}
    diffRegions = {}
    unionRegions = {}
    intervals1 = bedfiletools(bedlikefile1)
    intervals2 = bedfiletools(bedlikefile2)
    for chrom in intervals1:
        if chrom in intervals2:
            # collect regions before the first region of intervals
            if intervals2[chrom][0][0] < intervals1[chrom][0][0]:
                for i in range(len(intervals2[chrom])):
                    if intervals2[chrom][i][0] > intervals1[chrom][0][0] and intervals2[chrom][i - 1][0] < intervals1[chrom][0][0]:  # if i-1=-1 still in consider
                        break
                    elif intervals2[chrom][i][1] < intervals1[chrom][0][0]:  # ==will consider blow
                        unionRegions = collectRegion(unionRegions, chrom, (intervals2[chrom][i][0], intervals2[chrom][i][1]))
            else:
                print("intervals2[chrom][0][0]>=intervals1[chrom][0][0]")
            # start collect the other region  
            q_idx = 0
            while q_idx < len(intervals1[chrom]):
                q5 = intervals1[chrom][q_idx][0]
                q3 = intervals1[chrom][q_idx][1]
#             for q5,q3,optionfieldslist in intervals1[chrom]:
                low = 0
                high = len(intervals2[chrom]) - 1
                while low <= high:
                    mid = (low + high) >> 1
                    if intervals2[chrom][mid][0] < q5:
                        low = mid + 1
                    elif intervals2[chrom][mid][0] > q5:
                        high = mid - 1
                    else:
                        print(chrom, file=open("testout.txt", 'a'))
                        midcount=0
                        if intervals2[chrom][mid][1] < q3:
                            intersectionRegions=collectRegion(intersectionRegions, chrom,(intervals2[chrom][mid][0], intervals2[chrom][mid][1]))
                            mid += 1;midcount+=1
                            while mid < len(intervals2[chrom]) and intervals2[chrom][mid][0] <= intervals1[chrom][q_idx][1]:
                                
                                if intervals1[chrom][q_idx][1] >= intervals2[chrom][mid][1]:
                                    intersectionRegions = collectRegion(intersectionRegions, chrom, (intervals2[chrom][mid][0], intervals2[chrom][mid][1]))
#                                     unionRegions=collectRegion(unionRegions,chrom,(q5,intervals1[chrom][q_idx][1]))
                                    diffRegions = collectRegion(diffRegions, chrom, (intervals2[chrom][mid-1][1], intervals2[chrom][mid][0]))
                                    
                                    mid += 1
                                    midcount+=1
                                    continue
                                else:
                                    diffRegions = collectRegion(diffRegions, chrom, (intervals2[chrom][mid-1][1], intervals2[chrom][mid][0]))
                                    intersectionRegions = collectRegion(intersectionRegions, chrom, (intervals2[chrom][mid][0], intervals1[chrom][q_idx][1]))
                                    unionRegions = collectRegion(unionRegions, chrom, (intervals1[chrom][q_idx][0], intervals2[chrom][mid][1]))
                                    break
                            else:

                                diffRegions=collectRegion(diffRegions, chrom,(intervals2[chrom][mid-1][1],intervals1[chrom][q_idx][1]))

                                unionRegions = collectRegion(unionRegions, chrom, (intervals1[chrom][q_idx][0], intervals1[chrom][q_idx][1]))
                        else:  # intervals2[chrom][mid][1]>=q3
                            intersectionRegions = collectRegion(intersectionRegions,chrom,(intervals1[chrom][q_idx][0], intervals1[chrom][q_idx][1]))
                            q_idx+=1
                            while q_idx < len(intervals1[chrom]) and intervals1[chrom][q_idx][0] <= intervals2[chrom][mid][1]:
                                if intervals1[chrom][q_idx][1] >= intervals2[chrom][mid][1]:
                                    unionRegions = collectRegion(unionRegions, chrom, (q5, intervals1[chrom][q_idx][1]))
                                    diffRegions = collectRegion(diffRegions, chrom, (intervals2[chrom][mid][1], intervals1[chrom][q_idx][1]))
                                    intersectionRegions = collectRegion(intersectionRegions, chrom, (intervals1[chrom][q_idx][0], intervals2[chrom][mid][1]))
                                    break
                                else:
                                    intersectionRegions = collectRegion(intersectionRegions, chrom, (intervals1[chrom][q_idx][0], intervals1[chrom][q_idx][1]))
                                    q_idx += 1
#                                     print("continue",str(q_idx),str(intervals1[chrom][q_idx][0]),str(intervals1[chrom][q_idx][1]))
                            else:
                                unionRegions = collectRegion(unionRegions, chrom, (intervals2[chrom][mid][0], intervals2[chrom][mid][1]))
                                q_idx -= 1
                            
                        break
                else:
                    print("high:", str(high), "should < low:", str(low),chrom, file=open("testout.txt", 'a'))
                    high3 = intervals2[chrom][high][1];high5 = intervals2[chrom][high][0]
                    if high < 0:
                        high3 = -1;high5 = -2
                    lowcount = 0
                    if low >= len(intervals2[chrom]):
                        if intervals1[chrom][q_idx][0] < high3:
                            if intervals1[chrom][q_idx][1] < high3:
                                intersectionRegions = collectRegion(intersectionRegions, chrom, (intervals1[chrom][q_idx][0], intervals1[chrom][q_idx][1]))
                                
                            else:
                                intersectionRegions = collectRegion(intersectionRegions, chrom, (intervals1[chrom][q_idx][0], high3))
                                diffRegions = collectRegion(diffRegions, chrom, (high3, intervals1[chrom][q_idx][1]))
                        else:
                            diffRegions = collectRegion(diffRegions, chrom, (intervals1[chrom][q_idx][1], intervals1[chrom][q_idx][1]))
                    else:    
#                         while q_idx <len(intervals1[chrom]):
                        print(str(q_idx), str(low), str(high3), str(len(intervals1[chrom])), str(len(intervals2[chrom])), file=open("testout.txt", 'a'))
                        q5 = intervals1[chrom][q_idx][0]
                        q3 = intervals1[chrom][q_idx][1]
                        low5 = intervals2[chrom][low][0];low3 = intervals2[chrom][low][1]
                        if q5 >= high3:
                            while low < len(intervals2[chrom]):
                                low5 = intervals2[chrom][low][0];low3 = intervals2[chrom][low][1]
                                print(low3, low5)
                                if q3 > low5:
                                    if q3 > low3:
                                        intersectionRegions = collectRegion(intersectionRegions, chrom, (low5, low3))
                                        if lowcount == 0:
                                            diffRegions = collectRegion(diffRegions, chrom, (q5, low5))
                                        else:
                                            print("====================", file=open("testout.txt", 'a'))
                                            print(intervals2, file=open("testout.txt", 'a'))
                                            print(str(low), str(len(intervals2[chrom])), chrom, str(lowcount), file=open("testout.txt", 'a'))
                                            diffRegions = collectRegion(diffRegions, chrom, (intervals2[chrom][low - 1][1], low5))
                                        low += 1
                                        lowcount += 1
                                        print("q3>low5,q3>low3,low+=1", str(q5), str(q3), str(low5), file=aaa)
                                        continue
                                    else:  # out condition2
                                        intersectionRegions = collectRegion(intersectionRegions, chrom, (low5, q3))
                                        if lowcount == 0:
                                            diffRegions = collectRegion(diffRegions, chrom, (q5, low5))
                                        else:
                                            diffRegions = collectRegion(diffRegions, chrom, (intervals2[chrom][low - 1][1], low5))
                                        if high3 == q5:
                                            unionRegions = collectRegion(unionRegions, chrom, (high5, low3))
                                        else:
                                            unionRegions = collectRegion(unionRegions, chrom, (q5, low3))
                                        break
                                else:  # out condition1
                                    if lowcount == 0:
                                        diffRegions = collectRegion(diffRegions, chrom, (q5, q3))
                                    else:
                                        diffRegions = collectRegion(diffRegions, chrom, (intervals2[chrom][low - 1][1], q3))
                                    if high3 == q5:
                                        if q3 == low5:
                                            unionRegions = collectRegion(unionRegions, chrom, (high5, low3))
                                        else:
                                            unionRegions = collectRegion(unionRegions, chrom, (high5, q3))
                                    elif q3 == low5:
                                        unionRegions = collectRegion(unionRegions, chrom, (q5, low3))
                                    break
                            else:
                                # out condition1 ,but not intervel on the right side of q
                                if lowcount == 0:
                                    diffRegions=collectRegion(diffRegions, chrom, (q5, q3))
                                else:
                                    diffRegions = collectRegion(diffRegions, chrom, (intervals2[chrom][low - 1][1], q3))
                            lowcount = 0
#                             break
                        else:
                            while q3 < high3:
                                intersectionRegions = collectRegion(intersectionRegions, chrom, (q5, q3))
                                q_idx += 1
                                print("q3<high3", "q_idx+=1", str(q_idx), str(q5), str(q3), str(low5), file=aaa)
                                if q_idx >= len(intervals1[chrom]):
                                    unionRegions = collectRegion(unionRegions, chrom, (high5, high3))
#                                     nextchrom=True
                                    break
                                q5 = intervals1[chrom][q_idx][0]
                                q3 = intervals1[chrom][q_idx][1]
                                if q5 >= high3:
                                    q_idx -= 1
                                    break
                            else:                          
#                             if q3>high3:
                                while low < len(intervals2[chrom]):
                                    low5 = intervals2[chrom][low][0];low3 = intervals2[chrom][low][1]
                                    if q3 > low5:
                                        if q3 > low3:
                                            if lowcount == 0:
                                                intersectionRegions = collectRegion(intersectionRegions, chrom, (q5, high3))
                                                diffRegions=collectRegion(diffRegions,chrom,(high3,low5))
                                            else:
                                                intersectionRegions = collectRegion(intersectionRegions, chrom, (low5, low3))
                                                diffRegions = collectRegion(diffRegions, chrom, (intervals2[chrom][low - 1][1], low5))
                                            lowcount += 1
                                            print("q3>low5,q3>low3,low+=1", str(q5), str(q3), str(low5), file=aaa)
                                            low += 1
                                            continue
                                        else:  # out condition 4
                                            if lowcount == 0:
                                                print("condition4", str(lowcount))
                                                intersectionRegions = collectRegion(intersectionRegions, chrom, (q5, high3))
                                                intersectionRegions = collectRegion(intersectionRegions, chrom, (low5, q3))
                                                diffRegions = collectRegion(diffRegions, chrom, (high3, low5))
                                            else:
                                                intersectionRegions = collectRegion(intersectionRegions, chrom, (low5, q3))
                                                diffRegions = collectRegion(diffRegions, chrom, (intervals2[chrom][low - 1][1], low5))
                                            unionRegions = collectRegion(unionRegions, chrom, (high5, low3))
                                            break
                                    else:  # out condition 3
                                        intersectionRegions = collectRegion(intersectionRegions, chrom, (q5, high3))
                                        if lowcount==0:
                                            diffRegions = collectRegion(diffRegions, chrom, (high3, q3))
                                        else:
                                            diffRegions=collectRegion(diffRegions,chrom,(intervals2[chrom][low - 1][1],q3))
                                        if q3 == low5:
                                            unionRegions = collectRegion(unionRegions, chrom, (high5, low3))
                                        else:
                                            unionRegions = collectRegion(unionRegions, chrom, (high5, q3))
                                        break
                                else:
                                    if lowcount==0:
                                        diffRegions=collectRegion(diffRegions,chrom,(high3,q3))
                                    else:
                                        diffRegions=collectRegion(diffRegions,chrom,(intervals2[chrom][low - 1][1],q3))
                                lowcount = 0
                q_idx += 1
                print("overall", "q_idx+=1", str(q_idx), file=aaa)
    aaa.close()
    d_f = open("testd.txt", 'w')
    i_f = open("testi.txt", 'w')
    u_f = open("testu.txt", 'w')
    for chrom in diffRegions:
        for start, end in diffRegions[chrom]:
            print(chrom, str(start), str(end), file=d_f)
    for chrom in intersectionRegions:
        for start, end in intersectionRegions[chrom]:
            print(chrom, str(start), str(end), file=i_f)
    for chrom in unionRegions:
        for start, end in unionRegions[chrom]:
            print(chrom, str(start), str(end), file=u_f)
    d_f.close()
    i_f.close()
    u_f.close()
    return intersectionRegions, diffRegions
                                
def collectRegion(maplist, mapkey, maplistvalue):
    if maplistvalue[0] == maplistvalue[1]:
        return maplist
    if mapkey in maplist:
        maplist[mapkey].append(maplistvalue)
    else:
        maplist[mapkey] = [maplistvalue]
    returnmaplist = copy.deepcopy(maplist)
    return returnmaplist




def generateIndexByChrom(refFastaFileName, indexFileName, mapname=None,startchar=">",chrsignal=None):
    refFastaFile = open(refFastaFileName, 'r')
    refChromIndex = {}
    refline = refFastaFile.readline()
    while refline:
        if re.search(r'^['+startchar+']', refline) != None:
#            collist = re.split(r'\s+', refline)
            if mapname == "transcript:":
                currentChromNo = re.search(r'transcript:(.*?)\s+', refline).group(1).strip()
            else:
                if not chrsignal:
                    a = re.search(r'^'+startchar+'([^'+startchar+'|]+)', (re.split(r'\s+', refline))[0]).group(1).lower()
                else:
                    linelist=re.split(r'\s+', refline)
                    a=re.sub('[’!"#$%&\'()*+,-./:;<=>?@[\\]^_`{|}~]+',"", linelist[linelist.index(chrsignal)+1]).lower()
                currentChromNo=a.replace("chr","")
                currentChromNo=transform_roman_num2_alabo(currentChromNo)
                print(currentChromNo,type(currentChromNo))
            refChromIndex[currentChromNo] = int(refFastaFile.tell())  # from here is the sequence
        refline = refFastaFile.readline()
    pickle.dump(refChromIndex, open(indexFileName, 'wb'))
    refFastaFile.close()
def transform_roman_num2_alabo(one_str,changesignal=True):  
    ''''' 
    将罗马数字转化为阿拉伯数字 
    '''
    if re.search('^M{0,4}(CM|CD|D?C{0,3})(XC|XL|L?X{0,3})(IX|IV|V?I{0,3})$',one_str)!=None and changesignal:  
        define_dict={'I':1,'V':5,'X':10,'L':50,'C':100,'D':500,'M':1000}  
        if one_str=='0':  
            return 0  
        else:  
            res=0  
            for i in range(0,len(one_str)):  
                if i==0 or define_dict[one_str[i]]<=define_dict[one_str[i-1]]:  
                    res+=define_dict[one_str[i]]  
                else:  
                    res+=define_dict[one_str[i]]-2*define_dict[one_str[i-1]]  
            return str(res)
    else:
        return one_str
def generateFasterRefIndex(refFastaFileName, indexFileName,startchar=">",chrsignal=None,romanSignal=False):#,mapname=None
    refFastaFile = open(refFastaFileName, 'r')
    refChromIndex = {}
    refline = refFastaFile.readline()
    while refline:
        if re.search(r'^['+startchar+']', refline) != None:
            basecount=1
            m=1
            """
            chrsignal 可以代替原来 mapname == "transcript:"的功能，即 令  chrsignal= "transcript:" 
            """
#             if mapname == "transcript:":
#                 currentChromNo = re.search(r'transcript:(.*?)\s+', refline).group(1).strip()
#             else:
            if not chrsignal or chrsignal not in refline:
                a = re.search(r'^'+startchar+'([^'+startchar+'|]+)', (re.split(r'\s+', refline))[0]).group(1)
            else:
                linelist=re.split(r'\s+', refline)
                a=refline[refline.index(chrsignal)+len(chrsignal):].split()[0].strip('[’!"#$%&\'()*+,-./:;<=>?@[\\]^_`{|}~]+').strip()
                #a=re.sub('[’!"#$%&\'()*+,-./:;<=>?@[\\]^_`{|}~]+',"", refline[refline.index(chrsignal)+len(chrsignal):].split()[0])#for example chromosome 1,
            a=transform_roman_num2_alabo(a,romanSignal)
            currentChromNo=a.replace("chr","").replace("CHR", "")
            print(currentChromNo,type(currentChromNo))
            refChromIndex[currentChromNo] = [(basecount,int(refFastaFile.tell()))]# (no of base befor,cur file pos)
        else:
            basecount+=len(refline.strip())
            if basecount>=6000*m:
                refChromIndex[currentChromNo].append((basecount,int(refFastaFile.tell())))
                m+=1
        refline = refFastaFile.readline()
    pickle.dump(refChromIndex,open(indexFileName, 'wb'))
    refFastaFile.close()
def loadAnchorFile(anchorFile):
    anchorDATASTRUCTURE={}
    """
    {chr1:[(53353,53806,scaffold451,558997,558537,-),(57200,62371,scaffold451,553669,548504,-),(),,],chr2:[],,,,}
    """
    reverseAnchorDATASTRUCTURE={}
    """
    {scaffold451:{chr1:[0,1,2,,,,],chr2:[],,,},C17734302:{chr1:[idx,,,],chr2:[idx,,,],,,}}  idx is idx in the list of anchorDATASTRUCTURE[chr1] 
    """
    newanchorfilehandler=open(anchorFile,'r')
    print("chrNo\tstartpos\tendpos\tstrand",file=open("sexchromrecs.bed",'w'))
    for line in newanchorfilehandler:
        linelist=re.split(r"\s+",line.strip())
        if linelist[0].strip() in anchorDATASTRUCTURE:
            anchorDATASTRUCTURE[linelist[0].strip()].append((int(linelist[1].strip()),int(linelist[2].strip()),linelist[3].strip(),int(linelist[4].strip()),int(linelist[5].strip()),linelist[6].strip()))
        else:
            anchorDATASTRUCTURE[linelist[0].strip()]=[(int(linelist[1].strip()),int(linelist[2].strip()),linelist[3].strip(),int(linelist[4].strip()),int(linelist[5].strip()),linelist[6].strip())]
        #fill reverseAnchorDATASTRUCTURE
        if linelist[3].strip() in reverseAnchorDATASTRUCTURE:
            if linelist[0].strip() in reverseAnchorDATASTRUCTURE[linelist[3].strip()]:
                reverseAnchorDATASTRUCTURE[linelist[3].strip()][linelist[0].strip()].append(len(anchorDATASTRUCTURE[linelist[0].strip()])-1)
            else:
                reverseAnchorDATASTRUCTURE[linelist[3].strip()]={linelist[0].strip():[len(anchorDATASTRUCTURE[linelist[0].strip()])-1]}
        else:
            reverseAnchorDATASTRUCTURE[linelist[3].strip()]={linelist[0].strip():[len(anchorDATASTRUCTURE[linelist[0].strip()])-1]}
    for scaffold in reverseAnchorDATASTRUCTURE.keys():
        if "Z" in reverseAnchorDATASTRUCTURE[scaffold] or "z" in reverseAnchorDATASTRUCTURE[scaffold] or "W" in reverseAnchorDATASTRUCTURE[scaffold] or "w" in reverseAnchorDATASTRUCTURE[scaffold] or "X" in reverseAnchorDATASTRUCTURE[scaffold] or "x" in reverseAnchorDATASTRUCTURE[scaffold] or "Y" in reverseAnchorDATASTRUCTURE[scaffold] or "y" in reverseAnchorDATASTRUCTURE[scaffold]:
            for chrom,recs in reverseAnchorDATASTRUCTURE[scaffold].items():
                if "Z" ==chrom.upper() or "W" ==chrom.upper() or  "X" ==chrom.upper() or "Y" ==chrom.upper():
                    for idx in recs:
                        print("\t".join([str(e) for e in anchorDATASTRUCTURE[chrom][idx][2:]]) if anchorDATASTRUCTURE[chrom][idx][-1]=="+" else "\t".join((anchorDATASTRUCTURE[chrom][idx][2],str(anchorDATASTRUCTURE[chrom][idx][4]),str(anchorDATASTRUCTURE[chrom][idx][3]),anchorDATASTRUCTURE[chrom][idx][5])),file=open("sexchromrecs.bed",'a'))
    newanchorfilehandler.close()
    return anchorDATASTRUCTURE,reverseAnchorDATASTRUCTURE
def generateIndexByChromForFQ(refFastaFileName, indexFileName, mapname=None,startchar=">"):
    refFastaFile = open(refFastaFileName, 'r')
    refChromIndex = {}
    refline = refFastaFile.readline()
    while refline:
        if re.search(r'^['+startchar+']', refline) != None:
#            collist = re.split(r'\s+', refline)
            if mapname == "transcript:":
                currentChromNo = re.search(r'transcript:(.*?)\s+', refline).group(1).strip()
            else:
                currentChromNo = re.search(r'[^'+startchar+']+', (re.split(r'\s+', refline))[0]).group(0)
            refChromIndex[currentChromNo] = int(refFastaFile.tell())  # from here is the sequence
        if re.search(r"^[+]\s*$",refline)!=None:
            refChromIndex[currentChromNo]=(refChromIndex[currentChromNo],int(refFastaFile.tell()))
        refline = refFastaFile.readline()
    pickle.dump(refChromIndex, open(indexFileName, 'wb'))
    refFastaFile.close()

def getRefSeqBypos_faster(refFastahandle, fasterrefindex, currentChromNO, startpos, endpos, currentChromNOlen=None, seektuple=()):
    '''
    pos start at 1
    seektuple=(filepos,basesbeforefilepos)
    the refSeqMap has only one chromosome's sequence
    There is no restriction on refFastahander
    refindex must indexed by generateFasterRefIndex func
    return {'1': [0, 'A']}
    '''    
    betaTestNocheck_currentChromNOlen=True#i.e. read all until the next >
    refSeqMap = {}
    if startpos <= 0:
        startpos = 1
    print("getRefSeqBypos_faster currentChromNO:", currentChromNO, startpos, endpos)
    if currentChromNOlen != None and endpos > currentChromNOlen:
        endpos = currentChromNOlen

    filehandle = refFastahandle
    if not seektuple or seektuple[1] > startpos:
        refSeqMap[currentChromNO] = [startpos - 1]
        low=0;high=len(fasterrefindex[currentChromNO])-1
        while low<=high:
            mid=(low + high)>>1
            if fasterrefindex[currentChromNO][mid][0]<startpos:
                low=mid+1
            elif fasterrefindex[currentChromNO][mid][0]>startpos:
                high=mid-1
            else:
                perivouspos_idx=mid
#                 filehandle.seek(refindex[currentChromNO][mid][1])
                break
        else:
            if fasterrefindex[currentChromNO][high][0]>startpos:
                print("notice getRefSeqBypos_faster:",low,high,mid)
                high-=1
                perivouspos_idx=high
            else:
                perivouspos_idx=high
        filehandle.seek(fasterrefindex[currentChromNO][perivouspos_idx][1])
          # seekmap is empty so go to the first bases of the currentChromNO
        if fasterrefindex[currentChromNO][perivouspos_idx][0]<startpos:
            preseq = filehandle.read(startpos- fasterrefindex[currentChromNO][perivouspos_idx][0])
            dn = preseq.count('\n')
            while dn != 0:
                preseq = filehandle.read(dn)
                dn = preseq.count('\n')
        elif fasterrefindex[currentChromNO][perivouspos_idx][0]==startpos:
            pass
        else:
            print(fasterrefindex[currentChromNO][perivouspos_idx],startpos)
            print("getRefSeqBypos_faster ERROR error")
            
        # now filehander is right stay at the startpos
        """
        try
        """
        myseqline = filehandle.read(endpos - startpos + 1)
        """
        except :read to big space
        """
        myseqn = myseqline.count('\n')
#        if len(myseqline)>200:
#            print(myseqn)
#            exit(-1)
#        print("myseqline=",myseqline,"myseqn", myseqn)
        while myseqn != 0:  # fill the same number of \n with bases
            myseqline = myseqline.replace('\n', '')
            myseqline += filehandle.read(myseqn)
            myseqn = myseqline.count('\n')
            
#            print(currentChromNO,myseqline, myseqn)
        if myseqline.count('>') >= 1:
            
            print(currentChromNO, myseqline.index('>'),myseqline[myseqline.index('>')-10:myseqline.index('>')+10], myseqn)
            print("may be need chrlength")
            myseqline=myseqline[:myseqline.index('>')]
            if not betaTestNocheck_currentChromNOlen:
                exit(-1)
        refSeqMap[currentChromNO].extend(list(myseqline))
    else:
        filehandle.seek(seektuple[0])  # seekmap is not empty
        refSeqMap[currentChromNO] = [startpos - 1]
        preseq = filehandle.read(startpos - seektuple[1] - 1)
        dn = preseq.count('\n')
        while dn != 0:
            preseq = filehandle.read(dn)
            dn = preseq.count('\n')
        # now filehander is right stay at the startpos
        myseqline = filehandle.read(endpos - startpos + 1)
        myseqn = myseqline.count('\n')
        while myseqn != 0:  # fill the same number of \n with bases
            myseqline = myseqline.replace('\n', '')
            myseqline += filehandle.read(myseqn)
            myseqn = myseqline.count('\n')
        refSeqMap[currentChromNO].extend(list(myseqline))
    plus = myseqline.count('>')
    if plus != 0:
        print("getRefSeqBypos", currentChromNO, startpos, endpos)
        return -1
    
    return refSeqMap 
def getRefSeqBypos(refFastahandle, refindex, currentChromNO, startpos, endpos, currentChromNOlen=None, seektuple=()):
    '''
    pos start at 1
    seektuple=(filepos,basesbeforefilepos)
    the refSeqMap has only one chromosome's sequence
    There is no restriction on refFastahander
    '''    
    refSeqMap = {}
    if startpos <= 0:
        startpos = 1
    print("getRefSeqBypos", currentChromNO, startpos, endpos)
    if currentChromNOlen != None and endpos > currentChromNOlen:
        endpos = currentChromNOlen

    filehandle = refFastahandle
    if not seektuple or seektuple[1] > startpos:
        refSeqMap[currentChromNO] = [startpos - 1]
        filehandle.seek(refindex[currentChromNO])  # seekmap is empty so go to the first bases of the currentChromNO
        preseq = filehandle.read(startpos - 1)
        dn = preseq.count('\n')
        while dn != 0:
            preseq = filehandle.read(dn)
            dn = preseq.count('\n')
            
        # now filehander is right stay at the startpos
        myseqline = filehandle.read(endpos - startpos + 1)
        myseqn = myseqline.count('\n')
#        if len(myseqline)>200:
#            print(myseqn)
#            exit(-1)
#        print("myseqline=",myseqline,"myseqn", myseqn)
        while myseqn != 0:  # fill the same number of \n with bases
            myseqline = myseqline.replace('\n', '')
            myseqline += filehandle.read(myseqn)
            myseqn = myseqline.count('\n')
            
#            print(currentChromNO,myseqline, myseqn)
        if myseqline.count('>') >= 1:
            print(currentChromNO, myseqline, myseqn)
            exit(-1)
        refSeqMap[currentChromNO].extend(list(myseqline))
    else:
        filehandle.seek(seektuple[0])  # seekmap is not empty
        refSeqMap[currentChromNO] = [startpos - 1]
        preseq = filehandle.read(startpos - seektuple[1] - 1)
        dn = preseq.count('\n')
        while dn != 0:
            preseq = filehandle.read(dn)
            dn = preseq.count('\n')
        # now filehander is right stay at the startpos
        myseqline = filehandle.read(endpos - startpos + 1)
        myseqn = myseqline.count('\n')
        while myseqn != 0:  # fill the same number of \n with bases
            myseqline = myseqline.replace('\n', '')
            myseqline += filehandle.read(myseqn)
            myseqn = myseqline.count('\n')
        refSeqMap[currentChromNO].extend(list(myseqline))
    plus = myseqline.count('>')
    if plus != 0:
        print("getRefSeqBypos", currentChromNO, startpos, endpos)
        return -1
    
    return refSeqMap 
def getRefSeqMap(refFastafilehander, currentChromNO=None, preBaseTotal=0, linesOnce=500000, mapname=None):
    '''
    the refSeqMap has only one chromosome's sequence
    '''
    refSeqMap = {}
    print(refFastafilehander.tell(),currentChromNO)
    if currentChromNO == None:
        refline = refFastafilehander.readline() 
        print("getRefSeqMap", refline)
        if mapname == "transcript:":
            currentChromNO = re.search(r'transcript:(.*?)\s+', refline).group(1).strip()
        else:
            currentChromNO = re.search(r'[^>]+', (re.split(r'\s+', refline))[0]).group(0)
        refSeqMap[currentChromNO] = [preBaseTotal]  # preBaseTotal=0
        print("getRefSeqMap", currentChromNO)
    elif currentChromNO == "end of the reffile":
        return refSeqMap, currentChromNO, "end of the reffile"
    else:
        refSeqMap[currentChromNO] = [preBaseTotal]
#     for refline in refFastafilehander:
    
    while 1:
        refline = refFastafilehander.readline()
        if not refline:
            return refSeqMap, currentChromNO, "end of the reffile"
        if re.search(r'^[>]', refline) != None:
            collist = re.split(r'\s+', refline)
            print("getRefSeqMap","3", re.search(r'[^>]+', collist[0]).group(0))
#            refSeqMap[currentChromNO] = [0]
            nextChromNo = re.search(r'[^>]+', collist[0]).group(0)
            return refSeqMap, currentChromNO, nextChromNo  # clean the refSeqMap and report the current chromNO
        else:
            refSeqMap[currentChromNO].extend(list(refline.strip().lower()))
        linesOnce -= 1    
        if linesOnce == 0:
            break                
    else:
        return refSeqMap, currentChromNO, "end of the reffile"
    return refSeqMap, currentChromNO, currentChromNO
def getRefSeqBypos_fromFQ(refFastqhandle, refindex, currentChromNO, startpos, endpos, currentChromNOlen=None):
    '''
    pos start at 1
    seektuple=(filepos,basesbeforefilepos)
    the refSeqMap has only one chromosome's sequence
    There is no restriction on refFastahander
    '''    
    refSeqMap = {}
    if startpos <= 0:
        startpos = 1
    print("getRefSeqBypos", currentChromNO, startpos, endpos)
    if currentChromNOlen != None and endpos > currentChromNOlen:
        endpos = currentChromNOlen

    filehandle = refFastqhandle
    if True:
        refSeqMap[currentChromNO] = [startpos - 1]
        filehandle.seek(refindex[currentChromNO][0])  # seekmap is empty so go to the first bases of the currentChromNO
        preseq = filehandle.read(startpos - 1)
        dn = preseq.count('\n')
        while dn != 0:
            preseq = filehandle.read(dn)
            dn = preseq.count('\n')
            
        # now filehander is right stay at the startpos
        myseqline = filehandle.read(endpos - startpos + 1)
        myseqn = myseqline.count('\n')
#        if len(myseqline)>200:
#            print(myseqn)
#            exit(-1)
#        print("myseqline=",myseqline,"myseqn", myseqn)
        while myseqn != 0:  # fill the same number of \n with bases
            myseqline = myseqline.replace('\n', '')
            myseqline += filehandle.read(myseqn)
            myseqn = myseqline.count('\n')
            
#            print(currentChromNO,myseqline, myseqn)
        if myseqline.count('>') >= 1:
            print(currentChromNO, myseqline, myseqn)
            exit(-1)
        refSeqMap[currentChromNO].extend(list(myseqline))

    plus = myseqline.count('>')
    if plus != 0:
        print("getRefSeqBypos", currentChromNO, startpos, endpos)
        return -1
    refFastqhandle
    return refSeqMap   
#phylip format
class PhylipError(Exception):
    pass
def nunique_lengths(seq_of_seq):
    """
    Given a sequence of sequences, return the number of unique lengths.
    @param: a sequence of sequences
    @return: the number of unique sequence lengths
    """
    return len(set(len(seq) for seq in seq_of_seq))
     
def get_lines(f):
    """
    @param raw_lines: raw lines
    这从网上抄的程序，就离谱，改了好多bug 
    @return: a list of nonempty lines
    """
    lines=[]
    #process header line
    lines.append(f.readline().rstrip())#'\r\n'
    line=f.readline()
    #line[:10]+ #line[10:]
    line_len_first=1#'\r\n'
    #'\r\n'
    while line.split():
        lines.append(re.split(r"\s+",line)[0]+ re.search(r"\s+",line).group(0)+"".join(re.split(r"\s+",line.rstrip())[1:]).replace(" ","").replace("\t",""))#'\r\n'#line[:10]+line[10:].rstrip('\r\n').replace(" ","").replace("\t","")
        line=f.readline().strip();line_len_first+=1
        
    line_len=1    
    for line in f:
        if not line.split() and line_len_first==line_len:
            line_len=1
        else:
            lines[line_len]+=line.rstrip().replace(" ","").replace("\t","")#'\r\n'
            line_len+=1            
#     lines = [x.rstrip('\r\n') for x in raw_lines]
    return [x for x in lines if x]
def decode_phyliplines(raw_lines):
    """
    This parses lines of a non-interleaved phylip sequence file.
    @param raw_lines: raw lines of a non-interleaved phylip alignment file
    @return: headers, sequences
    """
    lines = get_lines(raw_lines)
    header_line, data_lines = lines[0], lines[1:]
    header_row = header_line.split()
    if len(header_row) != 2:
        raise PhylipError('the header should be a line with two integers')
    ntaxa_s, ncolumns_s = header_row
    try:
        ntaxa = int(ntaxa_s)
        ncolumns = int(ncolumns_s)
    except ValueError:
        raise PhylipError('the header should be a line with two integers')
    # check the number of data lines
    ntaxa_observed = len(data_lines)
    if ntaxa_observed != ntaxa:
        msg_a = 'the header says there are %d taxa' % ntaxa
        msg_b = 'but %d taxa were observed' % ntaxa_observed
        raise PhylipError(msg_a + msg_b)
    # all line lengths should be the same
    if nunique_lengths(data_lines) != 1:
        raise PhylipError('all data lines should be the same length')
    # break lines into taxa and data
    compound_data_rows = [[x[:10].strip(), x[10:].strip()] for x in data_lines]
#     compound_data_rows = [[re.split(r"\s+",x.strip())[0],re.split(r"\s+",x.strip())[1:]] for x in data_lines]
    headers, sequences = zip(*compound_data_rows)
    ncolumns_observed = len(sequences[0])
    if ncolumns_observed != ncolumns:
        msg_a = 'the header says there are %d alignment columns' % ncolumns
        msg_b = 'but %d alignment columns were observed' % ncolumns_observed
        raise PhylipError(msg_a + msg_b)
    maptoreturn={}
    for i in range(ntaxa_observed):
        maptoreturn[headers[i]]=sequences[i]
    return maptoreturn
def encode_phyliplines(headers, sequences,maxlen=10):
    """
    This creates the contents of a non-interleaved phylip sequence file.
    @param headers: some header strings
    @param sequences: some sequence strings
    """
    nrows = len(headers)
    ncols = len(sequences[0])
    out_lines = ['%d %d' % (nrows, ncols)]
    for h, seq in zip(headers, sequences):
        out_h = h[:maxlen].ljust(maxlen)
        out_lines.append(out_h + seq)
    return '\n'.join(out_lines)
#phylip format 

def mapFormatPrint(mapdata,filehander=sys.stdout):
    for m,k in mapdata:
        print(*m,":",*k,sep=":",file=filehander)
    
        

class GATK_depthfile():
    onecopy=None
    static_depthfileName=None
    static_allrecsforcurchrom_mapbypos=None
    def __init__(self, depthfileName, indexFileName,ismultplethreads=False):
        super().__init__()
        self.covfileidx = {}
        self.title = []
        self.depthfileName = depthfileName
        self.ismultplethreads=ismultplethreads
        
        if self.static_depthfileName==None:
            self.onecopy=True
            self.static_depthfileName=self.depthfileName
        elif self.onecopy and self.static_depthfileName==self.depthfileName:
            self.onecopy=True
        else:
            self.onecopy=False
        try:
            self.covfileidx = pickle.load(open(indexFileName, 'rb'))
        except IOError:
            GATK_depthfile.indexGATK_depthfile(depthfileName, indexFileName)
            self.covfileidx = pickle.load(open(indexFileName, 'rb'))
        self.title = self.covfileidx["title"]
        self.chromOrder=self.covfileidx["chromOrder"]
        self.depthfilefp = open(depthfileName, 'r')
        self.depthfilefp.readline()
        self.allrecsforcurchrom_mapbypos=None
        self.curchrom=None
    @staticmethod
    def indexGATK_depthfile(depthfileName, indexFileName):
        """
        {chrom:position_in_file_of_first_genomepos_of_this_chrom,chrom:position,,,,,,}
        """
        depthfile = open(depthfileName, 'r')
        covfileidx = {}
      
        currentChrom = None
        lastPosition = 0
        line = depthfile.readline()
        linelist = re.split(r"\s+", line.strip())
#        self.title = linelist
        print("title", line, linelist)
        covfileidx["title"] = linelist
        covfileidx["chromOrder"]=[]
        lastPosition = depthfile.tell()
        line = depthfile.readline()
        linelist = re.split(r"\s+", line)        
        while line:      
            linelist = re.split(r"\s+", line)
            if currentChrom != re.search(r"^([\w\W]*)[:]([\d]*)", linelist[0]).group(1):
                currentChrom = re.search(r"^([\w\W]*)[:]([\d]*)", linelist[0]).group(1)
                covfileidx["chromOrder"].append(currentChrom)
                covfileidx[currentChrom] = lastPosition
            lastPosition = depthfile.tell()
    
            line = depthfile.readline()
        pickle.dump(covfileidx, open(indexFileName, 'wb'))
        depthfile.close()
    def set_depthfilefp(self, targetchrom, targetloc, lastposoffilehandler=0):
        """
        set the self.depthfilefp to the line in the file where chrom==targetchrom pos==targetloc-1
        """
        
        if targetloc == 1:
            self.depthfilefp.seek(self.covfileidx[targetchrom])
            return "found"
        
        searchfp = open(self.depthfileName, 'r')
        searchfp.seek(lastposoffilehandler)
        linelist = re.split(r"\s+", searchfp.readline())
#        print("set_depthfilefp:",targetchrom,"search for targetloc:",targetloc,self.covfileidx[targetchrom],lastposoffilehandler,linelist)
        currentChrom = re.search(r"^([\w\W]*)[:]([\d]*)", linelist[0]).group(1)
        pos = int(re.search(r"^([\w\W]*)[:]([\d]*)", linelist[0]).group(2))
        if currentChrom == targetchrom and pos == targetloc:
            self.depthfilefp.seek(lastposoffilehandler)
            searchfp.close()
            return "found"
        if currentChrom != targetchrom or pos > targetloc:
            searchfp.seek(self.covfileidx[targetchrom])
            targetfpposition = self.covfileidx[targetchrom]
        else:
            targetfpposition = searchfp.tell()
        linelist = re.split(r"\s+", searchfp.readline())
        currentChrom = re.search(r"^([\w\W]*)[:]([\d]*)", linelist[0]).group(1)
        pos = int(re.search(r"^([\w\W]*)[:]([\d]*)", linelist[0]).group(2))
        # set the filehandler locate at the nearest location to the target location
        while currentChrom == targetchrom:
            if pos == targetloc:
                self.depthfilefp.seek(targetfpposition)
                searchfp.close()
                return "found"
            else:
                targetfpposition = searchfp.tell()
                linelist = re.split(r"\s+", self.depthfilefp.readline())
                currentChrom = re.search(r"^([\w\W]*)[:]([\d]*)", linelist[0]).group(1)
                pos = int(re.search(r"^([\w\W]*)[:]([\d]*)", linelist[0]).group(2))
        else:
            searchfp.close()
            return "didn't find"         
    def getnextposline(self):
        line = self.depthfilefp.readline()
        print("getnextposline", line)
        linelist = re.split(r"\s+", line)
        chrom = re.search(r"^([\w\W]*)[:]([\d]*)", linelist[0]).group(1)
        pos = int(re.search(r"^([\w\W]*)[:]([\d]*)", linelist[0]).group(2))
        return chrom, pos, linelist, self.depthfilefp.tell()
    def getdepthByPos_optimized(self, targetchr, targetloc):
        if self.curchrom!=targetchr:
            if self.ismultplethreads or (self.onecopy and self.static_allrecsforcurchrom_mapbypos==None) or not self.onecopy:#((first time only one obj )or not first time )and not multiplethreads. multimple copy
                self.depthfilefp.seek(self.covfileidx[targetchr])
                content=self.depthfilefp.read(self.covfileidx[self.chromOrder[self.chromOrder.index(targetchr.strip()) + 1]] - self.covfileidx[targetchr.strip()])
                contentlines=re.split(r"\n",content.strip())
                self.allrecsforcurchrom_mapbypos={}
                for line in contentlines:
                    chr_line=re.split(r":",line.strip())
                    linelist=re.split(r"\s+",chr_line[1])
                    self.allrecsforcurchrom_mapbypos[int(linelist[0])]=linelist
                self.curchrom=targetchr
                self.static_allrecsforcurchrom_mapbypos=self.allrecsforcurchrom_mapbypos
            elif self.onecopy and self.static_allrecsforcurchrom_mapbypos!=None:#only one copy
                self.allrecsforcurchrom_mapbypos=self.static_allrecsforcurchrom_mapbypos# 
    
                 
        if targetloc in self.allrecsforcurchrom_mapbypos:
#             print("getdepthByPos_optimized",self.allrecsforcurchrom_mapbypos[targetloc])
            return self.allrecsforcurchrom_mapbypos[targetloc]
        else:
            return ["0"]*len(self.title)
    def getdepthByPos(self, targetchr, targetloc, lastposoffilehandler=0):
         
        linelist = re.split(r"\s+", self.depthfilefp.readline())
        posoffilehandlerofnextchr=self.covfileidx[self.chromOrder[self.chromOrder.index(targetchr)+1]]
        if linelist[0] == "":  # read the last line of the depthfile
            self.depthfilefp.seek(self.covfileidx[targetchr])
        else:
            chrom = re.search(r"^([\w\W]*)[:]([\d]*)", linelist[0]).group(1)
            pos = int(re.search(r"^([\w\W]*)[:]([\d]*)", linelist[0]).group(2))
            if chrom == targetchr and pos == targetloc:
                return linelist
            if chrom != targetchr or pos > targetloc:
                self.depthfilefp.seek(self.covfileidx[targetchr])
            if chrom == targetchr and pos < targetloc and self.depthfilefp.tell()<lastposoffilehandler and lastposoffilehandler<posoffilehandlerofnextchr:
                self.depthfilefp.seek(lastposoffilehandler)
                pass  # use the lastposoffilehandler to set the filehanlder quickly
        linelist = re.split(r"\s+", self.depthfilefp.readline())
        chrom = re.search(r"^([\w\W]*)[:]([\d]*)", linelist[0]).group(1)
        pos = int(re.search(r"^([\w\W]*)[:]([\d]*)", linelist[0]).group(2))
        while chrom == targetchr:
            if pos == targetloc:
                return linelist
            linelist = re.split(r"\s+", self.depthfilefp.readline())
            chrom = re.search(r"^([\w\W]*)[:]([\d]*)", linelist[0]).group(1)
            pos = int(re.search(r"^([\w\W]*)[:]([\d]*)", linelist[0]).group(2))
        else:
            return ["0"]*len(self.title)
        
                    
    def closedepthfile(self):
        self.title.clear()
        self.covfileidx.clear()
        self.depthfilefp.close()
                
class FastQ_Util():
    def __init__(self):
        super().__init__()
    @staticmethod
    def generateIndexByChrom(FastQFileName, indexFileName):
        """
            to the consenus that produce by vcfutils.pl vcf2fq
            every line that start with one '@' and length of the line less than 20 is indexed
        """
        fasqfile = open(FastQFileName, 'r')
        refChromIndex = {}
        fqline = fasqfile.readline()
        while fqline:
            collist = re.split(r'\s+', fqline)
            if re.search(r'^[@][^@]+$', collist[0]) != None:                
                if len(collist[0]) > 20:  # may be fqline is located in the quality value block
                    fqline = fasqfile.readline()
                    continue
#                print(collist[0],fqline)
                currentChromNo = re.search(r'^[@]([^@]+)$', collist[0]).group(1).strip()
#                print(currentChromNo)
                refChromIndex[currentChromNo] = int(fasqfile.tell())  # from here is the sequence
            fqline = fasqfile.readline()
        pickle.dump(refChromIndex, open(indexFileName, 'wb'))
        fasqfile.close()
        
        
    @staticmethod
    def getConsenusSeqMap(fastQFileName, dbtools , tablename="chromosome", primaryID="chrID", bp_start=None, bp_end=None):
        '''
        the refSeqMap has only one chromosome's sequence
        '''
#        print(dbtools,fastQFileName,"inside FastQ_Util")
        fqfile = open(fastQFileName, 'r')
        sql = "select * from " + tablename
        seqMapByChrom = {}
        try:
            ChromIndexMap = pickle.load(open(fastQFileName + ".myindex", 'rb'))
        except IOError:
            FastQ_Util.generateIndexByChrom(fastQFileName, fastQFileName + ".myindex")
            ChromIndexMap = pickle.load(open(fastQFileName + ".myindex", 'rb'))
        
        totalChroms = dbtools.operateDB("select", "select count(*) from " + tablename)[0][0]
        
        print(totalChroms)
        currentchrID = dbtools.operateDB("select", sql + " limit 0,1")[0][0]
        seqMapByChrom[currentchrID] = ""
        for i in range(0, totalChroms, 20):
            currentsql = sql + " order by " + primaryID + " limit " + str(i) + ",20"
            result = dbtools.operateDB("select", currentsql)
            for row in result:
                currentchrID = row[0]
                if currentchrID in ChromIndexMap:
                    seqMapByChrom[currentchrID] = ""
                    fqfile.seek(ChromIndexMap[currentchrID])
                    line = fqfile.readline()
                    while line.strip() != "+":
                        seqMapByChrom[currentchrID] += line.strip()
#                        print(line.strip())
                        line = fqfile.readline()
        return seqMapByChrom

    
    
class Window():
    def __init__(self):
        super().__init__()
        self.winValueL = []  # [(startPos,lastPos,value),(),,,,,,]
    def forPhastConsFormat(self, L, L_End_Pos, windowWidth, Caculator, winStart=0):
        """
        without overlap
        L=[startpos,endpos,value]
        """
        self.winValueL = []
        currentIdx = 0
        while currentIdx != len(L):
            if L[currentIdx][0] >= winStart and L[currentIdx][1] <= (winStart + windowWidth):
#                 print(L[currentIdx][1] - L[currentIdx][0])
                Caculator.process(L[currentIdx], L[currentIdx][1] - L[currentIdx][0])
                if L[currentIdx][1] == (winStart + windowWidth):
                    value = Caculator.getResult()
                    self.winValueL.append((winStart, winStart + windowWidth, value))
                    winStart += windowWidth
                
            elif L[currentIdx][0] > winStart and L[currentIdx][0] < (winStart + windowWidth) and L[currentIdx][1] > (winStart + windowWidth):
                print("2")
                frontPartPosNum = winStart + windowWidth - L[currentIdx][0]
                rearPartPosNum = L[currentIdx][1] - (winStart + windowWidth)
                Caculator.process(L[currentIdx], frontPartPosNum)
                value = Caculator.getResult()
                self.winValueL.append((winStart, winStart + windowWidth, value))
                winStart += windowWidth
                Caculator.process(L[currentIdx], rearPartPosNum)
            elif L[currentIdx][0] <= winStart and L[currentIdx][1] > winStart and L[currentIdx][1] < (winStart + windowWidth):
                print("3")
                rearPartPosNum = L[currentIdx][1] - (winStart + windowWidth)
                Caculator.process(L[currentIdx], frontPartPosNum)
            elif (winStart + windowWidth) <= L[currentIdx][0]:
                print("4")
                while  (winStart + windowWidth) <= L[currentIdx][0]:
                    print("Util", winStart + windowWidth , L[currentIdx][0])
                    self.winValueL.append((winStart, winStart + windowWidth, Caculator.getResult()))
                    winStart += windowWidth
            elif L[currentIdx][1] == winStart:
                self.winValueL.append((winStart, winStart + windowWidth, Caculator.getResult()))
                winStart += windowWidth
            currentIdx += 1
        else:
            self.winValueL.append((winStart, winStart + windowWidth, Caculator.getResult()))
                
    def slidWindowOverlap(self, L, L_End_Pos, windowWidth, slideSize, Caculator,L_Start_Pos=0):
        print("L_End_Pos",L_End_Pos,L_Start_Pos)
        """
        require at least one SNP pos in L in sliding window
        window slide from L_Start_Pos to L_End_Pos
        L = [(pos,p1,p2,p3,A_base_idx),(pos,"a,b","c,d","e,f",0),(pos,"a,b","c,d","e,f",1),....] for D-statistics wihtout "no covered"
        or 
        L = [(pos, REF, ALT, INFO,FORMAT,sampleslist),(pos, REF, ALT, INFO,FORMAT,sampleslist),(),...........] for any score need one vcf,
        or 
        L = [(pos,REF,ALT,(INFO,FORMAT,sampleslist),(INFO,FORMAT,sampleslist)),(pos,REF,ALT,(INFO,FORMAT,sampleslist),(INFO,FORMAT,sampleslist)),(),...........] for any score need two or more vcf,for example two vcf's compare,eg. fst, one or multiple vcf caculate hp
        like the two situation upside,return a value
        or
        L = [(pos,samples1dp,samples2dp,samples3dp,,,),(pos,samples1dp,samples2dp,samples3dp,,,),(),(),......]#in this situation ,value formation like this ([sample1_pecentage,sample2_pecentage,,,],[sample1_average_depth,sample2_average_depth,,,])
        
        """
#         del self.winValueL[:]
        self.winValueL = []  # notice here
        nextIdx = -1  # always be -1 if windowWidth == slideSize
        currentIdx = 0
        winStart = L_Start_Pos
        FoundNextIdx = False
        firstComeInWin = True
        notjustforsnp = True
        for findfirstidx in range(len(L)):
            if L[findfirstidx][0]>winStart:
                currentIdx=findfirstidx
                break
        while currentIdx != len(L) and L[currentIdx][0]<= L_End_Pos:
#             print(L[currentIdx][0],L[currentIdx])
            if L[currentIdx][0] > winStart and  L[currentIdx][0] <= (winStart + windowWidth):
#                if notjustforsnp or (len(L[currentIdx][1])==1 and re.search(r'[^a-zA-Z]', L[currentIdx][2]) != None and len(L[currentIdx][2])==1):# it's not a snp? indel or cnv
                if firstComeInWin:
                    startPos = L[currentIdx][0]
                    firstComeInWin = False
                lastPos = L[currentIdx][0]
                Caculator.process(L[currentIdx])
                if FoundNextIdx == False and L[currentIdx][0] > (winStart + slideSize):  # always go to |currentIdx+=1|
                    nextIdx = currentIdx
                    FoundNextIdx = True
            else:
                noofsnps, value = Caculator.getResult()
                try:
                    self.winValueL.append((startPos, lastPos, noofsnps, value))
#                     print(startPos, lastPos, noofsnps, value)
                except:
                    print("no snp in first severl wins", len(L), currentIdx, value, L[currentIdx])
                    self.winValueL.append((0, 0, noofsnps, value))
                    winStart += slideSize
                    continue
#                 self.winValueL.append((0, 0, value))
                winStart += slideSize
                firstComeInWin = True
                
                FoundNextIdx = False
                if nextIdx == -1:
                    if slideSize >= windowWidth:
                        while not (L[currentIdx][0] > winStart and  L[currentIdx][0] <= (winStart + windowWidth)) and L[currentIdx][0] > winStart + windowWidth:
                            winStart += slideSize
                            noofsnps, value = Caculator.getResult()
                            self.winValueL.append((0, 0, noofsnps, value))
                        if L[currentIdx][0] < winStart:
                            while currentIdx != len(L):
                                if L[currentIdx][0] > winStart and L[currentIdx][0] <= (winStart + windowWidth):
                                    break
                                elif L[currentIdx][0] < winStart:
                                    winStart += slideSize
                                    noofsnps, value = Caculator.getResult()
                                    self.winValueL.append((0, 0, noofsnps, value))
                                currentIdx += 1
#                             self.winValueL.append((0,0,'NA'))
#                             winStart += slideSize
                    continue  # go to |if L[currentIdx][0] > winStart and L[currentIdx][0] < (winStart + windowWidth):| in upside block
                else:
                    currentIdx = nextIdx
                    nextIdx = -1
                    continue
                
            currentIdx += 1
        else:
            noofsnps, value = Caculator.getResult()
            try :
                self.winValueL.append((startPos, lastPos, noofsnps, value))
#                 print(startPos, lastPos, noofsnps, value)
            except UnboundLocalError:
                self.winValueL.append((0, 0, noofsnps, value))
#             if nextIdx!=-1:
#                 currentIdx = nextIdx
#                 nextIdx = -1
#                 while currentIdx != len(L):
#                     lastPos = L[currentIdx][0]
#                     Caculator.process(L[currentIdx])
#                     currentIdx += 1
#                 else:
#                     noofsnps, value = Caculator.getResult()
#                     try:
#                         self.winValueL.append((startPos, lastPos, noofsnps, value))
#                     except:
#                         self.winValueL.append((0, 0, noofsnps, value))
#            
        
        n = int((L_End_Pos-L_Start_Pos - (len(self.winValueL) * slideSize + windowWidth)) / slideSize) + 1
        for i in range(n):
            noofsnps, value = Caculator.getResult()
            self.winValueL.append((0, 0, noofsnps, value))
def distributionfuncdraft(intervalFileName,dataFileNames,col_to_bined1,col_to_bined2=0,col_to_mean=None):
    col_to_bined1-=1;col_to_bined2-=1;col_to_mean-=1
    intervalfile=open(intervalFileName,'r')
    intervalMap_count={}
    intervalMap_mean={}
    intervalMap_sum={}
    for line in intervalfile:
        linelist=re.split(r'\s+',line.strip())
        intervalMap_sum[float(linelist[0]),float(linelist[1])]=0
        intervalMap_mean[float(linelist[0]),float(linelist[1])]=[]
        intervalMap_count[float(linelist[0]),float(linelist[1])]=0
    intervalfile.close()
    for df in dataFileNames:
        print(df)
        datafile=open(df,'r')
        print(datafile.readline())
        for line in datafile:
            linelist=re.split(r'\s+',line.strip())
            if col_to_bined2!=0:
                value_to_bin=float(linelist[col_to_bined2])-float(linelist[col_to_bined1])
            else:
                value_to_bin=float(linelist[col_to_bined1])
            if col_to_mean!=None:
                value_to_mean=float(linelist[col_to_mean])
                
            for a,b in sorted(intervalMap_count.keys()):
                if value_to_bin>=a and value_to_bin<b:
                    if col_to_mean!=None:
                        intervalMap_mean[a,b].append(value_to_mean)
                    intervalMap_count[a,b]+=1
        #intervalMap_mean has two use one record value list one record mean value
        if col_to_mean!= None:
            for a,b in intervalMap_mean.keys():
                if len(intervalMap_mean[a,b])!=0:
                    intervalMap_sum[a,b]+=sum(intervalMap_mean[a,b])
                    intervalMap_mean[a,b]=[]
        datafile.close()
    if col_to_mean!= None:
        for a,b in intervalMap_count.keys():
            if intervalMap_count[a,b]==0:
                intervalMap_mean[a,b]="NA"
            else:
                intervalMap_mean[a,b]=intervalMap_sum[a,b]/intervalMap_count[a,b]
    return copy.deepcopy(intervalMap_count),copy.deepcopy(intervalMap_mean)


class BinDepth():
    def __init__(self, depthbinFileName):
        """
        chr:[(win0),(win1),(firstpos,endpos,sample1,sample2,,,),(40001,50000,passed,passed,passed,filtered,),,,]
        """
        self.speciesname, self.depthbinmap = BinDepth.readDepthfileintoaMap(depthbinFileName)
    @staticmethod
    def readDepthfileintoaMap(depthbinFileName):
        depthbinmap = {}
        depthfile = open(depthbinFileName, 'r')
        titlelist = re.split(r"\s+", depthfile.readline().strip())
        title = titlelist[3:]
        print("class BinDepth", titlelist)
        linelist = re.split(r'\s+', depthfile.readline().strip())
        depthbinmap[linelist[0].strip()] = [tuple(linelist[2:])]
        for line in depthfile:
            linelist = re.split(r'\s+', line.strip())
            if linelist[0] in depthbinmap:
                depthbinmap[linelist[0]].append(tuple(linelist[2:]))
            else:
                depthbinmap[linelist[0]] = [tuple(linelist[2:])]
        depthfile.close()
        return title, depthbinmap
        
     
        
        
        
        
        
