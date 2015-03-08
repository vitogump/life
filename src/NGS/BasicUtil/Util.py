# -*- coding: UTF-8 -*-
import copy
import random
import string,numpy
import re
import pickle
import os
import src.NGS.BasicUtil.DBManager as dbm
# from src.NGS.BasicUtil import *

'''
Created on 2013-6-30

@author: rui
'''
def bedfiletools(bedfilename, withtitle=False):
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
            print(linelist, file=open("testbedfile", 'a'))
            m[linelist[0].strip()] = [(int(linelist[1]), int(linelist[2]), linelist[3:])]
    f.close()
    for chrom in m.keys():
        m[chrom].sort(key=lambda listRec:listRec[0])
        print(chrom, file=open("testbedfile", 'a'))
        for i in m[chrom]:
            print(i, file=open("testbedfile", 'a'))
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
                        if intervals2[chrom][mid][1] < q3:
                            while mid < len(intervals2[chrom]) and intervals2[chrom][mid][0] <= intervals1[chrom][q_idx][1]:
                                
                                if intervals1[chrom][q_idx][1] >= intervals2[chrom][mid][1]:
#                                     unionRegions=collectRegion(unionRegions,chrom,(q5,intervals1[chrom][q_idx][1]))
                                    diffRegions = collectRegion(diffRegions, chrom, (intervals2[chrom][mid][1], intervals1[chrom][q_idx][1]))
                                    intersectionRegions = collectRegion(intersectionRegions, chrom, (intervals2[chrom][mid][0], intervals2[chrom][mid][1]))
                                    mid += 1
                                    continue
                                else:
                                    intersectionRegions = collectRegion(intersectionRegions, chrom, (intervals2[chrom][mid][0], intervals1[chrom][q_idx][1]))
                                    unionRegions = collectRegion(unionRegions, chrom, (intervals1[chrom][q_idx][0], intervals2[chrom][mid][1]))
                                    break
                            else:
                                unionRegions = collectRegion(unionRegions, chrom, (intervals1[chrom][q_idx][0], intervals1[chrom][q_idx][1]))
                        else:  # intervals2[chrom][mid][1]>=q3
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
                    print("high:", str(high), "should < low:", str(low), file=open("testout.txt", 'a'))
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
                                        diffRegions = collectRegion(diffRegions, chrom, (high3, q3))
                                        if q3 == low5:
                                            unionRegions = collectRegion(unionRegions, chrom, (high5, low3))
                                        else:
                                            unionRegions = collectRegion(unionRegions, chrom, (high5, q3))
                                        break
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
def complementary(seq):
    """
    ['tg', 'a', 't', 'g', 'c', 'acacacgatg', 'ctttttcccccccc', 'c', 'c', 'a', 'a', 'aaagagagagacagaaaaaggc', 'atatcgactg', 'catcga']
    reverse to 
    ['ac', 't', 'a', 'c', 'g', 'tgtgtgctac', 'gaaaaagggggggg', 'g', 'g', 't', 't', 'tttctctctctgtctttttccg', 'tatagctgac', 'gtagct']
    """
    newseq = []
    for i in range(0, len(seq)):
        if seq[i].lower() == 'a':
            newseq.insert(i, 't')
        elif seq[i].lower() == 't':
            newseq.insert(i, 'a')
        elif seq[i].lower() == 'c':
            newseq.insert(i, 'g')
        elif seq[i].lower() == 'g':
            newseq.insert(i, 'c')
        elif len(seq[i]) > 1:
            newseq.insert(i, complementary(seq[i]))
        else:
            newseq.insert(i, seq[i])
    if isinstance(seq, str):
        newseq = "".join(newseq)
    return newseq

def random_str(randomlength=8):
    a = list(string.ascii_letters)
    random.shuffle(a)
    return ''.join(a[:randomlength])

def generateIndexByChrom(refFastaFileName, indexFileName, mapname=None):
    refFastaFile = open(refFastaFileName, 'r')
    refChromIndex = {}
    refline = refFastaFile.readline()
    while refline:
        if re.search(r'^[>]', refline) != None:
#            collist = re.split(r'\s+', refline)
            if mapname == "transcript:":
                currentChromNo = re.search(r'transcript:(.*?)\s+', refline).group(1).strip()
            else:
                currentChromNo = re.search(r'[^>]+', (re.split(r'\s+', refline))[0]).group(0)
            refChromIndex[currentChromNo] = int(refFastaFile.tell())  # from here is the sequence
        refline = refFastaFile.readline()
    pickle.dump(refChromIndex, open(indexFileName, 'wb'))
    refFastaFile.close()

def getGtfMap(gtfFileName, elementTypes=["CDS", "stop_codon"]):
    """protein_codingMap={chromNo:[[transcript_id,strand,start,end,(feature, elemStart, elemEnd, frame),(),(),,,,,],
                        [transcript_id,strand,start,end,(),(),(),,,],[],,,,,,,],
               chromNo:[],,,,,,,,,,,,,,}
        chrtranscrpitididxMap{chromNo:{transcript_id:ttanscript_id_idx,transcript_id:ttanscript_id_idx,,,,,},
                                chromNo:{},chromNo:{},,,,}
    """
    gtfFileHandler = open(gtfFileName, 'r')
    protein_codingMap = {}
    chrtranscrpitididxMap = {}
#     gtfline = gtfFileHandler.readline()
    jumpout = False
    for getfirstcds in gtfFileHandler:
        if re.search(r"^#", getfirstcds) != None:
            print(getfirstcds)
            continue
        gtfColList = re.split(r'\s+', getfirstcds)
        chromNo = gtfColList[0].strip()
        protein_codingMap[chromNo] = []
        if "transcript_id" in gtfColList:
            transcript_id_idx = gtfColList.index("transcript_id") + 1
            gene_id = gtfColList.index("gene_id")
        else:
            continue
        print("transcript_id_idx", transcript_id_idx)
        transcript_id = re.search(r'\"(.*)\";', gtfColList[transcript_id_idx].strip()).group(1)
        countInChrom = 0        
        for elementType in elementTypes:
            if elementType == gtfColList[2].strip():
                jumpout = True
                protein_codingMap[chromNo] = [[transcript_id, gtfColList[6], int(gtfColList[3]), int(gtfColList[4]), (gtfColList[2], int(gtfColList[3]), int(gtfColList[4]), gtfColList[7])]]
        else:
            print(getfirstcds)
        if jumpout:
            break
    
    chrtranscrpitididxMap[chromNo] = {transcript_id:0}
    for gtfline in gtfFileHandler:
        gtfColList = re.split(r'\s+', gtfline)
        transcript_id = re.search(r'\"(.*)\";', gtfColList[transcript_id_idx].strip()).group(1)
        for elementType in elementTypes:
            if elementType == gtfColList[2].strip():
                break
        else:
            continue
        chromNo = gtfColList[0].strip()
        if chromNo in protein_codingMap:
            if transcript_id in chrtranscrpitididxMap[chromNo].keys():
                tanscript_id_idx = chrtranscrpitididxMap[chromNo][transcript_id]
                protein_codingMap[chromNo][tanscript_id_idx].append((gtfColList[2], int(gtfColList[3]), int(gtfColList[4]), gtfColList[7]))
                protein_codingMap[chromNo][tanscript_id_idx][2] = min(protein_codingMap[chromNo][tanscript_id_idx][2], int(gtfColList[3]))
                protein_codingMap[chromNo][tanscript_id_idx][3] = max(protein_codingMap[chromNo][tanscript_id_idx][3], int(gtfColList[4]))
            else:
                protein_codingMap[chromNo].append([transcript_id, gtfColList[6], int(gtfColList[3]), int(gtfColList[4]), (gtfColList[2], int(gtfColList[3]), int(gtfColList[4]), gtfColList[7])])
                chrtranscrpitididxMap[chromNo][transcript_id] = len(protein_codingMap[chromNo]) - 1
        else:
             protein_codingMap[chromNo] = [[transcript_id, gtfColList[6], int(gtfColList[3]), int(gtfColList[4]), (gtfColList[2], int(gtfColList[3]), int(gtfColList[4]), gtfColList[7])]]
             chrtranscrpitididxMap[chromNo] = {transcript_id:0}
    else:
        pass                 
    gtffilepath = re.search(r"^.*[/]", gtfFileName).group(0)
    testfile = open(gtffilepath + "protein_codingMap.sort.txt", 'w')    
    for chromNo in protein_codingMap.keys():
        protein_codingMap[chromNo].sort(key=lambda listRec:listRec[2])
        # 先按照转录本起始坐标排序，下面是对转录本内元件排序，不过是什么排序方法忘记了，仔细读一下吧
        for j in range(len(protein_codingMap[chromNo])):
            for t4_indx in range(5, len(protein_codingMap[chromNo][j])):
                t4_key = protein_codingMap[chromNo][j][t4_indx]
                t4_indxp = t4_indx - 1
                while t4_indxp >= 4 and protein_codingMap[chromNo][j][t4_indxp][1] > t4_key[1]:
                    protein_codingMap[chromNo][j][t4_indxp + 1] = protein_codingMap[chromNo][j][t4_indxp]
                    t4_indxp = t4_indxp - 1
                else:
                    protein_codingMap[chromNo][j][t4_indxp + 1] = t4_key
#            print(protein_codingMap[chromNo][j],protein_codingMap[chromNo][j][t4_indxp][1] > t4_key[1],t4_indxp)
        print(chromNo, "num of transcrpit:", len(protein_codingMap[chromNo]), file=testfile)
        for i in range(len(protein_codingMap[chromNo])):
#            print(protein_codingMap[chromNo][i][0],protein_codingMap[chromNo][i][1],protein_codingMap[chromNo][i][2],protein_codingMap[chromNo][i][3])
            for k in range(len(protein_codingMap[chromNo][i])):
                print(protein_codingMap[chromNo][i][k], file=testfile)
    testfile.close()
    gtfFileHandler.close()
    return protein_codingMap
def getNearestGenegroup(gtfList, pos):
    """
    input:for a chrom,contain all transcript of this chrom
    gtfList    =    [[transcript_id,strand,start,end,(feature, elemStart, elemEnd, frame),(),(),,,,,],
                        [transcript_id,strand,start,end,(),(),(),,,],[],,,,,,,]
    return: the first gene that after the pos and the genes contain in or overlap with or contact with this gene indirect
    geneOverlapList=[[transcript_id,strand,start,end,(feature, elemStart, elemEnd, frame),(),(),,,,,],
                        [transcript_id,strand,start,end,(),(),(),,,],[],,,,,,,]
    order by "start"
    """
    
    if gtfList == None:
        return []
    for i in range(len(gtfList)):
        if gtfList[i][2] > pos: 
            geneOverlapList = [gtfList[i]]
            break
    else:
        if pos > gtfList[-1][3]:
            print("what's wrong")
            return []
        else:
            exit(-1)

    furthest = gtfList[i][3]
    i += 1
    while len(gtfList) > i and furthest >= gtfList[i][2]:
        if gtfList[i][0] != geneOverlapList[-1][0]:  # this judgement maybe not need
            geneOverlapList.append(gtfList[i])
        furthest = max(furthest, gtfList[i][3])
        i += 1
    return furthest, geneOverlapList
def getSNPrecInCDS(codon_idx_in_tscptSeqAllCds, lenOftscpt, codon, codon_m, cds_frame_ONETRSCPT, gtftrscpt_list):
    for i in range(3):
        if codon[i] != codon_m[i]:
            snp_idx_in_tscptALLCds = i + codon_idx_in_tscptSeqAllCds
            ref_base = codon[i].upper()
            alt_base = codon_m[i].upper()
    if gtftrscpt_list[1] == "+":
        for cds_idx in sorted(cds_frame_ONETRSCPT.keys()):
            if cds_frame_ONETRSCPT[cds_idx][1] > snp_idx_in_tscptALLCds:
                cds_idx -= 1
                break
        snp_pos = snp_idx_in_tscptALLCds - cds_frame_ONETRSCPT[cds_idx][1] + gtftrscpt_list[cds_idx][1]
    elif gtftrscpt_list[1] == "-":
        for cds_idx in sorted(cds_frame_ONETRSCPT.keys()):
            if cds_frame_ONETRSCPT[cds_idx][1] > lenOftscpt - snp_idx_in_tscptALLCds - 1:
                cds_idx -= 1
                break
        snp_pos = lenOftscpt - snp_idx_in_tscptALLCds - 1 - cds_frame_ONETRSCPT[cds_idx][1] + gtftrscpt_list[cds_idx][1]
        ref_base = complementary(ref_base).upper()
        alt_base = complementary(alt_base).upper()
    return snp_pos, ref_base, alt_base
def getGeneGrouplist(gtfList):
    """
        geneGrouplist=[genegroup1,genegroup2,genegroup3,....]
                     =[[furthest,geneOverlapLists1],[furthest,geneOverlapLists2],[furthest,geneOverlapLists3],.....]
                     =[[furthest,[],[],[],[],....],[furthest,[],[],[],[],....],[furthest,[],[],[],[],....],.....]
                     =[[furthest,[transcript_id,strand,start,end,(),(),(),...],[transcript_id,strand,start,end,(),(),(),...],....],[furthest,[transcript_id,strand,start,end,(),(),(),...],[transcript_id,strand,start,end,(),(),(),...],....],[furthest,[],[],...],[],........]
    """
    newgtfList = copy.deepcopy(gtfList)
    geneGrouplist = []
    curpos = newgtfList[0][2] - 1
    while curpos < newgtfList[-1][2]:
        furthest, geneGroup = getNearestGenegroup(newgtfList, curpos)
        curpos = furthest + 1
        geneGroup.insert(0, furthest)
        geneGrouplist.append(geneGroup)
        
    return geneGrouplist
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
            print("getRefSeqMap", re.search(r'[^>]+', collist[0]).group(0))
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
class genes():
    def __init__(self, gtfList, pos, RefSeqList, minintervalbetweengenes_basesperfaline=60):
        super().__init__()
        self.lastgenesRearpos = 0
        self.minintervalbetweengenes_basesperfaline = minintervalbetweengenes_basesperfaline
        self.geneOverlapList = self.getNearestGeneOverlapList(gtfList, pos)
        self.tscptSeqAllCds = {}
        self.cds_frame = {}  # {transcript_id:{cdsidx:(frame,startpos of this cds),cdsidx:(),,,,,}}
        
        
        for gene in self.geneOverlapList:
            
            genename = gene[0]
            self.tscptSeqAllCds[genename] = []
            self.cds_frame[genename] = {}  # {cdsidx:(frame,startpos of this cds),cdsidx:(),,,,,}

            cdsidx = 3
            for feature, elemStart, elemEnd, frame in gene[4:]:
                cdsidx += 1
                if feature == 'CDS':
                    self.cds_frame[genename][cdsidx] = (int(frame), len(self.tscptSeqAllCds[genename]))
                    self.tscptSeqAllCds[genename] += RefSeqList[(elemStart - RefSeqList[0]):(elemEnd - RefSeqList[0] + 1)]  # ???如果不够呢
                elif  feature == "stop_codon":  # feature == 'start_codon' or
                    self.cds_frame[genename][cdsidx] = (int(frame), len(self.tscptSeqAllCds[genename]))
                    self.tscptSeqAllCds[genename] += RefSeqList[(elemStart - RefSeqList[0]):(elemEnd - RefSeqList[0] + 1)]

#         print(genename,elemStart - RefSeqList[0],elemEnd - RefSeqList[0] + 1,elemStart - RefSeqList[0],elemEnd - RefSeqList[0] + 1,RefSeqList[0],len(RefSeqList))
#         print(self.cds_frame)
    def getNearestGeneOverlapList(self, gtfList, pos):
        """
        input:for a chrom,contain all transcript of this chrom
        gtfList=[[transcript_id,strand,start,end,(feature, elemStart, elemEnd, frame),(),(),,,,,],
                            [transcript_id,strand,start,end,(),(),(),,,],[],,,,,,,]
        return: the first gene that after the pos and the genes contain in or overlap with or contact with this gene indirect
        geneOverlapList=[[transcript_id,strand,start,end,(feature, elemStart, elemEnd, frame),(),(),,,,,],
                            [transcript_id,strand,start,end,(),(),(),,,],[],,,,,,,]
        order by "start"
        """
        if gtfList == None:
            return []
        for i in range(len(gtfList)):
            if gtfList[i][2] >= pos and (pos == 1 or self.lastgenesRearpos + self.minintervalbetweengenes_basesperfaline <= gtfList[i][2]):  # the distance between two genes must longer than minintervalbetweengenes_basesperfaline except overlapped genes
                geneOverlapList = [gtfList[i]]
                break
        else:
            if pos > gtfList[-1][3]:
                return []
            else:
                print("getNearestGeneOverlapList", pos, gtfList)
                exit(-1)
 
        furthest = gtfList[i][3]
        i += 1
        while len(gtfList) > i and furthest >= gtfList[i][2]:
            if gtfList[i][0] != geneOverlapList[-1][0]:
                geneOverlapList.append(gtfList[i])
            furthest = max(furthest, gtfList[i][3])
            i += 1
        print("getNearestGeneOverlapList", i, geneOverlapList, pos)
        self.lastgenesRearpos = geneOverlapList[-1][3]
        return geneOverlapList
    def getgeneConsensus(self, RefSeqList, idx_RefSeq, VcfList, idx_vcf, depthfile):
        """
        RefSeqList didn't changed
        """
        CodonTable = {     'ttt': 'F', 'tct': 'S', 'tat': 'Y', 'tgt': 'C',
              'ttc': 'F', 'tcc': 'S', 'tac': 'Y', 'tgc': 'C',
              'tta': 'L', 'tca': 'S', 'taa': '*', 'tga': '*',
              'ttg': 'L', 'tcg': 'S', 'tag': '*', 'tgg': 'W',
              'ctt': 'L', 'cct': 'P', 'cat': 'H', 'cgt': 'R',
              'ctc': 'L', 'ccc': 'P', 'cac': 'H', 'cgc': 'R',
              'cta': 'L', 'cca': 'P', 'caa': 'Q', 'cga': 'R',
              'ctg': 'L', 'ccg': 'P', 'cag': 'Q', 'cgg': 'R',
              'att': 'I', 'act': 'T', 'aat': 'N', 'agt': 'S',
              'atc': 'I', 'acc': 'T', 'aac': 'N', 'agc': 'S',
              'ata': 'I', 'aca': 'T', 'aaa': 'K', 'aga': 'R',
              'atg': 'M', 'acg': 'T', 'aag': 'K', 'agg': 'R',
              'gtt': 'V', 'gct': 'A', 'gat': 'D', 'ggt': 'G',
              'gtc': 'V', 'gcc': 'A', 'gac': 'D', 'ggc': 'G',
              'gta': 'V', 'gca': 'A', 'gaa': 'E', 'gga': 'G',
              'gtg': 'V', 'gcg': 'A', 'gag': 'E', 'ggg': 'G'}
        curpos = RefSeqList[0] + idx_RefSeq
        tscptSeqAllCds_mut = {}
        ref_amino_seq = {}
        mutat_amino_seq = {}
        cns_append = ""
        originallen = {}  # just for test
#        indelatEdgeofCDS=open("indelatEdgeofCDS.txt",'w')
        for gene in self.geneOverlapList:
            genename = gene[0]
#            print(gene, self.cds_frame[genename], sep="\n", file=open("testgeneOverlapList.txt", 'a'))
            
            tscptSeqAllCds_mut[genename] = copy.deepcopy(self.tscptSeqAllCds[genename])
            originallen[genename] = len(tscptSeqAllCds_mut[genename])  # just for test
        
        while idx_vcf != -1 and idx_vcf != len(VcfList) and VcfList[idx_vcf][0] <= self.geneOverlapList[-1][3]:
            vcfpos = VcfList[idx_vcf][0];refalle = VcfList[idx_vcf][1];altalle = VcfList[idx_vcf][2]
            cns_append += ("".join(RefSeqList[idx_RefSeq:idx_RefSeq + (vcfpos - curpos)]))
            idx_RefSeq += (vcfpos - curpos)
            curpos = RefSeqList[0] + idx_RefSeq
            
            if re.search(r'[^a-zA-Z]', altalle) != None:  # contain ',' ie. multiple alle
                cns_append += ("".join(RefSeqList[idx_RefSeq:idx_RefSeq + len(refalle)]))
                idx_RefSeq += len(refalle)
                curpos = RefSeqList[0] + idx_RefSeq
                idx_vcf += 1
                continue
            cns_append += ("".join(RefSeqList[idx_RefSeq:idx_RefSeq + len(altalle)]))
            idx_RefSeq += len(refalle)  # here should still be refalle
            curpos = RefSeqList[0] + idx_RefSeq
            n_refbases = len(refalle);n_altbases = len(altalle)  # situation TAA     TA;     TTA     TTAAACTTCTATACTA;      C       T;    T       TATA;    ACG     A
# first for every variant making cns_append,and then substitute the seq in the cds seq,and finialy translate to protein 
            for gene in self.geneOverlapList:
                genename = gene[0]
                if gene[2] <= vcfpos and gene[3] >= vcfpos:
                    t4_indx = 3
                    for feature, elemStart, elemEnd, frame in gene[4:]:
                        t4_indx += 1
                        if feature == 'CDS' and vcfpos <= elemEnd and vcfpos >= elemStart:
                            if vcfpos + n_refbases - 1 > elemEnd or vcfpos + n_altbases - 1 > elemEnd:
                                print(VcfList[idx_vcf], genename, "indelatEdgeofCDS")
                                break
                            if n_refbases > n_altbases:  # situation TAA     TA;ACG     A
                                print(tscptSeqAllCds_mut[genename][(vcfpos - elemStart + self.cds_frame[genename][t4_indx][1]):(vcfpos - elemStart + self.cds_frame[genename][t4_indx][1] + n_refbases)])
                                tscptSeqAllCds_mut[genename][(vcfpos - elemStart + self.cds_frame[genename][t4_indx][1]):(vcfpos - elemStart + self.cds_frame[genename][t4_indx][1] + n_altbases)] = list(altalle)
                                tscptSeqAllCds_mut[genename][(vcfpos - elemStart + self.cds_frame[genename][t4_indx][1] + n_altbases):(vcfpos - elemStart + self.cds_frame[genename][t4_indx][1] + n_refbases)] = [' '] * (n_refbases - n_altbases)
                                print(tscptSeqAllCds_mut[genename][(vcfpos - elemStart + self.cds_frame[genename][t4_indx][1]):(vcfpos - elemStart + self.cds_frame[genename][t4_indx][1] + n_refbases)])
                            elif n_refbases < n_altbases:  # situation TTA     TTAAACTTCTATACTA;T       TATA;
                                if n_refbases == 1:
                                    tscptSeqAllCds_mut[genename][vcfpos - elemStart + self.cds_frame[genename][t4_indx][1]] = altalle
                                else:
                                    tscptSeqAllCds_mut[genename][(vcfpos - elemStart + self.cds_frame[genename][t4_indx][1]):(vcfpos - elemStart + self.cds_frame[genename][t4_indx][1] + n_refbases - 1)] = list(altalle[0:(n_refbases - 1)])
                                    tscptSeqAllCds_mut[genename][vcfpos - elemStart + self.cds_frame[genename][t4_indx][1] + n_refbases - 1] = altalle[(n_refbases - 1):]
                            else:  # n_refbases==n_altbases==1
                                try:
                                    tscptSeqAllCds_mut[genename][vcfpos - elemStart + self.cds_frame[genename][t4_indx][1]] = altalle
                                except IndexError:
                                    print(self.cds_frame)
                                    print(genename, vcfpos, t4_indx, altalle, elemStart, feature, len(tscptSeqAllCds_mut[genename]))
                                    exit(-1)
       
            idx_vcf += 1
# 该翻译蛋白了吧 还有 看看长度一样不  将最后一个vcf记录之后的序列加入一致序列字符串
        cns_append += "".join(RefSeqList[idx_RefSeq:idx_RefSeq + (self.geneOverlapList[-1][3] - curpos) + 1])
        for gene in self.geneOverlapList:
            
            genename = gene[0]
#            ref_amino_seq[genename] = []
            mutat_amino_seq[genename] = []
#            mutationTypeList=[]
            if originallen[genename] != len(tscptSeqAllCds_mut[genename]):  # just for test
                print(self.tscptSeqAllCds[genename])
                print(tscptSeqAllCds_mut[genename])
                print("Util getgeneConsensus: length of tscptSeqAllCds changed,so there is a indel in the rearpart of the trcptSeq", genename)

            tscptSeqAllCds_mut_str = "".join(filter(lambda e:e.strip() != "", tscptSeqAllCds_mut[genename]))           
            if gene[1] == '+':
                tscptSeqAllCds_mut[genename] = list(tscptSeqAllCds_mut_str)
                for i in range(self.cds_frame[genename][sorted(self.cds_frame[genename].keys())[0]][0], len(tscptSeqAllCds_mut_str), 3):
#                    codon = "".join(self.tscptSeqAllCds[genename][i:i + 3]).lower()
                    codon_m = tscptSeqAllCds_mut_str[i:i + 3].lower()
#                    ref_amino_seq[genename].append(CodonTable[codon])
                    try:
                        mutat_amino_seq[genename].append(CodonTable[codon_m])
                    except KeyError:
                        mutat_amino_seq[genename].append('X')
            else:  # strand == '-'
                tscptSeqAllCds_Revr_Cmplm = complementary(self.tscptSeqAllCds[genename])
                tscptSeqAllCds_Revr_Cmplm.reverse()                
                tscptSeqAllCds_mut_str = list(tscptSeqAllCds_mut_str)
                tscptSeqAllCds_mut_str_Revr_Cmplm = complementary(tscptSeqAllCds_mut_str)
                tscptSeqAllCds_mut_str_Revr_Cmplm.reverse()
                tscptSeqAllCds_mut_str_Revr_Cmplm = "".join(tscptSeqAllCds_mut_str_Revr_Cmplm)
                tscptSeqAllCds_mut[genename] = list(tscptSeqAllCds_mut_str_Revr_Cmplm)
#                 print(genename,sorted(self.cds_frame[genename].keys()),len(tscptSeqAllCds_mut_str_Revr_Cmplm),3)
                for i in range(self.cds_frame[genename][sorted(self.cds_frame[genename].keys())[-1]][0], len(tscptSeqAllCds_mut_str_Revr_Cmplm), 3):
#                    codon = "".join(self.tscptSeqAllCds[genename][i:i+3]).lower()
                    codon_m = tscptSeqAllCds_mut_str_Revr_Cmplm[i:i + 3].lower()
#                    ref_amino_seq[genename].append(CodonTable[codon])
                    try:
                        mutat_amino_seq[genename].append(CodonTable[codon_m])
                    except KeyError:
                        mutat_amino_seq[genename].append('X')
#        testfile = open("animo_acid.txt", 'w')
#        for gene in self.geneOverlapList:#for test only
#            genename = gene[0]
#            print(">" + genename + "\n", file=testfile)
#            print("".join(ref_amino_seq[genename]), file=testfile)
#            print("\n" + "".join(mutat_amino_seq[genename]) + "\n", file=testfile)
#        testfile.close()
        return  tscptSeqAllCds_mut, mutat_amino_seq, cns_append, idx_vcf

        

class GATK_depthfile():
    def __init__(self, depthfileName, indexFileName):
        super().__init__()
        self.covfileidx = {}
        self.title = []
        self.depthfileName = depthfileName
        try:
            self.covfileidx = pickle.load(open(indexFileName, 'rb'))
        except IOError:
            GATK_depthfile.indexGATK_depthfile(depthfileName, indexFileName)
            self.covfileidx = pickle.load(open(indexFileName, 'rb'))
        self.title = self.covfileidx["title"]
        self.depthfilefp = open(depthfileName, 'r')
        self.depthfilefp.readline()
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
        lastPosition = depthfile.tell()
        line = depthfile.readline()
        linelist = re.split(r"\s+", line)        
        while line:      
            linelist = re.split(r"\s+", line)
            if currentChrom != re.search(r"^([\w\W]*)[:]([\d]*)", linelist[0]).group(1):
                currentChrom = re.search(r"^([\w\W]*)[:]([\d]*)", linelist[0]).group(1)
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
    def getdepthByPos(self, targetchr, targetloc, lastposoffilehandler=0):
         
        linelist = re.split(r"\s+", self.depthfilefp.readline())
        if linelist[0] == "":  # read the last line of the depthfile
            self.depthfilefp.seek(self.covfileidx[targetchr])
        else:
            chrom = re.search(r"^([\w\W]*)[:]([\d]*)", linelist[0]).group(1)
            pos = int(re.search(r"^([\w\W]*)[:]([\d]*)", linelist[0]).group(2))
            if chrom == targetchr and pos == targetloc:
                return linelist
            if chrom != targetchr or pos > targetloc:
                self.depthfilefp.seek(self.covfileidx[targetchr])
            elif chrom == targetchr and pos < targetloc - 100:
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
            return []
        
                    
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
                
    def slidWindowOverlap(self, L, L_End_Pos, windowWidth, slideSize, Caculator):
        """
        L = [(pos,p1,p2,p3,A_base_idx),(pos,"a,b","c,d","e,f",0),(pos,"a,b","c,d","e,f",1),....] for D-statistics wihtout "no covered"
        or 
        L = [(pos, REF, ALT, INFO,FORMAT,sampleslist),(pos, REF, ALT, INFO,FORMAT,sampleslist),(),...........] for any score need one vcf,eg.  het
        or 
        L = [(pos,REF,ALT,(INFO,FORMAT,sampleslist),(INFO,FORMAT,sampleslist)),(pos,REF,ALT,(INFO,FORMAT,sampleslist),(INFO,FORMAT,sampleslist)),(),...........] for any score need two vcf's compare,eg. fst
        like the two situation upside,return a value
        or
        L = [(pos,samples1dp,samples2dp,samples3dp,,,),(pos,samples1dp,samples2dp,samples3dp,,,),(),(),......]#in this situation ,value formation like this ([sample1_pecentage,sample2_pecentage,,,],[sample1_average_depth,sample2_average_depth,,,])
        
        """
#         del self.winValueL[:]
        self.winValueL = []  # notice here
        nextIdx = -1  # always be -1 if windowWidth == slideSize
        currentIdx = 0
        winStart = 0
        FoundNextIdx = False
        firstComeInWin = True
        notjustforsnp = True
        while currentIdx != len(L):
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
                except:
                    print("no snp in first win", len(L), currentIdx, value, L[currentIdx])
                    self.winValueL.append((0, 0, 0, value))
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
            except UnboundLocalError:
                self.winValueL.append((0, 0, noofsnps, value))
#             
        
        n = int((L_End_Pos - (len(self.winValueL) * slideSize + windowWidth)) / slideSize) + 1
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
class WinInGenome():           
    def __init__(self, dbname, winFileName6Field, tableName=None):
        super().__init__()
        self.dbname = dbname
        self.chromOrder, self.windbtools, self.wintablewithoutNA, self.wintabletextvalueallwin = self.loadWinDataIntoDB(dbname, winFileName6Field, tableName)
        self.winContainTrscptMap = {}
    def loadWinDataIntoDB(self, dbname, winFileName6Field, tableNamewithoutNA=None):
        chromOrder = []
        if tableNamewithoutNA == None:
            tableNamewithoutNA = random_str()
        tableNametextValueForappendGeneName = tableNamewithoutNA + "textField"
        tempdbtools = dbm.DBTools("10.2.48.140", "root", "1234567", dbname)
        TABLES = {}
        TABLES[tableNamewithoutNA] = (
            "CREATE TABLE " + tableNamewithoutNA + " ("
            " `chrID` varchar(128) NOT NULL ,"
            " `winNo` varchar(128) NOT NULL,"
            " `bp_start` varchar(128) NOT NULL,"
            " `bp_end` varchar(128) NOT NULL,"
            " `snpcount` int(11) NOT NULL,"
            " `winvalue` double NOT NULL,"  #########why?
            " `zvalue` double NOT NULL,"  ##########
            " PRIMARY KEY (`chrID`,`winNo`)"
            ")"
            )
        TABLES[tableNametextValueForappendGeneName] = (
            "CREATE TABLE " + tableNametextValueForappendGeneName + " ("
            " `chrID` varchar(128) NOT NULL ,"
            " `winNo` varchar(128) NOT NULL,"
            " `bp_start` varchar(128) NOT NULL,"
            " `bp_end` varchar(128) NOT NULL,"
            " `snpcount` int(11) NOT NULL,"
            " `winvalue` text NOT NULL,"  ##############why?
            " `zvalue` text NOT NULL,"  ###############
            " PRIMARY KEY (`chrID`,`winNo`)"
            ")"
            )        
        print(TABLES)
        tempdbtools.create_table(TABLES)
        a = os.popen("awk '{print $1}' " + winFileName6Field + "|uniq")
        for chromNo in a:
            chromOrder.append(chromNo.strip())
        a.close()
        a = os.system("awk '$0!~/NA/ && NR!=1{print $0}' " + winFileName6Field + ">" + winFileName6Field + "_tmpfile")
        if a != 0:
            print("awk '$0!~/NA/ && NR!=1{print $0}' " + winFileName6Field + ">" + winFileName6Field + "_tmpfile" + ": failed")
            exit(-1)
        print("awk '$0!~/NA/ && NR!=1{print $0}' " + winFileName6Field + ">" + winFileName6Field + "_tmpfile" + ": ok")
        loaddatasql = "load data local infile '" + winFileName6Field + "_tmpfile' into table " + tableNamewithoutNA + " fields terminated by '\\t'"
        
        shellstatment = "mysql -uroot -p1234567 -D" + dbname.strip() + ' -e "' + loaddatasql + '"'
        
        a = os.system(shellstatment)
        if a != 0:
            print("Util : loadWinDataIntoDB func os.system return not 0")
            exit(-1)
        print(shellstatment + ":ok")
        os.system("rm " + winFileName6Field + "_tmpfile")
        
        loaddatasql = "load data local infile '" + winFileName6Field + "' into table " + tableNametextValueForappendGeneName + " fields terminated by '\\t'"
        
        shellstatment = "mysql -uroot -p1234567 -D" + dbname.strip() + ' -e "' + loaddatasql + '"'
        
        a = os.system(shellstatment)
        
        if a != 0:
            print(shellstatment + ":failed")
            exit(-1)
        print(shellstatment + ":ok")    
        tempdbtools.operateDB("delete", "delete from " + tableNametextValueForappendGeneName + " where chrID='chrNo' and winNo='winNo' and winvalue='winvalue' ")    

        return chromOrder, tempdbtools, tableNamewithoutNA, tableNametextValueForappendGeneName 
    def appendGeneName(self, TranscriptGenetable, genomedbtools, winwidth, slideSize, outfileName):
        outfile = open(outfileName, 'w')
        print("chrNo\twinNo\tfirstsnppos\tlastsnppos\tnoofsnps\twinvalue\tzvalue\tgeneName\ttrscptID", file=outfile)
        totalWins = self.windbtools.operateDB("select", "select count(*) from " + self.wintabletextvalueallwin)[0][0]
        allwins = self.windbtools.operateDB("select", "select * from " + self.wintabletextvalueallwin + " limit 1," + str(totalWins))
        self.windbtools.operateDB("callproc", "mysql_sp_add_column", data=(self.dbname, self.wintabletextvalueallwin, "geneName", "varchar(128)", "default null"))
        self.windbtools.operateDB("callproc", "mysql_sp_add_column", data=(self.dbname, self.wintabletextvalueallwin, "trscptID", "varchar(128)", "default null"))  
        for win in allwins:
            region = (win[0], int(win[1]) * slideSize, int(win[1]) * slideSize + winwidth, win[1], win[5])
            geneNames = ""
            trscptIDs = ""
            for rec in self.collectTrscptInWin(genomedbtools, TranscriptGenetable, None, region):
                print("rec", rec)
                trscptIDs += rec[0].strip() + ";"
                if rec[2].strip() != "":
                    geneNames += (rec[2].strip() + ";")
            self.windbtools.operateDB("update", "update " + self.wintabletextvalueallwin + " set geneName = '" + geneNames[0:-1] + "', trscptID= '" + trscptIDs[0:-1] + "' where chrID= '" + win[0] + "' and winNo=" + win[1])
        allwins = self.windbtools.operateDB("select", "select * from " + self.wintabletextvalueallwin)
        for win in allwins:
#                 win=tuple([a for a in win[0:4]]+['NA','NA']+[a for a in win[6:]])
            if win[-2] == "":
                if win[-1] == "":
                    print(*(win[:-2] + ("NA", "NA")), sep="\t", file=outfile)
                else:
                    print(*(win[:-2] + ("NA", win[-1])), sep="\t", file=outfile)
            else:
                print(*win, sep="\t", file=outfile)
        outfile.close()

    def collectTrscptInWin(self, genomedbtools, trscptableName, vcftable, region, upextend=0, downextend=0):
        """select region overlaped with the trscpt
        reture a list of trscpts [tp_generecord1,tp_generecord2,,,]
        """
        trscptlist = []
        transcripttable = trscptableName
        chrID = region[0]
        Region_start = region[1]
        Region_end = region[2]
        """
        region=(chrom,Region_start,Region_end,Nwin,extremeValue)
        
        """
        strandidx = 7
#         titlelist = [a[0].strip() for a in genomedbtools.operateDB("select", "select column_name  from information_schema.columns where table_schema='" + "ninglabvariantdata" + "' and table_name='" + trscptableName + "'")]
#         "select * from " + transcripttable + " where chrID='" + chrID + "' and trscpt_end_pos >= " + str(Region_start) + " and trscpt_start_pos <= " + str(Region_end)
        selectType1OverlapGenesql = "select * from " + transcripttable + " where chrID='" + chrID + "' and trscpt_start_pos >= " + str(Region_start) + " and trscpt_end_pos <= " + str(Region_end)
        selectType2OverlapGenesql = "select * from " + transcripttable + " where chrID='" + chrID + "' and trscpt_start_pos < " + str(Region_start) + " and trscpt_end_pos > " + str(Region_end)
        selectType3OverlapGenesql = "select * from " + transcripttable + " where chrID='" + chrID + "' and trscpt_start_pos < " + str(Region_start) + " and trscpt_end_pos > " + str(Region_start) + " and trscpt_end_pos < " + str(Region_end)
        selectType4OverlapGenesql = "select * from " + transcripttable + " where chrID='" + chrID + "' and trscpt_start_pos > " + str(Region_start) + " and trscpt_start_pos < " + str(Region_end) + " and trscpt_end_pos > " + str(Region_end)
        selectType5OverlapGenesql = "select * from " + transcripttable + " where chrID='" + chrID + "' and trscpt_start_pos < " + str(Region_start - upextend) + " and trscpt_end_pos > " + str(Region_start - upextend) + " and trscpt_end_pos < " + str(Region_end + downextend)
        selectType6OverlapGenesql = "select * from " + transcripttable + " where chrID='" + chrID + "' and trscpt_start_pos > " + str(Region_start - upextend) + " and trscpt_start_pos < " + str(Region_end + downextend) + " and trscpt_end_pos > " + str(Region_end + downextend)
        #1
        result = genomedbtools.operateDB("select", selectType1OverlapGenesql)
        for row in result:
            row += [1]
            trscptlist.append(row)
        #2
        result = genomedbtools.operateDB("select", selectType2OverlapGenesql)
        for row in result:
            row += [2]
            trscptlist.append(row)
        #3
        result = genomedbtools.operateDB("select", selectType3OverlapGenesql)
        for row in result:
            row += [3]
            trscptlist.append(row)
        #4
        result = genomedbtools.operateDB("select", selectType4OverlapGenesql)
        for row in result:
            row += [4]
            trscptlist.append(row)
        #5
        result = genomedbtools.operateDB("select", selectType5OverlapGenesql)
        for row in result:
            row += [5]
            trscptlist.append(row)
        #6
        result = genomedbtools.operateDB("select", selectType6OverlapGenesql)
        for row in result:
            row += [6]
            trscptlist.append(row)
        
        return trscptlist

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
        
     
        
        
        
        
        
