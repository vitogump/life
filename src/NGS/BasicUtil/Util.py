import copy
import random
import string
import re
import pickle
import os
import src.NGS.BasicUtil.DBManager as dbm
# from src.NGS.BasicUtil import *

'''
Created on 2013-6-30

@author: rui
'''
def complementary(seq):
    newseq = []
    for i in range(0, len(seq)):
        if seq[i].lower() == 'a':
            newseq.insert(i, 't')
        elif seq[i].lower() == 't':
            newseq.insert(i, 'a')
        elif seq[i].lower() == 'c':
            newseq.insert(i, 'g')
        else:
            newseq.insert(i, 'c')
    if isinstance(seq, str):
        newseq = "".join(newseq)
    return newseq

def random_str(randomlength=8):
    a = list(string.ascii_letters)
    random.shuffle(a)
    return ''.join(a[:randomlength])

def generateIndexByChrom(refFastaFileName, indexFileName):
    refFastaFile = open(refFastaFileName, 'r')
    refChromIndex = {}
    refline = refFastaFile.readline()
    while refline:
        if re.search(r'^[>]', refline) != None:
            collist = re.split(r'\s+', refline)
            currentChromNo = re.search(r'[^>]+', (re.split(r'\s+', refline))[0]).group(0)
            refChromIndex[currentChromNo] = int(refFastaFile.tell())  # from here is the sequence
        refline = refFastaFile.readline()
    pickle.dump(refChromIndex, open(indexFileName, 'wb'))
    refFastaFile.close()

def getGtfMap(gtfFileName):
    """protein_codingMap={chromNo:[[transcript_id,strand,start,end,(feature, elemStart, elemEnd, frame),(),(),,,,,],
                        [transcript_id,strand,start,end,(),(),(),,,],[],,,,,,,],
               chromNo:[],,,,,,,,,,,,,,}
        chrtranscrpitididxMap{chromNo:{transcript_id:ttanscript_id_idx,transcript_id:ttanscript_id_idx,,,,,},
                                chromNo:{},chromNo:{},,,,}
    """
    gtfFileHandler = open(gtfFileName, 'r')
    protein_codingMap = {}
    chrtranscrpitididxMap = {}
    gtfline = gtfFileHandler.readline()
    gtfColList = re.split(r'\s+', gtfline)
    chromNo = gtfColList[0].strip()
    protein_codingMap[chromNo] = []
    transcript_id = gtfColList[11]
    countInChrom = 0
    protein_codingMap[chromNo] = [[transcript_id, gtfColList[6], int(gtfColList[3]), int(gtfColList[4]), (gtfColList[2], int(gtfColList[3]), int(gtfColList[4]), gtfColList[7])]]
    chrtranscrpitididxMap[chromNo] = {transcript_id:0}
    for gtfline in gtfFileHandler:
        gtfColList = re.split(r'\s+', gtfline)
        transcript_id = gtfColList[11].strip()
        if "protein_coding" != gtfColList[1]:
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

    for chromNo in protein_codingMap.keys():
        protein_codingMap[chromNo].sort(key=lambda listRec:listRec[2])
    gtffilepath = re.search(r"^.*[/]", gtfFileName).group(0)
    testfile = open(gtffilepath + "protein_codingMap.sort.txt", 'w')    
    for chromNo in protein_codingMap.keys():
        protein_codingMap[chromNo].sort(key=lambda listRec:listRec[2])
        #先按照转录本起始坐标排序，下面是对转录本内元件排序，不过是什么排序方法忘记了，仔细读一下吧
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

def getRefSeqBypos(refFastahander, refindex, currentChromNO, startpos, endpos, seektuple=()):
    '''
    pos start at 1
    seektuple=(filepos,basesbeforefilepos)
    the refSeqMap has only one chromosome's sequence
    '''    
    refSeqMap = {}
    if startpos <= 0:
        startpos = 1
    print(currentChromNO, startpos, endpos)

    filehander = refFastahander
    if not seektuple or seektuple[1] > startpos:
        refSeqMap[currentChromNO] = [startpos - 1]
        filehander.seek(refindex[currentChromNO])  # seekmap is empty so go to the first bases of the currentChromNO
        preseq = filehander.read(startpos - 1)
        dn = preseq.count('\n')
        while dn != 0:
            preseq = filehander.read(dn)
            dn = preseq.count('\n')
            
        # now filehander is right stay at the startpos
        myseqline = filehander.read(endpos - startpos + 1)
        myseqn = myseqline.count('\n')
#        if len(myseqline)>200:
#            print(myseqn)
#            exit(-1)
#        print("myseqline=",myseqline,"myseqn", myseqn)
        while myseqn != 0:  # fill the same number of \n with bases
            myseqline = myseqline.replace('\n', '')
            myseqline += filehander.read(myseqn)
            myseqn = myseqline.count('\n')
            
#            print(currentChromNO,myseqline, myseqn)
            if myseqline.count('>') >= 1:
                exit(-1)
        refSeqMap[currentChromNO].extend(list(myseqline))
    else:
        filehander.seek(seektuple[0])  # seekmap is not empty
        refSeqMap[currentChromNO] = [startpos - 1]
        preseq = filehander.read(startpos - seektuple[1] - 1)
        dn = preseq.count('\n')
        while dn != 0:
            preseq = filehander.read(dn)
            dn = preseq.count('\n')
        # now filehander is right stay at the startpos
        myseqline = filehander.read(endpos - startpos + 1)
        myseqn = myseqline.count('\n')
        while myseqn != 0:  # fill the same number of \n with bases
            myseqline = myseqline.replace('\n', '')
            myseqline += filehander.read(myseqn)
            myseqn = myseqline.count('\n')
        refSeqMap[currentChromNO].extend(list(myseqline))
    plus = myseqline.count('>')
    if plus != 0:
        return -1
    
    return refSeqMap        


def getRefSeqMap(refFastafilehander, currentChromNO=None, preBaseTotal=0, linesOnce=500000):
    '''
    the refSeqMap has only one chromosome's sequence
    '''
    refSeqMap = {}
    if currentChromNO == None:
        refline = refFastafilehander.readline() 
        print("getRefSeqMap", refline)
        currentChromNO = re.search(r'[^>]+', (re.split(r'\s+', refline))[0]).group(0)
        refSeqMap[currentChromNO] = [preBaseTotal]  # preBaseTotal=0
        print("getRefSeqMap", currentChromNO)
    elif currentChromNO == "end of the reffile":
        return refSeqMap, currentChromNO, "end of the reffile"
    else:
        refSeqMap[currentChromNO] = [preBaseTotal]
    for refline in refFastafilehander:
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
    def __init__(self, gtfList, pos, RefSeqList):
        super().__init__()
        self.geneOverlapList = self.getNearestGeneOverlapList(gtfList, pos)
        self.tscptSeqAllCds = {}
        self.cds_frame = {}#{transcript_id:{cdsidx:(frame,startpos of this cds),cdsidx:(),,,,,}}
        for gene in self.geneOverlapList:
            
            genename = gene[0]
            self.tscptSeqAllCds[genename] = []
            self.cds_frame[genename] = {}#{cdsidx:(frame,startpos of this cds),cdsidx:(),,,,,}

            cdsidx = 3
            for feature, elemStart, elemEnd, frame in gene[4:]:
                cdsidx += 1
                if feature == 'CDS':
                    self.cds_frame[genename][cdsidx] = (int(frame), len(self.tscptSeqAllCds[genename]))
                    self.tscptSeqAllCds[genename] += RefSeqList[(elemStart - RefSeqList[0]):(elemEnd - RefSeqList[0] + 1)]#???如果不够呢
                elif  feature == "stop_codon":#feature == 'start_codon' or
                    self.cds_frame[genename][cdsidx] = (int(frame), len(self.tscptSeqAllCds[genename]))
                    self.tscptSeqAllCds[genename] += RefSeqList[(elemStart - RefSeqList[0]):(elemEnd - RefSeqList[0] + 1)]
            if genename=='"ENSAPLT00000005931";':
                print(genename,elemStart - RefSeqList[0],elemEnd - RefSeqList[0] + 1,elemStart - RefSeqList[0],elemEnd - RefSeqList[0] + 1,RefSeqList[0],len(RefSeqList))
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
            if gtfList[i][2]>=pos :
                geneOverlapList=[gtfList[i]]
                break
        else:
            if pos > gtfList[-1][3]:
                return []
            else:
                print(pos,gtfList)
                exit(-1)

#        high = len(gtfList) - 1
#        low = 0
#        mid = int((low + high) / 2)
#        while low < high:
#            mid = int((low + high) / 2)
#            if pos == (gtfList[mid][2]):
#                low = high#go to the else of the while block
#                high = mid
#            elif pos < (gtfList[mid][2]):
#                high = mid - 1
#            else:# pos > GtfMap[vcfChromNo][mid][2]:
#                low = mid + 1
#        else:
#            print("high:", high, "low:", low, "mid:", mid, file=open("testgetNearestGene.txt", 'a'))
#            if gtfList[high][3] >= pos and gtfList[high][2] <= pos:
#                geneOverlapList = [gtfList[high]];idx = high
#            elif gtfList[low][2] > pos:
#                geneOverlapList = [gtfList[low]];idx = low
#            elif low == high and low == 0:
#                geneOverlapList = [gtfList[0]];idx = 0
#            else:#out of end edge,so no gene after the pos,and returen a empty
#                return []
        print("getNearestGeneOverlapList", i, gtfList,pos)
        furthest = gtfList[i][3]
        i += 1
        while len(gtfList) > i and furthest >= gtfList[i][2]:
            if gtfList[i][0] != geneOverlapList[-1][0]:
                geneOverlapList.append(gtfList[i])
            furthest = max(furthest, gtfList[i][3])
            i += 1
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
        originallen = {}#just for test
        indelatEdgeofCDS=open("indelatEdgeofCDS.txt",'w')
        for gene in self.geneOverlapList:
            genename = gene[0]
            print(gene, self.cds_frame[genename], sep="\n", file=open("testgeneOverlapList.txt", 'a'))
            
            tscptSeqAllCds_mut[genename] = copy.copy(self.tscptSeqAllCds[genename])
            originallen[genename] = len(tscptSeqAllCds_mut[genename])#just for test
        
        while idx_vcf!=-1 and idx_vcf != len(VcfList) and VcfList[idx_vcf][0] <= self.geneOverlapList[-1][3]:
            vcfpos = VcfList[idx_vcf][0];refalle = VcfList[idx_vcf][1];altalle = VcfList[idx_vcf][2]
            cns_append += ("".join(RefSeqList[idx_RefSeq:idx_RefSeq + (vcfpos - curpos)]))
            idx_RefSeq += (vcfpos - curpos)
            curpos = RefSeqList[0] + idx_RefSeq
            
            if re.search(r'[^a-zA-Z]', altalle) != None:#contain ',' ie. multiple alle
                cns_append += ("".join(RefSeqList[idx_RefSeq:idx_RefSeq + len(refalle)]))
                idx_RefSeq += len(refalle)
                curpos = RefSeqList[0] + idx_RefSeq
                continue
            cns_append += ("".join(RefSeqList[idx_RefSeq:idx_RefSeq + len(altalle)]))
            idx_RefSeq += len(refalle)#here should still be refalle
            curpos = RefSeqList[0] + idx_RefSeq
            n_refbases = len(refalle);n_altbases = len(altalle)#situation TAA     TA;     TTA     TTAAACTTCTATACTA;      C       T;    T       TATA;    ACG     A
# first for every variant making cns_append,and then substitute the seq in the cds seq,and finialy translate to protein 
            for gene in self.geneOverlapList:
                genename = gene[0]
                if gene[2] <= vcfpos and gene[3] >= vcfpos:
                    t4_indx = 3
                    for feature, elemStart, elemEnd, frame in gene[4:]:
                        t4_indx += 1
                        if feature == 'CDS' and vcfpos <= elemEnd and vcfpos >= elemStart:
                            if vcfpos + n_refbases - 1 > elemEnd or vcfpos + n_altbases - 1 > elemEnd:
                                print(VcfList[idx_vcf],genename,file=indelatEdgeofCDS)
                                break
                            if n_refbases > n_altbases:#situation TAA     TA;ACG     A
                                print(tscptSeqAllCds_mut[genename][(vcfpos - elemStart + self.cds_frame[genename][t4_indx][1]):(vcfpos - elemStart + self.cds_frame[genename][t4_indx][1] + n_refbases)])
                                tscptSeqAllCds_mut[genename][(vcfpos - elemStart + self.cds_frame[genename][t4_indx][1]):(vcfpos - elemStart + self.cds_frame[genename][t4_indx][1] + n_altbases)] = list(altalle)
                                tscptSeqAllCds_mut[genename][(vcfpos - elemStart + self.cds_frame[genename][t4_indx][1] + n_altbases):(vcfpos - elemStart + self.cds_frame[genename][t4_indx][1] + n_refbases)] = [' '] * (n_refbases - n_altbases)
                                print(tscptSeqAllCds_mut[genename][(vcfpos - elemStart + self.cds_frame[genename][t4_indx][1]):(vcfpos - elemStart + self.cds_frame[genename][t4_indx][1] + n_refbases)])
                            elif n_refbases < n_altbases:#situation TTA     TTAAACTTCTATACTA;T       TATA;
                                if n_refbases == 1:
                                    tscptSeqAllCds_mut[genename][vcfpos - elemStart + self.cds_frame[genename][t4_indx][1]] = altalle
                                else:
                                    tscptSeqAllCds_mut[genename][(vcfpos - elemStart + self.cds_frame[genename][t4_indx][1]):(vcfpos - elemStart + self.cds_frame[genename][t4_indx][1] + n_refbases - 1)] = list(altalle[0:(n_refbases - 1)])
                                    tscptSeqAllCds_mut[genename][vcfpos - elemStart + self.cds_frame[genename][t4_indx][1] + n_refbases - 1] = altalle[(n_refbases - 1):]
                            else:#n_refbases==n_altbases==1
                                try:
                                    tscptSeqAllCds_mut[genename][vcfpos - elemStart + self.cds_frame[genename][t4_indx][1]] = altalle
                                except IndexError:
                                    print(genename, vcfpos, t4_indx, altalle, elemStart, feature, len(tscptSeqAllCds_mut[genename]))
                                    exit()
       
            idx_vcf += 1
#该翻译蛋白了吧 还有 看看长度一样不  将最后一个vcf记录之后的序列加入一致序列字符串
        cns_append += "".join(RefSeqList[idx_RefSeq:idx_RefSeq + (self.geneOverlapList[-1][3] - curpos) + 1])
        for gene in self.geneOverlapList:
            
            genename = gene[0]
#            ref_amino_seq[genename] = []
            mutat_amino_seq[genename] = []
#            mutationTypeList=[]
            if originallen[genename] != len(tscptSeqAllCds_mut[genename]):#just for test
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
            else:# strand == '-'
                tscptSeqAllCds_Revr_Cmplm = complementary(self.tscptSeqAllCds[genename])
                tscptSeqAllCds_Revr_Cmplm.reverse()                
                tscptSeqAllCds_mut_str = list(tscptSeqAllCds_mut_str)
                tscptSeqAllCds_mut_str_Revr_Cmplm = complementary(tscptSeqAllCds_mut_str)
                tscptSeqAllCds_mut_str_Revr_Cmplm.reverse()
                tscptSeqAllCds_mut_str_Revr_Cmplm = "".join(tscptSeqAllCds_mut_str_Revr_Cmplm)
                tscptSeqAllCds_mut[genename] = list(tscptSeqAllCds_mut_str_Revr_Cmplm)
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
        linelist = re.split(r"\s+", line)
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
        #set the filehandler locate at the nearest location to the target location
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
        chrom = re.search(r"^([\w\W]*)[:]([\d]*)", linelist[0]).group(1)
        pos = int(re.search(r"^([\w\W]*)[:]([\d]*)", linelist[0]).group(2))
        if chrom == targetchr and pos == targetloc:
            return linelist
        if chrom != targetchr and pos > targetloc:
            self.depthfilefp.seek(self.covfileidx[targetchr])
        elif chrom == targetchr and pos < targetloc - 100:
            pass#use the lastposoffilehandler to set the filehanlder quickly
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
        L = [(pos, REF, ALT, INFO),(),(),...........]
        """
        self.winValueL = []  # notice here
        nextIdx = -1
        currentIdx = 0
        winStart = 0
        FoundNextIdx = False
        firstComeInWin = True
        while currentIdx != len(L):
            if L[currentIdx][0] > winStart and L[currentIdx][0] <= (winStart + windowWidth):
                if re.search(r"INDEL", L[currentIdx][3]) == None and True:
                    if firstComeInWin:
                        startPos = L[currentIdx][0]
                        firstComeInWin = False
                    lastPos = L[currentIdx][0]
                    Caculator.process(L[currentIdx])
                if FoundNextIdx == False and L[currentIdx][0] > (winStart + slideSize):  # always go to |currentIdx+=1|
                    nextIdx = currentIdx
                    FoundNextIdx = True
            else:
                value = Caculator.getResult()
                try:
                    self.winValueL.append((startPos, lastPos, value))
                except UnboundLocalError:
                    self.winValueL.append((0, 0, value))
                winStart += slideSize
                firstComeInWin = True
                
                FoundNextIdx = False
                if nextIdx == -1:
                    if slideSize >= windowWidth:
                        while currentIdx != len(L):
                            if L[currentIdx][0] > winStart and L[currentIdx][0] <= (winStart + windowWidth):
                                break
                            currentIdx += 1
                    continue  # go to |if L[currentIdx][0] > winStart and L[currentIdx][0] < (winStart + windowWidth):| in upside block
                else:
                    currentIdx = nextIdx
                    nextIdx = -1
                    continue
                
            currentIdx += 1
        else:
            value = Caculator.getResult()
            try:
                self.winValueL.append((startPos, lastPos, value))
            except UnboundLocalError:
                self.winValueL.append((0, 0, value))
        
        n = int((L_End_Pos - (len(self.winValueL) * slideSize + windowWidth)) / slideSize) + 1
        for i in range(n):
            self.winValueL.append((0, 0, 'NA'))
        
class WinInGenome():           
    def __init__(self, dbname, winFileName6Field, tableName=None):
        super().__init__()
#        self.wintable="PoMuJVOr"
#        self.windbtools = dbm.DBTools("localhost", "root", "1234567", dbname)
        self.windbtools, self.wintable = self.loadWinDataIntoDB(dbname, winFileName6Field, tableName)
        self.winContainTrscptMap = {}
    def loadWinDataIntoDB(self, dbname, winFileName6Field, tableName=None):
        if tableName == None:
            tableName = random_str()
        tempdbtools = dbm.DBTools("localhost", "root", "1234567", dbname)
        TABLES = {}
        TABLES[tableName] = (
            "CREATE TABLE " + tableName + " ("
            " `chrID` varchar(128) NOT NULL ,"
            " `winNo` varchar(128) NOT NULL,"
            " `bp_start` varchar(128) NOT NULL,"
            " `bp_end` varchar(128) NOT NULL,"
            " `value` text NOT NULL,"
            " `zvalue` text NOT NULL,"
            " PRIMARY KEY (`chrID`,`winNo`)"
            ")"
            )
        
        tempdbtools.create_table(TABLES)
        loaddatasql = "load data local infile '" + winFileName6Field + "' into table " + tableName + " fields terminated by '\\t'"
        shellstatment = "mysql -uroot -p1234567 -D" + dbname.strip() + ' -e "' + loaddatasql + '"'
        print(shellstatment)
        a = os.system(shellstatment)
        if a != 0:
            print("loadWinDataIntaDB func os.system return not 0")
            exit(-1)
        print(a)
#        tempdbtools.load_file(tableName,"chrID","winNo","bp_start","bp_end","value","zvalue",fileName=winFileName6Field)
        return tempdbtools, tableName       
    def collectTrscptInWin(self, dbtools, trscptableName, vcftable, winRegion):
        transcripttable = trscptableName
        chrID = winRegion[0]
        winNo = int(winRegion[1])
        winWidth = int(winRegion[2])
        slideSize = int(winRegion[3])
        """
        winRegion=(chrID,winNo,winWidth,slideSize,zvalue)
        """

        selectsql = "select * from " + transcripttable + " where chrID='" + chrID + "' and end_pos >= " + str(winNo * slideSize) + " and start_pos <= " + str(winNo * slideSize + winWidth)
        result = dbtools.operateDB("select", selectsql)
        self.winContainTrscptMap[winRegion] = []
        for row in result:
            self.winContainTrscptMap[winRegion].append(row)

class Node(object):
    def __init__(self, val, p=0):
        self.data = val
        self.next = p

class LinkList(object):
    def __init__(self):
        self.head = 0

    def __getitem__(self, key):

        if self.is_empty():
            print('linklist is empty.')
            return

        elif key < 0  or key > self.getlength():
            print('the given key is error')
            return

        else:
            return self.getitem(key)



    def __setitem__(self, key, value):

        if self.is_empty():
            print('linklist is empty.')
            return

        elif key < 0  or key > self.getlength():
            print('the given key is error')
            return

        else:
            self.delete(key)
            return self.insert(key)

    def initlist(self, data):

        self.head = Node(data[0])

        p = self.head

        for i in data[1:]:
            node = Node(i)
            p.next = node
            p = p.next

    def getlength(self):

        p = self.head
        length = 0
        while p != 0:
            length += 1
            p = p.next

        return length

    def is_empty(self):

        if self.getlength() == 0:
            return True
        else:
            return False

    def clear(self):

        self.head = 0


    def append(self, item):

        q = Node(item)
        if self.head == 0:
            self.head = q
        else:
            p = self.head
            while p.next != 0:
                p = p.next
            p.next = q


    def getitem(self, index):

        if self.is_empty():
            print('Linklist is empty.')
            return
        j = 0
        p = self.head

        while p.next != 0 and j < index:
            p = p.next
            j += 1

        if j == index:
            return p.data

        else:

            print('target is not exist!')

    def insert(self, index, item):

        if self.is_empty() or index < 0 or index > self.getlength():
            print('Linklist is empty.')
            return

        if index == 0:
            q = Node(item, self.head)

            self.head = q

        p = self.head
        post = self.head
        j = 0
        while p.next != 0 and j < index:
            post = p
            p = p.next
            j += 1

        if index == j:
            q = Node(item, p)
            post.next = q
            q.next = p


    def delete(self, index):

        if self.is_empty() or index < 0 or index > self.getlength():
            print('Linklist is empty.')
            return

        if index == 0:
            self.head = self.head.next
            return
#            q = Node(item,self.head)
#
#            self.head = q

        p = self.head
        post = self.head
        j = 0
        while p.next != 0 and j < index:
            post = p
            p = p.next
            j += 1

        if index == j:
            post.next = p.next

    def index(self, value):

        if self.is_empty():
            print('Linklist is empty.')
            return

        p = self.head
        i = 0
        while p.next != 0 and not p.data == value:
            p = p.next
            i += 1

        if p.data == value:
            return i
        else:
            return -1        
        
        
        
        
        
