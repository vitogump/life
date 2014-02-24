import re, pickle, os
import random, string
import src.NGS.BasicUtil.DBManager as dbm
# from src.NGS.BasicUtil import *

'''
Created on 2013-6-30

@author: rui
'''
def complementary(seqlist):
    newseqlist = []
    for i in range(0, len(seqlist)):
        if seqlist[i].lower() == 'a':
            newseqlist.insert(i, 't')
        elif seqlist[i].lower() == 't':
            newseqlist.insert(i, 'a')
        elif seqlist[i].lower() == 'c':
            newseqlist.insert(i, 'g')
        elif seqlist[i].lower() == 'g':
            newseqlist.insert(i, 'c')
        else:
            newseqlist.insert(i, seqlist[i])
    return newseqlist
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

def getGtfMap(gtfFileHandler):
    """gtfMap={chromNo:[[transcript_id,strand,start,end,(feature, elemStart, elemEnd, frame),(),(),,,,,],
                        [transcript_id,strand,start,end,(),(),(),,,],[],,,,,,,],
               chromNo:[],,,,,,,,,,,,,,}
        chrtranscrpitididxMap{chromNo:{transcript_id:ttanscript_id_idx,transcript_id:ttanscript_id_idx,,,,,},
                                chromNo:{},chromNo:{},,,,}
    """
    gtfMap = {}
    chrtranscrpitididxMap={}
    gtfline = gtfFileHandler.readline()
    gtfColList = re.split(r'\s+', gtfline)
    chromNo = gtfColList[0].strip()
    gtfMap[chromNo] = []
    transcript_id = gtfColList[11]
    countInChrom = 0
    gtfMap[chromNo] = [[transcript_id, gtfColList[6], int(gtfColList[3]), int(gtfColList[4]), (gtfColList[2], int(gtfColList[3]), int(gtfColList[4]), gtfColList[7])]]
    chrtranscrpitididxMap[chromNo]={transcript_id:0}
    for gtfline in gtfFileHandler:
        gtfColList = re.split(r'\s+', gtfline)
        transcript_id = gtfColList[11].strip()
        chromNo = gtfColList[0].strip()
        if chromNo in gtfMap:
            if transcript_id in chrtranscrpitididxMap[chromNo].keys():
                tanscript_id_idx=chrtranscrpitididxMap[chromNo][transcript_id]
                gtfMap[chromNo][tanscript_id_idx].append((gtfColList[2], int(gtfColList[3]), int(gtfColList[4]), gtfColList[7]))
                gtfMap[chromNo][tanscript_id_idx][2]=min(gtfMap[chromNo][tanscript_id_idx][2],int(gtfColList[3]))
                gtfMap[chromNo][tanscript_id_idx][3]=max(gtfMap[chromNo][tanscript_id_idx][3],int(gtfColList[4]))
            else:
                gtfMap[chromNo].append([transcript_id, gtfColList[6], int(gtfColList[3]), int(gtfColList[4]), (gtfColList[2], int(gtfColList[3]), int(gtfColList[4]), gtfColList[7])])
                chrtranscrpitididxMap[chromNo][transcript_id]=len(gtfMap[chromNo])-1
        else:
             gtfMap[chromNo] = [[transcript_id, gtfColList[6], int(gtfColList[3]), int(gtfColList[4]), (gtfColList[2], int(gtfColList[3]), int(gtfColList[4]), gtfColList[7])]]
             chrtranscrpitididxMap[chromNo]={transcript_id:0}
    else:
        pass                 

    for chromNo in gtfMap.keys():
        gtfMap[chromNo].sort(key=lambda listRec:listRec[2])   
    testfile = open("gtfMap.sort.txt", 'w')    
    for chromNo in gtfMap.keys():
        gtfMap[chromNo].sort(key=lambda listRec:listRec[2])
        #先按照转录本起始坐标排序，下面是对转录本内元件排序，不过是什么排序方法忘记了，仔细读一下吧
        for j in range(len(gtfMap[chromNo])):
            for t4_indx in range(5, len(gtfMap[chromNo][j])):
                t4_key = gtfMap[chromNo][j][t4_indx]
                t4_indxp = t4_indx - 1
                while t4_indxp >= 4 and gtfMap[chromNo][j][t4_indxp][1] > t4_key[1]:
                    gtfMap[chromNo][j][t4_indxp + 1] = gtfMap[chromNo][j][t4_indxp]
                    t4_indxp = t4_indxp - 1
                else:
                    gtfMap[chromNo][j][t4_indxp + 1] = t4_key
#            print(gtfMap[chromNo][j],gtfMap[chromNo][j][t4_indxp][1] > t4_key[1],t4_indxp)
        print(chromNo,"num of transcrpit:",len(gtfMap[chromNo]),file=testfile)
        for i in range(len(gtfMap[chromNo])):
#            print(gtfMap[chromNo][i][0],gtfMap[chromNo][i][1],gtfMap[chromNo][i][2],gtfMap[chromNo][i][3])
            for k in range(len(gtfMap[chromNo][i])):
                print(gtfMap[chromNo][i][k], file=testfile)
    testfile.close()
    return gtfMap

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
#    try:
#        refindex = pickle.load(open(refFastaFileName + ".myindex", 'rb'))
#    except IOError:
#        generateIndexByChrom(refFastaFileName, refFastaFileName + ".myindex")
#        refindex = pickle.load(open(refFastaFileName + ".myindex", 'rb'))
#    filehander = open(refFastaFileName, 'r')
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
#    filehander.close()
    plus = myseqline.count('>')
    if plus != 0:
        return -1
    
    return refSeqMap        
#    filehander.close()

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
    def __init__(self,gtfList,pos):
        super.__init__()
        self.geneOverlapList=self.getNearestGeneOverlapList(gtfList, pos)
        self.tscptSeqAllCds = []
        self.cds_frame = {}#{cdsidx:(frame,startpos of this cds),cdsidx:(),,,,,}
        for gene in self.geneOverlapList:
            if gene[1]=="+":
                cdsidx=3
                for feature,elemStart,elemEnd,frame in gene[4:]:
                    cdsidx+=1
                    if feature == 'CDS':
                        cds_frame[cdsidx]=(int(frame), len(self.tscptSeqAllCds))
                        self.tscptSeqAllCds += RefSeqMap[]???如果不够呢
    def getNearestGeneOverlapList(self,gtfList,pos):
        """
        input:for a chrom,contain all transcript of this chrom
        gtfList=[[transcript_id,strand,start,end,(feature, elemStart, elemEnd, frame),(),(),,,,,],
                            [transcript_id,strand,start,end,(),(),(),,,],[],,,,,,,]
        return: the first gene that after the pos and the genes contain in or overlap with or contact with this gene indirect
        geneOverlapList=[[transcript_id,strand,start,end,(feature, elemStart, elemEnd, frame),(),(),,,,,],
                            [transcript_id,strand,start,end,(),(),(),,,],[],,,,,,,]
        order by "start"
        """
        high = len(gtfList) - 1
        low = 0
        while low < high:
            mid = int((low + high) / 2)
            if pos == (gtfList[mid][2]):
                low = high#go to the else of the while block
                high = mid
            elif pos < (gtfList[mid][2]):
                high = mid - 1
            else:# snpPos > GtfMap[vcfChromNo][mid][2]:
                low = mid + 1
        else:
            if gtfList[high][3]>=pos and gtfList[high][2]<=pos:
                geneOverlapList=[gtfList[high]];idx=high
            elif gtfList[low][2]>pos:
                geneOverlapList=[gtfList[low]];idx=low
            elif low ==high and low == 0:
                geneOverlapList=[gtfList[0]];idx=0
            else:#out of end edge,so no gene after the pos,and returen a empty
                return []
            furthest=gtfList[idx][3]
            idx+=1
            while furthest >= gtfList[idx][2]:
                geneOverlapList.append(gtfList[idx])
                furthest=max(furthest,gtfList[idx][3])
                idx+=1
        return geneOverlapList

class GATK_depthfile():
    def __init__(self, depthfileName, indexFileName):
        super.__init__()
        self.covfileidx = {}
        self.title=[]
        try:
            self.covfileidx = pickle.load(open(indexFileName, 'rb'))
        except IOError:
            self.indexGATK_depthfile(depthfileName, indexFileName)
            self.covfileidx = pickle.load(open(indexFileName, 'rb'))
        self.depthfilehandler = open(depthfileName, 'r')
    
    def indexGATK_depthfile(self, depthfileName, indexFileName):
        """
        {chrom:position_in_file_of_first_genomepos_of_this_chrom,chrom:position,,,,,,}
        """
        depthfile = open(depthfileName, 'r')
        covfileidx = {}
      
        currentChrom = None
        lastPosition = 0
        line = depthfile.readline()
        linelist = re.split(r"\s+", line)
        self.title=linelist
        print("title",line,linelist)
        while line:      
            linelist = re.split(r"\s+", line)
            if currentChrom != re.search(r"^([\w\W]*)[:]([\d]*)", linelist[0]).group(1):
                currentChrom = re.search(r"^([\w\W]*)[:]([\d]*)", linelist[0]).group(1)
                covfileidx[currentChrom] = lastPosition
            lastPosition = depthfile.tell()
    
            line = depthfile.readline()
        pickle.dump(covfileidx, open(indexFileName, 'wb'))
        depthfile.close()
    def set_depthfilehandler(self, locchrom, locingenome, lastposoffilehandler=0):
        """
        set the self.depthfilehandler to the line in the file where chrom==locchrom locingenome==locingenome
        """
        self.depthfilehandler.seek(lastposoffilehandler)
        linelist = re.split(r"\s+", self.depthfilehandler.readline())
        currentChrom = re.search(r"^([\w\W]*)[:]([\d]*)", linelist[0]).group(1)
        pos = re.search(r"^([\w\W]*)[:]([\d]*)", linelist[0]).group(2)
        if currentChrom == locchrom and pos <= locingenome:
            pass
        else:
            self.depthfilehandler.seek(self.covfileidx[locchrom])
            line = self.depthfilehandler.readline()
            currentChrom = re.search(r"^([\w\W]*)[:]([\d]*)", linelist[0]).group(1)
            pos = re.search(r"^([\w\W]*)[:]([\d]*)", linelist[0]).group(2)            
        #set the filehandler locate at the nearest location to the target location
        while currentChrom == locchrom:
            if pos == locingenome:
                return "found"
            line = self.depthfilehandler.readline()
            linelist = re.split(r"\s+", line)
            currentChrom = re.search(r"^([\w\W]*)[:]([\d]*)", linelist[0]).group(1)
            pos = re.search(r"^([\w\W]*)[:]([\d]*)", linelist[0]).group(2)
        else:
            return "didn't find"         
    def getnextposline(self):
        line = self.depthfilehandler.readline()
        linelist = re.split(r"\s+", line)
        chrom = re.search(r"^([\w\W]*)[:]([\d]*)", linelist[0]).group(1)
        pos = re.search(r"^([\w\W]*)[:]([\d]*)", linelist[0]).group(2)
        return chrom,pos,line,linelist
    def closedepthfile(self):
        self.title.clear()
        self.covfileidx.clear()
        self.depthfilehandler.close()
                
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
    def __init__(self,val,p=0):
        self.data = val
        self.next = p

class LinkList(object):
    def __init__(self):
        self.head = 0

    def __getitem__(self, key):

        if self.is_empty():
            print('linklist is empty.')
            return

        elif key <0  or key > self.getlength():
            print('the given key is error')
            return

        else:
            return self.getitem(key)



    def __setitem__(self, key, value):

        if self.is_empty():
            print('linklist is empty.')
            return

        elif key <0  or key > self.getlength():
            print('the given key is error')
            return

        else:
            self.delete(key)
            return self.insert(key)

    def initlist(self,data):

        self.head = Node(data[0])

        p = self.head

        for i in data[1:]:
            node = Node(i)
            p.next = node
            p = p.next

    def getlength(self):

        p =  self.head
        length = 0
        while p!=0:
            length+=1
            p = p.next

        return length

    def is_empty(self):

        if self.getlength() ==0:
            return True
        else:
            return False

    def clear(self):

        self.head = 0


    def append(self,item):

        q = Node(item)
        if self.head ==0:
            self.head = q
        else:
            p = self.head
            while p.next!=0:
                p = p.next
            p.next = q


    def getitem(self,index):

        if self.is_empty():
            print('Linklist is empty.')
            return
        j = 0
        p = self.head

        while p.next!=0 and j <index:
            p = p.next
            j+=1

        if j ==index:
            return p.data

        else:

            print('target is not exist!')

    def insert(self,index,item):

        if self.is_empty() or index<0 or index >self.getlength():
            print('Linklist is empty.')
            return

        if index ==0:
            q = Node(item,self.head)

            self.head = q

        p = self.head
        post  = self.head
        j = 0
        while p.next!=0 and j<index:
            post = p
            p = p.next
            j+=1

        if index ==j:
            q = Node(item,p)
            post.next = q
            q.next = p


    def delete(self,index):

        if self.is_empty() or index<0 or index >self.getlength():
            print('Linklist is empty.')
            return

        if index ==0:
            self.head=self.head.next
            return
#            q = Node(item,self.head)
#
#            self.head = q

        p = self.head
        post  = self.head
        j = 0
        while p.next!=0 and j<index:
            post = p
            p = p.next
            j+=1

        if index ==j:
            post.next = p.next

    def index(self,value):

        if self.is_empty():
            print('Linklist is empty.')
            return

        p = self.head
        i = 0
        while p.next!=0 and not p.data ==value:
            p = p.next
            i+=1

        if p.data == value:
            return i
        else:
            return -1        
        
        
        
        
        
