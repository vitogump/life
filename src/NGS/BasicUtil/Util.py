import re, pickle, os
import random, string
import src.NGS.BasicUtil.DBManager as dbm
#from src.NGS.BasicUtil import *

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
            newseqlist.insert(i,seqlist[i])
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
            refChromIndex[currentChromNo] = int(refFastaFile.tell())# from here is the sequence
        refline = refFastaFile.readline()
    pickle.dump(refChromIndex, open(indexFileName, 'wb'))
    refFastaFile.close()
def getRefSeqBypos(refFastahander,refindex, currentChromNO, startpos, endpos, seektuple=()):
    '''
    pos start at 1
    seektuple=(filepos,basesbeforefilepos)
    the refSeqMap has only one chromosome's sequence
    '''    
    refSeqMap = {}
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
        filehander.seek(refindex[currentChromNO])#seekmap is empty so go to the first bases of the currentChromNO
        preseq = filehander.read(startpos - 1)
        dn = preseq.count('\n')
        while dn != 0:
            preseq = filehander.read(dn)
            dn = preseq.count('\n')
            
        #now filehander is right stay at the startpos
        myseqline = filehander.read(endpos - startpos + 1)
        myseqn = myseqline.count('\n')
        if len(myseqline)>200:
            print(myseqn)
            exit(-1)
        print("myseqline=",myseqline,"myseqn", myseqn)
        while myseqn != 0:# fill the same number of \n with bases
            myseqline = myseqline.replace('\n', '')
            myseqline += filehander.read(myseqn)
            myseqn = myseqline.count('\n')
            
            print(currentChromNO,myseqline, myseqn)
            if myseqline.count('>') >= 1:
                exit(-1)
        refSeqMap[currentChromNO].extend(list(myseqline))
    else:
        filehander.seek(seektuple[0])#seekmap is not empty
        refSeqMap[currentChromNO] = [startpos - 1]
        preseq = filehander.read(startpos - seektuple[1] - 1)
        dn = preseq.count('\n')
        while dn != 0:
            preseq = filehander.read(dn)
            dn = preseq.count('\n')
        #now filehander is right stay at the startpos
        myseqline = filehander.read(endpos - startpos + 1)
        myseqn = myseqline.count('\n')
        while myseqn != 0:# fill the same number of \n with bases
            myseqline = myseqline.replace('\n', '')
            myseqline += filehander.read(myseqn)
            myseqn = myseqline.count('\n')
        refSeqMap[currentChromNO].extend(list(myseqline))
#    filehander.close()
    plus = myseqline.count('>')
    if plus != 0:
        return -1
    
    return refSeqMap

        
    filehander.close()
def getRefSeqMap(refFastafilehander, currentChromNO=None, preBaseTotal=0, linesOnce=500000):
    '''
    the refSeqMap has only one chromosome's sequence
    '''
    refSeqMap = {}
    if currentChromNO == None:
        refline = refFastafilehander.readline() 
        print("getRefSeqMap", refline)
        currentChromNO = re.search(r'[^>]+', (re.split(r'\s+', refline))[0]).group(0)
        refSeqMap[currentChromNO] = [preBaseTotal]#preBaseTotal=0
        print("getRefSeqMap", currentChromNO)
    else:
        refSeqMap[currentChromNO] = [preBaseTotal]
    for refline in refFastafilehander:
        if re.search(r'^[>]', refline) != None:
            collist = re.split(r'\s+', refline)
            print("getRefSeqMap", re.search(r'[^>]+', collist[0]).group(0))
#            refSeqMap[currentChromNO] = [0]
#            currentChromNO=re.search(r'[^>]+', collist[0]).group(0)
            return refSeqMap, currentChromNO#clean the refSeqMap and report the current chromNO
        else:
            refSeqMap[currentChromNO].extend(list(refline.strip().lower()))
        linesOnce -= 1    
        if linesOnce == 0:
            break                
    return refSeqMap, currentChromNO


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
                if len(collist[0]) > 20:# may be fqline is located in the quality value block
                    fqline = fasqfile.readline()
                    continue
#                print(collist[0],fqline)
                currentChromNo = re.search(r'^[@]([^@]+)$', collist[0]).group(1).strip()
#                print(currentChromNo)
                refChromIndex[currentChromNo] = int(fasqfile.tell())# from here is the sequence
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

        
        
        
        
        
        
