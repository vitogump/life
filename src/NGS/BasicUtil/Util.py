import re, pickle
#from src.NGS.BasicUtil import *

'''
Created on 2013-6-30

@author: rui
'''
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
    def getConsenusSeqMap(fastQFileName,dbtools ,tablename="chromosome",primaryID = "chrID",bp_start=None, bp_end=None):
        '''
        the refSeqMap has only one chromosome's sequence
        '''
#        print(dbtools,fastQFileName,"inside FastQ_Util")
        fqfile=open(fastQFileName,'r')
        sql="select * from "+tablename
        seqMapByChrom = {}
        try:
            ChromIndexMap = pickle.load(open(fastQFileName + ".myindex", 'rb'))
        except IOError:
            FastQ_Util.generateIndexByChrom(fastQFileName, fastQFileName + ".myindex")
            ChromIndexMap = pickle.load(open(fastQFileName + ".myindex", 'rb'))
        
        temresult = dbtools.operateDB("select","select count(*) from "+tablename)
        
        print(temresult)
        totalChroms = temresult[0][0]
        currentchrID=dbtools.operateDB("select",sql+" limit 0,1")[0][0]
        seqMapByChrom[currentchrID]=""
        for i in range(0,totalChroms-1,20):
            currentsql=sql+" order by "+primaryID+" limit "+str(i)+",20"
            result=dbtools.operateDB("select",currentsql)
            for row in result:
                currentchrID=row[0]
                if currentchrID in ChromIndexMap:
                    seqMapByChrom[currentchrID]=""
                    fqfile.seek(ChromIndexMap[currentchrID])
                    line=fqfile.readline()
                    while line.strip() !="+":
                        seqMapByChrom[currentchrID]+=line.strip()
#                        print(line.strip())
                        line=fqfile.readline()
        return seqMapByChrom
#        
#        if currentChromNO == None:
#            refline = fastQFileName.readline() 
#            print(refline)
#            currentChromNO = re.search(r'[^>]+', (re.split(r'\s+', refline))[0]).group(0)
#            refSeqMap[currentChromNO] = [preBaseTotal]#preBaseTotal=0
#            print(currentChromNO)
#        else:
#            refSeqMap[currentChromNO] = [preBaseTotal]
#        for refline in fastQFileName:
#            if re.search(r'^[>]', refline) != None:
#                collist = re.split(r'\s+', refline)
#                print(re.search(r'[^>]+', collist[0]).group(0))
#    #            refSeqMap[currentChromNO] = [0]
#                return refSeqMap,currentChromNO#clean the refSeqMap and report the current chromNO
#            else:
#                refSeqMap[currentChromNO].extend(list(refline.strip().lower()))
#            linesOnce -= 1    
#            if linesOnce == 0:
#                break                
#        return refSeqMap, currentChromNO
    
    
class Window():
    def __init__(self):
        super().__init__()
        self.winValueL = []  # [(startPos,lastPos,value),(),,,,,,]
    def slidWindowOverlap(self, L, windowWidth, slideSize, Caculator):
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
            
                
