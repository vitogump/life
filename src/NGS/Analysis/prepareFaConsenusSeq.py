import sys,NGS.BasicUtil.Util,pickle

import src.NGS.BasicUtil.DBManager as dbm

from NGS.BasicUtil import *
'''
Created on 2013-8-23

@author: liurui
'''
if len(sys.argv) < 2:
    print("python prepareFaConsenusSeq.py [cns1.fq] [cns2.fq] [cns3.fq] ...")
    exit(-1)
tablename='chromosome'
primaryID="chrID"
sql="select * from "+tablename
allchr=""
allseq=""
if __name__ == '__main__':
#    ChromIndexMap = pickle.load(open(fastQFileName + ".myindex", 'rb'))
    dbtools = dbm.DBTools("localhost","root","1234567","life_pilot")
    print(dbtools,"ssssssssssssss")
    for fastQFileName in sys.argv[1:]:
        print(fastQFileName)
        outfile=open(fastQFileName+".fa",'w')
        seqMapByChrom = Util.FastQ_Util.getConsenusSeqMap(fastQFileName, dbtools)

        totalChroms = dbtools.operateDB("select","select count(*) from "+tablename)[0][0]
    #    currentchrID=dbtools.operateDB("select",sql+" limit 0,1")[0][0]
    #    seqMapByChrom[currentchrID]=""
        for i in range(0,totalChroms,20):
            currentsql=sql+" order by "+primaryID+" limit "+str(i)+",20"
            result=dbtools.operateDB("select",currentsql)
            for row in result:
                currentchrID=row[0]
                if currentchrID in seqMapByChrom:
                    allchr+=currentchrID+"\n"
                    allseq+=seqMapByChrom[currentchrID]
                    seqMapByChrom.pop(currentchrID)
        for i in range(int(len(allseq)/100)):
            print(allseq[i*100:(i+1)*100],file=outfile)
        print(len(allseq)," ",i)
        if len(allseq)/100>int(len(allseq)/100):
            print(allseq[(i+1)*100:],file=outfile)             
        
        print(allchr,file=open(fastQFileName+"allchr.txt",'w'))
        outfile.close()
    dbtools.disconnect()
                