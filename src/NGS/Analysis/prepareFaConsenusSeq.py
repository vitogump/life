import sys,NGS.BasicUtil.Util,pickle
from NGS.BasicUtil import *
import src.NGS.BasicUtil.DBManager as dbm

'''
Created on 2013-8-23

@author: liurui
'''
tablename='chromosome'
primaryID="chrID"
fastQFileName=sys.argv[1]
outfile=open("testoutfile.cns.fa",'w')
sql="select * from "+tablename
allchr=""
if __name__ == '__main__':
#    ChromIndexMap = pickle.load(open(fastQFileName + ".myindex", 'rb'))
    dbtools = dbm.DBTools()
    seqMapByChrom = Util.FastQ_Util.getConsenusSeqMap(fastQFileName, dbtools)
    totalChroms = dbtools.operateDB("select","select count(*) from "+tablename)[0][0]
    currentchrID=dbtools.operateDB("select",sql+" limit 0,1")[0][0]
    seqMapByChrom[currentchrID]=""
    for i in range(0,totalChroms-1,20):
        currentsql=sql+" order by "+primaryID+" limit "+str(0)+",20"
        result=dbtools.operateDB("select",currentsql)
        for row in result:
            currentchrID=row[0]
            allchr+=currentchrID
            print(seqMapByChrom[currentchrID],file=outfile)
    outfile.close()
                