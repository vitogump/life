'''
Created on 2015-5-14

@author: liurui
'''
from optparse import OptionParser
from src.NGS.BasicUtil import Util
import src.NGS.BasicUtil.DBManager as dbm

extendlen=100000
parser = OptionParser()
parser.add_option("-m", "--minlength", dest="minlength")
parser.add_option("-g", "--gtffile", dest="gtffile",default=None, help="gtffile")
(options, args) = parser.parse_args()
if options.gtffile !=None:
    gtfMap,utrMap,allgeneSetMap = Util.getGtfMap(options.gtffile)
minlength = options.minlength
chromtable = Util.pekingduckchromtable#options.chromtablename
dbchromtools = dbm.DBTools(Util.ip, Util.username, Util.password, Util.genomeinfodbname)
sql = "select * from " + chromtable + " where chrlength>=" + minlength
primaryID = "chrID"
if __name__ == '__main__':
    totalChroms = dbchromtools.operateDB("select", "select count(*) from " + chromtable + " where chrlength>=" + minlength)[0][0]
    innergenecollectRegion={}
    for i in range(0, totalChroms, 20):
        currentsql = sql + " order by " + primaryID + " limit " + str(i) + ",20"
        result = dbchromtools.operateDB("select", currentsql)
        for row in result:
            seektuple = ()
            currentchrID = row[0]
            currentchrLen = int(row[1])
            if options.gtffile !=None:
                if currentchrID not in gtfMap:
                    continue
                GeneGrouplist = Util.getGeneGrouplist(gtfMap[currentchrID])
                for geneGroup in GeneGrouplist:
                    grouplist_idx=GeneGrouplist.index(geneGroup)
                    if len(geneGroup)==2:
                        
                        cdsidx = 3
                        for feature, elemStart, elemEnd, frame in geneGroup[1][4:]:
                            if feature == "CDS" or feature == "stop_codon":
                                cdsidx += 1
                                if cdsidx==4:
                                    pass
                                else:
                                    if geneGroup[1][cdsidx-1][2]+extendlen<geneGroup[1][cdsidx][1]-extendlen:
                                        innergenecollectRegion="suspend"
                