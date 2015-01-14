# -*- coding: utf-8 -*- 
'''
Created on 2014-12-1

@author: liurui
'''
from optparse import OptionParser
import pickle
import re

import NGS.BasicUtil.DBManager as dbm
from NGS.Service.Ancestralallele import AncestralAlleletabletools
from src.NGS.BasicUtil import Util


parser = OptionParser()

#"output data name is defined as 'inputdatapath folder name'+'is subfolder name'+'is subfolder name'+..."

parser.add_option("-m", "--mode", dest="mode", help="mode")
parser.add_option("-c", "--chromtable", dest="chromtable", help="it's the depth of the dir from the inputdatapath which the data file that need to be process in it,the depth of the inputdatapath is 0")
parser.add_option("-t","--toplevelsnptable",dest="toplevelsnptable",default="0",help="depth of the folder to output")
#mode 1
parser.add_option("-a", "--outgroupvcffile", dest="outgroupvcffile",help="mode1")
parser.add_option("-d", "--depthfile", dest="depthfile", help="mode1")
parser.add_option("-n", "--outgroupname", dest="outgroupname", help="外群在depth文件中的名字（也就是拦截层文件夹名）mode1")
#mode 2
parser.add_option("-1", "--ancenstralref", dest="ancenstralref",help="ancenstralref fa file mode2")
parser.add_option("-2", "--ref", dest="ref",help="ref fa file mode2")
parser.add_option("-q", "--quiet",
                  action="store_false", dest="verbose", default=True,
                  help="don't print status messages to stdout")
                                                                                                                                                          
(options, args) = parser.parse_args()
archicpopVcfFile=options.outgroupvcffile.strip()
chromtable=options.chromtable
toplevelsnptablename=options.toplevelsnptable
archicpopNameindepthFile=options.outgroupname
depthFile=options.depthfile


flanklen=30
if __name__ == '__main__':
    ancestralalleletabletools=AncestralAlleletabletools(database="ninglabvariantdata", ip="10.2.48.140", usrname="root", pw="1234567",dbgenome="genomebasicinfo")
    if options.mode.strip()=="1":
        ancestralalleletabletools.fillAncestral(archicpopVcfFile=archicpopVcfFile,depthFile=depthFile,archicpopNameindepthFile=archicpopNameindepthFile,chromtable=chromtable.strip(),toplevelsnptablename=toplevelsnptablename)
    elif options.mode.strip()=="2":
        originalspeciesref=options.ancenstralref
        colname=re.search(r'[^/]*$',originalspeciesref).group(0)
        colname=re.sub(r"[^\w^\d]","_",colname);colname=colname[:10]
        ancestralalleletabletools.dbvariant.operateDB("callproc", "mysql_sp_add_column", data=(ancestralalleletabletools.dbname, toplevelsnptablename, colname, "char(128)", "default null"))            
        OUTFILENAME="ducksnpflankseq.fa"
        outfile=open("ducksnpflankseq.fa",'w')
        duckrefhandler=open(options.ref,'r')
        dbtoolsforchrom = dbm.DBTools("10.2.48.140", "root", "1234567", "genomebasicinfo")
        try:
            duckrefindex = pickle.load(open(options.ref + ".myindex", 'rb'))
            originalspeciesindex = pickle.load(open(originalspeciesref + ".myindex", 'rb'))
        except IOError:
            Util.generateIndexByChrom(options.ref, options.ref + ".myindex")
            Util.generateIndexByChrom(originalspeciesref, originalspeciesref + ".myindex")
            duckrefindex = pickle.load(open(options.reference + ".myindex", 'rb'))
            originalspeciesindex = pickle.load(open(originalspeciesref + ".myindex", 'rb'))
        totalChroms = dbtoolsforchrom.operateDB("select","select count(*) from "+chromtable)[0][0]
        for i in range(0,totalChroms,20):
            currentsql="select * from " + chromtable+" order by chrlength limit "+str(i)+",20"
            result=dbtoolsforchrom.operateDB("select",currentsql)
            for row in result:
                currentchrID=row[0]
                currentchrLen=int(row[2])
                ancestralalleletabletools.getflankseqs(currentchrID,currentchrLen, 1+flanklen, currentchrLen, idxedreffilehandler=duckrefhandler,ancestralgenomename=colname, refindex=duckrefindex, flanklen=flanklen,outfile=outfile, tablename=toplevelsnptablename)
        outfile.close()
        duckrefhandler.close()