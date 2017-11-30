'''
Created on 2014-4-25

@author: liurui
'''
from NGS.BasicUtil import Util
from optparse import OptionParser
import NGS.BasicUtil.DBManager as dbm
import NGS.BasicUtil.DerivedalleleProcessor as DAP
# import NGS.RUtil.Make_Picture as MP

import pickle
"""
    the vcf file,every sample is a pop that has the same name in depthfile's title
"""
parser = OptionParser()
parser.add_option("-d", "--dbname", dest="dbname",# action="callback",type="string",callback=useoptionvalue_previous1,
                  help="write report to FILE")
parser.add_option("-1", "--mutiplepopsvcffile", dest="mutiplepopsvcffile",# action="callback",type="string",callback=useoptionvalue_previous2,
                  help="write report to FILE")
parser.add_option("-2", "--archicpopvcf", dest="archicpopvcf",# action="callback",type="string",callback=useoptionvalue_previous2,
                  help="write report to FILE")
parser.add_option("-3","--continuechrom",dest="continuechrom")
parser.add_option("-4","--continuepos",dest="continuepos")
# (options, args) = parser.parse_args()
parser.add_option("-D","--Depthfile",dest="Depthfile",help="default infile1_infile2")#
parser.add_option("-a","--ancenstryref",dest="ancenstryref")
parser.add_option("-r","--reference",dest="reference")
parser.add_option("-c","--chromtable",dest="chromtable")
parser.add_option("-f","--flanklen",dest="flanklen")
parser.add_option("-B","--pathtoblastn",dest="pathtoblastn")
parser.add_option("-b","--pathtoblastdb",dest="pathtoblastdb")
parser.add_option("-m","--minlength",dest="minlength")
parser.add_option("-p","--prefilename",dest="prefilename")
parser.add_option("-o","--ancenstryvcftable",dest="ancenstryvcftable")
parser.add_option("-q", "--quiet",
                  action="store_false", dest="verbose", default=True,
                  help="don't print status messages to stdout")
(options, args) = parser.parse_args()
allpop=args[:]
minlengthOfchrom=options.minlength
vcfFileName=options.mutiplepopsvcffile
DepthFileName=options.Depthfile
duckrefhandler=open(options.reference,'r')
originalspeciesref=options.ancenstryref
chromtable=options.chromtable
dbname=options.dbname
flanklen=int(options.flanklen)
pathtoblastn=options.pathtoblastn
pathtoblastdb=options.pathtoblastdb
ancestralsnptable=options.ancenstryvcftable
archicpopNameindepthFile="fanya"
archicpopVcfFile=options.archicpopvcf
if options.continuechrom!=None and options.continuepos!=None:
    continuechrom=options.continuechrom
    continuepos=int(options.continuepos)
else:
    continuechrom=None;continuepos=None
chromstable=options.chromtable
primaryID = "chrID"
OUTFILENAME="ducksnpflankseq.fa"
# outfile=open("ducksnpflankseq.fa",'w')
BlastOutFile="ducksnpflankseq.blast"
if __name__ == '__main__':
    aaa=DAP.MakeDerivedAlleletable(database=dbname,ip="10.2.48.96",usrname="root",pw="1234567")
#     aaa=DAP.MakeDerivedAlleletable(database=dbname,ip="10.2.48.140",usrname="root",pw="1234567")
#     ddd=MP.Dstistics_allpop(allpop)
#     ddd.caculateDofAllpossibleCombination(database=dbname,ip="10.2.48.140",usrname="root",pw="1234567", allpopssnptable="derived_alle_ref", chromstable=chromstable, winwidth=None, minlengthOfchrom=minlengthOfchrom, filenamepre=options.prefilename)
    dbtoolsforchrom = dbm.DBTools(Util.ip, Util.username, Util.password, Util.genomeinfodbname)
    try:
        duckrefindex = pickle.load(open(options.reference + ".myindex", 'rb'))
        originalspeciesindex = pickle.load(open(originalspeciesref + ".myindex", 'rb'))
    except IOError:
        Util.generateIndexByChrom(options.reference, options.reference + ".myindex")
        Util.generateIndexByChrom(originalspeciesref, originalspeciesref + ".myindex")
        duckrefindex = pickle.load(open(options.reference + ".myindex", 'rb'))
        originalspeciesindex = pickle.load(open(originalspeciesref + ".myindex", 'rb'))
#     aaa.createtable()
#     aaa.filldata(vcfFileName=vcfFileName,depthfileName=DepthFileName,continuechrom=continuechrom,continuepos=continuepos)
    aaa.fillarchicpop(archicpopVcfFile,DepthFileName,chromstable,archicpopNameindepthFile)
#     totalChroms = dbtoolsforchrom.operateDB("select","select count(*) from "+chromstable)[0][0]
#     for i in range(0,totalChroms,20):
#         currentsql="select * from " + chromstable+" order by chrlength limit "+str(i)+",20"
        result=dbtoolsforchrom.operateDB("select",currentsql)
#         for row in result:
#             currentchrID=row[0]
#             currentchrLen=int(row[2])
#             aaa.getflankseqs(currentchrID,currentchrLen, 1+flanklen, currentchrLen, idxedreffilehandler=duckrefhandler, refindex=duckrefindex, flanklen=flanklen,outfile=outfile, tablename="derived_alle_ref")
#     outfile.close()
#     duckrefhandler.close()
#     aaa.callblast(pathtoblastn,pathtoblastdb,OUTFILENAME,BlastOutFile)
#     aaa.extarctAncestryAlleleFromBlastOut(BlastOutFile,originalspeciesref,originalspeciesindex,tablename="derived_alle_ref",ancestralsnptable=ancestralsnptable)
    