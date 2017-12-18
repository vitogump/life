# -*- coding: UTF-8 -*-
'''
Created on 2014-12-1

@author: liurui
'''
from optparse import OptionParser
import pickle, re, sys,pysam

from NGS.BasicUtil import VCFutil
import NGS.BasicUtil.DBManager as dbm
from NGS.BasicUtil import Util
from src.NGS.Service.Ancestralallele import AncestralAlleletabletools


parser = OptionParser()

#"output data name is defined as 'inputdatapath folder name'+'is subfolder name'+'is subfolder name'+..."

parser.add_option("-m", "--mode", dest="mode", help="mode 1 according vcf,mode 2 ref compare")
parser.add_option("-c", "--chromlistfilename", dest="chromlistfilename", help="it's the depth of the dir from the inputdatapath which the data file that need to be process in it,the depth of the inputdatapath is 0")
parser.add_option("-t","--toplevelsnptable",dest="toplevelsnptable",default="0",help="depth of the folder to output")
#mode 1
parser.add_option("-C", "--vcfbamconfig", dest="vcfbamconfig",help="outgroupvcffile to make sure ancestral allele")

# parser.add_option("-n", "--outgroupname", dest="outgroupname", help="the name of outgroup in depth file,also the folder name of intercept level ,mode1")
#mode 2
parser.add_option("-1", "--ancenstralref", dest="ancenstralref",help="ancenstralref fa file mode2")
parser.add_option("-2", "--ref", dest="ref",help="ref fa file mode2")
parser.add_option("-f", "--flanklen", dest="flanklen",default='50',help="ref fa file mode2")
parser.add_option("-q", "--quiet",
                  action="store_false", dest="verbose", default=True,
                  help="don't print status messages to stdout")

(options,args)=parser.parse_args()
vcfnameKEY_vcfobj_pyBAMfilesVALUE={}
if options.ancenstralref==None:
    
    archicpopvcfbamconfig=options.vcfbamconfig.strip()
    vcfconfig=open(archicpopvcfbamconfig,"r")
    for line in vcfconfig:
        vcffilename_obj=re.search(r"vcffilename=(.*)",line.strip())
        if vcffilename_obj!=None:
            vcfname=vcffilename_obj.group(1).strip()
            vcfnameKEY_vcfobj_pyBAMfilesVALUE[vcfname]=[]
            vcfnameKEY_vcfobj_pyBAMfilesVALUE[vcfname].append(VCFutil.VCF_Data(vcfname))
        elif line.split():
            vcfnameKEY_vcfobj_pyBAMfilesVALUE[vcfname].append(pysam.Samfile(line.strip(),'rb'))
    vcfconfig.close()
toplevelsnptablename=options.toplevelsnptable




flanklen=int(options.flanklen.strip())
if __name__ == '__main__':
    ancestralalleletabletools=AncestralAlleletabletools(database=Util.vcfdbname, ip=Util.ip, usrname=Util.username, pw=Util.password,dbgenome=Util.genomeinfodbname)
    if options.mode.strip()=="1":
#         depthFile=options.depthfile
        chromlist=[]
        chromlistfile=open(options.chromlistfilename,"r")
        for chrrow in chromlistfile:
            chrrowlist=re.split(r'\s+',chrrow.strip())
            chromlist.append(chrrowlist[0].strip())
        ancestralalleletabletools.fillAncestral(archicpopVcfFile=vcfnameKEY_vcfobj_pyBAMfilesVALUE,chromlist=chromlist,toplevelsnptablename=toplevelsnptablename)
    elif options.mode.strip()=="2":
#         if options.depthfile!=None:
#             print(options.depthfile,"no need")
#         originalspeciesref=options.ancenstralref
#         colname=re.search(r'[^/]*$',originalspeciesref).group(0)
#         colname=re.sub(r"[^\w^\d]","_",colname);colname=colname[:10]
#         print(colname)
#         ancestralalleletabletools.dbvariant.operateDB("callproc", "mysql_sp_add_column", data=(ancestralalleletabletools.dbvariant, toplevelsnptablename, colname, "char(128)", "default null"))
        OUTFILENAME="ducksnpflankseq.fa"
        outfile=open(options.chromlistfilename+"snpflankseq.fa",'w')
        duckrefhandler=open(options.ref,'r')
        try:
            duckrefindex = pickle.load(open(options.ref + ".myfasteridx", 'rb'))
#             originalspeciesindex = pickle.load(open(originalspeciesref + ".myindex", 'rb'))
        except IOError:
            Util.generateFasterRefIndex(options.ref, options.ref + ".myfasteridx")
            duckrefindex = pickle.load(open(options.ref + ".myfasteridx", 'rb'))
            
#         try:
#             originalspeciesindex = pickle.load(open(originalspeciesref + ".myindex", 'rb'))
#         except IOError:
#             Util.generateIndexByChrom(originalspeciesref, originalspeciesref + ".myindex")
#             originalspeciesindex = pickle.load(open(originalspeciesref + ".myindex", 'rb'))
        chrom_lenlist=[]
        chromlistfile=open(options.chromlistfilename,"r")
        for chrrow in chromlistfile:
            chrrowlist=re.split(r'\s+',chrrow.strip())
            chrom_lenlist.append((chrrowlist[0].strip(),int(chrrowlist[1].strip())))
        for currentchrID,currentchrLen in chrom_lenlist:
            ancestralalleletabletools.getflankseqs(currentchrID,currentchrLen, flanklen, currentchrLen, idxedreffilehandler=duckrefhandler,ancestralgenomenameaddtotable=None, refindex=duckrefindex, flanklen=flanklen,outfile=outfile, tablename=toplevelsnptablename)
        outfile.close()
        duckrefhandler.close()
    chromlistfile.close()
    print("finished")