'''
Created on 2014-12-1

@author: liurui
'''
from optparse import OptionParser
import re
from NGS.BasicUtil import Util
from NGS.Service.Ancestralallele import AncestralAlleletabletools


parser = OptionParser()

#"output data name is defined as 'inputdatapath folder name'+'is subfolder name'+'is subfolder name'+..."
parser.add_option("-v", "--vcftablelist", dest="vcftablelist",action="append",default=[],help="")
# parser.add_option("-o", "--outputpath", dest="outputpath", help="outputpath")
parser.add_option("-c", "--chromlistfilename", dest="chromlistfilename", help="it's the depth of the dir from the inputdatapath which the data file that need to be process in it,the depth of the inputdatapath is 0")

parser.add_option("-t","--toplevelsnptable",dest="toplevelsnptable",default="ducksnp_toplevel",help="depth of the folder to output")
parser.add_option("-o","--outputtablename",dest="outputtablename",default="duckout",help="depth of the folder to output")

parser.add_option("-d", "--quiet",action="append", dest="depthfilenames",default=None)
                                                                                                                                                          
(options, args) = parser.parse_args()
vcftablenames=options.vcftablelist
print(vcftablenames)
toplevelsnptable=options.toplevelsnptable
outtable_filename=options.outputtablename
depthfilenames=options.depthfilenames
chromlist=[]
chromlistfile=open(options.chromlistfilename,"r")
for chrrow in chromlistfile:
    chrrowlist=re.split(r'\s+',chrrow.strip())
    chromlist.append(chrrowlist[0].strip())
print(chromlistfile,outtable_filename,toplevelsnptable,depthfilenames)
if __name__ == '__main__':
    atools=AncestralAlleletabletools(database=Util.vcfdbname, ip=Util.ip, usrname=Util.username, pw=Util.password,dbgenome=Util.genomeinfodbname)
    atools.leftjoinSelectedTables(chromlist,outtable_filename,vcftablenames,toplevelsnptable,depthfilenames=depthfilenames)