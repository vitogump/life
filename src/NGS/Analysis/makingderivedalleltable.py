'''
Created on 2014-4-25

@author: liurui
'''
import NGS.BasicUtil.DerivedalleleProcessor as DAP
from optparse import OptionParser

parser = OptionParser()
parser.add_option("-d", "--dbname", dest="dbname",# action="callback",type="string",callback=useoptionvalue_previous1,
                  help="write report to FILE")
parser.add_option("-v", "--vcffile", dest="vcffile",# action="callback",type="string",callback=useoptionvalue_previous2,
                  help="write report to FILE")
# (options, args) = parser.parse_args()
parser.add_option("-D","--Depthfile",dest="Depthfile",help="default infile1_infile2")#
parser.add_option("-s","--slidesize",dest="slidesize",help="default infile2_infile1")#
parser.add_option("-q", "--quiet",
                  action="store_false", dest="verbose", default=True,
                  help="don't print status messages to stdout")
(options, args) = parser.parse_args()
vcfFileName=options.vcffile
DepthFileName=options.Depthfile
if __name__ == '__main__':
    aaa=DAP.MakeDerivedAlleletable(database="life_pilot",ip="localhost",usrname="root",pw="1234567")
    aaa.createtable()
    aaa.filldata(vcfFileName=vcfFileName,depthfileName=DepthFileName)