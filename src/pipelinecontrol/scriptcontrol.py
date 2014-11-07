# -*- coding: UTF-8 -*-
'''
Created on 2014-11-7

@author: liurui
'''
from optparse import OptionParser
import re,os

parser = OptionParser()

parser.add_option("-p", "--inputdatapath", dest="inputdatapath", help="output data name is defined as 'inputdatapath folder name'+'is subfolder name'+'is subfolder name'+...")
parser.add_option("-c", "--cmdexample", dest="cmdexamplefile",help="oneline scriptexamplefile")
parser.add_option("-o", "--outputpath", dest="outputpath", help="outputpath")
parser.add_option("-n", "--newdatastream", dest="newdatastream",action="store_true", help="this will determine the way to define the output file name")


parser.add_option("-s", "--suffix", dest="suffixname", default="", help=".bam .bai .sam .vcf .blast and so on")
# parser.add_option("-b", "--winfileName", dest="winfileName", help="winfileName ")
# parser.add_option("-w","--winWidth",dest="winWidth",default="40000",help="win width ")
# parser.add_option("-X", "--winType", dest="winType", default="zvalue", help="winvalue or zvalue")
parser.add_option("-q", "--quiet",
                  action="store_false", dest="verbose", default=True,
                  help="don't print status messages to stdout")
                                                                                                                                                          
(options, args) = parser.parse_args()
isNewdatastream=options.newdatastream
if options.inputdatapath[-1]=="/":
    inputdatapath=options.inputdatapath
else:
    inputdatapath=options.inputdatapath+"/"
outputpath=options.outputpath
outsuffix=options.suffixname

cmdline=open(options.cmdexamplefile,'r').readline()
print(cmdline)
subtargets=re.findall(r"\${.*?}",cmdline)

if isNewdatastream:
    print("new data orignal")

if __name__ == '__main__':
    pops=os.listdir(path=inputdatapath)
    for dirpath,dirnames,filenames in os.walk(inputdatapath):
        print(dirpath,dirnames,filenames)
#     print(pops,inputdatapath)
#     for a in pops:
#         popdir=inputdatapath+a
#         indvds=os.listdir(path=popdir)
#         for b in indvds:
#             if os.path.isdir(popdir+"/"+b):
#                 print(popdir+"/"+b)
#             else:
#                 print("file")

