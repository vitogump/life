# -*- coding: UTF-8 -*-
'''
Created on 2014-11-7

@author: liurui
'''
from optparse import OptionParser
import re, os

from pipelinecontrol.Util import OperatorWithData_mode1, upTodownTravelDir, \
    OperatorWithData_mode2


parser = OptionParser()

parser.add_option("-p", "--inputdatapath", dest="inputdatapath", help="output data name is defined as 'inputdatapath folder name'+'is subfolder name'+'is subfolder name'+...")
parser.add_option("-c", "--cmdexample", dest="cmdexamplefile",help="oneline scriptexamplefile")
parser.add_option("-o", "--outputpath", dest="outputpath", help="outputpath")
parser.add_option("-d", "--datadepth", dest="datadepth", help="data need to be process in the depth from the inputdatapath,the depth of the inputdatapath is 0")


parser.add_option("-s", "--suffix", dest="suffixname", default="", help=".bam .bai .sam .vcf .blast and so on. note this is just used in the cmdline output parameter")
parser.add_option("-m", "--mode", dest="mode",
                  help="1 :means produce cmdline scripts for every terminal folder,the input data should be all the data files under the terminal folder. 2:use all selected data files as the input parameters in the only one cmdline script")
# parser.add_option("-w","--winWidth",dest="winWidth",default="40000",help="win width ")
# parser.add_option("-X", "--winType", dest="winType", default="zvalue", help="winvalue or zvalue")
parser.add_option("-q", "--quiet",
                  action="store_false", dest="verbose", default=True,
                  help="don't print status messages to stdout")
                                                                                                                                                          
(options, args) = parser.parse_args()
datadepth=int(options.datadepth)
inputdatapath=options.inputdatapath

outputpath=options.outputpath
outsuffix=options.suffixname
mode=options.mode

cmdline=open(options.cmdexamplefile,'r').readline()
print(cmdline)
subtargets=re.findall(r"\${.*?}",cmdline)
outputidx=subtargets.index("${output}")

dataInDirMap={}
targetdatasuffix=[]

for target in subtargets:
    c=re.search(r'\${(.*?)}',target).group(1)
    if "output"!=c:
        targetdatasuffix.append(c)
        
        
print("targetdatasuffix",targetdatasuffix)


operatorwithdata_mode1=OperatorWithData_mode1(cmdline)
operatorwithdata_mode2=OperatorWithData_mode2(cmdline)
if __name__ == '__main__':
    if mode==1:
        upTodownTravelDir(inputdatapath,OperatorWithData=operatorwithdata_mode1,maxdepth=datadepth)
        
    elif mode==2:
        upTodownTravelDir(inputdatapath,OperatorWithData=operatorwithdata_mode2,maxdepth=datadepth)
        finalcmdline=re.sub(r"-\wd+\s+\${.*?}"," ",operatorwithdata_mode2.cmdline)
    print("==============")
    cmdline=operatorwithdata_mode1.cmdline
#     popdirs=os.listdir(path=inputdatapath)
#     for a in popdirs:
#         if datadepth==3 and os.path.isdir(inputdatapath+"/"+a):
#             dataInDirMap[inputdatapath+"/"+a]={}
#             indvddirs=os.listdir(path=inputdatapath+"/"+a)
#             for b in indvddirs:
#                 if os.path.isdir(inputdatapath+"/"+a+"/"+b)==False:
#                     print(b+" :shouldn't accour here")
#                     continue
#                 dataInDirMap[inputdatapath+"/"+a][b]=[]
#                 datafiles=os.listdir(path=inputdatapath+"/"+a+"/"+b)
#                 print(inputdatapath+"/"+a+"/"+b,"datafiles",datafiles) 
#                 #ordered by the ${1} ${2} ... in cmdline
#                 for i in range(0,len(targetdatasuffix)):
#                     if targetdatasuffix[i]=="output":
#                         continue
#                     for datafilename in datafiles:
#                         #print("targetdatasuffix ",str(i),targetdatasuffix[i],datafilename,subtargets[i],re.search(r".*?"+targetdatasuffix[i],datafilename),".*?"+targetdatasuffix[i])
#                         if re.search(r".*?"+targetdatasuffix[i],datafilename)!=None:
#                             # sub was acted from the first to the rear most
#                             cmdline=re.sub(r"\${\s*"+targetdatasuffix[i]+"\s*}"," "+datafilename+" ",cmdline)
#         elif datadepth==2 and os.path.isdir(inputdatapath+"/"+a):
#             pass
#         elif datadepth==1:
#             datafiles=os.listdir(path=inputdatapath)
#             for datafilename in datafiles:
#                 if re.search(r".*?"+targetdatasuffix[0],datafilename)!=None:
#                     cmdline=re.sub(r"\${\s*"+targetdatasuffix+"\s*}"," "+datafilename+" ",cmdline)
#                     print(cmdline)
                    
                
    print(cmdline)
    print(dataInDirMap)        


