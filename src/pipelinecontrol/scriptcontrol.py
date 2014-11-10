# -*- coding: UTF-8 -*-
'''
Created on 2014-11-7

@author: liurui
'''
from optparse import OptionParser
import re, os

from src.pipelinecontrol.Util import OperatorWithData_mode1,upTodownTravelDir,OperatorWithData_mode2


parser = OptionParser()

parser.add_option("-p", "--inputdatapath", dest="inputdatapath", help="output data name is defined as 'inputdatapath folder name'+'is subfolder name'+'is subfolder name'+...")
parser.add_option("-c", "--cmdexample", dest="cmdexamplefile",help="oneline scriptexamplefile")
parser.add_option("-o", "--outputpath", dest="outputpath", help="outputpath")
parser.add_option("-d", "--datadepth", dest="datadepth", help="it's the depth of the dir from the inputdatapath which the data file that need to be process in it,the depth of the inputdatapath is 0")


parser.add_option("-s", "--suffix", dest="suffixname", default="", help="bam bai sam sorted.bam vcf blast and so on. note this is just used in the cmdline output parameter")
parser.add_option("-m", "--mode", dest="mode",
                  help="1 :means produce cmdline scripts for every terminal folder,the input data should be all the data files under the terminal folder. 2:use all selected data files as the input parameters in the only one cmdline script")
parser.add_option("-1","--depthoffoldertocopy",dest="depthoffoldertocopy",default="0",help="0 means don't creat folder in the output folder")
parser.add_option("-2", "--interceptdirs", dest="interceptdirs", default=[],action="append", help="winvalue or zvalue")
parser.add_option("-q", "--quiet",
                  action="store_false", dest="verbose", default=True,
                  help="don't print status messages to stdout")
                                                                                                                                                          
(options, args) = parser.parse_args()
datadepth=int(options.datadepth)
inputdatapath=options.inputdatapath

outputpath=options.outputpath
outsuffix=options.suffixname
mode=int(options.mode)
n_subdirs=int(options.depthoffoldertocopy)

cmdline=open(options.cmdexamplefile,'r').readline()

mode2_Interceptor=options.interceptdirs
print(mode2_Interceptor)

if __name__ == '__main__':
    if mode==1:
        #friendly ui interaction
        if options.depthoffoldertocopy==0:
            print("option -1 default 0")
        if len(mode2_Interceptor)>1:
            print("warning: option -2 is useless in mode 1")
        #progamma logic
        operatorwithdata_mode1=OperatorWithData_mode1(cmdline,outputpath,outsuffix,inputdatapath,n_subdirs=n_subdirs)
        upTodownTravelDir(inputdatapath,OperatorWithData=operatorwithdata_mode1,datadepth=datadepth)
        
    elif mode==2:
        if options.depthoffoldertocopy!=0:
            print("warning: option -1 is not used in mode 2")
            
        operatorwithdata_mode2=OperatorWithData_mode2(cmdline,outputpath,outsuffix)
        upTodownTravelDir(inputdatapath,OperatorWithData=operatorwithdata_mode2,datadepth=datadepth)
        #finalcmdline=re.sub(r"\${output}")
        finalcmdline=re.sub(r"[-\w\d]+[=\s]+\${.*?}"," ",operatorwithdata_mode2.newcmdline)
        print(finalcmdline,file=open("F:/work/pipelinecontrol/scripts/"+outsuffix+"_script.sh",'a'))
    print("==============")
#     cmdline=operatorwithdata_mode1.cmdline

                
    #print(cmdline,finalcmdline)
       


