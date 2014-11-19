# -*- coding: UTF-8 -*-
'''
Created on 2014-11-7

@author: liurui
'''
from optparse import OptionParser
import re, os

from src.pipelinecontrol.Util import OperatorWithData_mode1,upTodownTravelDir,OperatorWithData_mode2


parser = OptionParser()

#"output data name is defined as 'inputdatapath folder name'+'is subfolder name'+'is subfolder name'+..."
parser.add_option("-c", "--cmdexample", dest="cmdexamplefile",help="oneline scriptexamplefile")
# parser.add_option("-o", "--outputpath", dest="outputpath", help="outputpath")
parser.add_option("-d", "--datadepth", dest="datadepth", help="it's the depth of the dir from the inputdatapath which the data file that need to be process in it,the depth of the inputdatapath is 0")


parser.add_option("-s", "--scriptstorepath", dest="scriptstorepath", help="bam bai sam sorted.bam vcf blast and so on. note this is just used in the cmdline output parameter")
parser.add_option("-m", "--mode", dest="mode",
                  help="1 :means produce cmdline scripts for every terminal folder,the input data should be all the data files under the terminal folder. 2:use all selected data files as the input parameters in the only one cmdline script")
parser.add_option("-I","--Interceptor_depth",dest="Interceptor_depth",default="0",help="depth of the folder to output")
parser.add_option("-l", "--interceptdirs", dest="interceptdirs",action="append", default=[], help="winvalue or zvalue")
parser.add_option("-q", "--quiet",
                  action="store_false", dest="verbose", default=True,
                  help="don't print status messages to stdout")
                                                                                                                                                          
(options, args) = parser.parse_args()
datadepth=int(options.datadepth)


# outputpath=options.outputpath
scriptsstoredir=options.scriptstorepath
mode=int(options.mode)
Interceptor_depth=int(options.Interceptor_depth)

scriptcontent=open(options.cmdexamplefile,'r').read()

scriptcontext=re.search(r"([\s\S]*(\n)*)cmdline=.*",scriptcontent).group(1)

inputdatafilesrootpath=re.search(r"(\n)*inputdatafilesrootpath=\s*(.*)",scriptcontext).group(2)
scriptcmdline=re.search(r"(.*(\n)*)cmdline=\s*(.*)",scriptcontent).group(3)
print(scriptcontent,scriptcontext,inputdatafilesrootpath,scriptcmdline,sep="\n")

# outputpath=re.search(r"\${output=\s*([^\s^\|]*)\|suffix=(.*)}",scriptcmdline).group(1)
# outsuffix=re.search(r"\${output=\s*([^\s^\|]*)\|suffix=(.*)}",scriptcmdline).group(2)
# print(inputdatafilesrootpath,"outputpath=",outputpath,outsuffix)


interceptdirs=options.interceptdirs
print(interceptdirs)

if __name__ == '__main__':
    if mode==1:

        #progamma logic
        operatorwithdata_mode1=OperatorWithData_mode1(scriptcmdline,inputdatafilesrootpath,scriptcontext=scriptcontext,scriptsstoredir=scriptsstoredir,interceptdirs=interceptdirs)
        upTodownTravelDir(inputdatafilesrootpath,operatorwithdata_mode1,datadepth,Interceptor_depth)
        
    elif mode==2:
            
        operatorwithdata_mode2=OperatorWithData_mode2(scriptcmdline,outputpath,outsuffix,scriptsstoredir,interceptdirs)
        upTodownTravelDir(inputdatafilesrootpath,operatorwithdata_mode2,datadepth,Interceptor_depth)
        #finalcmdline=re.sub(r"\${output}")
        finalcmdline=re.sub(r"[-\w\d]+[=\s]+\${.*?}"," ",operatorwithdata_mode2.newcmdline)
        try:
            print(scriptcontext+finalcmdline,file=open("F:/work/pipelinecontrol/scripts/"+operatorwithdata_mode2.outsuffix+"Script.sh",'a'))
        except FileNotFoundError:
            print(scriptcontext+finalcmdline,file=open("F:/work/pipelinecontrol/scripts/"+operatorwithdata_mode2.outsuffix+"Script.sh",'w'))
    print("==============")
#     cmdline=operatorwithdata_mode1.cmdline

                
    #print(cmdline,finalcmdline)
       


