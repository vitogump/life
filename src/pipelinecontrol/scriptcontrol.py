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
parser.add_option("-c", "--cmdexample", dest="cmdtemplatefile",help="oneline scriptexamplefile")
# parser.add_option("-o", "--outputpath", dest="outputpath", help="outputpath")
parser.add_option("-d", "--datadepth", dest="datadepth", help="it's the depth of the dir from the inputdatapath which the data file that need to be process in it,the depth of the inputdatapath is 0")

parser.add_option("-s", "--scriptstorepath", dest="scriptstorepath", help="bam bai sam sorted.bam vcf blast and so on. note this is just used in the cmdline output parameter")
parser.add_option("-m", "--mode", dest="mode",
                  help="1 :means produce cmdline scripts for every terminal folder,the input data should be all the data files under the terminal folder. 2:use all selected data files as the input parameters in the only one cmdline script")
parser.add_option("-I","--Interceptor_depth",dest="Interceptor_depth",default="0",help="depth of the folder to output")
parser.add_option("-l", "--interceptdirs", dest="interceptdirs",action="append", default=[], help="winvalue or zvalue")
parser.add_option("-t","--dirsubtotaglen",dest="dirsubtotaglen",default=1,help="from dataupdir to -t step dir")
parser.add_option("-T","--directsubtag",dest="directsubtag",default=None,help="from dataupdir to -t step dir")
parser.add_option("-1","--collection_depth",dest="collection_depth",default="-1",help="depth of the folder to output")
parser.add_option("-2","--outfilepre",dest="outfilepre",default="mode2prename",help="depth of the folder to output")
parser.add_option("-q", "--quiet",
                  action="store_false", dest="verbose", default=True,
                  help="don't print status messages to stdout")
                                                                                                                                                          
(options, args) = parser.parse_args()
datadepth=int(options.datadepth)


# outputpath=options.outputpath
scriptsstoredir=options.scriptstorepath
if not os.path.exists(scriptsstoredir):
    os.makedirs(scriptsstoredir)
mode=int(options.mode)
Interceptor_depth=int(options.Interceptor_depth)



# outputpath=re.search(r"\${output=\s*([^\s^\|]*)\|suffix=(.*)}",scriptcmdline).group(1)
# outsuffix=re.search(r"\${output=\s*([^\s^\|]*)\|suffix=(.*)}",scriptcmdline).group(2)
# print(inputdatafilesrootpath,"outputpath=",outputpath,outsuffix)


interceptdirs=options.interceptdirs
dirsubtotaglen=int(options.dirsubtotaglen)
print(interceptdirs,"dirsubtotaglen:",str(dirsubtotaglen))

if __name__ == '__main__':
    if mode==1:
        if options.collection_depth!="-1":
            collection_depth=int(options.collection_depth)
            if collection_depth<Interceptor_depth:
                print("collection_depth<Interceptor_depth error")
#                 exit(-1)
        else:
            print("need -1 collection_depth")
            
        #progamma logic
        operatorwithdata_mode1=OperatorWithData_mode1(options.cmdtemplatefile,scriptsstoredir=scriptsstoredir,taglen=dirsubtotaglen,Dtag=options.directsubtag)
        if collection_depth<datadepth:
            operatorwithdata_mode1.datadepthequalcollectdepth=False
        print(collection_depth,datadepth,operatorwithdata_mode1.datadepthequalcollectdepth)
        upTodownTravelDir(operatorwithdata_mode1.inputdatapath,operatorwithdata_mode1,datadepth,Interceptor_depth,collection_depth=collection_depth,interceptdirs=interceptdirs,rootDirnotchange=operatorwithdata_mode1.inputdatapath,Interceptor_depth_notchange=Interceptor_depth)
        
    elif mode==2:
        outfilepre=options.outfilepre.strip() 
        operatorwithdata_mode2=OperatorWithData_mode2(options.cmdtemplatefile,scriptsstoredir,outfilepre)
        collection_depth=datadepth
        upTodownTravelDir(operatorwithdata_mode2.inputdatapath,operatorwithdata_mode2,datadepth,Interceptor_depth,collection_depth=collection_depth,interceptdirs=interceptdirs,rootDirnotchange=operatorwithdata_mode2.inputdatapath,Interceptor_depth_notchange=Interceptor_depth)
        #finalcmdline=re.sub(r"\${output}")
        finalcmdline=re.sub(operatorwithdata_mode2.outputfilenamewithoutoupfilesuffix+"."+operatorwithdata_mode2.outsuffix,operatorwithdata_mode2.outputfilenamewithoutoupfilesuffix+str(operatorwithdata_mode2.count)+"."+operatorwithdata_mode2.outsuffix,operatorwithdata_mode2.newcmdline)#add count info to outputfile
        finalcmdline=re.sub(r"[-\w\d]+[=\s]+\${.*?}"," ",finalcmdline)
        try:
            print(operatorwithdata_mode2.scriptinputdata[0:254]+"\n"+operatorwithdata_mode2.scriptoutputdata[0:-1]+"\n" +operatorwithdata_mode2.scriptcontext+finalcmdline,file=open(operatorwithdata_mode2.scriptsstoredir+"/"+re.search(r"[^/]*$",outfilepre).group(0)+str(operatorwithdata_mode2.count)+"_"+operatorwithdata_mode2.cmdtemplatefilename+"_"+operatorwithdata_mode2.suffixstr+"Script.sh",'a'))
        except FileNotFoundError:
            print(operatorwithdata_mode2.scriptinputdata[0:254]+"\n"+operatorwithdata_mode2.scriptoutputdata[0:-1]+"\n" +operatorwithdata_mode2.scriptcontext+finalcmdline,file=open(operatorwithdata_mode2.scriptsstoredir+"/"+re.search(r"[^/]*$",outfilepre).group(0)+str(operatorwithdata_mode2.count)+"_"+operatorwithdata_mode2.cmdtemplatefilename+"_"+operatorwithdata_mode2.suffixstr+"Script.sh",'w'))
    print("==============")
#     cmdline=operatorwithdata_mode1.cmdline

                
    #print(cmdline,finalcmdline)
       


