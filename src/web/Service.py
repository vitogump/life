'''
Created on 2014-11-17

@author: liurui
'''
import markdown2,time,datetime
import re,string,os,random
import src.web.dba as mydba
from src.pipelinecontrol.Util import OperatorWithData_webservice, upTodownTravelDir



def scriptproduce(datadepth,collectiondepth,scriptspath,inputdatapath,softwareconfig,parametersStr,inputList,outputList,lenOfdirtotag=1,taglist=[],selecteddepth=0,selecteddirs=[]):#selecteddepth=0 means check collectiondepth only
    inputstr=(" "+" ".join(taglist)+" ")
    inputstr+=" ".join(inputList)
    
    parametersStr,N=re.subn(r"\$\$\$\$",inputstr,parametersStr)
    if not os.path.exists(scriptspath):
        os.makedirs(scriptspath)
    outputStr=outputList[1]+" ${output="+outputList[0]+"|suffix="+outputList[2]+"}"
    print(datadepth,collectiondepth,scriptspath,inputdatapath,softwareconfig)
    cmdline=softwareconfig+" "+parametersStr+" "+outputStr
    print(cmdline)
    operatorwithdata=OperatorWithData_webservice(inputdatapath,cmdline,scriptspath,taglen=lenOfdirtotag)
    operatorwithdata.cmdtemplatefilename=softwareconfig+"Get"+outputList[2]
    upTodownTravelDir(inputdatapath,operatorwithdata,int(datadepth),int(selecteddepth),collection_depth=int(collectiondepth),interceptdirs=selecteddirs,rootDirnotchange=operatorwithdata.inputdatapath,Interceptor_depth_notchange=int(selecteddepth))
    return operatorwithdata.scriptsstorediruniq