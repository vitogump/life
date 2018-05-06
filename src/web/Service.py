'''
Created on 2014-11-17

@author: liurui
'''
import time,datetime
import re,string,os,random,markdown2
import src.web.dba as mydba
from src.pipelinecontrol.Util import OperatorWithData_webservice, upTodownTravelDir
from sqlalchemy.orm import session
from tabulate import tabulate
import src.web.dba as aaa
from src.web import entity
from src.web.dba import addJobs2jobstate

SLEEP_FOR_NEXT_TRY=5
def jobminitor(currentUstr):
#     
    session=aaa.getWebSession()
    l=[]
    whileStart=time.clock()
    while not l:
        whileEnd=time.clock ()
        if not currentUstr or int(whileEnd-whileStart)>=18:
            time.sleep(SLEEP_FOR_NEXT_TRY)
            l = session.query(entity.Jobs_recoder).all()
            break
        l=session.query(entity.Jobs_recoder).filter(entity.Jobs_recoder.foldername.like("%"+currentUstr+"%")).all()
        
    header=["*scriptname*","*scriptfolder*","*outputdata*","*starttime*","*finishtime*"," *state*","*outputinfo*"]
    mylist=[]
    
    for i in l:
        print("sssssssssssssssss",i.outputinfo)
#             mylist.append([("<br>"+i.scriptname),i.foldername[12:],("&nbsp;"+str(i.startdate)+"&nbsp;"),("&nbsp;"+str(i.finishdate)+"&nbsp;"),("&nbsp;"+str(i.state)),("""<input type="button" value="outputinfo" onclick="location.href='http://www.baidu.com'">""")])
        mylist.append([i.scriptname,i.foldername[12:],("&nbsp;"+str(i.outputdata)+"&nbsp;"),("&nbsp;"+str(i.startdate)+"&nbsp;"),("&nbsp;"+str(i.finishdate)+"&nbsp;"),("&nbsp;"+str(i.state)),("""<input type="button" value="outputinfo" onclick="location.href='http://www.baidu.com'">""")])
    print(mylist)
    print(header)
    print("======orgtbl=====================")
    text=tabulate(mylist,header,tablefmt="orgtbl")
    text=re.sub('\|\|[\-\+]+\|\|\n', '', text.replace("|", "||"))
    print(text)

    html=markdown2.markdown(text,extras=["wiki-tables"])
    html=html.replace("<tr", "<tr bgcolor='lightgrey'",1)

    
    return html
def random_uniqScriptDir(scriptspath,randomlength=8):
    a = list(string.ascii_letters)
    random.shuffle(a)
    ranUniscriptspath=(scriptspath.rstrip("/")+"/"+''.join(a[:randomlength]))
    if  os.path.exists(ranUniscriptspath):
        while True:
            random.shuffle(a)
            if  ''.join(a[:randomlength]) not in os.listdir(scriptspath): 
                ranUniscriptspath=(scriptspath.rstrip("/")+"/"+''.join(a[:randomlength]))
                break
    return ranUniscriptspath
def scriptproduce(datadepth,collectiondepth,scriptspath,inputdatapath,softwareconfig,parametersStr,inputList,outputList,lenOfdirtotag=0,taglist=[],selecteddepth=0,selecteddirs=[]):#selecteddepth=0 means check collectiondepth only
    
    inputstr=(" "+" ".join(taglist)+" ") if int(lenOfdirtotag)!=0 else ""
    inputstr+=" ".join(inputList)
    
    parametersStr,N=re.subn(r"\$\$\$\$",inputstr,parametersStr)
    if not os.path.exists(scriptspath):
        os.makedirs(scriptspath)
    outputStr=outputList[1]+" ${output="+outputList[0]+"|suffix="+outputList[2]+"}"
    print(datadepth,collectiondepth,scriptspath,inputdatapath,softwareconfig)
    cmdline=softwareconfig+" "+parametersStr+" "+outputStr
    print(cmdline)
    ranUniscriptspath=random_uniqScriptDir(scriptspath)
    os.makedirs(ranUniscriptspath)    
    operatorwithdata=OperatorWithData_webservice(inputdatapath,cmdline,ranUniscriptspath,taglen=lenOfdirtotag)
    operatorwithdata.cmdtemplatefilename=re.split(r'\s+',softwareconfig.strip())[0]+"Get"+outputList[2]
    if int(selecteddepth)==0:
        selecteddirs=[]
    upTodownTravelDir(inputdatapath,operatorwithdata,int(datadepth),int(selecteddepth),collection_depth=int(collectiondepth),interceptdirs=selecteddirs,rootDirnotchange=operatorwithdata.inputdatapath,Interceptor_depth_notchange=int(selecteddepth))
    return ranUniscriptspath