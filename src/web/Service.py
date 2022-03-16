'''
Created on 2014-11-17

@author: liurui
'''
import time,datetime,sys#,cv2
import re,string,os,random,markdown2
import src.web.dba as mydba
from src.pipelinecontrol.Util import OperatorWithData_webservice, upTodownTravelDir,longestCommonPrefix
from sqlalchemy.orm import session
from tabulate import tabulate
import src.web.dba as aaa
from src.web import entity
from src.web.dba import addJobs2jobstate

SLEEP_FOR_NEXT_TRY=5
scriptdir=entity.scriptdir 
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
        l=session.query(entity.Jobs_recoder).filter(entity.Jobs_recoder.foldername.like("%"+currentUstr[:-9]+"/"+currentUstr[-8:]+"%")).all()
        
    header=["*scriptname*","*scriptfolder*","*outputdata*","*starttime*","*finishtime*"," *state*","*outputinfo*"]
    mylist=[]
    
    for i in l:
        print("sssssssssssssssss",i.outputinfo)
#             mylist.append([("<br>"+i.scriptname),i.foldername[12:],("&nbsp;"+str(i.startdate)+"&nbsp;"),("&nbsp;"+str(i.finishdate)+"&nbsp;"),("&nbsp;"+str(i.state)),("""<input type="button" value="outputinfo" onclick="location.href='http://www.baidu.com'">""")])
        mylist.append([i.scriptname,i.foldername[12:],("&nbsp;"+str(i.outputdata)+"&nbsp;"),("&nbsp;"+str(i.startdate)+"&nbsp;"),("&nbsp;"+str(i.finishdate)+"&nbsp;"),("&nbsp;"+str(i.state)),("""<input type="button" value="outputinfo" onclick="location.href='http://www.baidu.com'">""")])
#     print(mylist)
    print(header)
    print("======orgtbl=====================")
    text=tabulate(mylist,header,tablefmt="orgtbl")
    text=re.sub('\|\|[\-\+]+\|\|\n', '', text.replace("|", "||"))
#     print(text)

    html=markdown2.markdown(text,extras=["wiki-tables"])
    html=html.replace("<tr", "<tr bgcolor='lightgrey'",1)

    
    return html
def random_uniqScriptDir(scriptspath,MSG,randomlength=8):
    a = list(string.ascii_letters)
    random.shuffle(a)
    ranUniscriptspath=(scriptspath.rstrip("/")+"/"+MSG+"/"+''.join(a[:randomlength]))
    if  os.path.exists(ranUniscriptspath):
        while True:
            random.shuffle(a)
            if  ''.join(a[:randomlength]) not in os.listdir(scriptspath): 
                ranUniscriptspath=(scriptspath.rstrip("/")+"/"+MSG+"/"+''.join(a[:randomlength]))
                break
    return ranUniscriptspath
def getrecords():
    print()
def genp():
    img=cv2.imread("lizerd.jpg")
    img=cv2.resize(img,(0,0),fx=0.5,fy=0.5)
    frame=cv2.imencode('.jpg', img)[1].tobytes()
    yield (b'--frame\r\n'b'Content-Type: image/jpeg\r\n\r\n' + frame + b'\r\n')

def gen(videoname):
    """Video streaming generator function."""
    cap = cv2.VideoCapture(videoname)
    # Read until video is completed
    while(cap.isOpened()):
      # Capture frame-by-frame
        ret, img = cap.read()
        if ret == True:
            img = cv2.resize(img, (0,0), fx=1, fy=1) 
            frame = cv2.imencode('.jpg', img)[1].tobytes()
            yield (b'--frame\r\n'b'Content-Type: image/jpeg\r\n\r\n' + frame + b'\r\n')
#             time.sleep(0.07)
        else: 
            break
#     while True:
#         frame = camera.get_frame()
#         yield (b'--frame\r\n'
#                b'Content-Type: image/jpeg\r\n\r\n' + frame + b'\r\n')
# class Camera(object):
#     def __init__(self):
#         self.frames = [open(f + '.jpg', 'rb').read() for f in ['1', '2', '3']]
# 
#     def get_frame(self):
#         return self.frames[int(time()) % 3]
def scriptproduce(datadepth,collectiondepth,scriptspath,inputdataroot,softwareconfig,parametersStr,inputList,outputList,bathchOfInPath,lenOfdirtotag=0,taglist=[],selecteddepth=0,selecteddirs=[]):#selecteddepth=0 means check collectiondepth only
    
    inputstr=(" "+" ".join(taglist)+" ") if int(lenOfdirtotag)!=0 else ""
    inputstr+=" ".join([pairIn[0]+" ${"+pairIn[1]+"}" for pairIn in inputList])

    if not os.path.exists(scriptspath[0]):
        os.makedirs(scriptspath[0])
#     exit(-1)
    if len(outputList)>1:
        outputStr= "  ".join([pairOut[0]+" ${output="+outputList[0]+"|suffix="+pairOut[1]+"}" for pairOut in outputList[1:]])#['', 'bam'] ["-o",'bam'] further of this func may need to adjust ,in the html page should like input not just one text seperate options and suffix
        osufix="".join([pairOut[1] for pairOut in outputList[1:]])
    else:
        outputStr=" "+outputList[0]# cp,mv command
        osufix=""
    if len(outputList[0])>1 and not os.path.exists(outputList[0]): os.makedirs(outputList[0])
    print(datadepth,collectiondepth,scriptspath,inputdataroot,softwareconfig)
    print(inputstr,parametersStr)
    parametersStr,Ni=re.subn(r"\$\$\$\$",inputstr,parametersStr)
    
    parametersStr= inputstr+" "+parametersStr if Ni==0 else parametersStr
    parametersStr,No = re.subn(r"\&\&\&\&",outputStr,parametersStr)
    cmdline=softwareconfig+"   "+parametersStr+" " if No!=0 else softwareconfig+"   "+parametersStr+" "+outputStr
    print("after Ni,No",parametersStr,Ni,No)
    ranUniscriptspath=random_uniqScriptDir(scriptspath[0],scriptspath[1])
    print(ranUniscriptspath)
    os.makedirs(ranUniscriptspath)    

    if int(selecteddepth)==0:
        selecteddirs=[]
    if inputdataroot.strip("")=="":
        compath=longestCommonPrefix(bathchOfInPath); compath=compath.rstrip(re.split(r""+os.sep,compath)[-1])
        for inputpath in bathchOfInPath:
            operatorwithdata=OperatorWithData_webservice(inputpath,inputList,cmdline,ranUniscriptspath,taglen=lenOfdirtotag)
            print("softwareconfig0",re.split(r'\s+',softwareconfig.strip())[0],outputList)
            operatorwithdata.cmdtemplatefilename=re.split(r'\s+|'+os.sep,re.sub("\s+","_",softwareconfig.strip()))[-1]+"Get"+osufix
#             newcmdline=operatorwithdata.process(inputpath.strip(), int(datadepth), int(collectiondepth),(int(collectiondepth),selecteddirs,int(selecteddepth)))
            creatDir=inputpath.strip().lstrip(compath);
            if creatDir.rfind("/")!=-1: updir=creatDir[creatDir.rfind("/")+1:];creatDir=creatDir.replace("/","").strip()
            operatorwithdata.taglen=0
            newcmdline=operatorwithdata.process(inputpath.strip(), 0, 0,(None,selecteddirs,int(selecteddepth)))
    else:
        operatorwithdata=OperatorWithData_webservice(inputdataroot,inputList,cmdline,ranUniscriptspath,taglen=lenOfdirtotag) 
        operatorwithdata.cmdtemplatefilename=re.split(r'\s+|'+os.sep,re.sub("\s+","_",softwareconfig.strip()))[-1]+"GET"+osufix 
        print("scriptproduce cmdline",operatorwithdata.cmdtemplatefilename)
        upTodownTravelDir(inputdataroot,operatorwithdata,int(datadepth),int(selecteddepth),collection_depth=int(collectiondepth),interceptdirs=selecteddirs,rootDirnotchange=operatorwithdata.inputdatapath,Interceptor_depth_notchange=int(selecteddepth))
    sys.stdout.flush();sys.stderr.flush()
    return ranUniscriptspath