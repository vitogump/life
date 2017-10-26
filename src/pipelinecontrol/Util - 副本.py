# -*- coding: utf-8 -*- 
'''
Created on 2014-11-8

@author: liurui
'''
from multiprocessing.dummy import Pool
import os, re
import time

import src.web.DBA as DBA


ISOTIMEFORMAT = '%Y-%m-%d %X'
def upTodownTravelDir(rootDir, OperatorWithData, datadepth=9999, Interceptor_depth=0,curdepth=0,collection_depth=0,interceptdirs=[],rootDirnotchange="",Interceptor_depth_notchange=0):
    """
        
    """

    if Interceptor_depth==0:
        print(rootDirnotchange,"++++++++++++  =0   +++++++++++++++++++++++++")
        if interceptdirs!=[] and re.search(r"" + rootDirnotchange + "(/.*?){" + str(Interceptor_depth_notchange-1) + "}[/]([^/]+)", rootDir)==None:
            print(rootDir,"is not the path a")
            return
        elif interceptdirs!=[] and  re.search(r"" + rootDirnotchange + "(/.*?){" + str(Interceptor_depth_notchange-1) + "}[/]([^/]+)", rootDir).group(2) not in interceptdirs:
            print(rootDir,interceptdirs,"is not the path b")
            return
        print("rootDir",rootDir,"datadepth",datadepth,"Interceptor_depth", Interceptor_depth,"curdepth", curdepth,"collection_depth",collection_depth)    
        if collection_depth == 0:
            print("OperatorWithData")
            newcmdline=OperatorWithData.process(rootDir, datadepth, curdepth)
        elif collection_depth > 0:
            for elem in os.listdir(path=rootDir):
                path = rootDir + "/" + elem
                if (not os.path.isdir(path)):
                    print(path,"is not the file")
                else:
                    upTodownTravelDir(path, OperatorWithData, datadepth, Interceptor_depth, curdepth+1,collection_depth-1,interceptdirs,rootDirnotchange,Interceptor_depth_notchange)
    elif Interceptor_depth>0:
        print("++++++++++  >0   ++++++++++")
        print("rootDir",rootDir,"datadepth",datadepth,"Interceptor_depth", Interceptor_depth,"curdepth", curdepth,"collection_depth",collection_depth)        
        for elem in os.listdir(path=rootDir):
            path = rootDir + "/" + elem
            if not os.path.isdir(path):
                pass
            else:
                upTodownTravelDir(path, OperatorWithData, datadepth, Interceptor_depth-1, curdepth+1,collection_depth-1,interceptdirs,rootDirnotchange,Interceptor_depth_notchange)

class OperatorWithData():
    def __init__(self, scriptsstoredir="F:/work/pipelinecontrol/scripts"):
        self.scriptsstoredir = scriptsstoredir + "/"
    def process(self, p, d):
        print(p, d)
# myprint=OperatorWithData()
class OperatorWithData_loadintodatabase(OperatorWithData):
    def __init__(self,inputdatapath,ancestralalleletabletools,interceptdirs,vcfsuffix,drop):
        self.inputdatapath=inputdatapath
        self.ancestralalleletabletools=ancestralalleletabletools
        self.interceptdirs=interceptdirs
        self.vcfsuffix=vcfsuffix.strip()
        self.drop=drop
    def process(self,curpath,datadepth,curdepth):
        if self.interceptdirs!=[] and re.search(r".*/([^/]+)$",curpath).group(1).strip() not in self.interceptdirs:
            return        
        lists =os.walk(curpath)
        for rootStr,dirs,files in lists:
            if len(re.split(r"/",rootStr))==len(re.split(r"/",self.inputdatapath))+datadepth:
                for datafilename in files:
                    if re.search(r".*?"+self.vcfsuffix+"$", datafilename) != None:
                        tablename=self.ancestralalleletabletools.createtable(rootStr + "/" +datafilename,drop=self.drop)
                        self.ancestralalleletabletools.filldata(rootStr + "/" +datafilename,tablename=tablename)
        return "OperatorWithData_loadintodatabase return"
class OperatorWithData_mode1(OperatorWithData):
    def __init__(self, cmdtemplatefile, scriptsstoredir,taglen=1):
        super().__init__(scriptsstoredir)
        self.taglen=int(taglen)
        self.cmdtemplatefilename=re.search(r"[^/]*$",cmdtemplatefile).group(0)
        scriptcontent=open(cmdtemplatefile,'r').read()
        
        self.scriptcontext=re.search(r"([\s\S]*(\n)*)cmdline=.*",scriptcontent).group(1)
        
        self.inputdatapath=re.search(r"(\n)*inputdatafilesrootpath=\s*(.*)",self.scriptcontext).group(2)
        self.cmdline=re.search(r"(.*(\n)*)cmdline=\s*(.*)",scriptcontent).group(3)
        print(scriptcontent,self.scriptcontext,self.inputdatapath,self.cmdline,sep="\n")
        self.outputlist=re.findall(r"\${output=\s*([^\s^\|]*)\|suffix=(.*?)}",self.cmdline)

#         self.interceptdirs=interceptdirs

    def process(self, curpath, datadepth, curdepth):
        print("mode1 process")
#         if self.interceptdirs!=[] and re.search(r".*/([^/]+)$",curpath).group(1).strip() not in self.interceptdirs:
#             return
        interceptdepth=curdepth
        newcmdline = self.cmdline
        subtargets = re.findall(r"\${.*?}", newcmdline)
        targetdatasuffix = []
        for target in subtargets[:]:
            c = re.search(r'\${(.*?)}', target).group(1)
            if re.search(r"output=.*", c) != None:
                print(target, subtargets)
                subtargets.remove(target)
                continue
            targetdatasuffix.append(c)
        tagname=""
        for i in range(self.taglen):
            tagname="/"+re.search(r".*/([^/]+)"+tagname+"$", curpath).group(1)+tagname
        tagname=re.sub(r"/","_",tagname[1:])
        updirname = re.search(r".*/([^/]+)$", curpath).group(1)
        newcmdline,no_of_tags=re.subn(r"\${tag}",tagname,newcmdline)

        pathToOutputdata_createdir = ""

        if curdepth<=datadepth:
            pathToOutputdata_createdir = re.search(r"" + self.inputdatapath + "((/.*?){" + str(interceptdepth) + "}[/])", curpath + "/").group(1) 
            #leftPathName_filenamepre = re.search(r"" + self.inputdatapath + pathToOutputdata_createdir + "(.*)", curpath + "/").group(1).replace("/", ".")
            for outputtuple in self.outputlist:
                if not os.path.exists(outputtuple[0] + pathToOutputdata_createdir):
                    os.makedirs(outputtuple[0] + pathToOutputdata_createdir)
        else:
            print(curdepth, datadepth, "OperatorWithData_mode1 error")
            exit(-1)
        for outputtuple in self.outputlist:
            outputpath=re.search(r"\${output=\s*("+outputtuple[0]+")\|suffix=("+outputtuple[1]+")}",newcmdline).group(1)
            outsuffix=re.search(r"\${output=\s*("+outputtuple[0]+")\|suffix=("+outputtuple[1]+")}",newcmdline).group(2)
            if outsuffix.strip()[-1]=="/":
                if outsuffix.strip()=="/":
                    outsuffix=""
                newcmdline = re.sub(r"\${output=\s*("+outputtuple[0]+")\|suffix=("+outputtuple[1]+")}", outputpath + pathToOutputdata_createdir + outsuffix, newcmdline)
                print(outputpath + pathToOutputdata_createdir )
                if not os.path.exists(outputpath + pathToOutputdata_createdir + outsuffix):
                    os.makedirs(outputpath + pathToOutputdata_createdir  + outsuffix)
            else:
                newcmdline = re.sub(r"\${output=\s*("+outputtuple[0]+")\|suffix=("+outputtuple[1]+")}", outputpath + pathToOutputdata_createdir + updirname + "myNtosub." + outsuffix, newcmdline)
        targetdata_count=0
        for i in range(0, len(targetdatasuffix)):
            lists =os.walk(curpath)    
            for rootStr,dirs,files in lists:
                if len(re.split(r"/",rootStr))==len(re.split(r"/",self.inputdatapath))+datadepth:# reach the depth that datafiles in it
                    print("reach the depth that datafiles in it",rootStr+"/",files,targetdatasuffix)
                    for datafilename in files:
                        if re.search(r".*?" + targetdatasuffix[i]+"$", datafilename) != None:
                            targetdata_count+=1
                            option_suffix_obj = re.search(r"([-\w\d]+[=\s]+)\${(\s*" + targetdatasuffix[i] + "\s*)}", newcmdline)  # for example "INPUT=${.bam} -i ${.sam}"
                            print("option_suffix_obj",option_suffix_obj,"make new cmdline:",newcmdline)
                            optionstr = option_suffix_obj.group(1)
                            suffixstr = option_suffix_obj.group(2)
                            newcmdline=re.sub(r"[-\w\d]+[=\s]+\${\s*" + targetdatasuffix[i] + "\s*}", optionstr  + rootStr + "/" + datafilename.strip() + " " + option_suffix_obj.group(0), newcmdline)                
        newcmdline = re.sub(r"[-\w\d]+[=\s]+\${.*?}", " ", newcmdline)                
                    # sub was acted from the first to the rear most
        print("pathToOutputdata_createdir", pathToOutputdata_createdir)
#         if len(targetdatasuffix)!=0 and targetdata_count!=len(targetdatasuffix)-no_of_tags:
#             print(targetdata_count,len(targetdatasuffix),no_of_tags)
#             print("targetdata_count!=len(targetdatasuffix)")
#             return newcmdline
        newcmdline=re.sub(r"myNtosub.",str(targetdata_count)+".",newcmdline)
        try:
            print(self.scriptcontext + newcmdline, file=open(self.scriptsstoredir + self.cmdtemplatefilename + "." + updirname + "Script.sh", "a"))
        except FileNotFoundError:
            print(self.scriptcontext + newcmdline, file=open(self.scriptsstoredir + self.cmdtemplatefilename + "." + updirname + "Script.sh", "w"))
        return newcmdline


class OperatorWithData_mode2(OperatorWithData):
    def __init__(self, cmdtemplatefile, scriptsstoredir,outfilepre):
        """
        interceptdirs=([subdir names list expected in the assigned depth])
        """
        super().__init__(scriptsstoredir)
        self.count=0
        self.cmdtemplatefilename=re.search(r"[^/]*$",cmdtemplatefile).group(0)
        scriptcontent=open(cmdtemplatefile,'r').read()
        
        self.scriptcontext=re.search(r"([\s\S]*(\n)*)cmdline=.*",scriptcontent).group(1)
        
        self.inputdatapath=re.search(r"(\n)*inputdatafilesrootpath=\s*(.*)",self.scriptcontext).group(2)
        self.cmdline=re.search(r"(.*(\n)*)cmdline=\s*(.*)",scriptcontent).group(3)
        print(scriptcontent,self.scriptcontext,self.inputdatapath,self.cmdline,sep="\n")
        self.outputlist=re.findall(r"\${output=\s*([^\s^\|]*)\|suffix=(.*?)}",self.cmdline)

        
#         self.interceptdirs=interceptdirs

        self.suffixstr=""
#         outputoptionstr = re.search(r"(-[\w\d]+[=\s]+)\${output=.*\|suffix=.*?}", self.cmdline).group(1)  # for example "OUTPUT=${output} -o ${output}"
        
        for outputtuple in self.outputlist:
            outputpath=re.search(r"\${output=\s*("+outputtuple[0]+")\|suffix=("+outputtuple[1]+")}",self.cmdline).group(1)
            outsuffix=re.search(r"\${output=\s*("+outputtuple[0]+")\|suffix=("+outputtuple[1]+")}",self.cmdline).group(2)
            if not os.path.exists(outputpath):
                os.makedirs(outputpath)
            if outsuffix=="/":
                outsuffix=""# dir
            self.outputfilenamewithoutoupfilesuffix=outputpath + "/" +outfilepre
            self.newcmdline = re.sub(r"\${output=\s*("+outputtuple[0]+")\|suffix=("+outputtuple[1]+")}", outputpath + "/" +outfilepre+"."+ outsuffix + " ", self.cmdline)
#         self.newcmdline = re.sub(r"[-\w\d]+[=\s]+\${output=.*\|suffix=.*?}", outputoptionstr + outputpath + "/" + suffix + " ", self.cmdline)
        print("OperatorWithData_mode2 __init__", self.newcmdline)
#         self.suffix = suffix
    def process(self, curpath, datadepth, curdepth):
        print("mode2 process")
        print(curpath, datadepth, curdepth)
#         if self.interceptdirs!=[] and re.search(r".*/([^/]+)$",curpath).group(1).strip() not in self.interceptdirs:
#             return
        newcmdline = self.newcmdline

        option_suffix_obj = re.search(r"(\s[-\w\d]+[=\s]+)\${(.*?)}", newcmdline)  # for example "INPUT=${.bam} -i ${.sam}"
        optionstr = option_suffix_obj.group(1)
        suffixstr = option_suffix_obj.group(2)
        print(optionstr)
        if re.search(r"(\s[-]+[\w\d]+[\s]+)", optionstr)==None and re.search(r"([\w\d]+\s*[=])", optionstr)==None:
            optionstr=""
        print("optionstr",optionstr,"suffixstr",suffixstr)
        self.suffixstr=suffixstr
        datafiles = os.listdir(path=curpath)
        print("OperatorWithData_mode2", datafiles)
        lists =os.walk(curpath)    
        for rootStr,dirs,files in lists:
            if len(re.split(r"/",rootStr))==len(re.split(r"/",self.inputdatapath))+datadepth:# reach the depth that datafiles in it
                for datafilename in files:
                    if re.search(r".*?" + suffixstr+"$", datafilename) != None:
                        self.count+=1
                        newcmdline = re.sub(r"[-\w\d]+[=\s]+\${.*?}", optionstr + " " + curpath + "/" + datafilename + " " + option_suffix_obj.group(0), newcmdline)


        self.newcmdline = newcmdline
        return newcmdline




class JobTracker():#for one dir
    def __init__(self, scriptDir, NumOfThread=8):
        self.scriptDir = scriptDir
        self.NumOfThread = int(NumOfThread)
def runashell(scriptDir,scriptname):
    session=DBA.getSession()
    session.execute("update jobsstate set state='1' where scriptname='"+scriptname+"' and foldername='"+scriptDir+"'")
    session.execute("update jobsstate set startdate='"+time.strftime(ISOTIMEFORMAT, time.localtime()) +"' where scriptname='"+scriptname+"' and foldername='"+scriptDir+"'")
    session.commit()
    scriptout=re.sub(r".sh$",".out",scriptname)
    a=os.system(scriptDir+"/"+scriptname+">>"+scriptDir+"/"+scriptout+" 2>&1")
#         logfile=open(self.scriptDir+"/"+scriptout,'r')
#         logtext=logfile.read()
#         logfile.close()
    if a!=0:
        session.execute("update jobsstate set state='-1' where scriptname='"+scriptname+"' and foldername='"+scriptDir+"'")
        session.commit()
        print("JobTracker "+scriptname+" runshell error")
        return#just exit this threads the python programma still go on
    else:
        #session.execute("update jobsstate set outputinfo='"+logtext+"' where scriptname='"+scriptname+"' and foldername='"+self.scriptDir+"'")
        session.execute("update jobsstate set state='2' where scriptname='"+scriptname+"' and foldername='"+scriptDir+"'")
        session.execute("update jobsstate set finishdate='"+time.strftime(ISOTIMEFORMAT, time.localtime()) +"' where scriptname='"+scriptname+"' and foldername='"+scriptDir+"'")
        session.commit()
    return
def callsh_updateDB(scriptDir,NumOfThread):
    pool=Pool(NumOfThread)
    scriptfiles = os.listdir(path=scriptDir)
    inputscriptfiles=[]
    for filename in scriptfiles:
        if re.search(r".*\.sh$", filename) == None:
            print("skip",filename)
            scriptfiles.remove(filename)
        inputscriptfiles.append((scriptDir,filename))
    DBA.addJobs2jobstate(scriptfiles,scriptDir)
    a = os.system("chmod +x " + scriptDir + "/*.sh")
#         if a!=0:
#             print("JobTracker chmod error")
#             exit(-1)
    pool.map(runashell,inputscriptfiles)
    pool.close()
    pool.join()
    