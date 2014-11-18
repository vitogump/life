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
class OperatorWithData():
    def __init__(self, scriptsstoredir="F:/work/pipelinecontrol/scripts"):
        self.scriptsstoredir = scriptsstoredir + "/"
    def process(self, p, d):
        print(p, d)
# myprint=OperatorWithData()
class OperatorWithData_mode1(OperatorWithData):
    def __init__(self, cmdline, outputpath, suffix, inputdatapath, scriptcontext, scriptsstoredir,interceptdirs=[]):
        super().__init__(scriptsstoredir)
        self.cmdline = cmdline
        self.outputpath = outputpath
        self.suffix = suffix
        self.inputdatapath = inputdatapath
        self.interceptdirs=interceptdirs
        self.scriptcontext = scriptcontext
    def process(self, curpath, datadepth, curdepth):
        print("mode1 process")
        if self.interceptdirs!=[] and re.search(r".*/([^/]+)$",curpath).group(1).strip() not in self.interceptdirs:
            return
        interceptdepth=curdepth
        newcmdline = self.cmdline
        subtargets = re.findall(r"\${.*?}", newcmdline)
        targetdatasuffix = []
        for target in subtargets:
            c = re.search(r'\${(.*?)}', target).group(1)
            if re.search(r"output=.*", c) != None:
                print(target, subtargets)
                subtargets.remove(target)
                continue
            targetdatasuffix.append(c)
        
        updirname = re.search(r".*/([^/]+)$", curpath).group(1)

        pathToOutputdata_createdir = ""

        if curdepth<=datadepth:
            pathToOutputdata_createdir = re.search(r"" + self.inputdatapath + "((/.*?){" + str(interceptdepth) + "}[/])", curpath + "/").group(1)
            
            #leftPathName_filenamepre = re.search(r"" + self.inputdatapath + pathToOutputdata_createdir + "(.*)", curpath + "/").group(1).replace("/", ".")
            if not os.path.exists(self.outputpath + pathToOutputdata_createdir):
                os.makedirs(self.outputpath + pathToOutputdata_createdir)
        else:
            print(curdepth, datadepth, "OperatorWithData_mode1 error")
            exit(-1)

        newcmdline = re.sub(r"\${output=.*\|suffix=.*}", self.outputpath + pathToOutputdata_createdir + updirname + "." + self.suffix, newcmdline)
        for i in range(0, len(targetdatasuffix)):
            lists =os.walk(curpath)    
            for rootStr,dirs,files in lists:
                if len(re.split(r"/",rootStr))==len(re.split(r"/",self.inputdatapath))+datadepth:# reach the depth that datafiles in it
                    print(rootStr+"/",files)
                    for datafilename in files:
                        if re.search(r".*?" + targetdatasuffix[i], datafilename) != None:
                            option_suffix_obj = re.search(r"([-\w\d]+[=\s]+)\${(\s*" + targetdatasuffix[i] + "\s*)}", newcmdline)  # for example "INPUT=${.bam} -i ${.sam}"
                            optionstr = option_suffix_obj.group(1)
                            suffixstr = option_suffix_obj.group(2)
                            newcmdline=re.sub(r"[-\w\d]+[=\s]+\${\s*" + targetdatasuffix[i] + "\s*}", optionstr + " " + curpath + "/" + datafilename.strip() + " " + option_suffix_obj.group(0), newcmdline)                
        newcmdline = re.sub(r"[-\w\d]+[=\s]+\${.*?}", " ", newcmdline)                
                    # sub was acted from the first to the rear most
        print("pathToOutputdata_createdir", pathToOutputdata_createdir)
        try:
            print(self.scriptcontext + newcmdline, file=open(self.scriptsstoredir + pathToOutputdata_createdir.replace("/", "_")[1:]  + updirname + "." + updirname + "Script.sh", "a"))
        except FileNotFoundError:
            print(self.scriptcontext + newcmdline, file=open(self.scriptsstoredir + pathToOutputdata_createdir.replace("/", "_")[1:]  + updirname + "." + updirname + "Script.sh", "w"))
        return newcmdline

def upTodownTravelDir(rootDir, OperatorWithData, datadepth=9999, Interceptor_depth=0, curdepth=0):
    """
        
    """
    print(rootDir)
    if Interceptor_depth == 0:
        # data files are under the curdepth
        newcmdline = OperatorWithData.process(rootDir, datadepth, curdepth)#curdepth== the beginning value of the Interceptor_depth
        print('rootDir', rootDir, newcmdline)
        return
    # now go into a deeper dir
    curdepth = curdepth + 1
    Interceptor_depth = Interceptor_depth - 1
    for elem in os.listdir(path=rootDir):
        path = rootDir + "/" + elem
        if not os.path.isdir(path):
            # this is a data file
            pass
            # print("data",path)
        else:
            # this is a folder
#             if len(mode2_Interceptor) > 1 and curdepth == mode2_Interceptor[0] and elem not in mode2_Interceptor[1:]:
#                 continue  # this is for mode 2 only
#             print("go into folder", path)
            upTodownTravelDir(path, OperatorWithData, datadepth, Interceptor_depth, curdepth)


class OperatorWithData_mode2(OperatorWithData):
    def __init__(self, cmdline, outputpath, suffix, scriptsstoredir, interceptdirs=[]):
        """
        interceptdirs=([subdir names list expected in the assigned depth])
        """
        super().__init__(scriptsstoredir)
        self.interceptdirs=interceptdirs
        outputoptionstr = re.search(r"([-\w\d]+[=\s]+)\${output=.*\|suffix=.*}", cmdline).group(1)  # for example "OUTPUT=${output} -o ${output}"

        self.newcmdline = re.sub(r"[-\w\d]+[=\s]+\${output=.*\|suffix=.*}", outputoptionstr + outputpath + "/" + suffix + " ", cmdline)
        print("OperatorWithData_mode2 __init__", self.newcmdline)
#         self.outputpath = outputpath
#         self.suffix = suffix
    def process(self, curpath, datadepth, curdepth):
        print("mode2 process")
        if self.interceptdirs!=[] and re.search(r".*/([^/]+)$",curpath).group(1).strip() not in self.interceptdirs:
            return
        newcmdline = self.newcmdline

        option_suffix_obj = re.search(r"([-\w\d]+[=\s]+)\${(.*?)}", newcmdline)  # for example "INPUT=${.bam} -i ${.sam}"
        optionstr = option_suffix_obj.group(1)
        suffixstr = option_suffix_obj.group(2)
        print(optionstr, suffixstr)
        
        datafiles = os.listdir(path=curpath)
        print("OperatorWithData_mode2", datafiles)
        lists =os.walk(curpath)    
        for rootStr,dirs,files in lists:
            if len(re.split(r"/",rootStr))==len(re.split(r"/",self.inputdatapath))+datadepth:# reach the depth that datafiles in it
                for datafilename in files:
                    if re.search(r".*?" + suffixstr, datafilename) != None:
                        newcmdline = re.sub(r"[-\w\d]+[=\s]+\${.*?}", optionstr + " " + curpath + "/" + datafilename + " " + option_suffix_obj.group(0), newcmdline)


        self.newcmdline = newcmdline
        return newcmdline




class JobTracker():#for one dir
    def __init__(self, scriptDir, NumOfThread=8):
        self.scriptDir = scriptDir
        self.NumOfThread = int(NumOfThread)
    def __runashell(self,scriptname):
        session=DBA.getSession()
        session.execute("update jobsstat set state='1' where scriptname='"+scriptname+"' and foldername='"+self.scriptDir+"'")
        session.execute("update jobsstat set startdate='"+time.strftime(ISOTIMEFORMAT, time.localtime()) +"' where scriptname='"+scriptname+"' and foldername='"+self.scriptDir+"'")
        scriptout=re.sub(r".sh$",".out",scriptname)
        a=os.system(self.scriptDir+"/"+scriptname+">>"+self.scriptDir+"/"+scriptout+" 2>&1")
        logfile=open(self.scriptDir+"/"+scriptout,'r')
        logtext=logfile.read()
        logfile.close()
        if a!=0:
            print("JobTracker . runshell error")
            exit(-1)
        else:
            session.execute("update jobsstat set outputinfo='"+logtext+"' where scriptname='"+scriptname+"' and foldername='"+self.scriptDir+"'")
            session.execute("update jobsstat set state='2' where scriptname='"+scriptname+"' and foldername='"+self.scriptDir+"'")
            session.execute("update jobsstat set finishdate='"+time.strftime(ISOTIMEFORMAT, time.localtime()) +"' where scriptname='"+scriptname+"' and foldername='"+self.scriptDir+"'")
    def callsh_updateDB(self):
        pool=Pool(self.NumOfThread)
        scriptfiles = os.listdir(path=self.scriptDir)
        for filename in scriptfiles:
            if re.search(r".*\.sh$", filename) == None:
                print("skip",filename)
                scriptfiles.remove(filename)
        DBA.addJobs(scriptfiles,self.scriptDir)
        a = os.system("chmod +x " + self.scriptDir + "/*.sh")
        if a!=0:
            print("JobTracker chmod error")
            exit(-1)
        pool.map(self.__runashell,scriptfiles)
        
#     def run(self):
#         scriptfiles = os.listdir(path=self.scriptDir)
#         print(scriptfiles)
#         a = os.system("chmod +x " + self.scriptDir + "/*.sh")
#         if a != 0:
#             print("JobTracker chmod error")
#             exit(-1)
#         
#         mypool=Pool(self.NumOfThread)
#         
#         if self.mode == "series":
#             for scriptfile in scriptfiles:
#                 if re.search(r".*\.sh$", scriptfile) == None:
#                     print("skip", scriptfile)
#                     continue
#                 print(scriptfile + "  " + time.strftime(ISOTIMEFORMAT, time.localtime()) + "\n\n", file=open(self.logfile, "a"))
#                 print(self.scriptDir + "/" + scriptfile + ">>" + self.logfile + " 2>&1")
#                 a = os.system(self.scriptDir + "/" + scriptfile + ">>" + self.logfile + " 2>&1")
#                 if a != 0:
#                     print("JobTracker run error" + scriptfile)
#         if self.mode == "parallel":
#             for scriptfile in scriptfiles:
#                 scriptfileout = re.sub(r"\.sh$", ".out", scriptfile)
#                 if re.search(r".*\.sh$", scriptfile) == None:
#                     print(scriptfileout, scriptfile)
#                     continue
#                 print("nohup " + self.scriptDir + "/" + scriptfile + ">" + self.scriptDir + "/" + scriptfileout + " 2>&1 &")
#                 a = os.system("nohup " + self.scriptDir + "/" + scriptfile + ">" + self.scriptDir + "/" + scriptfileout + " 2>&1 &")
#                 if a != 0:
#                     print("JobTracker run error" + scriptfile)