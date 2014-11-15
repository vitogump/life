# -*- coding: utf-8 -*- 
'''
Created on 2014-11-8

@author: liurui
'''
import os, re
import time


class OperatorWithData():
    def __init__(self, scriptsstoredir="F:/work/pipelinecontrol/scripts"):
        self.scriptsstoredir = scriptsstoredir + "/"
    def process(self, p, d):
        print(p, d)
# myprint=OperatorWithData()
class OperatorWithData_mode1(OperatorWithData):
    def __init__(self, cmdline, outputpath, suffix, inputdatapath, n_subdirs,scriptcontext,scriptsstoredir):
        super().__init__(scriptsstoredir)
        self.cmdline = cmdline
        self.outputpath = outputpath
        self.suffix = suffix
        self.inputdatapath = inputdatapath
        self.n_subdirs = n_subdirs
        self.scriptcontext=scriptcontext
    def process(self, curpath, curdepth):
        print("mode1 process")
        newcmdline = self.cmdline
        subtargets = re.findall(r"\${.*?}", newcmdline)
        targetdatasuffix = []
        for target in subtargets:
            c = re.search(r'\${(.*?)}', target).group(1)
            if re.search(r"output=.*",c) != None:
                print(target, subtargets)
                subtargets.remove(target)
                continue
            targetdatasuffix.append(c)
        
        datafiles = os.listdir(path=curpath)
        updirname = re.search(r".*/([^/]+)$", curpath).group(1)
        print("OperatorWithData_mode1", datafiles)
        pathToOutputdata_createdir = ""
        leftPathName_filenamepre=""
        if self.n_subdirs <= curdepth:
            pathToOutputdata_createdir = re.search(r"" + self.inputdatapath + "((/.*?){" + str(self.n_subdirs) + "}[/])", curpath + "/").group(1)
            leftPathName_filenamepre = re.search(r"" + self.inputdatapath + pathToOutputdata_createdir + "(.*)", curpath + "/").group(1).replace("/", ".")
            if not os.path.exists(self.outputpath + pathToOutputdata_createdir):
                os.makedirs(self.outputpath + pathToOutputdata_createdir)
        elif self.n_subdirs > curdepth:
            print(curdepth, self.n_subdirs, "OperatorWithData_mode1 error")
            exit(-1)

        newcmdline = re.sub(r"\${output=.*\|suffix=.*}", self.outputpath + pathToOutputdata_createdir + leftPathName_filenamepre + updirname + "." + self.suffix, newcmdline)
        for i in range(0, len(targetdatasuffix)):
            for datafilename in datafiles:

                if re.search(r".*?" + targetdatasuffix[i], datafilename) != None:
                    option_suffix_obj = re.search(r"([-\w\d]+[=\s]+)\${(\s*" + targetdatasuffix[i] + "\s*)}", newcmdline)  # for example "INPUT=${.bam} -i ${.sam}"
                    optionstr = option_suffix_obj.group(1)
                    suffixstr = option_suffix_obj.group(2)
                    newcmdline = re.sub(r"[-\w\d]+[=\s]+\${\s*" + targetdatasuffix[i] + "\s*}", optionstr + " " + curpath + "/" + datafilename.strip()+ " " + option_suffix_obj.group(0), newcmdline)
                    
        newcmdline=re.sub(r"[-\w\d]+[=\s]+\${.*?}"," ",newcmdline)                
                    # sub was acted from the first to the rear most
        print("pathToOutputdata_createdir",pathToOutputdata_createdir)
        try:
            print(self.scriptcontext+newcmdline, file=open(self.scriptsstoredir + pathToOutputdata_createdir.replace("/", "_")[1:] + leftPathName_filenamepre + updirname + "."+ updirname + "Script.sh", "a"))
        except FileNotFoundError:
            print(self.scriptcontext+newcmdline, file=open(self.scriptsstoredir + pathToOutputdata_createdir.replace("/", "_")[1:] + leftPathName_filenamepre + updirname + "."+ updirname + "Script.sh", "w"))
        return newcmdline
class OperatorWithData_mode2(OperatorWithData):
    def __init__(self, cmdline, outputpath, suffix,scriptsstoredir):
        super().__init__(scriptsstoredir)
        outputoptionstr = re.search(r"([-\w\d]+[=\s]+)\${output=.*\|suffix=.*}", cmdline).group(1)  # for example "OUTPUT=${output} -o ${output}"

        self.newcmdline = re.sub(r"[-\w\d]+[=\s]+\${output=.*\|suffix=.*}",outputoptionstr +outputpath+"/"+ suffix  + " ",cmdline)
        print("OperatorWithData_mode2 __init__",self.newcmdline)
#         self.outputpath = outputpath
#         self.suffix = suffix
    def process(self, curpath, curdepth):
        print("mode2 process")
        newcmdline = self.newcmdline

        option_suffix_obj = re.search(r"([-\w\d]+[=\s]+)\${(.*?)}", newcmdline)  # for example "INPUT=${.bam} -i ${.sam}"
        optionstr = option_suffix_obj.group(1)
        suffixstr = option_suffix_obj.group(2)
        print(optionstr, suffixstr)
        
        datafiles = os.listdir(path=curpath)
        print("OperatorWithData_mode2", datafiles)
        for datafilename in datafiles:
#             print("curfilename", curpath + "/" + datafilename)
            if re.search(r".*?" + suffixstr, datafilename) != None:
                newcmdline = re.sub(r"[-\w\d]+[=\s]+\${.*?}", optionstr + " " + curpath + "/" + datafilename + " " + option_suffix_obj.group(0), newcmdline)
#                 print(newcmdline)

        self.newcmdline = newcmdline
        return newcmdline


class JobTracker():
    def __init__(self,scriptDir,mode="series",logfile="/tmp/JobTrackerlife.log"):
        self.scriptDir=scriptDir
        self.mode=mode
        self.logfile=logfile

    def run(self):
        ISOTIMEFORMAT='%Y-%m-%d %X'
        scriptfiles=os.listdir(path=self.scriptDir)
        print(scriptfiles)
        a=os.system("chmod +x "+self.scriptDir+"/*.sh")
        if a!=0:
            print("JobTracker chmod error")
            exit(-1)
        if self.mode=="series":
            for scriptfile in scriptfiles:
                if re.search(r".*\.sh$",scriptfile)==None:
                    print("skip",scriptfile)
                    continue
                print(scriptfile+"  "+time.strftime(ISOTIMEFORMAT,time.localtime())+"\n\n",file=open(self.logfile,"a"))
                print(self.scriptDir+"/"+scriptfile+">>"+self.logfile+" 2>&1")
                a=os.system(self.scriptDir+"/"+scriptfile+">>"+self.logfile+" 2>&1")
                if a!=0:
                    print("JobTracker run error"+scriptfile)
        if self.mode=="parallel":
            for scriptfile in scriptfiles:
                scriptfileout=re.sub(r"\.sh$",".out",scriptfile)
                if re.search(r".*\.sh$",scriptfile)==None:
                    print(scriptfileout,scriptfile)
                    continue
                print("nohup "+self.scriptDir+"/"+scriptfile+">"+self.scriptDir+"/"+scriptfileout+ " 2>&1 &")
                a=os.system("nohup "+self.scriptDir+"/"+scriptfile+">"+self.scriptDir+"/"+scriptfileout+ " 2>&1 &")
                if a!=0:
                    print("JobTracker run error"+scriptfile)

def upTodownTravelDir(rootDir, OperatorWithData, datadepth=9999, curdepth=0, mode2_Interceptor=[]):
    """
        Interceptor=([subdir names list],depth of the names expected)
    """
    print(rootDir)
    if datadepth == 0:
        # data files are under the curdepth
        newcmdline = OperatorWithData.process(rootDir, curdepth)
        print('rootDir', rootDir, newcmdline)
        return
    # now go into a deeper dir
    curdepth = curdepth + 1
    datadepth = datadepth - 1
    for elem in os.listdir(path=rootDir):
        path = rootDir + "/" + elem
        if not os.path.isdir(path):
            # this is a data file
            pass
            # print("data",path)
        else:
            # this is a folder
            if len(mode2_Interceptor) > 1 and curdepth == mode2_Interceptor[0] and elem not in mode2_Interceptor[1:]:
                continue  # this is for mode 2 only
            print("go into folder", path)
            upTodownTravelDir(path, OperatorWithData, datadepth, curdepth)
