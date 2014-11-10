# -*- coding: utf-8 -*- 
'''
Created on 2014-11-8

@author: liurui
'''
import os, re

class OperatorWithData():
    def __init__(self, scriptsdir="F:/work/pipelinecontrol/scripts"):
        self.scriptsdir = scriptsdir + "/"
    def process(self, p, d):
        print(p, d)
# myprint=OperatorWithData()
class OperatorWithData_mode1(OperatorWithData):
    def __init__(self, cmdline, outputpath, suffix, inputdatapath, n_subdirs):
        super().__init__()
        self.cmdline = cmdline
        self.outputpath = outputpath
        self.suffix = suffix
        self.inputdatapath = inputdatapath
        self.n_subdirs = n_subdirs
    def process(self, curpath, curdepth):
        print("mode1 process")
        newcmdline = self.cmdline
        subtargets = re.findall(r"\${.*?}", newcmdline)
        targetdatasuffix = []
        for target in subtargets:
            c = re.search(r'\${(.*?)}', target).group(1)
            if "output" == c:
                print(target, subtargets)
                subtargets.remove(target)
                continue
            targetdatasuffix.append(c)
        
        datafiles = os.listdir(path=curpath)
        updirname = re.search(r".*/([^/]+)$", curpath).group(1)
        print("OperatorWithData_mode1", datafiles)
        pathToOutputdata_createdir = ""
        if self.inputdatapath != self.outputpath and self.n_subdirs <= curdepth:
            pathToOutputdata_createdir = re.search(r"" + self.inputdatapath + "((/.*?){" + str(self.n_subdirs) + "}[/])", curpath + "/").group(1)
            leftPathName_filenamepre = re.search(r"" + self.inputdatapath + pathToOutputdata_createdir + "(.*)", curpath + "/").group(1).replace("/", ".")
            if not os.path.exists(self.outputpath + pathToOutputdata_createdir):
                os.makedirs(self.outputpath + pathToOutputdata_createdir)
        elif self.n_subdirs > curdepth:
            print(curdepth, self.n_subdirs, "OperatorWithData_mode1 error")
            exit(-1)
        for i in range(0, len(targetdatasuffix)):
            for datafilename in datafiles:
                if re.search(r".*?" + targetdatasuffix[i], datafilename) != None:
                    newcmdline = re.sub(r"\${\s*" + targetdatasuffix[i] + "\s*}", " " + curpath + "/" + datafilename + " ", newcmdline)
                    if self.inputdatapath == self.outputpath:  # input data files and output data files are in the same dir. and this situation leftPathName_filenamepre==None 
                        newcmdline = re.sub(r"\${output}", self.outputpath + "/" + updirname + "." + self.suffix, newcmdline)
                    else:  # when curdepth ==  self.n_subdirs,the leftPathName_filenamepre contain the updirname
                        newcmdline = re.sub(r"\${output}", self.outputpath + pathToOutputdata_createdir + leftPathName_filenamepre + updirname + "." + self.suffix, newcmdline)
                    
                        
                    # sub was acted from the first to the rear most
        print(newcmdline, file=open(self.scriptsdir + pathToOutputdata_createdir.replace("/", ".")[1:] + updirname + "_script.sh", "a"))
        return newcmdline
class OperatorWithData_mode2(OperatorWithData):
    def __init__(self, cmdline, outputpath, suffix):
        super().__init__()
        optionstr = re.search(r"([-\w\d]+[=\s]+)\${output}", cmdline).group(1)  # for example "INPUT=${.bam} -i ${.sam}"

        self.newcmdline = re.sub(r"[-\w\d]+[=\s]+\${output}",optionstr +outputpath+"/"+ suffix  + " ",cmdline)
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
