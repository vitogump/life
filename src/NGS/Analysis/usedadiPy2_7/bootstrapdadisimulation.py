'''
Created on 2015-10-17

@author: liurui
'''
from optparse import OptionParser
from os.path import sys
import random
import re, os
import time

from numpy import array

from NGS.BasicUtil import Util


randomstr=None
parser = OptionParser()
parser.add_option("-n", "--popname", dest="popname",action="append",nargs=2,help="fs_file_name projection")
parser.add_option("-f", "--fsfile", dest="fsfile",help="fs file name ")
parser.add_option("-l", "--genomelengthwhichsnpfrom", dest="genomelengthwhichsnpfrom",help="fs file name ")
parser.add_option("-b", "--bootstrap", dest="bootstrap",default="100 10",nargs=2,help="times time ")
parser.add_option("-m", "--model", dest="model",help="1,model1 2,model2 ....")
parser.add_option("-p","--parameters",dest="parameters",action="append",nargs=4,help="""parametername lower upper""")
parser.add_option("-t", "--threads",dest="threads", default="4",help="don't print status messages to stdout")                                                                                                                
parser.add_option("-T", "--tag",dest="tag", default="TAG",help="don't print status messages to stdout")

(options, args) = parser.parse_args()
ll_param_MAPlist={}
ISOTIMEFORMAT='%Y-%m-%d %X'
pythonpath_pre="/home/liurui/softwarepakage/Python-2.7/python /home/liurui/life/src/NGS/Analysis/usedadiPy2_7/dadicode.py "
def call_system(commandline):
    a=os.system(commandline)
    return a
if __name__ == '__main__':
    popnamelist=[]
    projectionlist=[]
    for popname,projection in options.popname:
        popnamelist.append(popname)
        projectionlist.append(int(projection))
        pythonpath_pre=pythonpath_pre+" -n "+popname+" "+ projection
        
    namestr=""
    for name in popnamelist:
        namestr+=name        
    pythonpath_pre=pythonpath_pre+" -f "+options.fsfile+" -l "+options.genomelengthwhichsnpfrom+" -m "+options.model+" -T "+options.tag
    ll_param_MAPlist["likelihood"]=[]
    ll_param_MAPlist["theta"]=[]
    for n,v,l,u in options.parameters:
        ll_param_MAPlist[n]=[]
    print(options.bootstrap[0])
    for i in range(int(options.bootstrap[0])):
        print(time.strftime( ISOTIMEFORMAT, time.localtime() ))
        print("cycly ",i)
        
        paramslist=[]
        upper_boundlist=[]
        lower_boundlist=[]
        paramsname=[]
        pythonpath=pythonpath_pre
        for n,v,l,u in options.parameters:
            paramsname.append(n)
            #random initial value
            initvalue=random.gauss(float(v),0.01)
            while initvalue>float(u) or initvalue<float(l):
                initvalue=random.gauss(float(v),0.01)
            paramslist.append(float(initvalue))
            lower_boundlist.append(float(l))
            upper_boundlist.append(float(u))
#             ll_param_MAPlist[n].append()
        #produce command and run
            pythonpath=pythonpath+" -p "+n+" "+str(initvalue)+" "+l+" "+u+" "
        if randomstr!=None:
            os.system("rm "+namestr+options.tag+randomstr+".parameter")
        randomstr=Util.random_str()
        print(pythonpath+" -b "+randomstr+" "+str(70*int(options.bootstrap[1])))
        sys.stdout.flush()
        a=call_system(pythonpath+" -b "+randomstr+" "+str(70*int(options.bootstrap[1])))
        if a!=0:
            print("cycle",i,a,"wrong")
            continue
        #collection result
        inf=open(namestr+options.tag+randomstr+".parameter","r")
        for resultline in inf:
            linelist=re.split(r"\s+",resultline.strip())
            if len(linelist)>=2:
                name=linelist[0]
                convert_value=linelist[1]
                value=linelist[2]
                ll_param_MAPlist[name].append((convert_value,value))
        inf.close()
        print(ll_param_MAPlist)
##################
    of=open(options.fsfile[:20]+"_"+namestr+options.tag+".final_parameter","w")
    for a in sorted(ll_param_MAPlist.keys()):
        print(a,end="\t",file=of)
    else:
        print("",file=of)
    for i in range(len(ll_param_MAPlist["likelihood"])):
        for a in sorted(ll_param_MAPlist.keys()):
            print(ll_param_MAPlist[a][i][0],ll_param_MAPlist[a][i][1],end="\t",file=of)
        else:
            print("",file=of)
    of.close()
    
    