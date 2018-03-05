'''
Created on 2015-10-17

@author: liurui
'''
from numpy import array, savetxt, loadtxt, std
from optparse import OptionParser
from os.path import sys
import pickle
import random
import re, os
import time,math

from NGS.BasicUtil import Util
import matplotlib.pyplot as pyplot
import pylab
import matplotlib.font_manager as fm

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
pythonpath_pre="/home/liurui/software/Python-2.7/python /home/liurui/life/src/NGS/Analysis/usedadiPy2_7/dadicode.py "
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
    print(options.model)
    ll_param_MAPlist["likelihood"]=[]
    ll_param_MAPlist["theta"]=[]
    bootstraps=[]
    for n,v,l,u in options.parameters:
        ll_param_MAPlist[n]=[]
    print(options.bootstrap[0])
    residualarraylist=[]
    residualhistlist=[]
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
#             initvalue=random.gauss(float(v),0.01)
#             while initvalue>float(u) or initvalue<float(l):
#                 initvalue=random.gauss(float(v),0.01)
#             paramslist.append(float(initvalue))
            paramslist.append(float(v))
            lower_boundlist.append(float(l))
            upper_boundlist.append(float(u))
#             ll_param_MAPlist[n].append()
        #produce command and run
            pythonpath=pythonpath+" -p "+n+" "+str(v)+" "+l+" "+u+" "
#             pythonpath=pythonpath+" -p "+n+" "+str(initvalue)+" "+l+" "+u+" "
        if randomstr!=None:
            os.system("rm "+namestr+options.tag+options.model+randomstr+".parameter")
        randomstr=Util.random_str()
        print(pythonpath+" -b "+randomstr+" "+str(int(options.bootstrap[1])))
        sys.stdout.flush()
        a=call_system(pythonpath+" -b "+randomstr+" "+str(int(options.bootstrap[1])))
        if a!=0:
            print("cycle",i,a,"wrong")
            continue
        #collection result
        print(options.fsfile+namestr+options.tag+options.model+"array.pickle")
        u=pickle._Unpickler(open(options.fsfile+namestr+options.tag+options.model+"array.pickle","rb"))
        u.encoding='latin1'
        residualarray=u.load() #pickle.load(open(options.fsfile+namestr+options.tag+options.model+"array.pickle","rb"))
        u=pickle._Unpickler(open(options.fsfile+namestr+options.tag+options.model+"hist.pickle","rb"))
        u.encoding='latin1'
        residualhis=u.load()#pickle.load(open(options.fsfile+namestr+options.tag+options.model+"hist.pickle","rb"))
        bif=open(options.fsfile+namestr+options.tag+options.model+randomstr+"btstrap.temp",'r')
        btstrap=[]
        for e in bif:
            elist=re.split(r"\s+",e.strip())
            for ee in elist:
                btstrap.append(float(ee))
        bif.close()
        os.system("rm "+options.fsfile+namestr+options.tag+options.model+randomstr+"btstrap.temp")
            
#         u=pickle._Unpickler(open(options.fsfile+namestr+options.tag+options.model+randomstr+"btstrap.pickle","rb"))
#         u.encoding='latin1'
#         btstrap=u.load()
        residualarraylist.append(residualarray)
        residualhistlist.append(residualhis)
        bootstraps.append(btstrap)
        inf=open(namestr+options.tag+options.model+randomstr+".parameter","r")
        for resultline in inf:
            linelist=re.split(r"\s+",resultline.strip())
            if len(linelist)>=2:
                name=linelist[0]
                convert_value=linelist[1]
                value=linelist[2]
                ll_param_MAPlist[name].append((convert_value,value))
        inf.close()
#         os.system("rm "+namestr+options.tag+options.model+randomstr+".parameter")
        print(ll_param_MAPlist)
        print(bootstraps)
##################
    pickle.dump(residualarraylist,open(options.fsfile+namestr+options.tag+options.model+"arraylist.pickle",'wb'))
    pickle.dump(residualhistlist,open(options.fsfile+namestr+options.tag+options.model+"histlist.pickle",'wb'))
    of=open(options.fsfile[:20]+"_"+namestr+options.model+options.tag+".final_parameter","w")
    for a in sorted(ll_param_MAPlist.keys()):
        print(a+"convertvalue"+" "+a+"value",end="\t",file=of)
    else:
        print("",file=of)
    for i in range(len(ll_param_MAPlist["likelihood"])):
        for a in sorted(ll_param_MAPlist.keys()):
            print(ll_param_MAPlist[a][i][0],ll_param_MAPlist[a][i][1],end="\t",file=of)
        else:
            print("",file=of)
    of.close()

    print(bootstraps)
    bootstraps = array(bootstraps)
    savetxt(namestr+options.tag+options.model+'2Dboots.npy', bootstraps)
    bootstraps = loadtxt(namestr+options.tag+options.model+'2Dboots.npy')    
    sigma_boot = std(bootstraps, axis=0)[1:]
    
    
#     fig=pylab.figure(1)
#     fig.clear2
    a=math.ceil((len(options.parameters)-2)/2)
    for i in range(len(options.parameters)):
        fig = pyplot.figure(36, figsize=(18,18))
        fig.clear()
        ax = fig.add_subplot(1,1,1)
        ax.hist(bootstraps[:,i+2], bins=20, normed=True)# +2 (ll,theta)
        fig.savefig('hist'+namestr+options.tag+options.model+options.parameters[i][0]+'.png', dpi=300)
      

# ax.hist(bootstraps[:,1], bins=20, normed=True)
# lims = ax.axis()
def collectandcheckboundary(maplist,final_parameterF,totalmaplist=None,bootstraps_list=None):
    if totalmaplist==None :
        totalmaplist={"s_1":[]};bootstraps_list=[]
        for k in maplist.keys():
            totalmaplist[k]=[]
    else:
        bootstraps_list=[list(r) for r in bootstraps_list]
    j=0
    for i in range(len(maplist["likelihood"])):
        btstrap=[]
        if float(maplist["Ts"][i][0])> 1000000 or  float(maplist["s"][i][1])>0.9999 or float(maplist["s"][i][1])<0.1 or float(maplist["m12"][i][1])>19 or float(maplist["m21"][i][1])>19:# or float(maplist["nuW"][i][1])>19 or float(maplist["nuD"][i][1])>19 or float(maplist["nuW"][i][1])<9e-05 or float(maplist["nuD"][i][1])<9e-05 :#or
            print(i,maplist["Ts"][i],maplist["s"][i])#,maplist["m12"][i])
            continue
        Nref=float(maplist["theta"][i][1])/(4*9.97e-10*468232.1666268)
        for a in reversed(sorted(maplist.keys())):
            if a=='s':
                totalmaplist["s_1"]+=[((1-float(maplist[a][i][1]))*Nref,1-float(maplist[a][i][1]))]
                btstrap.append(float(totalmaplist["s_1"][j][0]));btstrap.append(float(totalmaplist["s_1"][j][1]))
                totalmaplist[a]+=[((float(maplist[a][i][1]))*Nref,float(maplist[a][i][1]))] 
            else:
                totalmaplist[a]+=[(float(maplist[a][i][0]),float(maplist[a][i][1]))]
            btstrap.append(float(totalmaplist[a][j][0]));btstrap.append(float(totalmaplist[a][j][1]))
        bootstraps_list.append(btstrap);j+=1
    f=open(final_parameterF,'w')
    for k in sorted(totalmaplist.keys()):
        print(k+"convert",k,sep="\t",end="\t",file=f)
    print("",file=f)
    for i in range(len(totalmaplist["likelihood"])):
        for k in sorted(totalmaplist.keys()):
            print(*totalmaplist[k][i],sep="\t",end="\t",file=f)
        print("",file=f)
    f.close()
    bootstraps_list = numpy.array(bootstraps_list)
    numpy.savetxt(final_parameterF[:-16]+'2Dboots.npy', bootstraps_list)
    i=0
    for k in reversed(sorted(totalmaplist.keys())):
        print(k)
        fig = pyplot.figure(36, figsize=(36,26))
        fig.clear()
        ax = fig.add_subplot(1,1,1)
        pyplot.xticks(fontsize=36)
        pyplot.yticks(fontsize=36)
        ax.hist(bootstraps_list[:,i*2], bins=20, normed=True,facecolor='g',alpha=1)
        ax.set_title(final_parameterF[-21:-16]+k, fontsize='x-large',size=36)
        fig.savefig(final_parameterF[:-16]+k+'hist.png', dpi=360)
        i+=1
    return totalmaplist,bootstraps_list