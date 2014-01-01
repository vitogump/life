# -*- coding: UTF-8 -*-
'''
Created on 2013-11-26

@author: rui
'''
import re
from NGS.BasicUtil import *
import NGS.BasicUtil.Util
from optparse import OptionParser
parser = OptionParser()
parser.add_option("-I", "--infile", dest="infilename",
                  help="write report to FILE", metavar="FILE")
parser.add_option("-W","--winwidth",dest="winwidth",help="windowidth")
parser.add_option("-q", "--quiet",
                  action="store_false", dest="verbose", default=True,
                  help="don't print status messages to stdout")
test=open("test.txt",'w')                                                                                                                                                             
(options, args) = parser.parse_args()
if __name__ == '__main__':
    phastConsfile=open(options.infilename,"r")
    L=[]
    firstline=re.split(r'\s+',phastConsfile.readline())
    currentchrom=firstline[0]
    winStart=int(firstline[1])
    L=[(winStart,int(firstline[2]),float(firstline[3]))]
    print(L)
    for line in phastConsfile:
        linelist=re.split(r'\s+', line)
        if currentchrom==linelist[0]:
            L.append((int(linelist[1]),int(linelist[2]),float(linelist[3])))
        
    caculate_phastConsValue = Caculators.Caculate_phastConsValue()
    win = Util.Window()
    print(int(options.winwidth),winStart)
    win.forPhastConsFormat(L=L, L_End_Pos=len(L),windowWidth= int(options.winwidth), Caculator=caculate_phastConsValue, winStart=winStart)
    for e in win.winValueL:
        print(*e,sep='\t',file=test)
    
    phastConsfile.close()
    test.close()