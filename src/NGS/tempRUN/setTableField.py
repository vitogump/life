'''
Created on 2015-4-18

@author: liurui
'''
import re
from optparse import OptionParser
from NGS.BasicUtil import Util
import src.NGS.BasicUtil.DBManager as dbm


parser = OptionParser()

parser.add_option("-T","--transcriptTablefile",dest="transcriptTablefile",help="scaffold per line")#

(options, args) = parser.parse_args()
f=open(options.transcriptTablefile,'r',encoding='utf-8')
titleline=f.readline()
titlelist=re.split(r",",titleline.strip())
FieldidxINFile=titlelist.index("Duck gene")

FieldINTable="geneID"
NameINTable="geneName"
NameINFile="Gene Symbol"
NameidxINFile=titlelist.index(NameINFile)
if __name__ == '__main__':
    dbtools = dbm.DBTools(Util.ip, Util.username, Util.password, Util.genomeinfodbname)
    for line in f:
        linelist=re.split(r",",line)
#         print("select * from "+Util.TranscriptGenetable+" where "+FieldINTable+"='"+linelist[FieldidxINFile]+"'")
        rec=dbtools.operateDB("select","select * from "+Util.TranscriptGenetable+" where "+FieldINTable+"='"+linelist[FieldidxINFile]+"'")
        if len(rec)==1 and rec[0][2]=="":
            print("update "+Util.TranscriptGenetable+" set "+NameINTable+"='"+linelist[NameidxINFile]+"' where "+FieldINTable+"='"+linelist[FieldidxINFile]+"'")
            dbtools.operateDB("update","update "+Util.TranscriptGenetable+" set "+NameINTable+"='"+linelist[NameidxINFile]+"' where "+FieldINTable+"='"+linelist[FieldidxINFile]+"'")
    f.close()
        