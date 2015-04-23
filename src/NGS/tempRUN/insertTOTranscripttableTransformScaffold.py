'''
Created on 2015-4-18

@author: liurui
'''
from optparse import OptionParser
from NGS.BasicUtil import Util
import src.NGS.BasicUtil.DBManager as dbm
import re
parser = OptionParser()

parser.add_option("-T","--transcriptTablefile",dest="transcriptTablefile",help="scaffold per line")#
parser.add_option("-m","--refchrNoMap",dest="refchrNoMap",help="scaffold per line")#

(options, args) = parser.parse_args()
ref_chrNO_Map={}
f=open(options.refchrNoMap,'r')
for line in f:
    linelist=re.split(r"\s+",line)
    ref_chrNO_Map[linelist[0].strip()]=linelist[1].strip()
f.close()
if __name__ == '__main__':
    f=open(options.transcriptTablefile,'r')
    genomedbtools = dbm.DBTools(Util.ip, Util.username, Util.password, Util.genomeinfodbname)
    for line in f:
        linelist=re.split(r"\s+",line.strip())
        if linelist[3].strip() in ref_chrNO_Map and ref_chrNO_Map[linelist[3].strip()]!="nofound":
            transcript_ID=linelist[0].strip()
            geneID=linelist[1].strip()
            geneName=linelist[2].strip()
            chrID=ref_chrNO_Map[linelist[3].strip()]
            GCcontent=float(linelist[4].strip())
            trscpt_start_pos=int(linelist[5])
            trscpt_end_pos=int(linelist[6])
            strand=linelist[7]
            genomedbtools.operateDB("insert","insert ignore into "+Util.TranscriptGenetable+"(transcript_ID,geneID,geneName,chrID,GCcontent,trscpt_start_pos,trscpt_end_pos,strand) value(%s,%s,%s,%s,%s,%s,%s,%s)",data=(transcript_ID,geneID,geneName,chrID,GCcontent,trscpt_start_pos,trscpt_end_pos,strand))
    f.close()
        
        
        