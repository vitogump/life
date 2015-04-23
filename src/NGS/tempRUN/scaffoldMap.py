'''
Created on 2015-4-17

@author: liurui
'''

import pickle,re
from optparse import OptionParser
from NGS.BasicUtil import Util
import src.NGS.BasicUtil.DBManager as dbm


parser = OptionParser()

parser.add_option("-r","--reffa",dest="reffa",help="scaffold per line")#
parser.add_option("-R","--reffa_linktoDB",dest="reffa_linktoDB",help="correspond to db")#

(options, args) = parser.parse_args()

reffa_linktoDB_Name=options.reffa_linktoDB.strip()
reffa_linktoDB_hanlder=open(reffa_linktoDB_Name,'r')
reffa_linktoDB_idxName = reffa_linktoDB_Name + ".myindex"


reffahanlder=open(options.reffa.strip(),'r')
outfile=open(options.reffa.strip()+"MAPTO"+reffa_linktoDB_Name,'w')
if __name__ == '__main__':
    dbtools = dbm.DBTools(Util.ip, Util.username, Util.password, Util.genomeinfodbname)
    try:
        refidxByChr = pickle.load(open(reffa_linktoDB_idxName, 'rb'))
    except IOError:
        Util.generateIndexByChrom(reffa_linktoDB_Name, reffa_linktoDB_idxName)
        refidxByChr = pickle.load(open(reffa_linktoDB_idxName, 'rb'))

    for onelineAscaffold in reffahanlder:
        onelineAscaffold=onelineAscaffold.lower()
        if re.search(r'^>',onelineAscaffold)!=None:
            current_scaffold=re.search(r'^>(.*)',onelineAscaffold).group(1).strip()
        else:
            current_len=len(onelineAscaffold.strip())
            
            selectedchr=dbtools.operateDB("select","select * from "+Util.pekingduckchromtable+" where chrlength="+str(current_len))
            print(current_scaffold,"select * from "+Util.pekingduckchromtable+" where chrlength="+str(current_len),"\n",selectedchr)
            if len(selectedchr)==0:
                print(current_scaffold,"nofound",sep="\t",file=outfile)
            elif len(selectedchr)==1:
                print(current_scaffold,selectedchr[0][0].strip(),sep="\t",file=outfile)
            elif len(selectedchr)>1:
                Found_breakUpperloop=False
                for chr_rec in selectedchr:
                    if Found_breakUpperloop:
                        break
                    reffa_linktoDB_hanlder.seek(refidxByChr[chr_rec[0].strip()])
                    posidx=0
                    for line in reffa_linktoDB_hanlder:
                        if re.search(r'^[>]', line) != None:
                            print(current_scaffold,chr_rec[0].strip(),sep="\t",file=outfile)
                            Found_breakUpperloop=True
                            break
                        posidx=onelineAscaffold.find(line.strip().lower(),posidx)
                        if posidx==-1:
                            break
                        else:
                            posidx=posidx+len(line.strip())-1
                else:
                    if not Found_breakUpperloop:
                        print(current_scaffold,"nofound",sep="\t",file=outfile)
    outfile.close()
    reffahanlder.close()