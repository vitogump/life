# -*- coding: UTF-8 -*-
'''
Created on 2014-4-3

@author: liurui
'''
from NGS.BasicUtil import *
from optparse import OptionParser
import pickle
import re
import sys, os
import time

parser = OptionParser()
parser.add_option("-m", "--homologousgene", dest="homologousgene",
                  help="homologous file")
parser.add_option("-p", "--proteincdspair", dest="proteincdspairfile", help="proteincdspairfile")
#parser.add_option("-o", "--outfileprename", dest="outfileprename", help="outfileprename")
parser.add_option("-c", "--configure", dest="configure")

parser.add_option("-q", "--quiet",
                  action="store_false", dest="verbose", default=True,
                  help="don't print status messages to stdout")
(options, args) = parser.parse_args()
homogenefile = open(options.homologousgene, 'r')
aa_cds_pair_file = open(options.proteincdspairfile, 'r')
configure = open(options.configure, 'r')
 
cline = configure.readline()
tempPath = re.search(r'tempdir=(.*)', cline).group(1).strip()
MuscleInputFileName = tempPath + "/muscleinseq.fa"
MuscleOutputFileName = tempPath + "/muscleinseq.afa"
pamlInputCDSFileName = tempPath + "/pamlinputfile.fa"
if __name__ == '__main__':
    aa_cds_filemap = {}
    """
        {species1:[aafafile,cdsfafile,aaindex,cdsfafile],species2:[,,,],,,,,}
    """
#load aa_cds_pair file and cdsfafile aafafile into memory
    for aa_cds_line in aa_cds_pair_file:
        speciesname = re.split(r':', aa_cds_line)[0].strip()
        aa_cds_pair = re.split(r':', aa_cds_line)[1]
        aa_cds_list = re.split(r';', aa_cds_pair)
        aafafileName = aa_cds_list[0].strip()
        cdsfafileName = aa_cds_list[1].strip()
        if aafafileName != None and cdsfafileName != None:
            aa_cds_filemap[speciesname] = [open(aafafileName, 'r'), open(cdsfafileName, 'r')]
            aaindex = {}
            cdsindex = {}
            try:
                aa_cds_filemap[speciesname].append(pickle.load(open(aafafileName + ".myindex", 'rb')))
                aa_cds_filemap[speciesname].append(pickle.load(open(cdsfafileName + ".myindex", 'rb')))
            except IOError:
                Util.generateIndexByChrom(aafafileName, aafafileName + ".myindex", "transcript:")
                Util.generateIndexByChrom(cdsfafileName, cdsfafileName + ".myindex", "transcript:")
                aa_cds_filemap[speciesname].append(pickle.load(open(aafafileName + ".myindex", 'rb')))
                aa_cds_filemap[speciesname].append(pickle.load(open(cdsfafileName + ".myindex", 'rb')))
            stat = os.system("rm " + aafafileName + ".myindex " + cdsfafileName + ".myindex")
            if stat != 0:
                print("rm " + aafafileName + ".myindex " + cdsfafileName + ".myindex" + " os.system return not 0")
                exit(-1)
            print(stat)
#make the fa file as the input of the file
    homotrscpttitle = re.split(r'~', homogenefile.readline())
    for homotrscptline in homogenefile:
        homotrscptlist = re.split(r'~', homotrscptline)
        i = 0
        muscleinfile = open(MuscleInputFileName, 'w')
        for trscpt in homotrscptlist:
            
            homotrscptlist[i] = re.search(r'transcript:(.*)', trscpt).group(1).strip()
            curspecies = homotrscpttitle[i].strip()
            print(curspecies, i, homotrscptlist)
            aa_cds_filemap[curspecies][0].seek(aa_cds_filemap[curspecies][2][homotrscptlist[i]])
            #homotrscptlist[i]==currentChromNO
            refSeqMap, currentChromNO, nextChromNO = Util.getRefSeqMap(aa_cds_filemap[curspecies][0], homotrscptlist[i], mapname="transcript:")
            if "".join(refSeqMap[homotrscptlist[i]][1:]).find("*") != "-1":
                print(">" + homotrscptlist[i], file=muscleinfile)
                print("".join(refSeqMap[homotrscptlist[i]][1:]), file=muscleinfile)
            i += 1
    muscleinfile.close()        
        
        
