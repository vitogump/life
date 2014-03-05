# -*- coding: UTF-8 -*-
import re, numpy, sys, pickle
'''
Created on 2013-6-30

@author: rui
'''
class VCF_Data():
    def __init__(self, vcffileName):
        super().__init__()
        self.VcfMap_AllChrom = {}
        self.VcfList_A_Chrom = []
        self.VcfIndexMap = {}
        try:
            self.VcfIndexMap = pickle.load(open(vcffileName + ".myindex", 'rb'))
        except:
            VCF_Data.indexVCF(VCFName=vcffileName, indexFileName=(vcffileName + ".myindex"))
            self.VcfIndexMap = pickle.load(open(vcffileName + ".myindex", 'rb'))
    @staticmethod
    def indexVCF(VCFName, indexFileName):
        """
        {chrom:position_in_file_of_first_SNP_of_this_chrom,chrom:position,,,,,,}
        """
        vcffile = open(VCFName, 'r')
        vcfChromIndex = {}
        line = vcffile.readline()
        
        while re.search(r'^##', line) != None:
            line = vcffile.readline()
        
        if re.search(r'^#', line) != None:
            vcfChromIndex["title"] = re.split(r'\s+', line)
        else:
            print("need title'#CHROM    POS    ID    REF    ALT    QUAL    FILTER    INFO    FORMAT'")
            exit(-1)        
        currentChrom = None
        lastPosition = vcffile.tell()
        print(line)
        line = vcffile.readline()
        print(line)
        while line:      
            linelist = re.split(r"\s+", line)
            if currentChrom != linelist[0]:
                currentChrom = linelist[0]
                vcfChromIndex[currentChrom] = lastPosition
            lastPosition = vcffile.tell()
    
            line = vcffile.readline()
        pickle.dump(vcfChromIndex, open(indexFileName, 'wb'))
        vcffile.close()
    def getVcfListByChrom(self, vcfFileName, chrom):
        """
            return a list that contain all vcf record of a chrom
        """
        self.VcfList_A_Chrom = []
        vcfFile = open(vcfFileName, 'r')
        try:
            print("getVcfListByChrom", self.VcfIndexMap[chrom], chrom)            
            vcfFile.seek(self.VcfIndexMap[chrom])
            line = vcfFile.readline()
        except KeyError:
            print(chrom + "didn't find in " + vcfFileName)
            return []
        while line and (re.split(r'\s+', line))[0] == chrom:
            linelist = re.split(r'\s+', line)
            chrom = linelist[0].strip()
            pos = int(linelist[1].strip())
            REF = linelist[3].strip()
            ALT = linelist[4].strip()
            INFO = linelist[7]
            self.VcfList_A_Chrom.append((pos, REF, ALT, INFO))
            line = vcfFile.readline()
        return self.VcfList_A_Chrom
        vcfFile.close()

    def getVcfMap(self, vcfFileName):
        """
        this func is from bio\test\posAroundGene\func.py ,and did some improvement,that is add  INFO = collist[7],and add INFO into
        read the vcffile into a map which keys are chrom,values are a list of tuple
        {chrNo:[(pos,REF,ALT,INFO),(pos,REF,ALT,INFO),,,,,],chrNo:[],,,,,,},the order of the tuples in the list,is according pos,
        you we can search a record by  binary chop search
        no matter self.VcfMap_AllChrom has or has not value,the value will be clean
        """
        vcfMap = {}
        vcfFile = open(vcfFileName, 'r')
        
        line = vcfFile.readline()
        while re.search(r'^##', line) != None:
    #        print(line)
            line = vcfFile.readline()
        if re.search(r'^#', line) != None:
            lineslist = vcfFile.readlines()
        else:
            print("need title'#CHROM    POS    ID    REF    ALT    QUAL    FILTER    INFO    FORMAT'\n" + line)
            exit(-1)
    #    print("pass")
        currentLine = 0
        totalRecs = len(lineslist)
        while currentLine != totalRecs:
    #        print(currentLine)
            collist = re.split(r'\s+', lineslist[currentLine])
            chrom = collist[0].strip()
            pos = int(collist[1].strip())
            REF = collist[3].strip()
            ALT = collist[4].strip()
            INFO = collist[7]
            if chrom in vcfMap:
                vcfMap[chrom].append((pos, REF, ALT, INFO))
            else:
                vcfMap[chrom] = [(pos, REF, ALT, INFO)]
            currentLine += 1
        vcfFile.close()
        self.VcfMap_AllChrom = vcfMap
#         for line in self.VcfMap["scaffold8"]:
#             print(line,file =open("vcfMapdata.txt",'a'))
