'''
Created on 2017年11月23日

@author: liurui
'''

import re,sys
if __name__ == '__main__':
    samplename=open(sys.argv[1],'r')
    filelistexistmark={"71S":[],"113S":[],"3S":[]}
    filenamehandler1=open(sys.argv[2],'r')
    filenamehandler2=open(sys.argv[3],'r')
    filenamehandler3=open(sys.argv[4],'r')
    fw=open("doessamplehasItsfile.txt",'w')
    filenamelist1=filenamehandler1.readlines()
    filenamelist2=filenamehandler2.readlines()
    filenamelist3=filenamehandler3.readlines()
    for accId in samplename:
        count=0
        count2=0
        print("Accession ID1:",accId.strip(),end="\t",file=fw)
        for line in filenamelist1:
            if re.search(r''+accId.strip(), line.strip())!=None:
                count+=1
                filelistexistmark["113S"].append(filenamelist1.index(line))
                if "2.fq.gz" in line.strip() and "adapter" not in line.strip() :
                    count2+=1
                print(line.strip(),end="\t",file=fw)
            else:
                pass
        else:
            
            if count==0:
                pass
                #print("sample does not exist a file",end="\t",file=fw)
            else:
                print("\nexist in folder 113S.files:",count,"sampledata:",count2,file=fw)
        if count!=0:
            continue
        for line in filenamelist2:
            if re.search(r''+accId.strip(), line.strip())!=None:
                count+=1
                filelistexistmark["3S"].append(filenamelist2.index(line))
                if "2.fq.gz" in line.strip() and "adapter" not in line.strip() :
                    count2+=1
                print(line.strip(),end="\t",file=fw)
            else:
                pass
        else:
            if count==0:
                pass
                #print("sample does not exist a file",end="\t",file=fw)
            else:
                print("\nexist in folder 3S.files:",count,"sampledata:",count2,file=fw)
        if count!=0:
            continue
        for line in filenamelist3:
            if re.search(r''+accId.strip(), line.strip())!=None:
                count+=1
                filelistexistmark["71S"].append(filenamelist3.index(line))
                if "2.fq.gz" in line.strip() and "adapter" not in line.strip() :
                    count2+=1
                print(line.strip(),end="\t",file=fw)
            else:
                pass
        else:
            if count==0:
                print("\n no sample data exist",file=fw)
                #print("sample does not exist a file",end="\t",file=fw)
            else:
                print("\nexist in folder 71S.files:",count,"sampledata:",count2,file=fw)
    fww=open("filenotinfilelist.txt",'w');fw.close();filenamehandler1.close();filenamehandler2.close();filenamehandler3.close()
    for k in filelistexistmark.keys():
        filelistexistmark[k].sort()
        print(filelistexistmark[k])
    for idx in range(len(filenamelist3)):
        if idx not in filelistexistmark["71S"]:
            print(filenamelist3[idx].strip(),file=fww)
    for idx in range(len(filenamelist2)):
        if idx not in filelistexistmark["3S"]:
            print(filenamelist2[idx].strip(),file=fww)
    for idx in range(len(filenamelist1)):
        if idx not in filelistexistmark["113S"]:
            print(filenamelist1[idx].strip(),file=fww)
    fww.close()