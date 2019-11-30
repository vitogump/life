# -*- coding: UTF-8 -*-
'''
Created on 2015-5-7

@author: liurui
this is for dadi
'''
import copy, random, re, os
from optparse import OptionParser
import sys
from src.NGS.BasicUtil import Util
import src.NGS.BasicUtil.DBManager as dbm


parser = OptionParser()
parser.add_option("-t","--toplevelsnptable",dest="toplevelsnptable",default="ducksnp_toplevel",help="depth of the folder to output")
parser.add_option("-m","--minlength",dest="minlength",help="require least chrom length")
parser.add_option("-q","--quantizing",dest="quantizing",action="append",help="select from some first param of -v ")
parser.add_option("-n","--noofindvds2quantizing",dest="noofindvds2quantizing",help="together with -q ")
parser.add_option("-C","--countsnpnumberfromvcf",dest="countsnpnumberfromvcf",default=None)
parser.add_option("-d","--snpperkb",dest="snpperkb")
parser.add_option("-j","--joinmanner",dest="joinmanner")
parser.add_option("-o","--outputfilename",dest="outputfilename")
parser.add_option("-v", "--vcftablelist", dest="vcftablelist",action="append",default=[],nargs=2,help="")
(options, args) = parser.parse_args()
joinmanner=" "+options.joinmanner+" join "
minlength=options.minlength;toplevelsnptable=options.toplevelsnptable;snpperkb=float(options.snpperkb);vcftableslist=options.vcftablelist#minAN=options.minAN
dadisnpfile=open(options.outputfilename+"dilutetodensity"+options.snpperkb.strip(),'w')
outgroupidx_in_topleveltable=[8,12];minoutgroupdepth=30#12 taihu 8 fanya
randomnamefile_recordchr=Util.random_str(8)
randomnamef_recchr=open(randomnamefile_recordchr,"w")
noofindvds2quantizing=int(options.noofindvds2quantizing)
allsnps=0
if __name__ == '__main__':
    print(sys.argv)
    genomedbtools = dbm.DBTools(Util.ip, Util.username, Util.password, Util.genomeinfodbname)
    dbvariantstools=dbm.DBTools(Util.ip, Util.username,Util.password, Util.vcfdbname)
    toplevelsnptable_titlelist=[a[0].strip() for a in dbvariantstools.operateDB("select", "select column_name  from information_schema.columns where table_schema='" + "ninglabvariantdata" + "' and table_name='" + toplevelsnptable + "'")]
    selectedchroms=genomedbtools.operateDB("select","select * from "+Util.pekingduckchromtable+" where chrlength>="+minlength)
    ######################## title print ##############################
#     print(Util.pekingduckchromtable[:9],toplevelsnptable_titlelist[outgroupidx_in_topleveltable[0]][:6]+toplevelsnptable_titlelist[outgroupidx_in_topleveltable[1]][:6],"Allele1",sep="\t",end="\t",file=dadisnpfile)
    print(Util.pekingduckchromtable[:9],toplevelsnptable_titlelist[outgroupidx_in_topleveltable[0]][:8],"Allele1",sep="\t",end="\t",file=dadisnpfile)
    vcftablesidxlist_toquantizing=[]
    idx=0
    for vcftable_name,minAN in vcftableslist:
        if options.quantizing!=None:
            for quantizingvcftable in options.quantizing:
                if quantizingvcftable.strip()==vcftable_name.strip():
                    vcftablesidxlist_toquantizing.append(idx)
        idx+=1
        popName=re.split(r'_',vcftable_name)[0]
        print(popName,end="\t",file=dadisnpfile)
    print("vcftablesidxlist_toquantizing",vcftablesidxlist_toquantizing)
    if options.quantizing!=None:
        print("quantizpool",end="\t",file=dadisnpfile)
    print("Allele2",end="\t",file=dadisnpfile)
    for vcftable_name,minAN in vcftableslist:
        popName=re.split(r'_',vcftable_name)[0]
        print(popName,end="\t",file=dadisnpfile)
    if options.quantizing!=None:
        print("quantizpool",end="\t",file=dadisnpfile)
    print("Gene\tPosition",file=dadisnpfile)
    ############               finish title print ##################################
    #######################    start make sqlstatment    ##################################
    totallength=0
    totallengthduilt=0
    allsnpcount=0
    totalduiltsnp=0
    for row in selectedchroms:
        duiltedsnpcountforcurrentchr=0
        currentchrID=row[0]
        currentchrLen=int(row[1])
        totallength+=currentchrLen
        sqlselectstatementpart="select t.*"
        sqlselectstatementpart_count="select count(*) "
#         sqlfromstatementpart=" from "+toplevelsnptable+" as t "
        for vcftable,minAN in vcftableslist:
            if re.search(r"indvd[^/]+",vcftable)!=None:
                sqlselectstatementpart=sqlselectstatementpart+","+vcftable.strip()+".AC as "+vcftable+"_AC,"+vcftable.strip()+".AN as "+vcftable+"_AN"
#                 sqlselectstatementpart_count=sqlselectstatementpart_count+","+vcftable.strip()+".AC as "+vcftable+"_AC,"+vcftable.strip()+".AN as "+vcftable+"_AN"
            elif re.search(r"pool[^/]+",vcftable)!=None:
                sqlselectstatementpart=sqlselectstatementpart+",floor("+vcftable.strip()+".AF*"+minAN+") as "+vcftable+"_AC,"+minAN+" as "+vcftable+"_AN"
#                 sqlselectstatementpart_count=sqlselectstatementpart_count+",floor("+vcftable.strip()+".AF*"+minAN+") as "+vcftable+"_AC,"+minAN+" as "+vcftable+"_AN"
        sqlselectstatementpart+=" from "+toplevelsnptable+" as t "
        sqlselectstatementpart_count+= " from "+vcftableslist[0][0]
        for vcftable,minAN in vcftableslist:
            sqlselectstatementpart=sqlselectstatementpart+joinmanner+vcftable.strip()+" using(chrID,snp_pos)"
#             sqlselectstatementpart_count=sqlselectstatementpart_count+" outer join "+vcftable.strip()+" using(chrID,snp_pos)"
        sqlselectstatementpart_count_left=sqlselectstatementpart_count
        sqlselectstatementpart_count_right=sqlselectstatementpart_count
        for vcftable,minAN in vcftableslist[1:]:
            sqlselectstatementpart_count_left=sqlselectstatementpart_count+" left join "+vcftable.strip()+" using(chrID,snp_pos)"
            sqlselectstatementpart_count_right=sqlselectstatementpart_count+" right join "+vcftable.strip()+" using(chrID,snp_pos)"
            sqlselectstatementpart_count=sqlselectstatementpart_count_left+" union "+sqlselectstatementpart_count_right
#         print("either fanya or weigeon are fixed")
#         sqlselectstatementpart=sqlselectstatementpart+" where chrID='"+currentchrID+"' and t.context is not null and length(t.ref_base)=length(t.alt_base)  and (t."+toplevelsnptable_titlelist[outgroupidx_in_topleveltable[0]+1] +" regexp '^0,[1234567890]+' or t."+toplevelsnptable_titlelist[outgroupidx_in_topleveltable[0]+1] +" regexp '[1234567890]+,0$' or t."+toplevelsnptable_titlelist[outgroupidx_in_topleveltable[1]+1]+" regexp '^0,[1234567890]+' or t."+toplevelsnptable_titlelist[outgroupidx_in_topleveltable[1]+1]+" regexp '[1234567890]+,0$')"
        sqlselectstatementpart=sqlselectstatementpart+" where chrID='"+currentchrID+"' and t.context is not null and length(t.ref_base)=length(t.alt_base)  and (t."+toplevelsnptable_titlelist[outgroupidx_in_topleveltable[0]+1] +" regexp '^0,[1234567890]+' or t."+toplevelsnptable_titlelist[outgroupidx_in_topleveltable[0]+1] +" regexp '[1234567890]+,0$')"# or t."+toplevelsnptable_titlelist[outgroupidx_in_topleveltable[1]+1]+" regexp '^0,[1234567890]+' or t."+toplevelsnptable_titlelist[outgroupidx_in_topleveltable[1]+1]+" regexp '[1234567890]+,0$')"
        sqlselectstatementpart_count_left=sqlselectstatementpart_count_left+" where chrID='"+currentchrID+"' and (length("+vcftableslist[0][0].strip()+".ref_base)=1 and length("+vcftableslist[0][0].strip()+".alt_base)=1) or ("+vcftableslist[0][0].strip()+".alt_base regexp '[ATCG],[ATCG]' and length("+vcftableslist[0][0].strip()+".ref_base)=1) "
        sqlselectstatementpart_count_right=sqlselectstatementpart_count_right+" where chrID='"+currentchrID+"' and (length("+vcftableslist[0][0].strip()+".ref_base)=1 and length("+vcftableslist[0][0].strip()+".alt_base)=1) or ("+vcftableslist[0][0].strip()+".alt_base regexp '[ATCG],[ATCG]' and length("+vcftableslist[0][0].strip()+".ref_base)=1) "
        for vcftable,minAN in vcftableslist:
            if re.search(r"indvd[^/]+",vcftable)!=None:
                sqlselectstatementpart=sqlselectstatementpart+" and ("+vcftable+".AN >="+minAN+" or "+vcftable+".AN is NULL) "
            elif re.search(r"pool[^/]+",vcftable)!=None:
                sqlselectstatementpart=sqlselectstatementpart+" and ("+vcftable+".DP >="+str(int(minAN)*1.5) +" or "+vcftable+".DP is NULL) "
        sqlselectstatementpart=sqlselectstatementpart+" and ( 0 "

        #variants exist in vcftablesidxlist_toquantizing
        for idx in vcftablesidxlist_toquantizing:
            if re.search(r"indvd[^/]+",vcftableslist[idx][0])!=None:
                sqlselectstatementpart=sqlselectstatementpart+" or ("+vcftableslist[idx][0]+".AF is not NULL and "+vcftableslist[idx][0]+".AF <1 )"

            elif re.search(r"pool[^/]+",vcftableslist[idx][0])!=None:
                sqlselectstatementpart=sqlselectstatementpart+" or ("+vcftableslist[idx][0]+".AF is not NULL and "+vcftableslist[idx][0]+".AF <1 )"

        sqlselectstatementpart+=") and (0 "
        #variants exist in non vcftablesidxlist_toquantizing
        idx=0
        for vcftable,minAN in vcftableslist:
            if idx in vcftablesidxlist_toquantizing:
                idx+=1
                continue
            idx+=1
            if re.search(r"indvd[^/]+",vcftable)!=None:
                sqlselectstatementpart=sqlselectstatementpart+" or ("+vcftable+".AF is not NULL and "+vcftable+".AF <1 )"
            elif re.search(r"pool[^/]+",vcftable)!=None:
                sqlselectstatementpart=sqlselectstatementpart+" or ("+vcftable+".AF is not NULL and "+vcftable+".AF <1 )"
        sqlselectstatementpart+=")"
#         sqlselectstatementpart_count+=")"
# make sql statement end start sql excute
        allsnpOfJoinTableinAchr=dbvariantstools.operateDB("select",sqlselectstatementpart)
        print(sqlselectstatementpart_count_left+" union "+sqlselectstatementpart_count_right)
        if options.countsnpnumberfromvcf==None:
            total_no_of_snpinAchrlist=dbvariantstools.operateDB("select",sqlselectstatementpart_count_left+" union "+sqlselectstatementpart_count_right)
            print(total_no_of_snpinAchrlist)
            total_no_of_snpinAchr=total_no_of_snpinAchrlist[0][0]+total_no_of_snpinAchrlist[1][0]
            allsnps=+total_no_of_snpinAchr
        else:
            total_no_of_snpinAchr=len(allsnpOfJoinTableinAchr)
        print(total_no_of_snpinAchr)
        if len(allsnpOfJoinTableinAchr)>int(currentchrLen*snpperkb)/1000:
            
            dilute = snpperkb*(currentchrLen) / (1000 * total_no_of_snpinAchr)
            totallengthduilt=totallengthduilt+currentchrLen*dilute
            allsnpOfJoinTableinAchr_sampled_idxlist=random.sample([j for j in range(len(allsnpOfJoinTableinAchr))],int(dilute*len(allsnpOfJoinTableinAchr))+1)
            allsnpOfJoinTableinAchr_sampled_idxlist.sort()
        elif len(allsnpOfJoinTableinAchr)>0:
            print("all snp in ",currentchrID)
            totallengthduilt+=currentchrLen
            allsnpOfJoinTableinAchr_sampled_idxlist=[j for j in range(len(allsnpOfJoinTableinAchr))]
        else:
            print("no snp in ",currentchrID)
            totallength-=currentchrLen
            continue
#         if options.countsnpnumberfromvcf!=None:
#             a=os.popen("grep -n "+currentchrID+" "+options.countsnpnumberfromvcf+"|wc -l")
#             snpforcurrentchr=int(a.readline().strip())
#             
#             allsnpcount+=(snpforcurrentchr-1)
#             a.close()
        print("allsnpcount",allsnpcount,"totallengthduilt",totallengthduilt,"totallength",totallength,sqlselectstatementpart)
        print(allsnpOfJoinTableinAchr_sampled_idxlist)
        for sampled_idx in allsnpOfJoinTableinAchr_sampled_idxlist:
            contextwithinspeces=allsnpOfJoinTableinAchr[sampled_idx][5].upper()
            contextoutgroup=copy.copy(contextwithinspeces)
#             firstoutgroupbase=allsnpOfJoinTableinAchr[sampled_idx][outgroupidx_in_topleveltable[0]].upper().strip();firstoutgroupdepthlist=re.split(r",",allsnpOfJoinTableinAchr[sampled_idx][outgroupidx_in_topleveltable[0]+1])
#             secondoutgroupbase=allsnpOfJoinTableinAchr[sampled_idx][outgroupidx_in_topleveltable[1]].upper().strip();secondoutgroupdepthlist=re.split(r",",allsnpOfJoinTableinAchr[sampled_idx][outgroupidx_in_topleveltable[1]+1])
# #             outgroupBase=re.split(r':',outgroup);outgroupdepth=re.split(r",",outgroupBase[1]);outgroupBase=outgroupBase[0].upper()
#             ALT=allsnpOfJoinTableinAchr[sampled_idx][4].upper().strip()
            sampled_idx_find_satisfied=sampled_idx
            #find the nearest correct rec
            direction=1
            continuesearch=-1#-1 continuesearch; 1 firstoutgroupbase; 2 secondgroupbase
            while continuesearch==-1:
                if sampled_idx_find_satisfied==len(allsnpOfJoinTableinAchr) or (sampled_idx!=allsnpOfJoinTableinAchr_sampled_idxlist[-1] and sampled_idx_find_satisfied==allsnpOfJoinTableinAchr_sampled_idxlist[allsnpOfJoinTableinAchr_sampled_idxlist.index(sampled_idx)+1] ) or sampled_idx_find_satisfied==-1 or (sampled_idx!=allsnpOfJoinTableinAchr_sampled_idxlist[0] and sampled_idx_find_satisfied==allsnpOfJoinTableinAchr_sampled_idxlist[allsnpOfJoinTableinAchr_sampled_idxlist.index(sampled_idx)-1]):
                    if direction==-1:
                        print("search snp out of range ",direction,sampled_idx,len(allsnpOfJoinTableinAchr_sampled_idxlist),sampled_idx_find_satisfied)
                        break
                    direction=-1
                    print("direction changed",direction)
                    sampled_idx_find_satisfied=sampled_idx#start again, but the opposit deriction
                else:
                    firstoutgroupbase=allsnpOfJoinTableinAchr[sampled_idx_find_satisfied][outgroupidx_in_topleveltable[1]].upper().strip();firstoutgroupdepthlist=re.split(r",",allsnpOfJoinTableinAchr[sampled_idx_find_satisfied][outgroupidx_in_topleveltable[1]+1])
                    secondoutgroupbase=allsnpOfJoinTableinAchr[sampled_idx_find_satisfied][outgroupidx_in_topleveltable[0]].upper().strip();secondoutgroupdepthlist=re.split(r",",allsnpOfJoinTableinAchr[sampled_idx_find_satisfied][outgroupidx_in_topleveltable[0]+1])
                    ALT=allsnpOfJoinTableinAchr[sampled_idx_find_satisfied][4].upper().strip()                    
#                 wigeondepthlist1=re.split(r",",snp[0][13])
#                 fanyadepthlist=re.split(r",",snp[0][9])
#                first is taihu goose;second is fanya
                if len(secondoutgroupdepthlist)==2 and (int(firstoutgroupdepthlist[0]) + int(firstoutgroupdepthlist[1])>=minoutgroupdepth or int(secondoutgroupdepthlist[0]) + int(secondoutgroupdepthlist[1])>=minoutgroupdepth) and ((firstoutgroupdepthlist[0].strip()=="0" and secondoutgroupdepthlist[0].strip()=="0") or (firstoutgroupdepthlist[1].strip()=="0" and secondoutgroupdepthlist[1].strip()=="0") ):
                    if firstoutgroupdepthlist[0].strip()=="0" and secondoutgroupdepthlist[0].strip()=="0":
                        A_base_idx=1;continuesearch=2
                    elif firstoutgroupdepthlist[1].strip()=="0" and secondoutgroupdepthlist[1].strip()=="0":
                        A_base_idx=0;continuesearch=2
                    else:
                        print("never get here")
                elif (int(secondoutgroupdepthlist[0]) + int(secondoutgroupdepthlist[1])>=minoutgroupdepth and (secondoutgroupdepthlist[0].strip()=="0" or secondoutgroupdepthlist[1].strip()=="0" )):# fanya  or (snp[0][7]=="no covered" and len(depthlist2)==2 and int(depthlist2[0]) + int(depthlist2[1])>=mindeptojudgefix and (depthlist2[1].strip()=="0" or depthlist2[0].strip()=="0")):
                    if  secondoutgroupdepthlist[0].strip()=="0":
                        A_base_idx=1;continuesearch=2
                    elif secondoutgroupdepthlist[1].strip()=="0":
                        A_base_idx=0;continuesearch=2
                else:
                    print("skip this snp,coverage not sufficent or different direction",allsnpOfJoinTableinAchr[sampled_idx_find_satisfied])
#                     print("skip snp",snp)
#                     A_base_idx=-1
                    sampled_idx_find_satisfied+=direction  
#                 if len(secondoutgroupdepthlist)==2 and int(secondoutgroupdepthlist[1])>=minoutgroupdepth and secondoutgroupdepthlist[0].strip()=="0":
#                     A_base_idx=1;continuesearch=2
#                 elif len(secondoutgroupdepthlist)==2 and int(secondoutgroupdepthlist[0])>=minoutgroupdepth and secondoutgroupdepthlist[1].strip()=="0":
#                     A_base_idx=0;continuesearch=2
# #                 if len(firstoutgroupdepthlist)==2 and len(secondoutgroupdepthlist)==2 and (int(firstoutgroupdepthlist[0])+int(firstoutgroupdepthlist[1])>=minoutgroupdepth or int(secondoutgroupdepthlist[0])+int(secondoutgroupdepthlist[1])>=minoutgroupdepth):
# #                     if (firstoutgroupdepthlist[0].strip()=="0" and secondoutgroupdepthlist[0].strip()=="0") :
# #                         A_base_idx=1;continuesearch=3
# #                     elif (firstoutgroupdepthlist[1].strip()=="0" and secondoutgroupdepthlist[1].strip()=="0"):
# #                         A_base_idx=0;continuesearch=3
# #                     else:
# #                         print("skip this snp ",allsnpOfJoinTableinAchr[sampled_idx_find_satisfied])
# #                         sampled_idx_find_satisfied+=direction
#                 else:
#                     print("skip this snp,coverage not sufficent or different direction",allsnpOfJoinTableinAchr[sampled_idx_find_satisfied])
#                     sampled_idx_find_satisfied+=direction
                #########
            else:#find
                duiltedsnpcountforcurrentchr+=1
                print("find",allsnpOfJoinTableinAchr[sampled_idx_find_satisfied])
                contextwithinspeces=allsnpOfJoinTableinAchr[sampled_idx_find_satisfied][5].upper()
                contextoutgroup=copy.copy(contextwithinspeces)
                if A_base_idx==1:
                    contextoutgroup=contextoutgroup[0]+ALT+contextoutgroup[2]
                REF=allsnpOfJoinTableinAchr[sampled_idx_find_satisfied][3].upper()
                postion=allsnpOfJoinTableinAchr[sampled_idx_find_satisfied][1]
                print(contextwithinspeces,contextoutgroup,REF,sep="\t",end="\t",file=dadisnpfile)
                k=0;l=0
                for vcftable_idx in range(len(vcftableslist)):
                    if allsnpOfJoinTableinAchr[sampled_idx_find_satisfied][len(toplevelsnptable_titlelist)+vcftable_idx*2]!=None:#AC IS NOT NULL
                        print(str(int(allsnpOfJoinTableinAchr[sampled_idx_find_satisfied][len(toplevelsnptable_titlelist)+vcftable_idx*2+1])-int(allsnpOfJoinTableinAchr[sampled_idx_find_satisfied][len(toplevelsnptable_titlelist)+vcftable_idx*2])),end="\t",file=dadisnpfile)#AN - AC
                        if options.quantizing!=None and vcftable_idx in vcftablesidxlist_toquantizing:
                            k+=int(round((int(allsnpOfJoinTableinAchr[sampled_idx_find_satisfied][len(toplevelsnptable_titlelist)+vcftable_idx*2+1])-int(allsnpOfJoinTableinAchr[sampled_idx_find_satisfied][len(toplevelsnptable_titlelist)+vcftable_idx*2]))*(round(noofindvds2quantizing/len(options.quantizing))/int(allsnpOfJoinTableinAchr[sampled_idx_find_satisfied][len(toplevelsnptable_titlelist)+vcftable_idx*2+1]))))#(AN-AC)*(n/x)*(AN)
                    else:#FIXED AS REF,the no of ref alle is vcftableslist[vcftable_idx][1] that is arbitrarily assigned,but dadi will projection it into 26,so any no greater than 26 is ok
                        print(vcftableslist[vcftable_idx][1],end="\t",file=dadisnpfile)
                        if options.quantizing!=None and vcftable_idx in vcftablesidxlist_toquantizing:
                            k+=int(round(int(vcftableslist[vcftable_idx][1])*(round(noofindvds2quantizing/len(options.quantizing))/int(vcftableslist[vcftable_idx][1]))))
                if options.quantizing!=None:
                    print(str(k),end="\t",file=dadisnpfile)
                        
                print(ALT,end="\t",file=dadisnpfile)
                for vcftable_idx in range(len(vcftableslist)):
                    if allsnpOfJoinTableinAchr[sampled_idx_find_satisfied][len(toplevelsnptable_titlelist)+vcftable_idx*2]!=None:
                        print(str(int(allsnpOfJoinTableinAchr[sampled_idx_find_satisfied][len(toplevelsnptable_titlelist)+vcftable_idx*2])),end="\t",sep="\t",file=dadisnpfile)#AC
                        if options.quantizing!=None and vcftable_idx in vcftablesidxlist_toquantizing:
                            l+=int(round(int(allsnpOfJoinTableinAchr[sampled_idx_find_satisfied][len(toplevelsnptable_titlelist)+vcftable_idx*2])*(round(noofindvds2quantizing/len(options.quantizing))/int(allsnpOfJoinTableinAchr[sampled_idx_find_satisfied][len(toplevelsnptable_titlelist)+vcftable_idx*2+1]))))  
                    else:
                        print("0",end="\t",file=dadisnpfile)
                        if options.quantizing!=None and vcftable_idx in vcftablesidxlist_toquantizing:
                            l+=0
                if options.quantizing!=None:
                    print(str(l),end="\t",file=dadisnpfile)
                print(currentchrID.replace(".","_"),postion,sep="\t",file=dadisnpfile)
        if options.countsnpnumberfromvcf!=None:
            totalduiltsnp+=duiltedsnpcountforcurrentchr
            if duiltedsnpcountforcurrentchr!=0:
                print(currentchrID,file=randomnamef_recchr)
            print(currentchrID,"totalduiltsnp",totalduiltsnp,"duiltedsnpcountforcurrentchr",duiltedsnpcountforcurrentchr)
        break
    exit()

    dadisnpfile.close()
    randomnamef_recchr.close()
    if options.countsnpnumberfromvcf!=None:
        a=os.popen("awk '$0~/#/ || length($4)==length($5) {print $0}' "+options.countsnpnumberfromvcf+" |grep -wFf "+randomnamefile_recordchr+"  - |grep ^[^#] - |less -S|wc -l")
        totalsnpforcount=int(a.readline().strip())
        print("should be none",a.readline())
        a.close()
    else:
        totalsnpforcount=allsnps
    
    
    dulttimes=totalduiltsnp/totalsnpforcount
    print(options.outputfilename)
    print("totalduiltsnp",totalduiltsnp,"totalsnpforcount",totalsnpforcount,totallength,"real totallengthduilt",dulttimes*totallength,"totallengthduilt",totallengthduilt,"finsih")
    os.system("rm "+randomnamefile_recordchr)