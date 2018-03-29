'''
Created on 2014-11-17

@author: liurui
'''
import datetime
import re
import time

import markdown2
from sqlalchemy.orm import session
from tabulate import tabulate

from src.web import entity
from src.web.dba import addJobs2jobstate
import src.web.dba as aaa


ISOTIMEFORMAT = '%Y-%m-%d %X'
if __name__ == '__main__':
#     ll=["addJobs","ddddd"]
#     addJobs(ll,"kkkk")
#     llll=["test","luowen"]
#     addJobs(llll,"lhomelcomlfuc")
#     
#     print("update jobsstat set startdate='"+time.strftime(ISOTIMEFORMAT, time.localtime()) +"' where id='1'")
#     print(datetime.datetime.now())
#     
#     results=[]
#     
#     file=open("F:\work\pyhtmlmarkdown\\tttt.txt",'r')
#     bowtieout=file.read()
#     file.close()
#     results.append(session.execute("update jobsstat set startdate='"+time.strftime(ISOTIMEFORMAT, time.localtime()) +"' where id='1'"))
#     results.append(session.execute("update jobsstat set outputinfo='"+bowtieout+"' where id = '1'"))
#     print("results",results)
    session=aaa.getWebSession()
    l = session.query(entity.Jobs_recoder).all()
    header=["*scriptname*","*foldername*","*starttime*","*finishtime*"," *state*","*outputinfo*"]
    mylist=[]
    
    for i in l:
        print("sssssssssssssssss",i.outputinfo)
#         mylist.append([("<br>"+i.scriptname),i.foldername[10:],("&nbsp;"+str(i.startdate)+"&nbsp;"),("&nbsp;"+str(i.finishdate)+"&nbsp;"),("&nbsp;"+str(i.state)),("""<input type="button" value="outputinfo" onclick="location.href='http://www.baidu.com'">""")])
        mylist.append([i.scriptname,i.foldername[12:],(str(i.startdate)),(str(i.finishdate)),(str(i.state)),("""<input type="button" value="outputinfo" onclick="location.href='http://www.baidu.com'">""")])
    print(mylist)
    print("======orgtbl=====================")
    print(tabulate(mylist,header,tablefmt="markdown"))    
    text=tabulate(mylist,header,tablefmt="orgtbl")
    text=re.sub('\|\|[\-\+]+\|\|\n', '', text.replace("|", "||"))
    print(text)

    html=markdown2.markdown(text,extras=["wiki-tables"])
    t="""
           <form name="myform" method="post" action="">
       <select name="software" onchange="change(this.value)">
           <option value="1">Bowtie2</option>
           <option value="2">SortSam.jar</option>
           <option value="3">MarkDuplicates.jar</option>
           <option value="4">GenomeAnalysisTK.jar</option>
       </select><br />
    """
    html=html.replace("<tr", "<tr bgcolor='lightgrey'",1)
    print(html,file=open("F:\work\pyhtmlmarkdown\mywikitable.html",'w'))
#    session=aaa.getSession()
 #   session.execute("update jobsstat set finishdate='"+time.strftime(ISOTIMEFORMAT, time.localtime()) +"' where id='4'")
    