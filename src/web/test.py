'''
Created on 2014-11-17

@author: liurui
'''
import datetime
import time

import markdown2
from sqlalchemy.orm import session
from tabulate import tabulate

from src.web import Entity
from src.web.DBA import addJobs
import src.web.DBA as aaa


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
    session=aaa.getSession()
    l = session.query(Entity.Jobstat).all()
    header=["*scriptname*","*foldername*","*starttime*","*finishtime*"," *state*","*outputinfo*"]
    mylist=[]
    
    for i in l:
        print("sssssssssssssssss",i.outputinfo)
        mylist.append([i.scriptname,i.foldername[10:],("&nbsp;"+str(i.startdate)+"&nbsp;"),("&nbsp;"+str(i.finishdate)+"&nbsp;"),("&nbsp;"+str(i.state)),("""<input type="button" value="outputinfo" onclick="location.href='http://www.baidu.com'">""")])
    print(mylist)
    print("======orgtbl=====================")
    print(tabulate(mylist,header,tablefmt="orgtbl"))    
    text=tabulate(mylist,header,tablefmt="markdown2")
    print(text)

    html=markdown2.markdown(text,extras=["wiki-tables"])
    print(html,file=open("F:\work\pyhtmlmarkdown\wikitable.html",'w'))
#    session=aaa.getSession()
 #   session.execute("update jobsstat set finishdate='"+time.strftime(ISOTIMEFORMAT, time.localtime()) +"' where id='4'")
    