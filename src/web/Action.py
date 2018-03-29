# -*- coding:utf8 -*-
'''
Created on 2014-11-17

@author: liurui
'''
import os, re, shutil, string
from time import sleep
from urllib.parse import quote, unquote

from bottle import route, run, template, get, post, request, static_file
import markdown2
from tabulate import tabulate

from src.web import entity
import src.web.dba as aaa




UPLOAD_BASE = "../../com"


@route('/monitorjobs/:urlpath#.+#')
def server_static(urlpath):
    session=aaa.getWebSession()
    l = session.query(entity.Jobs_recoder).all()
    header=["*scriptname*","*foldername*","*starttime*","*finishtime*"," *state*","*outputinfo*"]
    mylist=[]
    
    for i in l:
        mylist.append([("<br>"+i.scriptname),i.foldername[10:],("&nbsp;"+str(i.startdate)+"&nbsp;"),("&nbsp;"+str(i.finishdate)+"&nbsp;"),("&nbsp;"+str(i.state)),("""<input type="button" value="outputinfo" onclick="location.href='http://www.baidu.com'">""")])
    text=tabulate(mylist,header,tablefmt="markdown2")
    html=markdown2.markdown(text,extras=["wiki-tables"])
    print(html,file=open("../../com/temppage/wikitable.html",'w'))
    return static_file(urlpath,root='../../com/temppage')

run(host='localhost',port=8182,debug=True)
