# -*- coding:utf8 -*-
'''
Created on 2014-5-1

@author: liurui
'''
from bottle import route, run, template, get, post, request, static_file
from urllib.parse import quote, unquote
import os
import re
import shutil
# import src.web.DBA as dba
import string
UPLOAD_BASE = "../../classical_paper"
@get('/login')
def login_formc():
    return'''<form method = "POST" action="/login">
            <input name="name" type="text"/>
            <input name="password" type="password"/>
            <input type="submit"/>
            
            </form>
            '''
#@route('/static/<filename>')
#def server_static(filename):
#    return static_file(filename, root='statichtml')
#@route('/classicalpaper')
#def 
@route('/download/:urlpath#.+#')
def server_static(urlpath):
    print("sever_static")
    herfs=""""""
    catajudge=re.search(r'[^/]*$',urlpath).group(0)
    if urlpath.endswith("GBS"):
        path='../../toDownload/classical_paper/GBS/'
    elif urlpath.endswith("mRNA_microRNA"):
        path='../../toDownload/classical_paper/mRNA_microRNA/'
    elif urlpath.endswith("imprinting_Epigenetic"):
        path='../../toDownload/classical_paper/imprinting_epigenetic/'
    elif urlpath.endswith("selection"):
        path='../../toDownload/classical_paper/naturalselection_artificialselection/'
    elif urlpath.endswith("GeneMapping"):
        path='../../toDownload/classical_paper/geneMapping/'
    elif urlpath.endswith("aglorithms_in_biology"):
        path='../../toDownload/ppt_for_share/aglorithms_in_biology/'
    elif urlpath.endswith("freshtoread"):
        path='../../toDownload/ppt_for_share/freshtoread/'
    elif urlpath.endswith("genomics"):
        path='../../toDownload/ppt_for_share/genomics/'
    elif urlpath.endswith("renbin"):
        path='../../toDownload/ppt_for_share/renbin/'
    elif urlpath.endswith("Python_study"):
        path='../../toDownload/pythonstudy/'
    elif urlpath.endswith(".html"):
        print("sssss")
        return static_file(urlpath,root='../../index')
    else:
        print(unquote(urlpath))
        path = "../../"+re.search(r'catalog/(.*)',unquote(urlpath)).group(1)+"/"
        print(path)
#    else:#download file
#        path="../../"+re.search(r'^.*/',unquote(urlpath)).group(0)
#        print(urlpath,re.search(r'^.*/',urlpath).group(0),re.search(r'[^/]*$',urlpath).group(0))
#        return static_file(re.search(r'[^/]*$',urlpath).group(0), root='../../classical_paper/'+re.search(r'^.*/',urlpath).group(0),download=re.search(r'[^/]*$',urlpath).group(0))
    l=os.listdir(path=path)
    print(l,path)
    for a in l:
        if os.path.isdir(path+a) and not urlpath.endswith("rar"):
            url=quote('/download/catalog/'+re.search(r'\.\./\.\./(.*)',path+a).group(1))
            herfs+="""<a href="""+url+""" target="_self">"""+a+"""</a></p>
            """
        elif os.path.isfile(path+a):
            url=quote('/downloadfile/'+re.search(r'\.\./\.\./(.*)',path+a).group(1))
            herfs+="""<a href="""+url+""" target="_self">"""+a+"""</a></p>
            """
    return herfs      
#    return static_file(filename,root='../../index')
            
@post('/login')
def login_submit():
    name = request.forms.get("name")
    password = request.forms.get("password")
    if name =="liu" and password =="123":
        return"<p> login sucessed</p>"
    else:
        return"<p>login failed</p>"
@route('/hello/:name')
def greet(name='Stranger'):
    return 'Hello {},how are you?'.format(name)

@route('/downloadfile/:urlpath#.+#')
def send_static(urlpath):
    print("send_static")
    filename=re.search(r'[^/]*$',unquote(urlpath)).group(0)
    path="../../"+re.search(r'^.*/',unquote(urlpath)).group(0)
    print(path,filename,unquote(urlpath))
#    print(urlpath,re.search(r'^.*/',urlpath).group(0),re.search(r'[^/]*$',urlpath).group(0))
    return static_file(filename, root=path,download=filename)
#upload module

# @post('/upload')
# def do_upload():
#     
#     classname = request.forms.get('radio1')
#     filename=request.forms.get("filename")
#     papername=request.forms.get("papername")
#     data = request.files.get('data')
#     print("radio1","filename",filename,papername)
#     if filename and data.file and classname:
# #        raw = data.file.read() #当文件很大时，这个操作将十分危险
#         filename = data.filename
#         with open(os.path.join(UPLOAD_BASE, filename), 'wb') as f:
#             dba.addArticle(filename,classname)
#             shutil.copyfileobj(data.file, f, 8192)
#         return "Hello {}! You uploaded {} ( bytes).".format(classname, filename)
#     return "You missed a field"



run(host='192.168.0.14',port=5000)