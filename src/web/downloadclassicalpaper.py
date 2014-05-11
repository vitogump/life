'''
Created on 2014-5-1

@author: liurui
'''
from bottle import route, run, template, get, post, request, static_file
import os
import src.web.dba as dba
import shutil
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
@route('/classicalpaper')
def 
@route('/download/<filename>')
def server_static(filename):
    return static_file(filename,root='../../index')
            
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

@route('/classical_paper/:filename')
def send_static(filename):
    return static_file(filename, root='../../classical_paper',download=filename)
#upload module

@post('/upload')
def do_upload():
    
    classname = request.forms.get('radio1')
    filename=request.forms.get("filename")
    papername=request.forms.get("papername")
    data = request.files.get('data')
    print("radio1","filename",filename,papername)
    if filename and data.file and classname:
#        raw = data.file.read() #当文件很大时，这个操作将十分危险
        filename = data.filename
        with open(os.path.join(UPLOAD_BASE, filename), 'wb') as f:
            dba.addArticle(filename,classname)
            shutil.copyfileobj(data.file, f, 8192)
        return "Hello {}! You uploaded {} ( bytes).".format(classname, filename)
    return "You missed a field"



run(host='localhost',port=8080,debug=True)