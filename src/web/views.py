'''
Created on 2017年10月19日

@author: liurui
'''
from email.utils import unquote
import os
import re

from flask import request, jsonify, send_from_directory, abort, render_template
from werkzeug.routing import BaseConverter

from web import web
from web.forms import ParaForm,UsersForm


class RegexConverter(BaseConverter):
    def __init__(self, map, *args):
        self.map = map
        self.regex = args[0]
web.url_map.converters['regex'] = RegexConverter

@web.route('/mytest',methods=['GET','POST'])        
def testmyform():
    form=ParaForm()
    print("error",form.errors)
    if request.method=='POST'  :
        print("post",form.projectpath.data)
        if form.is_submitted():
            print("here")
            pass
#             return form.projectpath.data+"submitted"
        if form.validate_on_submit():
            print(request.form['projectpath'],"there\n",os.getcwd())
            return form.datadepth.data+"come on"

    return render_template('configsotware.html',form=form)
    

# @web.route('/login',methods=['GET','POST'])
# def login():
#     return render_template('hello.html')
# 
# @web.route('/hello',methods=['POST'])
# def hello():
#     ccc=request.form['textfield']
#     print(ccc,"there\n",os.getcwd())
#     return request.form['textfield']+"come on"

@web.route('/',methods=['GET','POST'])
def configsoftware():
    form=ParaForm()
    if form.validate_on_submit():
        print("configsoftware here",request.values)
        return form.datadepth.data+form.projectpath.data+"<br />"+form.outputpath.data+form.outputperfix.data+"<br />"+"".join(form.foldername.data)
    else:
        print("didn't validate")
        return render_template('commandtemplate.html',form=form)





@web.route('/ttt', methods=['GET', 'POST'])
def hello_world():

    form = UsersForm(users=[{}, {}, {}])
    form.users[0].experience.choices=[(1, 'One'), (2, 'Two')]
    form.users[1].experience.choices = [(1, 'Uno'), (2, 'Du')]
    form.users[2].experience.choices = [(0, 'Zero')]
    return render_template("dynamic.html", form=form)
    
@web.route('/downloadfile/:urlpath#.+#')
def send_static(urlpath):
    print("send_static")
    filename=re.search(r'[^/]*$',unquote(urlpath)).group(0)
    path="../../"+re.search(r'^.*/',unquote(urlpath)).group(0)
    print(path,filename,unquote(urlpath))
#    print(urlpath,re.search(r'^.*/',urlpath).group(0),re.search(r'[^/]*$',urlpath).group(0))
    return render_template(path+filename)