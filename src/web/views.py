'''
Created on 2017年10月19日

@author: liurui
'''
from web import web
from flask import request,jsonify,send_from_directory,abort,render_template
import os
from werkzeug.routing import BaseConverter
from web.forms import ParaForm


class RegexConverter(BaseConverter):
    def __init__(self, map, *args):
        self.map = map
        self.regex = args[0]
web.url_map.converters['regex'] = RegexConverter

@web.route('/mytest',methods=['GET','POST'])        
def testmyform():
    if not request.form['textfield']  :
        print(request.args.get('textfield'),"there\n",os.getcwd())
        return "come on"
    else:
        return render_template('configsotware.html')

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
        print("here")
    else:
        return render_template('processconfigure.html')
    
