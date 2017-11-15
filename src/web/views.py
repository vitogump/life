'''
Created on 2017年10月19日

@author: liurui
'''
from web import web
from flask import request,jsonify,send_from_directory,abort,render_template
import os
from werkzeug.routing import BaseConverter
from flask.ext.wtf import Form
from wtforms import StringField,SubmitField
from wtforms.validators import Required

class ParaForm(Form):
    projectpath=StringField("请输入项目路径：",validators=[Required()])
    datadepth=StringField("请输入所在层级：",validators=[Required()])
#     sftWname=StringField("请输入所在层级：",validators=[Required()])
    submit = SubmitField('Submit')

class RegexConverter(BaseConverter):
    def __init__(self, map, *args):
        self.map = map
        self.regex = args[0]
web.url_map.converters['regex'] = RegexConverter

@web.route('/',methods=['GET','POST'])
def configsoftware():
    form=ParaForm()
    if form.validate_on_submit():
        print("here")
    else:
        return render_template('processconfure.html')
        
    #else:
        print("there",os.getcwd())
        return send_from_directory(os.getcwd(),'..\\..\\toDownload\\pythonstudy\\MITjisuanjikexuebianchengdaolunlec01.mp4',as_attachment=True)