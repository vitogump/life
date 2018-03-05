# -*- coding: UTF-8 -*-
'''
Created on 2017年10月19日

@author: liurui
'''
from flask_wtf import FlaskForm
from flask_wtf.form import Form
from wtforms import FieldList, StringField, SubmitField, SelectField
from wtforms.fields.core import FormField
from wtforms.validators import Required

from web import web


class UserForm(FlaskForm):
# must inherit from wtforms.Form, not flask-WTForms'
# see http://stackoverflow.com/questions/15649027/wtforms-csrf-flask-fieldlist
    first_name = StringField('First Name')
    last_name = StringField('Last Name')

    experience = SelectField('Experience', coerce=int)
class UsersForm(Form):
    users = FieldList(FormField(UserForm), min_entries=2)
class ParaForm(FlaskForm):
    projectpath=StringField("请输入项目路径：",validators=[Required(message="根目录不能为空")])
    datadepth=StringField("数据所在在层级：",validators=[Required()])
    outputpath=StringField("输出路径：")
    outputperfix=StringField("输出后缀：")
    foldername=FieldList(StringField("org"),label='需过滤的文件夹名:',min_entries=2)
#     sftWname=StringField("请输入所在层级：",validators=[Required()])
#     sfware=
#     inputsuffix=
    
    submit = SubmitField('提交')