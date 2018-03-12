# -*- coding: UTF-8 -*-
'''
Created on 2017年10月19日

@author: liurui
'''

from flask_wtf import FlaskForm
from flask_wtf.form import Form
from wtforms import FieldList, StringField, SubmitField, SelectField,HiddenField
from wtforms.fields.core import FormField
from wtforms.validators import Required, DataRequired
import itertools
from web import web
from wtforms.utils import unset_value
from wtforms import validators as wtf_validators


class FileUploadForm(Form):
    pass

class StudentForm(Form):
    student_id = StringField('Student ID', validators = [DataRequired()])
    student_name = StringField('Student Name', validators = [DataRequired()])

class AddClassForm(Form):
    name = StringField('classname', validators=[DataRequired()])
    day = SelectField('classday', 
                      choices=[(1,"Monday"),(2,"Tuesday"),(3,"Wednesday"),(4,"Thursday"),(5,"Friday")],
                      coerce=int)

    students = FieldList(FormField(StudentForm), min_entries = 5) # show at least 5 blank fields by default

class UserForm(FlaskForm):
# must inherit from wtforms.Form, not flask-WTForms'
# see http://stackoverflow.com/questions/15649027/wtforms-csrf-flask-fieldlist
    first_name = StringField('First Name')
    last_name = StringField('Last Name')

    experience = SelectField('Experience', coerce=int)
class UsersForm(Form):
    users = FieldList(FormField(UserForm), min_entries=2)
_max_nb_entries = 100
_max_len_per_entry = 30  
_delimiter = "#;_"
class FieldListFromString(FieldList):
    """
    The idea here is to have a FieldList but to store the data in a string format instead of a list
    """
    def process(self, formdata, data=unset_value):
        self.entries = []
        if data is unset_value or not data:
            try:
                data = self.default()
            except TypeError:
                data = self.default
                
        ## Modification from classic FieldList
        if data and 0<len(data):
            data = data.split(_delimiter)
        else:
            data = []
        #
            
        self.object_data = data

        if formdata:
            indices = sorted(set(self._extract_indices(self.name, formdata)))
            if self.max_entries:
                indices = indices[:self.max_entries]

            idata = iter(data)
            for index in indices:
                try:
                    obj_data = next(idata)
                except StopIteration:
                    obj_data = unset_value
                self._add_entry(formdata, obj_data, index=index)
        else:
            for obj_data in data:
                self._add_entry(formdata, obj_data)

        while len(self.entries) < self.min_entries:
            self._add_entry(formdata)


    def populate_obj(self, obj, name):
        values = getattr(obj, name, None)
        try:
            ivalues = iter(values)
        except TypeError:
            ivalues = iter([])

        candidates = itertools.chain(ivalues, itertools.repeat(None))
        _fake = type(str('_fake'), (object, ), {})
        output = []
        for field, data in zip(self.entries, candidates):
            fake_obj = _fake()
            fake_obj.data = data
            field.populate_obj(fake_obj, 'data')
            output.append(fake_obj.data)
        
        ## Modification from classic FieldList
        setattr(obj, name, _delimiter.join(output))

class ParaForm(FlaskForm):
    projectpath=StringField("请输入项目路径：",validators=[Required(message="根目录不能为空")])
    datadepth=StringField("数据所在在层级：",validators=[Required()])
    outputpath=StringField("输出路径：")
    outputperfix=StringField("输出后缀：")
    tag=HiddenField("tag",id='foldtag')
    persons = FieldListFromString(StringField('需过滤的文件夹名:',default='',validators=[wtf_validators.Length(min=0, max=_max_len_per_entry)]),
                                  min_entries=1, max_entries=_max_nb_entries)
#     foldername=FieldList(StringField("org"),label='需过滤的文件夹名:',min_entries=2)
#     sftWname=StringField("请输入所在层级：",validators=[Required()])
#     sfware=
#     inputsuffix=
    
    submit = SubmitField('提交')