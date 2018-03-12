'''
Created on 2018-3-6

@author: Dr.liu
'''
# -*- coding: utf-8 -*-
"""
Example of a dynamic list (add, remove button handled on client side) with data stored in a classic String in database.
Works with Flask, Flask-WTForms, flask-SQLAlchemy (uses SQL lite for easy testing)
Inspired from https://gist.github.com/kageurufu/6813878
"""
import itertools
from flask import Flask, render_template_string
from flask import request
from flask_wtf import FlaskForm
from wtforms import FieldList, StringField
from wtforms import validators as wtf_validators
from wtforms.utils import unset_value
from flask_sqlalchemy import SQLAlchemy

app = Flask(__name__)
app.secret_key = 'TEST'

app.config['SQLALCHEMY_DATABASE_URI'] = r'sqlite:////test.db' # Will create a file in pwd (current path)
db = SQLAlchemy(app)

_max_nb_entries = 100
_max_len_per_entry = 30

class Crew(db.Model):
    __tablename__ = 'crew'
    id = db.Column(db.Integer, primary_key=True)
    persons = db.Column(db.String(_max_nb_entries*_max_len_per_entry),default="",nullable=False)

@app.before_request
def before_request():
    db.create_all()


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


     
class TestForm(FlaskForm):
    persons = FieldListFromString(StringField('Persons',default='',validators=[wtf_validators.Length(min=0, max=_max_len_per_entry)]),
                                  min_entries=1, max_entries=_max_nb_entries)


@app.route('/', methods=['POST', 'GET'])
def example():
    print('request.form',request.form)
    # this paragraph is trick for testing: (merging modification and add new)
    existing_crew = Crew.query.get(1)
    if existing_crew is None: 
        crew = Crew()
    else:
        crew = existing_crew
        
    form = TestForm(obj=crew)
    
    if form.validate_on_submit():
        form.populate_obj(crew)
        
        if existing_crew is None:# this trick for testing: (merging modification and add new)
            db.session.add(crew)
            
        db.session.commit()

    return render_template_string(
        """
<form method="post" action="{{url_for('example')}}">
  
    {{ form.hidden_tag() }}
    <div data-toggle="fieldset" id="persons-fieldset">
        <button type="button" data-toggle="fieldset-add-row" data-target="#persons-fieldset">+</button>
        
        <ul id="persons">
            {% for person in form.persons %}
            <li data-toggle="fieldset-entry" >
                <label for="persons-{{ loop.index0 }}">Person</label>
                <input id="persons-{{ loop.index0 }}" name="persons-{{ loop.index0 }}" type="text" value="{{ person.data }}">
                <button type="button" data-toggle="fieldset-remove-row" id="persons-{{loop.index0}}-remove">-</button>
            </li>
            {% endfor %}
        </ul>
    </div>
    
    <input type="submit"/>
</form>
<!-- just to display all potential errors -->
{% if form.errors %}
    <ul class="errors">
        {% for field_name, field_errors in form.errors|dictsort if field_errors %}
            {% for error in field_errors %}
                <li>{{ form[field_name].label }}: {{ error }}</li>
            {% endfor %}
        {% endfor %}
    </ul>
{% endif %}
<!-- The imports -->
<script src="//ajax.googleapis.com/ajax/libs/jquery/3.2.1/jquery.min.js"></script>
<script>
$(function() {
    $("div[data-toggle=fieldset]").each(function() {
        var $this = $(this);
            
        //Add new entry
        $this.find("button[data-toggle=fieldset-add-row]").click(function() {
            var target = $($(this).data("target"));
            //console.log(target);
            var oldrow = target.find("[data-toggle=fieldset-entry]:last");
            var row = oldrow.clone(true, true);
            var elem_id = row.find(":input")[0].id;
            var elem_prefix = elem_id.replace(/(.*)-(\d{1,4})/m, '$1')// max 4 digits for ids in list
            var elem_num = parseInt(elem_id.replace(/(.*)-(\d{1,4})/m, '$2')) + 1;
            //console.log(elem_prefix);
            //console.log(elem_num);
            row.children(':input').each(function() {
                var id = $(this).attr('id').replace(elem_prefix+'-' + (elem_num - 1), elem_prefix+'-' + (elem_num));
                $(this).attr('name', id).attr('id', id).val('').removeAttr("checked");
            });
            row.children('label').each(function() {
                var id = $(this).attr('for').replace(elem_prefix+'-' + (elem_num - 1), elem_prefix+'-' + (elem_num));
                $(this).attr('for', id);
            });
            row.show();
            oldrow.after(row);
        }); //End add new entry
        //Remove row
        $this.find("button[data-toggle=fieldset-remove-row]").click(function() {
            if($this.find("[data-toggle=fieldset-entry]").length > 1) {
                var thisRow = $(this).closest("[data-toggle=fieldset-entry]");
                thisRow.remove();
            }
        }); //End remove row
    });
});
</script>
""", form=form)


if __name__ == '__main__':
    app.run(host="localhost",port=5001,debug = True)