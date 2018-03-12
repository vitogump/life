from flask import Flask, render_template
from flask_wtf import Form
from wtforms import FormField, HiddenField, FieldList, IntegerField
from wtforms.validators import DataRequired

flask = Flask(__name__)
flask.secret_key = "I will not buy this record, it is scratched."

class InnerForm(Form):
    hidden = HiddenField(validators=[DataRequired()])
    integer = IntegerField(validators=[DataRequired()])

class OuterForm(Form):
    innerforms = FieldList(FormField(InnerForm))

@flask.route("/")
def index():
    form = OuterForm()
    for i in range(0,2):
        form.innerforms.append_entry(data={ 'hidden' : str(i), 'integer' : i})
    return render_template("demo.html", form=form)

flask.run("localhost", 5008, debug=True)