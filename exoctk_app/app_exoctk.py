from flask import Flask, render_template, request, redirect, make_response

app_exoctk = Flask(__name__)

@app_exoctk.route('/')
@app_exoctk.route('/index')

def index():
     return render_template('index.html')

# Page with a text box to take the SQL query
@app_exoctk.route('/ldc', methods=['GET', 'POST'])
def exoctk_ldc():
    return render_template('ldc.html')
