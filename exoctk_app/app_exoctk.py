from flask import Flask, render_template, request, redirect, make_response

app_exoctk = Flask(__name__)

import ExoCTK
import pandas
import numpy as np
import matplotlib.pyplot as plt
import bokeh
from bokeh.embed import components
from bokeh.models import Range1d, Label
from bokeh import mpl
from bokeh.core.properties import Override
from bokeh.models import Label
# from bokeh.plotting import figure
# from bokeh.models import ColumnDataSource, HoverTool, OpenURL, TapTool
# from bokeh.models.widgets import Panel, Tabs

# Redirect to the index
VERSION = ExoCTK.__version__
@app_exoctk.route('/')
@app_exoctk.route('/index')

# Load the Index page
def index():
    return render_template('index.html', version=VERSION)

# Load the LDC page
@app_exoctk.route('/ldc', methods=['GET', 'POST'])
def exoctk_ldc():
    return render_template('ldc.html', version=VERSION)
    
# Load the results page
@app_exoctk.route('/results', methods=['GET', 'POST'])
def exoctk_results():
    
    # Get the input from the form
    modeldir = request.form['modeldir']
    wave_rng = (float(request.form['wave_min']),float(request.form['wave_max']))
    Teff_rng = (float(request.form['teff_min']),float(request.form['teff_max']))
    logg_rng = (float(request.form['logg_min']),float(request.form['logg_max']))
    FeH_rng = (float(request.form['feh_min']),float(request.form['feh_max']))
    profile = request.form['profile']
    bandpass = request.form['bandpass']
    teff = int(request.form['teff'])
    logg = float(request.form['logg'])
    feh = float(request.form['feh'])
    
    # Make the model grid
    model_grid = ExoCTK.core.ModelGrid(modeldir)
    model_grid.customize(Teff_rng=Teff_rng, logg_rng=logg_rng, FeH_rng=FeH_rng, wave_rng=wave_rng)

    # Calculate the coefficients
    coeff, mu, r = ExoCTK.ldc.ldcfit.ldc(teff, logg, feh, model_grid, profile, bandpass=bandpass, plot=True)
    
    # Make some HTML for the input summary
    coeff_text = ', '.join(['c{} = {:.3f}'.format(n+1,c) for n,c in enumerate(coeff)])
    summary = """<h3>Input</h3>
    <table style="border-spacing: 10px; border-collapse: separate;">
    <tr>
        <td>Teff</td>
        <td>{}</td>
    </tr>
    <tr>
        <td>log(g)</td>
        <td>{}</td>
    </tr>
    <tr>
        <td>Fe/H</td>
        <td>{}</td>
    </tr>
    <tr>
        <td>Bandpass</td>
        <td>{}</td>
    </tr>
    <tr>
        <td>Profile</td>
        <td>{}</td>
    </tr>
    </table>
    <h3>Result</h3>
    <table style="border-spacing: 10px; border-collapse: separate;">
    <tr>
        <td>Coefficients</td>
        <td>{}</td>
    </tr>
    <tr>
        <td>mu</td>
        <td>{:.3f}</td>
    </tr>
    <tr>
        <td>Effective Radius</td>
        <td>{}</td>
    </tr>
    
    </table>
    """.format(teff, logg, feh, bandpass or 'None', profile, coeff_text, mu, int(r))
    
    # Make the matplotlib plot into a Bokeh plot
    p = mpl.to_bokeh()
    p.y_range = Range1d(0, 1)
    p.x_range = Range1d(0, 1)
    script, div = components(p)

    return render_template('results.html', summary=summary, script=script, plot=div, version=VERSION)

# Save the results to file
@app_exoctk.route('/savefile', methods=['POST'])
def exoctk_savefile():
    export_fmt = request.form['format']
    if export_fmt == 'votable':
        filename = 'exoctk_table.vot'
    else:
        filename = 'exoctk_table.txt'

    response = make_response(file_as_string)
    response.headers["Content-Disposition"] = "attachment; filename=%s" % filename
    return response