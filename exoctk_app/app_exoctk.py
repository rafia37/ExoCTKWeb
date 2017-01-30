from flask import Flask, render_template, request, redirect, make_response

app_exoctk = Flask(__name__)

import ExoCTK
import pandas
import numpy as np
from bokeh.plotting import figure
from bokeh.embed import components
from bokeh.models import ColumnDataSource, HoverTool, OpenURL, TapTool
from bokeh.models.widgets import Panel, Tabs

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
    
    # Get the input
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
    coeff, mu, r = ExoCTK.ldc.ldcfit.ldc(teff, logg, feh, model_grid, profile, bandpass=bandpass)
    
    # Plotting
    tools = "resize,crosshair,pan,wheel_zoom,box_zoom,reset"
    x = 'Wavelength (um)'
    y = 'Limb Darkening'
    
    # can specify plot_width if needed
    p = figure(tools=tools, title='Foo', x_axis_label=x, y_axis_label=y, plot_width=800)

    # Evaluate the limb darkening profile
    mu_vals = np.linspace(0, 1, 1000)
    ldfunc = ExoCTK.ldc.ldcfit.ld_profile(profile)
    ld_vals = ldfunc(mu_vals, *coeff)
    p.line(mu_vals, ld_vals, line_width=2)

    script, div = components(p)

    return render_template('results.html', coeff=coeff, mu=mu, r=r, script=script, plot=div, version=VERSION)

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
