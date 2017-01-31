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
    
    # Make the matplotlib plot into a Bokeh plot
    p = mpl.to_bokeh()
    p.y_range = Range1d(0, 1)
    p.x_range = Range1d(0, 1)
    
    # p.xaxis.axis_label = LatexLabel(text=plt.gca().xaxis.get_label().get_text(),
    #                      x_units='screen', y_units='screen', render_mode='css', text_font_size='16pt')
    # p.yaxis.axis_label = LatexLabel(text=plt.gca().yaxis.get_label().get_text(),
    #                      x_units='screen', y_units='screen', render_mode='css', text_font_size='16pt')

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

JS_CODE = """
import {Label, LabelView} from "models/annotations/label"

export class LatexLabelView extends LabelView
  render: () ->

    #--- Start of copied section from ``Label.render`` implementation

    ctx = @plot_view.canvas_view.ctx

    # Here because AngleSpec does units tranform and label doesn't support specs
    switch @model.angle_units
      when "rad" then angle = -1 * @model.angle
      when "deg" then angle = -1 * @model.angle * Math.PI/180.0

    if @model.x_units == "data"
      vx = @xmapper.map_to_target(@model.x)
    else
      vx = @model.x
    sx = @canvas.vx_to_sx(vx)

    if @model.y_units == "data"
      vy = @ymapper.map_to_target(@model.y)
    else
      vy = @model.y
    sy = @canvas.vy_to_sy(vy)

    if @model.panel?
      panel_offset = @_get_panel_offset()
      sx += panel_offset.x
      sy += panel_offset.y

    #--- End of copied section from ``Label.render`` implementation

    # ``katex`` is loaded into the global window at runtime
    # katex.renderToString returns a html ``span`` element
    latex = katex.renderToString(@model.text, {displayMode: true})

    # Must render as superpositioned div (not on canvas) so that KaTex
    # css can properly style the text
    @_css_text(ctx, latex, sx + @model.x_offset, sy - @model.y_offset, angle)

export class LatexLabel extends Label
  type: 'LatexLabel'
  default_view: LatexLabelView
"""

class LatexLabel(Label):
    """A subclass of the Bokeh built-in `Label` that supports rendering
    LaTex using the KaTex typesetting library.

    Only the render method of LabelView is overloaded to perform the
    text -> latex (via katex) conversion. Note: ``render_mode="canvas``
    isn't supported and certain DOM manipulation happens in the Label
    superclass implementation that requires explicitly setting
    `render_mode='css'`).
    """
    __javascript__ = ["https://cdnjs.cloudflare.com/ajax/libs/KaTeX/0.6.0/katex.min.js"]
    __css__ = ["https://cdnjs.cloudflare.com/ajax/libs/KaTeX/0.6.0/katex.min.css"]
    __implementation__ = JS_CODE