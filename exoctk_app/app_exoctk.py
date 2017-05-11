## -- IMPORTS
import os

import astropy.table as at
import matplotlib.pyplot as plt
import numpy as np

import bokeh
from bokeh import mpl
from bokeh.core.properties import Override
from bokeh.models import ColumnDataSource
from bokeh.embed import components
from bokeh.models import HoverTool
from bokeh.models import Label
from bokeh.models import Range1d
from bokeh.plotting import figure
from bokeh.plotting import show
from flask import Flask
from flask import render_template
from flask import request
from flask import redirect
from flask import make_response
from flask import current_app
from flask_cache import Cache

import ExoCTK

## -- FLASK SET UP (?)
app_exoctk = Flask(__name__)

# define the cache config keys, remember that it can be done in a settings file
app_exoctk.config['CACHE_TYPE'] = 'simple'

# register the cache instance and binds it on to your app
cache = Cache(app_exoctk)

# Redirect to the index
VERSION = ExoCTK.__version__
@app_exoctk.route('/')
@app_exoctk.route('/index')

## -- FUNCITONS

# Load the Index page
def index():
    return render_template('index.html')

# Load the LDC page
@app_exoctk.route('/ldc', methods=['GET', 'POST'])
def exoctk_ldc():
    # Get all the available filters
    filters = ExoCTK.core.filter_list()['bands']
    
    # Make HTML for filters
    filt_list = '\n'.join(['<option value="{0}"> {0}</option>'.format(b) for b in filters])
    
    return render_template('ldc.html', filters=filt_list)
    
# Load the LDC results page
@app_exoctk.route('/ldc_results', methods=['GET', 'POST'])
def exoctk_ldc_results():
        
    # Get the input from the form
    modeldir = request.form['modeldir']
    profiles = list(filter(None,[request.form.get(pf) for pf in ['uniform', 'linear', 'quadratic', 'square-root', 'logarithmic', 'exponential', '3-parameter', '4-parameter']]))
    bandpass = request.form['bandpass']
    teff = int(request.form['teff'])
    logg = float(request.form['logg'])
    feh = float(request.form['feh'])
    
    # Get models from local directory if necessary
    if not modeldir:
        modeldir = request.form['local_files']
    
    # Make the model grid, caching if necessary
    cached = cache.get(modeldir)
    if cached:
        model_grid = cached
        print('Fetching grid from cache:',modeldir)
    else:
        print('Not cached:',modeldir)
        model_grid = ExoCTK.core.ModelGrid(modeldir)
    
        if len(model_grid.data)>0:
            cache.set(modeldir, model_grid, timeout=300)
        
    if len(model_grid.data)==0:
        
        return render_template('ldc_error.html', teff=teff, logg=logg, feh=feh, \
                    band=bandpass or 'None', profile=profile, models=model_grid.path, \
                    script=script)
    
    # Trim the grid and calculate
    Teff_rng = find_closest(model_grid.Teff_vals,teff)
    logg_rng = find_closest(model_grid.logg_vals,logg)
    FeH_rng = find_closest(model_grid.FeH_vals,feh)
    model_grid.customize(Teff_rng=Teff_rng, logg_rng=logg_rng, FeH_rng=FeH_rng)
    
    # Calculate the coefficients for each profile
    coeff_list, mu_list, r_list, poly_list = [], [], [], []
    fig = plt.figure()
    for profile in profiles:
        coeff, mu, r = ExoCTK.ldc.ldcfit.ldc(teff, logg, feh, model_grid, profile, bandpass=bandpass, plot=fig)

        # LaTeX formatting
        r = '{:.4f} R_\odot'.format(r*1.4377979836321077e-11)
        mu = '{:.4f}'.format(mu)
        coeff_vals = ', '.join(['c{} = {:.3f}'.format(n+1,c) for n,c in enumerate(coeff)])
        poly = '\({}\)'.format(ExoCTK.ldc.ldcfit.ld_profile(profile, latex=True)).replace('*','\cdot').replace('\e','e')
        
        # Add the results to the lists
        coeff_list.append(coeff_vals)
        mu_list.append(mu)
        r_list.append(r)
        poly_list.append(poly)
        
    # Make the results into a list for easy printing
    table = at.Table([profiles, poly_list, coeff_list], names=['Profile','\(I(\mu)/I(\mu=1)\)','Coefficients'])
    table = '\n'.join(table.pformat(html=True)).replace('<table','<table class="results"')
    profile = ', '.join(profiles)

    # Make the matplotlib plot into a Bokeh plot
    p = mpl.to_bokeh()
    p.y_range = Range1d(0, 1)
    p.x_range = Range1d(0, 1)

    # # # Make LaTeX labels
    # # xlabel = LatexLabel(text="\mu", angle=0, angle_units='deg',
    # #                    x=250, y=-20, x_units='screen', y_units='screen',
    # #                    render_mode='css', text_font_size='16pt',
    # #                    background_fill_color='#ffffff')
    # # ylabel = LatexLabel(text="I(\mu)/I(\mu=1)", angle=90, angle_units='deg',
    # #                    x=-100, y=250, x_units='screen', y_units='screen',
    # #                    render_mode='css', text_font_size='16pt',
    # #                    background_fill_color='#ffffff')
    # # p.add_layout(xlabel)
    # # p.add_layout(ylabel)

    script, div = components(p)

    return render_template('ldc_results.html', teff=teff, logg=logg, feh=feh, \
                band=bandpass or 'None', profile=profile, mu=mu_list[0], \
                r=r_list[0], models=model_grid.path, table=table, script=script, plot=div)

# Load the LDC error page
@app_exoctk.route('/ldc_error', methods=['GET', 'POST'])
def exoctk_ldc_error():
    return render_template('ldc_error.html')

# Load the TOT page
@app_exoctk.route('/tot', methods=['GET', 'POST'])
def exoctk_tot():
    return render_template('tot.html')

# Load the TOR page
@app_exoctk.route('/tor', methods=['GET', 'POST'])
def exoctk_tor():
    return render_template('tor.html')

# Load the TOT results page
@app_exoctk.route('/tot_results', methods=['GET', 'POST'])
def exoctk_tot_results():
    
    exo_dict  = ExoCTK.tot.transit_obs.load_exo_dict()
    inst_dict = ExoCTK.tot.transit_obs.load_mode_dict('WFC3 G141')
    
    # Get the input from the form
    exo_dict['star']['hmag']      = float(request.form['hmag'])     # H-band magnitude of the system
    exo_dict['planet']['exopath'] = request.form['exopath']         # filename for model spectrum [wavelength, flux]
    exo_dict['planet']['w_unit']  = request.form['w_unit']          # wavelength unit (um or nm)
    exo_dict['planet']['f_unit']  = request.form['f_unit']          # flux unit (fp/f* or (rp/r*)^2)
    exo_dict['planet']['depth']   = 4.0e-3                          # Approximate transit/eclipse depth for plotting purposes
    exo_dict['planet']['i']       = float(request.form['i'])        # Orbital inclination in degrees
    exo_dict['planet']['ars']     = float(request.form['ars'])      # Semi-major axis in units of stellar radii (a/R*)
    exo_dict['planet']['period']  = float(request.form['period'])   # Orbital period in days
    
    # Detector and Observation inputs (Make these form inputs!)
    exo_dict['calculation'] = 'scale'
    inst_dict['configuration']['detector']['subarray']     = 'GRISM256'   # Subarray size, GRISM256 or GRISM512
    inst_dict['configuration']['detector']['nsamp']        = 10           # WFC3 NSAMP, 1..15
    inst_dict['configuration']['detector']['samp_seq']     = 'SPARS5'     # WFC3 SAMP-SEQ, SPARS5, SPARS10, or SPARS25
    exo_dict['observation']['transit_duration']            = 4170         # Full transit/eclipse duration in seconds
    exo_dict['observation']['norbits']                     = 4            # Number of HST orbits per visit
    exo_dict['observation']['noccultations']               = 5            # Number of transits/eclipses
    exo_dict['observation']['nchan']                       = 15           # Number of spectrophotometric channels
    exo_dict['observation']['scanDirection']               = 'Forward'    # Spatial scan direction, Forward or Round Trip
    exo_dict['observation']['schedulability']              = '30'         # % time HST can observe target (30 or 100)
    
    # Run PandExo
    deptherr, rms, ptsOrbit = ExoCTK.tot.transit_obs.run_pandexo(exo_dict, inst_dict, output_file='wasp43b_G141.p')
    
    # Plot the model spectrum with simpulated data and uncertainties
    specfile   = exo_dict['planet']['exopath']
    w_unit     = exo_dict['planet']['w_unit']
    grism      = inst_dict['configuration']['instrument']['disperser']
    nchan      = exo_dict['observation']['nchan']
    binspec = ExoCTK.tot.transit_obs.plot_PlanSpec(specfile, w_unit, grism, deptherr, nchan, smooth=10,
                     labels=['Model Spectrum', 'Simulated Obs.'])
    
    # Make the matplotlib plot into a Bokeh plot
    p = mpl.to_bokeh()
    plt.close()
    xmin, xmax = (1.125,1.650)
    ymin, ymax = (np.min(binspec)-2*deptherr, np.max(binspec)+2*deptherr)
    p.y_range = Range1d(ymin, ymax)
    p.x_range = Range1d(xmin, xmax)
    sim_script, sim_plot = components(p)
    
    # Plot the transit curves
    numorbits = exo_dict['observation']['norbits']
    depth     = exo_dict['planet']['depth']
    inc       = exo_dict['planet']['i']
    aRs       = exo_dict['planet']['ars']
    period    = exo_dict['planet']['period']
    windowSize= 20                                  # observation start window size in minutes
    minphase, maxphase = ExoCTK.tot.transit_obs.calc_StartWindow('eclipse', rms, ptsOrbit, numorbits, depth, inc, 
                                              aRs, period, windowSize, ecc=0, w=90.)
                                              
    # Make the matplotlib plot into a Bokeh plot
    p = mpl.to_bokeh()
    plt.close()
    obs_script, obs_plot = components(p)
    
    # Make some HTML for the input summary
    summary = """<h3>Input</h3>
    <table>
        <tr>
            <td>H-band magnitude of the system</td>
            <td>{}</td>
        </tr>
        <tr>
            <td>Orbital inclination [degrees]</td>
            <td>{}</td>
        </tr>
        <tr>
            <td>Semi-major axis [R*]</td>
            <td>{}</td>
        </tr>
        <tr>
            <td>Orbital period [days]</td>
            <td>{}</td>
        </tr>
    </table>
    
    <h3>Result</h3>
    <table>
        <tr>
            <td>Start observations between orbital phases:</td>
            <td>{:.4f} - {:.4f}</td>
        </tr>
    </table>
    """.format(exo_dict['star']['hmag'], exo_dict['planet']['i'], exo_dict['planet']['ars'], exo_dict['planet']['period'],
               minphase, maxphase)
    
    return render_template('tot_results.html', summary=summary, sim_script=sim_script, sim_plot=sim_plot, 
                           obs_script=obs_script, obs_plot=obs_plot)
    
def find_closest(A, a):
    """
    Find the two neighboring models for a given parameter
    
    Parameters
    ----------
    A: array-like
        The array to search
    a: float, int
        The value to search for
    
    Returns
    -------
    tuple
        The values to the left and right of 'a' in 'A'
    """
    idx = np.clip(A.searchsorted(a), 1, len(A)-1)
    return (A[idx-1], A[idx])

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
    __implementation__ = """
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
      
    if @model.orientation == 'v'
      angle = angle + 90

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

## -- RUN
if __name__ == '__main__':
    port = int(os.environ.get('PORT', 5000))
    app_exoctk.run(host='0.0.0.0', port=port, debug=True)
