from flask import Flask, render_template, request
from flask import redirect, make_response, current_app
from flask_cache import Cache

app_exoctk = Flask(__name__)

import ExoCTK
import numpy as np
import matplotlib.pyplot as plt
import astropy.table as at
import astropy.units as q
import bokeh
import os
from svo_filters import svo
from bokeh import mpl
from bokeh.core.properties import Override
from bokeh.embed import components
from bokeh.models import Range1d
from bokeh.models import Label
from bokeh.models import HoverTool
from bokeh.models import ColumnDataSource
from bokeh.plotting import figure, show
from bokeh.models.widgets import Panel, Tabs

# define the cache config keys, remember that it can be done in a settings file
app_exoctk.config['CACHE_TYPE'] = 'null'

# register the cache instance and binds it on to your app
# cache = Cache(app_exoctk)

# Nice colors for plotting
COLORS = ['blue', 'red', 'green', 'orange', 
          'cyan', 'magenta', 'pink', 'purple']

# Supported profiles
PROFILES = ['uniform', 'linear', 'quadratic', 
            'square-root', 'logarithmic', 'exponential', 
            '3-parameter', '4-parameter']

# Redirect to the index
VERSION = ExoCTK.__version__
@app_exoctk.route('/')
@app_exoctk.route('/index')

# Load the Index page
def index():
    return render_template('index.html')

# Load the LDC page
@app_exoctk.route('/ldc', methods=['GET', 'POST'])
def exoctk_ldc():
    # Get all the available filters
    filters = svo.filters()['Band']
    
    # Make HTML for filters
    filt_list = '\n'.join(['<option value="{0}"{1}> {0}</option>'.format(b,\
                ' checked' if b=='Kepler.K' else '') for b in filters])
    
    return render_template('ldc.html', filters=filt_list)
    
# Load the LDC results page
@app_exoctk.route('/ldc_results', methods=['GET', 'POST'])
def exoctk_ldc_results():
        
    # Get the input from the form
    modeldir = request.form['modeldir']
    profiles = list(filter(None,[request.form.get(pf) for pf in PROFILES]))
    bandpass = request.form['bandpass']
    teff = int(request.form['teff'])
    logg = float(request.form['logg'])
    feh = float(request.form['feh'])
    mu_min = float(request.form['mu_min'])
    n_bins = request.form.get('n_bins')
    n_channels = request.form.get('n_channels')
    wl_min = request.form.get('wave_min')
    wl_max = request.form.get('wave_max')
    
    # Get models from local directory if necessary
    if not modeldir:
        modeldir = request.form['local_files']
    
    # # Make the model grid, caching if necessary
    # cached = cache.get(modeldir)
    # if cached:
    #     model_grid = cached
    #     print('Fetching grid from cache:',modeldir)
    # else:
    #     print('Not cached:',modeldir)
    #     model_grid = ExoCTK.core.ModelGrid(modeldir, resolution=500)
    #
    #     if len(model_grid.data)>0:
    #         cache.set(modeldir, model_grid, timeout=300)
    
    model_grid = ExoCTK.core.ModelGrid(modeldir, resolution=500)
    
    # No data, redirect to the error page
    if len(model_grid.data)==0:
        
        return render_template('ldc_error.html', teff=teff, logg=logg, feh=feh, \
                    band=bandpass or 'None', profile=profile, models=model_grid.path, \
                    script=script)
                    
    # Trim the grid to the correct wavelength
    # to speed up calculations, if a bandpass is given
    min_max = model_grid.wave_rng
    if bandpass in svo.filters()['Band'] or bandpass=='tophat':
        kwargs = {'n_bins':int(n_bins)} if n_bins else \
                 {'n_channels':int(n_channels)} if n_channels else {}
                 
        if wl_min and wl_max:
            kwargs['wl_min'] = float(wl_min)*q.um
            kwargs['wl_max'] = float(wl_max)*q.um
            
        bandpass = svo.Filter(bandpass, **kwargs)
        min_max = (bandpass.WavelengthMin,bandpass.WavelengthMax)
        n_bins = bandpass.n_bins
        
        # Get the filter plot
        TOOLS = 'box_zoom,resize,reset'
        bk_plot = figure(tools=TOOLS, title=bandpass.filterID, plot_width=400, plot_height=300,
                         x_range=Range1d(bandpass.WavelengthMin,bandpass.WavelengthMax))
        if n_bins>1:
            bk_plot.line(*bandpass.raw, line_width=5, color='black', alpha=0.1)
            for i,(x,y) in enumerate(bandpass.rsr):
                bk_plot.line(x, y, color=(COLORS*5)[i])
        else:
            bk_plot.line(*bandpass.rsr)
        bk_plot.xaxis.axis_label = 'Wavelength [um]'
        bk_plot.yaxis.axis_label = 'Throughput'
        filt_script, filt_plot = components(bk_plot)
        
        plt.close()
        
    # Trim the grid to nearby grid points to speed up calculation
    full_rng = [model_grid.Teff_vals,model_grid.logg_vals,model_grid.FeH_vals]
    trim_rng = ExoCTK.core.find_closest(full_rng, [teff,logg,feh], 
                                        n=1, values=True)
                                        
    model_grid.customize(Teff_rng=trim_rng[0], logg_rng=trim_rng[1], 
                         FeH_rng=trim_rng[2], wave_rng=min_max)
                         
    # Calculate the coefficients for each profile
    grid_point = ExoCTK.ldc.ldcfit.ldc(teff, logg, feh, model_grid, profiles, 
                    mu_min=mu_min, bandpass=bandpass, plot=False, colors=COLORS,
                    save='output.txt')
                    
    # Draw the figure
    tabs = []
    for i in range(n_bins or 1):
        
        # PLot it
        TOOLS = 'box_zoom,box_select,crosshair,resize,reset,hover'
        fig = figure(tools=TOOLS, x_range=Range1d(0, 1), y_range=Range1d(0, 1),
                     plot_width=800, plot_height=400)
        ld_funcs = [ExoCTK.ldc.ldcfit.ld_profile(p) for p in profiles]
        ExoCTK.ldc.ldcplot.ld_plot(ld_funcs, grid_point, fig=fig, bin_idx=i)
                                    
        # Plot formatting
        fig.legend.location = 'bottom_right'
        fig.xaxis.axis_label = 'mu'
        fig.yaxis.axis_label = 'Intensity'
        
        tabs.append(Panel(child=fig, title=str(grid_point['centers'][0][i])))
                
    final = Tabs(tabs=tabs)
    
    # Get HTML
    script, div = components(final)
    
    # Read the file into a string and delete it
    with open('output.txt', 'r') as f:
        file_as_string = f.read()
    os.remove('output.txt')
    
    # # Format mu and r_eff vals
    # r_eff = '{:.4f} R_\odot'.format(grid_point['r_eff']*1.438e-11)
    # mu_eff = '{:.4f}'.format(0)
    r_eff = mu_eff = ''
    
    # Make a table for each profile with a row for each wavelength bin
    profile_tables = []
    for profile in profiles:
        
        # Make LaTeX for polynomials
        latex = ExoCTK.ldc.ldcfit.ld_profile(profile, latex=True)
        poly = '\({}\)'.format(latex).replace('*','\cdot').replace('\e','e')
        
        # Make the table into LaTeX
        table = grid_point[profile]['coeffs']
        table.rename_column('wavelength', '\(\lambda_\mbox{eff}\hspace{5px}(\mu m)\)')
        
        # Add the results to the lists
        html_table = '\n'.join(table.pformat(max_width=500, html=True))\
                     .replace('<table','<table class="results"')
        
        # Add the table title
        header = '<strong>{}</strong><br><p>\(I(\mu)/I(\mu=1)\) = {}</p>'.format(profile,poly)
        html_table = header+html_table
        
        profile_tables.append(html_table)
        
    return render_template('ldc_results.html', teff=teff, logg=logg, feh=feh, \
                band=bandpass.filterID or '-', mu=mu_eff, profile=', '.join(profiles), \
                r=r_eff, models=model_grid.path, table=profile_tables, \
                script=script, plot=div, file_as_string=repr(file_as_string), \
                filt_plot=filt_plot, filt_script=filt_script)

# Load the LDC error page
@app_exoctk.route('/ldc_error', methods=['GET', 'POST'])
def exoctk_ldc_error():
    return render_template('ldc_error.html')

# Load the TOT page
@app_exoctk.route('/tot', methods=['GET', 'POST'])
def exoctk_tot():
    return render_template('tot.html')

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

# Save the results to file
@app_exoctk.route('/download', methods=['POST'])
def exoctk_savefile():
    file_as_string = eval(request.form['file_as_string'])
    
    response = make_response(file_as_string)
    response.headers["Content-type"] = 'text; charset=utf-8'
    response.headers["Content-Disposition"] = "attachment; filename=ExoXTK_results.txt"
    return response

# @app_exoctk.route('/export', methods=['POST'])
# def exoctk_export():
#
#     # Get all the checked rows
#     checked = request.form
#
#     # Get column names
#     results = [list(eval(checked.get('cols')))]
#
#     for k in sorted(checked):
#         if k.isdigit():
#
#             # Convert string to list and strip HTML
#             vals = eval(checked[k])
#             for i,v in enumerate(vals):
#                 try:
#                     vals[i] = str(v).split('>')[1].split('<')[0]
#                 except:
#                     pass
#
#             results.append(vals)
#
#     # Make an array to export
#     results = np.array(results, dtype=str)
#     filename = 'ONCdb_results.txt'
#     np.savetxt(filename, results, delimiter='|', fmt='%s')
#
#     with open(filename, 'r') as f:
#         file_as_string = f.read()
#     os.remove(filename)  # Delete the file after it's read
#
#     response = make_response(str(file_as_string))
#     response.headers["Content-type"] = 'text; charset=utf-8'
#     response.headers["Content-Disposition"] = "attachment; filename={}".format(filename)
#     return response

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