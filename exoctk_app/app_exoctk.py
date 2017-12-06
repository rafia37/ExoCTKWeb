## -- IMPORTS
import glob
from datetime import datetime
import os
import shutil
from astropy.extern.six.moves import StringIO
import astropy.table as at
import astropy.units as q
import matplotlib.pyplot as plt
import numpy as np
import os
import pandas as pd
from sqlalchemy import *
import astropy.units as u 
import astropy.constants as constants

import bokeh
from bokeh import mpl
from bokeh.resources import INLINE
from bokeh.util.string import encode_utf8
from bokeh.core.properties import Override
from bokeh.embed import components
from bokeh.models import ColumnDataSource
from bokeh.models import HoverTool
from bokeh.models import Label
from bokeh.models import Range1d
from bokeh.models.widgets import Panel
from bokeh.models.widgets import Tabs
from bokeh.mpl import to_bokeh
from bokeh.plotting import figure
from bokeh.plotting import output_file
from bokeh.plotting import show
from bokeh.plotting import save

import flask
from flask import current_app
from flask import Flask
from flask import make_response
from flask import redirect
from flask import render_template
from flask import request
#from flask_cache import Cache

import ExoCTK
from ExoCTK.pal import exotransmit
from ExoCTK.tor.tor import create_tor_dict
from ExoCTK.tor.contam_tool.resolve import *
from ExoCTK.tor.contam_tool.visibilityPA import *

## -- FLASK SET UP (?)
app_exoctk = Flask(__name__)

# define the cache config keys, remember that it can be done in a settings file
app_exoctk.config['CACHE_TYPE'] = 'null'

exotransmit_dir = os.environ.get('EXOTRANSMIT_DIR')
modelgrid_dir = os.environ.get('MODELGRID_DIR')

#fortney grid
fortgrid_dir = os.environ.get('FORTGRID_DIR')

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
    filters = ExoCTK.svo.filters()['Band']
    
    # Make HTML for filters
    filt_list = '\n'.join(['<option value="{0}"{1}> {0}</option>'.format(b,\
                ' selected' if b=='Kepler.K' else '') for b in filters])
    
    return render_template('ldc.html', filters=filt_list)
    
# Load the LDC results page
@app_exoctk.route('/ldc_results', methods=['GET', 'POST'])
def exoctk_ldc_results():
        
    # Get the input from the form
    modeldir = request.form['modeldir']
    profiles = list(filter(None,[request.form.get(pf) for pf in PROFILES]))
    bandpass = request.form['bandpass']

    # protect against injection attempts
    bandpass = bandpass.replace('<', '&lt')
    profiles = [str(p).replace('<', '&lt') for p in profiles]

    # Get models from local directory if necessary
    if modeldir=='default':
        modeldir = modelgrid_dir
    
    try:
        teff = int(request.form['teff'])
        logg = float(request.form['logg'])
        feh = float(request.form['feh'])
        mu_min = float(request.form['mu_min'])
    except:
        teff = str(request.form['teff']).replace('<', '&lt')
        logg = str(request.form['logg']).replace('<', '&lt')
        feh = str(request.form['feh']).replace('<', '&lt')
        message = 'Could not calculate limb darkening with the above input parameters.'
        
        return render_template('ldc_error.html', teff=teff, logg=logg, feh=feh, \
                    band=bandpass or 'None', profile=', '.join(profiles), models=modeldir, \
                    message=message)
                    
    n_bins = request.form.get('n_bins')
    pixels_per_bin = request.form.get('pixels_per_bin')
    wl_min = request.form.get('wave_min')
    wl_max = request.form.get('wave_max')
    
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
    if not hasattr(model_grid, 'data'):
        message = 'Could not find a model grid to load. Please check your path.'
    
        return render_template('ldc_error.html', teff=teff, logg=logg, feh=feh, \
                    band=bandpass or 'None', profile=', '.join(profiles), models=model_grid.path, \
                    message=message)
        
    else:
        
        if len(model_grid.data)==0:
        
            message = 'Could not calculate limb darkening with the above input parameters.'
        
            return render_template('ldc_error.html', teff=teff, logg=logg, feh=feh, \
                        band=bandpass or 'None', profile=', '.join(profiles), models=model_grid.path, \
                        message=message)
                    
    # Trim the grid to the correct wavelength
    # to speed up calculations, if a bandpass is given
    min_max = model_grid.wave_rng
    if bandpass in ExoCTK.svo.filters()['Band'] or bandpass in ['tophat','NIRISS.GR700XD.1']:
        
        try:
            
            kwargs = {'n_bins':int(n_bins)} if n_bins else \
                     {'pixels_per_bin':int(pixels_per_bin)} if pixels_per_bin else {}
                 
            if wl_min and wl_max:
                kwargs['wl_min'] = float(wl_min)*q.um
                kwargs['wl_max'] = float(wl_max)*q.um
            
            # Manually create GR700XD filter
            if bandpass=='NIRISS.GR700XD.1':
                p = os.path.join(os.path.dirname(ExoCTK.__file__),'data/filters/NIRISS.GR700XD.1.txt')
                bandpass = ExoCTK.svo.Filter(bandpass, filter_directory=np.genfromtxt(p, unpack=True), **kwargs)
            else:
                bandpass = ExoCTK.svo.Filter(bandpass, **kwargs)
                
            min_max = (bandpass.WavelengthMin,bandpass.WavelengthMax)
            n_bins = bandpass.n_bins
            bp_name = bandpass.filterID
            
            # Get the filter plot
            TOOLS = 'box_zoom,resize,reset'
            bk_plot = figure(tools=TOOLS, title=bp_name, plot_width=400, plot_height=300,
                             x_range=Range1d(bandpass.WavelengthMin,bandpass.WavelengthMax))
                         
            bk_plot.line(*bandpass.raw, line_width=5, color='black', alpha=0.1)
            try:
                for i,(x,y) in enumerate(bandpass.rsr):
                    bk_plot.line(x, y, color=(COLORS*5)[i])
            except:
                bk_plot.line(*bandpass.raw)
            
            bk_plot.xaxis.axis_label = 'Wavelength [um]'
            bk_plot.yaxis.axis_label = 'Throughput'
            
            js_resources = INLINE.render_js()
            css_resources = INLINE.render_css()
            filt_script, filt_plot = components(bk_plot)
        
            plt.close()
        except:
            message = 'Insufficient filter information. Please complete the form and try again!'
        
            return render_template('ldc_error.html', teff=teff, logg=logg, feh=feh, \
                        band=bandpass or 'None', profile=', '.join(profiles), models=model_grid.path, \
                        message=message)
    else:
        bp_name = bandpass or '-'
        filt_plot = filt_script = ''
        
    # Trim the grid to nearby grid points to speed up calculation
    full_rng = [model_grid.Teff_vals,model_grid.logg_vals,model_grid.FeH_vals]
    trim_rng = ExoCTK.core.find_closest(full_rng, [teff,logg,feh], 
                                        n=1, values=True)
                                        
    if not trim_rng:
        
        message = 'Insufficient models grid points to calculate coefficients.'
        
        return render_template('ldc_error.html', teff=teff, logg=logg, feh=feh, \
                band=bp_name, profile=', '.join(profiles), models=model_grid.path,\
                message=message)
    
    elif not profiles:
        
        message = 'No limb darkening profiles have been selected. Please select at least one.'
        
        return render_template('ldc_error.html', teff=teff, logg=logg, feh=feh, \
                band=bp_name, profile=', '.join(profiles), models=model_grid.path,\
                message=message)
    
    else:
        
        try:
            model_grid.customize(Teff_rng=trim_rng[0], logg_rng=trim_rng[1], 
                         FeH_rng=trim_rng[2], wave_rng=min_max)
                         
        except:
            
            message = 'Insufficient wavelength coverage to calculate coefficients.'
        
            return render_template('ldc_error.html', teff=teff, logg=logg, feh=feh, \
                    band=bp_name, profile=', '.join(profiles), models=model_grid.path,\
                    message=message)
                         
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
                     .replace('<table','<table id="myTable" class="table table-striped table-hover"')
        
        # Add the table title
        header = '<br></br><strong>{}</strong><br><p>\(I(\mu)/I(\mu=1)\) = {}</p>'.format(profile,poly)
        html_table = header+html_table

        profile_tables.append(html_table)
        
    return render_template('ldc_results.html', teff=teff, logg=logg, feh=feh, \
                band=bp_name, mu=mu_eff, profile=', '.join(profiles), \
                r=r_eff, models=model_grid.path, table=profile_tables, \
                script=script, plot=div, file_as_string=repr(file_as_string), \
                filt_plot=filt_plot, filt_script=filt_script, js=js_resources, css=css_resources)

# Load the LDC error page
@app_exoctk.route('/ldc_error', methods=['GET', 'POST'])
def exoctk_ldc_error():
    return render_template('ldc_error.html')

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
    plt.savefig('static/plots/sim.png')
    sim_plot = mpl.to_bokeh(plt.gcf())
    output_file('test_sim.html')
    save(sim_plot)
    plt.close()
    xmin, xmax = (1.125,1.650)
    ymin, ymax = (np.min(binspec)-2*deptherr, np.max(binspec)+2*deptherr)
    sim_plot.y_range = Range1d(ymin, ymax)
    sim_plot.x_range = Range1d(xmin, 6)
#    sim_plot.y_axis_label = 'Wavelength (micron)'
    sim_script, sim_plot = components(sim_plot)
    
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
    plt.savefig('static/plots/obs.png')
    obs_plot = mpl.to_bokeh(plt.gcf())
    output_file('test_obs.html')
    save(obs_plot)
    plt.close()
    obs_script, obs_plot = components(obs_plot)
    
    exo_dict['minphase'] = round(minphase, 4)
    exo_dict['maxphase'] = round(maxphase, 4)

    sim = open('test_sim.html')
    lines = sim.readlines()
    sim_html = "".join(lines)
    sim.close()

    obs = open('test_obs.html')
    lines = obs.readlines()
    obs_html = "".join(lines)
    print(obs_html)
    html_dict = {'sim': sim_html, 'obs': obs_html}
    
    return render_template('tot_results.html', exo_dict=exo_dict, html_dict=html_dict)

# Load the WIP page
@app_exoctk.route('/wip')
def exoctk_wip():
    return render_template('wip.html')

# Load the TOR page
@app_exoctk.route('/tor', methods=['GET', 'POST'])
def exoctk_tor():

    ins = glob.glob('static/filter_dat/*')
    button_ins = [this_ins[18:] for this_ins in ins]
    
    return render_template('tor.html', button_ins=button_ins)

# Load the TOR results
@app_exoctk.route('/tor_results', methods=['GET', 'POST'])
def exoctk_tor_results():
    
    try:
        n_group = request.form['groups']
        mag = float(request.form['mag'])
        band = request.form['band']
        obs_time = float(request.form['T'])
#        temp = float(request.form['temp'])
        sat_max = float(request.form['sat_lvl'])
        sat_mode = request.form['sat']
#        throughput = request.form['tp']
        filt = request.form['filt']
        ins = request.form['ins']
        subarray = request.form['subarray']
        n_reset = int(request.form['n_reset'])
        infile = os.environ['tor_pandeia_path'] 

        tor_err = 0
        # Specific error catching
        if n_group.isdigit():
            n_group = int(n_group)
            if n_group <= 1:
               tor_err = 'Please try again with at least one group.'
        else:
            if n_group != 'optimize':
             tor_err = "You need to double check your group input. Please put the number of groups per integration or 'optimize' and we can calculate it for you."
        if (mag > 12) or (mag < 5.5):
            tor_err = 'Looks like we do not have useful approximations for your magnitude. Could you give us a number between 5.5-12 in a different band?'
        if obs_time <= 0:
            tor_err = 'You have a negative transit time -- I doubt that will be of much use to anyone.'
#    if temp <= 0:
#        tor_err = 'You have a negative temperature -- DOES THIS LOOK LIKE A MASER TO YOU?'
        if sat_max <= 0:
            tor_err = 'You put in an underwhelming saturation level. There is something to be said for being too careful...'
        if (sat_mode == 'well') and sat_max > 1:
            tor_err = 'You are saturating past the full well. Is that a good idea?'
        if n_reset < 1:
            tor_err = 'You have no or negative resets. That is not allowed!'

        if type(tor_err) == str:
            return render_template('tor_error.html', tor_err=tor_err)
    
    # Only create the dict if the form input looks okay.
        tor_output = create_tor_dict(float(obs_time), n_group, float(mag), str(band), str(filt), str(ins), str(subarray), str(sat_mode), float(sat_max), int(n_reset), infile)
        if type(tor_output) == dict:
            tor_dict = tor_output
            if tor_dict['n_group'] == 1:
                one_group_error = 'Be careful! This only predicts one group, and you may be in danger of oversaturating!'
            else:
                one_group_error = ""
            return render_template('tor_results.html', tor_dict=tor_dict, one_group_error=one_group_error)
        else:
            tor_err = tor_output
            return render_template('tor_error.html', tor_err=tor_err)
    
    except (IOError, KeyError) as e:
        if e == IOError:
            tor_err = 'One of you numbers is NOT a number! Please try again!'
        else:
            tor_err = 'Looks like you have mismatched your instrument/filter/subarray. Please try again.'
        return render_template('tor_error.html', tor_err=tor_err)
    except Exception as e:
        tor_err = 'This is not an error we anticipated, but the error caught was : ' + str(e)
        return render_template('tor_error.html', tor_err=tor_err)

# Load the TOR background
@app_exoctk.route('/tor_background')
def exoctk_tor_background():
    return render_template('tor_background.html')

# Load the Tor2 page
@app_exoctk.route('/tor2', methods = ['GET', 'POST'])
def exoctk_tor2():

    contamVars = {}
    if request.method == 'POST':
        tname = request.form['targetname']
        contamVars['tname'] = tname
        contamVars['ra'], contamVars['dec'] = request.form['ra'], request.form['dec']
        contamVars['binComp'] = request.form['bininfo']
        contamVars['PAmax'] = request.form['pamax']
        contamVars['PAmin'] = request.form['pamin']
        print(contamVars)

        if request.form['submit'] == 'Resolve Target':
            contamVars['ra'], contamVars['dec'] = resolve_target(tname)
    
        if request.form['submit'] == 'Calculate contamination':
            contamVars['contam'] = True
            png = 'results/visibilityPA-'+tname+'.png'
            contamVars['kelt-8'] = 'images/contamination-KELT-8_PA0-360.png'
            contamVars['tyc-55'] = 'images/contamination-TYC_5530-1795-1_PA0-360.png'
            contamVars['wasp-62'] = 'images/contamination-wasp-62_PA0-360.png'
            return render_template('tor2_results.html', contamVars = contamVars, png = png)
    
        if request.form['submit'] == 'Calculate visibility':
            contamVars['visPA'] = True
#             dir = 'static/results'
            png = calc_vis(contamVars['ra'], contamVars['dec'], tname)
            return render_template('tor2_results.html', contamVars = contamVars, png = png)

    return render_template('tor2.html', contamVars = contamVars)


# Load filter profiles pages
@app_exoctk.route('/filter_profile_<ins>')
def exoctk_filter_profile(ins):

    filt_imgs = glob.glob('static/filter_dat/' + ins + '/' + ins + '*')
    names = [filt_img[19+len(ins):-4] for filt_img in filt_imgs]
    print(filt_imgs)
    filt_imgs = ['../' + filt_img for filt_img in filt_imgs]

    return render_template('tor_filter_profile.html', names=names, filt_imgs=filt_imgs, ins=ins)
    
# Save the results to file
@app_exoctk.route('/download', methods=['POST'])
def exoctk_savefile():
    file_as_string = eval(request.form['file_as_string'])
    
    response = make_response(file_as_string)
    response.headers["Content-type"] = 'text; charset=utf-8'
    response.headers["Content-Disposition"] = "attachment; filename=ExoXTK_results.txt"
    return response

def exotransmit_run(eos, tp, g, R_p, R_s, P, Rayleigh):
    current_dir = os.path.abspath(os.curdir)
    now = datetime.now().isoformat()
    os.chdir(os.path.join(exotransmit_dir, 'runs'))
    os.mkdir(now)
    os.chdir(now)
    output_file = os.path.relpath('result.dat', start=exotransmit_dir)
    exotransmit.exotransmit(base_dir=exotransmit_dir,
                            EOS_file=os.path.join('/EOS', eos),
                            T_P_file=os.path.join('/T_P', tp),
                            g=g,
                            R_planet=R_p*69.911e6,
                            R_star=R_s*695.7e6,
                            P_cloud=P*100000,
                            Rayleigh=Rayleigh,
                            output_file='/'+output_file,
                            )
    x, y = np.loadtxt('result.dat', skiprows=2, unpack=True)
    os.chdir(current_dir)
    shutil.rmtree(os.path.join(exotransmit_dir, 'runs', now))
    return x, y

@app_exoctk.context_processor
def utility_processor():
    def process_eos(fname):
        params = fname[:-4].split('_')[1:]
        if len(params) == 1:
            params[0] = params[0] + ' Only Atmosphere'
            return params[0]
        else:
            if 'solar' in params[0]:
                params[0] = params[0].replace('p', '.').replace('X',
                                                                'x ').replace(
                    'solar', 'Solar Metallicity')
            if 'CtoO' in params:
                params[1] = params[1].replace('p', '.') + ' C:O'
                params.pop(2)
            if 'gas' in params:
                params[params.index('gas')] = 'Gas only'
            if 'cond' in params:
                if 'graphite' in params:
                    params[params.index('cond')] = 'w/ Condensation and rainout'
                    params[params.index('graphite')] = 'including graphite'
                else:
                    params[params.index(
                        'cond')] = 'w/ Condensation and rainout (no graphite)'

        return ', '.join(params)

    def process_tp(fname):
        return fname[4:-4] + ' Isothermal'

    return dict(process_eos=process_eos, process_tp=process_tp)

def _param_validation(args):
    invalid = {}
    eos = args.get('eos', 'eos_0p1Xsolar_cond.dat')
    try:
        str(eos)
        if eos not in os.listdir(os.path.join(exotransmit_dir,'EOS')):
            invalid['eos'] = "Invalid chemistry template file"
            eos = 'eos_0p1Xsolar_cond.dat'
    except ValueError:
        invalid['eos'] = "Chemistry template file must be a string"
        eos = 'eos_0p1Xsolar_cond.dat'

    tp = args.get('tp', 't_p_800K.dat')
    try:
        str(tp)
        if eos not in os.listdir(os.path.join(exotransmit_dir,'EOS')):
            invalid['tp'] = "Invalid TP file"
            tp = 't_p_800K.dat'
    except ValueError:
        invalid['tp'] = "TP file must be a string"
        tp = 't_p_800K.dat'

    try:
        g = float(args.get('g', 9.8))
        if g <= 0:
            invalid['g'] = "Surface gravity must be greater than zero"
            g = 9.8
    except ValueError:
        invalid['g'] = "Surface gravity must be a float"
        g = 9.8

    try:
        R_p = float(args.get('R_p', 0.0915))
        if R_p <= 0:
            invalid['R_p'] = "Planet Radius must be greater than zero"
            R_p = 0.0915
    except ValueError:
        invalid['R_p'] = "Planet Radius must be a float"
        R_p = 0.0915

    try:
        R_s = float(args.get('R_s', 1.0))
        if R_s <= 0:
            invalid['R_s'] = "Star Radius must be greater than zero"
            R_s = 1.0
    except ValueError:
        invalid['R_s'] = "Star Radius must be a float"
        R_s = 1.0

    try:
        P = float(args.get('P', 0.0))
        if P < 0:
            invalid['P'] = "Cloud Pressure must be greater than or equal to zero"
            P = 0.0
    except ValueError:
        invalid['P'] = "Cloud Pressure must be a float"
        P = 0.0

    try:
        Rayleigh = float(args.get('Rayleigh', 1.0))
        if Rayleigh <= 0:
            invalid['Rayleigh'] = "Rayleight augmentation must be greater than zero"
            Rayleigh = 1.0
    except ValueError:
        invalid['Rayleigh'] = "Rayleigh augmentation must be a float"
        Rayleigh = 1.0

    return eos, tp, g, R_p, R_s, P, Rayleigh, invalid

@app_exoctk.route('/exotransmit_portal', methods=['GET','POST'])
def exotransmit_portal():
    """
        Run Exo-Transmit taking inputs from the HTML form and plot the results
        """
    if exotransmit_dir is None:
        return render_template('tor_error.html', tor_err="There seems to be no directory in place for exo-transmit...")
    
    # Grab the inputs arguments from the URL
    args = flask.request.args
    # Get all the form arguments in the url with defaults
    eos, tp, g, R_p, R_s, P, Rayleigh, invalid = _param_validation(args)
    if invalid:
        return flask.render_template('exotransmit_validation.html', invalid=invalid)
    if args:
        x, y = exotransmit_run(eos, tp, g, R_p, R_s, P, Rayleigh)
    else:
        x, y = np.loadtxt(os.path.join(exotransmit_dir, 'Spectra/default.dat'), skiprows=2, unpack=True)
    
    tab = at.Table(data=[x/1e-6, y/100])
    fh = StringIO()
    tab.write(fh, format='ascii.no_header')
    table_string = fh.getvalue()

    fig = figure(plot_width=1100, plot_height=400, responsive=False)
    fig.line(x/1e-6, y/100, color='Black', line_width=0.5)
    fig.xaxis.axis_label = 'Wavelength (um)'
    fig.yaxis.axis_label = 'Transit Depth'
    
    js_resources = INLINE.render_js()
    css_resources = INLINE.render_css()
    
    script, div = components(fig)
    
    html = flask.render_template(
                                 'exotransmit_portal.html',
                                 plot_script=script,
                                 plot_div=div,
                                 js_resources=js_resources,
                                 css_resources=css_resources,
                                 eos_files=os.listdir(os.path.join(exotransmit_dir,'EOS')),
                                 tp_files=os.listdir(os.path.join(exotransmit_dir, 'T_P')),
                                 tp=tp,
                                 eos=eos,
                                 g=g,
                                 R_p=R_p,
                                 R_s=R_s,
                                 P=P,
                                 Rayleigh=Rayleigh,
                                 table_string=table_string
                                 )
    return encode_utf8(html)

def _param_fort_validation(args):
    temp = args.get('ptemp',1000)
    chem = args.get('pchem','noTiO')
    cloud = args.get('cloud','0')
    pmass = args.get('pmass','1.5')
    m_unit = args.get('m_unit','M_jup')
    reference_radius = args.get('refrad',1)
    r_unit = args.get('r_unit','R_jup')
    rstar = args.get('rstar',1)
    rstar_unit = args.get('rstar_unit','R_sun')
    return temp,chem,cloud,pmass,m_unit,reference_radius,r_unit,rstar,rstar_unit

@app_exoctk.route('/fortney_portal', methods=['GET','POST'])
def fortney_portal():
    """
    Pull up Forntey Grid plot the results and download
    """
    
    # Grab the inputs arguments from the URL
    args = flask.request.args
    
    temp,chem,cloud,pmass,m_unit,reference_radius,r_unit,rstar,rstar_unit = _param_fort_validation(args)
    #get sqlite database

    try:
        db = create_engine('sqlite:///'+os.environ.get('FORTGRID_DIR'))
        header= pd.read_sql_table('header',db)
    except:
        raise Exception('Fortney Grid File Path is incorrect, or not initialized')

    if args:
        rstar=float(rstar)
        rstar = (rstar*u.Unit(rstar_unit)).to(u.km)
        reference_radius = float(reference_radius)
        rplan = (reference_radius*u.Unit(r_unit)).to(u.km)

        #clouds 
        if cloud.find('flat') != -1: 
            flat = int(cloud[4:])
            ray = 0 
        elif cloud.find('ray') != -1: 
            ray = int(cloud[3:])
            flat = 0 
        elif int(cloud) == 0: 
            flat = 0 
            ray = 0     
        else:
            flat = 0 
            ray = 0 
            print('No cloud parameter not specified, default no clouds added')
        
        #chemistry 
        if chem == 'noTiO': 
            noTiO = True
        if chem == 'eqchem': 
            noTiO = False
            #grid does not allow clouds for cases with TiO
            flat = 0 
            ray = 0 

        fort_grav = 25.0*u.m/u.s/u.s

        temp = float(temp)
        df = header.loc[(header.gravity==fort_grav) & (header.temp==temp)
                           & (header.noTiO==noTiO) & (header.ray==ray) &
                           (header.flat==flat)]

        wave_planet=np.array(pd.read_sql_table(df['name'].values[0],db)['wavelength'])[::-1]
        r_lambda=np.array(pd.read_sql_table(df['name'].values[0],db)['radius'])*u.km
        z_lambda = r_lambda- (1.25*u.R_jup).to(u.km) #all fortney models have fixed 1.25 radii
        
        #scale with planetary mass 
        pmass=float(pmass)
        mass = (pmass*u.Unit(m_unit)).to(u.kg)
        gravity = constants.G*(mass)/(rplan.to(u.m))**2.0 #convert radius to m for gravity units
        #scale lambbda (this technically ignores the fact that scaleheight is altitude dependent)
        #therefore, it will not be valide for very very low gravities
        z_lambda = z_lambda*fort_grav/gravity
        
        #create new wavelength dependent R based on scaled ravity
        r_lambda = z_lambda + rplan
        #finally compute (rp/r*)^2
        flux_planet = np.array(r_lambda**2/rstar**2)

        x=wave_planet
        y=flux_planet
    else:   
        df= pd.read_sql_table('t1000g25_noTiO',db)
        x, y = df['wavelength'], df['radius']**2.0/7e5**2.0
    
    tab = at.Table(data=[x, y])
    fh = StringIO()
    tab.write(fh, format='ascii.no_header')
    table_string = fh.getvalue()

    fig = figure(plot_width=1100, plot_height=400, responsive=False)
    fig.line(x, 1e6*(y-np.mean(y)), color='Black', line_width=0.5)
    fig.xaxis.axis_label = 'Wavelength (um)'
    fig.yaxis.axis_label = 'Rel. Transit Depth (ppm)'
    
    js_resources = INLINE.render_js()
    css_resources = INLINE.render_css()
    
    script, div = components(fig)
    
    html = flask.render_template(
                                 'fortney_portal.html',
                                 plot_script=script,
                                 plot_div=div,
                                 js_resources=js_resources,
                                 css_resources=css_resources,
                                 temp=list(map(str, header.temp.unique())),
                                 table_string=table_string
                                 )
    return encode_utf8(html)

@app_exoctk.route('/fortney_result', methods=['POST'])
def save_fortney_result():
    table_string = flask.request.form['data_file']
    return flask.Response(table_string, mimetype="text/dat",
                          headers={"Content-disposition":
                          "attachment; filename=fortney.dat"})

@app_exoctk.route('/exotransmit_result', methods=['POST'])
def save_exotransmit_result():
    table_string = flask.request.form['data_file']
    return flask.Response(table_string, mimetype="text/dat",
                          headers={"Content-disposition":
                                      "attachment; filename=exotransmit.dat"})
                                      
## -- RUN
if __name__ == '__main__':
    port = int(os.environ.get('PORT', 5000))
    app_exoctk.run(host='0.0.0.0', port=port, debug=True)
