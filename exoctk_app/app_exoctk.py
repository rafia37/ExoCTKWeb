import flask
from flask import Flask, render_template, request, redirect, make_response, current_app
from flask_cache import Cache

app_exoctk = Flask(__name__)

import shutil
import os
from datetime import datetime

import ExoCTK
from ExoCTK.pal import exotransmit
import numpy as np
import matplotlib.pyplot as plt
import astropy.table as at
import astropy.units as q
from astropy.extern.six.moves import StringIO
import bokeh
import os
from bokeh import mpl
from bokeh.resources import INLINE
from bokeh.util.string import encode_utf8
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

exotransmit_dir = os.environ.get('EXOTRANSMIT_DIR')
modelgrid_dir = os.environ.get('MODELGRID_DIR')

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
    
    # Get models from local directory if necessary
    if modeldir=='default':
        modeldir = modelgrid_dir
    # elif not modeldir:
    #     modeldir = request.form['local_files']
    
    try:
        teff = int(request.form['teff'])
        logg = float(request.form['logg'])
        feh = float(request.form['feh'])
        mu_min = float(request.form['mu_min'])
    except:
        message = 'Could not calculate limb darkening with the above input parameters.'
        
        return render_template('ldc_error.html', teff=request.form['teff'], logg=request.form['logg'], feh=request.form['feh'], \
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
    if bandpass in ExoCTK.svo.filters()['Band'] or bandpass=='tophat':
        
        try:
            
            kwargs = {'n_bins':int(n_bins)} if n_bins else \
                     {'pixels_per_bin':int(pixels_per_bin)} if pixels_per_bin else {}
                 
            if wl_min and wl_max:
                kwargs['wl_min'] = float(wl_min)*q.um
                kwargs['wl_max'] = float(wl_max)*q.um
            
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
                     .replace('<table','<table class="results"')
        
        # Add the table title
        header = '<strong>{}</strong><br><p>\(I(\mu)/I(\mu=1)\) = {}</p>'.format(profile,poly)
        html_table = header+html_table
        
        profile_tables.append(html_table)
        
    return render_template('ldc_results.html', teff=teff, logg=logg, feh=feh, \
                band=bp_name, mu=mu_eff, profile=', '.join(profiles), \
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

@app_exoctk.route('/exotransmit', methods=['GET','POST'])
def exotransmit_page():
    """ 
    Run Exo-Transmit taking inputs from the HTML form and plot the results
    """
    if exotransmit_dir is None:
        flask.abort(404)
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
        'exotransmit.html',
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


@app_exoctk.route('/exotransmit_result', methods=['POST'])
def save_exotransmit_result():
    table_string = flask.request.form['data_file']
    return flask.Response(table_string, mimetype="text/dat",
                          headers={"Content-disposition":
                                      "attachment; filename=exotransmit.dat"})
