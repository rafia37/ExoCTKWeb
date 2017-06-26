## -- IMPORTS
import glob
from datetime import datetime
import os
import shutil

from astropy.extern.six.moves import StringIO
import astropy.table as at
import matplotlib.pyplot as plt
import numpy as np

import bokeh
from bokeh import mpl
from bokeh.core.properties import Override
from bokeh.embed import components
from bokeh.models import ColumnDataSource
from bokeh.models import HoverTool
from bokeh.models import Label
from bokeh.models import Range1d
from bokeh.mpl import to_bokeh
from bokeh.plotting import figure
from bokeh.plotting import output_file
from bokeh.plotting import show
from bokeh.plotting import save
from bokeh.resources import INLINE
from bokeh.util.string import encode_utf8
import flask
from flask import current_app
from flask import Flask
from flask import make_response
from flask import redirect
from flask import render_template
from flask import request
from flask_cache import Cache

import ExoCTK
from ExoCTK.pal import exotransmit
#from ExoCTK.tor.tor import create_tor_dict

from tor import create_tor_dict

#from ExoCTK.tor import create_tor_dict

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
        temp = float(request.form['temp'])
        sat_max = float(request.form['sat_lvl'])
        sat_mode = request.form['sat']
        throughput = request.form['tp']
        filt = request.form['filt']
        ins = request.form['ins']
        subarray = request.form['subarray']
        n_reset = int(request.form['n_reset'])
    
        tor_err = 0
        # Specific error catching
        if n_group.isdigit():
            if n_group <= 1:
                tor_err = 'Please try again with at least one group.'
        else:
            if n_group != 'optimize':
                tor_err = "You need to double check your group input. Please put the number of groups per integration or 'optimize' and we can calculate it for you."
        if (mag > 12) or (mag < 5.5):
            tor_err = 'Looks like we do not have useful approximations for your magnitude. Could you give us a number between 5.5-12 in a different band?'
        if obs_time <= 0:
            tor_err = 'You have a negative transit time -- I doubt that will be of much use to anyone.'
        if temp <= 0:
            tor_err = 'You have a negative temperature -- DOES THIS LOOK LIKE A MASER TO YOU?'
        if sat_max <= 0:
            tor_err = 'You put in an underwhelming saturation level. There is something to be said for being too careful...'
        if (sat_mode == 'well') and sat_max > 1:
            tor_err = 'You are saturating past the full well. Is that a good idea?'
        if n_reset < 1:
            tor_err = 'You have no or negative resets. That is not allowed!'

        if type(tor_err) == str:
            return render_template('tor_error.html', tor_err=tor_err)
        
        # Only create the dict if the form input looks okay.
        tor_output = create_tor_dict(obs_time, n_group, mag, band, temp, sat_max, sat_mode, throughput, filt, ins, subarray, n_reset)
        if type(tor_output) == dict:
            tor_dict = tor_output
            return render_template('tor_results.html', tor_dict=tor_dict)
        else:
            tor_err = tor_output
            return render_template('tor_error.html', tor_err=tor_err)
    
    except (ValueError, TypeError) as e:
        tor_err = 'One of you numbers is NOT a number! Please try again!'
        return render_template('tor_error.html', tor_err=tor_err)

# Load the TOR background
@app_exoctk.route('/tor_background')
def exoctk_tor_background():
    return render_template('tor_background.html')

# Load filter profiles pages
@app_exoctk.route('/filter_profile_<ins>')
def exoctk_filter_profile(ins):

    filt_imgs = glob.glob('static/filter_dat/' + ins + '/' + ins + '*')
    names = [filt_img[19+len(ins):-4] for filt_img in filt_imgs]
    print(filt_imgs)
    filt_imgs = ['../' + filt_img for filt_img in filt_imgs]

    return render_template('tor_filter_profile.html', names=names, filt_imgs=filt_imgs, ins=ins)
    

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

exotransmit_dir = os.environ.get('EXOTRANSMIT_DIR')

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
def exotransmit():
    """ 
    Run Exo-Transmit taking inputs from the HTML form and plot the results
    """
    if exotransmit_dir is None:
        return render_template('tor_error.html', tor_err="There seems to be no directory in place for exo-transmit...")
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
    """
    
    return render_template(
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
    """
    return encode_utf8(html)


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
