from flask import Flask, render_template, request, redirect, make_response

app_exoctk = Flask(__name__)

import ExoCTK
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
    
# Load the LDC results page
@app_exoctk.route('/ldc_results', methods=['GET', 'POST'])
def exoctk_ldc_results():
    
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

    return render_template('ldc_results.html', summary=summary, script=script, plot=div, version=VERSION)

# Load the TOT page
@app_exoctk.route('/tot', methods=['GET', 'POST'])
def exoctk_tot():
    return render_template('tot.html', version=VERSION)

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
                           obs_script=obs_script, obs_plot=obs_plot, version=VERSION)
                           
# tony darnell, scott lewis, ask to do one about Backyard Worlds?
# Write blog post

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