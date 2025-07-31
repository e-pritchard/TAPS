from astropy.io import fits
import astropy.units as u
from astropy import constants as const
import matplotlib.pyplot as plt
import numpy as np
import splat
import copy
import pandas as pd
import ucdmcmc
from astropy.utils.data import download_file
import copy
import io
import sys
import contextlib
%matplotlib inline

#need to add general code path to standards folder for others to use our code (see ucdmcmc code parameters section)

class spec(splat.Spectrum):
    def __init__(self, file, read_file: bool = True, flxtype = "flam"):
        self.file = file
        self.variance = []
        self.name_err = ""
        self.history = []
        self.instrument = "JWST-NIRSPEC-PRISM"

        if flxtype == "flam":
            self.flux_label = r'$f_{\lambda}\$'
            self.flux_unit = (10**-20)*u.erg*(u.cm**-2)*(u.s**-1)*(u.angstrom**-1)
        elif flxtype == "fnu":
            self.flux_label = r'$f_{\nu}\$'
            self.flux_unit = u.microjansky
            
        if read_file:
            self.readfile()


    def readfile(self):
        #THINKING OF UPDATING to specify flxtype in the readfile so that do not have to specify in multiple functions (also needs to be 
        #updated with julia's code that will reduce number of steps)
        #or have one default with a function to convert
        #Opens the fits file & select index 1 (SPECID) where the data we wish to access lives

        if "fits" in self.file:
            data = fits.open(self.file)[1].data
            self.wave = data['wave'] * u.micron
            self.flux = data['flux'] * u.microjansky #fnu
            self.noise = data["err"] * u.microjansky #err_fnu 
            
            if "/" in self.file:
                piecedir = self.file.split('/')
                pieceper = piecedir[-1].split('.')
                pieceundscr = pieceper[0].split('_')
                self.name = pieceundscr[0] + "_" + pieceundscr[2] + "_" + pieceundscr[3]
                self.name_err = "e_" + pieceundscr[3]
            else:
                pieceper = self.file.split('.')
                pieceundscr = pieceper[0].split('_')
                self.name = pieceundscr[0] + "_" + pieceundscr[2] + "_" + pieceundscr[3]
                self.name_err = "e_" + pieceundscr[3]


        elif "csv" in self.file:
            data = pd.read_csv(self.file)
            self.wave = data['wave'].values * u.micron
            self.flux = data['flux'].values * u.microjansky #fnu
            self.noise = data['unc'].values * u.microjansky #err_fnu

            if "/" in self.file:
                piecedir = self.file.split('/')
                pieceper = piecedir[-1].split('.')
                pieceundscr = pieceper[0].split('_')
                self.name = pieceundscr[1] + "_" + pieceundscr[2]
                self.name_err = "e_" + pieceundscr[2]
            else:
                pieceper = self.file.split('.')
                pieceundscr = pieceper[0].split('_')
                self.name = pieceundscr[1] + "_" + pieceundscr[2]
                self.name_err = "e_" + pieceundscr[2]
    
        if self.flux_label == r'$f_{\lambda}\$':
            self.flux = ((self.flux * const.c)/(self.wave**2)).to((10**-20)*u.erg*(u.cm**-2)*(u.s**-1)*(u.angstrom**-1))
            self.noise = ((self.noise * const.c)/(self.wave**2)).to((10**-20)*u.erg*(u.cm**-2)*(u.s**-1)*(u.angstrom**-1))

        self.variance = self.noise**2
    
    def plot(self):
    #plots the spectrum in either fnu or flam as specified
        if self.flux_label == r'$f_{\nu}\$':
            ylabel = r'$f_{\nu}\ [{\mu}Jy]$'
        elif self.flux_label == r'$f_{\lambda}\$':
            ylabel = r'$f_{\lambda}\ [10^{-20}ergs^{-1}cm^{-2}\AA^{-1}]$'
        y = self.flux
        e_y = self.noise

        plt.figure(figsize=(9,4))
        plt.plot(self.wave, y, color='#5d0eff', lw=1.2, label=self.name)
        plt.plot(self.wave, e_y, color = 'k', lw = 1.2, label=self.name_err) 
        plt.xticks(np.arange(0.5,5.5,step=0.5))
        plt.xlabel(r'$\lambda_{obs}\ [{\mu}m]$')
        plt.ylabel(ylabel)
        plt.grid()
        plt.legend()
        return plt.show()

#Everything commented out below are obselete functions and no longer in use (but are saved incase the code is useful

#        self.flam = ((self.flux * const.c)/(self.wave**2)).to((10**-20)*u.erg*(u.cm**-2)*(u.s**-1)*(u.angstrom**-1))
#        self.flam_err = ((self.noise * const.c)/(self.wave**2)).to((10**-20)*u.erg*(u.cm**-2)*(u.s**-1)*(u.angstrom**-1))

    #def normalize(self):
    #NEED TO UPDATE THIS (don't wan't it to actually update the spec object). may just write our own code here
    #this allows both flam and fnu to become normalized
        #self.flam_err = self.flam_err / np.nanmax(self.flam)
        #self.flam = self.flam / np.nanmax(self.flam)
    #splat's normalize applies to self.wave and self.flux (fnu) by default
        #return super().normalize()

#Below is a template for a convertflux method, which would update / change your current instance of spec if used:
#This method does not create a new instance, it only updates your current (Use fnutoflam and flamtofnu for new instances) 
        
#    def convertflux(spectrum):
#        if spectrum.flux_label == r"$f_{\nu}\$":
#            if not spectrum.flux_unit == u.microjansky:
#                spectrum.flux_unit = u.microjansky
#                spectrum.flux = spectrum.flux.to(u.microjansky)
#                spectrum.noise = spectrum.noise.to(u.microjansky)
#            spectrum.flux = ((spectrum.flux * const.c)/(spectrum.wave**2)).to((10**-20)*u.erg*(u.cm**-2)*(u.s**-1)*(u.angstrom**-1))
#            spectrum.noise = ((spectrum.noise * const.c)/(spectrum.wave**2)).to((10**-20)*u.erg*(u.cm**-2)*(u.s**-1)*(u.angstrom**-1))
#            spectrum.flux_unit = (10**-20)*u.erg*(u.cm**-2)*(u.s**-1)*(u.angstrom**-1)
#            spectrum.flux_label = r'$f_{\lambda}\$'
#            print(f"Flux has been converted to Flam with units {self.flux_unit}")
#    
#        elif spectrum.flux_label == r'$f_{\lambda}\$':
#            if not spectrum.flux_unit == (10**-20)*u.erg*(u.cm**-2)*(u.s**-1)*(u.angstrom**-1):
#                spectrum.flux_unit = (10**-20)*u.erg*(u.cm**-2)*(u.s**-1)*(u.angstrom**-1)
#                spectrum.flux = spectrum.flux.to((10**-20)*u.erg*(u.cm**-2)*(u.s**-1)*(u.angstrom**-1))
#            spectrum.flux = ((spectrum.wave**2)*(spectrum.flux))/(const.c).to(u.microjansky)
#            spectrum.noise = (((spectrum.wave**2)*(spectrum.noise))/(const.c)).to(u.microjansky)
#            spectrum.flux_unit = u.microjansky
#            spectrum.flux_label = r"$f_{\nu}\$"
#            print(f"Flux has been converted to Fnu with units {self.flux_unit}")
#        else:
#            print("something has gone wrong")

def fnutoflam(spectrum):
    newspec = copy.deepcopy(spectrum)
    if newspec.flux_label == r"$f_{\nu}\$":
        if not newspec.flux_unit == u.microjansky:
            newspec.flux_unit = u.microjansky
            newspec.flux = newspec.flux.to(u.microjansky)
            newspec.noise = newspec.noise.to(u.microjansky)
        newspec.flux = ((newspec.flux * const.c)/(newspec.wave**2)).to((10**-20)*u.erg*(u.cm**-2)*(u.s**-1)*(u.angstrom**-1))
        newspec.noise = ((newspec.noise * const.c)/(newspec.wave**2)).to((10**-20)*u.erg*(u.cm**-2)*(u.s**-1)*(u.angstrom**-1))
        newspec.flux_unit = (10**-20)*u.erg*(u.cm**-2)*(u.s**-1)*(u.angstrom**-1)
        newspec.flux_label = r'$f_{\lambda}\$'
        print(f"Flux have been converted to Flam with units {newspec.flux_unit}")
        return newspec
    else:
        print("Something has gone wrong")

def flamtofnu(spectrum):
    newspec = copy.deepcopy(spectrum)
    if newspec.flux_label == r'$f_{\lambda}\$':
        if not newspec.flux_unit == (10**-20)*u.erg*(u.cm**-2)*(u.s**-1)*(u.angstrom**-1):
            newspec.flux_unit = (10**-20)*u.erg*(u.cm**-2)*(u.s**-1)*(u.angstrom**-1)
            newspec.flux = newspec.flux.to((10**-20)*u.erg*(u.cm**-2)*(u.s**-1)*(u.angstrom**-1))
        newspec.flux = (((newspec.wave**2)*(newspec.flux))/(const.c)).to(u.microjansky)
        newspec.noise = (((newspec.wave**2)*(newspec.noise))/(const.c)).to(u.microjansky)
        newspec.flux_unit = u.microjansky
        newspec.flux_label = r"$f_{\nu}\$"
        print(f"Flux have been converted to Fnu with units {newspec.flux_unit}")
        return newspec
    else:
        print("Something has gone wrong")

def normalizespec(spectrum):
    if type(spectrum) == spec:
        output = spec(spectrum.file)
        output.noise = spectrum.noise / np.nanmax(spectrum.flux)
        output.flux = spectrum.flux / np.nanmax(spectrum.flux)
    else: 
        print("Need a SPAT.spec object!")

    return output   


def alpha(spec1, spec2):
    alphanum = 0 
    alphadenom = 0
    for i in range(len(spec1.flux)):
        if not np.isnan(spec1.flux[i].value) or not np.isnan(spec2.flux[i].value):
            alphanum += ((spec1.flux[i])*(spec2.flux[i])) / (spec1.noise[i]**2)
            alphadenom += (spec2.flux[i]**2) / (spec1.noise[i]**2)
            
    alpha = alphanum / alphadenom
    return alpha
            


def chisquare(spec1, spec2):
    #this is a bare bones chi
    #was modeled after sara's standard chi squared statistic

    chi_squared = 0
    alph = alpha(spec1, spec2)
    
    for i in range(len(spec1.flux)):
        if not np.isnan(spec1.flux[i].value) or not np.isnan(spec2.flux[i].value):
            chi_squared += ((spec1.flux[i] - (alph * spec2.flux[i])) / (spec1.noise[i]))**2     
            
    return float(chi_squared)
            #print(chi_squared)

flx = []
flx.append(griddata(standard1.wave.value, standard1.flux.value, spec1.wave.value, method='linear',rescale=True))

def classifystandard(specfile):
    spectrum = spec(specfile)
    standardset = '/Users/marylin/Desktop/UCSD/STARTastro/SPURS/NIRSpec_PRISM_standards/' #this needs to be general path directory
    chisquares = []
    alphas = []
    
    for standfile in os.listdir(standardset):
        standard = spec(standfile)
        alph = alpha(spectrum, standard)
        chisqur = chisquare(spectrum, standard)
    
        chisquares.append(chisqur)
        alphas.append(alph)
    
        
    chimin = np.min(chisquares)
    alphmin = np.min(alphas)
    
    minindex = np.argmin(chisquares)
    bestfit = standardset[minindex]

    chisqr_formatted = ("{:.1f}".format(chimin))
    alpha_formatted = ("{:.1f}".format(alphmin))
                       
    return f"$\chi^{2}$ = {chisqr_formatted}"
    return f"$\alpha$ = {alpha_formatted}"
    return "Best fit is" + bestfit


def compspec(spec1, spec2, err=True):

    #This function graphs two different spectra onto the same plot and calculates chi for a spectrum and a standard model
    #spec1 is intended as a source while spec2 is intended for a standard model

    chi_squared = chisquare(spec1, spec2)
    chisqr_formatted = ("{:.1f}".format(chi_squared))
    #chi = (chi_squared)**(1/2)
    
    #this simply plots the graphs inputed; as a visual for chi
    plt.figure(figsize=(9,4))
    plt.xlabel(r'$\lambda_{obs}\ [{\mu}m]$')
    plt.ylabel(r'$f_{\lambda}\ [10^{-20}ergs^{-1}cm^{-2}\AA^{-1}]$')
    plt.plot(spec1.wave, spec1.flux, label= spec1.id)
    plt.plot(spec2.wave, spec2.flux, label= spec2.id)
    if err==True:
        plt.plot(spec1.wave, spec1.noise, label= spec1.e_id)
    plt.plot([],[], label = f"$\chi^{2}$ = {chisqr_formatted}", alpha = 0)
    plt.legend(fontsize = "medium")
    plt.show()

    return chi_squared

def suppress_empty_figures(func, *args, **kwargs):
    """Runs a function and closes only truly empty figures (0 axes)."""
    before = set(plt.get_fignums())
    f = io.StringIO()
    with contextlib.redirect_stdout(f), contextlib.redirect_stderr(f):
        result = func(*args, **kwargs)
    after = set(plt.get_fignums())
    new_figs = after - before
    
    for fig_num in new_figs:
        fig = plt.figure(fig_num)
        if len(fig.axes) == 0:
            plt.close(fig_num)
    # If the fig axes is 0, it will not be plotted
    return result

def fit_models_to_sources(source, model_name):
    model, wave = ucdmcmc.getModelSet(model_name, 'JWST-NIRSPEC-PRISM')
    sp = spec(source)
    spsm = ucdmcmc.resample(sp, wave)

    ipar = suppress_empty_figures(ucdmcmc.fitGrid, spsm, model, file_prefix=sp.name + model_name + "Gridfit", output="allvalues", report = True)
    par = suppress_empty_figures(ucdmcmc.fitMCMC, spsm, model, p0=ipar, nstep=2500, file_prefix = sp.name + model_name + "_mcmc", output="all", verbose=False, report = True)

    plt.show() 

    return sp.name, model_name
