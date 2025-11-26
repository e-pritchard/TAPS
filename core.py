# TAPS: Toolkit for Analysis of PRISM Spectra
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
import os
import sys
import contextlib
from scipy.interpolate import griddata

#code parameters
CODE_PATH = os.path.dirname(os.path.abspath(__file__))+'/../'
MODEL_FOLDER = os.path.join(CODE_PATH,'NIRSpec_PRISM_standards/')
MODEL_FOLDER_NIR = os.path.join(CODE_PATH,'NIR_standards/')

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
        #Opens the fits file & select index 1 (SPECID) where the data we wish to access lives

        if "fits" in self.file:
            data = fits.open(self.file)[1].data
            self.wave = data['wave'] * u.micron
            self.flux = data['flux'] * u.microjansky #fnu #NORMALIZE
            self.noise = data["err"] * u.microjansky #err_fnu 
            self.variance = self.noise**2
            
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

            if self.flux_label == r'$f_{\lambda}\$':
                self.flux = ((self.flux * const.c)/(self.wave**2)).to((10**-20)*u.erg*(u.cm**-2)*(u.s**-1)*(u.angstrom**-1))
                self.noise = ((self.noise * const.c)/(self.wave**2)).to((10**-20)*u.erg*(u.cm**-2)*(u.s**-1)*(u.angstrom**-1))
                self.variance = self.noise**2


        elif "nirspec" in self.file and "csv" in self.file:
            data = pd.read_csv(self.file)
            self.wave = data['wave'].values * u.micron
            self.flux = data['flux'].values * ((10**-20)*u.erg*(u.cm**-2)*(u.s**-1)*(u.angstrom**-1)) #flam NORMALIZE BY DIVIDING BY MAX VALUE
            self.noise = data['unc'].values * ((10**-20)*u.erg*(u.cm**-2)*(u.s**-1)*(u.angstrom**-1)) #err_flam   
            self.variance = self.noise**2

            if "/" in self.file:
                piecedir = self.file.split('/')
                pieceundscr = piecedir[-1].split('_')
                self.name = pieceundscr[1] + "_" + pieceundscr[2]
                self.name_err = "e_" + pieceundscr[2]
            else:
                pieceundscr = self.file.split('_')
                self.name = pieceundscr[1] + "_" + pieceundscr[2]
                self.name_err = "e_" + pieceundscr[2]

        elif "classify" in self.file and "csv" in self.file:
            data = pd.read_csv(self.file)
            self.wave = data['#WAVELENGTH'].values * u.micron
            self.flux = data['FLUX'].values * ((10**-20)*u.erg*(u.cm**-2)*(u.s**-1)*(u.angstrom**-1)) #flam NORMALIZE BY DIVIDING BY MAX VALUE
            self.noise = data['UNCERTAINTY'].values * ((10**-20)*u.erg*(u.cm**-2)*(u.s**-1)*(u.angstrom**-1)) #err_flam   
            self.variance = self.noise**2

            if "/" in self.file:
                piecedir = self.file.split('/')
                pieceundscr = piecedir[-1].split('_')
                self.name = pieceundscr[1]
                self.name_err = "e_" + pieceundscr[1]
            else:
                pieceundscr = self.file.split('_')
                self.name = pieceundscr[1]
                self.name_err = "e_" + pieceundscr[1]

        elif "txt" in self.file:
            df_read = pd.read_csv(self.file)
            wrong_header = df_read.columns
            header_values = wrong_header.values
            header_split = header_values[0].split("  ")
            header_split[0] = float(header_split[0])
            header_split[1] = float(header_split[1])
            
            wave_list = []
            flux_list = []
            wave_list.append(header_split[0])
            flux_list.append(header_split[1])
            
            for i in range(len(df_read)):
                value_split = df_read.iloc[i,0].split("  ")
                value_split[0] = float(value_split[0])
                value_split[1] = float(value_split[1])
                wave_list.append(value_split[0])
                flux_list.append(value_split[1])
            df_dict = {"wave" : wave_list , "flux" : flux_list}
            
            data = pd.DataFrame(data=df_dict)
            self.wave = data['wave'].values * u.AA
            self.wave = self.wave.to(u.micron)
            self.flux = data['flux'].values * u.microjansky #I DO NOT KNOW THE REAL UNITS OF FLUX FOR THE 4 SOURCES!
            #-----------------------------------------------------------------------------------------------------
            self.noise = data['flux'].values * u.microjansky #err_fnu
            #---------------------------------------------------------------------------
            #BE WARE BE WARE BE WARE BE WARE
            #currently do not have uncertainty for luhman sources! THIS IS JUST SO SOURCES CAN BE READ
            self.variance = self.noise**2

            if self.flux_label == r'$f_{\lambda}\$':
                self.flux = ((self.flux * const.c)/(self.wave**2)).to((10**-20)*u.erg*(u.cm**-2)*(u.s**-1)*(u.angstrom**-1))
                self.noise = ((self.noise * const.c)/(self.wave**2)).to((10**-20)*u.erg*(u.cm**-2)*(u.s**-1)*(u.angstrom**-1))
                self.variance = self.noise**2
            
            piecedot = self.file.split('.')
            self.name = piecedot[0]

    
    
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


#creating standardset array
standardset = [] 
for standfile in os.listdir(MODEL_FOLDER):
    standard = spec(MODEL_FOLDER + standfile)
    standardset.append(standard)

standardset_NIR = []
for standfile in os.listdir(MODEL_FOLDER_NIR):
    standard = spec(MODEL_FOLDER_NIR + standfile)
    standardset_NIR.append(standard)

def interpolate(spectrum, stan):
    standard = copy.deepcopy(stan)
    stanflxint = griddata(standard.wave, standard.flux, spectrum.wave, method = 'linear', rescale = True)
    standard.flux = np.array(stanflxint) * ((10**-20)*u.erg*(u.cm**-2)*(u.s**-1)*(u.angstrom**-1))
    standard.wave = spectrum.wave 
    return standard
      

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
    if type(spectrum) == spec or type(spectrum) == TAPS.spec:
        output = copy.deepcopy(spectrum)
        output.noise = spectrum.noise / np.nanmax(spectrum.flux)
        output.flux = spectrum.flux / np.nanmax(spectrum.flux)
    else: 
        print("Need a TAPS.spec object!")

    return output   


def alpha(spec1, spec2):
    alphanum = 0 
    alphadenom = 0
    for i in range(len(spec1.flux)):
        if np.isfinite(spec1.flux[i].value) and np.isfinite(spec2.flux[i].value):
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
        if np.isfinite(spec1.flux[i].value) and np.isfinite(spec2.flux[i].value):
            chi_squared += ((spec1.flux[i] - (alph * spec2.flux[i])) / (spec1.noise[i]))**2     
            
    return float(chi_squared)
            #print(chi_squared)

def reducedchisquare(spec1, chisquare):
    dof = len(spec1.wave) - 1
    redchisqr = chisquare/dof
    return redchisqr
#ADD REDUCED CHI SQUARE

def trim(spectrum, rng):
    spec = copy.deepcopy(spectrum)
    mask = np.zeros(len(spec.wave))
    mask[np.where(np.logical_and(spec.wave.value > rng[0],spec.wave.value < rng[1]))] = 1
    spec.wave = spec.wave[mask == 1]
    spec.flux = spec.flux[mask == 1]
    spec.noise = spec.noise[mask == 1]
    spec.variance = spec.variance[mask == 1]
    return spec

def classifystandard(spectrum):
    specnorm = normalizespec(spectrum)  
    standnormlist = []
    chisquares = []
    alphas = []
    standardsetint = []
    #standnames = []
    #stanflxint = [] #list to hold interpolated standard flux
    
    for standard in standardset:
        stanint = interpolate(spectrum, standard)
        standnorm = normalizespec(stanint)
        alph = alpha(specnorm, standnorm)
        chisqur = chisquare(specnorm, standnorm)
        standnormlist.append(standnorm)
        chisquares.append(chisqur)
        alphas.append(alph)
   
        
    chimin = np.min(chisquares)
    redchisqr = reducedchisquare(spectrum, chimin)
    minindex = np.argmin(chisquares)
    bestfit = standnormlist[minindex]
    alphmin = alphas[minindex]
    bestfitname = bestfit.name

    chisqr_formatted = ("{:.1f}".format(chimin))
    alpha_formatted = ("{:.1f}".format(alphmin))

    compspec(specnorm, bestfit, alphmin, redchisqr)                
    return f"$\chi^{2}$ = {chisqr_formatted}" , f"$\alpha$ = {alpha_formatted}" , "Best fit is " + bestfitname


def classifystandard_NIR(spectrum): 
    #This iteration trims the inputted spectrum according to the standard compared
    #Trimming allows for a more accurate reduced chi squared 
    #This iteration plots the wave range 0.5-2.5 microns
    specnorm = normalizespec(spectrum)  
    standnormlist = []
    chisquares = []
    alphas = []
    standardsetint = []
    specset_trimmed = []
    #standnames = []
    #stanflxint = [] #list to hold interpolated standard flux
    
    for standard in standardset_NIR:
        #print(f"Standard's wave range before interpolation {standard.wave}")
        #print(f"Spectrum's wave range before trimming {len(specnorm.wave)}")
        stan_rng = [np.nanmin(standard.wave.value) - 0.01, np.nanmax(standard.wave.value) + 0.01]
        specnorm_trimmed = trim(specnorm, stan_rng)
        #print(f"Spectrum's wave range after trimming {len(specnorm_trimmed.wave)}")
        stanint = interpolate(specnorm_trimmed, standard)
        #print(f"Standard's wave range after interpolation {stanint.wave}")
        standnorm = normalizespec(stanint)
        alph = alpha(specnorm_trimmed, standnorm)
        chisqur = chisquare(specnorm_trimmed, standnorm)
        standnormlist.append(standnorm)
        specset_trimmed.append(specnorm_trimmed)
        chisquares.append(chisqur)
        alphas.append(alph)
   
        
    chimin = np.min(chisquares)
    minindex = np.argmin(chisquares)
    bestfit = standnormlist[minindex]
    bestfit_spec = specset_trimmed[minindex]
    redchisqr = reducedchisquare(bestfit_spec, chimin)
    alphmin = alphas[minindex]
    bestfitname = bestfit.name

    chisqr_formatted = ("{:.1f}".format(chimin))
    alpha_formatted = ("{:.1f}".format(alphmin))

    compspec(bestfit_spec, bestfit, alphmin, redchisqr)                
    return f"$\chi^{2}$ = {chisqr_formatted}" , f"$\alpha$ = {alpha_formatted}" , "Best fit is " + bestfitname
    #ADD PLOTTING OPTION TO CLASSIFY BY STANDARD

#------------------------------------------------------------------------------------------------

# def classifystandard_NIR_og(spectrum): 
    #This iteration does not take into account the usable range of the NIR standards
    #This iteration's reduced chi squared is less accurate as it includes unused points
    #This iteration plots the wave range 0.5-5.0 microns
#     specnorm = normalizespec(spectrum)  
#     standnormlist = []
#     chisquares = []
#     alphas = []
#     standardsetint = []
#     #standnames = []
#     #stanflxint = [] #list to hold interpolated standard flux
    
#     for standard in standardset_NIR:
#         stanint = interpolate(specnorm, standard)
#         standnorm = normalizespec(stanint)
#         alph = alpha(specnorm, standnorm)
#         chisqur = chisquare(specnorm, standnorm)
#         standnormlist.append(standnorm)
#         chisquares.append(chisqur)
#         alphas.append(alph)
   
        
#     chimin = np.min(chisquares)
#     redchisqr = reducedchisquare(spectrum, chimin)
#     minindex = np.argmin(chisquares)
#     bestfit = standnormlist[minindex]
#     alphmin = alphas[minindex]
#     bestfitname = bestfit.name

#     chisqr_formatted = ("{:.1f}".format(chimin))
#     alpha_formatted = ("{:.1f}".format(alphmin))

#     compspec(specnorm, bestfit, alphmin, redchisqr)                
#     return f"$\chi^{2}$ = {chisqr_formatted}" , f"$\alpha$ = {alpha_formatted}" , "Best fit is " + bestfitname
#     #ADD PLOTTING OPTION TO CLASSIFY BY STANDARD

#---------------------------------------------------------------------------------------------------------------------

def compspec(spec1, spec2, alpha=1, redchisqr = 1, err=True):

    #This function graphs two different spectra onto the same plot and calculates chi for a spectrum and a standard model
    #spec1 is intended as a source while spec2 is intended for a standard model

    # chi_squared = chisquare(spec1, spec2)
    # alpha_formatted = ("{:.1f}".format(alpha))
    # chisqr_formatted = ("{:.1f}".format(chisqr))
    #chi = (chi_squared)**(1/2)
    
     #this simply plots the graphs inputed; as a visual for chi

    fig, axs = plt.subplots(2, 1, gridspec_kw={'height_ratios': [5, 1]}, figsize=(9, 5), sharex=True)  # 2, 1 for two rows one column

    axs[0].tick_params(axis='y', labelsize=12)
    axs[0].set_ylabel(r'$f_{\lambda}\ [10^{-20}ergs^{-1}cm^{-2}\AA^{-1}]$', fontsize = 13)
    axs[0].axhline(0, color ='k', linestyle ='--', lw=1)
    axs[0].plot(spec1.wave, spec1.flux, label = spec1.name, color = 'black', lw=1.5)
    axs[0].plot(spec2.wave, alpha * spec2.flux, label = spec2.name.split('_')[0] + ' standard', color = '#37CDFA', lw = 4, alpha = 0.7)
    
    if err:
        #axs[0].plot(spec1.wave, spec1.noise, label = spec1.name_err, color = 'grey', lw = 1.2)
        axs[0].fill_between(spec1.wave.value, spec1.noise.value, -1*spec1.noise.value,
                    color='k', alpha=0.3) 
        
    axs[0].plot([], [], ' ', label=r"$\alpha = %.1f, \ \chi_{r}^2 = %.1f$" % (alpha, redchisqr))
    # axs[0].plot([], [], ' ', label=fr"$\alpha$ = {alpha_formatted}")
    # axs[0].plot([], [], ' ', label=fr"$\chi^2$ = {chisqr_formatted}")
    axs[0].legend(fontsize = 13)


    axs[1].tick_params(axis='both', labelsize=12)
    axs[1].set_xlabel(r'$\lambda_{obs}\ [{\mu}m]$', fontsize = 13)
    axs[1].set_ylabel(r'$\Delta$', fontsize = 13)
    axs[1].axhline(0, color ='k', linestyle ='--', lw=1)
    difference = spec1.flux - (alpha * spec2.flux)
    axs[1].plot(spec1.wave, difference, color = 'k', lw = 1.5)

    plt.show()



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



#PROBLEM 2

def x_HI(T, P):
    Z_rot = 0
    I = 1.67*10**(-41)*u.cm**2 *u.g
    j = 0
    while j <= 100:
        if (j%2) == 0:
            Z_rot += (2*j + 1)*np.exp(-(j*(j+1))*const.hbar**2 / (I*const.k_B*T*u.K).decompose())
            j += 1
        elif (j%2) != 0:
            Z_rot += 3*(2*j + 1)*np.exp(-j*(j+1)*const.hbar**2 / (I*const.k_B*T*u.K).decompose())
            j += 1
    
    Z_elec_H2 = 1 + np.exp(-4.52*u.eV/(const.k_B*T*u.K))
    
    Z_elec_HI = 0
    E_ion = 13.6*u.eV
    n = 1
    n_max = ((const.k_B*T*u.K)**(1/3) / (2*((P*u.erg/u.cm**3)**(1/3))*const.a0))**(1/2)
    while n <= (n_max//1):
        Z_elec_HI += 4*(n)*np.exp((-E_ion*(1-(1/n**2))/ (const.k_B*T*u.K)))
        n += 1
    
    m = 2*const.m_p
    Z_tra = (2*np.pi*m*const.k_B*T*u.K / const.h**2)**(3/2) * (const.k_B * T*u.K) / (P*u.erg / u.cm**3)
    
    Z_vib = 0
    n = 0
    w_vib = 3.8*10**14 * u.Hz
    while n <= 100:
        Z_vib += np.exp(-(n+0.5)*const.hbar*w_vib / (const.k_B * T*u.K))
        n+=1
    
    C = (np.pi*const.m_p *const.k_B *T* u.K/ const.h**2)**(3/2) * const.k_B * T*u.K * Z_elec_HI**2 * np.exp(-4.52*u.eV/ (const.k_B*T*u.K)) / ((P*u.erg/u.cm**3)*Z_elec_H2*Z_vib*Z_rot)
    
    x = (-C + (C**2 + 4*C)**(1/2)) / 2
    
    return x

temp_array = np.logspace(2, 5, 100)
pressure_array = np.logspace(0, 13, 100)

x_matrix = np.full([len(temp_array), len(pressure_array)], np.nan)

for i, P in enumerate(pressure_array):
    for j, T in enumerate(temp_array):
        x_matrix[i,j] = x_HI(T,P)

x
