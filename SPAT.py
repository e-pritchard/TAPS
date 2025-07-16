from astropy.io import fits
import astropy.units as u
from astropy import constants as const
import matplotlib.pyplot as plt
import numpy as np
import splat

class spec(splat.Spectrum):
    def __init__(self, file):
        self.flam = None
        self.flam_err = None
        self.file = file
        self.variance = []
        self.flux_unit = u.microjansky

    def readfile(self):
        #THINKING OF UPDATING to specify flxtype in the readfile so that do not have to specify in multiple functions (also needs to be 
        #updated with julia's code that will reduce number of steps)
        #Opens the fits file & select index 1 (SPECID) where the data we wish to access lives
        data = fits.open(self.file)[1].data
    
        #Within data is 10 columns, we wish to access the columns titled, "wave", "flux", "err"
        self.wave = data["wave"] * u.micron
        self.flux = data["flux"] * u.microjansky #fnu
        self.noise = data["err"] * u.microjansky #err_fnu
        self.variance = self.noise**2

        self.flam = ((self.flux * const.c)/(self.wave**2)).to((10**-20)*u.erg*(u.cm**-2)*(u.s**-1)*(u.angstrom**-1))
        self.flam_err = ((self.noise * const.c)/(self.wave**2)).to((10**-20)*u.erg*(u.cm**-2)*(u.s**-1)*(u.angstrom**-1))

    def plot(self, flxtype, name = None):
    #plots the spectrum in either fnu or flam as specified
        plt.figure(figsize=(9,4))

        if flxtype == "flam":
            y = self.flam
            e_y = self.flam_err
            ylabel = r'$f_{\lambda}\ [10^{-20}ergs^{-1}cm^{-2}\AA^{-1}]$'
            
        elif flxtype == "fnu":
            y = self.flux
            e_y = self.noise
            ylabel = r'$f_{\nu}\ [{\mu}Jy]$'


        if name:
            plt.plot(self.wave, y, color='#5d0eff', lw=1.2, label=name)
            plt.legend()
        
        else:
            plt.plot(self.wave, y, color='#5d0eff', lw=1.2)
            
        plt.plot(self.wave, e_y, color = 'k', lw = 1.2) 
        plt.xticks(np.arange(0.5,5.5,step=0.5))
        plt.xlabel(r'$\lambda_{obs}\ [{\mu}m]$')
        plt.ylabel(ylabel)
        plt.grid()
        return plt.show()

    def normalize(self):
    #NEED TO UPDATE THIS (don't wan't it to actually update the spec object). may just write our own code here
    #this allows both flam and fnu to become normalized
        self.flam_err = self.flam_err / np.nanmax(self.flam)
        self.flam = self.flam / np.nanmax(self.flam)
    #splat's normalize applies to self.wave and self.flux (fnu) by default
        return super().normalize()


def classifystandard(spec):
    #classifies by comparison to standard
    return splat.classifyByStandard(spec)

#need to add classification by index & template
    

def compspec(spec1, spec2, err=True):

    #This function graphs two different spectra onto the same plot
    #spec1 is intended as a source while spec2 is intended for a standard model
    #FOR NOW this is a bit bare bones but a good start!
    #need to add chi
    
    plt.figure(figsize=(9,4))
    plt.xlabel(r'$\lambda_{obs}\ [{\mu}m]$')
    plt.ylabel(r'$f_{\lambda}\ [10^{-20}ergs^{-1}cm^{-2}\AA^{-1}]$')
    plt.plot(spec1.wave, spec1.flam, label="Borg")
    plt.plot(spec2.wave, spec2.flam, label="RUBIES")
    if err==True:
        plt.plot(spec1.wave, spec1.err_flam, label="err_borg")
    plt.legend(fontsize = "medium")
    plt.show()
