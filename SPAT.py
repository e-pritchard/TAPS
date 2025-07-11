from astropy.io import fits
import astropy.units as u
from astropy import constants as const
import matplotlib.pyplot as plt
import numpy as np
import splat

class spec(splat.Spectrum):
    def __init__(self, file):
        self.fnu = None
        self.fnu_err = None
        self.flam = None
        self.flam_err = None
        self.file = file

    def readfile(self):
        #Opens the fits file & select index 1 (SPECID) where the data we wish to access lives
        data = fits.open(self.file)[1].data
    
        #Within data is 10 columns, we wish to access the columns titled, "wave", "flux", "err"
        self.wave = data["wave"] * u.micron
        self.fnu = data["flux"] * u.microjansky
        self.err_fnu = data["err"] * u.microjansky

        self.flam = ((self.fnu * const.c)/(self.wave**2)).to((10**-20)*u.erg*(u.cm**-2)*(u.s**-1)*(u.angstrom**-1))
        self.err_flam = ((self.err_fnu * const.c)/(self.wave**2)).to((10**-20)*u.erg*(u.cm**-2)*(u.s**-1)*(u.angstrom**-1))

    def plot(self, flxtype):
        plt.figure(figsize=(9,4))

        if flxtype == "flam":
            y = self.flam
            e_y = self.err_flam
            ylabel = r'$f_{\lambda}\ [10^{-20}ergs^{-1}cm^{-2}\AA^{-1}]$'
            
        elif flxtype == "fnu":
            y = self.fnu
            e_y = self.err_fnu
            ylabel = r'$f_{\nu}\ [{\mu}Jy]$'

        #else:
            #raise ExceptionType("Need a fluxtype")
            
        plt.plot(self.wave, y, color = '#5d0eff', lw = 1.2)
        plt.plot(self.wave, e_y, color = 'k', lw = 1.2) 
        plt.xticks(np.arange(0.5,5.5,step=0.5))
        plt.xlabel(r'$\lambda_{obs}\ [{\mu}m]$')
        plt.ylabel(ylabel)
        plt.grid()
        return plt.show()
