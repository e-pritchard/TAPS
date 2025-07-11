from astropy.io import fits
import astropy.units as u
from astropy import constants as const
import matplotlib.pyplot as plt
import numpy as np
import splat

class spectrum(splat.Spectrum):
    def __init__(self):
        self.fnu = None
        self.fnu_err = None
        self.flam = None
        self.flam_err = None

    def read_file(self, file_name):
        #Opens the fits file by:
        hdul = fits.open(file_name)
    
        #The data we wish to access is in index 1 (SPEC1D), so we can select index to it by:
        data = hdul[1].data
    
        #Within data is 10 columns, we wish to access the columns titled, "wave", "flux", "err"
        self.wave = data["wave"] * u.micron
        self.fnu = data["flux"] * u.microjansky
        self.fnu_err = data["err"] * u.microjansky
        
        lam_flux = self.fnu*(const.c) / self.wave**2
        lam_err = self.fnu_err*(const.c) / self.wave**2

        self.flam = lam_flux.to((1*10**(-20))*u.erg*(1/u.s)*(1/u.cm**2)*(1/u.AA))
        self.flam_err = lam_err.to((1*10**(-20))*u.erg*(1/u.s)*(1/u.cm**2)*(1/u.AA))

    def plot(self, fluxtype):
        plt.figure(figsize=(9,4))
        if fluxtype == "flam":
            plt.plot(self.wave, self.flam)
            plt.plot(self.wave, self.flam_err)
            plt.xlabel("Wavelength [${\mu}m$]")
            plt.ylabel("$f_{\lambda} [10^-20 ergs^-1 cm^-2 {\AA}^-1]$")
        elif fluxtype == "fnu":
            plt.plot(self.wave, self.fnu)
            plt.plot(self.wave, self.fnu_err)
            plt.xlabel("Wavelength [${\mu}m$]")
            plt.ylabel(r"$f_{\nu} [{\mu}Jy$]")
        else: 
            raise ExceptionType("Need a fluxtype")

        return plt.show()
