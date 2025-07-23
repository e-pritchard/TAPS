from astropy.io import fits
import astropy.units as u
from astropy import constants as const
import matplotlib.pyplot as plt
import numpy as np
import splat

class spec(splat.Spectrum):
    def __init__(self, file, read_file: bool = True):
        self.flam = None
        self.flam_err = None
        self.file = file
        self.variance = []
        self.flux_unit = u.microjansky
        self.id = ''
        self.e_id = ''

        if read_file:
            self.readfile()

    def readfile(self):
        #THINKING OF UPDATING to specify flxtype in the readfile so that do not have to specify in multiple functions (also needs to be 
        #updated with julia's code that will reduce number of steps)
        #or have one default with a function to convert
        #Opens the fits file & select index 1 (SPECID) where the data we wish to access lives
        data = fits.open(self.file)[1].data
    
        #Within data is 10 columns, we wish to access the columns titled, "wave", "flux", "err"
        self.wave = data["wave"] * u.micron
        self.flux = data["flux"] * u.microjansky #fnu
        self.noise = data["err"] * u.microjansky #err_fnu
        self.variance = self.noise**2

        self.flam = ((self.flux * const.c)/(self.wave**2)).to((10**-20)*u.erg*(u.cm**-2)*(u.s**-1)*(u.angstrom**-1))
        self.flam_err = ((self.noise * const.c)/(self.wave**2)).to((10**-20)*u.erg*(u.cm**-2)*(u.s**-1)*(u.angstrom**-1))

        #
        pieceper = self.file.split('.')
        pieceundscr = pieceper[0].split('_')
        self.id = pieceundscr[0] + "_" + pieceundscr[2] + "_" + pieceundscr[3]
        self.e_id = "e_" + pieceundscr[3]


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

    #def normalize(self):
    #NEED TO UPDATE THIS (don't wan't it to actually update the spec object). may just write our own code here
    #this allows both flam and fnu to become normalized
        #self.flam_err = self.flam_err / np.nanmax(self.flam)
        #self.flam = self.flam / np.nanmax(self.flam)
    #splat's normalize applies to self.wave and self.flux (fnu) by default
        #return super().normalize()

def normalizespec(spectrum):
    if type(spectrum) == spec:
        output = spec(spectrum.file)
        output.flam_err = spectrum.flam_err / np.nanmax(spectrum.flam)
        output.flam = spectrum.flam / np.nanmax(spectrum.flam)
        output.noise = spectrum.noise / np.nanmax(spectrum.flux)
        output.flux = spectrum.flux / np.nanmax(spectrum.flux)
    else: 
        print("Need a SPAT.spec object!")

    return output 

def classifystandard(spec):
    #classifies by comparison to standard
    #think this only works for fnu because splat.classifyByStandard refers to self.wave & self.noise
    return splat.classifyByStandard(spec)

#need to add classification by index & template


def chisquare(spec1, spec2):
    #this is a bare bones chi
    #was modeled after sara's standard chi squared statistic
    
    chi_squared = 0
    alphanum = 0 
    alphadenom = 0
    for i in range(len(spec1.flam)):
        if not np.isnan(spec1.flam[i].value) or not np.isnan(spec2.flam[i].value):
            alphanum += ((spec1.flam[i])*(spec2.flam[i])) / (spec1.flam_err[i]**2)
            alphadenom += (spec2.flam[i]**2) / (spec1.flam_err[i]**2)
            alpha = alphanum / alphadenom
            chi_squared += ((spec1.flam[i] - (alpha * spec2.flam[i])) / (spec1.flam_err[i]))**2
    return chi_squared
    #print(chi_squared)

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
    plt.plot(spec1.wave, spec1.flam, label= spec1.id)
    plt.plot(spec2.wave, spec2.flam, label= spec2.id)
    plt.plot([],[], label = f"$\chi^{2}$ = {chisqr_formatted}", alpha = 0)
    if err==True:
        plt.plot(spec1.wave, spec1.flam_err, label= spec1.e_id)
    plt.legend(fontsize = "medium")
    plt.show()

    return chi_squared
