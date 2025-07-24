from astropy.io import fits
import astropy.units as u
from astropy import constants as const
import matplotlib.pyplot as plt
import numpy as np
import splat

class spec(splat.Spectrum):
    def __init__(self, file, read_file: bool = True):
        self.file = file
        self.variance = []
        #self.flux_unit = u.microjansky
        self.flux_unit = ''
        #self.flux_label = r'$f_{\nu}\$'
        self.flux_label = ''
        self.id = ''
        self.e_id = ""
        self.name = ""
        self.history = []

        if read_file:
            self.readfile()

    #def readfile(self):
        #THINKING OF UPDATING to specify flxtype in the readfile so that do not have to specify in multiple functions (also needs to be 
        #updated with julia's code that will reduce number of steps)
        #or have one default with a function to convert
        #Opens the fits file & select index 1 (SPECID) where the data we wish to access lives
        #data = fits.open(self.file)[1].data
    
        #Within data is 10 columns, we wish to access the columns titled, "wave", "flux", "err"
        #self.wave = data["wave"] * u.micron
        #self.flux = data["flux"] * u.microjansky #fnu
        #self.noise = data["err"] * u.microjansky #err_fnu
        #self.variance = self.noise**2

        #self.flam = ((self.flux * const.c)/(self.wave**2)).to((10**-20)*u.erg*(u.cm**-2)*(u.s**-1)*(u.angstrom**-1))
        #self.flam_err = ((self.noise * const.c)/(self.wave**2)).to((10**-20)*u.erg*(u.cm**-2)*(u.s**-1)*(u.angstrom**-1))

        #pieceper = self.file.split('.')
        #pieceundscr = pieceper[0].split('_')
        #self.id = pieceundscr[0] + "_" + pieceundscr[2] + "_" + pieceundscr[3]
        #self.e_id = "e_" + pieceundscr[3]


     def readfile(self, flxtype):
        
        #Opens the fits file & select index 1 (SPECID) where the data we wish to access lives
        data = fits.open(self.file)[1].data
    
        #Within data is 10 columns, we wish to access the columns titled, "wave", "flux", "err"
        self.wave = data["wave"] * u.micron
        self.fnu = data["flux"] * u.microjansky #fnu
        self.fnu_err = data["err"] * u.microjansky #err_fnu

        self.flam = ((self.flux * const.c)/(self.wave**2)).to((10**-20)*u.erg*(u.cm**-2)*(u.s**-1)*(u.angstrom**-1))
        self.flam_err = ((self.noise * const.c)/(self.wave**2)).to((10**-20)*u.erg*(u.cm**-2)*(u.s**-1)*(u.angstrom**-1))

        pieceper = self.file.split('.')
        pieceundscr = pieceper[0].split('_')
        self.id = pieceundscr[0] + "_" + pieceundscr[2] + "_" + pieceundscr[3]
        self.e_id = "e_" + pieceundscr[3]

        if flxtype == "fnu":
            self.flux = self.fnu
            self.noise = self.fnu_err
            self.flux_unit = u.microjansky
            self.flux_label = r'$f_{\nu}\$'
            
            
        elif flxtype == "flam":
            self.flux = self.flam
            self.noise = self.flam_err
            self.flux_unit = ((10**-20)*u.erg*(u.cm**-2)*(u.s**-1)*(u.angstrom**-1))
            self.flux_label = r'$f_{\nu}\$'
            self.flux_label = r'$f_{\lambda}$'

         
        self.variance = self.noise**2
    

        


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

def convertflux(spectrum):
#PROCEED WITH CAUTION THIS IS NOT DONE LIKE AT ALLLLLLL
#This function converts from your current flux type (fnu or flam) to the other. 
#1st goal is to make this into to seperate functions in case you're already in the units you're attempting to convert to. 
#2nd goal is to initialize a new spec object which holds the new flux type 
    if spectrum.flux_label == r"$f_{\nu}\$":
        if not spectrum.flux_unit == u.microjansky:
            spectrum.flux_unit = u.microjansky
            spectrum.flux = spectrum.flux.to(u.microjansky)
            spectrum.noise = spectrum.noise.to(u.microjansky)
        spectrum.flux = ((spectrum.flux * const.c)/(spectrum.wave**2)).to((10**-20)*u.erg*(u.cm**-2)*(u.s**-1)*(u.angstrom**-1))
        spectrum.noise = ((spectrum.noise * const.c)/(spectrum.wave**2)).to((10**-20)*u.erg*(u.cm**-2)*(u.s**-1)*(u.angstrom**-1))
        spectrum.flux_unit = (10**-20)*u.erg*(u.cm**-2)*(u.s**-1)*(u.angstrom**-1)
        spectrum.flux_label = r'$f_{\lambda}\$'

    elif spectrum.flux_label == r'$f_{\lambda}\$':
        if not spectrum.flux_unit == (10**-20)*u.erg*(u.cm**-2)*(u.s**-1)*(u.angstrom**-1):
            spectrum.flux_unit = (10**-20)*u.erg*(u.cm**-2)*(u.s**-1)*(u.angstrom**-1)
            spectrum.flux = spectrum.flux.to((10**-20)*u.erg*(u.cm**-2)*(u.s**-1)*(u.angstrom**-1))
        spectrum.flux = ((spectrum.wave**2)*(spectrum.flux))/(const.c).to(u.microjansky)
        spectrum.noise = (((spectrum.wave**2)*(spectrum.noise))/(const.c)).to(u.microjansky)
        spectrum.flux_unit = u.microjansky
        spectrum.flux_label = r"$f_{\nu}\$"
    else:
        print("something has gone wrong")

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

#def classifystandard(spec):
    #classifies by comparison to standard
    #think this only works for fnu because splat.classifyByStandard refers to self.wave & self.noise
    #return splat.classifyByStandard(spec)

def classifystandard(spec, standardset):
    #classifies by comparison to standard
    chisqr = 10000
    for standard in standardset:
        newchisqr = chisquare(spec, standard)
        if newchisqr < chisqr:
            chisqr = newchisqur
            return chisqr

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
