{
 "cells": [
  {
   "cell_type": "code",
   "execution_count": 6,
   "id": "a51cb7ed-8266-4507-b78f-6c999e155448",
   "metadata": {},
   "outputs": [],
   "source": [
    "from astropy.io import fits\n",
    "import astropy.units as u\n",
    "from astropy import constants as const\n",
    "import matplotlib.pyplot as plt\n",
    "import numpy as np\n",
    "import splat"
   ]
  },
  {
   "cell_type": "code",
   "execution_count": 55,
   "id": "b7829cdf-2ec5-4556-ab4e-1f6d2dfa9c32",
   "metadata": {},
   "outputs": [],
   "source": [
    "class spectrum(splat.Spectrum):\n",
    "    def __init__(self):\n",
    "        self.fnu = None\n",
    "        self.fnu_err = None\n",
    "        self.flam = None\n",
    "        self.flam_err = None\n",
    "\n",
    "    def read_file(self, file_name):\n",
    "        #Opens the fits file by:\n",
    "        hdul = fits.open(file_name)\n",
    "    \n",
    "        #The data we wish to access is in index 1 (SPEC1D), so we can select index to it by:\n",
    "        data = hdul[1].data\n",
    "    \n",
    "        #Within data is 10 columns, we wish to access the columns titled, \"wave\", \"flux\", \"err\"\n",
    "        self.wave = data[\"wave\"] * u.micron\n",
    "        self.fnu = data[\"flux\"] * u.microjansky\n",
    "        self.fnu_err = data[\"err\"] * u.microjansky\n",
    "        \n",
    "        lam_flux = self.fnu*(const.c) / self.wave**2\n",
    "        lam_err = self.fnu_err*(const.c) / self.wave**2\n",
    "\n",
    "        self.flam = lam_flux.to((1*10**(-20))*u.erg*(1/u.s)*(1/u.cm**2)*(1/u.AA))\n",
    "        self.flam_err = lam_err.to((1*10**(-20))*u.erg*(1/u.s)*(1/u.cm**2)*(1/u.AA))\n",
    "\n",
    "    def plot(self, fluxtype):\n",
    "        plt.figure(figsize=(9,4))\n",
    "        if fluxtype == \"flam\n",
    "        plt.plot(self.wave, self.flam)\n",
    "        plt.plot(self.wave, self.flam_err)\n",
    "        plt.xlabel(\"Wavelength [${\\mu}m$]\")\n",
    "        plt.ylabel(\"$f_{\\lambda} [10^-20 ergs^-1 cm^-2 {\\AA}^-1]$\")\n",
    "        return plt.show()\n",
    "        "
   ]
  }
 ],
 "metadata": {
  "kernelspec": {
   "display_name": "Python 3 (ipykernel)",
   "language": "python",
   "name": "python3"
  },
  "language_info": {
   "codemirror_mode": {
    "name": "ipython",
    "version": 3
   },
   "file_extension": ".py",
   "mimetype": "text/x-python",
   "name": "python",
   "nbconvert_exporter": "python",
   "pygments_lexer": "ipython3",
   "version": "3.11.13"
  }
 },
 "nbformat": 4,
 "nbformat_minor": 5
}
