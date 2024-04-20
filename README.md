# owl spectroscopic library
  
Library to calculate optical emission spectra from atoms as well as line shapes.

## Installation
Owl is on [pypi](https://pypi.org/project/owlspec/) and can simply be installed as any python package:
```
pip install owlspec
```

For manual installation, you can also simply drop the `owlspec` folder in the folder where your script is located. However, you will still need to install the dependencies (numpy, scipy, mendeleev, astroquery, roman).

## Capabilities
The library has two core capabilities:
1. Fetch information about transitions and levels from NIST.
2. Calculating the broadening of emission lines for a selection of cases.
3. Calculating complete PLTE spectra.

### Levels and transitions data
Owl can be used to access information from the NIST atomic spectra database[^NIST], which we access using [astroquery](https://astroquery.readthedocs.io/en/latest/). That does mean _owl requires internet access to function_. To access information about a transition, create a transition object by supplying the emitter name in spectroscopic notation (e.g. "Ar I" for argon neutrals or "Fe IV" for triply charged iron ions) and the wavelength in nm:
```python
transition = owl.emitter.transition("Ar I", 751.5)
```
The object then contains the information as attributes. For example, the energy of the upper level of the transition can be accessed as `transition.upperE`, the Einstein coefficient for emission as `transition.Aik` and so on. Possible attributes are: `name`, `charge`, `spec_name`, `wl`, `upperE`, `lowerE`, `upperl`, `lowerl`,`Aik`,`upperg`, `lowerg`. Additionally, the levels corresponding to the transition can be accessed as `transition.upper` and `transition.lower`. The levels contain `J` and the Lande g factor `G`.

Alternatively, the energy levels can also be accessed directly by specifying the energy (in eV):
```python
level = owl.emitter.level("Ar I", 11.623)
```
The level contains attributes `E`, `J`, `l` and `G` and `conf`. 

By itself, the level and transition data is not very useful, but would be helpful in, for example, assembling a collisional radiative model. However, the level class contains a method that might be of immediate use: `get_lifetime()` which calculates the radiative lifetime of the level from all transitions listed in NIST. The result is returned in units of seconds.

### Line broadening calculations

Owl supports a range of different line broadening mechanisms that are automatically activated when providing information about the physical situation surrounding the emitters. For example, specifying an electron density will automatically switch on Stark broadening calculations, if they are available for the emitter. We currently support Doppler, Stark, Zeeman, van der Waals and instrumental broadening. Self-resonance broadening is not yet supported. Multiple broadening mechanisms are combined by numerical convolution of the individual profiles.

#### Stark broadening
Since Stark broadening is different for each transition and each emitter species and no generalized theory is available, owl only supports a few selected transitions.

Owl uses the Stark broadening calculations of Gigosos et al for H alpha, beta and gamma, as well as He 447.1 nm and He 492.2 nm. The precalculated tables are downloaded from the publishers when the respective functions are first executed. The space in between the calculated datapoints is interpolated as advised by the respective authors.

Additionally, we have support for O 777 nm, Ar 810.369 nm and Ar 738.398 nm based on Griems tabulated constants[^1]. Extending the library with more data from Griems calculations for other transitions would be trivial, but doing everything would be a tremendous amount of busywork. Thus, if you need support for any specific line, let me know and I would be happy to put it in.

| Transition | Wavelength | Source                                                                                                                                         |
| ---------- | ---------- | ---------------------------------------------------------------------------------------------------------------------------------------------- |
| Hα         | 656.3 nm   | [Gigosos et al 2003 _Spectrochimica Acta Part B: Atomic Spectroscopy_ **58** 1489–504](https://doi.org/10.1016/S0584-8547(03)00097-1)          |
| Hβ         | 486.1 nm   | [Gigosos et al 2003 _Spectrochimica Acta Part B: Atomic Spectroscopy_ **58** 1489–504](https://doi.org/10.1016/S0584-8547(03)00097-1)          |
| Hγ         | 434.0 nm   | [Gigosos et al 2003 _Spectrochimica Acta Part B: Atomic Spectroscopy_ **58** 1489–504](https://doi.org/10.1016/S0584-8547(03)00097-1)          |
| He I       | 447.1 nm   | [Gigosos et al 2009 _A&A_ **503** 293–9](https://doi.org/10.1051/0004-6361/200912243)                                                          |
| He I       | 492.2 nm   | [Lara N et al 2012  _A&A_ **542** A75](https://dx.doi.org/10.1051/0004-6361/201219123)                                                         |
| Ar I       | 810.369 nm | [H.R. Griem: Spectral Line Broadening by Plasmas](https://shop.elsevier.com/books/spectral-line-broadening-by-plasmas/griem/978-0-12-302850-1) |
| Ar I       | 738.398 nm | [H.R. Griem: Spectral Line Broadening by Plasmas](https://shop.elsevier.com/books/spectral-line-broadening-by-plasmas/griem/978-0-12-302850-1) |

[^1]: H.R. Griem: Spectral Line Broadening by Plasmas

#### Zeeman splitting
Splitting due to magnetic fields is calculated analytically following the books of Cowan[^Cowan] as well as Condon and Shortley[^CS]. The necessary base data (as the Lande g factors) are obtained from the NIST atomic spectra database [^NIST]. Transitions for which the base data is missing from NIST are not currently supported. Polarization of the different components is calculated, but is currently not exposed to the user. Let me know if you need that feature.

[^Cowan]: Cowan R D 1981 _The theory of atomic structure and spectra_ (Berkeley: University of California Press)
[^CS]: Condon E U and Shortley G H 1979 _The theory of atomic spectra_ (Cambridge: Univ. Pr)
[^NIST]: https://www.nist.gov/pml/atomic-spectra-database

#### Doppler broadening
Doppler broadening currently only supports Maxwell VDFs.

#### van der Waals broadening
Van der Waals broadening is calculated using the equations from Konjević [^Konj]. The polarizability of atoms surrounding the emitters is automatically obtained using [mendeleev](https://mendeleev.readthedocs.io/en/stable/). 

[^Konj]: Konjević  1999 _Physics Reports_ **316** 339–401

#### Instrumental broadening
Instrumental broadening is included either using a pseudo Voigt function with user-specified width and shape parameter or by passing a function that takes the x-axis and the central wavelength position as arguments. The example for Zeeman broadening below shows how to do that.

### Spectra simulations

The library queries the NIST atomic spectra database[^NIST] to automatically obtain base data about transitions and levels. This can be used to calculate complete spectra. For example, to obtain the spectrum of chromium neutrals between 300 nm and 500 nm in partial local thermal equilibrium at 3 eV, all you need to do is:

```python
spec = owl.spectrum('Cr I', wl_range=[300,500])
y = spec.get_LTE_spectrum(x, Te=3, width=0.2, mu=0.5, norm=True) 
```

This capabilities can be used to quickly identify unknown lines in measurements or fit measured spectra to obtain an excitation temperature in cases where Boltzmann plots are difficult due to insufficient resolution.


## Examples

The following examples show most of what owl is currently capable of. All examples should be completely copy-pasteable, if you have the library and matplotlib installed.

### H beta Stark broadening

```python
#!/bin/python3
import numpy as np
import matplotlib.pyplot as plt
import owlspec as owl

cw = central_wavelength = 486
x = np.linspace(cw-0.5, cw+0.5, 3000)

transition = owl.emitter.transition("H I", cw) # Hydrogen atom emission, H beta
# H I emitters are 6000 K hot, and sourrounded by argon neutrals and ions 
# and 6000 K hot electrons
line = owl.emission_line(transition, cw, pert="Ar I", T=6000, Te=6000)
# setting T siwtches on Doppler broadening
# Te and the pert will be used for Stark broadening below

plt.figure() # Show influence of Stark broadening
x = np.linspace(cw-5, cw+7.5, 3000)
y1 = line.get_profile(x, ne=5e21, N=1e20)
y2 = line.get_profile(x, ne=2e22, N=1e20)
y3 = line.get_profile(x, ne=5e22, N=1e20)

plt.plot(x,y1/np.max(y1), label=r"n$_e$ = $5 \times 10^{21}$ m$^{-3}$")
plt.plot(x,y2/np.max(y2), label=r"n$_e$ = $2 \times 10^{22}$ m$^{-3}$")
plt.plot(x,y3/np.max(y3), label=r"n$_e$ = $5 \times 10^{22}$ m$^{-3}$")
plt.legend()
plt.xlabel("wavelength / nm")
plt.ylabel("normalized intensity")
```    

<img src="https://github.com/mimurrayy/owl/assets/3911345/2f150413-2e24-49b3-978f-5cc08828cee5" width=400px>


### van der Waals and Doppler broadening of a helium line

```python
#!/bin/python3
import numpy as np
import matplotlib.pyplot as plt
import owlspec as owl

plt.figure() # Show influence of van der Waals broadening
cw = central_wavelength = 501.567
transition = owl.emitter.transition("He I", cw) # Helium atom emission
x = np.linspace(cw-0.1, cw+0.15, 3000)

# He I emitters are 500 K hot, and sourrounded by argon neutrals
line = owl.emission_line(transition, cw, pert="Ar I", T=500)

y1 = line.get_profile(x, N=1e24) # N is the neutral density of the perturber (Ar)
y2 = line.get_profile(x, N=5e24) # setting N switches on vdW broadening
y3 = line.get_profile(x, N=1e25) # T (defined as 500 above) siwtches on Doppler
plt.plot(x,y1/np.max(y1), label=r"n$_g$ = $1 \times 10^{24}$ m$^{-3}$")
plt.plot(x,y2/np.max(y2), label=r"n$_g$ = $5 \times 10^{24}$ m$^{-3}$")
plt.plot(x,y3/np.max(y3), label=r"n$_g$ = $1 \times 10^{25}$ m$^{-3}$")
plt.legend()
plt.xlabel("wavelength / nm")
plt.ylabel("normalized intensity")


plt.figure() # Show influence of Doppler broadening
y1 = line.get_profile(x, T=1000) # Doppler broadening
y2 = line.get_profile(x, T=10000)

plt.plot(x,y1/np.max(y1), label="T$_g$ = 1000 K")
plt.plot(x,y2/np.max(y2), label="T$_g$ = 10000 K")
plt.legend()
plt.xlabel("wavelength / nm")
plt.ylabel("normalized intensity")
```
<img src="https://github.com/mimurrayy/owl/assets/3911345/0aab3190-e863-4644-aab7-266676493d01" width=400px><img src="https://github.com/mimurrayy/owl/assets/3911345/942fc7d7-ba98-4238-a8ac-7f225a89fa24" width=400px>


### Zeeman splitting of Ar I 706.722 nm


```python
#!/bin/python3
import numpy as np
import matplotlib.pyplot as plt
import owlspec as owl

def instr(x, xc): # Function to use for instrumental broadening
    return owl.util.psd_voigt(x, xc, w=0.003, mu=0.2)

plt.figure() # Show Zeeman splitting
cw = central_wavelength = 706.722
x = np.linspace(cw-0.05, cw+0.05, 3000)
transition = owl.emitter.transition("Ar I", cw)
line = owl.emission_line(transition, cw, instr_func=instr)


y1 = line.get_profile(x, B = 0.2) # Now only Zeeman and Doppler
y2 = line.get_profile(x, B = 0.8)

plt.plot(x,y1/np.max(y1), label="B = 0.1 T")
plt.plot(x,y2/np.max(y2), label="B = 0.8 T")
plt.legend()
plt.xlabel("wavelength / nm")
plt.ylabel("normalized intensity")
```

<img src="https://github.com/mimurrayy/owl/assets/3911345/8eb81589-e755-48ed-8be5-d0960c9637df" width=400px>


### Stick spectrum for the identification of emission lines
In order to identify unknown lines in a measured spectrum, it's useful to compare the measurement to lines reported in NIST using their reported relative intensities. To this end, `owl.spectrum.get_ident_spectrum()` can be used, which marks the line position and intensity with a single dot, which is much faster than calculating realistic line shapes. Intensities according to PLTE are also supported, by using `owl.spectrum.get_ident_spectrum_LTE(Te)`, instead.

```python
#!/bin/python3
import numpy as np
import matplotlib.pyplot as plt
import owlspec as owl

plt.figure() # Ident spectrum

 # Select all Ti lines between 600 and 900 nm listed by NIST
spec = owl.spectrum('Ti I', wl_range=[310,415])
spec2 = owl.spectrum('Ar II', wl_range=[310,415])

# simulate spectra with a given line braodening (w, mu) and normalized to the maximum
# min_int and min_Aik allow filtering for NIST-reported relative intensity and Aik value
x,y = spec.get_ident_spectrum(min_int=200, min_Aik=1e6) 
x2,y2 = spec2.get_ident_spectrum() 
plt.plot(x,y, '-', label="Ti I")
plt.plot(x2,y2, '-', label="Ar II")
plt.legend()
plt.xlabel("wavelength / nm")
plt.ylabel("relative intensity (NIST)")
```

<img src="https://github.com/mimurrayy/owl/assets/3911345/c57b40c2-d596-40c9-bf70-101847030393" width=400px>


### Simulation of a complete PLTE spectrum

The example below shows how to calculate a complete LTE spectrum. The function is fast enough that it can be used to fit measured data to obtain an excitation temperature.

```python
#!/bin/python3
import numpy as np
import matplotlib.pyplot as plt
import owlspec as owl

plt.figure() # PLTE Spectrum simulation
x = np.linspace(310,415,10000)

 # Select all Ti lines between 600 and 900 nm listed by NIST
spec = owl.spectrum('Ti I', wl_range=[310,415])

# simulate spectra with a given line braodening (w, mu) and normalized to the maximum
y = spec.get_LTE_spectrum(x, Te=3, width=0.2, mu=0.5, norm=True) 
y2 = spec.get_LTE_spectrum(x, Te=1, width=0.2, mu=0.5, norm=True)
plt.plot(x,y, '-', label="T$_e$ = 3 eV")
plt.plot(x,y2, '-', label="T$_e$ = 1 eV")
plt.legend()
plt.xlabel("wavelength / nm")
plt.ylabel("normalized intensity")
```

<img src="https://github.com/mimurrayy/owl/assets/3911345/99354512-8a80-411b-93be-2e4d89b42ce9" width=400px>

