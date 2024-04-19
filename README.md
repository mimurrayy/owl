# owl spectroscopic library
  
Library to calculate optical emission spectra from atoms as well as line shapes.

## Installation
Owl is on [https://pypi.org/project/owlspec/](pypi) and can simply be installed as any python package:
```
pip install owlspec
```

For manual installation, you can also simply drop the `owlspec` folder in the folder where your script is located. However, you will still need to install the dependencies.

## Capabilities
The library has two core capabilities:
1. Fetching information about emission lines and levels from NIST and calculating complete PLTE spectra.
2. Calculating the broadening of emission lines for a selection of cases.

What can be done based on these core capabilities is shown in the examples below.

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

transition1 = owl.emitter.transition("H I", cw) # Hydrogen atom emission, H beta
# H I emitters are 6000 K hot, and sourrounded by argon neutrals and ions 
# and 6000 K hot electrons
line1 = owl.emission_line(transition1, cw, pert="Ar I", T=6000, Te=6000)
# setting T siwtches on Doppler broadening
# Te and the pert will be used for Stark broadening below

plt.figure() # Show influence of Stark broadening
x = np.linspace(cw-5, cw+7.5, 3000)
y1 = line1.get_profile(x, ne=5e21, N=1e20)
y2 = line1.get_profile(x, ne=2e22, N=1e20)
y3 = line1.get_profile(x, ne=5e22, N=1e20)

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
transition1 = owl.emitter.transition("He I", cw) # Helium atom emission
x = np.linspace(cw-0.1, cw+0.15, 3000)

# He I emitters are 500 K hot, and sourrounded by argon neutrals
line = owl.emission_line(transition1, cw, pert="Ar I", T=500)

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

plt.figure() # Show Zeeman splitting
cw = central_wavelength = 706.722
x = np.linspace(cw-0.05, cw+0.05, 3000)
transition = owl.emitter.transition("Ar I", cw)
line = owl.emission_line(transition, cw, T=2000)

y1 = line.get_profile(x, B = 0.2) # Now only Zeeman and Doppler
y2 = line.get_profile(x, B = 0.8)

plt.plot(x,y1/np.max(y1), label="B = 0.1 T")
plt.plot(x,y2/np.max(y2), label="B = 0.8 T")
plt.legend()
plt.xlabel("wavelength / nm")
plt.ylabel("normalized intensity")
```

<img src="https://github.com/mimurrayy/owl/assets/3911345/8eb81589-e755-48ed-8be5-d0960c9637df" width=400px>


### Stick spectrum for the intefication of emission lines

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
x,y = spec.get_ident_spectrum(min_int=200, min_Aik=1e6) # filtered for relative intensity and Aik value
x2,y2 = spec2.get_ident_spectrum() 
plt.plot(x,y, '-', label="Ti I")
plt.plot(x2,y2, '-', label="Ar II")
plt.legend()
plt.xlabel("wavelength / nm")
plt.ylabel("relative intensity (NIST)")
```

<img src="https://github.com/mimurrayy/owl/assets/3911345/c57b40c2-d596-40c9-bf70-101847030393" width=400px>


### Simulation of complete PLTE spectrum

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

