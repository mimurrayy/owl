# owl
  
Library to calculate optical emission line shapes.

## Usage
The library has currently four main classes: 

- spectrometer
- particle
- transition
- emission_line

## Example
```python
#!/bin/python3
import matplotlib.pyplot as plt
import numpy as np

import sys
import os.path
home = os.path.expanduser('~')
sys.path.append(home + '/.local/lib/')
import owl

cw = central_wavelength = 453.4
pgs = owl.spectrometer.pgs(cw,3)
x = pgs.x
ti_atom = owl.emitter.particle("TiI")
ti_ion = owl.emitter.particle("TiII")
transition1 = owl.emitter.transition(ti_atom, 453.324)
transition2 = owl.emitter.transition(ti_ion, 453.396)
transition3 = owl.emitter.transition(ti_atom, 453.478)
line1 = owl.emission_line(pgs,transition1,453.337,B=0.1)
line2 = owl.emission_line(pgs,transition2,453.410,B=0.1)
line3 = owl.emission_line(pgs,transition3,453.491,B=0.1)
emission_lines = [line1,line2,line3]
x = pgs.make_x_scale(cw,order=3)
y = line1.get_profile(x, A=1, T=1000)

plt.plot(x,y)
```    
