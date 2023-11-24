# owl
  
Library to calculate optical emission line shapes.

## Usage
The library has currently three main classes: 

- spectrometer
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

cw = central_wavelength = 453.37
x = np.linspace(cw-0.1, cw+0.1, 1000)

transition1 = owl.emitter.transition("Ti I", 453.324)
transition2 = owl.emitter.transition("Ti II", 453.396)

line1 = owl.emission_line(transition1,453.337,B=0.5)
line2 = owl.emission_line(transition2,453.410,B=0.5)

y1 = line1.get_profile(x, A=1, T=5000)
y2 = line2.get_profile(x, A=1, T=5000)

plt.plot(x,y1)
plt.plot(x,y2)
```    
