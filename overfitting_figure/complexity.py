import os
import numpy as np
import matplotlib.pyplot as pp
from xkcdplot import XKCDify

x = np.linspace(0.1, 5)
training = 1/x + 1.0
test = 1/x + 1.4 + 0.2*(x-1)**2

pp.plot(x+0.1, test, label='test error')
pp.plot(x, training, label='training error')

#pp.legend()
#pp.ylabel('Error')
#pp.xlabel('Model Complexity')
XKCDify(pp.gca())






name = __file__ + '.png'
pp.savefig(name)
os.system('open %s' % name)

