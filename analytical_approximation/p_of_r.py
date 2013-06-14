from __future__ import division
import numpy as np
import matplotlib.pyplot as pp

r = np.linspace(0,1)
for N in range(3,10,2):
    p = N / (2**(N/2 + 1)) * (1/r**2 - 1/2)**(-N/2 - 1)
    np.log(1/r**2 - 0.5)
    pp.plot(r,p, lw=2, label='N=%d' % N)

pp.ylabel('p(r)')
pp.xlabel('r')
pp.legend(loc=2)
#pp.show()
pp.savefig('p_of_r.png')

