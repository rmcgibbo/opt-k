import os
import pandas as pd
import numpy as np
import statsmodels.formula.api as smf
import matplotlib
import matplotlib.pyplot as pp

# Set pandas matplotlib style
matplotlib.rcParams.update(pd.tools.plotting.mpl_stylesheet)
matplotlib.rcParams.update({'font.family':'sans-serif'})

pp.figure(figsize=(8,4))

pp.subplot(121)

pd.DataFrame(data=dict(x=[0,1], y=[0,1.2])).plot(x='x', y='y', style='.', ms=15, label='Data')
pp.plot([1], [1], label='Robust')
pp.plot([1], [1], label='Overfit')

pp.legend()

pp.subplot(122)
pp.plot([1], [1], label='Training Error')
pp.plot([1], [1], label='Generalization Error')


pp.xlabel('Model Complexity', size=16)
pp.ylabel('Error', size=16)
pp.legend()

pp.gcf().subplots_adjust(bottom=0.15)
pp.savefig('text.png')
#os.system('open text.png')
