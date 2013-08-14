import os
import pandas as pd
import numpy as np
import statsmodels.formula.api as smf
import matplotlib
import matplotlib.pyplot as pp

# Set pandas matplotlib style
matplotlib.rcParams.update(pd.tools.plotting.mpl_stylesheet)
matplotlib.rcParams.update({'font.family':'sans-serif'})

# 3rd degree polynomial
coeff = [0.4, -0.7, 0.3, 0]
sigma = 0.01

x = np.linspace(0, 1, 12)
y = np.polyval(coeff, x) + sigma*np.random.randn(len(x))


pp.figure(figsize=(8,4))#figsize=(3,3))
pp.subplot(121)

df = pd.DataFrame(data=dict(x=x, y=y))
for i in range(2, 11):
    df['x%d' % i] = df['x']**i
df.plot(x='x', y='y', style='.', label='data', ms=15)

m3 = smf.ols('y ~ x + x2 + x3', data=df).fit()
m6 = smf.ols('y ~ x + x2 + x3 + x4 + x5 + x6', data=df).fit()

df2 = pd.DataFrame(data=dict(x=np.linspace(0,1,1000)))
for i in range(2, 11):
    df2['x%d' % i] = df2['x']**i

df2['y_pred_3'] = m3.predict(df2)
df2['y_pred_6'] = m6.predict(df2)
df2.plot(x='x', y='y_pred_3', label='3rd degree')
df2.plot(x='x', y='y_pred_6', label='6th degree')
pp.xlabel('Model Inputs', size=16)
pp.ylabel('Predictions', size=16)
pp.legend(loc=1)
pp.gcf().subplots_adjust(bottom=0.15)

pp.subplot(122)


name = __file__ + '.png'
pp.savefig(name)
os.system('open %s' % name)

