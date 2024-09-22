import pandas as pd
import numpy as np
import scipy as sc
import keras as k
import tensorflow as tf
import matplotlib.pyplot as plt
import statsmodels.api as sm
from FE_Library import yearfrac

# Load Time Series
quotes = pd.read_csv("DatasetPythonAss3.csv")
quotes["Date"] = pd.to_datetime(quotes["Date"], format="%d/%m/%Y")
quotes = quotes.set_index("Date")
AAPL = quotes["AAPL"]
SPX = quotes["SPX"]

# Plot Time Series
# get dates from the table
AAPL_dates = AAPL.index
SPX_dates = SPX.index
plt.subplot(2,1,1)
# plot the time series
# use a subplot as the values have different orders of magnitude
plt.plot(AAPL_dates, AAPL, label = 'AAPL')
plt.title('Time Series: AAPL')
plt.subplot(2,1,2)
plt.plot(SPX_dates, SPX, label = 'SPX')
plt.title('Time Series: SPX')


# Compute Log_returns
# transform the arrays in numpy arrays in order to do matrix operations
AAPL = np.array(AAPL)
SPX = np.array(SPX)
# compute the logReturns
returnsAAPL = AAPL[1:]/AAPL[:-1]
returnsSPX = SPX[1:]/SPX[:-1]
logRetAAPL = np.log(returnsAAPL)
logRetSPX = np.log(returnsSPX)
# plot the logReturns
plt.figure()
plt.plot(AAPL_dates[1:], logRetAAPL, label = 'logReturnAAPL')
plt.plot(SPX_dates[1:], logRetSPX, label = 'logReturnSPX')
plt.legend()
plt.title('logReturns')

# Regressions
# find the fitting first-degree polynomial and get the slope
slope, intercept, r, p, err = sc.stats.linregress(logRetSPX, logRetAAPL)
plt.figure()
# plot the fitting polynomial
plt.scatter(logRetSPX, logRetAAPL)
plt.plot(logRetSPX, slope*logRetSPX+intercept)
plt.title('Linear Regression')
plt.xlabel('logReturnSPX')
plt.ylabel('logReturnAAPL')
print('slope: ', slope)

# YearFrac
Act365 = 3
# use the given yearfrac function
yf = yearfrac(AAPL_dates[0], AAPL_dates[-1], 3)
print('yearfrac: ', yf)

# Interpolate
y = [1, 2, 3.5, 4, 2]
x = [0, 1, 2, 3, 4]
# find the linear interpolating function
f = sc.interpolate.interp1d(x, y)
# print the interpolated value
print('interpolated value: ', f(2.7))

# Simulation
# set the seed
np.random.seed(2)
# get the standard Normal random variable
standardGaussian = np.random.randn()
print('Standard Normal random variable: ', standardGaussian)
n = 10000
# create increasing array
xG = np.arange(n) + 1
# initialization of the array
yG = np.zeros(n)
# populate a random variable array of increasing dimensions
# and then compute variances for those arrays
for i in xG:
    yrand = np.random.randn(xG[i-1])
    yG[i-1] = np.var(yrand)
# plot the variances
plt.figure()
plt.plot(xG, yG)
plt.title('Gaussian Variance plot')
plt.xlabel('Number of iterations')
plt.ylim(0.8,1.2)

# Find the quantile
quantile = np.quantile(yrand, 0.9)
# evaluate the Gaussian CDF at quantile 0.9
value_quantile = sc.stats.norm.cdf(quantile)
print('Quantile of 0.9 computed: ', quantile)
print('Gaussian CDF evaluated at quantile: ', value_quantile)

# Minimization
def f(x):
    return (x[0] - 3)**2 + (x[1] - 7)**2

# Analitically
# compute the gradient of the function and find the stationary point [3,7],
# then, deriving again, verify that it's the minimum

# Numerically
# using the minimizing function
minimum = sc.optimize.minimize(f, [0, 0])
print('Minimum : [x y] = ', minimum.x, ', value of the function = ', minimum.fun)
plt.show()
