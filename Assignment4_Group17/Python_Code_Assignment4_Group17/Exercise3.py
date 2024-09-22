# ----- EXERCISE 3 ------
import pandas as pd
import numpy as np
import datetime as dt
from datetime import datetime
import matplotlib.pyplot as plt
import Utilities as u
import FE_Library as FEL
import random
import scipy as sp
import time


start_time = time.time()
# Setting the seed:
np.random.seed(0)
# Read Bootstrap results:
bootstrap = pd.read_csv("BootstrapResults.csv")
# Read CDS Bootstrap results for ISP:
ISPSurvProb = pd.read_csv("ISPSurvProb.csv")
# Option Cash flow dates:
optionDates = ["2009-02-19","2010-02-19","2011-02-21","2012-02-20","2013-02-19","2014-02-19","2015-02-19"]
# Discount at option dates calculation:
discounts = u.DiscountsAtOptionDates(optionDates,bootstrap)
# Forward discounts:
fwdDiscounts = discounts[1:]/discounts[:-1]
# Add settlement date for computation:
optionDates.insert(0,"2008-02-19")
date_objects = [datetime.strptime(date, '%Y-%m-%d') for date in optionDates]
# Calculation year fraction between consecutive dates
year_frac = np.zeros(len(optionDates)-1)
for i in range(len(optionDates)-1):
    year_frac[i] = FEL.yearfrac(pd.to_datetime(optionDates[i]),pd.to_datetime(optionDates[i+1]), basis=6)
year_frac = np.array(year_frac)
# Notional of the option:
notional = 30000000
# Participation coefficient:
L = 0.99
# Underlying volatility:
sigma = 0.2
# Recovery value:
pi = .40
# Forward starting Call initialization:
call_prices = np.zeros(len(year_frac))
# Forward starting Call computation:
for i in range(len(call_prices)):
    # Formula for Call options with Fwd discount:
    call_prices[i] = u.CallPriceComputationGeman(1,1/L,fwdDiscounts[i],0,sigma,year_frac[i])
cliquetCleanPrice = call_prices.sum()*L*notional

# Number of Monte-Carlo simulations:
N = 100000
# Underlying initialization:
S = np.zeros((N, len(optionDates)))
S[:, 0] = 1
# Gaussian random variable initialization, used as diffusion in the simulation of GBM:
Z = np.random.randn(N, len(optionDates))
# Payoff initialization:
Payoff = np.zeros((N, len(year_frac)))
# Forward rate calculation from forward discounts:
fwdRate = -np.log(fwdDiscounts) / year_frac
# Option date for loop:
for i in range(len(fwdDiscounts)):
    # GBM simulation
    S[:, i + 1] = S[:, i] * np.exp((fwdRate[i] - sigma ** 2 / 2) * year_frac[i] + sigma * np.sqrt(year_frac[i]) * Z[:, i])
    # Option payoff at option date i
    Payoff[:, i] = np.maximum(L * S[:, i + 1] - S[:, i], 0)
# Discounted payoff for each forward start Call:
DiscountedPayoff = Payoff * discounts[1:].reshape(1, -1)
# Discounted payoff for Cliquet option:
DiscountedPayoff = np.sum(DiscountedPayoff,axis=1)
# Monte-Carlo Price and standard deviation:
cliquetCleanPriceMC, cliquetMCstd = sp.stats.norm.fit(DiscountedPayoff * notional)
# Quantile for 95% CI:
z = sp.stats.norm.ppf(0.975)
# Cliquet lower bound:
cliquetBid = cliquetCleanPriceMC - z * cliquetMCstd/np.sqrt(N)
# Cliquet upper bound:
cliquetAsk = cliquetCleanPriceMC + z * cliquetMCstd/np.sqrt(N)
# Survival probability:
surv = ISPSurvProb['SurvivalProbabilities'].values
# Survival probability at settlement date:
surv = np.insert(surv,0,1)
# Cliquet price with counterparty risk:
cliquetDirtyPrice = L*notional*((np.dot(call_prices,surv[1:].T)) + pi * np.dot(call_prices[::-1].cumsum(), (surv[:-1] - surv[1:]).T))

# Monte-Carlo for dirty price:
# Simulation of default:
samples = np.random.uniform(0, 1, N)
# Initialization of the vector  of defaults:
defaultTime = np.zeros(N)
for i in range(N):
    defaultTime[i] = np.argmax(samples[i] > surv)
defaultIndeces = np.nonzero(defaultTime)
for idx in defaultIndeces:
    time_idx = defaultTime[idx][0]  # Get the first element of the array
    Payoff[idx, int(time_idx):]= 0
    nettingValue = 0
    SDefault = S[idx,int(time_idx)]
    for it in range(len(fwdDiscounts) - int(time_idx)-1):
        # Formula for call options with Fwd discount:
        nettingValue += u.CallPriceComputationGeman(SDefault, SDefault / L, fwdDiscounts[it+1+int(time_idx)], 0, sigma, year_frac[it +1 +int(time_idx)])
    Payoff[idx, int(time_idx)] = pi * nettingValue
# Discounted payoff for each forward start Call:
DiscountedPayoff = Payoff * discounts[1:].reshape(1, -1)
# Discounted payoff for Cliquet option:
DiscountedPayoff = np.sum(DiscountedPayoff,axis=1)
# Monte-Carlo Price and standard deviation:
cliquetDirtyPriceMC, cliquetDirtyMCstd = sp.stats.norm.fit(DiscountedPayoff*notional)
# Quantile for 95% CI:
z = sp.stats.norm.ppf(0.975)  # Quantile for 95% CI
# Cliquet lower bound:
cliquetDirtyBid = cliquetDirtyPriceMC - z * cliquetDirtyMCstd/np.sqrt(N)
# Cliquet upper bound:
cliquetDirtyAsk = cliquetDirtyPriceMC + z * cliquetDirtyMCstd/np.sqrt(N)
# Clock end and computational time:
end_time = time.time()
computational_time = end_time - start_time

# Print the results:
print('\n\n----------- EXERCISE 3 -----------')
print("Cliquet clean price via closed formula: ", round(cliquetCleanPrice,2), "€")
print("Cliquet clean price via MonteCarlo :", round(cliquetCleanPriceMC,2), "€")
print("Bid: ", round(cliquetBid,2), "€        Ask: ", round(cliquetAsk,2), "€")
print("Cliquet dirty price via closed formula: ", round(cliquetDirtyPrice,2),"€")
print("Cliquet dirty price via MonteCarlo :", round(cliquetDirtyPriceMC,2), "€")
print("Bid: ", round(cliquetDirtyBid,2), "€        Ask: ", round(cliquetDirtyAsk,2),"€")
print("MC simulation time: ", computational_time)