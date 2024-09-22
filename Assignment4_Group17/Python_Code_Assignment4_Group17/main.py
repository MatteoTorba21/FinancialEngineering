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

# Import dataset with market data:
dataset = pd.read_csv('EUROSTOXX50_Dataset.csv')
dataset['Date'] = pd.to_datetime(dataset['Date'])
dataset.set_index('Date', inplace = True)

# Import dataset with Names and Tickers:
titles = pd.read_csv('_indexes.csv')
titles = titles.reindex(columns=['0','Name','Ticker'])
titles.set_index('Name', inplace=True)
titles.drop('0', axis= 1, inplace=True)


# ----- EXERCISE 0 ------

# Vector of the firms names:
Names_0 = ['Adidas', 'Allianz', 'Munich Re', "L'Oréal"]
# Fill the dataset: add previous day value in case of missing share price
dataset = dataset.ffill()
# Define today's date:
settlement = dt.datetime.strptime('2020-02-20', '%Y-%m-%d').date()
# Define the start date:
startDate = settlement - dt.timedelta(days = 365*5+1) # Add 1 for the leap year 2016
# Build a matrix with the values of the stocks for the firms and the period considered
dataArray_0 = u.getTablenp(dataset, titles, settlement, startDate, Names_0)

# Significance level:
alpha_0 = 0.99
# Vector of weights of the portfolio:
weights_0 = np.array([1/len(Names_0)]*len(Names_0))
# Notional of the portfolio:
Notional_0 = 15*10**6
# To have daily values of the risk measure set:
delta_0 = 1
# Compute the log returns of the stocks:
logReturns_0 = u.logReturnsComputation(dataArray_0)
# Degrees of freedom of the t-Student distribution:
nu = 4
# Clock start:
start_time = time.time()
# Compute the VaR and the ES with tha analytical approach with t-Student:
VaR_0, ES_0 = u.AnalyticaltStudentMeasures(alpha_0, weights_0, Notional_0, delta_0, logReturns_0, nu)
# Clock end and computational time:
end_time = time.time()
computational_time = end_time - start_time

# Print the results:
print('----------- EXERCISE 0 -----------')
print("Equally weighted equity portfolio with Adidas, Allianz, Munich Re and L’Oreal")
print('\nCheck for the t-student distribution of the loss')
u.verifytstudentdistributionhypothesis(weights_0,logReturns_0,nu)
print('\nVaR and the ES via analytic approach with t-Student:')
print('VaR: ', VaR_0,'     ES:', ES_0)
print('\nComputational time: ', computational_time)



# ----- EXERCISE 1 ------

# Significance level:
alpha_1 = 0.95
# Define today's date:
settlement = dt.datetime.strptime('2019-03-20', '%Y-%m-%d').date()
# Define the start date:
startDate = settlement - dt.timedelta(days = 365*5+1)

# QUESTION A:
# To have daily values of the risk measure set:
delta_1A = 1
# Vectors with the number of shares of Total, AXA, Sanofi and Volkswagen respectively:
nshares = {'TotalEnergies': 25000, 'AXA': 20000, 'Sanofi': 20000, 'Volkswagen Group': 10000}
# Vector of the firm Names:
Names_1A = nshares.keys()
# Vector of indices:
valsNames = titles.loc[nshares.keys(), 'Ticker']
# Value of shares for each stock:
weightsDirty_1A = [nshares[n]*dataset.loc[[settlement],titles.loc[[n], 'Ticker'].values[0]].values[0] for n in Names_1A]
# Build a matrix with the values of the stocks for the firms and the period considered:
dataArray_1A = u.getTablenp(dataset, titles, settlement, startDate, Names_1A)
# Compute portfolio value:
portfolioValue_1A = sum(weightsDirty_1A)
# Compute portfolio weights:
portfolioWeights_1A = np.array(weightsDirty_1A)/portfolioValue_1A
# Compute the log returns of the stocks:
logReturns_1A = u.logReturnsComputation(dataArray_1A)
# Clock start:
start_time = time.time()
# Compute the Var and the ES via a Historical Simulation approach:
VaR_1A_HS, ES_1A_HS = u.HSMeasurements(logReturns_1A, alpha_1, portfolioWeights_1A, portfolioValue_1A, delta_1A)
# Clock end and computational time:
end_time = time.time()
computational_time_HS = end_time - start_time
# Setting the seed:
np.random.seed(1)
# Number of Bootstrap samples:
numberOfSamples = 200
# Selection of the random samples:
samples = u.bootstrapStatistical(numberOfSamples, logReturns_1A)
# Clock start:
start_time = time.time()
# Compute the VaR and the ES via a statistical Bootstrap using Historical simulation:
VaR_1A_SB, ES_1A_SB = u.HSMeasurements(samples, alpha_1, portfolioWeights_1A, portfolioValue_1A, delta_1A)
# Clock end and computational time:
end_time = time.time()
computational_time_Bootstap = end_time - start_time

# Plausibility check of the VaR:
VaR_1A_PC = u.plausibilityCheck(logReturns_1A, portfolioWeights_1A, alpha_1, portfolioValue_1A, delta_1A)

# Print the results:
print('\n\n----------- EXERCISE 1 -----------')
print('\nQUESTION A:')
print("Portfolio with:25000 shares of Total,20000 shares of AXA,20000 shares of Sanofi and 10000 shares of Volkswagen")
print('\nVaR and the ES via Historical Simulation:')
print('VaR: ', VaR_1A_HS,'     ES:', ES_1A_HS)
print('\nVaR and the ES via Statistical Bootstrap using Historical Simulation:')
print('VaR_SB: ', VaR_1A_SB,'     ES_SB:', ES_1A_SB)
print('\nPlausibility check of the VaR in the question A: ', VaR_1A_PC)
print('\nComputational time for full HS: ', computational_time_HS)
print('Computational time for Bootstrap HS: ', computational_time_Bootstap)


# Question B:
# Vector of the firms names:
Names_1B = ['Adidas', 'Airbus', 'BBVA', 'BMW', 'Deutsche Telekom']
# Build a matrix with tha values of the stocks for the firms and the period considered:
dataArray_1B = u.getTablenp(dataset, titles, settlement, startDate, Names_1B)
# Compute the log returns of the stocks:
logReturns_1B = u.logReturnsComputation(dataArray_1B)
# Lambda for the weights in the Weighted Historical Simulation:
Lambda_1B = 0.95
# Vector of weights of the portfolio:
weights_1B = np.array([1/len(Names_1B)]*len(Names_1B))
# Notional of the portfolio:
portfolioValue_1B = 1
# To have daily values of the risk measure set:
delta_1B = 1
# Clock start:
start_time = time.time()
# Compute the Var and the ES via a Weighted Historical simulation approach:
[VaR_1B_WHS, ES_1B_WHS] = u.WHSMeasurements(logReturns_1B, alpha_1, Lambda_1B, weights_1B, portfolioValue_1B, delta_1B)
# Clock end and computational time:
end_time = time.time()
computational_time = end_time - start_time
# Plausibility check of the VaR:
VaR_1B_PC = u.plausibilityCheck(logReturns_1B, weights_1B, alpha_1, portfolioValue_1B, delta_1B)

print('\n\nQUESTION B:')
print("Equally weighted equity portfolio with Adidas, Airbus, BBVA, BMW and Deutsche Telekom")
print('\nVaR and the ES via Weighted Historical Simulation:')
print('VaR_WHS: ', VaR_1B_WHS,'     ES_WHS:', ES_1B_WHS)
print('\nPlausibility check of the VaR in the question B: ', VaR_1B_PC)
print('\nComputational time: ', computational_time)


# Question C:
# Vector of the firm Names:
Names_1C = titles.index[0:19].values
Names_1C = Names_1C[Names_1C != 'Adyen']
# Vector of weights of the portfolio:
weights_1C = np.array([1/len(Names_1C)]*len(Names_1C))
# Build a matrix with the values of the stocks for the firms and the period considered:
dataArray_1C = u.getTablenp(dataset, titles, settlement, startDate, Names_1C)
# Compute the log returns of the stocks:
logReturns_1C = u.logReturnsComputation(dataArray_1C)
# Compute the covariance matrix of logReturns for 1 year:
yearlyCovariance = np.cov(logReturns_1C.T)*256
# Compute the mean vector of logReturns computed for 1 year:
yearlyMeanReturns = -np.mean(logReturns_1C, axis=0)*256
# Time interval expressed in years in order to compute 10 daily values of the risk measure set:
H = 10/256
# Notional of the portfolio:
portfolioValue_1C = 1
# Vector of number of principal components:
numberOfPrincipalComponents = np.arange(5) + 1
# # Clock start:
start_time = time.time()
# Compute the Var and the ES via a Gaussian parametric approach:
result_1C = [u.PrincCompAnalysis(yearlyCovariance, yearlyMeanReturns, weights_1C, H, alpha_1, i, portfolioValue_1C) for i in numberOfPrincipalComponents]
# Clock end and computational time:
end_time = time.time()
computational_time = end_time - start_time

# Plausibility check of the VaR:
VaR_1C_PC = u.plausibilityCheck(logReturns_1C, weights_1C, alpha_1, portfolioValue_1C, H*256)

# Vector of number of principal components used for plotting:
AllComponents = np.arange(18)+1
# Plot the VaR and the ES computed via the PCA for different numbers of principal components
u.plotPCA(yearlyCovariance, yearlyMeanReturns, weights_1C, H, alpha_1, AllComponents, portfolioValue_1C)
# Plot the explained variance ratio of each principal component
u.plotPCA_Explain(yearlyCovariance)

# Print the results:
print('\n\nQUESTION C:')
print("Equally weighted equity portfolio with 18 companies")
print('\nVaR and the ES via Gaussian parametric approach (PCA):')
[print('Number of Components: ',i, '   VaR: ',result_1C[i-1][0],'   ES: ', result_1C[i-1][1]) for i in range(1,6)]
print('\nPlausibility check of the VaR in the question C: ', VaR_1C_PC)
print('\nComputational time: ', computational_time)



# ----- EXERCISE 2 ------

# Define today's date:
settlement_2 = dt.datetime.strptime('2017-01-16', '%Y-%m-%d').date()
# Define the star date:
startDate_2 = settlement_2 - dt.timedelta(days = 365*2+1)
# Vector of the firm name:
Names_2 = ['BMW']
# Build a matrix with the values of the stocks of the firm in the period considered:
dataArray_2 = u.getTablenp(dataset, titles, settlement_2, startDate_2, Names_2)
# Compute the log returns of the stocks:
logReturns_2 = u.logReturnsComputation(dataArray_2)
# set the random seed
random.seed(42)
# Number of simulations:
Nsim_2 = 100000
# Notional value of the stocks in the portfolio:
Notional_2 = 1186680
# Strike of the call options:
strike_2 = 25
# Volatility:
volatility_2 = 0.154
# Dividend yield:
dividend_2 = 0.031
# Fixed interest rate for the period:
rate_2 = 0.005
# Significance level:
alpha_2 = 0.95
# Lambda for the weights in the Weighted Historical Simulation:
Lambda_2 = 0.95
# Compute the  time interval in years for which compute the VaR:
riskMeasureTimeIntervalInYears_2 = FEL.yearfrac(pd.to_datetime('2017-01-16'),pd.to_datetime('2017-01-16')+pd.Timedelta(days=10),3)
# Compute the Time to Maturity in years for the Call;
timeToMaturityInYears_2 = FEL.yearfrac(pd.to_datetime('2017-01-16'),pd.to_datetime('2017-04-18'),3)
# Number of days for each year:
NumberOfDaysPerYears_2 = 365
Ticker2 = titles.loc[Names_2, 'Ticker'].values
stockPrice_2 = dataset.loc['2017-01-16'][Ticker2].values[0]
# Compute the number of shares:
numberOfShares_2 = Notional_2/stockPrice_2
# Number of Calls equal to the number of shares:
numberOfCalls_2 = numberOfShares_2
start_time = time.time()
# Compute the Var and the ES via a Gaussian parametric approach:
VaR_2_FMC = u.FullMonteCarloVaR(logReturns_2, numberOfShares_2, numberOfCalls_2, stockPrice_2, strike_2, rate_2, dividend_2,
volatility_2, timeToMaturityInYears_2, riskMeasureTimeIntervalInYears_2, alpha_2,NumberOfDaysPerYears_2, Lambda_2, Nsim_2)
# Clock end and computational time:
end_time = time.time()
computational_time_FMC = end_time - start_time
start_time = time.time()
# Compute the Var and the ES via a Gaussian parametric approach:
VaR_2_DN = u.DeltaNormalVaR(logReturns_2, numberOfShares_2, numberOfCalls_2, stockPrice_2, strike_2, rate_2, dividend_2,
volatility_2, timeToMaturityInYears_2, riskMeasureTimeIntervalInYears_2, alpha_2,NumberOfDaysPerYears_2, Lambda_2, Nsim_2)
# Clock end and computational time:
end_time = time.time()
computational_time_Delta = end_time - start_time

# Print the results:
print('\n\n----------- EXERCISE 2 -----------')
print('\nPortfolio with', numberOfShares_2, 'shares of BMW and',numberOfCalls_2, 'short Call Options')
print('\nVaR via Full Monte-Carlo approach:')
print('VaR_FMC: ', VaR_2_FMC)
print('\nVaR via Delta Normal approach:')
print('VaR_DN: ', VaR_2_DN[0])
print('\nVaR via Delta-Gamma approach:')
print('VaR_DN: ', VaR_2_DN[1])
print('\nComputational time for Full MC: ', computational_time_FMC)
print('Computational time for Delta Normal: ', computational_time_Delta)
