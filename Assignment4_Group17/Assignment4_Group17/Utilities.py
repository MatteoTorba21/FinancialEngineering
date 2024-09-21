import numpy as np
import scipy as sp
import pandas as pd
import datetime as dt
import matplotlib.pyplot as plt
from scipy.stats import kstest, t, shapiro, norm
from datetime import datetime
import FE_Library as FEL

def logReturnsComputation(returns):
    """ Function to compute the logReturns
        INPUT:     returns:     vector of returns
        OUTPUT:    logReturns:  vector of logReturns
        """
    # Compute returns from given data:
    Returns = returns[1:, :] / returns[0:-1, :]
    # Compute logReturns:
    logReturns = np.log(Returns)
    return logReturns

def AnalyticaltStudentMeasures(alpha, weights, portfolioValue, riskMeasureTimeIntervalInDay, returns, nu):
    """ Function to compute the VaR and the ES of a portfolio via Variance-covariance method t-student (with 4 degrees of freedom) parametric approach
        INPUT:     alpha:                             significance level
                   weights:                           vector of weights of the assets in the portfolio
                   portfolioValue:                    notional of the portfolio
                   riskMeasureTimeIntervalInDay:      time interval in days in which compute the VaR and ES
                   returns:                           vector of returns
                   nu:                                degrees of freedom of the t-Student distribution
        OUTPUT:    VaR:                               Value of Risk of the portfolio with Variance/Covariance method with t-Student quantile
                   ES:                                Expected Shortfall of the portfolio with Variance/Covariance method with t-Student quantile
        """
    # Compute mean vector:
    means_vector = np.mean(returns, axis=0)
    # Compute portfolio mean:
    ptf_mean = -np.dot(np.array(weights),means_vector.reshape(4,1))[0]
    # Compute variance matrix:
    variance_matrix = np.cov(returns.T)
    # Compute portfolio variance:
    ptf_covariance = np.dot(weights, np.dot(variance_matrix,weights.T))
    # Compute the standard VaR for a t-Student:
    VaR_std = sp.stats.t.ppf(alpha, nu)
    # Compute the standard ES for a t-Student:
    ES_std = ((nu + sp.stats.t.ppf(alpha, nu) ** 2) / (nu - 1)) * ((sp.stats.t.pdf(sp.stats.t.ppf(alpha, nu), df=nu)) / (1 - alpha))
    # Compute the VaR:
    VaR = (riskMeasureTimeIntervalInDay * ptf_mean + np.sqrt(riskMeasureTimeIntervalInDay * ptf_covariance) * VaR_std) * portfolioValue
    # Compute the ES:
    ES = (riskMeasureTimeIntervalInDay * ptf_mean + np.sqrt(riskMeasureTimeIntervalInDay * ptf_covariance) * ES_std) * portfolioValue
    return VaR, ES

def verifytstudentdistributionhypothesis(weights,returns,df):
    """ Function to check the t-student hypothesis for the distribution of the loss of the portfolio, and plot the loss against the PDF of a t-student distribution
        INPUT:     weights:             vector of weights of the assets in the portfolio
                   returns:             vector of log-returns of the assets in the portfolio
                   df:                  degrees of freedom of the t-Student
        """
    # Compute the loss of the linear portfolio (without scaling for the portfolio value):
    L = - (np.dot(weights, returns[1:, :].T))
    # Portfolio mean:
    ptf_mean = -np.dot(np.array(weights), np.mean(returns, axis=0).reshape(4, 1))[0]
    # Portfolio standard deviation:
    ptf_variance = np.sqrt(np.dot(weights, np.dot(np.cov(returns.T), weights.T)))
    # Calculate the t-statistic with 5 degrees of freedom:
    t_statistic, _ = sp.stats.ttest_1samp(L, ptf_mean)
    # Compute the p-value:
    p_value = sp.stats.t.sf(abs(t_statistic), df=5)
    print("t-statistic:", t_statistic)
    print("p-value:", p_value)
    # Interpretation:
    if p_value > 0.05:
        print("Fail to reject the null hypothesis. Data is consistent with t-distribution.")
    else:
        print("Reject the null hypothesis. Data is unlikely from the t-distribution.")
    # Graphical representation:
    # Histogram of the losses:
    plt.hist(L, bins=30, density=True, alpha=0.5, color='blue', label='Histogram of Losses')
    # Plot the probability density function (PDF) of the t-distribution:
    x = np.linspace(min(L), max(L), 1000)
    pdf = t.pdf(x, df, loc=ptf_mean, scale=ptf_variance)
    plt.plot(x, pdf, color='red', linestyle='--', label='t-Distribution PDF')
    # Add labels and legend:
    plt.xlabel('Losses')
    plt.ylabel('Probability Density')
    plt.title('Histogram of Losses vs. t-Distribution PDF')
    plt.legend()
    # Show plot:
    plt.show()
    return

def verifyGaussiandistributionhypothesis(weights, returns):
    """ Function to verify that the Gaussian distribution hypothesis are verified
        INPUT:     weights:           vector of weights of the portfolio
                   returns:           vector of log returns
        """
    # Compute the loss of the linear portfolio (without scaling for the portfolio value):
    L = - (np.dot(weights, returns[1:, :].T))
    # Shapiro test:
    statistic, pValue = shapiro(L)
    print("Statistic:", statistic)
    print("P-value:%.10f" % pValue)
    # Portfolio mean:
    ptf_mean = np.dot(np.array(weights), np.mean(returns, axis=0).reshape(18, 1))[0]
    # Portfolio variance:
    ptf_variance = np.sqrt(np.dot(weights, np.dot(np.cov(returns.T), weights.T)))
    # Plot the histogram of the losses:
    plt.hist(L, bins=30, density=True, alpha=0.5, color='blue', label='Histogram of Losses')
    # Plot the probability density function (PDF) of the Gaussian:
    x = np.linspace(min(L), max(L), 1000)
    pdf = sp.stats.norm.pdf(x, ptf_mean, ptf_variance)
    plt.plot(x, pdf, color='red', linestyle='--', label='Gaussian PDF')
    # Add labels and legend:
    plt.xlabel('Losses')
    plt.ylabel('Probability Density')
    plt.title('Histogram of Losses vs. t-Distribution PDF')
    plt.legend()
    # Show plot:
    plt.show()
    return

def getTablenp (dataset, titles, settlement, startDate, Names):
    """ Compute Stock Values as numpy Matrix, keeping the ones from "startDate" to "settlement" only for "Names" titles
        INPUT:      dataset:            table with all the stock values
                    titles:             table with titles names and tickers
                    settlement:         settlement date
                    startDate:          first date to consider
                    Names:              names of the titles to keep
        OUTPUT:     dataPar.values:     numpy array of the wanted stock values
        """
    # Get the tickers needed:
    valsNames = titles.loc[Names, 'Ticker']
    # Get the stock values only for the dates needed:
    dataCut = dataset.loc[startDate:settlement, :]
    # Get only the stock values as an array:
    dataPar = dataCut[valsNames.values]
    return dataPar.values

def HSMeasurements(returns, alpha, weights, portfolioValue, riskMeasureTimeIntervalInDay):
    """ Historical Simulation approach
        INPUT:     returns:                          vector of logReturns
                   alpha:                            significance level
                   weights:                          vector of weights
                   portfolioValue:                   portfolio Notional
                   riskMeasureTimeIntervalInDay:     time interval in days for which compute the VaR and the ES
        OUTPUT:    VaR:                              Value at Risk computed with the Historical Simulation approach
                   ES:                               Expected Shortfall computed with the Historical Simulation approach
        """
    # Compute the losses' array under the frozen portfolio assumption:
    L = -portfolioValue*(np.dot(weights,returns[1:, :].T))
    # Length (in days) of the time window considered:
    n = len(returns)
    # Sort the losses' array to make it decreasing:
    LOrdered = sorted(L, reverse = True)
    # Get the biggest integer before n*(1-alpha):
    istar = int(np.floor(n * (1 - alpha)))
    # Get the VaR and rescale it (the istar-th element of the ordered vector is in position istar-1, since the first element is in position 0):
    VaR = np.sqrt(riskMeasureTimeIntervalInDay)*LOrdered[istar-1]
    # Compute the ES and rescale it:
    ES = np.sqrt(riskMeasureTimeIntervalInDay)*np.mean(LOrdered[:istar])
    return VaR, ES

def bootstrapStatistical(numberOfSamplesToBootstrap, returns):
    """ Generate "numberOfSamplesToBootstrap" samples randomly picked from the log returns
        INPUT:      numberOfSamplesToBootstrap:         number of samples that we need
                    returns:                            vector of log returns
        OUTPUT:     samples:                            the wanted samples' vector
    """
    # Compute all the random indices needed:
    random = np.random.randint(1, 256*5, numberOfSamplesToBootstrap)
    # Get the samples of random logReturns:
    samples = returns[random, :]
    return samples

def WHSMeasurements(returns, alpha, Lambda, weights, portfolioValue, riskMeasureTimeIntervalInDay):
    """ Weighted Historical Simulation approach
        INPUT:      returns:                         vector of logReturns
                    alpha:                           significance level
                    Lambda:                          value of lambda used to compute the weights
                    weights:                         vector of weights
                    portfolioValue:                  portfolio Notional
                    riskMeasureTimeIntervalInDay:    time interval in days for which compute the VaR and the ES
        OUTPUT:     VaR:                             Value at Risk computed with the Weighted Historical Simulation approach
                    ES:                              Expected Shortfall computed with the Weighted Historical Simulation approach
        """
    # Length of the logReturn vector:
    n = len(returns)
    # Compute the normalization factor:
    C = (1-Lambda)/(1-Lambda**n)
    # Compute the weights' array, decreasing in the past:
    w = C * Lambda ** np.arange(n-1,-1,-1)
    # Compute the losses' array:
    L = -portfolioValue * np.dot(returns, weights)
    # Indices of the descending losses (np.argsrot applied to -L to have a descending sequence):
    desc_indices = np.argsort(-L)
    # Loss sequence ordered in descendent way:
    desc_L = L[desc_indices]
    # Weight sequence ordered in descendent way (in correspondence with the losses):
    desc_w = w[desc_indices]
    # Vector in which the i^{th} element is the sum of the weights of the biggest i losses:
    sum_w = np.cumsum(desc_w)
    # Position of the Loss vector corresponding to the Loss to be considered to compute the VaR:
    i_star = np.searchsorted(sum_w, 1-alpha, side='right')-1
    # Compute the VaR:
    VaR = np.sqrt(riskMeasureTimeIntervalInDay)*desc_L[i_star]
    # Compute the ES:
    ES = np.sqrt(riskMeasureTimeIntervalInDay)*np.sum(desc_w[:i_star+1]*desc_L[:i_star+1])/np.sum(desc_w[:i_star+1])
    return VaR,ES

def PrincCompAnalysis(yearlyCovariance, yearlyMeanReturns, weights, H, alpha, numberOfPrincipalComponents, portfolioValue):
    """ Gaussian parametric approach (PCA)
        INPUT:     yearlyCovariance:               covariance matrix of logReturns for 1 year
                   yearlyMeanReturns:              mean vector of logReturns for 1 year
                   weights:                        vector of weights
                   H:                              time interval in years for which compute the VaR and the ES
                   alpha:                          significance level
                   numberOfPrincipalComponents:    vector with the number of principal components
                   portfolioValue:                 portfolio Notional
        OUTPUT:    VaR:                            Value at Risk computed with the Gaussian parametric approach
                   ES:                             Expected Shortfall computed with the Gaussian parametric approach
        """
    # Compute eigenvalues and the matrix of eigenvectors:
    eigenval, eigenvecMat = np.linalg.eig(yearlyCovariance)
    # Order the eigenvalues matrix:
    indeces = np.argsort(-eigenval)
    # Order the eigenvectors matrix in the same way of the eigenvalues' one:
    eigenvecMat = eigenvecMat[:,indeces]
    # Get the first n indices of the decreasing array of eigenvalues, where n is the number of principal components:
    indeces_to_consider = np.argsort(-eigenval)[0:numberOfPrincipalComponents]
    # Create the diagonal matrix composed by the n smallest eigenvalues, where n is the number of principal components:
    Lambda = np.diag(eigenval[indeces_to_consider])
    # Compute the weights and the mean of the portfolio projected on the principal components:
    weightsHat = np.dot(eigenvecMat.T, weights)
    MeanHat = np.dot(eigenvecMat.T, yearlyMeanReturns)
    # Compute the mean of the reduced portfolio:
    MeanReduced = np.dot(MeanHat[indeces_to_consider].T, weightsHat[indeces_to_consider])
    # Compute the volatility of the reduced portfolio:
    SigmaReduced = np.dot(weightsHat[indeces_to_consider], np.dot(Lambda, weightsHat[indeces_to_consider]).T)
    # Compute the standardized VaR considering a gaussian distribution:
    VaR_Std = sp.stats.norm.ppf(alpha)
    # Compute the standardized ES considering a gaussian distribution:
    ES_Std = sp.stats.norm.pdf(sp.stats.norm.ppf(alpha))/(1 - alpha)
    # Compute the VaR for the reduced portfolio:
    VaR = (H*MeanReduced + np.sqrt(H*SigmaReduced)*VaR_Std)*portfolioValue
    # Compute the ES for the reduced portfolio:
    ES = (H*MeanReduced + np.sqrt(H*SigmaReduced)*ES_Std)*portfolioValue
    return VaR, ES

def plotPCA_Explain (yearlyCovariance):
    """ Plot cumulated explained variance ratios
        INPUT:      yearlyCovariance:           covariance matrix of logReturns for 1 year
    """
    # Compute eigenvalues and eigenvectors matrix:
    eigenval, eigenvecMat = np.linalg.eig(yearlyCovariance)
    # Order the eigenvalues indices:
    indices = np.argsort(-eigenval)
    # Compute the Explained Variance Ratios:
    explained_variance_ratio = eigenval[indices] / np.sum(eigenval[indices])
    # Compute the cumulated Explained Variance Ratios:
    cumulated_explained_variance_ratio = np.cumsum(explained_variance_ratio)
    # Plot cumulated explained variance ratios:
    plt.figure()
    plt.bar(range(1, len(cumulated_explained_variance_ratio) + 1), cumulated_explained_variance_ratio, alpha=0.5, align='center')
    plt.xticks(range(1, len(explained_variance_ratio) + 1))
    plt.grid(True)
    # Add labels and legend
    plt.xlabel('Principal Component')
    plt.ylabel('Cumulative Explained Variance Ratio')
    plt.title('Explained Variance Ratio')
    # Show plot:
    plt.show()
    return

def plotPCA(yearlyCovariance, yearlyMeanReturns, weights, H, alpha, numberOfPrincipalComponents, portfolioValue):
    """ Plot of the VaR computed via the PCA (Gaussian parametric approach) for different number of principal components
        INPUT:     yearlyCovariance:               covariance matrix of logReturns for 1 year
                   yearlyMeanReturns:              mean vector of logReturns for 1 year
                   weights:                        vector of weights
                   H:                              time interval in years for which compute the VaR and the ES
                   alpha:                          significance level
                   numberOfPrincipalComponents:    vector with the number of principal components
                   portfolioValue:                 portfolio Notional
        """
    # Compute the VaR and the ES via PCA for different numbers of principal components:
    PCAs = [PrincCompAnalysis(yearlyCovariance, yearlyMeanReturns, weights, H, alpha, i, portfolioValue) for
              i in numberOfPrincipalComponents]
    # Save the results:
    VaRs = [x for x,y in PCAs]
    ESs = [y for x,y in PCAs]
    # Compute the gaussian VaR and ES:
    VaR, ES = AnalyticalGaussianVaRandES(alpha, yearlyCovariance, yearlyMeanReturns, H, weights)
    # Plot the VaR results for the PCA:
    plt.figure()
    plt.plot(np.arange(1, len(PCAs) + 1), VaRs, linestyle='--', marker='o',label='PCA')
    plt.plot(np.arange(1, len(PCAs) + 1),VaR*np.ones((len(PCAs))), color='r', linestyle='--', label='Gaussian')
    # Set the x axis:
    plt.xticks(range(1, len(PCAs) + 1))
    plt.grid(True)
    # Add labels and legend:
    plt.xlabel('Number of principal components')
    plt.ylabel('VaR')
    plt.title('PCA approximation for VaR')
    plt.legend()
    # Show plot:
    plt.show()
    # Plot the ES results for the PCA:
    plt.figure()
    plt.plot(np.arange(1, len(PCAs) + 1), ESs, linestyle='--', marker='o', label='PCA')
    plt.plot(np.arange(1, len(PCAs) + 1), ES * np.ones((len(PCAs))), color='r', linestyle='--', label='Gaussian')
    plt.xticks(range(1, len(PCAs) + 1))
    plt.grid(True)
    # Add labels and legend:
    plt.xlabel('Number of principal components')
    plt.ylabel('ES')
    plt.title('PCA approximation for ES')
    plt.legend()
    # Show plot:
    plt.show()
    return

def AnalyticalGaussianVaRandES(alpha, yearlyCovariance, yearlyMeanReturns, riskMeasureTimeIntervalInYears, weights):
    """ Function to compute the VaR and the ES of a portfolio via Variance-covariance method t-student parametric approach
        INPUT:      alpha:                             significance level
                    weights:                           vector of weights of the assets in the portfolio
                    portfolioValue:                    notional of the portfolio
                    riskMeasureTimeIntervalInYears:    time interval in years in which compute the VaR and ES
        OUTPUT:     VaR:                               Value at risk of the portfolio computed with Variance/Covariance method with Gaussian quantile
                    ES:                                Expected Shortfall of the portfolio computed with Variance/Covariance method with Gaussian quantile
        """
    # Compute mean vector:
    Mean_ptf = np.dot(yearlyMeanReturns.T, weights)
    # Compute portfolio variance:
    Sigma_ptf = np.dot(weights, np.dot(yearlyCovariance, weights).T)
    # Compute the standard VaR for Gaussian:
    VaR_std = sp.stats.norm.ppf(alpha)
    # Compute the VaR:
    VaR = riskMeasureTimeIntervalInYears * Mean_ptf + np.sqrt(riskMeasureTimeIntervalInYears * Sigma_ptf) * VaR_std
    # Compute the ES:
    ES = riskMeasureTimeIntervalInYears * Mean_ptf + (np.sqrt(riskMeasureTimeIntervalInYears * Sigma_ptf) * norm.pdf(norm.ppf(alpha)) / (1 - alpha))
    return VaR, ES

def plausibilityCheck(returns, portfolioWeights, alpha, portfolioValue, riskMeasureTimeIntervalInDay):
    """ Check the order of magnitude of the Value at Risk
        INPUT:     returns:                          vector of logReturns
                   portfolioWeights:                 vector of portfolio weights
                   alpha:                            significance level
                   portfolioValue:                   portfolio Notional
                   riskMeasureTimeIntervalInDay:     time interval in days for which compute the VaR and the ES
        OUTPUT:    VaR:                              Value at Risk computed with the Plausibility check
        """
    # Matrices of losses:
    L = -returns
    # Compute correlation matrix:
    C = np.corrcoef(returns.T)
    # Compute lower and upper bounds:
    l, u = [np.percentile(L, 100*(star), axis = 0) for star in [alpha, 1-alpha]]
    # Compute the Signed-VaR:
    sVaR = portfolioWeights*(abs(l) + abs(u)) /2
    # Compute the VaR:
    VaR = np.sqrt(riskMeasureTimeIntervalInDay) * portfolioValue * np.sqrt(np.dot(np.dot(sVaR, C), sVaR))
    return VaR

def CallPriceComputation(stockPrice,strike,rate,dividend,volatility,timeToMaturityInYears):
    """ Computation of a Call Option price with dividends
        INPUT:      stockPrice:                value of the stock
                    strike:                    strike price
                    rate:                      fixed rate
                    dividend:                  dividend yield
                    volatility:                volatility
                    timeToMaturityInYears:     time to maturity computed in years
        OUTPUT:     CallPrice:                 price of a Call option
        """
    # Compute d1:
    d1 = (np.log(stockPrice/strike)+(rate-dividend+volatility**2/2)*timeToMaturityInYears)/(volatility*np.sqrt(timeToMaturityInYears))
    # Compute d2:
    d2 = d1-volatility*np.sqrt(timeToMaturityInYears)
    # Compute the price of the Call with dividends:
    CallPrice = stockPrice*np.exp(-dividend*timeToMaturityInYears)*sp.stats.norm.cdf(d1)-strike * np.exp(-rate*timeToMaturityInYears)*sp.stats.norm.cdf(d2)
    return CallPrice

def CallPriceComputationGeman(stockPrice,strike,discount,dividend,volatility,timeToMaturityInYears):
    """ Computation of a Call Option price with dividends via Geman-El Karroui-Rochet model
        INPUT:      stockPrice:                value of the stock
                    strike:                    strike price
                    discount:                  discount factor
                    dividend:                  dividend yield
                    volatility:                volatility
                    timeToMaturityInYears:     time to maturity computed in years
        OUTPUT:     CallPrice:                 price of a Call option
        """
    # Compute d1:
    d1 = (np.log(stockPrice/(strike*discount))+(dividend+volatility**2/2)*timeToMaturityInYears)/(volatility*np.sqrt(timeToMaturityInYears))
    # Compute d2:
    d2 = d1-volatility*np.sqrt(timeToMaturityInYears)
    # Compute the price of the Call with dividends via Geman-El Karroui-Rochet:
    CallPrice = stockPrice*np.exp(-dividend*timeToMaturityInYears)*sp.stats.norm.cdf(d1)-strike * discount*sp.stats.norm.cdf(d2)
    return CallPrice

def deltaCallComputation(stockPrice,strike,rate,dividend,volatility,timeToMaturityInYears):
    """ Computation of the Delta of a Call Option
        INPUT:      stockPrice:                value of the stock
                    strike:                    strike price
                    rate:                      fixed rate
                    dividend:                  dividend yield
                    volatility:                volatility
                    timeToMaturityInYears:     time to maturity computed in years
        OUTPUT:     CallPrice:                 price of a Call option
        """
    # Compute d1:
    d1 = (np.log(stockPrice/strike)+(rate-dividend+volatility**2/2)*timeToMaturityInYears)/(volatility*np.sqrt(timeToMaturityInYears))
    # Compute the delta of a Call option:
    deltaCall = np.exp(-dividend*timeToMaturityInYears)*sp.stats.norm.cdf(d1)
    return deltaCall

def gammaCallComputation(stockPrice,strike,rate,dividend,volatility,timeToMaturityInYears):
    """ Computation of the Gamma of a Call Option
        INPUT:     stockPrice:                value of the stock
                   strike:                    strike price
                   rate:                      fixed rate
                   dividend:                  dividend yield
                   volatility:                volatility
                   timeToMaturityInYears:     Time to Maturity computed in years
        OUTPUT:    CallPrice:                 price of a Call option
        """
    # Compute d1:
    d1 = (np.log(stockPrice/strike)+(rate-dividend+volatility**2/2)*timeToMaturityInYears)/(volatility*np.sqrt(timeToMaturityInYears))
    # Compute the gamma of a Call option:
    gammaCall = np.exp(-dividend*timeToMaturityInYears)*sp.stats.norm.pdf(d1)/(stockPrice*volatility*np.sqrt(timeToMaturityInYears))
    return gammaCall

def FullMonteCarloVaR(logReturns, numberOfShares, numberOfCalls, stockPrice, strike, rate, dividend,
                      volatility, timeToMaturityInYears, riskMeasureTimeIntervalInYears, alpha, NumberOfDaysPerYears, Lambda, Nsim):
    """ Full Monte-Carlo method using Weighted Historical Simulation approach
        INPUT:     logReturns:                        vector of logReturns
                   numberOfShares:                    number of shares of the stock
                   numberOfCalls:                     number of Calls
                   stockPrice:                        value of the stock
                   strike:                            strike price
                   rate:                              fixed rate
                   dividend:                          dividend yield
                   volatility:                        volatility
                   timeToMaturityInYears:             Time to Maturity computed in years
                   riskMeasureTimeIntervalInYears:    time interval in years for which compute the VaR
                   alpha:                             significance leve
                   NumberOfDaysPerYears               number of days in each year
                   Lambda:                            value of lambda used to compute the weights
                   Nsim:                              number of simulations
        OUTPUT:    VaR:                               Value at Risk computed via Full Monte-Carlo method using Weighted Historical Simulations
        """
    # Length of the logReturn vector:
    n = len(logReturns)
    # Compute the normalization factor:
    C = (1 - Lambda) / (1 - Lambda**n)
    # Compute the weights' array, decreasing in the past:
    w = C * Lambda ** np.arange(n-1,-1,-1)
    # Increasing array of indices, starting from 0, as long as the logReturns:
    indices = np.arange(len(logReturns))
    # Sampling using weights as probabilities:
    random = np.random.choice(indices, size=(Nsim,int(riskMeasureTimeIntervalInYears*NumberOfDaysPerYears)), p = w)
    # Get the samples:
    samples = np.empty_like(random, dtype=logReturns.dtype)
    # Initialization of the vector of weights:
    weights = np.zeros(random.shape)
    # Iterate over each element of random and assign values to samples and weights:
    for i in range(random.shape[0]):
        for j in range(random.shape[1]):
            samples[i, j] = logReturns[random[i, j]]
            weights[i][j] = w[random[i][j]]
    # Select elements from logReturns using the randomly generated indices in random:
    samples = np.sum(samples, axis=1)
    # Select elements from w using the same randomly generated indices in random:
    weights = np.mean(weights, axis=1)
    # Normalize the weights in order to sum up to 1:
    weights = weights/weights.sum()
    # Compute the stock price in t + riskMeasureTimeIntervalInDays:
    stockPrice_deltat_simulation = stockPrice * np.exp(samples)
    # Compute the Call price in t:
    callPrice = CallPriceComputation(stockPrice,strike,rate,dividend,volatility,timeToMaturityInYears)
    # Compute the Call price in t + riskMeasureTimeIntervalInDays:
    callPrice_deltat_simulation = np.array([CallPriceComputation(s10d,strike,rate,dividend,volatility,timeToMaturityInYears-riskMeasureTimeIntervalInYears) for s10d in stockPrice_deltat_simulation])
    # Compute the loss of the Call of the portfolio:
    LCall = numberOfCalls * (callPrice_deltat_simulation - callPrice)
    # Compute the loss of the stock of the portfolio:
    LStock = numberOfShares * (stockPrice - stockPrice_deltat_simulation)
    # Compute the loss of the portfolio:
    L = np.array(LCall + LStock)
    # Indices of the descending losses (np.argsrot applied to -L to have a descending sequence):
    desc_indices = np.argsort(-L.flatten())
    # Loss sequence ordered in descendent way:
    desc_L = L[desc_indices]
    # Weight sequence ordered in descendent way (in correspondence with the losses):
    desc_w = weights[desc_indices]
    # Vector in which the i^{th} element is the sum of the weights of the biggest i losses:
    sum_w = np.cumsum(desc_w)
    # Position of the Loss vector corresponding to the Loss to be considered to compute the VaR:
    i_star = np.searchsorted(sum_w, 1-alpha, side='right') - 1
    # Compute the VaR:
    VaR = desc_L[i_star]
    return VaR

def DeltaNormalVaR(logReturns, numberOfShares, numberOfCalls, stockPrice, strike, rate, dividend,
volatility, timeToMaturityInYears, riskMeasureTimeIntervalInYears, alpha, NumberOfDaysPerYears, Lambda, Nsim):
    """ Delta Normal method using Weighted Historical Simulation approach
        INPUT:    logReturns:                       vector of logReturns
                  numberOfShares:                   number of shares of the stock
                  numberOfCalls:                    number of Calls
                  stockPrice:                       value of the stock
                  strike:                           strike price
                  rate:                             fixed rate
                  dividend:                         dividend yield
                  volatility:                       volatility
                  timeToMaturityInYears:            Time to Maturity computed in years
                  riskMeasureTimeIntervalInYears:   time interval in years for which compute the VaR
                  alpha:                            significance level
                  NumberOfDaysPerYears:             number of days in each year
                  Lambda:                           value of lambda used to compute the weights
                  Nsim:                             number of simulations
        OUTPUT:   VaR: Value at Risk computed via Full Monte-Carlo method using Weighted Historical Simulations
        """
    # Length of the logReturn vector:
    n = len(logReturns)
    # Compute the normalization factor:
    C = (1 - Lambda) / (1 - Lambda ** n)
    # Compute the weights' array, decreasing in the past:
    w = C * Lambda ** np.arange(n-1,-1,-1)
    # Increasing array of indices, starting from 0, as long as the logReturns:
    indices = np.arange(len(logReturns))
    # Generate an array of random integers between 1 and the length of logReturns (exclusive) with a size of Nsim:
    random = np.random.choice(indices, size=(Nsim,int(riskMeasureTimeIntervalInYears*NumberOfDaysPerYears)), p = w)
    # Get the samples:
    samples = np.empty_like(random, dtype=logReturns.dtype)
    # Initialization of the vector of weights:
    weights = np.zeros(random.shape)
    # Iterate over each element of random and assign values to samples and weights:
    for i in range(random.shape[0]):
        for j in range(random.shape[1]):
            samples[i, j] = logReturns[random[i, j]]
            weights[i][j] = w[random[i][j]]
    # Select elements from logReturns using the randomly generated indices in random:
    samples = np.sum(samples, axis=1)
    # Select elements from w using the same randomly generated indices in random:
    weights = np.mean(weights, axis=1)
    # Normalize the weights in order to sum up to 1:
    weights = weights/weights.sum()
    # Calculate the delta of the Call option:
    deltaCall = deltaCallComputation(stockPrice,strike,rate,dividend,volatility,timeToMaturityInYears)
    # Compute the loss of the Call of the portfolio:
    LCallDelta = deltaCall * stockPrice * samples * numberOfCalls
    # Compute the loss of the stock of the portfolio:
    LStockNeutral = - samples * stockPrice * numberOfShares
    # Compute the loss of the portfolio considering the Delta:
    L = np.array(LCallDelta + LStockNeutral)
    # Indices of the descending losses (np.argsort is applied to -L in order to have a descending sequence):
    desc_indices = np.argsort(-L.flatten())
    # Loss sequence ordered in descendent way:
    desc_L = L[desc_indices]
    # Weight sequence ordered in descendent way (in correspondence with the losses):
    desc_w = weights[desc_indices]
    # Vector in which the i^{th} element is the sum of the weights of the biggest i losses:
    sum_w = np.cumsum(desc_w)
    # Position of the Loss vector corresponding to the Loss to be considered to compute the VaR:
    i_star = np.searchsorted(sum_w, 1-alpha, side='right') - 1
    # Compute the VaR considering the Delta:
    VaR_delta = desc_L[i_star]
    # Calculate the Gamma of the call option:
    gammaCall = gammaCallComputation(stockPrice, strike, rate, dividend, volatility, timeToMaturityInYears)
    # Compute the loss of the call of the portfolio considering the Gamma:
    LCallGamma = (deltaCall * stockPrice * samples + 1 / 2 * gammaCall * (stockPrice * samples) ** 2) * numberOfCalls
    # Compute the loss of the portfolio considering the Gamma:
    L_2 = np.array(LCallGamma + LStockNeutral)
    # Indices of the descending losses (np.argsort applied to -L to have a descending sequence):
    desc_indices2 = np.argsort(-L_2.flatten())
    # Loss sequence ordered in descendent way:
    desc_L2 = L_2[desc_indices2]
    # Weight sequence ordered in descendent way (in correspondence with the losses):
    desc_w2 = weights[desc_indices2]
    # Vector in which the i^{th} element is the sum of the weights of the biggest i losses:
    sum_w2 = np.cumsum(desc_w2)
    # Position of the Loss vector corresponding to the Loss to be considered to compute the VaR:
    i_star = np.searchsorted(sum_w2, 1 - alpha, side='right') - 1
    # Compute the VaR  considering the Gamma:
    VaR_gamma = desc_L2[i_star]
    return VaR_delta, VaR_gamma

def interpolate_rates(data, new_date):
  """ Performs linear interpolation on Discount and ZeroRate to estimate Discount, ZeroRate, Rate, and Return Rate for a new date
      INPUT:   data:         a pandas DataFrame containing columns 'Date', 'Discount', and 'ZeroRate'
               new_date:     the date for which to estimate the rates (format: YYYY-MM-DD)
      OUTPUT: 
           A dictionary containing interpolated values for Discount, ZeroRate, Rate, and Return Rate
      """
  # Convert 'Date' column to datetime format:
  data['Date'] = pd.to_datetime(data['Date'])
  # Sort data by date:
  data = data.sort_values(by='Date')
  # Find bracketing dates for the new date:
  bracketing_dates = data[data['Date'].dt.strftime('%Y-%m-%d') <= new_date].tail(2)
  if len(bracketing_dates) < 2:
    raise ValueError("New date is outside the range of data")
  # Extract data for bracketing dates:
  earlier_date, later_date = bracketing_dates['Date'].tolist()
  discount_earlier, discount_later = bracketing_dates['Discount'].tolist()
  zerorate_earlier, zerorate_later = bracketing_dates['ZeroRate'].tolist()
  # Calculate weights:
  time_diff_total = (later_date - earlier_date).days
  weight_earlier = (later_date - pd.to_datetime(new_date)).days / time_diff_total
  weight_later = 1 - weight_earlier
  # Perform linear interpolation for ZeroRate:
  zerorate_new = weight_earlier * zerorate_earlier + weight_later * zerorate_later
  # Calculate Discount from interpolated ZeroRate (exponential compounding):
  # Replace '2008-02-19' with your actual reference date if needed:
  ref_date = pd.to_datetime('2008-02-19')
  time_diff_new = (pd.to_datetime(new_date) - ref_date).days / 365
  discount_new = np.exp(-zerorate_new * time_diff_new)
  return {  "Discount": discount_new,
            "ZeroRate": zerorate_new}

def DiscountsAtOptionDates(optionDates,bootstrap):
    """ Compute discounts for optionDates
        INPUT:      optionDates:        vector of dates for which to compute discounts
        OUTPUT:     discounts:          vector of discount factors
    """
    # Generate list of discount factors interpolating rates for each date in optionDates:
    discounts = [interpolate_rates(bootstrap.copy(), date)['Discount'] for date in optionDates]
    # Insert discount factor 1 at beginning:
    discounts.insert(0, 1)
    # Convert to numpy array:
    discounts = np.array(discounts)
    return discounts



