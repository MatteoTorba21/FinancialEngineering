function [s,CI,P] = CDSspread(dates,discounts,TTM,R,datesCDS,intensities,N)
% Function to compute the spread of a Credit Default Swap
%
% INPUTS:
% dates:         Vector of dates obtained from the bootstrap
% discounts:     Vector of discounts obtained from the bootstrap
% TTM:           Time to maturity expressed in years
% R:             Recovery rate of the first obligor
% datesCDS:      Vector of dates of the Credit Default Swap
% intensities:   Vector of the intensities of the obligor
% N:             Number of Monte-Carlo simulations
%
% OUTPUTS:
% s:             Spread of the annual fee
% CI:            Confidence interval
% P:             Default percentage
%

% Set the seed to use the same random numbers at each simulation:
rng(321)           
% Convenction 30/360
SwapDayCount = 6;
% Yearfrac convenction: ACT/365
Act365 = 3;

% Generation of gaussian random variables with mean 0 and covariance function Sigma:
x = randn(1,N);
% Uniform distributed variables:
u = normcdf(x);
% Dates od the CDS up to Maturity:
datesCDS = datesCDS(1:TTM);
% Initialization of the default time vector:
tau=zeros(N,1);

% Time to reset dates:
deltasStart = yearfrac(dates(1), datesCDS, Act365);

% Computation of the default times by reversing the probability:
for i = 1:N
    survFunct1 = @(t) (t<TTM)*exp(-min(1,max(0,t-[0;deltasStart(1:end-1)]))'*intensities(1:TTM))-u(1,i);
    tau(i) = fzero(survFunct1,0);
end

% Remove the cases in which the default did not occur:
idx = find(tau<TTM-1e-6); % (TTM-1e-6 chose to avoid errors in the approximation)
result = tau(idx);

% Discounts and deltas in the CDS dates
DFCDS = InterpDFviaRates(dates, discounts, datesCDS);
deltas = yearfrac([dates(1); datesCDS(1:end-1)], datesCDS, SwapDayCount);

% Default percentage:
P = length(result)/N;

% Contingent and fee leg initialization:
contingentLeg = zeros(N,1);
% iIn fee leg sum the payments when no default occurs up to FTD maturity:
feeLeg = (DFCDS'*deltas) * ones(N,1);

% Cases in which defaults occur:
for ii=1:length(result)
    % Index of first reset date after default occurs:
    ContingentResetDate = find(deltasStart > result(ii),1);
    % Discount in default event:
    discount_contingentleg = InterpDFviaRates(dates, discounts, dates(1)+365*result(ii));
    % Vector of discounts of reset dates before default:
    discounts_feeleg = DFCDS(1:ContingentResetDate-1);
    % Vector of deltas of reset dates before default:
    deltas_feeleg = deltas(1:ContingentResetDate-1);
    % Accrual delta calculation:
    if ContingentResetDate == 1
        accrualDelta = result(ii);
    else
        accrualDelta = result(ii) - sum(deltas(1:ContingentResetDate-1));
    end
    % Contingent and fee legs update:
    contingentLeg(ii) =(1-R)*discount_contingentleg;
    feeLeg(ii) = (discounts_feeleg'*deltas_feeleg) + accrualDelta*discount_contingentleg;
end

% Equivalent spread:
s = sum(contingentLeg)/sum(feeLeg);

% Fieller confidence interval calculation by using delta method as in the Franz paper:
alpha = 0.05;
t = tinv(1 - alpha/2, N-1);
VC = cov(contingentLeg,feeLeg);

meanf = mean(feeLeg);
meanc = mean(contingentLeg);
sc = var(contingentLeg)/N;
sf = var(feeLeg)/N;
cv = VC(1,2)/N;

LL = ((meanf*meanc-t^2*cv) - sqrt( (meanf*meanc-t^2*cv)^2 -...
    (meanf^2 - sf * t^2) * (meanc^2 - sc * t^2)) )/ (meanf^2 - t^2*sf);
UL = ((meanf*meanc-t^2*cv) + sqrt( (meanf*meanc-t^2*cv)^2 -...
    (meanf^2 - sf * t^2) * (meanc^2 - sc * t^2)) )/ (meanf^2 - t^2*sf);

CI = [LL UL];

end % function CDSspread