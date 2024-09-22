function [s, CI, P] = firstToDefault2NamesNoAccrual(dates,discounts,TTM,R1,R2,datesCDS,intensities1,intensities2,rho,N)
% Function to compute the price of a First to Default derivative written
% considering two obligors by considering the accrual term
%
% INPUTS:
% dates:            Vector of dates obtained from the bootstrap
% discounts:        Vector of discounts obtained from the bootstrap
% TTM:              Time to maturity expressed in years
% R1:               Recovery rate of the first obligor
% R2:               Recovery rate of the second obligor
% datesCDS:         Vector of dates of the Credit Default Swap
% intensities1:     Vector of the approximated intensities of the first obligor
% intensities2:     Vector of the approximated intensities of the second obligor
% rho:              Default correlation of the two obligors
% N:                Number of Monte-Carlo simulations
%
% OUTPUTS:
% s:                Spread of the annual fee
% CI:               Confidence interval
% P:                Default matrix:
%                      P(1,1): no default percentage
%                      P(1,2): only name 1 defaults percentage
%                      P(2:1): only name 2 defaults percentage
%                      P(2,2): both defaults percentage
%

warning('off')
% Set the seed to use the same random numbers at each simulation:
rng(321)
% Convenction 30/360
SwapDayCount = 6;
% Yearfrac convenction: ACT/365
Act365 = 3;

% Correlation matrix:
Sigma = [1 rho; rho 1];
% Lower triangular matrix such that A*A'=Sigma:
A = chol(Sigma)';
% Generation of gaussian random variables with mean 0 and covariance function Sigma:
x = A*randn(2,N);
% Uniform distributed variables:
u = normcdf(x);
% Dates od the CDS up to Maturity:
datesCDS = datesCDS(1:TTM);
% Initialization of the default time vectors:
tau1=zeros(N,1);
tau2=zeros(N,1);

% Time to reset dates:
deltasStart = yearfrac(dates(1), datesCDS, Act365);

% Computation of the default times by reversing the probability:
for i = 1:N
    survFunct1 = @(t)  (t<TTM)*exp(-min(1,max(0,t-[0;deltasStart(1:end-1)]))'*intensities1(1:TTM))-u(1,i);
    tau1(i) = fzero(survFunct1,0);
    survFunct2 = @(t)  (t<TTM)*exp(-min(1,max(0,t-[0;deltasStart(1:end-1)]))'*intensities2(1:TTM))-u(2,i);
    tau2(i) = fzero(survFunct2,0);
end

% Calculate the minimum between tau1 and tau2
min_values = min([tau1, tau2], [], 2);

% Initialize the Nx2 matrix to store the minimum times and a flag equal to
% 1 or 2 which indicates the first one defaulted:
N = length(tau1);
result_matrix = zeros(N, 2);

% Assign values to the matrix:
result_matrix(:, 1) = min_values; % Time of the first to default
result_matrix(:, 2) = 1 + (min(tau1, tau2) == tau2); % Flag to indicate which is the first to default

% Remove the cases in which the default did not occur:
idx = find(result_matrix(:,1)<TTM-1e-6); % (TTM-1e-6 chose to avoid errors in the approximation)
result_matrix = result_matrix(idx,:);

% Discounts and deltas in the CDS dates:
DFCDS = InterpDFviaRates(dates, discounts, datesCDS);
deltas = yearfrac([dates(1); datesCDS(1:end-1)], datesCDS, SwapDayCount);

% Vector of the recovery rates:
R = [R1; R2];

% P calculations:
% no default percentage
percentageOfNoDefault = 1 - length(result_matrix(:,1))/N;
% only name 1 defaults percentage
percentageOfOnlyIssuer1Default = sum((tau1<4-1e-6).*(tau2>4-1e-6))/N;
% only name 2 defaults percentage
percentageOfOnlyIssuer2Default = sum((tau1>4-1e-6).*(tau2<4-1e-6))/N;
% both defaults percentage
percentageOfTwoDefault = 1-percentageOfNoDefault-percentageOfOnlyIssuer1Default-...
    percentageOfOnlyIssuer2Default;
P = [percentageOfNoDefault percentageOfOnlyIssuer1Default
    percentageOfOnlyIssuer2Default percentageOfTwoDefault]*100;

% Contingent and fee leg initialization
contingentLeg = zeros(N,1);
% In fee leg sum the payments when no default occurs up to FTD maturity
feeLeg = (DFCDS'*deltas) * ones(N,1);

% Cases in which defaults occur:
for ii=1:length(result_matrix(:,1))
    % Index of first reset date after default occurs:
    ContingentResetDate = find(deltasStart > result_matrix(ii,1),1);
    % Discount in default event:
    discount_contingentleg = InterpDFviaRates(dates, discounts, dates(1)+365*result_matrix(ii,1));
    % Vector of discounts of reset dates before default:
    discounts_feeleg = DFCDS(1:result_matrix(ii,1)-1);
    % Vector of deltas of reset dates before default:
    deltas_feeleg = deltas(1:result_matrix(ii,1)-1);
    % Contingent and fee legs update:
    if ContingentResetDate == 1
        feeLeg(ii) = 0;
    else
        feeLeg(ii) = (discounts_feeleg'*deltas_feeleg);
    end
    contingentLeg(ii) =(1-R(result_matrix(ii,2)))*discount_contingentleg;
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

end % function firstToDefault2NamesNoAccrual