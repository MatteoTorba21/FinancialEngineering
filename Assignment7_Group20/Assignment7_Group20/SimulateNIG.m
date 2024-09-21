function Asset = SimulateNIG(deltas,sigma,k,nu,Nsim,F0,discounts,div,lnL)
% This function simulates the asset under the assumption of
% INPUTS
% deltas:           vector of the yearfractions of dates, where the asset is simulated
% sigma:            NIG average volatility parameter
% k:                NIG vol-of-vol parameter
% nu:               NIG skewness parameter
% F0:               Forward price with delivery at the end of deltas' date
% discouts:         Discounts at simulation dates
% div:              Dividend yield
% lnL:              log of Laplace exponent function
%
% OUTPUT
% Asset:            Asset matrix simulation (rows are simulations, columns are the value of the simulated asset in the date)
n = length(deltas);
rng(126)
% Uniform sampling for Inverse gaussian sampling
U = rand(Nsim,n);
% Standard normal sampling
Z = randn(Nsim,n);

% log-forward prices simulation with AV
F = zeros(Nsim,n+1);
F_minus = zeros(Nsim,n+1);
F(:,1) = F0;
F_minus(:,1) = F0;

% Asset matrices initializzation
S = zeros(Nsim,n);
S_minus = zeros(Nsim,n);

% Forward discounts and rates computation
fwDiscounts = discounts(2:end)./discounts(1:end-1);
fwRates = -log(fwDiscounts)./deltas;

for ii =1:n
    % Inverse gaussian sampling via icdf
    G = icdf('InverseGaussian',U(:,ii),1,deltas(ii)/k);
    % log-forward prices computation
    fMC = sigma * sqrt(G)*sqrt(deltas(ii)) .* Z(:,ii) - (.5 + nu) * sigma^2 * G *deltas(ii) - lnL(nu);
    fMC_minus = - sigma * sqrt(G)*sqrt(deltas(ii)) .* Z(:,ii) - (.5 + nu) * sigma^2 * G *deltas(ii) - lnL(nu);
    % Simulated forward prices
    F(:,ii+1)= F(:,ii) .* exp(fMC);
    F_minus(:,ii+1)= F_minus(:,ii) .* exp(fMC_minus);

    % Obtain spot from forward price
    if isempty(fwRates(ii+1:end))
        S(:,ii) = F(:,ii+1);
        S_minus(:,ii) = F_minus(:,ii+1);
    else
        S(:,ii) = F(:,ii+1)*exp(-sum(fwRates(ii+1:end)-div).*deltas(ii+1:end));
        S_minus(:,ii) = F_minus(:,ii+1)*exp(-sum(fwRates(ii+1:end)-div).*deltas(ii+1:end));
    end
end

% return AV simulated spot prices as a unique matrix
Asset=[S;S_minus];
end
