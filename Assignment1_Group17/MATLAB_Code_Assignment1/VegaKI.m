function vega=VegaKI(F0,K,KI,B,T,sigma,N,flagNum)
% Vega calculation with CRR, Monte-Carlo and ClosedFormula. For CRR and
% Monte-Carlo the central difference approximetion method is used.
%
% INPUT:
% F0:       forward price vector
% K:        strike
% KI        barrier
% B:        discount factor
% T:        time-to-maturity
% sigma:    volatility
% N:        either number of time steps or Monte-Carlo simulations
% flagNum:  1 CRR, 2 Monte Carlo, 3 exact
%
%OUTPUT:
% vega:     vector of values of vega for different forwards expressed in
%           the currency

epsilon = sigma*0.1; % increment for derivative approximation
vega = zeros(size(F0));

for ii = 1:length(F0)
switch (flagNum)
    case 1  % CRR
        vega(ii) = (EuropeanOptionKICRR(F0(ii),K, KI,B,T,sigma + epsilon,N) - ...
         EuropeanOptionKICRR(F0(ii),K, KI,B,T,sigma - epsilon,N))/(2*epsilon) * 0.01;
    case 2  % MonteCarlo
        vega(ii) = (EuropeanOptionKIMC(F0(ii),K, KI,B,T,sigma + epsilon,N) - ...
         EuropeanOptionKIMC(F0(ii),K, KI,B,T,sigma - epsilon,N))/(2*epsilon) * 0.01;
    case 3  % closedFormula
        d2 = log(F0(ii)/KI)/(sigma*sqrt(T)) - sigma * sqrt(T)/2;
        d1 = d2 + sigma * sqrt(T);
        vegaDigital = (KI - K ) * normpdf(d2) * (- log(F0(ii)/KI)/(sqrt(T) * sigma^2 ) - sqrt(T) / 2);
        vegaCall = B * F0(ii) * sqrt(T) * normpdf(d1);
        vega(ii) = (vegaCall + vegaDigital)*0.01;
        %The following is numerical estimation for vega using the black formula
        % vega = (EuropeanCallKIClosed(F0,K, KI,B,T,sigma + epsilon) - EuropeanCallKIClosed(F0,K, KI,B,T,sigma - epsilon))/(2*epsilon); 
    otherwise
end
end

return
end % function VegaKI