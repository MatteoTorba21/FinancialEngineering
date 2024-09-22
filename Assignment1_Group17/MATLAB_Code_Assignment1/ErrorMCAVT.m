function MCAVTerror = ErrorMCAVT(F0,K,B,T,sigma,flag)
% MonteCarlo error computed with Antithetic Variable Technique
%
%INPUT
% F0:          forward price
% K:           strike
% B:           discount factor
% T:           time-to-maturity
% sigma:       volatility
% flag:        1 call, -1 put
%
%OUTPUT
% MCAVTerror:  row vector of errors computed with Antithetic Variables Technique

m = 1:20;
M = 2.^m;
MCAVTerror = zeros(1,length(M));
if flag==1 % Call case
    for ii = 1 : length(M)
        % sample from normal
        g = randn(1,M(ii)/2);
        % discounted payoff simulations with g and -g
        callPrices = B * max(F0*exp(-(sigma^2)*T/2+sigma*sqrt(T)*g)-K,0);
        callPricesAV = B * max(F0*exp(-(sigma^2)*T/2+sigma*sqrt(T)*(-g))-K,0);
        MCAVTerror(ii) = std((callPrices+callPricesAV)/2)/sqrt(M(ii));
    end
end
if flag==-1 % Put case 
    for ii = 1 : length(M)
        % sample from normal
        g = randn(1,M(ii)/2);
        % discounted payoff simulations with g and -g
        putPrices = B * max(K-F0*exp(-(sigma^2)*T/2+sigma*sqrt*g),0);
        putPricesAV = B * max(K-F0*exp(-(sigma^2)*T/2+sigma*sqrt*(-g)),0);
        MCAVTerror(ii) = std((putPrices+putPricesAV)/2)/sqrt(M(ii));
    end
end

end % function ErrorMCAVT