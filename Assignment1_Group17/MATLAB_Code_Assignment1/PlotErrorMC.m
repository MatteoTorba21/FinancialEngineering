function [M,stdEstim]=PlotErrorMC(F0,K,B,T,sigma)
% compute and plot the Monte-Carlo error
%INPUT
% F0:        forward price
% K:         strike
% B:         discount factor
% T:         time-to-maturity
% sigma:     volatility
%
%OUTPUT
% M:         row vector of number of simulation for which errorMC is calcutaed 
% stdEstim:  row vector of errors computed with Monte-Carlo method

% we are considering only the Call case 
m = 1:20;
M = 2.^m;
% As error for Monte-Carlo method we consider the unbiased
% standard deviation of the Monte-Carlo price
stdEstim = zeros(1,length(M));
for ii = 1 : length(M)
    g = randn(1,M(ii));
    % the discounted payoff of a Call is B*max(F0*exp(-(sigma^2)*T/2+sigma*sqrt(T))-K)
    callPrices = B * max(F0*exp(-(sigma^2)*T/2+sigma*sqrt(T)*g)-K,0);
    stdEstim(ii) = std(callPrices)/sqrt(M(ii));
end

% verify that the errorMC for a call rescaled with M has the same trend of 1/sqrt(M)
figure()
grid on
loglog(M,stdEstim,'r','LineWidth',2);
hold on 
loglog(M,1./sqrt(M),'b','LineWidth',2);
xlabel('M simulations');
ylabel('Error');
legend('errorMC','1/sqrt(M)');
hold off

end %function PlotErrorMC
