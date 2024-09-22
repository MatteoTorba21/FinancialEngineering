function [M,errorCRR]=PlotErrorCRR(F0,K,B,T,sigma)
% the function calculates the errors between Black and CRR price for a set
% of default number of steps
%
%INPUT
% F0:        forward price
% K:         strike
% B:         discount factor
% T:         time-to-maturity
% sigma:     volatility
%
%OUTPUT
% M:         row vector of intervals in CRR tree for which the errorCRR is calculated
% errorCRR:  row vector of errors computed with CRR tree

m = 1:10;
M = 2.^m;
errorCRR = zeros(1,length(M));
flag = 1; % Call case
% As error for the CRR we consider the difference in absolute value between
% the exact value (computed with the closed formula) and the value computed
% with CRR
closedFormulaPrice = EuropeanOptionClosed(F0,K,B,T,sigma,flag);
for ii = 1 : length(M)
   errorCRR(ii) = abs(closedFormulaPrice - EuropeanOptionCRR(F0,K,B,T,sigma,M(ii),flag));
end

% we want to verify that the errorCRR for a call rescaled with M
% has the same trend of 1/M 
figure()
grid on;
loglog(M,errorCRR,'r','LineWidth',2);
hold on
loglog(M,1./M,'b','LineWidth',2);
xlabel('Number of intervals M');
ylabel('Error');
legend('errorCRR','1/M');
hold off

end % function PlotErrorCRR 