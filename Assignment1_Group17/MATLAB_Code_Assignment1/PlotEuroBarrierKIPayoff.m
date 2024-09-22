function [] = PlotEuroBarrierKIPayoff(K,KI,flag)
% the function plots the payoff for a European Barrier option
%INPUT
% K:     strike
% KI:    barrier
% flag:  1 call, -1 put

N = 1000;
x = linspace(0,max(KI,K)*1.5,N);

figure();
switch (flag)
    case 1
        payoff = max((x-K).*(x>KI),0);
    case -1
        payoff = max((K-x).*(x>KI),0);
end
plot(x,payoff,"LineWidth", 5);
xlabel('F(T,T)');
ylabel('Payoff');
hold on
switch (flag)
    case 1
        payoffVanilla = max((x-KI),0);
        payoffDigital = (x>KI)*(KI-K);
    case -1
        payoffVanilla = max((KI-x),0);
        payoffDigital = (x<KI)*(KI-K);
end
plot(x,payoffVanilla,"LineWidth",2);
plot(x,payoffDigital,"LineWidth",2);
switch (flag)
    case 1
        legend('European Barrier Payoff','European Call Payoff (K = barrier)', '#(Barrier - K) Digital Payoff')
    case -1
        legend('European Barrier Payoff','European Put Payoff (K = barrier)', '#(K - Barrier) Digital Payoff')
end

end % function PlotEuroBarrierKIPayoff