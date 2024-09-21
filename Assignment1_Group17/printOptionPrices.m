function [] = printOptionPrices(OptionPrice, Notional)
% print option prices obtained with closed formula, CRR and Monte Carlo
%
%INPUT
% OptionPrice: row vector of the option prices computed with the three methods
% Notional:    falg variable which indicates whether our prices include the
%               mupltiplication for the notional value (Notional==1) or not (Notional==0)

if length(OptionPrice) ~= 3
    error('Wrong input size')
end

if Notional==0
for ii = 1:3
    switch (ii)
        case 1
            disp(['Price with closed formula for one contract: ', num2str(OptionPrice(ii))]);
        case 2
            disp(['Price with CRR for one contract: ', num2str(OptionPrice(ii))]);
        case 3
            disp(['Price with Monte Carlo for one contract: ', num2str(OptionPrice(ii))]);
    end
end
end

if Notional==1
for ii = 1:3
    switch (ii)
        case 1
            disp(['Price with closed formula considering the notional: ', num2str(OptionPrice(ii))]);
        case 2
            disp(['Price with CRR considering the notional: ', num2str(OptionPrice(ii))]);
        case 3
            disp(['Price with Monte Carlo considering the notional: ', num2str(OptionPrice(ii))]);
    end
end
end

end % function OptionPrices
