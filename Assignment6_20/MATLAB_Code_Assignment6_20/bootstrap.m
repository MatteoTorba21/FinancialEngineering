function [dates, discounts] = bootstrap(datesSet, ratesSet)
% Bootstrap function given dates and rates of IB deposits, STIR futures and IR Swaps contracts.
% For the first part of the curve depos rates are used, for the second part the STIR futures are used and for 
% the third part the swaps are used. In each part the most liquid contracts are chosen.
%
% INPUT:        
% datesSet:   Struct of dates composed of settlement date, depo dates,
%             future settlement and expiry dates and swap dates
% ratesSet:   Struct of rates composed of depo, future and swap bid/ask
%             
% OUTPUT:
% dates:      Sates at which discount factors are computed
% discounts:  Discount factors computed via bootstrap
%

ACT_360 = 2; % yearfrac ACT_360 (day-count depo and future convenction)
European_30_360 = 6;  % yearfrac 30/360 European (day-count swap convenction)

% Save dates and rates in variables easier to handle than the struct
depoDate = datesSet.depos;
futureDate = datesSet.futures;
swapDate = datesSet.swaps;

depoRate = ratesSet.depos;
futureRate = ratesSet.futures;
swapRate = ratesSet.swaps;

% Initial date and initial discount:
dates = datesSet.settlement;
discounts = 1;

% Short end -> depos
mid = find(depoDate > futureDate(1,1), 1); % to find the last depo to consider for the bootstrap
dates = [dates; depoDate(1:mid)]; % adding the dates needed to dates
depoRate = (depoRate(1:mid,1) + depoRate(1:mid,2))./2; % computing the rate
% Compute the discounts obtained from the depos:
discounts = [discounts; 1./(1+yearfrac(dates(1)*ones(mid,1), dates(2:mid+1), ACT_360).*depoRate)]; 

% Middle area -> futures
futureDate = futureDate(1:7,:); % to consider only the first 7 futures (which are the most liquid)
futureRate = (futureRate(1:7,1) + futureRate(1:7,2))./2; % mean rate of the first 7 futures

for ii = 1 : 7
    % Search the position of the settlement settlement date of the iith future in the dates vector
    % to understand if interpolation or extrapolation is needed:
    Dindex = find(dates >= futureDate(ii,1), 1);
    if isempty(Dindex) % settlement date of the Future is greater than the biggest date in dates
        % Extrapolate the discount at the settlement date:
        Discount_at_Future_settlement = InterpDFviaRates(dates, discounts, futureDate(ii,1));

    elseif futureDate(ii,1) ~= dates(Dindex) % Settlement date of the future is between two dates in the vector dates
        % Intrapolate the discount at the settlement date:
        Discount_at_Future_settlement = InterpDFviaRates(dates, discounts, futureDate(ii,1));

    else % No need for interpolation nor extrapolation case, the discount at settlement date is known:
        Discount_at_Future_settlement = discounts(Dindex);
    end

    % Compute the forward discount 
    FWD_Discount = 1./(1+yearfrac(futureDate(ii,1),futureDate(ii,2),ACT_360).*futureRate(ii));
    % Discount at future expiry:
    Discount_at_Future_expiry = Discount_at_Future_settlement*FWD_Discount;
       
    % Update the discounts vector:
    discounts = [discounts; Discount_at_Future_expiry];   
    % Update dates:
    dates = [dates; futureDate(ii,2)];
end

% Long end -> swaps
firstSwapExpiry = swapDate(1); % expiry of the first swap (not returned in outputs)
firstSwapDF = InterpDFviaRates(dates,discounts,firstSwapExpiry); % discount factor for the first swap (not returned in outputs)

swapDate = [dates(1); swapDate]; % settlement and all swaps maturities (for computations)
swapDF = [firstSwapDF;zeros(length(swapDate)-2,1)]; % preallocate the memory for the discounts for all swaps expiries

for ii = 2:length(swapDF)
    % Find next DF by inverting par rate swap formula
   newDF = (1-mean(swapRate(ii,:))*yearfrac(swapDate(1:ii-1),swapDate(2:ii),European_30_360)'*...
       swapDF(1:ii-1))/(1+yearfrac(swapDate(ii),swapDate(ii+1),European_30_360)*mean(swapRate(ii,:)));
   % Update in the auxiliary vector
   swapDF(ii) = newDF; 
   % Update output vectors of dates and discounts
   dates = [dates; swapDate(ii+1)];
   discounts = [discounts; newDF];
end

end % function bootstrap