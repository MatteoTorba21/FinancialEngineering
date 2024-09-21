function [BermudanPrice,EuropeanSwaptionPrices] = buildTrinomialTree(T_start,K,Maturity,N,a,sigma,dates,discounts,fwdD_matrix,deltas,flag_disp)
% Function to price a Bermudan Swaption and the European Swaption for all
% the years in which the exercise of the Bermudan is possible via a
% trinomial tree, in the Hull-White framework.
%
% INPUTS:
% T_start:         Start time to enter in the swaption 
% K:               Swaption's strike 
% Maturity:        Maturity of the underlying swap 
% N:               Number of time steps in a year 
% a:               Hull-White parameter a 
% sigma:           Hull-White parameter σ 
% dates:           Dates from the bootstrap 
% discounts:       Discounts from the bootstrap 
% fwdD_matrix:     Forward discounts matrix 
% deltas:          Vector of yearly year fractions 
% flag_disp:       Flag variable to activate/deactivate the display of
%                  early exercises each year:
%                  -If flag==1 display the early exercises
%                  -If flag==0 do not display the early exercises
%
% OUTPUT:
% BermudanPrice:            Price of the Bermudan Swaption
% EuropeanSwaptionPrices:   Vector of European Swaptions' prices

% Convenctions for yearfrac function:
ACT365 = 3;
% Time step
dt = 1/N;          

% delta x parameters:
mu_hat = 1-exp(-a*dt);
sigma_hat = sigma*sqrt((1-exp(-2*a*dt))/(2*a));

% Maximum value for l:
l_max = min(ceil((1-sqrt(2/3))/mu_hat),floor(sqrt(2/3)/mu_hat));

dx = sqrt(3)*sigma_hat; % dx
l = (l_max:-1:-l_max)'; % l vector
x = l*dx;               % x vector

% Vector of probabilities
pu = 1/2*(1/3-mu_hat*l+mu_hat^2*l.^2);
pm = 2/3-mu_hat^2*l.^2;
pd = 1/2*(1/3+mu_hat*l+mu_hat^2*l.^2);
% Modify the values in the case in which lmax steps up have been made
pu(1) = 1/2*(7/3-3*l(1)*mu_hat+(l(1)*mu_hat)^2);
pm(1) = (-1/3+2*l(1)*mu_hat-(l(1)*mu_hat)^2);
pd(1) = 1/2*(1/3-l(1)*mu_hat+(l(1)*mu_hat)^2);
% Modify the values in the case in which lmax steps down have been made
pu(end) = 1/2*(1/3+l(end)*mu_hat+(l(end)*mu_hat)^2);
pm(end) = (-1/3-2*l(end)*mu_hat-(l(end)*mu_hat)^2);
pd(end) = 1/2*(7/3+3*l(end)*mu_hat+(l(end)*mu_hat)^2);

% Define the starting time t0
t0 = 0;
% Times vector
t = (0:dt:Maturity-1)';

% Function handles to compute the integral value in the exponent:
sigma_integrand = @(u,t) sigma*(1-exp(-a*(t-u)))/a;
sigma_difference = @(u,t,tau) sigma_integrand(u,t+tau).^2-sigma_integrand(u,t).^2;

% Compute the discounts in 0 via interpolation on the bootstrap results:
discounts_grid = InterpDFviaRates(yearfrac(dates(1),dates,ACT365),discounts,t);
% Compute the forward discounts in t=0 from the bootsrap results
fwd_discounts_grid = discounts_grid(2:end)./discounts_grid(1:end-1);
fwd_discounts_grid = [1; fwd_discounts_grid];
% Ausiliary variable:
jj = length(fwd_discounts_grid);
% Value in the exponent of the approximated ratio of discounts:
sigma_star_hat = sigma/a * sqrt( dt - 2*(1-exp(-a*dt))/a + (1-exp(-2*a*dt))/(2*a) );

% Continuation value vector initialization
continuation_value = zeros(length(x),1);
% Initialization of the payoff matrix
Payoff_matrix = zeros(length(l),Maturity-T_start);
% Initialization of the European payoff matrix
Payoff_EU = zeros(length(l),Maturity-T_start);

for ii=Maturity-1:-dt:0
    % Integral in the forward discount expression
    integral_value = integral(@(u) sigma_difference(u,ii,dt),t0,ii+dt,"AbsTol",1e-12,"RelTol",1e-12);
    % Coefficient in the forward discount expression
    coeff_x = sigma_integrand(t0,dt) / sigma;
    % Fowrard discounts grid: for the value T_alpha=ii, a matrix of
    % discounts for different x values is built:
    B_grid = exp(-x.*coeff_x-1/2*integral_value).*fwd_discounts_grid(jj);

    % Compute the ratios between the stochastic discounts and the discounts:
    DiscountRatios_up = exp(-1/2*sigma_star_hat^2 - sigma_star_hat/sigma_hat * (exp(-a*dt)*dx + mu_hat*x)); 
    DiscountRatios_mid = exp(-1/2*sigma_star_hat^2 - sigma_star_hat/sigma_hat * mu_hat*x);
    DiscountRatios_down = exp(-1/2*sigma_star_hat^2 - sigma_star_hat/sigma_hat * (-exp(-a*dt)*dx + mu_hat*x));
    % Also in the 'extremities' cases:
    DiscountRatios_upabove = exp(-1/2*sigma_star_hat^2 - sigma_star_hat/sigma_hat * (mu_hat*x(1)));
    DiscountRatios_upbelow = exp(-1/2*sigma_star_hat^2 - sigma_star_hat/sigma_hat * (exp(-a*dt)*2*dx + mu_hat*x(end)));
    DiscountRatios_midabove = exp(-1/2*sigma_star_hat^2 - sigma_star_hat/sigma_hat * (-exp(-a*dt)*dx + mu_hat*x(1)));
    DiscountRatios_midbelow = exp(-1/2*sigma_star_hat^2 - sigma_star_hat/sigma_hat * (exp(-a*dt)*dx + mu_hat*x(end)));
    DiscountRatios_downabove = exp(-1/2*sigma_star_hat^2 - sigma_star_hat/sigma_hat * (-exp(-a*dt)*2*dx + mu_hat*x(1)));
    DiscountRatios_downbelow = exp(-1/2*sigma_star_hat^2 - sigma_star_hat/sigma_hat * (mu_hat*x(end)));

    % Save the old continuation value:   
    old_continuation_value = continuation_value;
    % Update the continuation value in different scenarios:
    continuation_up = zeros(length(continuation_value),1);
    continuation_up(2:end-1) = old_continuation_value(1:end-2).*pu(2:end-1).*DiscountRatios_up(1:end-2).*B_grid(2:end-1);
    continuation_mid = zeros(length(continuation_value),1);
    continuation_mid(2:end-1) = old_continuation_value(2:end-1).*pm(2:end-1).*DiscountRatios_mid(2:end-1).*B_grid(2:end-1);
    continuation_down = zeros(length(continuation_value),1);
    continuation_down(2:end-1) = old_continuation_value(3:end).*pd(2:end-1).*DiscountRatios_down(3:end).*B_grid(2:end-1);
    % In the 'extremities' cases:
    continuation_up(1) = old_continuation_value(1)*pu(1)*DiscountRatios_upabove*B_grid(1);
    continuation_up(end) = old_continuation_value(end-2)*pu(end)*DiscountRatios_upbelow*B_grid(end);
    continuation_mid(1) = old_continuation_value(2)*pm(1)*DiscountRatios_midabove*B_grid(1);
    continuation_mid(end) = old_continuation_value(end-1)*pm(end)*DiscountRatios_midbelow*B_grid(end);
    continuation_down(1) = old_continuation_value(3)*pd(1)*DiscountRatios_downabove*B_grid(1);
    continuation_down(end) = old_continuation_value(end)*pd(end)*DiscountRatios_downbelow*B_grid(end);
    % Update the discounted continuation value:
    continuation_value = continuation_up + continuation_mid + continuation_down;
    
    % Discount the european payoffs: 
    old_european_value = Payoff_EU;
    % Update the european payoff in different scenarios:
    european_up = zeros(size(Payoff_EU));
    european_up(2:end-1,:) = old_european_value(1:end-2,:).*pu(2:end-1).*DiscountRatios_up(1:end-2).*B_grid(2:end-1);
    european_mid = zeros(size(Payoff_EU));
    european_mid(2:end-1,:) = old_european_value(2:end-1,:).*pm(2:end-1).*DiscountRatios_mid(2:end-1).*B_grid(2:end-1);
    european_down = zeros(size(Payoff_EU));
    european_down(2:end-1,:) = old_european_value(3:end,:).*pd(2:end-1).*DiscountRatios_down(3:end).*B_grid(2:end-1);
    % In the 'extremities' cases:
    european_up(1,:) = old_european_value(1,:)*pu(1)*DiscountRatios_upabove*B_grid(1);
    european_up(end,:) = old_european_value(end-2,:)*pu(end)*DiscountRatios_upbelow*B_grid(end);
    european_mid(1,:) = old_european_value(2,:)*pm(1)*DiscountRatios_midabove*B_grid(1);
    european_mid(end,:) = old_european_value(end-1,:)*pm(end)*DiscountRatios_midbelow*B_grid(end);
    european_down(1,:) = old_european_value(3,:)*pd(1)*DiscountRatios_downabove*B_grid(1);
    european_down(end,:) = old_european_value(end,:)*pd(end)*DiscountRatios_downbelow*B_grid(end);
    % Update the discounted european payoffs:
    Payoff_EU = european_up + european_mid + european_down;
    
    % Check for early exercise:
    if mod(ii,1)==0 && ii>=T_start 
        % Initialization of the forward discounts matrix
        B_matrix = zeros(length(x),(Maturity-ii));
        
        % Compute the forward discounts matrix
        for kk=1:Maturity-ii
            integral_value = integral(@(u) sigma_difference(u,ii,kk),t0,ii,"AbsTol",1e-12,"RelTol",1e-12);
            coeff_x = sigma_integrand(t0,kk) / sigma;
            B_matrix(:,kk) = exp(-x.*coeff_x-1/2*integral_value).*fwdD_matrix(ii-1+kk,ii-1);
        end
        % Compute the BPV:
        BPV = B_matrix*deltas(ii-1:end);
        % Compute the Swap rate:
        SwapRate = (1-B_matrix(:,end))./BPV;
        % Instert in the payoff matrix the current payoff:
        Payoff_matrix(:,ii) = max(SwapRate-K,0).*BPV;
        
        % Update the European payoff matrix:
        Payoff_EU(:,ii-T_start+1) = Payoff_matrix(:,ii);
        
        if flag_disp==1
            % Print the number of early exercises in the current year:
            fprintf('At year %d: %d early exercise \n', ii ,sum( continuation_value<Payoff_matrix(:,ii)))
        end
        % Update the continuation value:
        continuation_value = max(continuation_value,Payoff_matrix(:,ii));
    end
    % Update the auxiliary variable:
    jj = jj-1;
end

% Update the results' variables (take the value in the initial position of
% the tree):
BermudanPrice = continuation_value(l_max+1);
EuropeanSwaptionPrices = Payoff_EU(l_max+1,:);
end