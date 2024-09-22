function Put_Prices = JamshidianPrices(a,sigma,Maturity,t0,T_start,fwdD_matrix,K,deltas,discountsyearly,flag_plot)
% Function to compute the prices of the Put options on coupon bond
% (corresponding to the European Swaption prices) via the Jamshidian
% formula for different expiries.
%
% INPUTS:
% a:                Hull-White parameter a
% sigma:            Hull-White parameter σ
% Maturity:         Maturity of the underlying
% t0:               Time in which compute the price
% T_start:          Start time to enter in the swaption
% fwdD_matrix:      Forward discounts matrix
% K:                Strike
% deltas:           Vector of yearly year fractions
% discountsyearly:  Vector of yearly discounts
% flag_plot:        Flag variable to determine whether the plot for the
%                   "graphical resolution" of the equation to find X* must
%                   be showed or not
%
% OUPUT:
% Put_Prices:       Vector of the prices of the Put options with expiries
%                   from T_start to Maturity-1

% Function handles to compute the integral value in the exponent:
sigma_integrand = @(u,t) sigma*(1-exp(-a*(t-u)))/a;
sigma_difference = @(u,t,tau) sigma_integrand(u,t+tau).^2-sigma_integrand(u,t).^2;
integral_value = @(t_alpha) integral(@(u) sigma_difference(u,t_alpha*ones(Maturity-t_alpha,1),(t_alpha+1:Maturity)'-t_alpha),t0,t_alpha,'ArrayValued', true,"AbsTol",1e-12,"RelTol",1e-12);

% Strike of the put coupon bond option equal to the swaption:
K_CB = 1;
% Function handle to compute the forward discounts:
fwd_discount = @(x,t_alpha) fwdD_matrix(t_alpha:end,t_alpha-T_start+1).* ...
    exp(-x/sigma .* (sigma_integrand(t_alpha*ones(Maturity-t_alpha,1),(t_alpha+1:Maturity)')-sigma_integrand(t_alpha,t_alpha))...
   -1/2*integral_value(t_alpha));
% Coupons vector:
coupons = [K*deltas(1:end-1); 1+K*deltas(end)];

% Initialization of the vectors:
x_star = zeros(Maturity-T_start,1);
Call_Prices = zeros(Maturity-T_start,1);
Put_Prices = zeros(Maturity-T_start,1);

for ii=T_start:Maturity-1
    % Numerical procedure to compute X*:
    Price = @(x) coupons(ii-T_start+1:end)'*fwd_discount(x,ii);
    eqn = @(x) Price(x)-K_CB;
    options = optimoptions('fsolve', 'TolFun', 1e-10, 'TolX', 1e-12,'Display', 'off');
    x_star(ii-T_start+1) = fsolve(eqn,.0,options);
    % Strikes computed from the X* value found:
    K = fwd_discount(x_star(ii-T_start+1),ii);
    % Maturities vector:
    Maturities = (ii+1:Maturity)'; 
    % Zero coupon call prices via the formula in the Gaussian HJM framework:
    ZC_Call_Prices = ZC_Call_Gaussian_HJM(sigma_integrand,discountsyearly(ii-T_start+1),fwdD_matrix(ii:end,ii-T_start+1),ii,Maturities,K,0);
    % Call option price:
    Call_Prices(ii-T_start+1) = coupons(ii-T_start+1:end)'*ZC_Call_Prices;
    % Put option price computed via Put-Call parity:
    Put_Prices(ii-T_start+1) = Call_Prices(ii-T_start+1)-discountsyearly(ii-T_start+2:end)'*coupons(ii-T_start+1:end)+K_CB*discountsyearly(ii-T_start+1);
end

if flag_plot==1 % Plot the "Graphical resolution" 
    % Vectors for the x axis:
    xx = [-0.5:0.001:0.5];
    x_stop = [-0.5:0.001:x_star];
    % Calculate y axis values:
    yy = Price(xx);
    
    figure();
    plot(xx, yy, 'LineWidth', 2); % Plot the Price function
    hold on;
    grid on;
    plot(x_stop, K_CB * ones(size(x_stop)), 'r--', 'LineWidth', 3); % Plot the horizontal line
    plot([x_star, x_star], [0, K_CB], 'r--', 'LineWidth', 2); % Plot a vertical line at x_star
    hold off;
    xlabel('X*', FontSize=24)
    ylabel('Prices', FontSize=24)
    title('Graphical resolution to find X*',FontSize=24)
end

end