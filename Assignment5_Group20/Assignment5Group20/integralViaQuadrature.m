function [integral_value] = integralViaQuadrature(f,x_values)
% Function which computes the integral in the Lewis formula for the price
% of a call option via the Matlab function integral
% 
% INPUTS:
% f:            function to be integrated         
% x_values:     values of moneyness in which the function has to be evaluated          
% OUTPUT:
% I:            value of the integral

% Define the function I(x) which computes the integral with respect to u
I = @(x) integral(@(u) f(u, x), -Inf, Inf, "AbsTol",1e-12,"RelTol",1e-12); 
% Compute the function I(x) for a specific value of x
integral_value = (arrayfun(@(x) I(x), x_values));

end