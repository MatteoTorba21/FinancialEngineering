function Q = Transition_Matrix(IG_h_curve, HY_h_curve)
% Computation of Transition Matrix, trhough the use of the homogeneity
% property of the Markov chain
%
% INPUT:
% IG_h_curve:       Hazard rates vector for the Investment Grade Bonds
% HY_h_curve:       Hazard rates vector for the High Yield Bonds
%
% OUTPUT:
% Q:                Transition Matrix

% Default probability matrixes: [maturities, default probabilities]
IG_p_default = [1, 1 - exp(-IG_h_curve(1,2)); 2, 1 - exp(-sum(IG_h_curve(:,2)))];
HY_p_default = [1, 1 - exp(-HY_h_curve(1,2)); 2, 1 - exp(-sum(HY_h_curve(:,2)))];

syms x y %variables to find

% Transition Matrix with variables
Q1 = [1-x-IG_p_default(1,2), x, IG_p_default(1,2);
      y, 1-y-HY_p_default(1,2), HY_p_default(1,2);
      0, 0, 1];

% From the Chapman-Kolmogorov equation, formulate a system
% equalizing the default probabilities at two years and the ones computed
% through the product of the Transition Matrix (from t = 0 to t = 1y) by
% itself
sol = solve([Q1(1,:)*(Q1(:,3)) - IG_p_default(2,2), Q1(2,:)*(Q1(:,3)) - HY_p_default(2,2)], [x, y]);

xx = double(sol.x);
yy = double(sol.y);

Q = [1-xx-IG_p_default(1,2), xx, IG_p_default(1,2);
     yy, 1-yy-HY_p_default(1,2), HY_p_default(1,2);
     0, 0, 1];

end % function Transition_Matrix

