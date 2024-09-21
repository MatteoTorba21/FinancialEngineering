%% ES 3: PRICING
%% FFT
clear;
load("cSelect20230131_B.mat")
daycount365 = 3;    % yearfrac Act/365

% Exercise data:
settlementDate = datenum('19-Feb-2008');  % settlement date
k = 1;                                    % vol-of-vol
sigma = 20/100;                           % average volatility
nu = 3;                                   % volatility skew

endDate = addtodate(settlementDate, 1, 'year'); % final date
% Check whether it is a business date and eventually take the next business
% date according to the 'follow' convenction:
Holidays = [datenum('21-Mar-2008') datenum('24-Mar-2008') datenum('15-Ago-2008') datenum('24-Dec-2008') datenum('25-Dec-2008') datenum('26-Dec-2008') datenum('31-Dec-2008') datenum('01-Jan-2009') datenum('10-Apr-2009')];
if ~isbusday(endDate,Holidays)
    endDate = busdate(FLDates(nonBusinessDays), 'follow', Holidays);
end

TTM = yearfrac(settlementDate,endDate,daycount365);  % time to maturity
S0 = cSelect.reference;            % spot price
x_values = [-0.25:0.01:0.25];      % required moneynesses
dividend = cSelect.dividend;       % dividend
discount_TTM = 0.961345775982266;  % discount at maturity (previously computed)

F0 = S0*exp(-dividend*TTM)/discount_TTM;  % fwd price
K = F0*exp(-x_values);                    % Required strikes from moneyness

% NIG alpha parameter
alpha = .5;
% Logarithm of Laplace transform function
lnL = @(w) TTM./k*(1-alpha)/alpha*(1-(1+(w.*k*sigma^2)/(1-alpha)).^alpha);
% Characteristic function
Phi = @(x) exp(-1i*x.*lnL(nu) + lnL((x.^2 + 1i*(1+2*nu).*x)./2));
% Function to be Fourier transformed
f = @(x) 1/(2*pi).*Phi(-x-1i/2).*1./(x.^2 + 1/4);

% Fast Fourier Transform
tic
% Number of intervals (to use FFT it must be a power of 2)
M = 15;
N = 2^M;
% x inizialization (x is moneyness vector)
dx = 0.0025;
x1 = -dx*(N-1)/2;
x = [x1:dx:-x1];
% u inizialization (u is fourier space variable)
du = 2*pi/(dx*N);
u1 = -du*(N-1)/2;
u = [u1:du:-u1];
% Discrete Fourier transform using respective function
dft = integralViaFFT(f,M,dx);
% Call prices from integral computation
C = S0*exp(-dividend*TTM)*(1-exp(-x/2).*dft);
% Interpolation of the results for the required moneyness values:
C_fft = interp1(x,C,x_values);
computational_time_FFT1 = toc;
% Plot the results
figure(1)
plot(x_values,C_fft);
grid on
hold on


%% Quadrature
% Integrand function
f_u = @(u,x) 1/(2*pi).*Phi(-u-1i/2).*1./(u.^2 + 1/4).*exp(-1i*u.*x);
tic     % start timer
% Integral computation via respective function
integral_value = integralViaQuadrature(f_u,x_values);

% Call prices from integral computation
CQuadrature = real(S0*exp(-dividend*TTM)*(1-exp(-x_values/2).*integral_value));
% Computational time upload
computational_time_quadrature = toc;
% Plot of the results
plot(x_values,CQuadrature,'--');

%% Montecarlo
tic
Nsim = 10000;    % Number of simulations
rng(126)
% Uniform sampling for Inverse gaussian sampling
U = rand(Nsim,1);
% Standard normal sampling
Z = randn(Nsim,1);
% Inverse gaussian sampling via icdf
G = icdf('InverseGaussian',U,1,TTM/k);
% log-forward prices simulation with AV
fMC = sigma * sqrt(G)*sqrt(TTM) .* Z - (.5 + nu) * sigma^2 * G *TTM - lnL(nu);
fMC_minus = - sigma * sqrt(G)*sqrt(TTM) .* Z - (.5 + nu) * sigma^2 * G *TTM - lnL(nu);
% forward prices simulation values
F1 = F0 * exp(fMC);
F1_minus = F0 * exp(fMC_minus); 
% Payoff computations for the call option
F_matrix = repmat(F1, 1, length(K));
F_minus_matrix = repmat(F1_minus, 1, length(K));
K_matrix = repmat(K, Nsim, 1);
Payoff = max(F_matrix-K_matrix,0);
Payoff_minus = max(F_minus_matrix-K_matrix,0);
% Discounted payoffs
discountedPayoff = discount_TTM * Payoff;
discountedPayoff_minus = discount_TTM * Payoff_minus;
% MC prices and respective confidence intervals
[Prices,~, MCCI] = normfit((discountedPayoff+discountedPayoff_minus)/2);
computational_time_MC = toc;
% Plot of the results
plot(x_values,Prices,'--');
plot(x_values,MCCI(1,:),'--');
plot(x_values,MCCI(2,:),'--');
legend('FFT','Quadrature','MC','MC bid','MC ask', Location='northwest',FontSize=20)
xlabel('Moneyness', FontSize=18)
ylabel('Call price',FontSize=18)
title('Prices of the call options for different moneyness values for α=1/2', FontSize=22)
hold off

% Compute the errors
FFT_error = max(abs(C_fft-CQuadrature));
MC_error = max(abs(Prices-CQuadrature));
% Display the errors
disp('--- 3. EXERCISE: PRICING --- ')
disp(' ')
fprintf(['Errors using the quadrature as a benchmark and computational times:\nFFT error = %.10f ' ...
    '  \nMC error = %.10f \nComputational time FFT = %.10f \nComputational time MC = %.10f \n\n'],FFT_error,MC_error,computational_time_FFT1,computational_time_MC)

%% FFT Computational time and error for different values of M
% Using Integral as benchmark, we compared fft method with different values
% of M
maxM = 24;
computational_time = zeros(maxM,1);
error = zeros(maxM,1);
for M = [1:maxM]
    % FFT parameters
    N = 2^M;
    dx = 0.0025;
    x1 = -dx*(N-1)/2;
    x = [x1:dx:-x1];
    du = 2*pi/(dx*N);
    u1 = -du*(N-1)/2;
    u = [u1:du:-u1];
    tic
    dft = integralViaFFT(f,M,dx);
    % Error upload (using quadrature method as a benchmark)
    % error(M) = max(abs(integral_value-interp1(x,dft,x_values)));
    C = S0*exp(-dividend*TTM)*(1-exp(-x/2).*dft);
    % Interpolation of the results for the required moneyness values:
    C_fft_iter = interp1(x,C,x_values);
    % Computational time upload
    computational_time(M) = toc;
    % Error among prices upload (using quadrature method as a benchmark)
    error(M) = max(abs(C_fft_iter-CQuadrature));
end

% Execution time plot
figure()
loglog(2.^[1:maxM], computational_time,"*b",'LineWidth',4)
grid on
hold on
loglog(linspace(2^1,2^maxM,100), computational_time_quadrature*ones(1,100),"--",'LineWidth',2)
loglog(2.^[1:maxM], 2.^[1:maxM].*log(2.^[1:maxM]),'LineWidth',4)
loglog(2.^[1:maxM], computational_time,"--b",'LineWidth',1)
legend('Computational Time FFT', 'Computational Time Quadrature', 'N log(N)',Location='northwest', Fontsize=24)
title('FFT Computational time for diffferent values of N',FontSize=28)
xlabel('Number of intervals',FontSize=28)
ylabel('Time',FontSize=28)
hold off

% Error plot
figure()
loglog(2.^[1:maxM], error,"*b",'LineWidth',4)
grid on
hold on
loglog(2.^[1:maxM], error,"--b",'LineWidth',1)
title('FFT error for different values of N using Quadrature method as benchmark',FontSize=28)
xlabel('Number of intervals',FontSize=28)
ylabel('Error',FontSize=28)
hold off

%% FFT Computational time and error for different values of dx
% Using Integral as benchmark, we compared fft method with different values
% of dx

coefficients = [0.125; 0.25; 0.5; 1; 2; 5; 10]; % 0.125; 0.25;
dx = 0.0025.*coefficients;
computational_time_dx = zeros(length(coefficients),1);
error_dx = zeros(length(coefficients),1);
for ii=1:length(coefficients)
    % FFT parameters
    dx_iter = dx(ii);
    M = 15;
    N = 2^M;
    x1 = -dx_iter*(N-1)/2;
    x = [x1:dx_iter:-x1];
    du = 2*pi/(dx_iter*N);
    u1 = -du*(N-1)/2;
    u = [u1:du:-u1];
    tic
    dft = integralViaFFT(f,M,dx_iter);
    % Error upload (using quadrature method as a benchmark)
    % error_dx(ii) = max(abs(integral_value-interp1(x,dft,x_values)));
    C = S0*exp(-dividend*TTM)*(1-exp(-x/2).*dft);
    % Interpolation of the results for the required moneyness values:
    C_fft_iter = interp1(x,C,x_values);
    % Computational time upload
    computational_time_dx(ii) = toc;
    % Error among prices upload (using quadrature method as a benchmark)
    error_dx(ii) = max(abs(C_fft_iter-CQuadrature));
end

% Execution time plot
figure()
loglog(dx, computational_time_dx,"*b",'LineWidth',4)
grid on
hold on
loglog(dx, computational_time_quadrature*ones(1,length(dx)),"--",'LineWidth',2)
loglog(dx, computational_time_dx,"--b",'LineWidth',1)
legend('Computational Time FFT', 'Computational Time Quadrature',Location='northwest',FontSize=24)
title('FFT Computational time for diffferent values of dx',FontSize=28)
ylim([10^(-3) 1])
xlabel('Step of the grid (dx)',FontSize=28)
ylabel('Time',FontSize=28)
hold off

% Error plot
figure()
loglog(dx, error_dx,"*b",'LineWidth',4)
grid on
hold on
loglog(dx, error_dx,"--b",'LineWidth',1)
title('FFT error for different values of dx using Quadrature method as benchmark',FontSize=28)
xlabel('Step of the grid (dx)',FontSize=28)
ylabel('Error',FontSize=28)
hold off

%% facultative : FFT
% Model parameters and fft computations
alpha = 2/3;
% Logarithm of Laplace transform function
lnL = @(w) TTM./k*(1-alpha)/alpha*(1-(1+(w.*k*sigma^2)/(1-alpha)).^alpha);
% Characteristic function
Phi = @(x) exp(-1i*x.*lnL(nu) + lnL((x.^2 + 1i*(1+2*nu).*x)./2));
% Function to be Fourier transformed
f = @(x) 1/(2*pi).*Phi(-x-1i/2).*1./(x.^2 + 1/4);

% Fast Fourier Transform
% Number of intervals (to use FFT it must be a power of 2)
M = 15;
N = 2^M;
% x inizialization (x is moneyness vector)
dx = 0.0025;
x1 = -dx*(N-1)/2;
x = [x1:dx:-x1];

% Discrete Fourier transform using respective function
dft = integralViaFFT(f,M,dx);
% Call prices from integral computation
C_fac = S0*exp(-dividend*TTM)*(1-exp(-x/2).*dft);
% Interpolation of the results for the required moneyness values:
C_fft_fac = interp1(x,C_fac,x_values);

% Plot the results
figure()
plot(x_values,C_fft_fac);
grid on
hold on

%% facultative : quadrature method, different alpha
% Integrand function
f_u = @(u,x) 1/(2*pi).*Phi(-u-1i/2).*1./(u.^2 + 1/4).*exp(-1i*u.*x);
% Integral computation via respective function
I = integralViaQuadrature(f_u,x_values);
% Call prices from integral computation
CQuadrature_fac = real(S0*exp(-dividend*TTM)*(1-exp(-x_values/2).*I));
% Plot the results
plot(x_values,CQuadrature_fac,'--');

legend('FFT','Quadrature',Location='northwest',FontSize=20)
xlabel('Moneyness', FontSize=18)
ylabel('Call price',FontSize=18)
title('Prices of the call options for different moneyness values for α=2/3', FontSize=22)
hold off

% Compute the error
FFT_error_fac = max(abs(C_fft_fac-CQuadrature_fac));
% Display the error
disp('Facultative point with α=2/3: ')
fprintf('Errors using the quadrature as a benchmark:\nFFT error = %.10f  \n\n', FFT_error_fac)

% Plot the FFT errors for different values of alpha
figure()
plot(x_values,C_fft-CQuadrature)
hold on
plot(x_values,C_fft_fac-CQuadrature_fac)
legend('Error of FFT method with α=1/2','Error of FFT method with α=2/3')
hold off

% Plot the prices with α=1/2 VS the prices with α=2/3
figure()
plot(x_values,C_fft);
grid on
hold on
plot(x_values,C_fft_fac,'--');
legend('Call prices via FFT method with α=1/2','Call prices via FFT method with α=2/3',Location='northwest',FontSize=20)
xlabel('Moneyness', FontSize=18)
ylabel('Call price',FontSize=18)
title('Prices of the call options for different moneyness values for α=1/2 and α=2/3', FontSize=22)
hold off
hold off

% Plot the call prices difference for the two values of α
figure()
plot(x_values,C_fft-C_fft_fac, LineWidth=3);
hold on
plot(x_values,zeros(length(x_values),1),'--k', LineWidth=1);
grid on
xlabel('Moneyness', FontSize=18)
ylabel('Call price difference',FontSize=18)
title('Difference between the prices for different moneyness values for α=1/2 and α=2/3', FontSize=22)
hold off

% Max percentage error:
alpha_error_percentage = 100*max(abs((C_fft_fac-C_fft)./C_fft));
fprintf('Max percentage error for the FFT computed prices for different alphas: error = %.10f  \n\n', alpha_error_percentage)


%% FFT Computational time and error for different values of M for α=2/3
% Using Integral as benchmark, we compared fft method with different values
% of M
maxM = 24;
computational_time_23 = zeros(maxM,1);
error_23 = zeros(maxM,1);
for M = [1:maxM]
    % FFT parameters
    N = 2^M;
    dx = 0.0025;
    x1 = -dx*(N-1)/2;
    x = [x1:dx:-x1];
    du = 2*pi/(dx*N);
    u1 = -du*(N-1)/2;
    u = [u1:du:-u1];
    tic
    dft = integralViaFFT(f,M,dx);
    % Error upload (using quadrature method as a benchmark)
    C = S0*exp(-dividend*TTM)*(1-exp(-x/2).*dft);
    % Interpolation of the results for the required moneyness values:
    C_fft_iter_23 = interp1(x,C,x_values);
    % Computational time upload
    computational_time_23(M) = toc;
    % Error among prices upload (using quadrature method as a benchmark)
    error_23(M) = max(abs(C_fft_iter_23-CQuadrature));
end

% Error plot
figure()
loglog(2.^[1:maxM], error,"*b",'LineWidth',4)
grid on
hold on
loglog(2.^[1:maxM], error_23,"*r",'LineWidth',4)
loglog(2.^[1:maxM], error_23,"--r",'LineWidth',1)
loglog(2.^[1:maxM], error,"--b",'LineWidth',1)
title('FFT error wrt Quadrature method varying N for α=1/2 and α=2/3',FontSize=28)
xlabel('Number of intervals',FontSize=28)
ylabel('Error',FontSize=28)
legend('Error with α=1/2','Error with α=2/3',FontSize=20)
hold off
