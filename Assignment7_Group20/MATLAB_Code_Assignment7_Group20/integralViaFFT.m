function I = integralViaFFT(f,M,dx)
% Function which computes the integral in the Lewis formula for the price
% of a call option via the Fast Fourier Transform
% 
% INPUTS:
% f:            function to be Fourier transformed         
% M:            power of 2 to compute the number of intervals          
% dx:           step of the moneyness grid
% OUTPUT:
% I:            value of the integral

% FFT parameters:
N = 2^M;
x1 = -dx*(N-1)/2;
x = [x1:dx:-x1];
du = 2*pi/(dx*N);
u1 = -du*(N-1)/2;
u = [u1:du:-u1];

% Compute the the Fast Fourier Transform:
FFT = fft(f(u).*exp(-1i.*x1.*du*[0:N-1]));
% Compute the integral by multiplying for the prefactor:
I = real(du * exp( -1i .* u1 .* x) .* FFT);

end


