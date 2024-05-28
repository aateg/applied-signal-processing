%% 1. Interference Elimination
% Guilherme Fernandes Goncalves Silva
% Num. USP 10297272
clear;
clc;

%% 1. Variables definition
M = 2; % adaptative filter coeficients
N = 1000; % number of points
n = 1:N;
phi_v = (2*pi)*rand(N,1); % U ~ [0, 2pi]
phi_u = phi_v;


s = sqrt(0.01) .* randn(N, 1); % s(n) - white noise (0, 0.01)
x = sin(2*pi*n/10 + pi/6 + phi_v')' ; % x(n) - interference
u = 5*sin(2*pi*n/10 + phi_u')'; % u(n) - corr signal with interference

d = s + x; % d(n) - 
%e = d - y; % e(n) - estimation error

%% 1a. 
% Autocorrelation Matrix (R) calculation
r = xcorr(u, M-1, 'biased');
ru = r(M:end);
R = toeplitz(ru);

% cross-correlation (p) vector calculation
rdu = xcorr(d, u, M-1, 'biased');
p = rdu(M:end);

% optimal coeficients vector (w0)
w0 = R \ p; % inv(R) * p

% MMSE
J_min = var(d) - p'*w0;

% plot resposta em frequencia na freq da interferencia
% compare com x e u

%% 1b. 
% Calcule a faixa de valores do passo de adaptacao mu
% que garante a convergencia do algoritmo Steepest Descent

[V, D] = eig(R);

%% 1c. Algoritmo LMS
mu = 0.03;
N_iter = 500;


%% 1d.

%% 1e.

%% 1f. 