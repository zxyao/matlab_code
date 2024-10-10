% Parameter setting
N = 10; % Number of antennas for STAR-RIS
Pmax = -10; % Transmit power limitation
sigma2 = 1; % Noise power

% Channel matrix initialization
G = sqrt(1/2) * (randn(N, N) + 1i * randn(N, N));
hI = sqrt(1/2) * (randn(1, N) + 1i * randn(1, N));
hO = sqrt(1/2) * (randn(1, N) + 1i * randn(1, N));
hE1 = sqrt(1/2) * (randn(1, N) + 1i * randn(1, N));
hE2 = sqrt(1/2) * (randn(1, N) + 1i * randn(1, N));

% Initialize variables
beta_t = ones(1, N);
beta_r = zeros(1, N);
theta_t = zeros(1, N);
theta_r = zeros(1, N);

for i = 1:1000
    % Optimization of transmission and reflection coefficients
    beta_t = optimize_beta_t_and_beta_r(beta_t, beta_r, hI, hO, hE1, hE2, G, sigma2, Pmax);
    beta_r = 1 - beta_t;
    % Phase optimization
    theta_t = optimize_theta_t_and_theta_r(theta_t, theta_r, beta_t, hI, hO, hE1, hE2, G, sigma2, Pmax);
    theta_r = theta_t + pi/2;
    % Transmission power constraints
    W = beta_t * beta_t' + beta_r * beta_r';
    if trace(W) > Pmax
        beta_t = beta_t * Pmax / trace(W);
        beta_r = 1 - beta_t;
    end
end

% output
disp(['The optimal beta_t and beta_r are: ', num2str(beta_t)]);
disp(['The optimal theta_t and theta_r are: ', num2str(theta_t)]);

% Optimized function of transmission and reflection coefficients
function beta_t = optimize_beta_t_and_beta_r(beta_t, beta_r, hI, hO, hE1, hE2, G, sigma2, Pmax)
    beta_t = beta_t + 0.01 .* (hI * beta_t - hO * beta_r - hE1 * beta_t - hE2 * beta_r - G * (beta_t * beta_t'));
    beta_t = beta_t / norm(beta_t);
end

% Optimized function of phase
function theta_t = optimize_theta_t_and_theta_r(theta_t, theta_r, beta_t, hI, hO, hE1, hE2, G, sigma2, Pmax)
    theta_t = theta_t + 0.01 .* (hI * exp(1i * theta_t) - hO * exp(1i * theta_r) - hE1 * exp(1i * theta_t) - hE2 * exp(1i * theta_r) - G * beta_t * exp(1i * theta_t));
    theta_t = theta_t / norm(theta_t);
end
