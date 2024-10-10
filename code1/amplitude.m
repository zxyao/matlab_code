function [beta_tt, beta_rr] = amplitude(L, theta_t, theta_r)

beta_tt = zeros(L,1); beta_rr = zeros(L,1);
% theta_t = -theta_t + rho*lambda_t; theta_t = theta_t'*diag(qtt);
% theta_r = -theta_r + rho*lambda_r; theta_r = theta_r'*diag(qrr);

for n = 1:L
    a = abs(theta_t(n,n))*cos(angle(theta_t(n,n)));
    b = abs(theta_r(n,n))*cos(angle(theta_r(n,n)));
    
    phi = sign(b) * acos(a / sqrt(a^2 + b^2));

    if phi >= -pi && phi <-1/2*pi
        w = -1/2*pi - phi;
    elseif phi >= -1/2*pi && phi <= 1/4*pi
        w = 0;
    else
        w = 1/2*pi;
    end
    beta_tt(n) = sin(w);
    beta_rr(n) = cos(w);
w
end
end