
function  [W, theta_t, theta_r, theta_tt, theta_rr] = algorithm_BCD(M,N,L, theta_t, theta_r, theta_tt, theta_rr, lambda_t, lambda_r, rho,H_RB,H_mB,H_nB,H_mR,H_Rn,H_mn,H_nn,gamma_const,alpha,rho_opt,s_m,t_const1,t_const2,P_m,P_t)
%The BCD algorithm for solve the augmented Lagrangian problem

sum_sm_pre = 0;
% phas_diff_all = [];
% sum_rate_all = [];

for i = 1:40
    % update auxiliary variables
%     [Rm, Sm] = update_weights( theta_t, theta_r,H_RB,H_mB,H_nB,h_mR,h_Rn,H_mn,H_nn,gamma_const,alpha,M,rho_opt,s_m,t_const1,t_const2);
    % update beamformers
%     [W] = update_W(para, omega, upsilon, theta_t, theta_r, G, h);
    % update original STAR coefficients
    H_Rn
    h_rn=H_Rn';hmb=H_mB';
    [theta_t, theta_r] = update_theta(L, N, theta_tt, theta_rr, lambda_t, lambda_r,rho,hmb,H_nB,H_mR,h_rn,H_mn,H_nn,gamma_const,alpha,M,rho_opt,s_m,t_const1,t_const2,P_m,P_t);
    % update auxiliary STAR coefficients
    [theta_tt, theta_rr] = update_theta_aux(M,N,L, theta_t, theta_r, theta_tt, theta_rr, lambda_t, lambda_r, rho);
    
%     [gamma, ~] = SINR(para, W, theta_t, theta_r, G, h);
%     sum_rate = sum(log2(1 + gamma));
%     disp(['Inner loop - ' num2str(i) ', Sum rate - ' num2str(sum_rate)]);
%     sum_rate_all = [sum_rate_all, sum_rate];
%     phas_diff_all = [phas_diff_all abs(angle(theta_r) - angle(theta_t))];
    Rm=caculate-rate(P_m,P_t,alpha,theta_r,W_sigma,M,N,H_nB,H_mB,h_RB,H_mR);
    sum_sm=sum(caculate_S(rho_opt,t_const1,t_const2,lambda,M,R,Ab_Q,s_m));
    % check convergence
    reduction = abs(sum_sm - sum_sm_pre) / sum_sm;
    if reduction < 1e-3 % convergence
        break; 
    end
    sum_sm_pre = sum_sm;

end


end

%% update auxiliary variables
function [Rm, Sm] = update_weights(P_m,P_t,alpha,Theta_r,W_sigma,M,N,H_nB,H_mB,h_RB,H_mR)
    [Rm] = caculate_rate(P_m,P_t,alpha,Theta_r,W_sigma,M,N,H_nB,H_mB,h_RB,H_mR);
    [Sm] = caculate_S(rho,t_const1,t_const2,lambda,M,R,Ab_Q,s_m);
end

%% update beamformers
% function [W] = update_W(para, omega, upsilon, theta_t, theta_r, G, h)
%     Theta_t = diag(theta_t); ht = G' * Theta_t' * h; 
%     Theta_r = diag(theta_r); hr = G' * Theta_r' * h; 
%     cvx_begin quiet
%         % optimization variables
%         variable W(para.M, para.K) complex
% 
%         % constraint
%         square_pos(norm(W,'fro')) <= para.Pt;
%         
%         % calculate objective function
%         obj = 0;
%         for k = 1:para.K
%             if k <= para.K/2
%                 hk = ht(:,k); 
%             else
%                 hk = hr(:,k);
%             end   
%             wk = W(:,k);
%             ek = abs(upsilon(k))^2 * (square_pos(norm(sqrtm(hk*hk')*W,'fro')) + 1) - 2*real( conj(upsilon(k))*hk'*wk ) + 1;
%             obj = obj + omega(k)*ek;
%         end
% 
% 
%         minimize(obj);
%     cvx_end   
% end

%% update original STAR coefficients
function [theta_t, theta_r] = update_theta(L, N, theta_tt, theta_rr, lambda_t, lambda_r,rho,hmb,H_nB,H_mR,h_rn,H_mn,H_nn,gamma_const,alpha,M,rho_opt,s_m,t_const1,t_const2,P_m1,P_t1)

    cvx_begin quiet
        % optimization variables
        variable theta_t(L, 1) complex
        variable theta_r(L, 1) complex
        expression P_t(N,1)
        expression P_m(M,1)
        expression fenzi(N,1)
        expression hnb(N,1)
        expression fenmu(N,1)
        expression R(M,1)
        expression SINR(M,1)
        expression fenmu(M,1)
        P_t=P_t1;
        P_m=P_m1;
        % constraint
        for i = 1:L
            theta_t(i)*conj(theta_t(i)) + theta_r(i)*conj(theta_r(i)) <= 1;
        end
        % calculate objective function
        Theta_t = diag(theta_t); 
        Theta_r = diag(theta_r); 
        
        %constraint gamma
        h_mn=sum(H_mn,1);

              for j=1:N
                  for i=1:M
                       fenzi(j)=P_t(j)*pow_pos(abs(H_nn(j)+h_rn(j,:)*Theta_t*H_mR(:,i)),2);
                  end
              end
       for i=1:N
           for j=1:M
               fenmu(i)=fenmu(i)+alpha(j,i).*P_m(j).*abs(h_mn(i).^2);
           end
       end
       for i=1:N
            gamma_const.*fenmu(i)<=fenzi(i);
       end
       
       %caculate 
       lambda=0.7;
       Ab_Q=0.8;

           for i=1:N
               hnb(i)=sum(H_nB(:,i));
           end
           for i=1:M
               for j=1:N
                   fenmu(i)=fenmu(i)+alpha(i,j)*P_t(j)*abs(hnb(j))^2;
               end
               fenmu(i)=fenmu(i)+W_sigma;
           end
          h_rb=sum(h_RB,1);
           h_mb=sum(hmb,2);
        for i=1:M
            SINR(i)=((P_m(i).*abs(h_mb(i)+h_rb*Theta_r*H_mR(:,i)).^2))./fenmu(i);
            R(i)=log2(1+SINR(i));
        end
       for i=1:M
            const=(1-lambda)*Ab_Q*rho_opt(i)+lambda*Ab_Q;
            e1=exp((1-rho_opt(i)).*t_const1(i));
            e2=exp(rho_opt(i).*(s_m(i).*inv_pos(R(i))+t_const2(i)));
            S_m(i)=log(e1+e2)./const;
        end
       obj=sum(S_m);
        
        penalty = sum_square_abs(theta_tt - theta_t + rho*lambda_t) + sum_square_abs(theta_rr - theta_r + rho*lambda_r);
        obj = obj + 1/(2*rho) * penalty;
        minimize(obj);
    cvx_end  

end

%% update auxiliary STAR coefficients
function [theta_tt, theta_rr] = update_theta_aux(M,N,L, theta_t, theta_r, theta_tt, theta_rr, lambda_t, lambda_r, rho)
beta_tt = abs(theta_tt); 
beta_rr = abs(theta_rr); 

% update phase shift
[qtt, qrr] = update_phase(M,N,L, theta_t, theta_r, beta_tt, beta_rr, lambda_t, lambda_r, rho);

% update amplitude
[beta_tt, beta_rr] = update_amplitude(M,N,L, theta_t, theta_r, qtt, qrr, lambda_t, lambda_r, rho);

% combine amplitudes ad phase shifts
theta_tt = diag(beta_tt)*qtt;
theta_rr = diag(beta_rr)*qrr;

end

%% update phase shifts
function [qtt, qrr] = update_phase(M,N,L, theta_t, theta_r, beta_tt, beta_rr, lambda_t, lambda_r, rho)
qtt = zeros(L,1); qrr = zeros(L,1);

theta_t = -theta_t + rho*lambda_t; theta_t = theta_t'*diag(beta_tt);
theta_r = -theta_r + rho*lambda_r; theta_r = theta_r'*diag(beta_rr);

for n = L
    
    phi_p = theta_t(n) + 1i*theta_r(n);
    phi_m = theta_t(n) - 1i*theta_r(n);

    qtt_1 = exp(1i*(pi - angle(phi_p))); qrr_1 =  1i*qtt_1;
    qtt_2 = exp(1i*(pi - angle(phi_m))); qrr_2 =  -1i*qtt_2;


    o1 = real(theta_t(n)*qtt_1) + real(theta_r(n)*qrr_1);
    o2 = real(theta_t(n)*qtt_2) + real(theta_r(n)*qrr_2);

    if o1 < o2
        qtt(n) = qtt_1; 
        qrr(n) = qrr_1;
    else
        qtt(n) = qtt_2; 
        qrr(n) = qrr_2;
    end
end
qtt
end

%% update amplitudes
function [beta_tt, beta_rr] = update_amplitude(M,N,L, theta_t, theta_r, qtt, qrr, lambda_t, lambda_r, rho)

beta_tt = zeros(L,1); beta_rr = zeros(L,1);
theta_t = -theta_t + rho*lambda_t; theta_t = theta_t'*diag(qtt);
theta_r = -theta_r + rho*lambda_r; theta_r = theta_r'*diag(qrr);

for n = 1:L
    a = abs(theta_t(n))*cos(angle(theta_t(n)));
    b = abs(theta_r(n))*cos(angle(theta_r(n)));
    
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

end
beta_tt
end