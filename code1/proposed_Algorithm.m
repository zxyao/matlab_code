function [s,R_m,Rm,SINR,sum_capacity,sum_Sm,rho_opt,P_m,P_t,alpha,Theta,Theta1,W_sigma]=proposed_Algorithm(M,N,K,L,h_mR,h_Rn,h_RB,h_nn,h_mB,h_nB,h_mn,theta_init,alpha,lambda,Ab_Q,s_m,P_t,t_const1,t_const2,P_max,Psi_solution)
    %% Initialize rho and set index r
    s=zeros(M,1);
    s_last=zeros(M,1);
    rate=zeros(M,1);
    r=1;
    eqsi=0.001;
    %% Initialization rate R_m
    R_m=zeros(M,1);
    sum_capacity=zeros(14,1);
    sum_Sm=zeros(14,1);
    for i=1:M
        R_m(i)=10;
    end
    %Initialization omiga
    omiga=8;
    W=1;
    gamma_th=2.06; %bps/Hz
    P_outage=0.01;
    gamma_const=(gamma_th+1/omiga*log(1/P_outage-1));
    
    % Realigning channels into 2D arrays
    HmB=reshape(h_mB,K,M);
    HmR=reshape(h_mR,L,M);
    HnB=reshape(h_nB,K,N);
    HRn=reshape(h_Rn,L,N);
    
    H_mB=reshape(h_mB,K,M).*10^4;
    H_mR=reshape(h_mR,L,M).*10^4;
    H_nB=reshape(h_nB,K,N).*10^4;
    H_Rn=reshape(h_Rn,L,N).*10^4;
    H_RB=h_RB.*10^2;
    H_mn=h_mn.*10^4;
    H_nn=h_nn.*10^4;
    Theta=diag(theta_init); % Initial Phase Shift 
    save('Hmn.mat','H_mn');
    %参数设置
    W_sigma=dBm_W(-110); %-110dBm for % interference
    %% Quadratic
    delta=1e-6;
    IterMax=20;
    r_=(1-lambda)*Ab_Q;
    r_c=lambda*Ab_Q;
    
    %% 交替迭代s
%     r<15
    Theta1=0;
    while (sum(s)-sum(s_last))/sum(s_last)>eqsi || r==1
        %for given P,alpha,phi caculate rho
        [isFeasible,obj_opt,rho_opt,F_log,F_i] = LFP_Quadratic(r_,r_c,delta,IterMax,t_const1,t_const2,R_m,M,s_m);
%         rho_opt=randi([0,1],1,M);
        save('rho.mat','rho_opt');
        %for given rho,alpha,phi
        [P_m,P_t]=cvx_power(alpha,rho_opt,M,N,P_t,Theta,H_mB,H_nB,H_mR,H_RB,H_nn,gamma_const,L,H_Rn,t_const1,t_const2,s_m,P_max);
        P_m
        P_t
        %for given rho,alpha,phi
        alpha=caculate_alpha(rho_opt,P_m,P_t,Theta,M,N,gamma_const,H_mR,H_RB,H_nB,H_Rn,H_mB,H_nn,W_sigma,L,t_const1,t_const2,s_m);
        alpha
        
        %for given rho,P,alpha
        [Theta,Theta1]=SDR(N,P_m,P_t,L,H_RB,H_mB,H_nB,h_mR,h_Rn,H_mn,H_nn,gamma_const,alpha,M,rho_opt,s_m,t_const1,t_const2);
        R_m=caculate_rate(P_m,P_t,alpha,Theta,W_sigma,M,N,H_nB,H_mB,H_RB,H_mR);
        
        sum_capacity(r)=sum(R_m);
        
        %caculate result
        r=r+1;
%         R_m=caculate_rate(P_m,P_t,alpha,Theta,W_sigma,M,N,H_nB,H_mB,h_RB,H_mR);
        s_last=s;
        s=caculate_S(rho_opt,t_const1,t_const2,lambda,M,R_m,Ab_Q,s_m);
        sum_Sm(r)=sum(s);
    end
    [beta_tt, beta_rr] = amplitude(L, Theta1, Theta);
    Theta1=Theta1.*diag(beta_tt).*diag(Psi_solution);
    [Rm,SINR]=caculate_rate_V2V(P_m,P_t,alpha,Theta1,W_sigma,M,N,H_mn,H_nn,H_Rn,H_mR);
end