% % generate simulate settings
M=10;N=5;r=500;
[s_m,c_m,f_m,t_const1,t_const2,lambda,Ab_Q,bs_vue,ris_vue,bs_ris,wavelen,K,r]=parameter_settings(M,N);
[bs_pos,ris_pos,CV_points,V2V_points,V2V_dist]=scene(r,M);
% %RIS Number of elements is 30
L=30;
Nx=5;
Ny=6;
% Initialize phase shift
theta_init=exp(1j.*rand(Nx*Ny,1).*2.*pi);
alpha=init_alpha(M,N);
% 
delta=1e-6;
IterMax=20;
r_=(1-lambda)*Ab_Q;
r_c=lambda*Ab_Q;
P_max=dBm_W(23);
[h_mR,h_Rn,h_RB,h_nn,h_mB,h_nB,h_mn]=init(L,Nx,Ny,M,N,K,CV_points,V2V_points,V2V_dist,ris_pos,bs_pos,lambda,Ab_Q,IterMax,delta);
load('Psi.mat');
[s1,rate1,rate,SINR,sum_capacity,sum_Sm,rho_opt,P_m,P_t,alpha,Theta,Theta1] = proposed_Algorithm(M,N,K,L,h_mR,h_Rn,h_RB,h_nn,h_mB,h_nB,h_mn,theta_init,alpha,lambda,Ab_Q,s_m,P_t,t_const1,t_const2,P_max,Psi_solution);
