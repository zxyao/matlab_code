function [s_m,c_m,f_m,t_const1,t_const2,lambda,Ab_Q,bs_vue,ris_vue,bs_ris,wavelen,K,r]=parameter_settings(M,N)
    %% Initial scene radius design
    r=400;
    %% Initialize the number of CVs, V2Vs
%     M=10;
%     N=5;
    %% Task J latency
    c_m=zeros(M,1);
    f_m=zeros(M,1);
    s_m=zeros(M,1);
    F=50000;
    R=100; %bps/Hz(Mean value calculated)
    for i=1:M
        s_m(i)=randi([10,20]);
%         s_m(i)=10;
        c_m(i)=randi([5000,6000]);
        f_m(i)=randi([1000,5000]);
        
    end
    t_const1=c_m./f_m;
    t_const2=(c_m./F);

    %% Reasoning accuracy
    lambda=0.8;
    Ab_Q=0.9;
    Am_Q=lambda*Ab_Q;

    %% 卸载功率
    Poutage=0.01;
    R_th=2;

    %% 路径损耗指数
    bs_vue=3;
    ris_vue=2.2;
    bs_ris=2.5;

    %% RIS元素块数
    Nx=5;
    Ny=10;
    L=50;

    wavelen=0.001;

    %BS天线数
    K=32;
end