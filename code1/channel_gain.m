function [h_mR,h_Rn,h_RB,h_nn,h_mB,h_nB,h_mn]=channel_gain(Nx,Ny,M,N,K,CV_points,V2V_points,V2V_dist,ris_pos,bs_pos)
    %% parameter settings
    % path loss exponent
    a_bm=3;
    a_rm=2.2;
    a_br=2.5;
    a_Rn=2.3;

    %Rician factors
    k_rician=3;
    k_const1=k_rician/(1+k_rician);
    k_const2=1/(1+k_rician);

    rho=10^(-20/10); %-20dB
    %load dist
    % CV_points=load('CV_location.mat');
    % V2V_points=load('V2V_location.mat');
    % V2V_dist=load('V2V_dist.mat');

    wavelen=0.001;
    L=Nx*Ny;
    % channel gain
 
    %% CV->RICS
    h_mR=zeros(L,1,M);
    dist_mR=zeros(M,1);
    CV_location=CV_points;
    V2V_location=V2V_points;
    for i=1:M
        dist_mR(i)=pdist([CV_points(i,:);ris_pos]);
    end
    for i=1:M
        h_mR(:,:,i)=nlos(rho,dist_mR(i),a_rm)*(sqrt(k_const1)*mR_los(Nx,Ny,M,ris_pos,CV_location,dist_mR,wavelen)+sqrt(k_const2)*Rayleigh_model());
    end

    %% RICS->V2V
    h_Rn=zeros(L,1,N);
    dist_Rn=zeros(N,1);
    for i=1:5
        dist_Rn(i)=pdist([ris_pos;V2V_points(i,1:3)]);
    end
    save('dist_Rn.mat','dist_Rn');
    for i=1:N
        h_Rn(:,:,i)=nlos(rho,dist_Rn(i),a_Rn)*(sqrt(k_const1)*Rn_los(Nx,Ny,N,ris_pos,V2V_location,dist_Rn,wavelen)+sqrt(k_const2)*Rayleigh_model());
    end

    %% RICS->BS NLos+Los
    h_RB=zeros(K,L);
    dist_RB=pdist([bs_pos;ris_pos]);
    for i=1:K
       h_RB(i,:)=nlos(rho,dist_RB,a_br)*(sqrt(k_const1)*RB_los(Nx,Ny,ris_pos,bs_pos,wavelen,dist_RB)+sqrt(k_const2)*Rayleigh_model());
    end
  

    %% V2V NLos
    h_nn=zeros(N,1);
    for i=1:N
        h_nn(i)=nlos(rho,V2V_dist(i),a_Rn)*Rayleigh_model();
    end

    %% direct link between BS and CVs Los+NLos
    h_mB=zeros(K,1,M);
    dist_mB=zeros(M,1);
    for i=1:M
        dist_mB(i)=pdist([bs_pos;CV_points(i,1:3)]);
    end
    for i=1:M
        h_mB(:,:,i)=nlos(rho,dist_mB(i),a_bm)*(sqrt(k_const1)*mB_los(M,bs_pos,CV_location,dist_mB,wavelen,K)+sqrt(k_const2)*Rayleigh_model());
    end

    %% BS -> V2V(Tx) Los+NLos
    h_nB=zeros(K,1,N);
    dist_nB=zeros(N,1);
    for i=1:N
        dist_nB(i)=pdist([bs_pos;V2V_points(i,4:6)]);
    end
    for i=1:N
        h_nB(:,:,i)=nlos(rho,dist_nB(i),a_bm)*(sqrt(k_const1)*nB_los(N,bs_pos,CV_location,dist_nB,wavelen,K)+sqrt(k_const2)*Rayleigh_model());
    end

    %% V2V(Rx)-CV NLos
    h_mn=zeros(M,N);
    dist_mn=zeros(M,N);
    for i=1:M
        for j=1:N
            dist_mn(i,j)=pdist([CV_points(i,1:3);V2V_points(j,1:3)]);
        end
    end
    for i=1:M
        for j=1:N
            h_mn(i,j)=nlos(rho,dist_mn(i,j),a_Rn)*Rayleigh_model();
        end
    end

end