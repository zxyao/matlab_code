function [h_mR,h_Rn,h_RB,h_nn,h_mB,h_nB,h_mn]=init(L,Nx,Ny,M,N,K,CV_points,V2V_points,V2V_dist,ris_pos,bs_pos,lambda,Ab_Q,IterMax,delta)
    %% gain the location of cars and set parameters
    % initialization parameters

    %% compute the pathloss according to the cars' locations and generate the channel coefficients
    [h_mR,h_Rn,h_RB,h_nn,h_mB,h_nB,h_mn]=channel_gain(Nx,Ny,M,N,K,CV_points,V2V_points,V2V_dist,ris_pos,bs_pos);

    %% update parameters alternately
    %% Quadratic
    delta=1e-6;
    IterMax=20;
    r_=(1-lambda)*Ab_Q;
    r_c=lambda*Ab_Q;
    %% alternate optimization of transmission power P
    %% proposed_Algorithm
end
