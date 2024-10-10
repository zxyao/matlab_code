function [hlos_RB]=RB_los(Nx,Ny,ris_pos,bs_pos,wavelen,dist_RB)
    %RICS is planar and modeled according to UPA and the base station antenna is modeled according to ULA
    % Calculation of antenna response
    projection_xie=sqrt((ris_pos(1)-bs_pos(1))^2+(ris_pos(2)-bs_pos(2))^2);% First calculate the hypotenuse projected onto the xoy plane
    xie=sqrt(projection_xie+(ris_pos(3)-bs_pos(3))^2);% Calculate the length of the hypotenuse in space
    sin_theta=(ris_pos(2)-bs_pos(2))/projection_xie;
    sin_phi=(ris_pos(3)-bs_pos(3))/xie;
    cos_phi=projection_xie/xie;
    at_phi_theta=array_response(sin_theta,sin_phi,cos_phi,Nx,Ny);
    
    hlos_RB=exp(-1i*2*pi*dist_RB*wavelen^(-1)).*at_phi_theta;
end