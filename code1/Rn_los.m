function [hlos_Rn] = Rn_los(Nx,Ny,N,ris_pos,V2V_location,dist_Rn,wavelen)
    %RICS is planar and modeled according to UPA and the base station antenna is modeled according to ULA
    for k=1:N
        h_nR=exp(-1i*2*pi*dist_Rn(k)*wavelen^(-1));
        % Calculation of antenna response
        projection_xie=sqrt((V2V_location(k,1)-ris_pos(1)).^2+(V2V_location(k,2)-ris_pos(2)).^2);% First calculate the hypotenuse projected onto the xoy plane
        xie=sqrt(projection_xie+(ris_pos(3)-V2V_location(k,3)).^2);% Calculate the length of the hypotenuse in space
        sin_theta=(ris_pos(2)-V2V_location(k,2))./projection_xie;
        sin_phi=(ris_pos(3)-V2V_location(k,3))./xie;
        cos_phi=projection_xie./xie;
        ar_phi_theta=array_response(sin_theta,sin_phi,cos_phi,Nx,Ny);

        hlos_Rn=h_nR.*ar_phi_theta;
    end
end

