function [hlos_mR] = mR_los(Nx,Ny,M,ris_pos,CV_location,dist_mR,wavelen)
    %RICS is planar and modeled according to UPA and the base station antenna is modeled according to ULA
    for k=1:M
        h_mR=exp(-1i*2*pi*dist_mR(k)*wavelen^(-1));
        % Calculation of antenna response
        projection_xie=sqrt((CV_location(k,1)-ris_pos(1)).^2+(CV_location(k,2)-ris_pos(2)).^2);% First calculate the hypotenuse projected onto the xoy plane
        xie=sqrt(projection_xie+(ris_pos(3)-CV_location(k,3)).^2);% Calculate the length of the hypotenuse in space
        sin_theta=(ris_pos(2)-CV_location(k,2))./projection_xie;
        sin_phi=(ris_pos(3)-CV_location(k,3))./xie;
        cos_phi=projection_xie./xie;
        ar_phi_theta=array_response(sin_theta,sin_phi,cos_phi,Nx,Ny);

        hlos_mR=h_mR.*ar_phi_theta;
    end
end

