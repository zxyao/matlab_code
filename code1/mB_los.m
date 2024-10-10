function [hlos_mB] = mB_los(M,bs_pos,CV_location,dist_mB,wavelen,K)
    %RICS is planar and modeled according to UPA and the base station antenna is modeled according to ULA
    for k=1:M
        h_mB=exp(-1i*2*pi*dist_mB(k)*wavelen^(-1));
        % Calculation of antenna response
        projection_xie=sqrt((CV_location(k,1)-bs_pos(1)).^2+(CV_location(k,2)-bs_pos(2)).^2);% First calculate the hypotenuse projected onto the xoy plane
        xie=sqrt(projection_xie+(bs_pos(3)-CV_location(k,3)).^2);% Calculate the length of the hypotenuse in space
        sin_theta=(bs_pos(2)-CV_location(k,2))./projection_xie;
        sin_phi=(bs_pos(3)-CV_location(k,3))./xie;
        cos_phi=projection_xie./xie;
        ar_phi_theta=array_response_bs(sin_theta,sin_phi,cos_phi,K);

        hlos_mB=h_mB.*ar_phi_theta;
    end
end

