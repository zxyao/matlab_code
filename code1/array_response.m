function [y] = array_response(sin_theta,sin_phi,cos_phi,Nx,Ny)
%RICS is deployed in the xoz plane
    y=zeros(Nx*Ny,1);
        for m=0:Nx-1
            for n=0:Ny-1
                y(m*Ny+n+1)=exp( -1i * pi * (n*sin_theta*sin_phi + m*cos_phi) );
            end
        end
end

