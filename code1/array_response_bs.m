function [y] = array_response_bs(sin_theta,sin_phi,cos_phi,K)
%RICS is deployed in the xoz plane
    y=zeros(K,1);
        for n=1:K
            y(n)=exp( -1i * pi * (n*sin_theta*sin_phi + n*cos_phi) );
        end
end

