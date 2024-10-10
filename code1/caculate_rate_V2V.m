% function [R,SINR] = caculate_rate_V2V(P_m,P_t,alpha,Theta,W_sigma,M,N,H_mn,h_nn,H_Rn,H_mR)
%     R=zeros(N,1);
%     SINR=zeros(N,1);
%     fenzi=zeros(N,1);
%        hnb=zeros(N,1);  
%        h_mn=sum(H_mn,1);
%       h_rn=H_Rn';
% %        for i=1:M
%            for j=1:N
%                for i=1:M
%                     fenzi(j)=P_t(j)*abs(h_nn(j)+h_rn(j,:)*Theta*H_mR(:,i))^2;
%                end
%            end
% %            fenzi(i)=fenzi(i)+W_sigma;
% %        end
%       fenmu=zeros(N,1); 
%     for i=1:N
%         for j=1:M
%             fenmu(i)=fenmu(i)+alpha(j,i).*P_m(j).*abs(h_mn(i).^2);
%         end
%     end
%     for i=1:N
%         SINR(i)=fenzi(i)./fenmu(i);
%         R(i)=log2(1+SINR(i));
%     end
% %     SINR
% end
function [R,SINR] = caculate_rate_V2V(P_m,P_t,alpha,Theta,W_sigma,M,N,H_mn,h_nn,H_Rn,H_mR)
    R=zeros(N,1);
    SINR=zeros(N,1);
    fenzi=zeros(N,1);
       hnb=zeros(N,1);
%        for i=1:M
           for j=1:N
               fenzi(j)=P_t(j)*abs(h_nn(j))^2;
           end
%            fenzi(i)=fenzi(i)+W_sigma;
%        end
%       h_mn=sum(H_mn,1);
      h_rn=H_Rn';
      fenmu=zeros(N,1); 
    for i=1:N
        for j=1:M
            a=h_rn(i,:)*Theta*H_mR(:,j)*0.1;
            fenmu(i)=fenmu(i)+alpha(j,i)*0.1.*P_m(j).*(abs(H_mn(j,i)+a).^2);
        end
    end
    for i=1:N
        SINR(i)=fenzi(i)./fenmu(i);
        R(i)=log2(1+SINR(i));
    end
%     SINR
end