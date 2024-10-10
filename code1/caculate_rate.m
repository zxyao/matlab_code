function [R] = caculate_rate(P_m,P_t,alpha,Theta,W_sigma,M,N,H_nB,H_mB,h_RB,H_mR)
    R=zeros(M,1);
    SINR=zeros(M,1);
    fenmu=zeros(M,1);
       hnb=zeros(N,1);
       for i=1:N
           hnb(i)=sum(H_nB(:,i));
       end
       for i=1:M
           for j=1:N
               fenmu(i)=fenmu(i)+alpha(i,j)*P_t(j)*abs(hnb(j))^2;
           end
           fenmu(i)=fenmu(i)+W_sigma;
       end
      h_rb=sum(h_RB,1);
       hmb=H_mB';
       h_mb=sum(hmb,2);
    for i=1:M
        SINR(i)=((P_m(i).*abs(h_mb(i)+h_rb*Theta*H_mR(:,i)).^2))./fenmu(i);
        R(i)=log2(1+SINR(i));
    end
end