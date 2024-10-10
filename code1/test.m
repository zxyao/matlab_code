sum_Pt=zeros(M,1);
hnB=H_nB';
     hmb=H_mB';
     hmb=sum(hmb,2);
     for i=1:L
        hrb=h_RB(i,:);
     end
     for i=1:N
         hnb(i)=sum(H_nB(:,i));
     end
 alpha=init_alpha(M,N);
 P_t=zeros(N,1);
     array=randi([1,23],1,5);
     for i=1:N
         P_t(i)=dBm_W(array(i));
     end
      W_sigma=dBm_W(-110);
      Theta=diag(theta_init);
      P_m=zeros(M,1);
     array=randi([1,23],1,M);
     for i=1:M
         P_m(i)=dBm_W(array(i));
     end
     R_1=zeros(M,1);
     R_2=zeros(M,1);
     for i=1:M
        for j=1:N
             sum_Pt(i)=sum_Pt(i)+alpha(i,j)*P_t(j)*abs(hnb(j))^2;
        end
        sum_Pt(i)=sum_Pt(i)+W_sigma+1;
        R_1(i)=log(sum_Pt(i))./log(2);
        sum_h(i)=abs(hmb(i)+hrb*Theta*H_mR(:,i)).^2;
        R_2(i)=log(P_m(i)*sum_h(i)+sum_Pt(i)+W_sigma)./log(2);
%         R(i)=R_1(i)-R_2(i);  
%         R(i)=log(1+(P_m(i).*sum_h(i))*inv_pos(sum_Pt(i)))./log(2);
     end
      Pt=zeros(N,1);
     array=randi([1,10],1,N);
     for i=1:N
         Pt(i)=dBm_W(array(i));
     end
     sumPt=zeros(M,1);
     R1=zeros(M,1);
      % The derivative of the % expansion at a given point Pt
    for i=1:M
        for j=1:N
             sumPt(i)=sumPt(i)+alpha(i,j)*Pt(j)*abs(hnb(j))^2;
        end
        R1(i)=log(sumPt(i)+W_sigma)/log(2);
       % R(i)=log(1+P_m(i).*sum_h(i))*inv_pos(R1(i))./log(2);
    end
    R1_expand=zeros(M,1);
    alpha1=sum(alpha,2);
    aa=(R1.*log(2));
    R1_expand=alpha1./aa;
    R_approx = R1 + R1_expand .* (sum(P_t) - sum(Pt));
    R=R_2-R_approx;
    
sum_left=zeros(N,1);
hRn=H_Rn';
    for i=1:N
        for j=1:M
            sum_left(i)=sum_left(i)+(alpha(j,i)*P_m(j)*(1+abs(hRn(i,:)*Theta*H_mR(:,j))^2))*10^(-2);
        end
%         sum_left(i)=sum_left(i)+W_sigma;
    end