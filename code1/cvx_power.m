function [P_m,P_t]=cvx_power(alpha,rho1,M,N,Pt,Theta,H_mB,H_nB,H_mR,H_RB,H_nn,gamma_const,L,H_Rn,t_const11,t_const22,s_m1,P_max)
%     Pt_max=dBm_W(23);
%     Pm_max=dBm_W(23);
%     sum_Pt=zeros(M,1);
%     R_1=zeros(M,1);
%     R_2=zeros(M,1);
%     R=zeros(M,1);
%     sum_h=zeros(M,1);
    %% CVX求解
%    cvx_begin quiet
%     variable P_t(N,1) nonnegative
%     variable P_m(M,1) nonnegative
%     variable sm
%     expression sum_Pt(M,1)
%     expression sumPt(M,1)
%     expression R_1(M,1) 
%     expression R1(M,1) 
%     expression R_2(M,1) 
%     expression R(M,1) 
%     expression sum_h(M,1)
%     expression alpha1(M,1)
%     expression R1_expand(M,1)
%     expression R_approx(M,1)
%     expression sum_left(N,1)
%     expression S_m(M,1)
%     expression Pt_max(N,1)
%     expression Pm_max(M,1)
%     expression gamma_const1(N,1)
%     expression const(M,1)
%     expression t_const1(M,1)
%     expression t_const2(M,1)
%     expression s_m(M,1)
%     expression rho(M,1)
%     expression e1(M,1)
%     expression e2(M,1)
%     %%   
%     W_sigma=dBm_W(-110); %干扰项
%     lambda=0.7;
%     Ab_Q=0.8;
%     t_const1=t_const11;
%     t_const2=t_const22;
%     s_m=s_m1;
%     rho=rho1;
%     for i=1:N
%         Pt_max(i)=dBm_W(23);
%     end
%     for i=1:M
%         Pm_max(i)=P_max;
%     end
%     for i=1:N
%         gamma_const1(i)=gamma_const;
%     end
%     %% 计算R_m
%      hnB=H_nB';
%      hmb=H_mB';
%      hmb=sum(hmb,2);
%      hrb=sum(H_RB,1);
%      for i=1:N
%          hnb(i)=sum(H_nB(:,i));
%      end
%     for i=1:M
%         for j=1:N
%              sum_Pt(i)=sum_Pt(i)+alpha(i,j)*P_t(j)*abs(hnb(j))^2;
%         end
%         sum_Pt(i)=sum_Pt(i)+W_sigma;
% %         R_1(i)=log(sum_Pt(i))./log(2);
%         sum_h(i)=abs(hmb(i)+hrb*Theta*H_mR(:,i)).^2;
%         R_2(i)=log(P_m(i)*sum_h(i)+sum_Pt(i)+W_sigma)./log(2);
% %         R(i)=R_1(i)-R_2(i);  %向量
% %         R(i)=log(1+(P_m(i).*sum_h(i))*inv_pos(sum_Pt(i)))./log(2);
%     end
%     %在给定点Pt处展开的导数
%     for i=1:M
%         for j=1:N
%              sumPt(i)=sumPt(i)+alpha(i,j)*Pt(j)*abs(hnb(j))^2;
%         end
%         R1(i)=log(sumPt(i)+W_sigma)/log(2);
%        % R(i)=log(1+P_m(i).*sum_h(i))*inv_pos(R1(i))./log(2);
%     end
%     
%     %对R_1展开 对P_t求导数
%     alpha1=sum(alpha,2);
%     R1_expand=alpha1./(R1.*log(2));
%     
%     %计算展开式
%     R_approx = R1 + R1_expand .* (sum(P_t) - sum(Pt));
%     R=R_2-R_approx;
    
%     %% 处理约束
%     hRn=H_Rn';
%     for i=1:N
%         for j=1:M
%             sum_left(i)=sum_left(i)+(alpha(j,i)*P_m(j)*(1+abs(3*hRn(i,:)*Theta*H_mR(:,j))^2))*10^(-3);
%         end
% %         sum_left(i)=sum_left(i)+W_sigma;
%     end
    %% 计算安全系数S_m
%     for i=1:M
%         const=(1-lambda)*Ab_Q*rho(i)+lambda*Ab_Q;
%         e1=exp((1-rho(i)).*t_const1(i));
%         e2=exp(rho(i).*(s_m(i).*inv_pos(R(i))+t_const2(i)));
%         S_m(i)=log(e1+e2)./const;
%     end
%         minimize sum(S_m)
%         subject to 
%             0<=P_t<=Pt_max;
%             0<=P_m<=Pm_max;       
%             (gamma_const1.*sum_left)<=P_t;
%     cvx_end
    for i=1:N
        P_t(i)=dBm_W(25);
    end
    for i=1:M
        P_m(i)=P_max;
    end
end