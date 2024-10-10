function [Alp] = caculate_alpha(rho,P_m1,P_t1,Theta,M,N,gamma_const,H_mR,H_RB,H_nB,H_Rn,H_mB,H_nn,W_sigma,L,t_const11,t_const22,s_m1)
    cvx_begin quiet
    variable Alp(M,N) nonnegative
    expression R_1(M,1)
    expression R_2(M,1)
    expression R(M,1)
    expression sum_h(M,1)
    expression sum_left(N,1)
    expression sum_right(N,1)
    expression e1(M,1)
    expression e2(M,1)
    expression S_m(M,1)
    expression sum_Pt(M,1)
    expression P_t(N,1)
    expression P_m(M,1)
    %%
    lambda=0.7;
    Ab_Q=0.8;
    P_t=P_t1;
    P_m=P_m1;
    t_const1=t_const11;
    t_const2=t_const22;
    s_m=s_m1;
    %% caculate alpha with known theta
     hnB=H_nB';
     hmb=H_mB';
     hmb=sum(hmb,2);
     hrb=sum(H_RB,1);
     for i=1:N
         hnb(i)=sum(H_nB(:,i));
     end
    for i=1:M
        for j=1:N
             sum_Pt(i)=sum_Pt(i)+Alp(i,j)*P_t(j)*abs(hnb(j))^2;
        end
        sum_Pt(i)=sum_Pt(i)+W_sigma;
        R_1(i)=log(sum_Pt(i))./log(2);
        sum_h(i)=abs(hmb(i)+hrb*Theta*H_mR(:,i)).^2;
        R_2(i)=log(P_m(i)*sum_h(i)+sum_Pt(i)+W_sigma)./log(2);
%         R(i)=R_1(i)-R_2(i);  %向量
    end
    %对R_1展开
%     R1_prime=sp.diff(R_1,alpha); %对alpha求一阶导数
%     %计算展开式
%     taylor_expanssion=R_1.sub(alpha,alp)+R1_prime.subs(alpha,alp)*(alpha-alp);
%     R=taylor_expanssion-R2;
    %% 对约束进行处理
    hRn=H_Rn';
    for i=1:N
        for j=1:M
            sum_left(i)=sum_left(i)+Alp(j,i)*P_m(j);
        end
        sum_left(i)=sum_left(i)+W_sigma;
    end
    for i=1:N
        for j=1:M
            sum_right(i)=1+abs(3*hRn(i,:)*Theta*H_mR(:,j))^2;
        end
    end
     %% 计算安全系数之和
    for i=1:M
        const=(1-lambda)*Ab_Q*rho(i)+lambda*Ab_Q;
        e1=exp((1-rho(i)).*t_const1(i));
        if rho(i)==0
            e2=1;
        else
            e2=exp(rho(i).*(s_m(i).*inv_pos(R(i))+t_const2(i)));
        end
        S_m(i)=log(e1+e2)./const;
    end
    P_t=P_t';
    minimize sum(S_m)
    subject to
        P_t.*sum_right>=sum_left.*gamma_const;
        0<=sum(Alp,1)<=1;
        0<=sum(Alp,1)<=1;
    cvx_end
end

