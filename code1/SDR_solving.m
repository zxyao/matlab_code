function [v, lower_bound] = SDR_solving(N,P_m,P_t,L,h_RB,H_mB,H_nB,h_mR,h_Rn,h_mn,h_nn,gamma_const,alpha,M,rho1,s_m1,t_const11,t_const22)
    W_sigma=dBm_W(-110);
    LL = 1000; 
    %% Construct auxiliary matrices H, V
    for i=1:L
        h_rb(i)=sum(h_RB(:,i));
    end
    h=zeros(L,1,M);
    H=zeros(L+1,L+1,M);
    for i=1:M
        h_mb(i)=sum(H_mB(:,i));
    end
    for k=1:M
        h(:,:,k)=h_rb'.*h_mR(:,:,k);
        H(:,:,k)=P_m(k).*[h(:,:,k)*h(:,:,k)',h(:,:,k)*h_mb(k)'; h_mb(k)*h(:,:,k)',h_mb(k)'*h_mb(k)]; %辅助变量H
    end
     %% 
    ans_h=zeros(L+1,L+1);
    for k=1:M
        ans_h=ans_h+H(:,:,k);
    end
    
    %% Construct the auxiliary matrix of constraints
    for i=1:N
        h1=zeros(L,1,N);
        W=zeros(L+1,L+1,N);
        for j=1:M 
            h1(:,:,i)=h_Rn(:,:,i).*h_mR(:,:,j);
            W(:,:,i)=P_t(i).*[h1(:,:,i)*h1(:,:,i)',h1(:,:,i)*h_nn(i)'; h_nn(1)*h1(:,:,i)',h_nn(i)*h_nn(i)']; %辅助变量W
        end
    end
    ans_h1=zeros(L+1,L+1);
    for k=1:N
        ans_h1=ans_h1+W(:,:,k);
    end
    %% Calculate the constant part of the constraint
    sum_ph=zeros(N,1);
    for i=1:N
        for j=1:M
            sum_ph(i)=sum_ph(i)+alpha(j,i)*P_m(j)*abs(h_mn(j,i))^2;
        end
    end
    sum_right=sum_ph.*gamma_const;
    
    cvx_begin sdp quiet
        variable V_r(L+1,L+1) hermitian
        variable V_t(L+1,L+1) hermitian
        expression sum_left(N,1)
        expression fenmu(M,1)
        expression hnb(N,1)
        expression S_m(M,1)
        expression const
        expression e1
        expression e2
        expression t_const1(M,1)
        expression t_const2(M,1)
        expression s_m(M,1)
        expression rho(M,1)
        expression H1(L+1,L+1,M)
        expression W1(L+1,L+1,N)
        expression sumH(L+1,L+1)
        expression sumW(L+1,L+1)
        expression R(M,1)
        %% 
        lambda=0.7;  
        Ab_Q=0.8;
        t_const1=t_const11;
        t_const2=t_const22;
        s_m=s_m1;
        rho=rho1;
        H1=H;
        W1=W;
        sumH=ans_h;
        sumW=ans_h1;
        for i=1:N
            for j=1:M
                sum_left(i)=sum_left(i)+real(trace(V_t*ans_h1)); 
            end
        end
       
       %% Calculation rate
       for i=1:N
           hnb(i)=sum(H_nB(:,i));
       end
       for i=1:M
           for j=1:N
               fenmu(i)=fenmu(i)+alpha(i,j)*P_t(j)*abs(hnb(j))^2;
           end
           fenmu(i)=fenmu(i)+W_sigma;
       end
       for i=1:M
           R(i)=log(1+real(trace(V_r*H1(:,:,i)))./fenmu(i))./log(2);
       end
       
       %% calculate sum sm
       
       for i=1:M
           const=(1-lambda)*Ab_Q*rho(i)+lambda*Ab_Q;
           e1=exp((1-rho(i)).*t_const1(i));
           e2=exp(rho(i).*(s_m(i).*inv_pos(R(i))+t_const2(i)));
           S_m(i)=log(e1+e2)./const;
       end
       
        minimize sum(S_m);
        subject to
%            diag(V_t)+diag(V_r)==1.8;
            for i=1:L
                V_t(i,i)==1.3*(1-V_r(i,i));
            end
%             diag(V_r)>=0;
%             diag(V_t)>=0;
            V_r>=0;
            V_t>=0;
            V_r(L+1,L+1)==1;
            sum_left>=sum_right;
    cvx_end
    
    lower_bound = cvx_optval;
      % guass
        %% method 1
        max_F = 0;
        max_v = 0;
        sum_Sm=0;
        [U, Sigma] = eig(V_r);
        for l = 1 : LL
            r = sqrt(2) / 2 * (randn(L+1, 1) + 1j * randn(L+1, 1));
            vr = U * Sigma^(0.5) * r;
            for i=1:M
                sum_Sm=sum_Sm+vr'*H(i)*vr;
            end
            if sum_Sm > max_F
                max_v = vr;
                max_F = sum_Sm;
            end
        end

        vr = exp(1j * angle(max_v / max_v(end)));
        vr = vr(1 : L);
        vr=diag(vr);
    %     v' * h(:,:,k)*h(:,:,k)' * v;
    V_t
    V_r
    max_F = 0;
        max_v = 0;
        sum_Sm=0;
        [U, Sigma] = eig(V_t);
        for l = 1 : LL
            r = sqrt(2) / 2 * (randn(L+1, 1) + 1j * randn(L+1, 1));
            vt = U * Sigma^(0.5) * r;
            for i=1:M
                sum_Sm=sum_Sm+vt'*H(i)*vt;
            end
            if sum_Sm > max_F
                max_v = vt;
                max_F = sum_Sm;
            end
        end

        vt = exp(1j * angle(max_v / max_v(end)));
        vt = vt(1 : L);
        vt=diag(vt);

end