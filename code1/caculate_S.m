function [s] = caculate_S(rho,t_const1,t_const2,lambda,M,R,Ab_Q,s_m)
%caculate sum safety coefficient
    s=zeros(M,1);
    for i=1:M
        k=Ab_Q.*(lambda+rho(i).*(1-lambda))
        s(i)=k./max((1-rho(i))*t_const1(i),rho(i)*(s_m(i)/R(i)+t_const2(i)))
    end
end

