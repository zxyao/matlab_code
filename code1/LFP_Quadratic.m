function [isFeasible,obj_opt,rho_opt,F_log,F_i] = LFP_Quadratic(r_,r_c,delta,IterMax,t_const1,t_const2,R_m,M,s_m)
% Input parameter description:
% delta small positive number, convergence threshold
% IterMax Maximum number of iterations
% M dimension
% Output parameter description:
% isFeasible
% obj_opt Optimal objective function value
% tau_opt Optimal solution
% F_log Record the value of F for each iteration.
% Use CVX to solve the initial solution A Linear Programming Problem
cvx_begin quiet
    variable rho(M,1) nonnegative
    subject to
        0<=rho<=1;
cvx_end
% Initialization data
rho_i=rho;
F_i=(r_.*rho+r_c)./(max((1-rho).*t_const1,rho.*(s_m./R_m+t_const2)));
F_log=[F_i];
F_i=sum(F_i);
isFeasible=true;

for j=1:IterMax
    y=sum(sqrt(r_.*rho_i+r_c)./(max((1-rho_i).*t_const1(),rho_i.*(s_m./R_m+t_const2()))));
    cvx_clear
    cvx_begin
        variable rho(M,1) nonnegative
        maximize sum((2.*y.*sqrt(r_.*rho+r_c)-y.^2.*(max((1-rho).*t_const1,rho.*(s_m./R_m+t_const2))) ))
        
        subject to
            0<=rho<=1;
    cvx_end
    rho_i_1=rho;
    F_i_1=(r_.*rho_i_1+r_c)./(max((1-rho_i_1).*t_const1,rho_i_1.*(s_m./R_m+t_const2)));
    F_log=[F_log;F_i_1];
    
   % If not feasible or not bounded, just exit the loop.
   if strcmp(cvx_status,'Infeasible') || strcmp(cvx_status,'Unbounded') || isnan(cvx_optval) || isinf(cvx_optval)
        isFeasible = false;
        break
   end
   
   % Exit the loop if it converges, if not update Q
   if (sum(F_i_1) - sum(F_i)) <= delta
       break
   else
       rho_i = rho_i_1;
       F_i = F_i_1;
   end
   
end

obj_opt = F_i_1;
rho_opt = rho_i_1;
% rho_opt=randi([0,0],1,M);
end

