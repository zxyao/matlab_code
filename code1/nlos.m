function [loss]=nlos(rho,dist,a)
    loss=sqrt(rho*dist^(-a));
end
