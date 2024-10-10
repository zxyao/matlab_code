% using gradient descent method
function [Psi_solution, Psi_history, obj_history] = gradient_descent(Psi_initial, h_mn, h1, h2)
    learning_rate = 0.01;
    tol = 1e-6;
    max_iter = 1000;
    l=30;
    
    Psi = Psi_initial
    Psi_history = [];
    obj_history = [];
    
    for i = 1:max_iter
        grad = gradient(Psi, h_mn, h1, h2);
        Psi = Psi - learning_rate * grad;
        
%         Psi(Psi < 1) = 1;
         Psi(Psi <= 0) = 1e-6;
        % checking for convergence
        if norm(grad) < tol
            break;
        end
        
         % collecting the current state
        Psi_history(end+1,:) = Psi;
        obj_history(end+1,:) = objective_function(Psi, h_mn, h1, h2);
    end
    
    Psi_solution = Psi;
end