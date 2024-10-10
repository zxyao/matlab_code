function [h] = Rayleigth_model()
%% Complex Gaussian random variables with h zero mean and unit variance - the Rayleigh model
    %s=sqrt(var/2)*(randn(1,K) +j*randn(1,K)) 
    %where s denotes the complex Gaussian matrix, var denotes the power (i.e., variance), and K denotes the number of samples (var is 1 in this example)
    h = sqrt(1/2)*(randn(1,1)+1i*randn(1,1));
end