
function result = objective_function(Psi, h_mn, h1, h2)
    h_mn_re = real(h_mn);
    h_mn_im = imag(h_mn);
    h1_re = real(h1);
    h1_im = imag(h1);
    h2_re = real(h2);
    h2_im = imag(h2);
    M=10;result=0;
    for i=1:M
        result = result + (h_mn_re(i) + Psi * h1_re').^2 + (h_mn_im(i) + Psi * h2_im').^2;
    end
end
