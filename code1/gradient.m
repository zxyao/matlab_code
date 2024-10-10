% calculate the gradient of the objective function
function grad = gradient(Psi, h_mn, h1, h2)
    h_mn_re = real(h_mn);
    h_mn_im = imag(h_mn);
    h1_re = real(h1);
    h1_im = imag(h1);
    h2_re = real(h2);
    h2_im = imag(h2);
    M=10;
    L=30;
    grad=zeros(1,L);
    for i=1
        grad(1,:) = grad(1,:)+2 * (h_mn_re(i) + Psi * h1_re') * h1_re + 2 * (h_mn_im(i) + Psi * h2_im') *h2_im;
    end
end
