function R = chap_eigen_scan(P,D,Vx_vec,E,Afun)
% Sweep forward speed, compute generalized eigenvalues

n=length(Vx_vec);
lambda = zeros(n,10);
freq_hz = zeros(n,10);
zeta = zeros(n,10);

for i=1:n
    Vx=Vx_vec(i);
    A = Afun(Vx,0);
    [V,Dlam] = eig(A,E);
    lam = diag(Dlam);
    lambda(i,:) = lam;
    freq_hz(i,:) = abs(imag(lam))/(2*pi);
    zeta(i,:) = -real(lam)./sqrt(real(lam).^2+imag(lam).^2);
end

R.v = Vx_vec; R.lambda=lambda;
R.freq_hz=freq_hz; R.zeta=zeta;
end
