function [s,z] = noise(sys,n)

s(1)    = randn(1);
w       = sqrt(sys.W)*randn(1,n);
if sys.L==0
v       = sqrt(sys.V)*randn(1,n);
else
    v=w;
end

e = sqrt(sys.Psi)*randn(1,n);
shat(1) = 0;
for i=1:n
    s(i+1)      = sys.F*s(i) + w(i);
    z(i)        = sys.H*s(i) + v(i);
    shat(i+1) = sys.F*shat(i) + sys.Kp*(z(i) - sys.H*shat(i));
end

% mean((z - sys.H*shat(1:end-1)).^2)
% 
% disp(z-sys.H*s(1:end-1));