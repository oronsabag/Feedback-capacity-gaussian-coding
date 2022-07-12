clear all
close all
clc

beta = 0.5;
sys.F = beta;
sys.H = beta;
sys.G = 1;
sys.V= 1;
sys.W = 1;
sys.L = 1;
[sys.P,eig_P,K_Temp] = dare(sys.F',sys.H',sys.G*sys.W*sys.G',sys.V,sys.G*sys.L);
sys.Psi     = sys.H*sys.P*sys.H' + sys.V;
sys.Kp      = (sys.F*sys.P*sys.H' + sys.G*sys.L)*sys.Psi^(-1);
sys.Pow     = 1;

sys.CVX = CVX_cap(sys);
sys.dec.PsiY = sys.CVX.Pi + sys.H*sys.CVX.Sig*sys.H'+ 2*sys.CVX.Gamma*sys.H'+ sys.Psi;
sys.dec.KY = (sys.F*sys.CVX.Gamma + sys.F*sys.CVX.Sig*sys.H + 1)*inv(sys.dec.PsiY);

% load('sys.mat');
n = 50;

HAT_SIG(1)  = sys.Psi;
F_bar(1) = 1;
COV_Z1_0   = sys.Kp*sys.Psi*sys.Kp';
for i=1:n
    PSIY(i)      = (sys.CVX.Gamma/sys.CVX.Sig + sys.H)*HAT_SIG(i)*(sys.CVX.Gamma/sys.CVX.Sig + sys.H)' + 1;
    Ky(i)        = (sys.F*HAT_SIG(i)*(sys.CVX.Gamma/sys.CVX.Sig + sys.H) + 1) * inv(PSIY(i));
    HAT_SIG(i+1) = sys.F*HAT_SIG(i)*sys.F' + 1 - Ky(i)*PSIY(i)*Ky(i)';
    F_bar(i+1)   = (sys.F - sys.Kp*(sys.CVX.Gamma/sys.CVX.Sig + sys.H))*F_bar(i);
    H(i)         = F_bar(i)*(sys.CVX.Gamma/sys.CVX.Sig + sys.H);
    if i>1
        COV_Z1(i)  = Var_sm(COV_Z1(i-1),H(i),PSIY(i));
    else
        COV_Z1(i)  = Var_sm(COV_Z1_0,H(i),PSIY(i));
    end
end

EXP = 10000;
for j=1:EXP
    [s,z] = noise(sys,n);
    
    %Transmission
    sh(1)       = sys.Kp*z(1);
    shh(1)      = 0;
    shh2(1)      = 0;
    zp = 0;
    s(1) = randn(1);
    sh(1) = s(1);
    for i=1:n
        x(i)         = (sys.CVX.Gamma/sys.CVX.Sig)*(sh(i)-shh(i));
        y(i)         = x(i) + sys.H*sh(i) + z(i) - sys.H*sh(i);
        sh(i+1)      = shat_enc(sh(i),z(i),sys);
        e(i) = z(i) - sys.H*sh(i);
        ytilde(i) = (sys.CVX.Gamma/sys.CVX.Sig+sys.H)*(sh(i)-shh(i)) + e(i);
        shh(i+1)     = shat_dec(shh(i),Ky(i),y(i),sys);
        shh2(i+1)    = shat_dec(shh2(i),sys.dec.KY,y(i),sys);
    end
    xi(1) = sh(1);
    for i=1:n
        if i>1
            z0(i)        = z0(i-1) + COV_Z1(i-1)*H(i)*inv(PSIY(i))*(y(i) - sys.H*shh(i));
        else
            z0(i)        = zp + COV_Z1_0*H(i)*inv(PSIY(i))*(y(i)-sys.H*shh(i));
        end
    end
    e      = z - sys.H*sh(1:end-1);
    MSE(j) = (z0(end)-sh(1))^2;
end
plot(z0);
sh(1)
formula = -(1/(2*n))*log2(COV_Z1(end))
emp = log2(1/mean(MSE))/(2*n)

