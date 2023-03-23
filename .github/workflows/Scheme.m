clear all
close all
clc

%State-space parameters
beta = 1.3;
sys.F = beta;
sys.H = beta;
sys.G = 1;
sys.V= 1;
sys.W = 1;
sys.L = 1;

%Power parameter
sys.P     = 1;

%%Kalman filter for the encoder
[sys.SIG,eig_P,K_Temp] = dare(sys.F',sys.H',sys.G*sys.W*sys.G',sys.V,sys.G*sys.L);
sys.Psi     = sys.H*sys.SIG*sys.H' + sys.V;
sys.Kp      = (sys.F*sys.SIG*sys.H' + sys.G*sys.L)*inv(sys.Psi);

%%This function computes the capacity
sys.CVX = CVX_cap(sys);
% sys.dec.PsiY = sys.CVX.Pi + sys.H*sys.CVX.Sig*sys.H'+ 2*sys.CVX.Gamma*sys.H'+ sys.Psi;
% sys.dec.KY = (sys.F*sys.CVX.Gamma + sys.F*sys.CVX.Sig*sys.H + 1)*inv(sys.dec.PsiY);

% load('sys.mat');
n = 40; %blocklength

KF_DEC.HAT_SIG(1)  = sys.Kp*sys.Psi*sys.Kp';
for i=1:n
    KF_DEC.PSIY(i)      = (sys.CVX.Gamma/sys.CVX.HAT_Sig + sys.H)*KF_DEC.HAT_SIG(i)*(sys.CVX.Gamma/sys.CVX.HAT_Sig + sys.H)' + 1;
    KF_DEC.Ky(i)        = (sys.F*KF_DEC.HAT_SIG(i)*(sys.CVX.Gamma/sys.CVX.HAT_Sig + sys.H) + 1) * inv(KF_DEC.PSIY(i));
    KF_DEC.HAT_SIG(i+1) = sys.F*KF_DEC.HAT_SIG(i)*sys.F' + 1 - KF_DEC.Ky(i)*KF_DEC.PSIY(i)*KF_DEC.Ky(i)';
end

COV_Z1_0   = sys.Psi;
cl_loop(1) = 1; % We use it compute the power of the closed loop: F-Kp*(Lambda*Gamma*inv(HAT_SIG)+H))
for i=1:n
       cl_loop(i+1)   = (sys.F - sys.Kp*(sys.CVX.Gamma/sys.CVX.HAT_Sig + sys.H))*cl_loop(i);
    kappa(i)         = (sys.CVX.Gamma/sys.CVX.HAT_Sig + sys.H)*cl_loop(i)*sys.Kp;
    if i>1
        COV_Z1(i)  = Var_sm(COV_Z1(i-1),kappa(i),KF_DEC.PSIY(i));
    else
        COV_Z1(i)  = Var_sm(COV_Z1_0,kappa(i),KF_DEC.PSIY(i));
    end
 
end

%%Random experiment
EXP = 10000;
for j=1:EXP
    %Noise and state generation
    [s,z] = noise(sys,n);
    
    %Initialization step
    sh(1)       = sys.Kp*z(1);
    sh(1)       = sys.Kp*randn(1);
    shh(1)      = 0;
%     shh2(1)      = 0;
    s(1) = randn(1);
%     sh(1) = s(1);
    for i=1:n
        x(i)         = (sys.CVX.Gamma/sys.CVX.HAT_Sig)*(sh(i)-shh(i));
        y(i)         = x(i) + z(i);
        sh(i+1)      = KF_enc(sh(i),z(i),sys);
        e(i) = z(i) - sys.H*sh(i);
%         ytilde(i) = (sys.CVX.Gamma/sys.CVX.HAT_Sig+sys.H)*(sh(i)-shh(i)) + e(i);
        shh(i+1)     = KF_dec(shh(i),KF_DEC.Ky(i),y(i),sys);
%         shh2(i+1)    = shat_dec(shh2(i),sys.dec.KY,y(i),sys);
    end
    
    %%Smoothing formula based on y1,...,yn
    z_hat = 0;
    for i=1:n
        if i>1
            z0(i)        = z0(i-1) + COV_Z1(i-1)*kappa(i)*inv(KF_DEC.PSIY(i))*(y(i) - sys.H*shh(i));
        else
            z0(i)        = z_hat + COV_Z1_0*kappa(i)*inv(KF_DEC.PSIY(i))*(y(i)-sys.H*shh(i));
        end
    end
%     e      = z - sys.H*sh(1:end-1);
    MSE(j) = (z0(end)-sh(1))^2;
end
plot(z0);
sh(1)
formula = -(1/(2*n))*log2(COV_Z1(end))
emp = log2(1/mean(MSE))/(2*n)

