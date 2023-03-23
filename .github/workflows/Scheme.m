clear all
close all
clc

%State-space
beta = 0.5;
sys.F = beta;
sys.H = beta;
sys.G = 1;
sys.V= 1;
sys.W = 1;
sys.L = 1;

%Power
sys.P     = 1;

%%Constants for the KF-ENC
[sys.SIG,eig_P,K_Temp] = dare(sys.F',sys.H',sys.G*sys.W*sys.G',sys.V,sys.G*sys.L);
sys.Psi     = sys.H*sys.SIG*sys.H' + sys.V;
sys.Kp      = (sys.F*sys.SIG*sys.H' + sys.G*sys.L)*inv(sys.Psi);

%This function computes the capacity
sys.CVX = CVX_cap(sys);

%Two examples can be loaded at this point to skip CVX
% load('Ex_AR_05.mat');

n = 60; %blocklength

%%Constants for KF-DEC
KF_DEC.HAT_SIG(1)  = sys.Kp*sys.Psi*sys.Kp';
for i=1:n
    KF_DEC.PSIY(i)      = (sys.CVX.Gamma*inv(sys.CVX.HAT_Sig) + sys.H)*KF_DEC.HAT_SIG(i)*(sys.CVX.Gamma*inv(sys.CVX.HAT_Sig) + sys.H)' + sys.Psi;
    KF_DEC.Ky(i)        = (sys.F*KF_DEC.HAT_SIG(i)*(sys.CVX.Gamma*inv(sys.CVX.HAT_Sig) + sys.H) + sys.Kp*sys.Psi) * inv(KF_DEC.PSIY(i));
    KF_DEC.HAT_SIG(i+1) = sys.F*KF_DEC.HAT_SIG(i)*sys.F' + sys.Kp*sys.Psi*sys.Kp' - KF_DEC.Ky(i)*KF_DEC.PSIY(i)*KF_DEC.Ky(i)';
end


%%Constants for smoothing
COV_Z_0   = sys.Psi;
cl_loop(1) = 1; % We use it compute the power of the closed loop: F-Kp*(Lambda*Gamma*inv(HAT_SIG)+H))
for i=1:n
    kappa(i)            = (sys.CVX.Gamma*inv(sys.CVX.HAT_Sig) + sys.H)*cl_loop(i)*sys.Kp;
    cl_loop(i+1)        = (sys.F - sys.Kp*(sys.CVX.Gamma*inv(sys.CVX.HAT_Sig) + sys.H))*cl_loop(i);
    if i>1
        COV_Z0(i)       = Var_sm(COV_Z0(i-1),kappa(i),KF_DEC.PSIY(i));
    else
        COV_Z0(i)       = Var_sm(COV_Z_0,kappa(i),KF_DEC.PSIY(i));
    end
end

%%Random experiment
EXP = 5000;
for j=1:EXP
    %Noise generation
    [s,z] = noise(sys,n);
    
    %Message transmission
    z_0         = sqrt(sys.Psi)*randn(1); %This corresponds to the z_0 in the paper.
    
    %Initialization
    sh(1)       = sys.Kp*z_0;
    shh(1)      = 0;
    
    %Refinement stage
    for i=1:n
        x(i)         = (sys.CVX.Gamma/sys.CVX.HAT_Sig)*(sh(i)-shh(i));
        y(i)         = x(i) + z(i);
        sh(i+1)      = KF_enc(sh(i),z(i),sys);
        e(i)         = z(i) - sys.H*sh(i);
        shh(i+1)     = KF_dec(shh(i),KF_DEC.Ky(i),y(i),sys);
    end
    
    %%Smoothing based on y1,...,yn
    z_hat = 0;
    for i=1:n
        if i>1
            z0(i)        = z0(i-1) + COV_Z0(i-1)*kappa(i)*inv(KF_DEC.PSIY(i))*(y(i) - sys.H*shh(i));
        else
            z0(i)        = z_hat + COV_Z_0*kappa(i)*inv(KF_DEC.PSIY(i))*(y(i)-sys.H*shh(i));
        end
    end
    %Error at the j'th trial
    MSE_z0(j)   = (z0(end)-z_0)^2;
end

plot(z0);
hold on
plot([1:n],linspace(z_0,z_0,n));

%%Comparison of analytical and numerical rates
%Rate based on MSE formula in Lemma 1 (Analytcal)
MSE_formula = -(1/(2*(n+1)))*log2(COV_Z0(end))
%Rate based on numerical reaults
Empirical_rate = log2(1/mean(MSE_z0))/(2*(n+1))
%The capacity
Capacity = sys.CVX.Capacity
