function CVX = CVX_cap(sys)

d=size(sys.F,1);

cvx_begin sdp
variable Pi(1,1) symmetric
variable Gamma(1,d)
variable Sig(d,d) symmetric
Pi>=0;
Sig>=0;
[Pi Gamma; Gamma' Sig] > 0;
trace(Pi)<=sys.P;
[sys.F*Sig*sys.F' + sys.Kp*sys.Psi*sys.Kp' - Sig , sys.F*Gamma' + sys.F*Sig*sys.H' + sys.Kp*sys.Psi;
(sys.F*Gamma' + sys.F*Sig*sys.H' + sys.Kp*sys.Psi')' ,  Pi + sys.H*Sig*sys.H' + Gamma*sys.H' + sys.H*Gamma'+ sys.Psi] >= 0;
maximize log_det(Pi + sys.H*Sig*sys.H' + Gamma*sys.H' + sys.H*Gamma' + sys.Psi) - log_det(sys.Psi);
cvx_end

CVX.Capacity = 0.5*cvx_optval/log(2);
CVX.Gamma = Gamma;
CVX.HAT_Sig = Sig;
CVX.Pi = Pi;
CVX.PsiY = Pi + sys.H*CVX.HAT_Sig*sys.H' + Gamma*sys.H' + sys.H*Gamma' + sys.Psi;
CVX.KY = (sys.F*CVX.Gamma + sys.F*CVX.HAT_Sig*sys.H + sys.Kp*sys.Psi)*inv(CVX.PsiY);

M = num2str(Pi-Gamma*inv(Sig)*Gamma);
disp(sprintf('M is equal to %s',M));
