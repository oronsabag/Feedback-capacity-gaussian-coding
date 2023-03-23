function shat = KF_enc(shat_p,z,sys)

   shat = sys.F*shat_p + sys.Kp*(z - sys.H*shat_p);