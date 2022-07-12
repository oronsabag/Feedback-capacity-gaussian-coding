function shh = shat_dec(shh_p,Ky,y,sys)

    shh = sys.F*shh_p + Ky*(y - sys.H*shh_p);
