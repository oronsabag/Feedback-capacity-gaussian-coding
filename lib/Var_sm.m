function V_Z = Var_sm(V_p,H,Psi)

V_Z = V_p - V_p*H'*inv(Psi)*H*V_p;
