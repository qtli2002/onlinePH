function [temp1,temp2]=T_temp_new(delta,sigma,m,n,p,seed)
mu = 0;
rng(p+seed+m)
S = normrnd(mu,sigma,[n 1]);
temp1 = delta'*S;
temp2 = sum(delta);
end