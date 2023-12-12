function [lambda_best]=BIC_0(c0_min,c0_max,Z,T_tilde,delta,n,p,t)
%t = 0.3;
c = 1;
C_n = c*log(log(p));
bic =ones(1,100)*10000;
iter=0;
for i=c0_min:3:c0_max
%for i=exp(c0_min:0.8:c0_max)
    iter = iter+1;
    [beta_lambda]=LASSO_0(Z,T_tilde,delta,n,i*log(p)/n,p,t);
    Y_0 = T_tilde >= T_tilde';
    C1=(beta_lambda'*Z'*delta)/n;
    C2=(delta'*log(Y_0'*exp(Z*beta_lambda)))/n;
    l_0 = p-sum(beta_lambda==0);
    %1/n log likelihood function: C1-C2
    bic_temp = C2-C1+C_n*(log(n)/n)*l_0;
    bic(iter)=bic_temp;
end
index = find(bic==min(bic));
lambda_best = (c0_min+3*(index(1)-1))*log(p)/n;
%lambda_best = exp(c0_min+0.8*(index(1)-1))*log(p)/n;
end


