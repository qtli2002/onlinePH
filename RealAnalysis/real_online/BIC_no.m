function [Beta_best,lambda_best]=BIC_no(c0_min,c0_max,Z,T_tilde,delta,n,N_m,p,t,H_sum,Beta_last,beta_temp)
c = 1;
C_n = c*log(log(p));
bic =ones(1,100)*10000;
iter=0;
for i=(c0_min):3:c0_max
    iter=iter+1;
    [beta_lambda]=LASSO(Z,T_tilde,delta,n,N_m,i*log(p)/N_m,p,t,H_sum,Beta_last,beta_temp);
    Beta_all(iter,:) = beta_lambda';
    Y_0 = T_tilde >= T_tilde';
    C1=(beta_lambda'*Z'*delta)/n;
    C2=(delta'*log(Y_0'*exp(Z*beta_lambda)))/n;
    l_0 = p-sum(beta_lambda==0);
    %1/n log likelihood function: C1-C2
    l_m = C2-C1;
    bic_temp =(((beta_lambda-Beta_last)'*H_sum*(beta_lambda-Beta_last))/2+l_m*n)/N_m...
        +C_n*(log(N_m)/(N_m))*l_0;
    bic(iter)=bic_temp;
end
index = find(bic==min(bic));
lambda_best = (c0_min+3*(index(1)-1))*log(p)/(N_m);
Beta_best = Beta_all(index(1),:)';
end


