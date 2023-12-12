function [index,beta_all,lambda_best]=BIC_0(c0_min,c0_max,Z,T_tilde,delta,n,p,t)
batch = 93;
beta_all = zeros([6 p]);
Y_0 = T_tilde >= T_tilde';
c = 1;
C_n = c*log(log(p));
bic = ones(1,100)*10000;
iter0 = 0;
for i=((c0_min+3):3:c0_max)*batch
    iter0 = iter0+1;
    display(iter0)
    [beta_lambda]=LASSO_0_SGD(Z,T_tilde,delta,n,i*log(p)/n,p,t);
    beta_all(iter0,:) = beta_lambda';
    %Compute BIC where BIC(lambda) = L(beta_lambda) + C_n*log(n)/n*l_0(beta_lambda)
    C1=(beta_lambda'*Z'*delta)/n;
    C2_sum = 0;
    for j = 1:batch
        Z_temp = Z(1+(j-1)*2300:j*2300,:);
        delta_temp = delta(1+(j-1)*2300:j*2300,:);
        T_tilde_temp = T_tilde(1+(j-1)*2300:j*2300,:);
        Y_0_temp = T_tilde_temp >= T_tilde_temp';
        C2_sum = C2_sum + ...
            delta_temp'*log(Y_0_temp'*exp(Z_temp*beta_lambda));
    end
    C2 = C2_sum/n;
    l_0 = p-sum(beta_lambda==0);
    bic_temp = C2-C1+C_n*(log(n)/n)*l_0;
    bic(iter0)=bic_temp;
end
index = find(bic==min(bic));
lambda_best = (c0_min+3+3*(index(1)-1))*batch*log(p)/n;
end


