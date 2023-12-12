function [beta]=LASSO_0(Z,T_tilde,delta,n,lambda,p,t,beta_temp)
old=digits(10);
Y_0 = (T_tilde >= T_tilde');
df1 = mean(Z.*repmat(delta,1,p))';
df2 = ones([p 1]);
X = ones([p 1]);
%first iteration before while
temp1 = Y_0.*repmat(exp(Z*beta_temp),1,n);
temp2 = mean(Y_0.*repmat(exp(Z*beta_temp),1,n));
for i =1:p
    df2(i) = mean(mean(repmat(Z(:,i),1,n).*temp1)./temp2.*delta');
end
df = df2 - df1;
X = beta_temp - t*df;
beta = zeros(p,1)...
    + (X>repmat(lambda*t,p,1)).*(X-repmat(lambda*t,p,1))...
    + (X<(-repmat(lambda*t,p,1))).*(X+repmat(lambda*t,p,1));
iter=0;
while vpa(norm(beta-beta_temp,inf)) > 10^(-3)
    iter=iter+1;
    %display([iter,n,lambda*n/log(p)])
    beta_temp = beta;
    temp1 = Y_0.*repmat(exp(Z*beta),1,n);
    temp2 = mean(Y_0.*repmat(exp(Z*beta),1,n));
    for i =1:p
        df2(i) = mean(mean(repmat(Z(:,i),1,n).*temp1)./temp2.*delta');
    end
    df = df2 - df1;
    X = beta_temp - t*df;
    beta = zeros(p,1)...
    + (X>repmat(lambda*t,p,1)).*(X-repmat(lambda*t,p,1))...
    + (X<(-repmat(lambda*t,p,1))).*(X+repmat(lambda*t,p,1));
end
end

