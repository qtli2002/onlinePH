function [beta]=LASSO_0_SGD(Z0,T_tilde0,delta0,n,lambda,p,t)
old=digits(10);
k = 7000;
seed = 0;
rnd_index = zeros([1 k]);
X = ones([p 1]);
%initialize beta
beta_temp = zeros([p 1]);
%first iteration before while
rng(seed)
rnd_index = randperm(n,k);
Z = Z0(rnd_index,:); delta = delta0(rnd_index,:); T_tilde = T_tilde0(rnd_index,:);
Y_0 = (T_tilde >= T_tilde');
df1 = mean(Z.*repmat(delta,1,p))';
df2 = ones([p 1]);
temp1 = Y_0.*repmat(exp(Z*beta_temp),1,k);
temp2 = mean(Y_0.*repmat(exp(Z*beta_temp),1,k));
for i =1:p
    df2(i) = mean(mean(repmat(Z(:,i),1,k).*temp1)./temp2.*delta');
end
df = df2 - df1;
X = beta_temp - t*df;
beta = zeros(p,1)...
    + (X>repmat(lambda*t,p,1)).*(X-repmat(lambda*t,p,1))...
    + (X<(-repmat(lambda*t,p,1))).*(X+repmat(lambda*t,p,1));
%iteration using ISTA
iter=0;
while vpa(norm(beta-beta_temp,inf)) > 10^(-3)
    seed = seed+1;
    iter=iter+1;
    display(iter);
    rng(seed)
    rnd_index = randperm(n,k);
    Z = Z0(rnd_index,:); delta = delta0(rnd_index,:); T_tilde = T_tilde0(rnd_index,:);
    Y_0 = (T_tilde >= T_tilde');    
    beta_temp = beta;  
    df1 = mean(Z.*repmat(delta,1,p))';
    temp1 = Y_0.*repmat(exp(Z*beta),1,k);
    temp2 = mean(Y_0.*repmat(exp(Z*beta),1,k));
    for i =1:p
        df2(i) = mean(mean(repmat(Z(:,i),1,k).*temp1)./temp2.*delta');
    end
    df = df2 - df1;
    X = beta_temp - t*df;
    beta = zeros(p,1)...
    + (X>repmat(lambda*t,p,1)).*(X-repmat(lambda*t,p,1))...
    + (X<(-repmat(lambda*t,p,1))).*(X+repmat(lambda*t,p,1));
end
end

