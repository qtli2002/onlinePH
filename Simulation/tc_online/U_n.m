function [df]=U_n(Z,T_tilde,delta,n,p,m,H_sum,Beta_last_real,Beta_best)
old=digits(10);
%p = 200;
Y_0 = (T_tilde >= T_tilde');
df1 = mean(Z.*repmat(delta,1,p))';
df2 = ones([p 1]);
temp1 = Y_0.*repmat(exp(Z*Beta_best),1,n);
temp2 = mean(Y_0.*repmat(exp(Z*Beta_best),1,n));
for i =1:p
    df2(i) = mean(mean(repmat(Z(:,i),1,n).*temp1)./temp2.*delta');
end
l_df = df2 - df1;
df = (H_sum*(Beta_best-Beta_last_real)+l_df)/m;
end

