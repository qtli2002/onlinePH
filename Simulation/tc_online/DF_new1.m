function [df,df2,df_each,df_times,H]=DF_new1(beta,Z,T_tilde,delta,n,p)
Y_0 = (T_tilde >= T_tilde');
temp1 = Y_0.*repmat(exp(Z*beta),1,n);
temp2 = mean(Y_0.*repmat(exp(Z*beta),1,n));
%df
df1 = mean(Z.*repmat(delta,1,p))';
df2 = ones([p 1]);
for i =1:p
    df2(i) = mean(mean(repmat(Z(:,i),1,n).*temp1)./temp2.*delta');
end
df = df2 - df1;
df_times1 = Z'.*repmat(delta',p,1);
df_times2 = zeros(p,n);
for i =1:p
    df_times2(i,:) = mean(repmat(Z(:,i),1,n).*temp1)./temp2.*delta';
end
df_times = (df_times2-df_times1)*((df_times2-df_times1)');
df_each = df_times2-df_times1;
%Hessian
H = zeros(p,p);
for i =1:n
    S_2 = (Z'.*repmat(Y_0(:,i)',p,1).*repmat(exp((Z*beta))',p,1))*Z;
    S_0 = Y_0(:,i)'*exp(Z*beta);
    S_1 = (Z'.*repmat(Y_0(:,i)',p,1))*exp(Z*beta);
    H_each = delta(i,1)*(S_2/S_0-S_1*(S_1')/(S_0*S_0));
    H = H+H_each;
end
H = H/n;
end