function[Z,T_tilde,delta] = g_data(n,p,beta_0,a)
mu = zeros(1,p);
sigma = ones(p)*0.2;
temp_s1 = repmat(1:p,p,1);
temp_s2 = repmat((1:p)',1,p);
sigma = sigma.^abs(temp_s1-temp_s2);

Z = mvnrnd(mu,sigma,n);%'
C = unifrnd(0,a,1,n);
S_T = rand(1,n);
T = -log(S_T)./exp(beta_0'*Z'); 
delta = (T <= C)';
T_tilde = delta .* T' + (1-delta).*C';
end