%tic;
B = 200;
t = 0.02;
p = 200;
a = 8.081;
seed = 7;
n = 200;
batch = 16;
c0_min = 1;
c0_max = 20;
Beta_true = [0.4;0.8;0.4;zeros(p-3,1)];
Beta_all = zeros(B,p);
l1 = 0;
l2 = 0;
for iter = 1:B
    %generate data
    Z = cell(batch,1);
    T_tilde = cell(batch,1);
    delta = cell(batch,1);
    Z1 = [];
    T_tilde1 =[];
    delta1 = [];
    rng(seed+p+iter)
    for i = 1:batch
        [Z_0,T_tilde_0,delta_0] = g_data(n,a,Beta_true,p);
        Z{i,1} = Z_0;
        T_tilde{i,1} = T_tilde_0;
        delta{i,1} = delta_0;
        Z1 = [Z1;Z{i,1}];
        T_tilde1 =[T_tilde1;T_tilde{i,1}];
        delta1 = [delta1;delta{i,1}];
    end
    %estimation
    lambda_best_1 = BIC_0(c0_min,c0_max,Z1,T_tilde1,delta1,n*batch,p,t);
    [Beta_best_1] = LASSO_0(Z1,T_tilde1,delta1,n*batch,lambda_best_1(1),p,t);
    Beta_all(iter,:) = Beta_best_1';
    l1_temp = norm(Beta_best_1-Beta_true,1);
    l2_temp = norm(Beta_best_1-Beta_true,2);
    l1 = l1+l1_temp;
    l2 = l2+l2_temp;
    iter
end
l1_mean = l1/B;
l2_mean = l2/B;
Beta_mean = mean(Beta_all);
%toc;

save 1Beta_all.mat Beta_all -v6
save 1Beta_mean.mat Beta_mean -v6
save 1delta.mat delta -v6
save 1l1_mean.mat l1_mean -v6
save 1l2_mean.mat l2_mean -v6
save 1T_tilde.mat T_tilde -v6
save 1Z.mat Z -v6