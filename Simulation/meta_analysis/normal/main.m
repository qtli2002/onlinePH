B = 200;
t = 0.02;
p = 200;
a = 7.248;
seed = 2;
n = 200;
batch = 16;

c0_min = 1;
c0_max = 20;
Beta_true = [0.4;0.8;0.4;zeros(p-3,1)];

Beta_offline_store = cell(B,1);
Beta_final_store = cell(B,1);

bias_offline_mean = zeros(batch,p);
bias_final_mean = zeros(batch,p);
l1_offline_mean = zeros(batch,1);
l1_final_mean = zeros(batch,1);
l2_offline_mean = zeros(batch,1);
l2_final_mean = zeros(batch,1);
for iter = 1:B
    Beta_offline_all = zeros(batch,p);
    Beta_final_all = zeros(batch,p);
    Z = cell(batch,1);
    T_tilde = cell(batch,1);
    delta = cell(batch,1);
    Z1 = [];
    T_tilde1 =[];
    delta1 = [];
    rng(seed+p+iter)
    for i = 1:batch
        [Z_0,T_tilde_0,delta_0] = g_data(n,p,Beta_true,a);
        Z{i,1} = Z_0;
        T_tilde{i,1} = T_tilde_0;
        delta{i,1} = delta_0;
        Z1 = [Z1;Z{i,1}];
        T_tilde1 =[T_tilde1;T_tilde{i,1}];
        delta1 = [delta1;delta{i,1}];
    end
    %first batch
    m = 1;
    N_m = length(T_tilde{1,1});
    lambda_best_1 = BIC_0(c0_min,c0_max,Z{1,1},T_tilde{1,1},delta{1,1},n,p,t);
    [Beta_offline_1] = LASSO_0(Z{1,1},T_tilde{1,1},delta{1,1},n,lambda_best_1(1),p,t);
    Beta_final_1 = Beta_offline_1;
    Beta_final_last = Beta_final_1;
    Beta_offline_all(1,:) = Beta_offline_1';
    Beta_final_all(1,:) = Beta_final_1';
    
    bias_offline_mean(m,:) = bias_offline_mean(m,:)+abs(Beta_offline_1 - Beta_true)';
    bias_final_mean(m,:) = bias_final_mean(m,:)+abs(Beta_final_1 - Beta_true)';
    l1_offline_mean(m,1) = l1_offline_mean(m,1) + norm(Beta_offline_1 - Beta_true,1);
    l1_final_mean(m,1) = l1_final_mean(m,1) + norm(Beta_final_1 - Beta_true,1);
    l2_offline_mean(m,1) = l2_offline_mean(m,1) + norm(Beta_offline_1 - Beta_true,2);
    l2_final_mean(m,1) = l2_final_mean(m,1) + norm(Beta_final_1 - Beta_true,2);
    for m=2:batch
        n_m = length(T_tilde{m,1});
        N_m = N_m + n_m;
        N_last = N_m - n_m;
        lambda_best_m = BIC_0(c0_min,c0_max,Z{m,1},T_tilde{m,1},delta{m,1},n,p,t);
        [Beta_offline_m] = LASSO_0(Z{m,1},T_tilde{m,1},delta{m,1},n,lambda_best_m(1),p,t);
        Beta_final_m = N_last/N_m*Beta_final_last + n_m/N_m*Beta_offline_m;
        Beta__final_last = Beta_final_m;
        Beta_offline_all(m,:) = Beta_offline_m';
        Beta_final_all(m,:) = Beta_final_m';

        bias_offline_mean(m,:) = bias_offline_mean(m,:)+abs(Beta_offline_m - Beta_true)';
        bias_final_mean(m,:) = bias_final_mean(m,:)+abs(Beta_final_m - Beta_true)';
        l1_offline_mean(m,1) = l1_offline_mean(m,1) + norm(Beta_offline_m - Beta_true,1);
        l1_final_mean(m,1) = l1_final_mean(m,1) + norm(Beta_final_m - Beta_true,1);
        l2_offline_mean(m,1) = l2_offline_mean(m,1) + norm(Beta_offline_m - Beta_true,2);
        l2_final_mean(m,1) = l2_final_mean(m,1) + norm(Beta_final_m - Beta_true,2);
        [iter,m]
    end
    Beta_offline_store{iter,1} = Beta_offline_all;
    Beta_final_store{iter,1} = Beta_final_all;
end

bias_offline_mean = bias_offline_mean/B;
bias_final_mean = bias_final_mean/B;
l1_offline_mean = l1_offline_mean/B;
l1_final_mean = l1_final_mean/B;
l2_offline_mean = l2_offline_mean/B;
l2_final_mean = l2_final_mean/B;

save 1Beta_offline_store.mat Beta_offline_store
save 1Beta_final_store.mat Beta_final_store
save 1bias_offline_mean.mat bias_offline_mean
save 1bias_final_mean.mat bias_final_mean
save 1l1_offline_mean.mat l1_offline_mean
save 1l1_final_mean.mat l1_final_mean
save 1l2_offline_mean.mat l2_offline_mean
save 1l2_final_mean.mat l2_final_mean
save 1Z.mat Z
save 1T_tilde.mat T_tilde
save 1delta.mat delta