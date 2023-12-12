censor_time = 18;
p = 161;
n = 2300;
batch = 93;
t = 0.5;
c0_min = 1;
c0_max = 20;
Z_0 = load('D:\1mat\real_data_analysis\SEERdata\dataset_covariate.csv');
[max, index] = max(Z_0,[],1);
[min, index] = min(Z_0,[],1);
Z = (Z_0-min)./(max-min);
T = load('D:\1mat\real_data_analysis\SEERdata\T.csv');
delta = T <= censor_time;
T_tilde = delta .* T + (1-delta).* censor_time;

[best_index,beta_all,lambda_best_1] = BIC_0(c0_min,c0_max,Z,T_tilde,delta,n*batch,p,t);
Beta_best_1 = beta_all(best_index(1),:)';

save 1best_index.mat best_index -v6
save 1beta_all.mat beta_all -v6
save 1Beta_best_1.mat Beta_best_1 -v6

