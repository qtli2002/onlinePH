censor_time = 18;
p = 161;
N = 213900;
n_1 = 20000;
batch = 84;
n_each = floor((N-n_1)/(batch-1));
t = 0.5;
c0_min = 1;
c0_max = 20;
t_alpha = icdf('norm',0.975,0,1);
z_alpha = icdf('norm',0.95,0,1);
Beta_all = zeros([batch,p]);
beta_temp_all = zeros([batch,p]);
H_sum_all = cell(1,batch);

%load data
load('D:\1mat\real_data_analysis\output\offline_output\1Beta_best_1.mat')
beta_temp = Beta_best_1;
Z_0 = load('D:\1mat\real_data_analysis\SEERdata\dataset_covariate.csv');
[max, index] = max(Z_0,[],1);
[min, index] = min(Z_0,[],1);
Z = (Z_0-min)./(max-min);
T = load('D:\1mat\real_data_analysis\SEERdata\T.csv');
delta = T <= censor_time;
T_tilde = delta .* T + (1-delta).* censor_time;

%the first batch
N_m = n_1;
Z_1 = Z(1:n_1,:);
T_tilde_1 = T_tilde(1:n_1,:);
delta_1 = delta(1:n_1,:);
beta_temp_all(1,:) = beta_temp';
[Beta_best_1,lambda_best_1] = BIC_0(c0_min,c0_max,Z_1,T_tilde_1,delta_1,n_1,p,t,beta_temp);
[df_1,df2_1,df_each_1,df_times_1,H_1] = DF_new1(Beta_best_1,Z_1,T_tilde_1,delta_1,n_1,p);

%other batches
H_sum = H_1*n_1;
Beta_last = Beta_best_1;
Beta_all(1,:) = Beta_last';
H_sum_all{1,1} = H_sum;
k = 0;
for m=2:batch
    n_m = n_each;
    if(m==batch)
        n_m = n_each + rem(N-n_1,batch-1);
    end
    N_m = N_m + n_m; N_last = N_m - n_m;
    Z_m = Z((N_last+1):N_m,:);
    T_tilde_m = T_tilde((N_last+1):N_m,:);
    delta_m = delta((N_last+1):N_m,:);
    Beta_last_real = Beta_last;
    beta_temp_all(m,:) = beta_temp';
    [Beta_best_2,lambda_best_2] = BIC_no(c0_min,c0_max,Z_m,T_tilde_m,delta_m,n_m,N_m,p,t,H_sum,Beta_last,beta_temp);
    [df_2,df2_2,df_each_2,df_times_2,H_2] = DF_new1(Beta_best_2,Z_m,T_tilde_m,delta_m,n_m,p);
    H_sum = H_sum+H_2*n_m;
    Beta_last = Beta_best_2;%After ending, Beta_last = Beta_hat_final
    Beta_all(m,:) = Beta_last';
    H_sum_all{1,m} = H_sum;
    m
end

save 1Beta_all.mat Beta_all -v6
save 1beta_temp_all.mat beta_temp_all -v6
save 1H_sum_all.mat H_sum_all -v6
save 1df_1.mat df_1 -v6
save 1df2_1.mat df2_1 -v6
save 1df_each_1.mat df_each_1 -v6
save 1df_times_1.mat df_times_1 -v6
save 1H_1.mat H_1 -v6
