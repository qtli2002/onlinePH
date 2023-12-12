censor_time = 18;
p = 161;
seed = 2023;
N = 213900;
n_1 = 1000;
batch = 84;
n_each = floor((N-n_1)/(batch-1));
t = 0.5;
c0_min = 1;
c0_max = 20;
sigma0 = sqrt(0);
sigma1 = sqrt(1);
ir = batch - 14;
T_all_s0 = zeros([ir,4]);
VAR_all_s0 = zeros([ir,4]);
p_all_s0 = zeros([ir,4]);
T_all_s1 = zeros([ir,4]);
VAR_all_s1 = zeros([ir,4]);
p_all_s1 = zeros([ir,4]);

% load data
load('D:\1mat\real_data_analysis\output\online_output\n1000\1Beta_all.mat');
load('D:\1mat\real_data_analysis\output\online_output\n1000\inference_result\1omega_all.mat')
load('D:\1mat\real_data_analysis\output\online_output\n1000\inference_result\1sig_index_all.mat')
sig_index = [1:p;sig_index_all];

Z_0 = load('D:\1mat\real_data_analysis\SEERdata\dataset_covariate.csv');
[max,index] = max(Z_0,[],1);
[min,index] = min(Z_0,[],1);
Z_t = (Z_0-min)./(max-min);
T = load('D:\1mat\real_data_analysis\SEERdata\T.csv');
delta_t = T <= censor_time;
T_tilde_t = delta_t .* T + (1-delta_t).* censor_time;

Z = cell(batch,1);
T_tilde = cell(batch,1);
delta = cell(batch,1);

G_1 = sig_index(1,sig_index_all(ir,:)==1);
G_2 = sig_index(1,sig_index_all(ir,:)==0);

%the first batch
N_m = n_1;
Z_1 = Z_t(1:n_1,:);
T_tilde_1 = T_tilde_t(1:n_1,:);
delta_1 = delta_t(1:n_1,:);
Z{1,1} = Z_1;
T_tilde{1,1} = T_tilde_1;
delta{1,1} = delta_1;
m = 1;
Beta_best_1 = Beta_all(1,:)';
load('D:\1mat\real_data_analysis\output\online_output\n1000\1df_1.mat')
load('D:\1mat\real_data_analysis\output\online_output\n1000\1df2_1.mat')
load('D:\1mat\real_data_analysis\output\online_output\n1000\1df_each_1.mat')
load('D:\1mat\real_data_analysis\output\online_output\n1000\1df_times_1.mat')
load('D:\1mat\real_data_analysis\output\online_output\n1000\1H_1.mat')
[temp1_1s0,temp2_1s0] = T_temp_new(delta_1,sigma0,m,n_1,p,seed);
[temp1_1s1,temp2_1s1] = T_temp_new(delta_1,sigma1,m,n_1,p,seed);
[s,Cov0,Cov_est] = Cov_threshold_new(T_tilde_t,Z_t,m,N_m,p,seed);
[temp4_1,temp5_1] = T_temp4_new(G_1,G_2,Beta_best_1,Z_1,T_tilde_1,delta_1,Cov_est,df2_1,n_1,p);

H_sum = H_1*n_1;
H_Beta = H_1*Beta_best_1*n_1;
Beta_last = Beta_best_1;
df_sum = df_1*n_1;
df_times_sum = df_times_1;
T_temp1s0_sum = temp1_1s0;
T_temp2s0_sum = temp2_1s0;
T_temp1s1_sum = temp1_1s1;
T_temp2s1_sum = temp2_1s1;
Z_times_sum = Z_1'*Z_1;
T_temp4_sum = temp4_1;
T_temp5_sum = temp5_1;

k = 0;
for m= 2:batch
    n_m = n_each;
    if(m==batch)
        n_m = n_each + rem(N-n_1,batch-1);
    end
    N_m = N_m + n_m; N_last = N_m - n_m;
    Z{m,1} = Z_t((N_last+1):N_m,:);
    T_tilde{m,1} = T_tilde_t((N_last+1):N_m,:);
    delta{m,1} = delta_t((N_last+1):N_m,:);
    Beta_best_2 = Beta_all(m,:)';
    Beta_last_real = Beta_last;
    [df_2,df2_2,df_each_2,df_times_2,H_2] = DF_new1(Beta_best_2,Z{m,1},T_tilde{m,1},delta{m,1},n_m,p);
    [temp1_2s0,temp2_2s0] = T_temp_new(delta{m,1},sigma0,m,n_m,p,seed);
    [temp1_2s1,temp2_2s1] = T_temp_new(delta{m,1},sigma1,m,n_m,p,seed);
    [s,Cov0,Cov_est] = Cov_threshold_new(T_tilde_t,Z_t,m,N_m,p,seed);
    [temp4_2,temp5_2] = T_temp4_new(G_1,G_2,Beta_best_2,Z{m,1},T_tilde{m,1},delta{m,1},Cov_est,df2_2,n_m,p);
    H_sum = H_sum+H_2*n_m;
    df_sum = df_sum+df_2*n_m;
    df_times_sum = df_times_sum+df_times_2;
    H_Beta = H_Beta+H_2*Beta_best_2*n_m;
    T_temp1s0_sum = T_temp1s0_sum+temp1_2s0;
    T_temp2s0_sum = T_temp2s0_sum+temp2_2s0;
    T_temp1s1_sum = T_temp1s1_sum+temp1_2s1;
    T_temp2s1_sum = T_temp2s1_sum+temp2_2s1;
    Z_times_sum = Z_times_sum+Z{m,1}'*Z{m,1};
    T_temp4_sum = T_temp4_sum+temp4_2;
    T_temp5_sum = T_temp5_sum+temp5_2;
    Beta_last = Beta_best_2;%After ending, Beta_last = Beta_hat_final
    if(m>=15)
        k = k+1;
        H_sum_less = H_sum-H_2*n_m;
        H_Beta_less = H_Beta-H_2*Beta_best_2*n_m;
        df_sum_less = df_sum-df_2*n_m;
        H_ave = H_sum/N_m;
        omega = omega_all{1,k};
        Beta_d = Beta_last+omega'*(H_Beta_less-H_sum_less*Beta_last-df_sum)/N_m;
        var = omega'*df_times_sum*omega/(N_m);
        st = sqrt(diag(var));
        U = U_n(Z{m,1},T_tilde{m,1},delta{m,1},n_m,p,N_m,H_sum_less,Beta_last_real,Beta_last);
        T_temp3 = omega*(H_Beta_less-H_sum_less*Beta_last_real-df_sum_less)/N_m;
        %s0
        %test1
        [T_1s0,VAR_1s0,T_2s0,VAR_2s0] = ...
                test1_new(G_1,G_2,Beta_last,omega,U,T_temp1s0_sum,T_temp2s0_sum,T_temp3,df_times_sum,N_m,sigma0);
        p1_t1_s0 = 1-cdf('Normal',T_1s0*sqrt(N_m)/sqrt(VAR_1s0),0,1);
        p2_t1_s0 = 1-cdf('Normal',T_2s0*sqrt(N_m)/sqrt(VAR_2s0),0,1);
        %test2
        [T_4s0,VAR_4s0,T_5s0,VAR_5s0]=...          
                test2_new(G_1,G_2,Beta_last,omega,U,T_temp1s0_sum,T_temp2s0_sum,T_temp3,df_times_sum,N_m,sigma0,Cov0,Cov_est,T_temp4_sum,T_temp5_sum);
    
        p1_t2_s0 = 1-cdf('Normal',T_4s0*sqrt(N_m)/sqrt(VAR_4s0),0,1);
        p2_t2_s0 = 1-cdf('Normal',T_5s0*sqrt(N_m)/sqrt(VAR_5s0),0,1);
        T_all_s0(k,:) = [T_1s0,T_2s0,T_4s0,T_5s0];
        VAR_all_s0(k,:) = [VAR_1s0,VAR_2s0,VAR_4s0,VAR_5s0];
        p_all_s0(k,:) = [p1_t1_s0,p2_t1_s0,p1_t2_s0,p2_t2_s0];
        %s1
        %test1
        [T_1s1,VAR_1s1,T_2s1,VAR_2s1] = ...
                test1_new(G_1,G_2,Beta_last,omega,U,T_temp1s1_sum,T_temp2s1_sum,T_temp3,df_times_sum,N_m,sigma1);
        p1_t1_s1 = 1-cdf('Normal',T_1s1*sqrt(N_m)/sqrt(VAR_1s1),0,1);
        p2_t1_s1 = 1-cdf('Normal',T_2s1*sqrt(N_m)/sqrt(VAR_2s1),0,1);
        %test2
        [T_4s1,VAR_4s1,T_5s1,VAR_5s1]=...          
                test2_new(G_1,G_2,Beta_last,omega,U,T_temp1s1_sum,T_temp2s1_sum,T_temp3,df_times_sum,N_m,sigma1,Cov0,Cov_est,T_temp4_sum,T_temp5_sum);
    
        p1_t2_s1 = 1-cdf('Normal',T_4s1*sqrt(N_m)/sqrt(VAR_4s1),0,1);
        p2_t2_s1 = 1-cdf('Normal',T_5s1*sqrt(N_m)/sqrt(VAR_5s1),0,1);
        T_all_s1(k,:) = [T_1s1,T_2s1,T_4s1,T_5s1];
        VAR_all_s1(k,:) = [VAR_1s1,VAR_2s1,VAR_4s1,VAR_5s1];
        p_all_s1(k,:) = [p1_t1_s1,p2_t1_s1,p1_t2_s1,p2_t2_s1];
    end
    m
end

save 1T_all_s0.mat T_all_s0 -v6
save 1VAR_all_s0.mat VAR_all_s0 -v6
save 1p_all_s0.mat p_all_s0 -v6
save 1temp1_1s0.mat temp1_1s0 -v6
save 1temp2_1s0.mat temp2_1s0 -v6
save 1T_all_s1.mat T_all_s1 -v6
save 1VAR_all_s1.mat VAR_all_s1 -v6
save 1p_all_s1.mat p_all_s1 -v6
save 1temp1_1s1.mat temp1_1s1 -v6
save 1temp2_1s1.mat temp2_1s1 -v6
save 1temp4_1.mat temp4_1 -v6
save 1temp5_1.mat temp5_1 -v6

