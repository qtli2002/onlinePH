censor_time = 18;
p = 161;
N = 213900;
n_1 = 1000;
batch = 84;
n_each = floor((N-n_1)/(batch-1));
t = 0.5;
c0_min = 1;
c0_max = 20;
t_alpha = icdf('norm',0.975,0,1);
z_alpha = icdf('norm',0.95,0,1);
ir = batch - 14;
Beta_d_all = zeros([ir,p]);
omega_all = cell(1,ir);
var_all = cell(1,ir);
st_all = zeros([ir,p]);
sig_index_all = zeros([ir p]);
for i=1:ir
    omega_all{1,i} = zeros([p p]);
    var_all{1,i} = zeros([p p]);
end
% load data
load('D:\1mat\real_data_analysis\output\online_output\n1000\1Beta_all.mat');
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

%the first batch
N_m = n_1;
Z_1 = Z_t(1:n_1,:);
T_tilde_1 = T_tilde_t(1:n_1,:);
delta_1 = delta_t(1:n_1,:);
Z{1,1} = Z_1;
T_tilde{1,1} = T_tilde_1;
delta{1,1} = delta_1;

Beta_best_1 = Beta_all(1,:)';
load('D:\1mat\real_data_analysis\output\online_output\n1000\1df_1.mat')
load('D:\1mat\real_data_analysis\output\online_output\n1000\1df_each_1.mat')
load('D:\1mat\real_data_analysis\output\online_output\n1000\1df_times_1.mat')
load('D:\1mat\real_data_analysis\output\online_output\n1000\1H_1.mat')
H_sum = H_1*n_1;
H_Beta = H_1*Beta_best_1*n_1;
Beta_last = Beta_best_1;
df_sum = df_1*n_1;
df_times_sum = df_times_1;

k = 0;
for m=2:batch
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
    H_sum = H_sum+H_2*n_m;
    df_sum = df_sum+df_2*n_m;
    df_times_sum = df_times_sum+df_times_2;
    H_Beta = H_Beta+H_2*Beta_best_2*n_m;
    Beta_last = Beta_best_2;
    if(m>=15)
        k = k+1;
        H_sum_less = H_sum-H_2*n_m;
        H_Beta_less = H_Beta-H_2*Beta_best_2*n_m;
        df_sum_less = df_sum-df_2*n_m;
        H_ave = H_sum/N_m;
        clear omega
        delete("data_H_sum\omega.mat")
        save H_ave.mat H_ave -v6
        Rpath = 'D:\Program Files\R\R-4.2.1\bin';
        filepath = 'D:\1mat\real_data_analysis\code\real_online_inference\data_H_sum\omega_compute.R';
        for ei = 1:100
            ei
            try
                RunRcode(filepath,Rpath);
                load('D:\1mat\real_data_analysis\code\real_online_inference\data_H_sum\omega.mat')
            end
            if exist('omega','var')
                break
            end
        end
        Beta_d = Beta_last+omega'*(H_Beta_less-H_sum_less*Beta_last-df_sum)/N_m;
        Beta_d_all(k,:) = Beta_d';
        omega_all{1,k} = omega;
        var = omega'*df_times_sum*omega/N_m;
        var_all{1,k} = var;
        st = sqrt(diag(var));
        st_all(k,:) = st';
        sig_index = 1-((((Beta_d-t_alpha/sqrt(N_m)*st)<0)...
            +(0<(Beta_d+t_alpha/sqrt(N_m)*st)))==2);
        sig_index_all(k,:) = sig_index';
    end
    m  
end


save 1Beta_d_all.mat Beta_d_all -v6
save 1omega_all.mat omega_all -v6
save 1var_all.mat var_all -v6
save 1st_all.mat st_all -v6
save 1sig_index_all.mat sig_index_all -v6
