B = 200;
t = 0.02;
p = 400;
a = 16.649;
seed = 40;
n = 200;
batch = 16;
c0_min = 1;
c0_max = 20;
t_alpha = icdf('norm',0.975,0,1);
z_alpha = icdf('norm',0.95,0,1);
sigma0 = sqrt(0);
sigma1 = sqrt(1);
Beta_true = [0.4;0.8;0.4;zeros(p-3,1)];
Beta_all = cell(1,4);
omega_all = cell(B,4);
Beta_d_all = cell(1,4);
Beta_mean = zeros(4,p);
Beta_d_mean = zeros(4,p);
Beta_store = cell(1,batch);
df_store = cell(1,batch);
H_store = cell(B,batch);
df_times_store = cell(B,batch);
Cov_est_all = cell(B,batch);
G12s0_store = cell(B,1);
G12s1_store = cell(B,1);
G3_store = cell(B,1);
G4_store = cell(B,batch);
l1 = zeros(1,4);
l2 = zeros(1,4);
l1_d = zeros(1,4);
l2_d = zeros(1,4);
st_all = cell(1,4);
cp_sum = zeros(p,4);
reject_H0_1s0 = zeros(3,4);
reject_H0_2s0 = zeros(3,4);
reject_H0_1s1 = zeros(3,4);
reject_H0_2s1 = zeros(3,4);
T_all = cell(B,1);
VAR_all = cell(B,1);
for i = 1:4
    Beta_all{1,i} = zeros(B,p);
    Beta_d_all{1,i} = zeros(B,p);
    st_all{1,i} = zeros(B,p);
end
for i = 1:batch
    Beta_store{1,i} = zeros(B,p);
    df_store{1,i} = zeros(B,p);
end
for iter = 1:B
    [G_1,G_2,G_3] = G_index(p);
    %generate data
    G12s0_store_each = zeros(batch,2);
    G12s1_store_each = zeros(batch,2);
    G3_store_each = zeros(batch,3);
    Z = cell(batch,1);
    T_tilde = cell(batch,1);
    delta = cell(batch,1);
    Z1 = [];
    T_tilde1 =[];
    delta1 = [];
    T_all_each = zeros(12,4);
    VAR_all_each = zeros(12,4);
    rng((seed+p)*iter)
    for i = 1:batch
        [Z_0,T_tilde_0,delta_0] = g_data_tc(n,p,Beta_true,a);
        Z{i,1} = Z_0;
        T_tilde{i,1} = T_tilde_0;
        delta{i,1} = delta_0;
        Z1 = [Z1;Z{i,1}];
        T_tilde1 =[T_tilde1;T_tilde{i,1}];
        delta1 = [delta1;delta{i,1}];
    end
    %the first batch
    m = 1;
    lambda_best_1 = BIC_0(c0_min,c0_max,Z{1,1},T_tilde{1,1},delta{1,1},n,p,t);
    [Beta_best_1] = LASSO_0(Z{1,1},T_tilde{1,1},delta{1,1},n,lambda_best_1(1),p,t);
    [df_1,df2_1,df_each_1,df_times_1,H_1] = DF_new1(Beta_best_1,Z{1,1},T_tilde{1,1},delta{1,1},n,p);
    [temp1_1s0,temp2_1s0] = T_temp_new(delta{1,1},sigma0,m,n,p,seed,iter);
    [temp1_1s1,temp2_1s1] = T_temp_new(delta{1,1},sigma1,m,n,p,seed,iter);
    [s,Cov0,Cov_est] = Cov_threshold_new(T_tilde1,Z1,m,n,p,seed,iter);
    [temp4_1,temp5_1] = T_temp4_new(G_1,G_2,G_3,Beta_best_1,Z{m,1},T_tilde{m,1},delta{m,1},Cov_est,df2_1,n,p);

    %other batches
    H_sum = H_1;
    H_Beta = H_1*Beta_best_1;
    Beta_last = Beta_best_1;
    df_sum = df_1;
    df_times_sum = df_times_1;
    T_temp1s0_sum = temp1_1s0;
    T_temp2s0_sum = temp2_1s0;
    T_temp1s1_sum = temp1_1s1;
    T_temp2s1_sum = temp2_1s1;
    Z_times_sum = Z{1,1}'*Z{1,1};
    T_temp4_sum = temp4_1;
    T_temp5_sum = temp5_1;
    Beta_store{1,1}(iter,:) = Beta_last';
    df_store{1,1}(iter,:) = df_1';
    H_store{iter,1} = H_1;
    df_times_store{iter,1} = df_times_1;
    Cov_est_all{iter,1} = Cov_est;
    G12s0_store_each(m,:) = [temp1_1s0,temp2_1s0];
    G12s1_store_each(m,:) = [temp1_1s1,temp2_1s1];
    G3_store_each(m,:) = temp4_1;
    G4_store{iter,m} = temp5_1;
    for m=2:batch
        Beta_last_real = Beta_last;
        [lambda_best_2] = BIC_no(c0_min,c0_max,Z{m,1},T_tilde{m,1},delta{m,1},n,p,t,m,H_sum,Beta_last);
        [Beta_best_2] = LASSO(Z{m,1},T_tilde{m,1},delta{m,1},n,lambda_best_2(1),p,t,m,H_sum,Beta_last);
        [df_2,df2_2,df_each_2,df_times_2,H_2] = DF_new1(Beta_best_2,Z{m,1},T_tilde{m,1},delta{m,1},n,p);
        [temp1_2s0,temp2_2s0] = T_temp_new(delta{m,1},sigma0,m,n,p,seed,iter);
        [temp1_2s1,temp2_2s1] = T_temp_new(delta{m,1},sigma1,m,n,p,seed,iter);
        [s,Cov0,Cov_est] = Cov_threshold_new(T_tilde1,Z1,m,n,p,seed,iter);
        [temp4_2,temp5_2] = T_temp4_new(G_1,G_2,G_3,Beta_best_2,Z{m,1},T_tilde{m,1},delta{m,1},Cov_est,df2_2,n,p);
       
        H_sum = H_sum+H_2;
        df_sum = df_sum+df_2;
        df_times_sum = df_times_sum+df_times_2;
        H_Beta = H_Beta+H_2*Beta_best_2;
        T_temp1s0_sum = T_temp1s0_sum+temp1_2s0;
        T_temp2s0_sum = T_temp2s0_sum+temp2_2s0;
        T_temp1s1_sum = T_temp1s1_sum+temp1_2s1;
        T_temp2s1_sum = T_temp2s1_sum+temp2_2s1;
        Z_times_sum = Z_times_sum+Z{m,1}'*Z{m,1};
        T_temp4_sum = T_temp4_sum+temp4_2;
        T_temp5_sum = T_temp5_sum+temp5_2;
        Beta_last = Beta_best_2;%After ending, Beta_last = Beta_hat_final
        Beta_store{1,m}(iter,:) = Beta_last';
        df_store{1,m}(iter,:) = df_2';
        H_store{iter,m} = H_2;
        df_times_store{iter,m} = df_times_2;
        Cov_est_all{iter,m} = Cov_est;
        G12s0_store_each(m,:) = [temp1_2s0,temp2_2s0];
        G12s1_store_each(m,:) = [temp1_2s1,temp2_2s1];
        G3_store_each(m,:) = temp4_2;
        G4_store{iter,m} = temp5_2;
        if(m==batch)
            G12s0_store{iter,1} = G12s0_store_each;
            G12s1_store{iter,1} = G12s1_store_each;
            G3_store{iter,1} = G3_store_each;
        end
        %save 4 batches
        if(m==batch/4||m==batch/2||m==batch-batch/4||m==batch)
            k = m/4;
            H_sum_less = H_sum-H_2;
            H_Beta_less = H_Beta-H_2*Beta_best_2;
            df_sum_less = df_sum-df_2;
            Beta_all{1,k}(iter,:) = Beta_last';
            H_ave = H_sum/m;
            clear omega
            delete("data_H_sum\omega.mat")
            save H_ave.mat H_ave -v6
            Rpath = 'D:\Program Files\R\R-4.2.1\bin';
            filepath = 'D:\1mat\simulation\code\tc_online\data_H_sum\omega_compute.R';
            for ei = 1:100
                ei
                try
                    RunRcode(filepath,Rpath);
                    load('D:\1mat\simulation\code\tc_online\data_H_sum\omega.mat')
                end
                if exist('omega','var')
                    break
                end
            end

            omega_all{iter,k} = omega;
            Beta_d = Beta_last+omega'*(H_Beta_less-H_sum_less*Beta_last-df_sum)/m;
            Beta_d_all{1,k}(iter,:) = Beta_d';
            %
            l1_temp = norm(Beta_last-Beta_true,1);
            l2_temp = norm(Beta_last-Beta_true,2);
            l1(1,k) = l1(1,k)+l1_temp;
            l2(1,k) = l2(1,k)+l2_temp;
            l1_d_temp = norm(Beta_d-Beta_true,1);
            l2_d_temp = norm(Beta_d-Beta_true,2);
            l1_d(1,k) = l1_d(1,k)+l1_d_temp;
            l2_d(1,k) = l2_d(1,k)+l2_d_temp;
            var = omega'*df_times_sum*omega/(m*n);
            st = sqrt(diag(var));
            st_all{1,k}(iter,:) = st';
            cp_num = ((((Beta_d-t_alpha/sqrt(n*m)*st)<Beta_true)...
                +(Beta_true<(Beta_d+t_alpha/sqrt(n*m)*st)))==2);
            cp_sum(:,k) = cp_sum(:,k)+cp_num;
            
            U = U_n(Z{m,1},T_tilde{m,1},delta{m,1},n,p,m,H_sum_less,Beta_last_real,Beta_last);
            T_temp3 = omega*(H_Beta_less-H_sum_less*Beta_last_real-df_sum_less)/m;
            %s0
            %test1
            [T_1s0,VAR_1s0,T_2s0,VAR_2s0,T_3s0,VAR_3s0] = ...
                test1_new(G_1,G_2,G_3,Beta_last,omega,U,T_temp1s0_sum,T_temp2s0_sum,T_temp3,df_times_sum,n*m,sigma0);
            reject_H0_1s0(:,k) = reject_H0_1s0(:,k)+...
                [T_1s0>(sqrt(VAR_1s0)*z_alpha/sqrt(n*m));T_2s0>(sqrt(VAR_2s0)*z_alpha/sqrt(n*m));T_3s0>(sqrt(VAR_3s0)*z_alpha/sqrt(n*m))];
            %test2
            [T_4s0,VAR_4s0,T_5s0,VAR_5s0,T_6s0,VAR_6s0]=...          
                test2_new(G_1,G_2,G_3,Beta_last,omega,U,T_temp1s0_sum,T_temp2s0_sum,T_temp3,df_times_sum,n*m,sigma0,Cov0,Cov_est,T_temp4_sum,T_temp5_sum);
            reject_H0_2s0(:,k) = reject_H0_2s0(:,k)+...
                [T_4s0>(sqrt(VAR_4s0)*z_alpha/sqrt(n*m));T_5s0>(sqrt(VAR_5s0)*z_alpha/sqrt(n*m));T_6s0>(sqrt(VAR_6s0)*z_alpha/sqrt(n*m))];
            %s1
            %test1
            [T_1s1,VAR_1s1,T_2s1,VAR_2s1,T_3s1,VAR_3s1] = ...
                test1_new(G_1,G_2,G_3,Beta_last,omega,U,T_temp1s1_sum,T_temp2s1_sum,T_temp3,df_times_sum,n*m,sigma1);
            reject_H0_1s1(:,k) = reject_H0_1s1(:,k)+...
                [T_1s1>(sqrt(VAR_1s1)*z_alpha/sqrt(n*m));T_2s1>(sqrt(VAR_2s1)*z_alpha/sqrt(n*m));T_3s1>(sqrt(VAR_3s1)*z_alpha/sqrt(n*m))];
            %test2
            [T_4s1,VAR_4s1,T_5s1,VAR_5s1,T_6s1,VAR_6s1]=...          
                test2_new(G_1,G_2,G_3,Beta_last,omega,U,T_temp1s1_sum,T_temp2s1_sum,T_temp3,df_times_sum,n*m,sigma1,Cov0,Cov_est,T_temp4_sum,T_temp5_sum);
            reject_H0_2s1(:,k) = reject_H0_2s1(:,k)+...
                [T_4s1>(sqrt(VAR_4s1)*z_alpha/sqrt(n*m));T_5s1>(sqrt(VAR_5s1)*z_alpha/sqrt(n*m));T_6s1>(sqrt(VAR_6s1)*z_alpha/sqrt(n*m))];
            T_all_each(:,k) = [T_1s0;T_2s0;T_3s0;T_4s0;T_5s0;T_6s0;T_1s1;T_2s1;T_3s1;T_4s1;T_5s1;T_6s1];
            VAR_all_each(:,k) = [VAR_1s0;VAR_2s0;VAR_3s0;VAR_4s0;VAR_5s0;VAR_6s1;VAR_1s1;VAR_2s1;VAR_3s1;VAR_4s1;VAR_5s1;VAR_6s1];
        end
        m;
    end
    T_all{iter,1} = T_all_each;
    VAR_all{iter,1} = VAR_all_each;
    iter
end
st_mean = zeros(4,p);
st_empirical = zeros(4,p);
for i = 1:4
    Beta_mean(i,:) = mean(Beta_all{1,i});
    Beta_d_mean(i,:) = mean(Beta_d_all{1,i});
    st_mean(i,:) = mean(st_all{1,i});
    st_empirical(i,:) = std(Beta_d_all{1,i});
end
cp = cp_sum/B;
cp_mean = mean(cp);
cp_1 = mean(cp(1:3,:));
cp_2 = mean(cp(4:p,:));
l1_mean = l1/B;
l2_mean = l2/B;
l1_d_mean = l1_d/B;
l2_d_mean = l2_d/B;
reject_H0_1s0_mean = reject_H0_1s0/B;
reject_H0_2s0_mean = reject_H0_2s0/B;
reject_H0_1s1_mean = reject_H0_1s1/B;
reject_H0_2s1_mean = reject_H0_2s1/B;

% online
save 1Beta_all.mat Beta_all -v6 %store
save 1Beta_d_all.mat Beta_d_all -v6
save 1Beta_mean.mat Beta_mean -v6
save 1Beta_d_mean.mat Beta_d_mean -v6
save 1Beta_store.mat Beta_store -v6
save 1cp.mat cp -v6
save 1cp_1.mat cp_1 -v6
save 1cp_2.mat cp_2 -v6
save 1cp_mean.mat cp_mean -v6
save 1delta.mat delta -v6
save 1df_store.mat df_store -v6 %store
save 1H_store.mat H_store -v6 %store
save 1l1_d_mean.mat l1_d_mean -v6
save 1l2_d_mean.mat l2_d_mean -v6
save 1l1_mean.mat l1_mean -v6
save 1l2_mean.mat l2_mean -v6
save 1reject_H0_1s0_mean.mat reject_H0_1s0_mean -v6
save 1reject_H0_2s0_mean.mat reject_H0_2s0_mean -v6
save 1reject_H0_1s1_mean.mat reject_H0_1s1_mean -v6
save 1reject_H0_2s1_mean.mat reject_H0_2s1_mean -v6
save 1st_all.mat st_all -v6
save 1st_empirical.mat st_empirical -v6
save 1st_mean.mat st_mean -v6
save 1T_tilde.mat T_tilde -v6
save 1Z.mat Z -v6
save 1omega_all.mat omega_all -v6 %store
save 1df_times_store.mat df_times_store -v6 %store
save 1Cov_est_all.mat Cov_est_all -v6 %store
save 1G123_store.mat G12s0_store -v6 %store
save 1G4_store.mat G12s1_store -v6 %store
save 1G123_store.mat G3_store -v6 %store
save 1G4_store.mat G4_store -v6 %store
save 1T_all.mat T_all -v6
save 1VAR_all.mat VAR_all -v6
