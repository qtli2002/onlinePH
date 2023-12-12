function [T_4,VAR_4,T_5,VAR_5]=test2_new(G_1,G_2,Beta_best,omega,U,T_temp1_sum,T_temp2_sum,T_temp3,df_times_sum,N,sigma,Cov0,Cov_est,T_temp4_sum,T_temp5_sum)
%Beta_G_1: {1,2,3}
Beta_G_1 = Beta_best;
Beta_G_1(G_1) = 0;
Cov_G_1 = Cov0; Cov_G_1(:,G_1) = 0; Cov_G_1(G_1,:) = 0;
Sigma_tilde_1 = Cov_est; Sigma_tilde_1(:,G_1) = 0; Sigma_tilde_1(G_1,:) = 0;
%Beta_G_2: {1,0,0,49*other}
Beta_G_2 = Beta_best;
Beta_G_2(G_2) = 0;
Cov_G_2 = Cov0; Cov_G_2(:,G_2) = 0; Cov_G_2(G_2,:) = 0;
Sigma_tilde_2 = Cov_est; Sigma_tilde_2(:,G_2) = 0; Sigma_tilde_2(G_2,:) = 0;


T_4 = Beta_G_1'*Cov_G_1*Beta_G_1-2*Beta_G_1'*Sigma_tilde_1*omega*U...
    +(2/N)*T_temp1_sum+2*Beta_G_1'*Sigma_tilde_1*T_temp3;
VAR_4 = 4*Beta_G_1'*Sigma_tilde_1*omega*(df_times_sum/N)*omega*(Beta_G_1'*Sigma_tilde_1)'...
    +4*(sigma^2)*T_temp2_sum/N...
    +T_temp4_sum(1,1)/N-4*Beta_G_1'*Sigma_tilde_1*omega*T_temp5_sum(:,1)/N;
T_5 = Beta_G_2'*Cov_G_2*Beta_G_2-2*Beta_G_2'*Sigma_tilde_2*omega*U...
    +(2/N)*T_temp1_sum+2*Beta_G_2'*Sigma_tilde_2*T_temp3;
VAR_5 = 4*Beta_G_2'*Sigma_tilde_2*omega*(df_times_sum/N)*omega*(Beta_G_2'*Sigma_tilde_2)'...
    +4*(sigma^2)*T_temp2_sum/N...
    +T_temp4_sum(1,2)/N-4*Beta_G_2'*Sigma_tilde_2*omega*T_temp5_sum(:,2)/N;