function [T_1,VAR_1,T_2,VAR_2]=test1_new(G_1,G_2,Beta_best,omega,U,T_temp1_sum,T_temp2_sum,T_temp3,df_times_sum,N,sigma)

%Beta_G_1
Beta_G_1 = Beta_best;
Beta_G_1(G_1) = 0;
%Beta_G_2
Beta_G_2 = Beta_best;
Beta_G_2(G_2) = 0;


T_1 = Beta_G_1'*Beta_G_1-2*Beta_G_1'*omega*U+(2/N)*T_temp1_sum+2*Beta_G_1'*T_temp3;
VAR_1 = 4*Beta_G_1'*omega*(df_times_sum/N)*omega*Beta_G_1+4*(sigma^2)*T_temp2_sum/N;
T_2 = Beta_G_2'*Beta_G_2-2*Beta_G_2'*omega*U+(2/N)*T_temp1_sum+2*Beta_G_2'*T_temp3;
VAR_2 = 4*Beta_G_2'*omega*(df_times_sum/N)*omega*Beta_G_2+4*(sigma^2)*T_temp2_sum/N;
end