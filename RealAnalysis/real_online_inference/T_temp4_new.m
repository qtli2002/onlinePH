function [temp4,temp5] = T_temp4_new(G_1,G_2,Beta_best,Z_j,T_tilde_j,delta_j,Cov_est_j,df2_j,n,p)
%Beta_G_1
Beta_G_1 = Beta_best;
Beta_G_1(G_1) = 0;
Z_G_1 = Z_j;
Z_G_1(:,G_1) = 0;
Sigma_tilde_1 = Cov_est_j; Sigma_tilde_1(:,G_1) = 0; Sigma_tilde_1(G_1,:) = 0;
%Beta_G_2
Beta_G_2 = Beta_best;
Beta_G_2(G_2) = 0;
Z_G_2 = Z_j;
Z_G_2(:,G_2) = 0;
Sigma_tilde_2 = Cov_est_j; Sigma_tilde_2(:,G_2) = 0; Sigma_tilde_2(G_2,:) = 0;
temp4_G1 = 0; temp4_G2 = 0;
temp5_G1 = zeros(p,1);temp5_G2 = zeros(p,1);
for i=1:n
    temp_G1 =  Beta_G_1'*(T_tilde_j(i,:)*Z_G_1(i,:)'*Z_G_1(i,:)-Sigma_tilde_1)*Beta_G_1;
    temp_G2 =  Beta_G_2'*(T_tilde_j(i,:)*Z_G_2(i,:)'*Z_G_2(i,:)-Sigma_tilde_2)*Beta_G_2;
    temp4_G1 = temp4_G1+temp_G1^2;
    temp4_G2 = temp4_G2+temp_G2^2;
    temp5_G1 = temp5_G1+(Z_G_1(i,:)'-df2_j)*delta_j(i,:)*temp_G1;
    temp5_G2 = temp5_G2+(Z_G_2(i,:)'-df2_j)*delta_j(i,:)*temp_G2;
end
temp4 = [temp4_G1,temp4_G2];
temp5 = [temp5_G1,temp5_G2];