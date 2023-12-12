function [s,Cov0,Cov_est]=Cov_threshold_new(T_tilde,Z,batch,N_m,p,seed)
n1 = round(N_m*(1-1/log(N_m)));
n2 = N_m-n1;
T_tilde0 = T_tilde(1:N_m,:);
Z0 = Z(1:N_m,:);
N = 10;
s_min = log(p)/(N_m);
rng(batch+p+seed)
R_s_all = [];
for s_each=s_min:0.001:(s_min+0.01)
    R_s = 0;
    for i = 1:N
        n1_index = randperm(N_m,n1);
        n2_index = setdiff(1:(N_m),n1_index);
        T_tilde_n1 = T_tilde0(n1_index,:);
        Z_n1 = Z0(n1_index,:);
        T_tilde_n2 = T_tilde0(n2_index,:);
        Z_n2 = Z0(n2_index,:);
        Cov_n2 = Z_n2'*(Z_n2.*repmat(T_tilde_n2,1,p))/n2;
        Cov_n1 = Z_n1'*(Z_n1.*repmat(T_tilde_n1,1,p))/n1;
        T_n1 = zeros(p)+(abs(Cov_n1)>=s_each).*Cov_n1;
        r_s = trace((T_n1-Cov_n2)*((T_n1-Cov_n2)'));
        R_s = R_s+r_s;
    end
    R_s = R_s/N;
    R_s_all = [R_s_all,R_s];
end
index = find(R_s_all==min(R_s_all));
s = s_min+0.001*(index(1)-1);
Cov0 = Z0'*(Z0.*repmat(T_tilde0,1,p))/(N_m);
Cov_est = zeros(p)+(abs(Cov0)>=s).*(Cov0);    
