library(R.matlab)
library(fastclime)
nl<-100 
path <- ("D:/1mat/real_data_analysis/code/real_online_inference/")
pathname <- file.path(path, "H_ave.mat")
H<-readMat(pathname)
H_compute<-H$H.ave

set.seed(2023)
out1 = fastclime(H_compute,lambda.min = 1e-6,nlambda =2000)
set.seed(2023)
out2 = fastclime.selector(out1$lambdamtx, out1$icovlist,1e-5)

omega<-out2[["icov"]]
writeMat("D:/1mat/real_data_analysis/code/real_online_inference/data_H_sum/omega.mat",omega = omega)