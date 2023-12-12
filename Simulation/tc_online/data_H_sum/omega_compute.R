library(R.matlab)
library(fastclime)
nl<-100 
#nlambda
p<-400
path <- ("D:/1mat/1tc_online_1/")
pathname <- file.path(path, "H_ave.mat")
H<-readMat(pathname)
H_compute<-H$H.ave
#out1 = fastclime(H_compute,lambda.min = 1e-6,nlambda = 800)
#out2 = fastclime.selector(out1$lambdamtx, out1$icovlist,1e-5)
out1 = fastclime(H_compute,lambda.min = 1e-5,nlambda = 120)
out2 = fastclime.selector(out1$lambdamtx, out1$icovlist,1e-4)
omega<-out2[["icov"]]
writeMat("D:/1mat/1tc_online_1/data_H_sum/omega.mat",omega = omega)