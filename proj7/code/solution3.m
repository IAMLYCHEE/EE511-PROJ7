clear 
close all

clear 
%generate samples
%test1 : same weight well seperated spherical
%generate the first population
w = [0.5,0.5];
mu1 = [2;5];
mu2 = [0;1];
cov1 = [3,0;0,1/2];
cov2 = [1,0;0,2];
GMM2d2popuEM(mu1,mu2,cov1,cov2,w,300)

clear
%test2 : ellipsoidal covariance 
w = [0.5,0.5];
mu1 = [2;5];
mu2 = [0;1];
cov1 = [3,1;1,1/2];
cov2 = [1,1/2;1/2,2];
GMM2d2popuEM(mu1,mu2,cov1,cov2,w,300)

clear
%test3 : different weight
w = [0.2,0.8];
mu1 = [2;5];
mu2 = [0;1];
cov1 = [3,1;1,1/2];
cov2 = [1,1/2;1/2,2];
GMM2d2popuEM(mu1,mu2,cov1,cov2,w,300)


