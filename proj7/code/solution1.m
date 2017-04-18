clear
sampleAmount = 300;
cov = [3,-1,1;-1,5,3;1,3,4];
mu = [1;2;3];
sample = multiGauGenerator(cov,mu,sampleAmount);
scatter3(sample(:,1),sample(:,2),sample(:,3))
title('sample distribution')