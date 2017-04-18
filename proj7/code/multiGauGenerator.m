function sample = multiGauGenerator(cov,mu,sampleAmount)
%use:  sample = multiGauGenerator(cov,mu,sampleAmount)
%generate multivariate gaussian random samples with given mean and
%covariance. 

%reshape the input
mu = reshape(mu,[length(mu),1]);
a = chol(cov,'lower');
% a * a';
aSize = size(cov);
sample = zeros(sampleAmount,aSize(1));
for i = 1 : sampleAmount
    x = a * randn(aSize(2),1) + mu;
    sample(i,:) = x';
end

