function GMM2d2popuEM(mu1,mu2,cov1,cov2,w,sampleAmount)
%GMM2d2popuEM(mu1,mu2,cov1,cov2,w,sampleAmount)


%generate samples
sample1 = multiGauGenerator(cov1,mu1,sampleAmount*w(1));
sample2 = multiGauGenerator(cov2,mu2,sampleAmount*w(2));
%sample plot
s1 = scatter(sample1(:,1),sample1(:,2));
s1.LineWidth = 0.5;
s1.MarkerEdgeColor = 'b';
s1.MarkerFaceColor = [0,0.5,0.5];
hold on
s2 = scatter(sample2(:,1),sample2(:,2));
s2.LineWidth = 0.5;
s2.MarkerEdgeColor = 'r';
s2.MarkerFaceColor = [0.5,0.5,0.0];
title(strcat(num2str(sampleAmount),' samples'));
hold off

disp('true GMM')
mu1
mu2
cov1
cov2
%use these samples perform EM test
plot2d2subpopuGMM(mu1,mu2,cov1,cov2,w,0,0);



%EM algorithm
%initialization
k = 2;
w = [1/k,1/k];
sample = [sample1;sample2];%we combine the samples so we donot know which sample comes from which cluster
%for the initial guess of mean ,we take the k-means algorithm
% mu1 = [mean(sample(1:sampleAmount*w(1),1)),mean(sample(1:sampleAmount*w(1),2))];
% mu2 = [mean(sample(sampleAmount*w(1)+1:sampleAmount,1)),mean(sample(sampleAmount*w(1)+1:sampleAmount,2))];
[y,C] = kmeans(sample,k);
mu1 = C(1,:);
mu2 = C(2,:);
mu = zeros(2,1,k);
mu(:,:,1) = mu1';
mu(:,:,2) = mu2';
covES = zeros(2,2,k);
covES(:,:,1) = eye(2);
covES(:,:,2) = eye(2);
%compute initial likelihood
sumLH = 0;
for i = 1 : sampleAmount
    sumLH = sumLH + log(w(1) * mvnpdf(sample(i,:),mu(:,:,1)',covES(:,:,1)) + ...
                w(2) * mvnpdf(sample(i,:),mu(:,:,2)',covES(:,:,2))); 
end
L = sumLH/sampleAmount;
L_new = 0;
ebsilon = 0.001;
%record data
disp('initial guess')
mu
covES
%plot the initial guess
plot2d2subpopuGMM(mu(:,:,1),mu(:,:,2),covES(:,:,1),covES(:,:,2),w,1,0);

%EM
stepCount = 1;
gamma = zeros(sampleAmount,k);
n = zeros(1,2);
while abs(L_new-L) >  ebsilon
%E-step
    L = L_new;
    for j = 1 : k
        for i = 1 : sampleAmount
            sumLHPoint = 0;
            for l = 1 : k
                sumLHPoint = sumLHPoint + w(l) * mvnpdf(sample(i,:),mu(:,:,l)',covES(:,:,l));
            end
            gamma(i,j)= w(j) * mvnpdf(sample(i,:),mu(:,:,j)',covES(:,:,j))/sumLHPoint;
        end
    end
    n(1) = sum(gamma(:,1));
    n(2) = sum(gamma(:,2));
    
%M-step
    w = n/sampleAmount;
    for j = 1 : k
        mu(:,:,j) = (gamma(:,j)' * sample /n(j) )';
        sumCov = zeros(2,2);
        for i = 1 : sampleAmount
            sumCov  = sumCov + gamma(i,j) * (sample(i,:)' - mu(:,:,j))*...
                         (sample(i,:)' - mu(:,:,j))';
        end
        covES(:,:,j)=sumCov / n(j);
    end
%the new likelihood
    
    stepCount = stepCount + 1;
    plot2d2subpopuGMM(mu(:,:,1),mu(:,:,2),covES(:,:,1),covES(:,:,2),w,stepCount,0);
    disp(strcat('step:',num2str(stepCount)))
    mu
    covES
    
%Convergence check
    sumLH = 0;
    for i = 1 : sampleAmount
        sumLH = sumLH + log(w(1) * mvnpdf(sample(i,:),mu(:,:,1)',covES(:,:,1)) + ...
                w(2) * mvnpdf(sample(i,:),mu(:,:,2)',covES(:,:,2))); 
    end
    L_new = sumLH/sampleAmount;
end