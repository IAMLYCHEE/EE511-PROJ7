function GMM2d2popuEMcore(sample,option)
sampleAmount = size(sample,1);
k = 2;
w = [1/k,1/k];
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
plot2d2subpopuGMM(mu(:,:,1),mu(:,:,2),covES(:,:,1),covES(:,:,2),w,1,option);

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
    plot2d2subpopuGMM(mu(:,:,1),mu(:,:,2),covES(:,:,1),covES(:,:,2),w,stepCount,option);
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