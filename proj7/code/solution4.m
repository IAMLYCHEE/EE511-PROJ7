clear
clc
close all
duration = load('duration.mat');
duration = duration.duration;
waiting = load('waiting.mat');
waiting = waiting.waiting;
scatter(duration,waiting);
title('without kmeans clustering')
xlabel('duration')
ylabel('waiting')


%cascade data points
samples = [duration,waiting];
%shuffle
samples = samples(randperm(size(samples,1)),:);
%kmeans 
[y,C] = kmeans(samples,2);

%figure the clusters
figure
hold on
scatter(samples(y == 1,1),samples(y==1,2),'x')
scatter(samples(y == 2,1),samples(y==2,2),'o')
plot(C(1,1),C(1,2),'rx','LineWidth',2)
plot(C(2,1),C(2,2), 'ro','LineWidth',2)
legend('Points of cluster 1','Points of cluster 2')
title('Data Points with Labels by K-means Clustering')
xlabel('duration')
ylabel('waiting')
hold off

figure
GMM2d2popuEMcore(samples,0)
hold on
scatter3(samples(y == 1,1),samples(y==1,2),ones(length(samples(y==1,1)),1),'rx')
scatter3(samples(y == 2,1),samples(y==2,2),ones(length(samples(y==2,1)),1),'ro')
plot(C(1,1),C(1,2),'rx','LineWidth',2)
plot(C(2,1),C(2,2), 'ro','LineWidth',2)
% legend('Points of cluster 1','Points of cluster 2')
title('Data Points with Labels by K-means Clustering & EM GMM estimation')
xlabel('duration')
ylabel('waiting')
hold off
