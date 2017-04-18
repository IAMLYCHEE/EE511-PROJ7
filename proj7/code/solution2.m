sampleAmount = 100000;
sample = zeros(sampleAmount,1);
for i = 1 : sampleAmount
    if rand(1) < 0.4
        sample(i) = randn(1) - 1;
    else
        sample(i) = randn(1) + 1;
    end
end
histogram(sample,'BinLimits',[-4,4],'Normalization','probability')
hold on
%theoretical
i = 1;
for t = -4: 0.05 : 4
    c = sqrt(2 * pi);
    p1 = exp(-((t+1)^2)/2);
    p2 = exp(-((t-1)^2)/2);
    ft(i) = 0.4 / c * p1 + 0.6 / c * p2; 
    i = i+1;
end

plot(-4 : 0.05 : 4, ft/9);
legend('histogram','theoretical')
title('mixture gaussian(100000 sample)')


    