function plot2d2subpopuGMM(mu1,mu2,cov1,cov2,w,k,option)
%plot2d2subpopuGMM(mu1,mu2,cov1,cov2,w,k,option)
if option == 0
    figure
else
    hold on
    pause(0.5)
end
mu1 = reshape(mu1,[2,1]);
mu2 = reshape(mu2,[2,1]);
alpha = 4.7; %scale the plot
if mu1(1) < mu2(1)
    x1rangeLB = mu1(1)- alpha * sqrt(cov1(1,1));
    x1rangeUB = mu2(1)+ alpha * sqrt(cov2(1,1));
else
    x1rangeLB = mu2(1)- alpha * cov2(1,1);
    x1rangeUB = mu1(1)+ alpha * cov1(1,1);
end
x1 = x1rangeLB:0.05:x1rangeUB; 

if mu1(2) < mu2(2)
    x2rangeLB = mu1(2)- alpha * sqrt(cov1(2,2));
    x2rangeUB = mu2(2)+ alpha * sqrt(cov2(2,2));
else
    x2rangeLB = mu2(2)- alpha * sqrt(cov2(2,2));
    x2rangeUB = mu1(2)+ alpha * sqrt(cov1(2,2));
end 
x2 = x2rangeLB: 0.05 : x2rangeUB;

[X1,X2] = meshgrid(x1,x2);
f = w(1) * mvnpdf([X1(:) X2(:)],mu1',cov1)+w(2) * mvnpdf([X1(:) X2(:)],mu2',cov2);
f = reshape(f,length(x2),length(x1));
h = surf(x1,x2,f);
if k == 0
    title(' true GMM density');
else
    if k == 1
        title('initial guess');
    else
        title(strcat(num2str(k),' step estimation'));
    end
end
axis tight
view(0,90)
set(h,'LineStyle','none')