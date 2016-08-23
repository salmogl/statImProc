clc;
clear all;
close all;

% This script run experiment of 3 different kernel density estimation: 
%                  (1) fixed kernel. (equivalent to rows normalization)
%                  (2) normalized kernel (equivalent to columns and rows normalization)
%                  (3) addaptive kernel (equivalent to Addaptive KDE)
%==========================================================================
% AUTHOR        Almog Lahav
% INSTITUTION   Technion
% DATE          23th August 2016
%==========================================================================
    
N = 4000; % samples
sigma = [ 1 0.1 0.2 0.3];
mu = [-3 3 4 5];

% Mu = kron(mu,ones(1,N));
% Sigma = kron(sigma,ones(1,N));
class = randi(length(mu),[N 1]);
y = zeros(1,N);
for i = 1:length(mu)
    classSize = length(find(class == i));
    y(class == i) = mu(i) + sigma(i).*randn(1,classSize) ;
end
    
% y = Mu + Sigma.*randn(1,length(mu)*N) ;
[val,idxSort] = sort(y);

[p,yi] = ksdensity(y);
U = -2*log(p);

realP = zeros(1,length(y));
realY = linspace(-6,6,length(y));
for i = 1:length(mu)
    realP = realP + (2*pi*sigma(i)^2)^-0.5*exp(-(realY-mu(i)).^2/(2*sigma(i)^2));
end
realP = realP/length(mu);


figure;
subplot(221)
plot(yi,p)
hold on 
plot(realY,realP,'r')
xlabel('y')
ylabel('p(y)')
title('Density Using "ksdensity"')

args.knn = 50;
args.eps = 0.1;
args.norm = 'rows';
p_noNorm = kde(y,args);

args.norm = 'columns';
p_norm = kde(y,args);

args.norm = 'adaptive';
p_adapt = kde(y,args);

subplot(222)
plot(y(idxSort),p_noNorm)
hold on 
plot(realY,realP,'r')
title('D')

subplot(223)
plot(y(idxSort),p_norm)
hold on 
plot(realY,realP,'r')
title('Dr')

subplot(224)
plot(y(idxSort),p_adapt)
hold on 
plot(realY,realP,'r')
title('Adaptive')

args.norm = 'rows';
args.eps = 1;
p_noNorm = kde(y,args);

args.norm = 'columns';
p_norm = kde(y,args);

args.norm = 'adaptive';
p_adapt = kde(y,args);

figure
subplot(121)
plot(y(idxSort),p_noNorm)
hold on 
plot(realY,realP,'r')
h=legend('$\hat{P}_0(x)$','$f(x)$');
set(h,'Interpreter','latex','Fontsize',15)
title('The Density Estimator $\hat{P}_0(x)$ For $\alpha=1$','Interpreter','latex','Fontsize',20)
xlabel('$x$','Interpreter','latex','Fontsize',20)

subplot(122)
plot(y(idxSort),p_norm)
hold on 
plot(realY,realP,'r')
h=legend('$\hat{P}_1(x)$','$f(x)$');
set(h,'Interpreter','latex','Fontsize',15)
title('The Density Estimator $\hat{P}_1(x)$ For $\alpha=1$','Interpreter','latex','Fontsize',20)
xlabel('$x$','Interpreter','latex','Fontsize',20)


figure
plot(y(idxSort),p_adapt)
hold on 
plot(realY,realP,'r')
title('Adaptive')


realP0 = zeros(1,length(y));
for i = 1:length(mu)
    realP0 = realP0 + (2*pi*sigma(i)^2)^-0.5*exp(-(y-mu(i)).^2/(2*sigma(i)^2));
end
realP0 = realP0/length(mu);
realP0 = realP0(idxSort);


eps = linspace(0.005,0.6,5);
iter = 1;

err_norm = zeros(length(eps),iter);
err_noNorm = zeros(length(eps),iter);
err_adapt = zeros(length(eps),iter);
mErr_norm = zeros(length(eps),iter);
mErr_noNorm = zeros(length(eps),iter);
mErr_adapt = zeros(length(eps),iter);


for i = 1:length(eps)
    for j = 1:iter
        
        args.eps = eps(i);
        args.norm = 'rows';
        p_noNorm = kde(y,args);
        args.norm = 'columns';
        p_norm = kde(y,args);
        args.norm = 'adaptive';
        p_adapt = kde(y,args);
        
        err_norm(i,j) = mean((p_norm-realP0).^2); 
        err_noNorm(i,j) = mean((p_noNorm-realP0).^2);
        err_adapt(i,j) = mean((p_adapt-realP0).^2);
        mErr_norm(i,j) = max((p_norm-realP0).^2); 
        mErr_noNorm(i,j) = max((p_noNorm-realP0).^2);
        mErr_adapt(i,j) = max((p_adapt-realP0).^2);

        
    end
end

AIMSE_norm = mean(err_norm,2);
AIMSE_noNorm = mean(err_noNorm,2);
% AIMSE_adapt = mean(err_adapt,2);

MMSE_norm = mean(mErr_norm,2);
MMSE_noNorm = mean(mErr_noNorm,2);
% MMSE_adapt = mean(mErr_adapt,2);

figure;
subplot(121)
plot(eps,AIMSE_norm,'r')
hold on
plot(eps,AIMSE_noNorm,'b')
% plot(eps,AIMSE_adapt,'k')
h=legend('$\hat{P}_1$','$\hat{P}_0$','$\hat{P}_A$');
set(h,'Interpreter','latex','Fontsize',15)
title('AIMSE','Interpreter','latex','Fontsize',15)
xlabel('$\alpha$','Interpreter','latex','Fontsize',15)
ylabel('Error','Interpreter','latex','Fontsize',15)
xlim([0 0.6])

subplot(122)
plot(eps,MMSE_norm,'r')
hold on
plot(eps,MMSE_noNorm,'b')
% plot(eps,MMSE_adapt,'k')
h=legend('$\hat{P}_1$','$\hat{P}_0$','$\hat{P}_A$');
set(h,'Interpreter','latex','Fontsize',15)
title('MMSE','Interpreter','latex','Fontsize',15)
xlabel('$\alpha$','Interpreter','latex','Fontsize',15)
ylabel('Error','Interpreter','latex','Fontsize',15)
xlim([0 0.6])




