clc;
clear all;
close all;

% This script plot the potential estimated by kernel density estimation of:
% (1) white Gaussian noise; (2) step function with white Gaussian noise.
%==========================================================================
% AUTHOR        Almog Lahav
% INSTITUTION   Technion
% DATE          23th August 2016
%
%
% SCRIPT PARAMETERS
%   runSim      0 = plot the potential calculated in previous simulation
%               1 = calculate the potential and then plot it
%==========================================================================


runSim = 0;

if(~runSim)
    
    open('figures/1well.fig')
    title('(a)')
    open('figures/2well.fig')
    title('(b)')

    
else
    
N = 2000;
sigma = 1;

y = sigma*randn(N,1);
[p,yi] = ksdensity(y);
U = -2*log(p);

figure;
plot(yi,U)
xlabel('y')
ylabel('U(y)')

xa = [.23 .33];
ya = [.6 .4];
annotation('arrow',xa,ya)
xa = 1 - xa + 0.04;
annotation('arrow',xa,ya)

sigma = 0.4;
y1 = [ones(round(N/2),1); -ones(round(N/2),1)] + sigma*randn(N,1);
[p,yi] = ksdensity(y1);
U = -2*log(p);

figure;
plot(yi,U)
xlabel('y')
ylabel('U(y)')

xa = [.28 .35];
ya = [.5 .3];
annotation('arrow',xa,ya)
xa = 1-xa+0.04;
annotation('arrow',xa,ya)
xa = [.495 .455];
ya = [.26 .21];
annotation('arrow',xa,ya)
xa = 1-xa+0.04;
annotation('arrow',xa,ya)

end