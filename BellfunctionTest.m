%% Description
% Bell functions test
%% Specifiactions
% MAT-files required: none
%% Authors
% Authors: Ammar Mian, Ph.D., Signal processing at Supelec SONDRA
% Email address: ammar.mian@centralesupelec.fr
%% Copyright
% Copyright 2017 CentraleSupelec
%
% Licensed under the Apache License, Version 2.0 (the "License");
% you may not use this file except in compliance with the License.
% You may obtain a copy of the License at
%
%     http://www.apache.org/licenses/LICENSE-2.0
%
% Unless required by applicable law or agreed to in writing, software
% distributed under the License is distributed on an "AS IS" BASIS,
% WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
% See the License for the specific language governing permissions and
% limitations under the License.
%% ------------- BEGIN CODE --------------
    close all
    clc
    clear variables
    addpath('..', './Hyperimages')
%% Plot gbell function
    x = [-10:0.1:10];
    a=3;
    c=0;
    figure('Position', [100,100,700,300])
    linestyles = {'--', ':', '-.', '-'};
    markerStyle = {'+', 'o', '*', '.', 's', 'd', '^', 'p', 'h'};
    i=0;
    lg = {};
    for b=[1,10,100,inf]
        i = i+1;
        plot(x, gbellmfunction(x, [a,b,c]), 'k', 'linestyle', linestyles{i})
        hold on
%         plot(x, gbellmfunction(x, [a,b,c+2*a]), 'b', 'linewidth', 1)
%         plot(x, gbellmfunction(x, [a,b,c-2*a]), 'r', 'linewidth', 1)
        lg = cat(2,lg, sprintf('$b=%d$', b));
    end
    h = legend(lg);
    set(h, 'interpreter', 'latex')
%     grid
%     grid minor
    xlabel('$x$', 'interpreter', 'latex')
    ylabel(sprintf('$g^{\\mathrm{Bell}}_{%d, b, %d}(x)$', a, c), 'interpreter', 'latex')
%%     
% % %% Derivative of Bell function
% % 
% %     x = [-10:0.1:10];
% %     a=3;
% %     c=5;
% %     i=0;
% %     lg = {};
% %     figure('Position', [100,100,700,300])
% %     for b=10^10
% %         i = i+1;
% %         plot(x, dgbell(x, a,b,c), 'linewidth', 2)
% %         hold on
% %         lg = cat(2,lg, sprintf('$b=%d$', b));
% %     end
% %     h = legend(lg);
% %     set(h, 'interpreter', 'latex')
% %     grid
% %     grid minor
% %     xlabel('$x$', 'interpreter', 'latex')
% %     ylabel(sprintf('$\\frac{\\mathrm{d}g^{\\mathrm{Bell}}_{%d, b, %d}(x)}{\\mathrm{d}x}$', a, c), 'interpreter', 'latex')
% %     
% %% Bell function in context of SAR geometry
%     
%     c = 3e8;
%     f0 = 9e9;
%     B = 5e9;
%     k_0 = 2*f0/c;
%     k_B = 2*B/c;
%     k = linspace(k_0-k_B,k_0+k_B,10000);
%     R = 3;
%     j = 1;
% %   m = 0;
%     lg = {};
%     figure('Position', [100,100,700,300])
%     for m=0:R^j-1
% 
%         for d1=[inf]
%             i = i+1;
%             plot(k, Hk(k, j, m, d1, R, k_0, k_B), 'linewidth', 2)
%             hold on
%             lg = cat(2,lg, sprintf('$d_1=%d$', d1));
%         end
% %         if m==0
% %             h = legend(lg);
% %             set(h, 'interpreter', 'latex')
% %         end
%         grid
%         grid minor
%         xlabel('$\kappa$', 'interpreter', 'latex')
%     end
% % 
% %     
% % %% Derivative Bell function in context of SAR geometry
% % 
% %     c = 3e8;
% %     f0 = 9.6e9;
% %     B = 5e9;
% %     k_0 = 2*f0/c;
% %     k_B = 2*B/c;
% %     k = linspace(k_0-k_B/2,k_0+k_B/2,10000);
% %     R = 10;
% %     j = 1;
% %     m = 0;
% %     lg = {};
% %     figure('Position', [100,100,700,300])
% %     for d1=[1,10, inf]
% %         i = i+1;
% %         plot(k, dHk(k, j, m, d1, R, k_0, k_B), 'linewidth', 2)
% %         hold on
% %         lg = cat(2,lg, sprintf('$d_1=%d$', d1));
% %     end
% %     h = legend(lg);
% %     set(h, 'interpreter', 'latex')
% %     grid
% %     grid minor
% %     xlabel('$\kappa$', 'interpreter', 'latex')
% %     
% 
%     
%     
%% Redundancy Bell function in context of SAR geometry waveabs
    
    c = 3e8;
    f0 = 9.6e9;
    B = 640e6;
    k_0 = 2*f0/c;
    k_B = 2*B/c;
    k = linspace(k_0-k_B/2,k_0+k_B/2,10000);
    R = 2;
    j = 1;
    lg = {};
    figure('Position', [100,100,700,300])
    i = 0;
    for d1=[1,3,10,100,inf]
        i = i+1;
        Q = 0;
        for m=0:R^j-1
            Q = Q +  Hk(k, j, m, d1, R, k_0, k_B).^2;
        end
        pl = plot(k,Q,'linestyle', '--', 'marker',markerStyle{i});
        pl.MarkerIndices = floor(linspace(1,10000,30));
        hold on
        lg = cat(2,lg, sprintf('$d_1=%d$', d1));
    end
    h = legend(lg);
    set(h, 'interpreter', 'latex', 'fontsize', 12)
    grid
    grid minor
    xlabel('$\kappa$', 'interpreter', 'latex', 'fontsize', 18)
    ylabel('$Q_\kappa$', 'interpreter', 'latex', 'fontsize', 18)
    
    
%% Redundancy Bell function in context of SAR geometry theta
    
    theta_B = 0.25;
    theta = linspace(-theta_B, theta_B, 10000);
    L = 2;
    j = 1;
    lg = {};
    figure('Position', [100,100,700,300])
    i = 0;
    for d2=[1,3,10,100,inf]
        i = i+1;
        Q = 0;
        for n=0:L^j-1
            Q = Q + Htheta(theta, j, n, d2,L, theta_B).^2;
        end
        pl = plot(theta,Q,'linestyle', '--', 'marker',markerStyle{i});
        pl.MarkerIndices = floor(linspace(1,10000,30));
        hold on
        lg = cat(2,lg, sprintf('$d_2=%d$', d2));
    end
    h = legend(lg);
    set(h, 'interpreter', 'latex', 'fontsize', 12)
    grid
    grid minor
    xlabel('$\theta$', 'interpreter', 'latex', 'fontsize', 18)
    ylabel('$Q_{\theta}$', 'interpreter', 'latex', 'fontsize', 18)
    
    
%% Conjoint function SDMS
    c = 3e8;
    f0 = 9.6e9;
    B = 640*10^6;
    theta_B = 0.25;
    theta = linspace(-theta_B, theta_B, 1000);
    k_0 = 2*f0/c;
    k_B = 2*B/c;
    k = linspace(k_0-k_B/2,k_0+k_B/2,1000);
    [K,Theta] = meshgrid(k,theta);
    R = 5;
    L = 5;
    j = 1;
    
    maxQ = -inf;
    minQ = inf;
    d1_vector=[0.5,1,3,10,50];
    i = 1;
    Q = cell(length(d1_vector),1);
    for d1=d1_vector
        Q{i} = 0;
        for m=0:R^j-1
            for n=0:L^j-1
                Q{i} = Q{i} + (Hk(K, j, m, d1, R, k_0, k_B).^2).*(Htheta(Theta, j, n, d1,L, theta_B).^2);
            end
        end
        maxQ = max(maxQ, max(max(Q{i})));
        minQ = min(minQ, min(min(Q{i})));
        i = i+1;
    end
 %%   
    i = 1;
    for d1=d1_vector
        figure
        imagesc(k,theta,Q{i}, [minQ, maxQ])
        colorbar
        title(sprintf('$d_1=d_2=%d$',d1), 'interpreter', 'latex')
        xlabel('$\kappa$', 'interpreter', 'latex', 'fontsize', 18)
        ylabel('$\theta$', 'interpreter', 'latex', 'fontsize', 18)
        colormap(jet(100))
        i = i+1;
    end    
    
%% Functions
function y = dgbell(x, a, b, c)
    y = -2*b/a*((x-c)/a).^(2*b-1) .* gbellmfunction(x, [a,b,c]);
end


function y = Hk(x, j, m, d1, R, k_0, k_B)
    a = k_B / (2 * R^j)  ;
    b = d1;
    c = k_0 - k_B/2 + (2*m+1)*a;
    y = gbellmfunction(x, [a,b,c]).*(abs(x-k_0)<=k_B/2);
end

function y = Htheta(x, j, n, d2, L, theta_B)
    a = theta_B / (L^j)  ;
    b = d2;
    c =  - theta_B + (2*n+1)*a;
    y = gbellmfunction(x, [a,b,c]);
end

function y = dHk(x, j, m, d1, R, k_0, k_B)
    a = k_B / (2 * R^j)  ;
    b = d1;
    c = k_0 - k_B/2 + (2*m+1)*a;
    y = R^(j/2)*dgbell(x, a,b,c);
end