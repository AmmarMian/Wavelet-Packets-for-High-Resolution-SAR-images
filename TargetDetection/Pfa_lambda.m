%% Description
% Simulation of Pfa-lambda
%% Specifiactions
% Other m-files required: target creation functions
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
    clear variables
    clc
    close all
    set(0,'defaulttextinterpreter','latex')

%% Reading Image
    addpath('..', '../Hyperimages', '../Target Dictionary', '../Detectors')
    imagePath = '../../../../../Data/SANDIA/MiniSAR20050519p0010image002.mat';
    [Image, Header] = readSandia(imagePath, true, [1, 1, 500, 500]);
%     [Image, Header] = readSDMS('../../../../..\Data\SDMS', 'FP0121', true, [800, 800, 1201, 1201]);
%     Image = Image(:,:,3);
%% Decomposition
    Decomposition = struct;
    Decomposition.method = 'Exponential';
    Decomposition.J = 1;
    Decomposition.numberOfBands = 5;
    Decomposition.numberOfLooks = 5;
    N_k = Decomposition.numberOfBands^Decomposition.J;
    N_theta = Decomposition.numberOfLooks^Decomposition.J;
    Decomposition.d1 = inf;
    Decomposition.d2 = inf;
    Decomposition.visualise = false;
    [ImageFTheta, k, theta] = decomposeFTheta(Image, Header, Decomposition);  
    ImageFTheta = ImageFTheta(1:N_theta:end, 1:N_k:end,:);

%% Detectors
    DetectorsStruct = initTargetDetectors();
    p = randn(N_k*N_theta,1) + 1i*randn(N_k*N_theta,1);
    p = p/(p'*p);
    Results = computeTestOnImage(ImageFTheta, p, DetectorsStruct);
    
%% Plotting
    
    close all
    for d=1:DetectorsStruct.number
        figure('Position', [100,100,350,300])
        Nsimu = length(Results{d}(:));
        Pfa = (Nsimu:-1:1)/ Nsimu;
        lambda = sort(Results{d}(:));
        loglog(lambda, Pfa, 'k', 'linewidth', 2, 'linestyle', '--')
        hold on
        N = sum(DetectorsStruct.mask(:));
        m = N_k*N_theta;
        loglog(lambda, DetectorsStruct.list(d).theo(lambda, N,m), 'k', 'linewidth', 2, 'linestyle', '-')
%         title(DetectorsStruct.list(d).name, 'interpreter', 'latex')
        h = legend('Experimental', 'Theoritical');
        set(h, 'interpreter', 'latex', 'location', 'southwest')
        grid
        grid minor
        xlabel('$\lambda$', 'interpreter', 'latex')
        ylabel('$P_{FA}$', 'interpreter', 'latex')
        ylim([10/Nsimu, 1])
        textString = {sprintf('$j=%d$, $R=%d$, $L=%d$', Decomposition.J, Decomposition.numberOfBands, Decomposition.numberOfLooks), sprintf('$K=%d$', sum(DetectorsStruct.mask(:)))};
    annotation('textbox',...
    [0.62, 0.60, 0.1, 0.1], ...
    'Units','characters',...
    'String', textString,...
    'backgroundcolor', 'w', ...
    'FitBoxToText','on',...
    'interpreter', 'latex');
    end
    
%% ------------- END CODE ----------------
