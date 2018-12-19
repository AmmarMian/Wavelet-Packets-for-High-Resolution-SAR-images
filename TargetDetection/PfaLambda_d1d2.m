%% Description
% Simulation of Pfa-lambda for d1 and d2 that vary
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
    
%% Simulation Params
    Decomposition = struct;
    Decomposition.method = 'Exponential';
    Decomposition.J = 1;
    Decomposition.numberOfBands = 5;
    Decomposition.numberOfLooks = 5;
    Decomposition.visualise = false;
    N_k = Decomposition.numberOfBands^Decomposition.J;
    N_theta = Decomposition.numberOfLooks^Decomposition.J;
    DetectorsStruct = initTargetDetectors();
    p = randn(N_k*N_theta,1) + 1i*randn(N_k*N_theta,1);
    p = p/(p'*p);
    
    d1_vec = [0.1,0.5,1,10,inf];
    d2_vec = d1_vec;
%% Reading Image
    addpath('..', '../Hyperimages', '../Target Dictionary', '../Detectors')
    imagePath = '../../../../Data/SANDIA/MiniSAR20050519p0010image002.mat';
    [Image, Header] = readSandia(imagePath, true, [1, 1, 400, 400]);
%     [Image, Header] = readSDMS('C:\Users\mian_amm\Dropbox\Thèse\Data\SDMS', 'FP0121', false, [100, 500, 1100, 1500]);
%     Image = Image(:,:,3);


%% Pfa-lambda
    lambda = cell(length(d1_vec),DetectorsStruct.number);
    for i_d1 = 1:length(d1_vec)
        Decomposition.d1 = d1_vec(i_d1);
        Decomposition.d2 = d2_vec(i_d1);
        % Decomposition
        [ImageFTheta, k, theta] = decomposeFTheta(Image, Header, Decomposition);  
        ImageFTheta = ImageFTheta(1:N_theta:end, 1:N_k:end,:);
        % Detection
        Results = computeTestOnImage(ImageFTheta, p, DetectorsStruct);
        for d=1:DetectorsStruct.number
            lambda{i_d1, d} = sort(Results{d}(:));
            Nsimu = length(lambda{i_d1, d});
        end
    end
    
%% Plotting
    close all
    numberOfPoints = 15;
    Pfa = (Nsimu:-1:1)/ Nsimu;
    index = floor(logspace(1,log10(Nsimu),numberOfPoints));
    index = unique(floor(index(:)));
    markers = {'s', 'd', 'o', '*', '+'};
    for d=1:DetectorsStruct.number
        figure('Position', [100, 100, 700, 300]);
        lg = {};
        for i_d1 = 1:length(d1_vec)
            pl = loglog(lambda{i_d1, d}, Pfa, 'k', 'linestyle', '--', 'marker', markers{i_d1}, 'markersize', 10);
            pl.MarkerIndices = Nsimu-index+1;
            hold on
            lg = cat(2, lg, sprintf('$d_1=%d$, $d_2=%d$', d1_vec(i_d1), d2_vec(i_d1))); 
        end
        title(DetectorsStruct.list(d).name, 'interpreter', 'latex')
        h = legend(lg);
        set(h, 'interpreter', 'latex', 'location', 'southwest')
        grid
        grid minor
        xlabel('$\lambda$', 'interpreter', 'latex')
        ylabel('$P_{FA}$', 'interpreter', 'latex')
    end
    
%% ------------- END CODE ----------------
