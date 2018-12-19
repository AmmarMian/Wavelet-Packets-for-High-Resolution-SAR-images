%% Description
% Visualisation of decompositon according to d1_d2
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
%     [Image, Header] = readSDMS('C:\Users\mian_amm\Dropbox\Thèse\Data\SDMS', 'FP0120', true, [100, 500, 1100, 1500]);
%     Image = Image(:,:,1);
    
%% Target Creation

%     N_k = 5;
%     N_theta = 5;
%     steeringVector1 = randn(N_theta,N_k) + 1i*randn(N_theta,N_k);
% %     steeringVector = ones(N_theta,N_k);
% %     steeringVector = reshape([1:13,12:-1:1], N_theta,N_k);
%     Target1 = struct;
%     Target1.position = [381, 336];
%     Target1.Args = {steeringVector1, N_k, N_theta};
%     Target1.SNR_dB = 20;
%     Target1.function = @steeringVectorTarget;
%     
%     steeringVector2 = randn(N_theta,N_k) + 1i*randn(N_theta,N_k);
%     Target2 = struct;
%     Target2.position = [261, 261];
%     Target2.Args = {steeringVector2, N_k, N_theta};
%     Target2.SNR_dB = 70;
%     Target2.function = @steeringVectorTarget;
    
    
    gaussTarget = struct;
    gaussTarget.position = [251, 251];
    mean = [112, 0];
    scale = [1, 0.01];
    gaussTarget.Args = {mean, scale};
    gaussTarget.SNR_dB = 30;
    gaussTarget.function = @gaussianTarget;
%     
%     sTarget = struct;
%     sTarget.position = [250, 300];
%     mean = [116, 0];
%     scale = [0.5, 0.01];
%     sTarget.Args = {mean, scale};
%     sTarget.SNR_dB = 10;
%     sTarget.function = @sincTarget;

    TargetStruct = struct;
    TargetStruct.list = [gaussTarget];
    TargetStruct.number = length(TargetStruct.list);
    TargetStruct.Mx = 20;
    TargetStruct.My = 20;
    
    ImageTargets = createTargets(Image, Header, TargetStruct);


%% Decomposition
    d1_vec = [3,10,100,inf];
    for i_d1 = 1:length(d1_vec)
        Decomposition = struct;
        Decomposition.method = 'Exponential';
        Decomposition.J = 1;
        Decomposition.numberOfBands = 2;
        Decomposition.numberOfLooks = 2;
        N_k = Decomposition.numberOfBands^Decomposition.J;
        N_theta = Decomposition.numberOfLooks^Decomposition.J;
        Decomposition.d1 = d1_vec(i_d1);
        Decomposition.d2 = d1_vec(i_d1);
        Decomposition.visualise = true;
        Decomposition.dyn_dB = 50;
        [ImageFTheta, k, theta, handleFTheta, handleXY, handleWavelet] = decomposeFTheta(ImageTargets, Header, Decomposition);
        mkdir(sprintf('Decomposition_d1_d2/d_1 %d/', d1_vec(i_d1)))
        figure(handleFTheta);
        save2pdf(sprintf('Decomposition_d1_d2/d_1 %d/FTheta N_k %d N_theta %d', d1_vec(i_d1), N_k, N_theta)) 
        figure(handleXY);
        save2pdf(sprintf('Decomposition_d1_d2/d_1 %d/XY N_k %d N_theta %d', d1_vec(i_d1), N_k, N_theta)) 
        figure(handleWavelet);
        save2pdf(sprintf('Decomposition_d1_d2/d_1 %d/Wavelet N_k %d N_theta %d', d1_vec(i_d1), N_k, N_theta)) 
    end

%% ------------- END CODE ----------------
