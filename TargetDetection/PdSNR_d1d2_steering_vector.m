%% Description
% Simulation of Pd-SNR with variation of d1 and d2
%% Specifiactions
% Other m-files required: target creation functions, detectors,
% decomposeFTheta
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
    addpath('..', '../Hyperimages', '../Target Dictionary', '../Detectors')

%% Simulations parameters

    % M-C Trials
    numberOfTrials = 1000;
    numberOfTrials_steering_vector = 100;
    SNR = 0;
    Pfa = 0.001;
    targetGuard = 70; % Guard for position of target on the image
    
    % Decomposition
    Decomposition = struct;
    Decomposition.method = 'Exponential';
    Decomposition.J = 1;
    Decomposition.numberOfBands = 5;
    Decomposition.numberOfLooks = 5;
    Decomposition.visualise = false;
    Decomposition.dyn_dB = 70;
    N_k = Decomposition.numberOfBands^Decomposition.J;
    N_theta = Decomposition.numberOfLooks^Decomposition.J;
    m = N_k*N_theta;
    d1_vec = [10, inf];
    d2_vec = d1_vec;
    
    % Detectors
    DetectorsStruct = initTargetDetectors();
    [My, Mx] = size(DetectorsStruct.mask); % Mask
    
%% Reading Image
 
    imagePath = '../../../../../Data/SANDIA/MiniSAR20050519p0010image002.mat';
    [Image, Header] = readSandia(imagePath, true, [1, 1, 501, 501]);
%     [Image, Header] = readSDMS('C:\Users\mian_amm\Dropbox\Thèse\Data\SDMS', 'FP0121', true, [800, 800, 1201, 1201]);
%     Image = Image(:,:,3);
    [Ny, Nx] = size(Image);
    
%% Pfa-lambda to select threshold
    Decomposition.d1 = inf;
    Decomposition.d2 = inf;
    [ImageFTheta, k, theta] = decomposeFTheta(Image, Header, Decomposition);  
    ImageFTheta = ImageFTheta(1:N_theta:end, 1:N_k:end,:);
    p = randn(N_k*N_theta,1) + 1i*randn(N_k*N_theta,1);
    p = p/(p'*p);
    Results = computeTestOnImage(ImageFTheta, p, DetectorsStruct);
    thresholds = zeros(DetectorsStruct.number,1);
    for d=1:DetectorsStruct.number
        Nsimu = length(Results{d}(:));
        Pfa_vector = (Nsimu:-1:1)/ Nsimu;
        lambda = sort(Results{d}(:));
        thresh = lambda(abs(Pfa - Pfa_vector) == min(abs(Pfa - Pfa_vector)));
        thresholds(d) = thresh(1);
    end
    
%% Pd-SNR
    Pd_comparison_mean = zeros(DetectorsStruct.number, length(d1_vec), numberOfTrials_steering_vector);
    h = waitbar(0,'Please wait for Pd-SNR...');
    for i_steering_vector=1:numberOfTrials_steering_vector
        p = randn(N_theta,N_k) + 1i*randn(N_theta,N_k);
        for i_d1 = 1:length(d1_vec)
            Decomposition.d1 = d1_vec(i_d1);
            Decomposition.d2 = d2_vec(i_d1);

            Pd_temp = cell(numberOfTrials,1);

            parfor trial=1:numberOfTrials % M-C Trials

                svTarget = struct;
                svTarget.position = [N_k * randi([floor(targetGuard/N_k), ...
                            floor((Ny - targetGuard)/N_k)],1) + 1, ...
                            N_theta * randi([floor(targetGuard/N_theta), ...
                            floor((Nx - targetGuard)/N_theta)],1) + 1];
                svTarget.Args = {p, N_k, N_theta};
                svTarget.SNR_dB = SNR;
                svTarget.function = @steeringVectorTarget;

                % Creation of Target
                TargetStruct = struct;
                TargetStruct.list = [svTarget];
                TargetStruct.number = length(TargetStruct.list);
                TargetStruct.Mx = 20; % For Noise Estimation
                TargetStruct.My = 20;
                ImageTargets = createTargets(Image, Header, TargetStruct);

                % Decomposition in (f, theta)
                [ImageFTheta,~,~] = decomposeFTheta(ImageTargets + Image, Header, Decomposition);
                ImageFTheta = ImageFTheta(1:N_theta:end, 1:N_k:end,:);

                % Detection test on target pixel
                NI = find(reshape(DetectorsStruct.mask, Mx*My,1) == 1);
                ii = (svTarget.position(2) - 1) / N_k + 1 - floor(My/2);
                jj = (svTarget.position(1) - 1) / N_theta + 1 - floor(Mx/2);
                Data = struct;
                Data.X = ImageFTheta(ii:ii+(My-1), jj:jj+(Mx-1), :);
                Data.X = reshape(Data.X, Mx*My, m).';
                Data.X = Data.X(:, NI);
                Data.y = reshape(ImageFTheta(ii+floor((My-1)/2), jj+floor((Mx-1)/2), :),m,1);
                Data.p = p.';
                Data.p = Data.p(:);
                Data.p = Data.p/sqrt(Data.p'*Data.p);
                for d=1:DetectorsStruct.number
                    detector = DetectorsStruct.list(d);
                    Pd_temp{trial}(d) = detector.detect(Data, detector.Estimation) >= thresholds(d);
                end

            end
            Pd_temp = cell2mat(Pd_temp);
            % Averaging trials
            for d=1:DetectorsStruct.number
               Pd_comparison_mean(d, i_d1, i_steering_vector) = mean(Pd_temp(:, d)); 
            end
        end

       percent = i_steering_vector / numberOfTrials_steering_vector;
       waitbar(percent, h);  
    end
    close(h)

%% ------------- END CODE ----------------
