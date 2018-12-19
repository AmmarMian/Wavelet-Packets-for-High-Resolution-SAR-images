%% Description
% Simulation of Pd-SNR
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
    SNR_vec = [-15:2:25];
    Pfa = 0.001;
    targetGuard = 50; % Guard for position of target on the image
    
    % Decomposition
    Decomposition = struct;
    Decomposition.method = 'Exponential';
    Decomposition.J = 1;
    Decomposition.numberOfBands = 5;
    Decomposition.numberOfLooks = 5;
    N_k = Decomposition.numberOfBands^Decomposition.J;
    N_theta = Decomposition.numberOfLooks^Decomposition.J;
    Decomposition.d1 = inf;
    Decomposition.d2 = inf;
    m = N_k*N_theta;
    
    % Detectors
    DetectorsStruct = initTargetDetectors();
    [My, Mx] = size(DetectorsStruct.mask); % Mask
    
%% Reading Image

    imagePath = '../../../../Data/SANDIA/MiniSAR20050519p0010image002.mat';
    [Image, Header] = readSandia(imagePath, true, [1, 1, 200, 200]);
%     [Image, Header] = readSDMS('C:\Users\mian_amm\Dropbox\Thèse\Data\SDMS', 'FP0121', false, [100, 500, 1100, 1500]);
%     Image = Image(:,:,3);
    [Ny, Nx] = size(Image);
    
%% Pfa-lambda to select threshold

%     [ImageFTheta, k, theta] = decomposeFTheta(Image, Header, Decomposition);  
%     ImageFTheta = ImageFTheta(1:N_theta:end, 1:N_k:end,:);
%     p = randn(N_k*N_theta,1) + 1i*randn(N_k*N_theta,1);
%     p = p/(p'*p);
%     Results = computeTestOnImage(ImageFTheta, p, DetectorsStruct);
%     thresholds = zeros(DetectorsStruct.number,1);
%     for d=1:DetectorsStruct.number
%         Nsimu = length(Results{d}(:));
%         Pfa_vector = (Nsimu:-1:1)/ Nsimu;
%         lambda = sort(Results{d}(:));
%         thresh = lambda(abs(Pfa - Pfa_vector) == min(abs(Pfa - Pfa_vector)));
%         thresholds(d) = thresh(1);
%     end
    thresholds = [1.358426864904436e+02;0.342878292230499]; % For Pfa=10^-3
    
%% Pd-SNR
    
    Pd_vec = zeros(length(SNR_vec),DetectorsStruct.number);
    h = waitbar(0,'Please wait for Pd-SNR...');
    for i_SNR = 1:length(SNR_vec) % Loop on SNR

       SNR = SNR_vec(i_SNR);
       Pd_temp = zeros(numberOfTrials,DetectorsStruct.number);
       
       for trial=1:numberOfTrials
           
            % Target description: p =random, posiiton = random
            p = randn(N_theta,N_k) + 1i*randn(N_theta,N_k);
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
            Data.p = p(:)/(p(:)'*p(:));
            for d=1:DetectorsStruct.number
                detector = DetectorsStruct.list(d);
                Pd_temp(trial, d) = detector.detect(Data, detector.Estimation) >= thresholds(d);
            end

       end
       
       % Averaging trials
        for d=1:DetectorsStruct.number
            Pd_vec(i_SNR, d) = mean(Pd_temp(:, d)); 
        end
        
       percent = i_SNR / length(SNR_vec);
       waitbar(percent, h);  
       
    end
    close(h);
%% Plotting
    
    figure('Position', [100, 100, 700, 300]);
    markerstyle = {'d', 's', '*', '+'};
    lg = {};
    for d=1:DetectorsStruct.number
        plot(SNR_vec, Pd_vec(:, d), 'k', 'linestyle', '--', 'marker', markerstyle{d}, 'markersize', 10)
        lg = cat(2,lg, DetectorsStruct.list(d).name);
        hold on
    end
    h = legend(lg);
    set(h, 'interpreter', 'latex', 'fontsize', 15, 'location', 'northwest')
    xlabel('SNR (dB)', 'interpreter', 'latex', 'fontsize', 18)
    ylabel('$P_D$', 'interpreter', 'latex', 'fontsize', 18)
    grid
    grid minor
    textString = {sprintf('$d_1=%d$, $d_2=%d$', Decomposition.d1, Decomposition.d2) ...
        ,sprintf('$j=%d$, $m=%d$, $n=%d$, $K=%d$', Decomposition.J, Decomposition.numberOfBands, Decomposition.numberOfLooks, sum(DetectorsStruct.mask(:)))...,
        sprintf('$P_{FA}=%.3f$, %d Trials', Pfa, numberOfTrials)};
    annotation('textbox',...
    [0.15, 0.57, 0.1, 0.1], ...
    'Units','characters',...
    'String', textString,...
    'backgroundcolor', 'w', ...
    'FitBoxToText','on',...
    'interpreter', 'latex');
    savefig(sprintf('Results/PdSNR %d_trials Pfa_%d j_%d m_%d n_%d d1_%d d2_%d K_%d.fig', ...
                    numberOfTrials, Pfa, Decomposition.J, Decomposition.numberOfBands, ...
                    Decomposition.numberOfLooks, Decomposition.d1, Decomposition.d2, ...
                    sum(DetectorsStruct.mask(:))));
        
%% ------------- END CODE ----------------
