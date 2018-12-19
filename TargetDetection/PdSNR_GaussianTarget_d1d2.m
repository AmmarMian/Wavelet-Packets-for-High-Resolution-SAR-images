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
    SNR_vec = [0:2:50];
    Pfa = 0.001;
    targetGuard = 70; % Guard for position of target on the image
    
    % Decomposition
    Decomposition = struct;
    Decomposition.method = 'Exponential';
    Decomposition.J = 1;
    Decomposition.numberOfBands = 7;
    Decomposition.numberOfLooks = 7;
    Decomposition.visualise = false;
    N_k = Decomposition.numberOfBands^Decomposition.J;
    N_theta = Decomposition.numberOfLooks^Decomposition.J;
    m = N_k*N_theta;
    d1_vec = [10, 100, inf];
    d2_vec = d1_vec;
    
    % Detectors
    DetectorsStruct = initTargetDetectors();
    [My, Mx] = size(DetectorsStruct.mask); % Mask
    
%% Reading Image

    imagePath = '../../../../../Data/SANDIA/MiniSAR20050519p0010image002.mat';
    [Image, Header] = readSandia(imagePath, true, [1, 1, 700, 700]);
%     [Image, Header] = readSDMS('C:\Users\mian_amm\Dropbox\Thèse\Data\SDMS', 'FP0121', false, [100, 500, 1100, 1500]);
%     Image = Image(:,:,3);
    [Ny, Nx] = size(Image);
    
%% Pfa-lambda to select threshold
%     Decomposition.d1 = inf;
%     Decomposition.d2 = inf;
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
    thresholds = [2.202285750664898e+02;0.319801970816285;4.761946266981409e+06]; % For Pfa=10^-3
    
%% Pd-SNR
    Pd_comparison = cell(DetectorsStruct.number, length(d1_vec));
    p = randn(N_theta,N_k) + 1i*randn(N_theta,N_k);
    for i_d1 = 1:length(d1_vec)
        Decomposition.d1 = d1_vec(i_d1);
        Decomposition.d2 = d2_vec(i_d1);

        Pd_vec = zeros(length(SNR_vec),DetectorsStruct.number);
        h = waitbar(0,'Please wait for Pd-SNR...');
        for i_SNR = 1:length(SNR_vec) % Loop on SNR

           SNR = SNR_vec(i_SNR);
           Pd_temp = cell(numberOfTrials,1);

           parfor trial=1:numberOfTrials % M-C Trials

                % Target description
                svTarget = struct;
                svTarget.position = [N_k * randi([floor(targetGuard/N_k), ...
                            floor((Ny - targetGuard)/N_k)],1) + 1, ...
                            N_theta * randi([floor(targetGuard/N_theta), ...
                            floor((Nx - targetGuard)/N_theta)],1) + 1];
                meanG = [112, 0];
                scaleG = [1, 0.01];
                svTarget.Args = {meanG, scaleG};
                svTarget.SNR_dB = SNR;
                svTarget.function = @gaussianTarget;
                TargetStruct = struct;
                TargetStruct.list = [svTarget];
                TargetStruct.number = length(TargetStruct.list);
                TargetStruct.Mx = 20; % For Noise Estimation
                TargetStruct.My = 20;
                ImageTargets = createTargets(Image, Header, TargetStruct);
                p = spectreToSteeringVector(fft2(ImageTargets), Header.kx, Header.ky, Decomposition);

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
               Pd_vec(i_SNR, d) = mean(Pd_temp(:, d)); 
           end

           percent = i_SNR / length(SNR_vec);
           waitbar(percent, h);  

        end
        close(h);

        % Plotting
        figure('Position', [100, 100, 700, 300]);
        markerstyle = {'d', 's', '*', '+', '.', '<', '>', 'v', '^', 'pentagram', 'o'};
        lg = {};
        for d=1:DetectorsStruct.number
            Pd_comparison{d, i_d1} = Pd_vec(:, d);
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
        textString = {sprintf('$d_1=%.2f$, $d_2=%.2f$', Decomposition.d1, Decomposition.d2) ...
            ,sprintf('$j=%d$, $R=%d$, $L=%d$, $K=%d$', Decomposition.J, Decomposition.numberOfBands, Decomposition.numberOfLooks, sum(DetectorsStruct.mask(:)))...,
            sprintf('$P_{FA}=%.3f$, %d Trials', Pfa, numberOfTrials)};
        annotation('textbox',...
        [0.15, 0.57, 0.1, 0.1], ...
        'Units','characters',...
        'String', textString,...
        'backgroundcolor', 'w', ...
        'FitBoxToText','on',...
        'interpreter', 'latex');
%         savefig(sprintf('Results/PdSNR %d_trials Pfa_%d j_%d m_%d n_%d d1_%.2f d2_%.2f K_%d.fig', ...
%                         numberOfTrials, Pfa, Decomposition.J, Decomposition.numberOfBands, ...
%                         Decomposition.numberOfLooks, Decomposition.d1, Decomposition.d2, ...
%                         sum(DetectorsStruct.mask(:))));
    end

%% Comparison
linestyles = {'--', ':', '-','-.'};
colors = {'y', 'm', 'c' 'r', 'g', 'b', 'k'};
for d=1:DetectorsStruct.number
    figure('Position', [100, 100, 700, 400])
    lg = {};
    for i_d1 = 1:length(d1_vec)
        lg = cat(2, lg, sprintf('%s $d_1=%d$, $d_2=%d$', DetectorsStruct.list(d).name, d1_vec(i_d1), d2_vec(i_d1)));
        Pd =  Pd_comparison{d, i_d1};
        plot(SNR_vec, Pd, 'k', 'linestyle', '--', 'marker', markerstyle{i_d1}, 'markersize', 10)
        hold on
    end
    h = legend(lg);
    set(h, 'interpreter', 'latex', 'location', 'southeast')
    grid
    grid minor
    set(h, 'interpreter', 'latex', 'location', 'northwest')
    xlabel('SNR (dB)', 'interpreter', 'latex', 'fontsize', 18)
    ylabel('$P_D$', 'interpreter', 'latex', 'fontsize', 18)
    textString = {sprintf('$d_1=%.2f$, $d_2=%.2f$', Decomposition.d1, Decomposition.d2) ...
        ,sprintf('$j=%d$, $R=%d$, $L=%d$, $K=%d$', Decomposition.J, Decomposition.numberOfBands, Decomposition.numberOfLooks, sum(DetectorsStruct.mask(:)))...,
        sprintf('$P_{FA}=%.3f$, %d Trials', Pfa, numberOfTrials)};
    annotation('textbox',...
    [0.15, 0.6, 0.1, 0.1], ...
    'Units','characters',...
    'String', textString,...
    'backgroundcolor', 'w', ...
    'FitBoxToText','on',...
    'interpreter', 'latex');
    savefig(sprintf('Results/Gaussian/PdSNR Comparison detector %d %d_trials Pfa_%d j_%d m_%d n_%d d1_%.2f d2_%.2f K_%d.fig', ...
                d, numberOfTrials, Pfa, Decomposition.J, Decomposition.numberOfBands, ...
                Decomposition.numberOfLooks, Decomposition.d1, Decomposition.d2, ...
                sum(DetectorsStruct.mask(:))));
%     axes('position',[.65 .175 .25 .25])
%     for i_d1 = 1:length(d1_vec)
%         Pd =  Pd_comparison{d, i_d1};
%         plot(SNR_vec, Pd, 'k', 'linestyle', '--', 'marker', markerstyle{i_d1}, 'markersize', 10)
%         hold on
%         grid
%         grid minor
% %         ylim([0.16, 0.24])
% %         xlim([4, 6])
%         xlabel('SNR (dB)', 'interpreter', 'latex')
%         ylabel('$P_D$', 'interpreter', 'latex')
%     end
end

%% ------------- END CODE ----------------
