%% Description
% Initialisation of Target structures
%% Specifiactions
% Other m-files required: Detectors function
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

    gaussTarget = struct;
    gaussTarget.position = [50, 30];
    mean = [110, 1.555];
    scale = [0.5, 0.001];
    gaussTarget.Args = {mean, scale};
    gaussTarget.SNR_dB = 70;
    gaussTarget.function = @gaussianTarget;
    
    rectTarget = struct;
    rectTarget.position = [20, 70];
    mean = [112, 1.6];
    scale = [5, 0.01];
    rectTarget.Args = {mean, scale};
    rectTarget.SNR_dB = 30;
    rectTarget.function = @rectangularTarget;
    
    TargetStruct = struct;
    TargetStruct.list = [rectTarget, gaussTarget];
    TargetStruct.number = length(TargetStruct.list);
    TargetStruct.Mx = 20;
    TargetStruct.My = 20;

%% ------------- END CODE ----------------
