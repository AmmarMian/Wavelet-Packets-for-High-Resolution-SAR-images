%% Description
% Initialisation of detectors structures for target detection
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
function DetectorsStruct = initTargetDetectors()

    ANMFTylerStruct = struct;
    ANMFTylerStruct.name = '$\hat{\Lambda}_{\mathrm{ANMF-Tyler}}$';
    ANMFTylerStruct.detect = @ANMFFunction;
    ANMFTylerStruct.Estimation = struct;
    ANMFTylerStruct.Estimation.Args = [0.0001, 50];
    ANMFTylerStruct.Estimation.MxEstimator = @Tyler;
    ANMFTylerStruct.threshold = 0.3566;
    ANMFTylerStruct.theo = @(lambda, N, m) exp(((m/(m+1)*N)-m+1)*log(1-lambda)+log(mhygfx(lambda,(m/(m+1)*N)-m+2,(m/(m+1)*N)-m+1,(m/(m+1)*N)+1)));
    
    
    MahalanobisStruct = struct;
    MahalanobisStruct.name = '$\hat{\Lambda}_{\mathrm{Mahalanobis}}$';
    MahalanobisStruct.detect = @MahalanobisDistance;
    MahalanobisStruct.Estimation = struct;
    MahalanobisStruct.Estimation.Args = [0.0001, 50];
    MahalanobisStruct.Estimation.MxEstimator = @Tyler;
    MahalanobisStruct.threshold = 0;
    MahalanobisStruct.theo = 0;
    
    ANMFSCMStruct = struct;
    ANMFSCMStruct.name = '$\hat{\Lambda}_{\mathrm{ANMF-SCM}}$';
    ANMFSCMStruct.detect = @ANMFFunction;
    ANMFSCMStruct.Estimation = struct;
    ANMFSCMStruct.Estimation.Args = 0;
    ANMFSCMStruct.Estimation.MxEstimator = @SCM;
    ANMFSCMStruct.threshold = 0.3566;
    ANMFSCMStruct.theo = @(lambda, N, m) exp((N-m+1)*log(1-lambda)+log(mhygfx(lambda,N-m+2,N-m+1,N+1)));
    
    AMFStruct = struct;
    AMFStruct.name = '$\hat{\Lambda}_{\mathrm{AMF}}$';
    AMFStruct.detect = @AMFFunction;
    AMFStruct.Estimation = struct;
    AMFStruct.Estimation.Args = 0;
    AMFStruct.Estimation.EstimateMx = true;
    AMFStruct.Estimation.MxEstimator = @SCM;
    AMFStruct.threshold = 487.4567;
    AMFStruct.theo = @(lambda, N, m) mhygfx(-lambda/N, N-m+1, N-m+2, N+1);
    
    DetectorsStruct = struct;
    DetectorsStruct.list = [AMFStruct, ANMFTylerStruct];
    DetectorsStruct.number = length(DetectorsStruct.list);
    DetectorsStruct.mask = ones(13);
    [My, Mx] = size(DetectorsStruct.mask);
    guard = 4;
    DetectorsStruct.mask(floor(My/2)+1-guard:floor(My/2)+1+guard,floor(Mx/2)+1-guard:floor(Mx/2)+1+guard) = 0;
    

end
%% ------------- END CODE ----------------
