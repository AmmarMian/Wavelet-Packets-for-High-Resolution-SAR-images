%% Description
% Initialisation of detectors structures
%% Specifiactions
% Other m-files required: none
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
function Results = computeTestOnImage(Image, p, DetectorsStruct)
% ImageTargets - Create several targets
% Syntax:  Results = computeTestOnImage(Image, p, DetectorsStruct)
%
% Inputs:
%    - Image = the image on which we add the targets
%    - p = steering vector
%    - DetectorsStruct - structure for detectors
%
% Outputs:
%    Results - results of detection test in a cell
%

    mask = DetectorsStruct.mask;
    [Ny, Nx, m] = size(Image);
    [My, Mx] = size(mask);
    II = 1 : Ny - My + 1;
    JJ = 1 : Nx - Mx + 1;
    NI = find(reshape(mask, Mx*My,1) == 1);
    
    % Results container
    Results = cell(DetectorsStruct.number,1);
    for d=1:DetectorsStruct.number
        Results{d} = zeros(Ny - My + 1, Nx - Mx + 1);
    end
    Data = struct;
    Data.p = p;
    
    h = waitbar(0,'Please wait...');
        for ii = II
            waitbar((ii-II(1))/length(II), h);
            for jj = JJ

                Data.X = Image(ii:ii+(My-1), jj:jj+(Mx-1), :);
                Data.X = reshape(Data.X, Mx*My, m).';
                Data.X = Data.X(:, NI);
                Data.y = reshape(Image(ii+floor((My-1)/2), jj+floor((Mx-1)/2), :),m,1);

                for d=1:DetectorsStruct.number
                    detector = DetectorsStruct.list(d);
                    Results{d}(ii,jj) = detector.detect(Data, detector.Estimation);
                end
                              
            end
        end
    close(h)
end
%% ------------- END CODE ----------------
