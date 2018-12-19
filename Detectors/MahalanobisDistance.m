%% Description
% Mahalanobis distance
%% Specifiactions
% Other m-files required: Tyler.m, SCM.m
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
function ratio = MahalanobisDistance(Data, Estimation)
% ANMFFunction: Compute ANMF detection Test
% Syntax:  ratio = ANMFFunction(Data, Estimates)
% Inputs:
%    - Data = A struct with fields:
%       * y = test pixel
%       * X = secondary obersertions.
%       * m - size of vectors
%       * N - Number Of secodnary observations.
%    - Estimation = A struct with fields:
%       * MxEstimator = pointer to the function or estimation of Mx
%       * ArgsEstimator = if Tyler estimator: [Tol, Itermax] with
%       Tol = tolerance for Tyler estimation,
%       Itermax = Number of maximum iterations for Tyler estimation.
% Outputs:
%    ratio - the test ratio
    

    Mx = Estimation.MxEstimator(Data.X, Estimation.Args);
    iMx = inv(Mx);
    y = Data.y;
    ratio = abs(y'*iMx*y);

end
%% ------------- END CODE ----------------