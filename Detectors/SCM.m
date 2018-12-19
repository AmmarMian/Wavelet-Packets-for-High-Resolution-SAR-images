%% Description
% SCM estimator for covariance matrices
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
function R = SCM(x, Args)
% SCM: Compute tyler fixed point estimate of covariance matrix
% Syntax:  R = SCM(x, Args)
% Inputs:
%    x - Data for estimation
%    Args - Unused
%
% Outputs:
%    R - The estimate

    [p,N] = size(x);
    R = x*x'/N;

end
%% ------------- END CODE ----------------