%% Description
% Tyler fixed point estimator for covariance matrices
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
function R = Tyler(x, Args)
% Tyler: Compute Tyler fixed point estimate of covariance matrix
% Syntax:  R = Tyler(x, Args)
% Inputs:
%    x - Data for estimation
%    Args - [Tol, Itermax] with
%       Tol = tolerance for Tyler estimation,
%       Itermax = Number of maximum iterations for Tyler estimation.
% Outputs:
%    R - The estimate

    % Variables
    tol = Args(1);
    nitermax = Args(2);
    [p,N] = size(x);
    err = inf;
    % Init
    R = eye(p);
    niter = 0;
    % Iterations
    while err>tol && niter<nitermax
        % Compute next iteration
        v = chol(R)' \ x;
        a = mean(v.*conj(v));
        y = x./ sqrt( a(ones(p,1),:) );
        Rnew =  y*y' / N;
        Rnew = p * Rnew / sum(diag(Rnew)); % Normalize by the trace
        
        % Criterion
        err = norm(Rnew-R,'fro')/norm(R,'fro');
        R=Rnew;
        niter = niter+1;
    end

    if niter == nitermax 
        disp('Tyler: did not converge')
    end

end
%% ------------- END CODE ----------------