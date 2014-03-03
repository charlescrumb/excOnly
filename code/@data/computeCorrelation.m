function [C, T] = computeCorrelation(this, varargin)
%[C, T] = computeCorrelation(this, varargin)
% COMPUTECORRELATION calls the Matlab native function to compute correlation.
% Also returns the execution time in seconds.
%
% USAGE:
%    mydata = data(F, P, N);
%    [C, T] = computeCorrelation(mydata);
%
%    optional arguments ('key' followed by its value): 
%    'debug', true/false
%    'discretize', true/false

%==========================================================================
% Package: ChaLearn Connectomics Challenge Sample Code
% Source: http://connectomics.chalearn.org
% Author: Isabelle Guyon
% Date: Dec 2013
% Last modified: NA
% Contact: causality@chalearn.org
% License: GPL v3 see http://www.gnu.org/licenses/
%==========================================================================

C=[];
if isempty(this.F), return; end 

%%% Assign default values
params.debug = false;
params.discretize = true;
params = parse_pv_pairs(params,varargin); 

if params.discretize, discretize(this); end 

fprintf('Computing correlation...\n');
tic; 
if params.discretize
    C= corrcoef(this.D); 
else
	C= corrcoef(this.F); 
end
% Suppress self connections
C(eye(size(C))==1)=0;
T=toc;
toc
fprintf('...done\n');

this.S.Correlation=C;
this.T.Correlation=T;

end
