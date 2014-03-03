function [Y, map] = gmat_display(X, num, noshow, h)
%[Y, map] = gmat_display(X, num, noshow, h)
% Display matrix X in gray levels after rescaling.
% Inputs:
% X -- matrix to be displayed (real neumbers).
% Optional:
% num -- number of gray levels (default 256).
% noshow -- a flag. If 1, the figure is not displayed.
% Returns:
% Y -- rescaled matrix X (in uint8, to make an image).

%==========================================================================
% Package: ChaLearn Connectomics Challenge Sample Code
% Source: http://connectomics.chalearn.org
% Author: Isabelle Guyon
% Date: April 2002
% Last modified: 01/30/2014
% Contact: causality@chalearn.org
% License: GPL v3 see http://www.gnu.org/licenses/
%==========================================================================

if (nargin<2||isempty(num)), num=256; end
if (nargin<3||isempty(noshow)), noshow=0; end
if (nargin<4||isempty(h)), h=figure; end

map=gray(num);
if isempty(X), Y=[]; return, end;

Y = inormalize(full(X));
Y = uint8(Y*(num-1));
if(~noshow) 
    figure(h);
    colormap(map);
    image(Y); 
end



