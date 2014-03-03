function browse(this, h)
% browse(this, h)
% Graphical display of fluorescence values of pairs

%==========================================================================
% Package: ChaLearn Connectomics Challenge Sample Code
% Source: http://connectomics.chalearn.org
% Authors: Javier Orlandi, Mehreen Saeed, Isabelle Guyon
% Date: Jan 2014
% Last modified: 01/29/2014
% Contact: causality@chalearn.org
% License: GPL v3 see http://www.gnu.org/licenses/
%==========================================================================

if nargin<2 || isempty(h), h=figure('Position', [0 245 560 560]); else figure(h); end

N=length(this);
n=-2;
gg=-3;
g=-4;
m=-5;
b=-6;
bb=-7;
p=-1;
e=0;
n1=1;
n2=2;
while 1
    cla;
    show_pair(this, n1, n2, [], h);
    dir=' ? ';
    if ~isempty(this.N)
        if this.N(n1, n2) >0
            if this.N(n2, n1) >0
                dir=' <-> ';
            else
                dir=' -> ';
            end
        else
            if this.N(n2, n1) >0
                dir=' <- ';
            else
                dir=' | ';
            end
        end
    end
    title(['neuron' num2str(n1) dir 'neuron' num2str(n2) ' (of ' num2str(N) ' neurons)'], 'FontSize', 16);


    ans = input('Neuron numbers [n1 n2] (or n for next, p for previous, e exit)? ');
    idx=ans;
    if isempty(idx), idx=n; end
    if length(idx)>1,
        idx=Inf;
    end
    switch idx
    case {n, gg, g, m, b, bb}
        n2=n2+1;
        if n2>N, n1=n1+1; n2=1; end
        if n1>N, n1=1; end
    case p
        n2=n2-1;
        if n2<1, n2=N; n1=n1-1; end
        if n1<1, n1=N; end
    case e	
        break
    otherwise
        n1=ans(1);
        n2=ans(2);
    end
    
end


