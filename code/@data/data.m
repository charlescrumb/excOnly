%==========================================================================
% data structure             
%==========================================================================
% D = data(F, P, N) 
% D = data(F, P)
% D = data(F)
% D = data(other_data)
% 
% Creates a data structure holding Fluorescence data (F) a matrix (txn)
% of n times series of length t, corresponding to the activities of n
% neurons.
% P contained the coordinates of the neurons as a matrix (nx2).
% N contains the network connections N(i,j) is the connection strngth 
% of the connection neuron i -> neuron j. N(i, j)<=0 means no connection or
% disable connection.
%
% USAGE:
% mydata=data('data', 'mockvalid');
% score(mydata);
% roc(mydata); % or auc(mydata)
% write_score(mydata, 'results');
%
% Warning: a data object is a handle. To make a physical copy of D0, use
% D=data(D0); not D=D0; 

%==========================================================================
% Package: ChaLearn Connectomics Challenge Sample Code
% Source: http://connectomics.chalearn.org
% Author: Isabelle Guyon 
% Date: 01/29/2014
% Last modified: 01/30/2014
% Contact: causality@chalearn.org
% License: GPL v3 see http://www.gnu.org/licenses/
%==========================================================================

classdef data < handle
	properties (SetAccess = public)
        F=[];   % Fluorescence (matrix t x n)
        P=[];   % Positions (matrix n x 2)
        N=[];   % Network (matrix n x n)
        D=[];   % Discretized time series (matrix t x n)
        G=[];   % Global conditioning levels of discretization (matrix t x 1)
        S=[];   % Scores/predicted connection strengths for various methods.
        T=[];   % Time to compute the scores.
        A=[];   % AUC for computed scores.
        name='';% Base name of the network.
    end
    properties (SetAccess = private)
        rediscretize = false; % If true, rediscretize whenever the discretize method is called
                              % If false, discretize only the first time
                              % the method is called. In that case, if new
                              % settings are chosen, they are ignored and
                              % the first settings are used.
        rescore = false;      % If true, recomputes the scores when score is called.
        reauc = false;        % If true, recomputes the auc when auc is called.
                              % call clear_D, clear_S, clear_A to clear the
                              % discretized values, the scores, or the AUC.
    end
    methods
        %%%%%%%%%%%%%%%%%%%
        %%% CONSTRUCTOR %%%
        %%%%%%%%%%%%%%%%%%%
        function this = data(F, P, N) 
            %this = data
            %this = data(F, P, N)  
            %this = data(F, P) 
            %this = data(F)
            %this = data(other_data)
            %this = data(dirname, basename)
            % The last one loads data from file:
            % dirname/fluorescence_basename.txt
            % dirname/fluorescence_networkPositions.txt
            % dirname/fluorescence_network.txt
            if nargin<1, return; end
            if isa(F, 'data') 
                f=intersect(fields(obj), fields(this));
                for k=1:length(f)
                    this.(f{k})=obj.(f{k});
                end
            elseif ischar(F)
                this.load(F, P);
            else
                this.F=F;
                if nargin>1
                    this.P=P;
                end
                if nargin>2 
                    this.N=N;
                end
            end
        end     
        
        function F=get_F(this, num)
            %P=get_F(this, num)
            % Get the fluorescence values.
            % Get the neuron num only, if num is given.
            % If num is [] or not given, get all the values.
            if isempty(this.F), F=[]; return; end
            % Find the pattern number
            if nargin<2 || isempty(num)
                num=1:size(this.F, 2);
            else
                if length(num)==1 && (num<1 || num>size(this.F, 2)), F=[]; return; end
            end          
            F=this.F(:,num);
        end
        
        function P=get_P(this, num)
            %P=get_P(this, num)
            % Get the position values.
            % Get the neuron num only, if num is given.
            % If num is [] or not given, get all the values.
            if isempty(this.P), P=[]; return; end
            % Find the pattern number
            if nargin<2 || isempty(num)
                num=1:size(this.P, 1);
            else
                if length(num)==1 && (num<1 || num>size(this.P, 1)), P=[]; return; end
            end          
            P=this.P(num,:);
        end
        
        function net=network(this)
            %net=network(this)
            % get network as a structure with fields X, Y, and RS
            net.X=this.P(:,1);
            net.Y=this.P(:,2);
            net.RS=this.N;
        end
        
        function N=get_N(this)
            %N=get_N(this)
            % Get the network architecture.
            N=this.N;
        end
        
        function C=get_children(this, num)
            %C=get_children(this, num)
            % Get the children of neuron i (immediate descendents).
            C=this.N(num,:);
        end
        
        function C=get_parents(this, num)
            %C=get_parents(this, num)
            % Get the parents of neuron i (immediate ascendents).
            C=this.N(:,num);
        end              
        
        function [D, G]=get_D(this, num)
            %[D, G]=get_F(this, num)
            % Get the discretized fluorescence values and the corresponding
            % conditioning levels.
            % Get the neuron num only, if num is given.
            % If num is [] or not given, get all the values.
            if isempty(this.D), D=[]; G=[]; return; end
            % Find the neuron number
            if nargin<2 || isempty(num)
                num=1:size(this.D, 2);
            else
                if length(num)==1 && (num<1 || num>size(this.D, 2)), D=[]; G=[]; return; end
            end          
            D=this.D(:,num);
            G=this.G(num);
        end
        
        function clear_D(this)
            this.D=[];
            this.G=[];
        end
        
        function clear_S(this)
            this.S=[];
            this.T=[];
        end

        function clear_A(this)
            this.A=[];
        end                
        
        function n=length(this)
            %n=length(this)
            % number of neurons
            n=size(this.F, 2);
        end
        
        function n=sampnum(this)
            %n=sampnum(this)
            % number of samples
            n=size(this.F, 1);
        end
            
        function save_as_kaggle(this, filename)
            %save_as_kaggle(this, filename)
            % save the network architecture is Kaggle format
            writeNetworkScoresInCSV(filename,this.N,this.name);
        end
        
        function h=show_network(this, h)
            %h=show_network(this, h)
            % display the network connectivity as a matrix
            if nargin<2,
                h=figure('Name', 'Network connections');
            else 
                figure(h);
            end
            gmat_display(full(this.N),[],[],h);
        end
        
        function h=show_pair(this, neuron_i, neuron_j, range, h)
            %h=show_pair(this, neuron_i, neuron_j, range)
            % show the time series of a pair of neurons in range "range"
            if nargin<5,
                h=figure('Name', 'visualizePair: Compare original and discretized');
            else
                figure(h);
            end
            if nargin<2, neuron_i=1; end
            if nargin<3, neuron_j=2; end
            if nargin<4 || isempty(range) , range = 1:1000; end
            if isempty(this.F), return; end
            if ~isempty(this.D), subplot(2,1,1); end
            visualizePair(this.network, this.F, neuron_i, neuron_j, 'offset', 0.5,'samples',range);
            if ~isempty(this.D),
                subplot(2,1,2);
                visualizePair(this.network, this.D, neuron_i, neuron_j, 'offset', 0.5,'samples',range);
            end
        end
        
        function h=show_fluo(this, range, h)
            if nargin<3,
                h=figure('Name', 'visualizeFluorescenceTimeSeries: All series slightly offset');
            else
                figure(h);
            end
            if nargin<2, range = 1:1000; end
            if isempty(this.F), return; end        
            visualizeFluorescenceTimeSeries(this.F, 'type', 'single','offset',0.02, 'neuronList', 1:size(this.F,2), 'samples',range);
        end
        
        function nm=get_name(this)
            nm=this.name;
        end
        
        function load(this, dirname, basename)
            %load(this, dirname, basename)
            % Load data from files:
            % dirname/fluorescence_basename.txt
            % dirname/fluorescence_networkPositions.txt
            % dirname/fluorescence_network.txt
            this.name=basename;
            fluorescenceFile = [dirname filesep 'fluorescence_' basename  '.txt'];
            this.F = load_data(fluorescenceFile);
            networkFile = [dirname filesep 'network_' basename '.txt'];
            if exist(networkFile,'file')
                fprintf('Loading %s\n', networkFile);
                this.N = readNetworkScores(networkFile);
            end
            networkPositionsFile = [dirname filesep 'networkPositions_' basename '.txt'];
            fprintf('Loading %s\n', networkPositionsFile);
            this.P = load(networkPositionsFile);
        end
        
        function score(this, scorename, varargin)
            %score(this, scorename)
            %score(this, scorename, varargin)
            % Computes all the scores or just scorename, if 2nd arg given.
            % Stores them in the structure this.S.
            % Optional arguments ('key' followed by its value): 
            % 'debug', true/false
            % 'discretize', true/false.
            if nargin<2, 
                mymethods = methods(this);
                mymethods = mymethods(strmatch('compute', mymethods));
            else
                mymethods = {['compute' scorename]};
            end
            for k=1:length(mymethods)
                mname=mymethods{k}(8:end);
                if ~isfield(this.S, mname) || this.rescore 
                    if nargin>2
                        [this.S.(mname),this.T.(mname)]=feval(mymethods{k}, this, varargin{:});
                    else
                        [this.S.(mname),this.T.(mname)]=this.(mymethods{k});
                    end
                end
            end
        end
        
        function filenames=write_score(this, dirname, scorename)
            %filenames=write_score(this, dirname, scorename)
            %write_score(this, dirname)
            % Writes the score of name scorename or all scores if scorename is not
            % given, to a file [dirname]/[this.name]_[scorename].csv 
            % Also log the AUC, if we have the network connection truth
            % values in [dirname]/logfile.txt 
            if isempty(this.S), return; end
            if nargin<3 || isempty(scorename)
                scorenames=fields(this.S);
            else
                scorenames={scorename};
            end
            auc(this);
            logfile = [dirname filesep 'logfile.txt'];
            flog=fopen(logfile, 'a');
            the_date = datestr(now, 'yyyymmddTHHMMSS');
            filenames=cell(length(scorenames),1);
            for k=1:length(scorenames)
                filenames{k}=[dirname filesep scorenames{k} '_' this.name '_' the_date '.csv'];
                writeNetworkScoresInCSV(filenames{k},this.S.(scorenames{k}),this.name);
                fprintf(flog, '\n%s\t%s\t%s\t%5.4f\t%5.4f', ...
                    the_date, scorenames{k}, this.name, this.A.(scorenames{k}), this.T.(scorenames{k}));
            end
            fclose(flog);
        end
        
        function h=show_score(this, scorename, h)
            if nargin<2 || isempty(scorename)
                score=this.N;
                scorename = 'original network';
            else
                score=this.S.(scorename);
            end
            if nargin<3,
                h=figure('Name', ['Connection strength for: ' scorename]);
            else
                figure(h);
            end
            gmat_display(score,[],[],h);
            axis square
        end 
        
        function area=auc(this, scorename)
            %area=auc(this, scorename)
            %area=auc(this)
            % Compute the area under the ROC curve for the score scorename
            % or for all the scores, if scorename is not given.
            % The results are  stored in this.A
            % See also: roc.
            % Important: for efficiency reason, the AUC's already computed
            % are not recomputed.
            area=[];
            if isempty(this.S), return; end
            if nargin<2 || isempty(scorename)
                scorenames=fields(this.S);
            else
                scorenames={scorename};
            end
            area=[];
            for k=1:length(scorenames)
                if isempty(this.N)
                    this.A.(scorenames{k})=NaN;
                elseif ~isfield(this.A, scorenames{k}) || this.reauc
                    this.A.(scorenames{k})=computeAUC(this.N, this.S.(scorenames{k}));
                end
            end
            area=this.A;
        end
        
        function area=roc(this, scorename)
            %area=roc(this, scorename)
            %area=roc(this)
            % Compute the ROC curve for the score scorename
            % or for all the scores, if scorename is not given.
            % The AUC is returned and stored in this.A
            % See also: auc.
            area=[];
            if isempty(this.N), return; end
            if isempty(this.S), return; end
            if nargin<2 || isempty(scorename)
                scorenames=fields(this.S);
            else
                scorenames={scorename};
            end
            area=[];
            for k=1:length(scorenames)
                area.(scorenames{k})=computeROC(this.N, this.S.(scorenames{k}));
            end
            this.A=area;
        end
        
    end %methods
end % classdef
