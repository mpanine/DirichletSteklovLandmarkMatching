function [SpQR,QROptions] = gsqr(A, Options)

% function SpQR = gsqr(A, [Options])
%
% GSQR computes a partial QR decomposition of the matrix A via modified
% Gram-Schmidt orthogonalization.  A rank revealing QR decomposition is a
% factorization
%
%      A*P = | Q11 Q12 | * | R11 R12 |
%                          |  0  R22 |
%
% where P is a permutation, Q = [Q11 Q12] is orthogonal, and R22 is a matrix
% with small L^2 norm.  Under these assumptions, the contribution to the
% matrix product from the term Q12*R22 can be safely ignored (introducing
% a controllable error).
%
% GSQR computes the matrices Q11 and R11 in the RRQR decomposition, which
% yields the approximate factorization
%
%      [A(:,SpQR.Idxs) A(:,SpQR.DIdxs)] = [Q11*R11 Q11*R12]
%
% where R12 = Q11'*A(:,SpQR.DIdxs).
%
%
% In:
%   A        = Input matrix
%   Options  = A structure specifying algorithm settings.
%
%       QRroutine:  any of {'gsqr','suitesparse'}. Default: 'gsqr'. Note that
%                   'suitesparse' does not necessarily give a Rank-Revealing Factorization.
%
%   (1) Stopping conditions
%
%     StopPrecision   : stop when the l^2 norm of all remaining columns falls below this value
%     StopN           : specifies the maximum number of iterations
%     StopDensity     : stop when the density of all remaining columns becomes
%                       larger than this. Exception: no column could be picked, in which case returns one column.
%     StopRelJump     : stop at the last time that the diagonal entry of R_11 jumps down by a relative 
%                       amount greater than this, w.r.t. the previous diagional entry
%
%   (2) Thresholding
%
%     QThreshold       : threshold for the entries of the matrix Q. 
%                        For QRroutine=='gsqr', it is applied to each new column of Q computed.
%                        For QRroutine=='suitesparse', it is applied to the final Q.
%
%   (3) Other Options
%
%     Quiet           : when true, all output other than error messages is suppressed. Default: 0.
%     PivotingLayers  : labels for pivoting layers. Not used if QRroutine=='suitesparse'. Default: only one layer.
%     NesDis          : use nested dissection to create the pivoting layers. Not used if QRroutine=='suitesparse'.
%     ClusterMinSize  : size of small cluster for nesdis. Not used if QRroutine=='suitesparse'.
%     MinFcnPerCluster: minimum number of functions to be picked per cluser. Not used if QRroutine=='suitesparse'. Default: 0.
%
%   (4) Diagnostics:
%
%     RunDiagnostics  : runs diagnostics on numerical accuracy of procedure. Typically very expensive ( O(N^3) ).
%
%
%   The default options are as follows:
%
%     Options.StopPrecision      = eps;
%     Options.StopN              = Inf;
%     Options.StopDensity        = 1.0;
%     Options.StopRelJump        = 0.0;
%     Options.QThreshold         = 0.0;
%     Options.Quiet              = false;
%     Options.Symmetric          = false;
%
% OUT: 
%   SpQR structure containing the following fields:
%       Idxs     = List of chosen columns, in a rank-revealing fashion
%       DIdxs    = List of discarded columns, in no specific order
%       Q        = The matrix Q11 in the RRQR decomposition
%       R11      = The matrix R11 in the RRQR decomposition
%       R12      = The matrix R12 in the RRQR decomposition
%
% Dependences:
%   none
%
% Version History:
%   jcb        2/2006         cell array version completed; elimanted RI
%                             computation (for now)
%   mm         4/2007         added parameters, fixed a few minor bugs
%   mm         12/2008
%
%


SpQR.Idxs   = [];
SpQR.DIdxs  = [];
SpQR.Q      = [];
SpQR.R11    = [];
SpQR.R12    = [];
SpQR.NesDis = [];

Idxs = []; DIdxs = []; Q=[]; R=[];

if isempty(A),
    return;
end

% Setup some basic options
[M N] = size(A);

if ~exist('Options')    Options = []; end;

if ~isfield(Options, 'QRroutine')       Options.QRroutine = 'gsqr';         end;
if ~isfield(Options, 'StopN')           Options.StopN = Inf;                end;
if ~isfield(Options, 'StopPrecision')   Options.StopPrecision = eps;        end;
if ~isfield(Options, 'StopDensity')     Options.StopDensity = 1.0;          end;
if ~isfield(Options, 'StopRelJump')     Options.StopRelJump = 0.0;          end;
if ~isfield(Options, 'QThreshold')      Options.QThreshold = 0.0;           end;
if ~isfield(Options, 'Quiet')           Options.Quiet = false;              end;
if ~isfield(Options, 'PivotingLayers')  Options.PivotingLayers = ones(1,N); end;
if ~isfield(Options, 'NesDis')          Options.NesDis = false;             end;
if ~isfield(Options, 'ClusterMinSize')  Options.ClusterMinSize = [];        end;
if ~isfield(Options, 'MinFcnPerCluster') Options.MinFcnPerCluster = 0;      end;
if ~isfield(Options, 'NumProjections')  Options.NumProjections = 2;         end;


T = cputime;

% Prepare display
if ~Options.Quiet
    if ~isfield(Options, 'dwtree')
        fprintf('gsqr.m: ');
    end;
    fprintf('00000');
end;



switch lower(Options.QRroutine),
    case 'suitesparse'
        % Compute the QR decomposition. A(:,P)=QR
        [Q,R,P] = spqr(A,struct('ordering','metis','permutation','vector','tol',Options.StopPrecision,'econ',0));

        NumChosen = size(Q,2);
        
        % Threshold Q if required
        if Options.QThreshold>0.0,
            Q = Q.*(abs(Q)>Options.QThreshold);
        end;        
        
        % Sort the QR decomposition based on the size of the diagonal of R
        lDiagR = abs(diag(R));
        [lTmp,lSortIdxs] = sort(full(abs(lDiagR)),1,'descend');

        % Truncate if stopping based on relative jump in the diagonal of R is requested
        if Options.StopRelJump>0,       %% THIS WON'T WORK, TODO
            lLastBigJump = NumChosen;
            for j = 2:NumChosen,
                if lDiagR(lSortIdxs(j)) < Options.StopRelJump * lDiagR(lSortIdxs(j-1)),
                    lLastBigJump = j-1;
                end;
            end;
            if lLastBigJump < NumChosen,
                Chosen(Idxs(lLastBigJump+1:NumChosen)) = false;
                NumChosen = lLastBigJump;
            end;
        end;

        % Sort list of indices of chosen columns
        Idxs  = lSortIdxs(1:NumChosen);

        % Reorder Q and R
        Q = Q(:,lSortIdxs(1:NumChosen));
        R = R(lSortIdxs(1:NumChosen),lSortIdxs(1:NumChosen));
        
        Idxs = P(Idxs);
        DIdxs = setdiff(1:size(A,2),Idxs);
        
    case 'gsqr'
        % Run nested dissection if required
        if (isstruct(Options.NesDis)) || (Options.NesDis==1),
            if ~isstruct(Options.NesDis),
                % Set the minimum size of clusters
                if isempty(Options.ClusterMinSize),
                    % Find size of smallest chunks in nested dissection
                    lSupports = sum(abs(A)>1e-4,1);
                    lSupportSize = min(lSupports(find(lSupports>0)));
                    Options.ClusterMinSize = round(N/40); %min([lSupportSize*10]);
                    %Options.ClusterMinSize,
                end;
                % Perform hierarchical nested dissection
                [SpQR.NesDis.Perm SpQR.NesDis.cp SpQR.NesDis.cmember]=nesdis(A,'col',[Options.ClusterMinSize]);
                % Build the hierarchical tree corresponding to the nest dissection
                SpQR.NesDis.Tree = BuildTreeFromSeparatorTree( SpQR.NesDis.cp );
            else
                SpQR.NesDis = Options.NesDis;
            end;
            % Thicken up the boundaries between the partitions
            SpQR.NesDis.NewMember = ProcessTreeFromNesDis( SpQR.NesDis.Tree,SpQR.NesDis.cmember,A );
            Options.PivotingLayers = SpQR.NesDis.NewMember;
        end;

        lUniquePivotingLayers = unique(Options.PivotingLayers);
        lNPivotingLayers = length(lUniquePivotingLayers);

        % Allocate memory.
        Q = cell(1,N); Qt = cell(1,N); R = cell(1,N);

        Chosen = zeros(1,N);
        Norms = zeros(1,N);
        Idxs = zeros(1,N);
        NumChosen = 0;
        ips = zeros(N,1);

        % Compute norms
        for j=1:N
            Norms(j) = norm(A(:,j));
            if Norms(j)==0,
                Norms(j) = -1;
            end;
        end
        
        lCurColIdx = 1;
        ltreeChildrenIdxs = [];
        At = A';
        
        for p = 1:lNPivotingLayers,
            if lCurColIdx > Options.StopN,
                break;
            end;
            lLocalNumChosen = 0;
            % Select current allowed columns for pivoting
            lCurPivotingCols = find( Options.PivotingLayers==lUniquePivotingLayers(p) );
            if ~isempty(SpQR.NesDis),
                % Find (j,k) position of current pivoting layer in tree
                for ltreej = length(SpQR.NesDis.Tree):-1:1,
                    ltreek = find( SpQR.NesDis.Tree{ltreej}.Idxs==lUniquePivotingLayers(p) );
                    if ~isempty(ltreek),
                        break;
                    end;
                end;
                % Only indices of points in the children components need to be checked
                ltreeChildrenIdxs = GetLeafIdxs( SpQR.NesDis.Tree, SpQR.NesDis.NewMember , ltreej, ltreek );
            end;
            
            % Select pivot column and orthogonalize
            for j = 1:length(lCurPivotingCols),
                % From among the remaining columns, choose the one with maximum L^2 norm
                [HeapNorm ChosenColumn] = max(Norms(lCurPivotingCols));
                if HeapNorm<=0, break; end;
                ChosenColumn        = lCurPivotingCols(ChosenColumn);
                % Move the chosen column into the cell array Q
                Q{lCurColIdx}       = A(:,ChosenColumn);
                Qt{lCurColIdx}      = Q{lCurColIdx}';
                ComputedNorm        = norm(Q{lCurColIdx});

                % Project the chosen column into the orthogonal complement of the column space of Q
                % Index of columns with which the inner products should be computed
                lNIpIdxs = lCurColIdx-1;
                lIpIdxs = 1:lNIpIdxs;
                if ~isempty(ltreeChildrenIdxs),
                    [lTmp,lintersectIdxs] = intersect( Idxs(1:lCurColIdx-1),ltreeChildrenIdxs );
                    lIpIdxs  = lIpIdxs(lintersectIdxs);
                    lNIpIdxs = length(lIpIdxs);
                end;
                if lNIpIdxs>0,
                    % Compute inner products
                    ips = zeros(1,lNIpIdxs);
                    % Project onto already picked columns
                    for i=1:Options.NumProjections
                        for r=1:lNIpIdxs
                            ip = Qt{lIpIdxs(r)}*Q{lCurColIdx};           % ip = sum(Q{r}.*Q{j});
                            ips(r) = ip;
                            Q{lCurColIdx} = Q{lCurColIdx} - ip*Q{lIpIdxs(r)};
                        end
                        % Update R
                        if i==1,
                            R{lCurColIdx} = sparse(lIpIdxs, ones(lNIpIdxs, 1), ips, N, 1, lNIpIdxs+1);
                        else
                            R{lCurColIdx} = R{lCurColIdx} - sparse(lIpIdxs, ones(lNIpIdxs, 1), ips, N, 1, lNIpIdxs);
                        end
                    end;
                else
                    R{lCurColIdx} = spalloc(N,1,1);
                end;

                % Check for stopping condition
                ChosenNorm = norm(Q{lCurColIdx});
                if (ChosenNorm < Options.StopPrecision)
                    if ((lLocalNumChosen>=min(Options.MinFcnPerCluster,length(lCurPivotingCols))) || ...
                            (ChosenNorm<=0)) && (NumChosen>1), break; end;
                end;

                % Normalize the new column
                Q{lCurColIdx}   = Q{lCurColIdx}/ChosenNorm;
                Qt{lCurColIdx}  = Q{lCurColIdx}';
                R{lCurColIdx}(lCurColIdx) = ChosenNorm;

                % Threshold Q
                if Options.QThreshold>0.0,
                    Q{lCurColIdx} = Q{lCurColIdx} .* (abs(Q{lCurColIdx}) > Options.QThreshold);
                end;

                % Stop if density condition is given and is satisfied
                if(nnz(Q{lCurColIdx})/M > Options.StopDensity),
                    if (lLocalNumChosen>=min(Options.MinFcnPerCluster,length(lCurPivotingCols))) && (NumChosen>1), break; end;
                end;

                % Mark column as chosen and update counters
                Chosen(ChosenColumn)    = true;
                Norms(ChosenColumn)     = -1.0;
                NumChosen               = NumChosen+1;
                lLocalNumChosen         = lLocalNumChosen+1;
                Idxs(NumChosen)         = ChosenColumn;

                % Update the norms of the overlapping columns
                valididxs = find(Norms > 0.0);

                if ~isempty(valididxs)
                    % Compute inner products between the new column of Q and the remaining columns of A
                    ips = At(valididxs,:)*Q{lCurColIdx};

                    for i=1:length(valididxs),
                        reali = valididxs(i);
                        Norms(reali) = sqrt(abs(Norms(reali)^2 - ips(i)^2));
                    end;
                end;

                % Display
                if ~Options.Quiet && rand(1,1) < .25,  fprintf('\b\b\b\b\b%05d', j);  end;

                lCurColIdx = lCurColIdx+1;
                if lCurColIdx > Options.StopN,
                    break;
                end;
            end;            
        end;

        % Sort the columns based on the diagonal of R
        lDiagR = [];
        for l = NumChosen:-1:1,
            lDiagR(l) = R{l}(l);
        end;
        [lTmp,lSortIdxs] = sort(full(abs(lDiagR)),2,'descend');

        % Truncate if stopping based on relative jump in the diagonal of R is requested
        if Options.StopRelJump>0,
            lLastBigJump = NumChosen;
            for j = 2:NumChosen,
                if lDiagR(lSortIdxs(j)) < Options.StopRelJump * lDiagR(lSortIdxs(j-1)),
                    lLastBigJump = j;
                end;
            end;
            if lLastBigJump < NumChosen,
                Chosen(Idxs(lLastBigJump+1:NumChosen)) = false;
                NumChosen = lLastBigJump;
            end;
        end;

        % Truncate list of indices, find list of columns not chosen
        Idxs  = Idxs(lSortIdxs(1:NumChosen));
        DIdxs = sort(find(~Chosen));

        % Convert Q,R to sparse matrices and reorder
        Q = sparse([Q{lSortIdxs(1:NumChosen)}]);
        R = sparse([R{lSortIdxs(1:NumChosen)}]);
        R = R(lSortIdxs(1:NumChosen), :);

        if ~Options.Quiet
            if isfield(Options, 'dwtree')   % special dwtree mode
                fprintf('\b\b\b\b\b');
            else
                fprintf('\b\b\b\b\bdone (%d fcns chosen, %g secs)\n', NumChosen, cputime-T);
            end;
        end;
        
    otherwise
        error(sprintf('gsqr: Unknown QRroutine %s.',Options.QRroutine));
end;

if isempty(Q), Q=0; end;

SpQR.Idxs   = Idxs;
SpQR.DIdxs  = DIdxs;
SpQR.Q      = Q;
SpQR.R11    = R;
SpQR.R12    = Q'*A(:,SpQR.DIdxs);

QROptions = Options;

return;