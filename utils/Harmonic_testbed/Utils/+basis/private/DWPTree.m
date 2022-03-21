function Tree = DWPTree (T, Levels, Precision, Options)

%
% function Tree = DWPTree (T, Levels, Precision, Options)
%
% DWPTREE generates a "diffusion wavelet packet tree" for a given diffusion
% operator.  A diffusion wavelet packet tree stores the bases and operators
% associated with the diffusion wavelet packet transform.
%
% DWPTree returns a cell array represented a tree of subspaces.  Each node
% in the cell array represents a single subspace.
%
% IN:
%
%    Diffusion = diffusion operator represented in the delta basis
%    Levels    = maximum number of levels for the tree
%    Precision = precision for the calculations
%    Options   = a structure for passing algorithm options; the following
%                fields are recognized:
%
%     Construction options:
%      [Wavelets]           : when true, compute the wavelets. Default: true.
%      [WaveletPackets]     : when true, compute wavelets and wavelet packets. Default is false.
%      [GS]                 : string with the name of the QR routine to use; the gs function must have the 
%                               same syntax as gsqr.m. Default: 'gsqr'.
%      [SplitFcn]           : handle to a user specified function for deciding where to
%                               split the wavelet packet spaces
%      [StopFcns]           : don't split spaces with fewer functions than this. Default: 1.
%      [StopRelJump]        : split at a jump in the approximated singular values. This attempts at discovering relevant
%                               scales, corresponding to steps in the eigenvalues of the diffusion.
%      [AdaptiveScales]     : structure with the following fields. Default: 'None'.
%                   Type    : {'None','Fixed','AdaptSparsity'}
%                   Param   : if Type=='Fixed', the power to be taken at each scale
%                             if Type=='AdaptSparsity', sparsity factor w.r.t to standard power not to be 
%                               exceeded from one scale to next.
%
%     Numerics
%      [GSOptions]          : structure containing options to pass to gsqr. Default: [].
%      [OpThreshold]        : Threshold for the entries of Op at each level. Default: 0.
%      [ExtBasisThreshold]  : Threshold for the entries of ExtBasis at each level. Default: 0.
%      [BiorthThresold]     : Thresold for Basis, in order to create the dual basis. Default: 0.
%      [Symm]               : Whether T is symmetric or not. Used only for numerical accuracy purposes. Default: false.
%
%     Extra service options:
%      [ExtendBases]        : when true, DWPTree compute representations of each basis
%                               with respect to the delta basis
%      [Diagnostics]        : compute and display diagnostic information while running. Computationally 
%                               expensive. Default: 0.
%
% OUT:
%
%    Tree      = a cell array representing a diffusion wavelet packet tree
%       Tree{j,k} contains the wavelet packet subspace at frequency 2^j and of index k. k=1 corresponds to
%       scaling functions, k=1 to waveletes, and the other indices to the index of dyadically generated
%       wavelet packet subspaces.
%       The fields of Tree{j,1} are (notation as in the "Diffusion Wavelets" paper):
%           Op      : [T]_{\Phi_{j-1}}^{\Phi_j}
%           Basis   : [\Phi_j]_{\Phi_{j-1}}
%           T{1}    : [T^2]_{\Phi_j}^{\Phi_j}
%           Idxs    : indices of the columns of [T]_{\Phi_{j-1}}^{\Phi_{j-1}} used to produce Basis
%           ExtIdxs : indices of the columns of [T]_{\Phi_{0}}^{\Phi_{0}} used to produce Basis
%       Depending on the options above, it may also contain:
%           ExtBasis: [\Phi_j]_{\Phi_0}
%
%
% Dependences:
%    gsqr.m, fprintf.m
%
% Version History:
%   jcb        2/2006         Original version
%   mm         4/2007         Added a few options, bug fixes, nonsymmetric case...
%   mm         5/2007         Added a few more options
%   mm         12/2008
%
%
% (c) Copyright Yale University, 2006, James C Bremer Jr.
% (c) Copyright Duke University, 2007, Mauro Maggioni
%
%
% EXAMPLES:
%   N=512;e=ones(N,1);T=spdiags([e/4 e/2 e/4],-1:1,N,N);
%   Tree = DWPTree (T, 12, 1e-6, struct('Wavelets',false,'OpThreshold',1e-6,'Symm',true,'GSOptions',struct('StopDensity',1,'Threshold',1e-6)));
%   figure;plot(Tree{4,1}.ExtBasis(:,10))
%   figure;plot(Tree{4,1}.ExtBasis(:,10));hold on;Tp=T^15;plot(Tp(:,Tree{4,1}.ExtIdxs(10)),'r');
%


% Handle options and defaults
if ( nargin<4 )                             Options = [];                           end;
if ~isfield(Options, 'ExtendBases')         Options.ExtendBases = true;             end;
if ~isfield(Options, 'StopFcns')            Options.StopFcns = 1;                   end;
if ~isfield(Options, 'StopRelJump')         Options.StopRelJump = 0.0;              end;
if ~isfield(Options, 'GS')                  Options.GS = @gsqr;                     end;
if ~isfield(Options, 'GSOptions')           Options.GSOptions  = [];                end;
if ~isfield(Options, 'Wavelets')            Options.Wavelets = true;                end;
if ~isfield(Options, 'WaveletPackets')      Options.WaveletPackets = false;         end;
if ~isfield(Options, 'SplitFcn')            Options.SplitFcn = @DefaultSplitFcn;    end;
if ~isfield(Options, 'OpThreshold')         Options.OpThreshold = Precision;        end;
if ~isfield(Options, 'ThresThreshold')      Options.ThresThreshold = Precision;     end;
if ~isfield(Options, 'ExtBasisThreshold')   Options.BasisThreshold = 0;             end;
if ~isfield(Options, 'Symm')                Options.Symm = false;                   end;
if ~isfield(Options, 'ScalePower')          Options.ScalePower = 2;                 end;
if ~isfield(Options, 'AdaptiveScales')      Options.AdaptiveScales.Type = 'None';   end;
if ~isfield(Options, 'BiorthThreshold')     Options.BiorthThreshold = 0.0;          end;
if ~isfield(Options, 'ExtBasisThreshold')   Options.ExtBasisThreshold = 0.0;        end;
if ~isfield(Options, 'Diagnostics')         Options.Diagnostics = false;            end;


% process input
N = size(T,1);

% Initialize the tree structure
if Options.WaveletPackets
    MaxIndex = min([2^floor(Levels/2),N]);  % The maximum possible index of any node in the tree
elseif Options.Wavelets
    MaxIndex = 2;
else
    MaxIndex = 1;
end;
Tree = cell(Levels, MaxIndex);

% Save original T
T_0 = T;

% QROptions specifies parameters for the sparse QR process
QROptions               = Options.GSOptions;
QROptions.StopPrecision = Precision;
QROptions.dwtree        = 1;  % special formatting for gsqr.m output

if Options.Diagnostics,
    Diagnostics.Compr_TandOp_infty(1)= 0;
    Diagnostics.Dens_Basis(1)        = 1/N;
    Diagnostics.Dens_Op(1)           = MatrixDensity(T);
    Diagnostics.Spars_Basis(1)       = N;
    Diagnostics.Spars_Op(1)          = nnz(T);
    Diagnostics.Compr_TandTExt_infty(1) = 0;
    Diagnostics.Dens_Tj(1)           = MatrixDensity(T);
    Diagnostics.Spars_Tj(1)          = nnz(T);
end;

%
% Go through the different scales
%
for j=1:Levels
    % Stop if there's only a small number of functions at this level
    if size(T,2) <= Options.StopFcns, Tree = Tree(1:(j-1),:); break; end;

    % Initialize the tree node corresponding to V_{j}
    Tree{j,1} = struct('Level', j, 'Index',1, 'Basis', [], 'Op', [], 'T', [], 'ExtBasis', [], 'Freq', [(Precision^(2^(-j+1))) 1.0]);

    % Print out the name of the node and its index in the tree structure
    fprintf(sprintf('%s', DWNodeName(j,1)), 15);

    % Perform sparse QR on the columns of the "current" T operator
    fprintf('gsqr: '); TIME = cputime;
    SpQR = feval(Options.GS, T, QROptions);

    % Put the sparse QR information in the tree
    Tree{j,1}.Basis     = SpQR.Q;
    Tree{j,1}.Op        = [SpQR.R11,SpQR.R12];
    Tree{j,1}.OOp       = Tree{j,1}.Op;
    Tree{j,1}.Idxs      = SpQR.Idxs;
    Tree{j,1}.NesDis    = SpQR.NesDis;
    if Options.BiorthThreshold>0.0,
        % Threshold basis functions and mark for biorthogonality        
        BasisThres          = Tree{j,1}.Basis.*(abs(Tree{j,1}.Basis)>Options.BiorthThreshold);
        Tree{j,1}.Op        = (BasisThres'*((BasisThres*BasisThres')\Tree{j,1}.Basis)*Tree{j,1}.Op);
        Tree{j,1}.Basis     = BasisThres;
        Tree{j,1}.Biorth    = true;
    else
        Tree{j,1}.Biorth    = false;        
    end;

    if j>1, Tree{j,1}.ExtIdxs = Tree{j-1,1}.ExtIdxs(SpQR.Idxs);
    else    Tree{j,1}.ExtIdxs = SpQR.Idxs;  end;
    
    % Compute the inverse permutation
    [lTemp,lIdxs]   = sort([SpQR.Idxs,SpQR.DIdxs]);
    Tree{j,1}.Op    = Tree{j,1}.Op(:,lIdxs);
    if Options.OpThreshold > 0.0,
        Tree{j,1}.Op  = Tree{j,1}.Op.*(abs(Tree{j,1}.Op)>Options.OpThreshold);
        Tree{j,1}.OOp = Tree{j,1}.OOp.*(abs(Tree{j,1}.OOp)>Options.OpThreshold);
    end;

    if Options.Diagnostics,
        Diagnostics.Compr_TandOp_infty(j+1)= max(max(abs(T(:,:)-Tree{j,1}.Basis*Tree{j,1}.Op)));
        Diagnostics.Dens_Basis(j+1)        = MatrixDensity(Tree{j,1}.Basis);
        Diagnostics.Dens_Op(j+1)           = MatrixDensity(Tree{j,1}.Op);
        Diagnostics.Spars_Basis(j+1)       = nnz(Tree{j,1}.Basis);
        Diagnostics.Spars_Op(j+1)          = nnz(Tree{j,1}.Op);
    end;
    
    N = size(T,2);                  % Number of functions at previous level
    K = size(Tree{j,1}.Basis,2);    % Number of functions at this level
    fprintf(sprintf('%4d fcns, %4.2f secs', K, cputime-TIME), 25);

    % Compute representations of the dyadic powers T, T^2, T^4 in our new basis
    % Needed only for a wavelet packet transform
    fprintf('T reps: '); TIME = cputime;
    if Options.WaveletPackets,
        Tree{j,1}.T    = cell(j,1);
    elseif Options.Wavelets
        Tree{j,1}.T    = cell(2);
    else
        Tree{j,1}.T    = cell(1);
    end;
    
    % Compute the representation of T in the basis at the current scale
    if ~isstruct(Tree{j,1}.Basis),
        try
        Tree{j,1}.T{1} = Tree{j,1}.Op*Tree{j,1}.Basis;      % MM: This is the RQ product Tree{j,1}.Basis'*T*Tree{j,1}.Basis;
        catch
            1,
        end;
    else
        Tree{j,1}.T{1} = spqr_qmult(Tree{j,1}.Basis,spqr_qmult(Tree{j,1}.Basis,T,3),0);
    end;
    
    if Options.Diagnostics,
        Diagnostics.Compr_TandTExt_infty(j+1) = max(max(T-Tree{j,1}.Basis*Tree{j,1}.T{1}*Tree{j,1}.Basis'));
    end;
    
    % Compute power of T{1} as requested
    switch lower(Options.AdaptiveScales.Type),
        case 'none'            
            Tree{j,1}.Time = 1;
        case 'fixed'
            if j==1,
                if ~isfield(Options.AdaptiveScales,'Param'),
                    Options.AdaptiveScales.Param = 2;
                end;
            end;
            % Raise Op to the requested power
            Tree{j,1}.Op = Tree{j,1}.Basis'*(T^Options.AdaptiveScales.Param);
            Tree{j,1}.Time = Options.AdaptiveScales.Param;
        case 'adaptsparsity'
            if j==1,
                if ~isfield(Options.AdaptiveScales,'Param'),
                    Options.AdaptiveScales.Param = 1;
                end;
            end;
            Tp = T;
            for lTime = 1:10,
                nnz(Tree{j,1}.Op),nnz(Tp),
                if nnz(Tree{j,1}.Op)<Options.AdaptiveScales.Param*nnz(Tp),
                    break;
                else
                    Tree{j,1}.Op = Tree{j,1}.Basis'*(T^lTime);
                    Tree{j,1}.Op = (Tree{j,1}.Op).*(abs(Tree{j,1}.Op)>Options.OpThreshold);
                    Tp = T^lTime;
                    Tp = (Tp).*(abs(Tp)>Options.OpThreshold);
                end;
            end;            
            Tree{j,1}.Time = lTime;
        case 'relativejump'
            while true,
                SpQR= feval(Options.GS, T, QROptions);
                if (length(SpQR.Idxs)>size(T,1)-3) || (size(T,1)<8),
                    break;
                else
                    % Estimate the power of T that will send the first step below the gap below threshold
                    lTime = GetPowerForSplitBelowThreshold(SpQR.R11, QROptions.StopPrecision, Options.StopRelJump);
                    T = T^lTime;
                end;
            end;
    end;
    
    if (Options.Symm) && (~Tree{j,1}.Biorth),                                        % Iron out numerical wrinkles
        Tree{j,1}.T{1} = 0.5*(Tree{j,1}.T{1}+Tree{j,1}.T{1}');
    end;
    if Options.OpThreshold > 0.0
        Tree{j,1}.T{1} = (Tree{j,1}.T{1}).*(abs(Tree{j,1}.T{1})>Options.OpThreshold);
    end;
    if Options.WaveletPackets
        for r=1:j-1     % Representation of various powers of T for the packet subspaces
            Tree{j,1}.T{r+1} = Tree{j,1}.Basis'*Tree{j-1,1}.T{r}*Tree{j,1}.Basis;
        end;
    end;
    fprintf(sprintf('%4.2f secs', cputime-TIME), 14);    
    
    % Square the operator T for the next level
    fprintf('T^2: '); TIME=cputime;
    if Options.ScalePower == 2,
        if Options.Symm,
            T = Tree{j,1}.Op*Tree{j,1}.Op';            % Compute and store T^2 w.r.t. the latest bass. No need to enforce symmetry apparently.
        else
            T = Tree{j,1}.T{1}^2;                      % Compute and store T^2 on the new basis, general (non-symmetric) case, but less precise in the symmetric case
        end;
    else
        T = Tree{j,1}.T{1};                            % Compute and store T on the new basis, general (non-symmetric) case, but less precise in the symmetric case
    end;
    % Threshold T
    if Options.OpThreshold>0,
        T = T.*(abs(T)>Options.OpThreshold);
    end;
    
    if Options.Diagnostics,
        Diagnostics.Dens_Tj(j+1)     = MatrixDensity(T);
        Diagnostics.Spars_Tj(j+1)    = nnz(T);
    end;
    
    fprintf(sprintf('%4.2f secs\n', cputime-TIME), 14); %fprintf('freq: [%g %g]', Tree{j,1}.Freq); %fprintf('dens: %.4f\n',nnz(Tree{j,1}.T{1})/numel(Tree{j,1}.T{1}));

    % Compute representation of current basis w.r.t. delta-basis
    if Options.ExtendBases
        if j==1
            Tree{j,1}.ExtBasis = Tree{j,1}.Basis;
        else
            if ~isstruct(Tree{j,1}.Basis),
                Tree{j,1}.ExtBasis = Tree{j-1,1}.ExtBasis*Tree{j,1}.Basis;
            else
                Tree{j,1}.ExtBasis = spqr_qmult(Tree{j,1}.Basis,Tree{j-1,1}.ExtBasis,3);
            end;
            if Options.ExtBasisThreshold>0.0,
                Tree{j,1}.ExtBasis = Tree{j,1}.ExtBasis.*(abs(Tree{j,1}.ExtBasis)>Options.ExtBasisThreshold);
            end;
        end;
    end;
    
    % Initialize the node for the wavelet space W_j
    Tree{j,2} = struct('Level', j, 'Index', 2, 'Basis', [], 'Op', [], 'T', [], 'ExtBasis', [], 'Freq', []);

    if j==1 Tree{j,2}.Freq   = [0 Tree{j,1}.Freq(1)];
    else    Tree{j,2}.Freq   = [Tree{j-1,2}.Freq(2) Tree{j,1}.Freq(1)]; end;

    if Options.Diagnostics,
        % Plot diagnostic information
        figure(1);subplot(1,3,1);
        plot(1:(j+1),log10(abs([Diagnostics.Compr_TandOp_infty;Diagnostics.Compr_TandTExt_infty])+eps)');grid on;axis tight;
        title('Compression accuracy');
        legend({'||T-Op||_\infty','||T_0^{2^j}-T_{j,ext}||_\infty'},'Location','Best');
        subplot(1,3,2);
        semilogy(1:(j+1),[Diagnostics.Spars_Basis;Diagnostics.Spars_Op;Diagnostics.Spars_Tj]');grid on;axis tight;
        legend({'Basis','Op','T_j'},'Location','Best'); title('Sparsity');
        subplot(1,3,3);
        plot(1:(j+1),[Diagnostics.Dens_Basis;Diagnostics.Dens_Op;Diagnostics.Dens_Tj]');grid on;axis tight;
        legend({'Basis','Op','T_j'},'Location','Best'); title('Density');
        drawnow;
    end;
    
    % Compute the wavelets
    if (Options.Wavelets) || (Options.WaveletPackets),
        % Print out node name and position in tree
        fprintf(sprintf('%s', DWNodeName(j,2)), 15);

        % SpQR on the columns of I-Q*Q' in order to choose basis for the wavelet space
        fprintf('gsqr: '); TIME = cputime;
        W                   = speye(N) - Tree{j,1}.Basis*Tree{j,1}.Basis';
        % SpQR options
        QRwOptions          = Options.GSOptions;
        QRwOptions.StopN    = N-K;                    % Choose the "right" number of columns
        if QRwOptions.StopN > 0,
            QRwOptions.dwtree   = 1;                      % Forces special formatting for gsqr.m output
            SpQR                = feval(Options.GS, W, QRwOptions);
            Tree{j,2}.Basis     = SpQR.Q;
            Tree{j,2}.Idxs      = SpQR.Idxs;
            % Store the "extended indices"
            if j>1, Tree{j,2}.ExtIdxs = Tree{j-1,1}.ExtIdxs(Tree{j,2}.Idxs);  end;

            fprintf(sprintf('%4d fcns, %4.2f secs', N-K, cputime-TIME), 25);

            % Compute representation of T for the wavelet space as well
            if Options.WaveletPackets
                fprintf('T reps: '); TIME = cputime;
                for r = 1:j-1
                    Tree{j,2}.T{r} = Tree{j,2}.Basis'*Tree{j-1,1}.T{r}*Tree{j,2}.Basis;     % representation of T w.r.t. Sigma_1
                    if Options.Symm, % Iron out numerical wrinkles
                        Tree{j,2}.T{r} = 0.5*(Tree{j,2}.T{r}+Tree{j,2}.T{r}');
                    end;
                end
                fprintf(sprintf('%4.2f secs', cputime-TIME), 19);
            end;
            
            % Compute "extended bases" (representations of the bases at this level w.r.t. delta basis)
            if Options.ExtendBases,
                if j==1, Tree{j,2}.ExtBasis = Tree{j,2}.Basis;
                else     Tree{j,2}.ExtBasis = Tree{j-1,1}.ExtBasis*Tree{j,2}.Basis; end;
            end;
        else
            Tree{j,2}.Basis     = [];
            Tree{j,2}.Idxs      = [];
            Tree{j,2}.ExtIdxs   = [];
        end;

        % We don't have to square T
        fprintf('', 14);
        fprintf('freq: [%g %g]', Tree{j,2}.Freq);
        fprintf('\n');
    end;

    % Split any wavelet nodes at the preceeding level
    if Options.WaveletPackets,
        for node = 2: 2^ceil((j-2)/2)
            if (~isempty(Tree{j-1,node})) && (~isempty(Tree{j-1,node}.T))
                N = size(Tree{j-1,node}.T{1}, 2);
                % Check to see if the node is worth splitting
                if N > Options.StopFcns
                    % Node is the index of the parent node at the preceeding level
                    % LeftNode is the index of its left child
                    % RightNode is the index of its right child
                    LeftNode  = 2*node-1;
                    RightNode = 2*node;
                    % Determine what precision we should use for splitting the function
                    SplitFreq = feval(Options.SplitFcn, Tree{j-1,node}, Precision);
                    Tree{j,LeftNode}.Freq  = [Tree{j-1, node}.Freq(1) SplitFreq];
                    Tree{j,RightNode}.Freq = [SplitFreq Tree{j-1, node}.Freq(2)];
                    % TWPower adjusts for the fact that we are processing T^(2^TPower)
                    TPower = (2)^(length(Tree{j-1,node}.T)-1);
                    % SpQR on the columns of T
                    fprintf(sprintf('%s', DWNodeName(j,LeftNode)), 15);
                    fprintf('gsqr: ');TIME = cputime;
                    % Parameters for SpQR
                    QRwOptions                  = Options.GSOptions;
                    QRwOptions.StopPrecision    = SplitFreq^(TPower);
                    QRwOptions.dwtree           = 1;
                    SpQR                        = feval(Options.GS,Tree{j-1,node}.T{1}, QRwOptions);
                    Tree{j,LeftNode}.Op         = SpQR.R11;
                    Tree{j,LeftNode}.Basis      = SpQR.Q;
                    Tree{j,LeftNode}.Idxs       = SpQR.Idxs;
                    
                    % Compute ExtIdxs
                    try
                    if j>1, Tree{j,LeftNode}.ExtIdxs = Tree{j-1,node}.ExtIdxs(Tree{j,LeftNode}.Idxs); end;
                    catch
                        1,
                    end;
                    K = size(Tree{j,LeftNode}.Basis,2);
                    fprintf(sprintf('%4d fcns, %4.2f secs', K, cputime-TIME), 25);

                    % Compute T reps
                    fprintf('T reps: '); TIME = cputime;
                    Tree{j,LeftNode}.T = cell(length(Tree{j-1,node}.T)-1,1);

                    for r=1:length(Tree{j-1,node}.T)-1,
                        Tree{j,LeftNode}.T{r} = Tree{j,LeftNode}.Basis'*Tree{j-1,node}.T{r+1}*Tree{j,LeftNode}.Basis;
                        if Options.Symm, % Iron out numerical wrinkles
                            Tree{j,LeftNode}.T{r} = 0.5*(Tree{j,LeftNode}.T{r}+Tree{j,LeftNode}.T{r}');
                        end;                        
                    end;                    
                    % We don't have to square T
                    fprintf(sprintf('%4.2f secs', cputime-TIME), 19); fprintf('', 14);
                    fprintf('freq: [%g %g] \n', Tree{j,LeftNode}.Freq);
                    
                    fprintf(sprintf('%s', DWNodeName(j,RightNode)),15);
                    fprintf('gsqr: '); TIME = cputime;

                    % Projection onto wavelet packet subspace, and its SpQR
                    W                        = speye(N) - Tree{j,LeftNode}.Basis*Tree{j,LeftNode}.Basis';
                    QRwOptions               = Options.GSOptions;
                    QRwOptions.StopN         = N-K;
                    if QRwOptions.StopN > 0,
                        QRwOptions.dwtree        = 1;
                        SpQR                     = feval(Options.GS, W, QRwOptions);
                        Tree{j,RightNode}.Basis  = SpQR.Q;
                        Tree{j,RightNode}.Idxs   = SpQR.Idxs;
                        if j>1, Tree{j,RightNode}.ExtIdxs = Tree{j-1,1}.ExtIdxs(Tree{j,RightNode}.Idxs); end;
                        fprintf(sprintf('%4d fcns, %4.2f secs', K, cputime-TIME), 25);

                        % Compute T reps
                        fprintf('T reps: '); TIME = cputime;

                        if ~isempty(Tree{j,RightNode}.Basis)
                            Tree{j,RightNode}.T = cell(length(Tree{j-1,node}.T)-1,1);
                            for r=1:length(Tree{j-1,node}.T)-1
                                Tree{j,RightNode}.T{r} = Tree{j,RightNode}.Basis'*Tree{j-1,node}.T{r+1}*Tree{j,RightNode}.Basis;
                                if Options.Symm, % Iron out numerical wrinkles
                                    Tree{j,RightNode}.T{r} = 0.5*(Tree{j,RightNode}.T{r}+Tree{j,RightNode}.T{r}');
                                end;
                            end;
                        else
                            Tree{j,RightNode}.T = {};
                        end;

                        % Compute the "extended basis" if needed
                        if Options.ExtendBases
                            Tree{j,LeftNode}.ExtBasis  = Tree{j-1,node}.ExtBasis*Tree{j,LeftNode }.Basis;
                            Tree{j,RightNode}.ExtBasis = Tree{j-1,node}.ExtBasis*Tree{j,RightNode}.Basis;
                        end;
                    else
                        Tree{j,RightNode}.Basis     = [];
                        Tree{j,RightNode}.Idxs      = [];
                        Tree{j,RightNode}.ExtIdxs   = [];
                        Tree{j,RightNode}.T         = cell(length(Tree{j-1,node}.T)-1,1);
                    end;
                    
                    fprintf(sprintf('%4.2f secs', cputime-TIME), 19); fprintf('', 14);
                    fprintf('freq: [%g %g]', Tree{j,LeftNode}.Freq);
                    fprintf('\n');
                end;
            end;
        end;
    end;
end;

fprintf('\n');

return;


function SplitFreq = DefaultSplitFcn(Node, Precision)
% function [SplitFreq NFcns] = DefaultSplitFcn(Node, Precision)
%
% Averages the endpoints of the approximate frequency range of the node
% Node.
%
% In:
%    Node       = Node to split
%    Precision  = Precision of the process
% Out:
%    SplitFreq  = Frequency value to split at.
%

SplitFreq = sqrt(Node.Freq(1)*Node.Freq(2));


return;



function vPower = GetPowerForSplitBelowThreshold(cR, cPrecision, cRelJump)

lN = size(cR,2);

lLastBigJump = lN;

for j = 2:lLastBigJump,
    if cR(j,j) < cRelJump*cR(j-1,j-1),
        lLastBigJump = j;
    end;
end;

vPower = ceil(log(cPrecision)/log(cR(lLastBigJump,lLastBigJump)));

return;