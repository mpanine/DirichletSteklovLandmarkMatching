function [Phi] = compressed_manifold_mode(X,T,W,A,numEigs,mu)
% |-------------------------------------------------------|
% |CODE FROM https://github.com/matansel/geometric_L1_norm|
% |-------------------------------------------------------|
A=diag(A);
nv = length(X);

opt.mu=mu/nv; %sparsity constant, mu=100 by default

opt.xi = 1e5; %orthonormality constant
opt.nIter = 50; %Inner iterations of IRLS
opt.nEV=numEigs; % #eigenvectors to compute

opt.epsilonStart = 1e-4;%also defined by default
opt.epsilonTresh = 1e-7;%also defined by default
% opt.epsilonStart = 1e-1;% for faster computations
% opt.epsilonTresh = 1e-2;% for faster computations
opt.epsilonStep = 2; %epsilon=epsilon/opt.epsilonStep at each iteration

opt.A = A; %Your column Vector of Vertex Area
opt.W = W; %Laplacian weights (sparse)
opt.verbose = 1;
%Critical Step!
%Column Vector of weights such that opt.func(epsilon,ww,V1,A) = diag(Vi^2*da)
%Matan's l1 weights:
% opt.func = @(epsilon,ww,V1,A)(abs(ww)./((sqrt(epsilon+V1.^2)).*A));
opt.func = @(epsilon,ww,V1,A)(max(0,(ww)./(sign(V1+epsilon).*(sqrt(epsilon+V1.^2)).*A)));

%Yoni's l1 weights:
% opt.func = @(epsilon,ww,V1,A)(1./(2*sqrt(epsilon+(V1).^2)));

%%

mu=opt.mu;
nIter = opt.nIter;
nEV=opt.nEV;
xi = opt.xi;
epsilonStart = 1e-3;
if isfield(opt,'epsilonStart')
  epsilonStart = opt.epsilonStart;
end
epsilonTresh = 1e-5;
if isfield(opt,'epsilonTresh')
  epsilonTresh = opt.epsilonTresh;
end
epsilonStep =2;
if isfield(opt,'epsilonStep')
  epsilonStep = opt.epsilonStep;
end
A = opt.A;W = opt.W+1e-8*speye(length(A));
Am = sparse(1:length(A), 1:length(A), A);

% goalFuncTot = cell(1,nIter);
% goalFunc = zeros(1,nIter);
% Potential = zeros(nv,nEV);
Phi= zeros(nv,nEV);
D = zeros(1,nEV);

time = 0;

tic;fprintf('Computing CMM...');
for jj = 1:nEV
  
%   tic
  V = sparse(1:length(A), 1:length(A), 1);
  %     V = Am;
%   disp(['CMM EigenVec #',num2str(jj)])
  
  epsilon = epsilonStart;
  
  for ii=1:nIter
    op.issym = 1;
    op.isreal = 1;
    
    if jj~=1
      [V1,D1] = eigs((mu)*V+Phitr, Am, 1,'sm',op);
      
      [~,p] = chol((mu*Am)*V+Phitr);
%       if p>0
%         disp('CMM Error: Matrix is not PSD!!!!');
%       end;
    else
      [~,p] = chol(-W+mu*V);
%       if p>0
%         disp('CMM Error: Matrix is not PSD!!!!');
%       end;
      [V1,D1] = eigs(-W+mu*V, Am, 1,'sm',op);
    end
    
%     goalFunc(ii) = -V1'*(W)*V1 +mu*sum(abs(A.*V1));
    
    ww = (l1_weight(T,X,V1))';
    
    %Here V=V^2 of the functional
    %         Vtmp = opt.func(epsilon,ww,V1,A);
    %    matan:     Vtmp = (max(0,(ww)./(sign(V1+epsilon).*(sqrt(epsilon+V1.^2)).*A)));
    Vtmp =(1./(2*sqrt(epsilon+(V1).^2)));
    V = sparse(1:length(A), 1:length(A),Vtmp );
    
    if epsilon > epsilonTresh
      epsilon = epsilon/epsilonStep;
    end
  end
%   time = time+toc;
  
%   goalFuncTot{jj} = goalFunc;
%   potential = full(diag(V));
  %Enforce Sparsity (can remove but increase computation time)
  V1(abs(V1)< 1e-6) = 0;
  Phi(:,jj) = V1;
  D(jj) = D1;
  
  Phit = Phi(:,1:jj)*(Phi(:,1:jj))';
  if length(find(Phit)) >size(V1,1)^2/3
    Phitr = -W + (xi*Am*(Phit)*Am);
  else
    Phitr = -W + sparse(xi*Am*(Phit)*Am);
  end
%   Potential(:,jj) = potential;
  
end
fprintf(sprintf('done\n'));toc;

% optResult.goalFuncTot = goalFuncTot;

% optResult.time = time;
% [obj_fun,const_viol] = compressed_modes_objective(T,X,-W,diag(A),opt.mu,Phi);


end