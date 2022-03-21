function[Phi,objfn,iter]=acmm(W,A,mu,K,varargin)
% 
% Accelerated calculation of compressed manifold modes
%
% Calculates K compressed manifold modes of an operator of the 
% form A^{-1}W (such as a Laplace-Beltrami operator) with respect to the
% parameter mu.
%
% Correspondence address:
% Kevin Houston, School of Mathematics, University of Leeds, Leeds, LS2 9JT, United Kingdom
% k.houston@leeds.ac.uk
% http://www.maths.leeds.ac.uk/~khouston/
%
% Creative commons license: BY-SA
% (You are free to use this in a commercial program but I would like you to notify me.) 
%
% Inputs:
% W : Weights matrix of the operator. Eg, the cotan operator.
% A : Mass matrix (or Area matrix) of the operator. Eg, Voronoi areas.
%     This matrix does not need to diagonal (also known as lumped).
% mu : The compression parameter. Low mu means no compression and hence 
%      support of the modes is large. High mu is high compression meaning
%      small support.
% K : Number of modes to calculate.
% Phi_init : Optional input for an initialization for Phi. Phi is randomly
%            generated if this option is not used.
% maxiter : Optional input to limit the number of iterations. Set at 10000
%
% Outputs:
% Phi : Matrix with K columns and number of rows equal to number of rows
%       and columns of W. Contains the modes as columns.
% objfn : The value of the objective function.
% iter : Number of iterations taken
%
%
% References
% Code is implementation of material in Kevin Houston, Compressed Manifold Modes: Fast Calculation and Natural Ordering, preprint University of Leeds April 2015
%
% [Boyd] Stephen Boyd, Neal Parikh, Eric Chu, Borja Peleato, and Jonathan Eckstein. Distributed optimization and statistical learning via the alternating direction method of multipliers. Found. Trends Mach. Learn., 3(1):1–122, January 2011.
% [Neu] Thomas Neumann, Kiran Varanasi, Christian Theobalt, Marcus Magnor, and Markus Wacker. Compressed manifold modes for mesh processing. Computer Graphics Forum (Proc. of Symposium on Geometry Processing SGP), 33(5):1– 10, July 2014.

%
% EXAMPLES
%
%path(path, 'toolbox_graph/');
%path(path, 'meshes/');
%
%FilePath = 'meshes/';
%FileName = 'Lshape';K=15; mu=0.01;
%[vertices,faces] = read_mesh([FilePath FileName '.ply']);  % read ply.
%[vertices,faces] = removeRedundantPoints(vertices,faces);
%[vertices,faces] = patchslim(vertices,faces);
%[W,A,vertices,faces] = FEM(vertices,faces,'linear'); % Construct LBO from the mesh
%
%[Phi_init,eigval] = eigs(W,A,K,'sm'); % Initialization based on eigenfunctions
%[Phi,ob]=acmm(W,A,mu,K,Phi_init,200); % Use at most 200 iterations
%ob; % Display the objective function
%[Phi,lambda]=sortcmm(W,D,mu,Phi); % Sort the modes into ascending compressed eigenvalue
%
%fn=eigenfunction(:,5)% Take fifth mode
%figure;
%set(gcf,'color',[1,1,1]);
%trisurf(faces,vertices(:,1),vertices(:,2),vertices(:,3),fn);
%
%
%Phi_init=rand(size(W,1),K);%/100000; % Random initialization
%[Phi2,obj,iter]=acmm(W,A,mu,K,Phi_init); % Use the above random initialization
%obj; % give the objective funtion
%iter; % give the number of iterations used
%
%Phi3=acmm(W,A,mu,K); % Use a random initialization of numbers between 0 and 1, max number of iterations is 10000


%
% ---- Start of code ----
%
N=size(W,1); % Number of vertices / Size of operator, ie W is NxN matrix

% Deal with extra arguments
if nargin>6
    disp('Too many inputs');
    return
end

Phi_init=[0];
maxiter=0;

if nargin==5 
    if size(varargin{1})==[1 1]
        maxiter=varargin{1};
    else
        Phi_init=varargin{1};
    end
end

if nargin==6
    if size(varargin{1})==[1 1]
        maxiter=varargin{1};
        Phi_init=varargin{2};
    else
        Phi_init=varargin{1};
        maxiter=varargin{2};
    end
end

if size(maxiter)~=[1 1]
    disp('Maximum iterations should be a number'); return;
end
if maxiter==0 
    maxiter=10000;
end

if Phi_init==0
    Phi=(0.5*rand(size(W,1),K))+0.5;%-rand(size(W,1),K); % Create a random initialization
else
    Phi=Phi_init;
end

 
if (size(Phi,1)~=N) | (size(Phi,2)~=K)
    disp('The initialization matrix is wrong size or possibly was not really a matrix'); return;
end

% Begin main algorithm

rho=1; % This is the penalty parameter for the ADMM method. See [Boyd] section 2.3

auto_adjust_penalty=true; % False will remove the autoadjust. True is best as the algorithm converges significantly faster
rho_adjust=2.0; % The factor by which rho changes 
rho_adjust_sensitivity=5; % Factor difference between rnorm and snorm
tol_abs=1.0e-8; % Tolerance for stopping 
tol_rel=1.0e-6; % Tolerance for stopping

E=Phi;
S=Phi;
E_old=E;
S_old=S;
E_hat=E;
S_hat=S;
UE=zeros(N,size(Phi,2)); 
US=zeros(N,size(Phi,2));
UE_old=UE;
US_old=US;
UE_hat=UE;
US_hat=US;

c=1;
eta=0.999;

c_old=eta*c;
alpha_old=1;
alpha=1;

Hsolve=(rho*speye(N))+2*W; % Not -2*W as our W is pos semi def not neg semi def as in Neu

%tracesparse_old=0;
iter=1;
rnorm=0;
snorm=0;

while (iter<=maxiter)
    % Update Phi 
    % Equations (10),  (14-16) in [Neu]
    Y=S_hat-US_hat+E_hat-UE_hat; 
    dy=transpose(Y)*A*Y;
    [V,W1,VT]=svd(dy);
    Wsqrtinv=diag(1./sqrt(diag(W1)));
    Phi=Y*V*Wsqrtinv*transpose(VT); 
  
    % Update S 
    % Equation (12) in [Neu]
    S_old=S;
    S=sign(Phi+US_hat).*max(zeros(N,K),bsxfun(@minus,abs(Phi+US_hat),mu/rho)); % Eqn 18 in [Neu], S_{ij}<-prox etc 
    
    % Update E 
    % Equations (11) and (17) in [Neu]
    E_old=E;
    rhs=rho.*(Phi+UE_hat);
    E=Hsolve\rhs;
    
    % Update UE, US
    UE=UE_hat+(Phi-E);
    US=US_hat+(Phi-S);
    
    %
    % ---- Acceleration phase ----
    %
    c=rho*(sum(sum((UE-UE_hat).^2+(US-US_hat).^2)))+rho*(sum(sum((E-E_hat).^2+(S-S_hat).^2)));
    
    % Compute primal and dual residual
    snorm = sqrt(sum(sum(((rho*(S-S_old)).^2))+sum(sum((rho*(E-E_old)).^2))));
    rnorm = sqrt(sum(sum(((Phi-S).^2)+(Phi-E).^2)));
    
     if (c/c_old)>=eta
        alpha_old=1;
        %disp('Restarted!');
     end
        alpha=(1+sqrt(1+4*alpha_old^2))/2; 
        scale_alpha=(alpha_old-1)/alpha; 
        
        E_hat=E+scale_alpha*(E-E_old);
        S_hat=S+scale_alpha*(S-S_old);
        UE_hat=UE+scale_alpha*(UE-UE_old);
        US_hat=US+scale_alpha*(US-US_old);    
        UE_old=UE;
        US_old=US;
    
        alpha_old=alpha;
        c_old=c;
        
        if auto_adjust_penalty==true 
            if rnorm> rho_adjust_sensitivity*snorm
                rho = rho*rho_adjust;
                Hsolve=(rho*speye(N))+2*W; 
                UE_old = UE_old/rho_adjust;
                US_old = US_old/rho_adjust;
                UE_hat = UE_hat/rho_adjust;
                US_hat = US_hat/rho_adjust;
                %disp(['Changing rho up to ',num2str(rho)]);
            else
                if snorm>rho_adjust_sensitivity*rnorm
                rho = rho/rho_adjust;
                Hsolve=(rho*speye(N))+2*W;
                UE_old = UE_old*rho_adjust;
                US_old = US_old*rho_adjust;
                UE_hat = UE_hat*rho_adjust;
                US_hat = US_hat*rho_adjust;
                %disp(['Changing rho down to ',num2str(rho)]);
               
                end
            end
        end
     
 eps_pri=sqrt(2*N)*tol_abs+tol_rel*(max(sqrt(sum(sum(Phi.^2))),sqrt(sum(sum(E.^2))+sum(sum(S.^2)))));
 eps_dual = sqrt(N)*tol_abs+tol_rel*sqrt(sum(sum((rho*UE).^2+(rho*US).^2)));
    
    if rnorm<eps_pri & snorm<eps_dual
        %disp('Converged')
         objfn=trace(transpose(E)*W*E)+mu*norm(Phi,1);
        break
    end
   iter=iter+1;
end
objfn=trace(transpose(E)*W*E)+mu*norm(Phi,1);
iter=iter-1;