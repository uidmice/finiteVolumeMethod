
function [r1,r2] = two2d (varargin)
%TWO2D: numeric simulation of the dynamics of two species in 2D
%   [r1,r2] = two2d(r1_0,r2_0, l, W, dt, T) 
%   Input:
%       r1_0,r2_0:  initial density of the two species, NxN matrix
%       l:          domain[-l,l]x[-l,l]
%       W:          interacting potentials. 
%                   W is a structure with fields 'W11','W12','W22' and 'W21':
%       dt:         time step
%       T:          simulation time. Total #iterations = T/dt
%   Output: 
%       r1,r2:      densities of the two species at t = T, NxN matrix
%
%   Optional parm:
%   [r1,r2] = two2d(.. ,H): H is a symbolic function for internal energy as a
%   function of the density. Default H(r) = 0
%
%   [r1,r2] = two2d(.. ,V) optionally sets the environmental confinement
%   potential V, which is a NxN matrix. Default: 0
%
%   [r1,r2] = two2d(.. ,e) sets the diffusion coefficient for some e > 0.
%   Default e = 0
%
%   [r1,r2] = two2d(.. ,'v') or [r1,r2] = two2d(.. ,'V') enables visual display
%   during the simulation. Default disabled.
%   
%   [r1,r2] = two2d(.. ,'solver') where 'solver' sets the numeric method used
%   for ODE. Possible options: 'Euler', 'SSPRK3'. Default 'SSPRK3'.


narginchk(6,11);

[r1_0,r2_0, l, W11,W12,W22,W21, H, V, dt, T, e,v,s] = arguments(varargin{:});
N = length(r1_0);
dx = 2*l/N;
K11 = getKernal(W11, dx, N);
K12 = getKernal(W12, dx, N);
K21 = getKernal(W21, dx, N);
K22 = getKernal(W22, dx, N);
dH = diff(H);
L = @(X) drho2(X, K11,K12,K21,K22, dH, V, dx, e);
nt = ceil(T/dt);
R = zeros([N,N,2]);
R(:,:,1) = r1_0;
R(:,:,2) = r2_0;

X = -l+dx/2:dx:l-dx/2;
X = meshgrid(X);
Y = X.';

for i = 1:nt
    if s=='s'
        R = SSPRK3(R,L,dt);
    else
        R = Euler(R,L,dt);
    end 
    
    if v
        draw(X,Y,R(:,:,1),R(:,:,2),[-l l -l l 0 Inf]);
        title(sprintf('t = %.3f', i*dt));
        drawnow;
    end
end

r1 = R(:,:,1);
r2 = R(:,:,2);

end

function [K] = getKernal (W, dx, N)

%{
input:
    W:              interaction potential function
    dx:             space size
    N:              number of points 
output:
    Wij:            help matrix for calculation W(i, j)
%}

   xi = 0:N-1;
   xi = repmat(xi, 1,N)*dx;
   yi = 0:N-1;
   yi = repelem(yi, N)*dx;
   
   K = double(W(xi, yi));
   
   t = logical(isnan(K)+isinf(K));
   if any(t)
       in = find(double(t));
       a = 0.5773503*dx/2;
       if length(in)==1
           K(in) =(W(xi(in)+a,yi(in)+a)+W(xi(in)-a,yi(in)+a)+W(xi(in)+a,yi(in)-a)+W(xi(in)-a,yi(in)-a))/4;
       else
           for i = range (length(in))
               ind = in(i);
               K(ind) =(W(xi(ind)+a,yi(ind)+a)+W(xi(ind)-a,yi(ind)+a)+W(xi(ind)+a,yi(ind)-a)+W(xi(ind)-a,yi(ind)-a))/4;
           end
       end
   end
   K = reshape(K, N, N);
   
   hf = flip(K(:,2:end),2);
   vf = flip(K(2:end,:));
   f = flip(vf(:,2:end),2);
   K = [f, vf;hf,K];
end

function [r1_0,r2_0,l, W11,W12,W22,W21, H, V, dt, T, e,v, s] = arguments(varargin)
r1_0 = varargin{1};
r2_0 = varargin{2};
l = varargin{3};
W = varargin{4};
dt = varargin{5};
T = varargin{6};

W11 = W.W11;
W12 = W.W12;
W22 = W.W22;
W21 = W.W21;

N = length(r1_0);
syms x H(x)

% Optional input defaults
e = 0;                      % No diffusion term
v = false;                  % No visible update during simulation
V = zeros(N, N);            % No confinement potential
H(x) = 0;                   % No internal energy
s = 's';

% Loop over optional arguments
for i = 7:nargin
    a = varargin{i};
    if ischar(a) && isscalar(a) && lower(a) == 'v'
        v = true;
    elseif isreal(a) && isscalar(a) && isfinite(a)
        e = a;
    elseif isa(a,'symfun')
        H = a;
    elseif isstring(a) && ~strcmp(a,'Euler')
        s = 'e';
    elseif ismatrix(a) && length(a)==N
        V = a;
    elseif ismatrix(a) && length(a)~=N
        error("The shape of the input confinement energy matrix is incompatible");
    else
        error('arguments:nonsense','Failed to interpret argument #%d.',i);
    end
                
end
end






