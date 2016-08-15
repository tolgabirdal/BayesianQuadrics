%
% Generate vertices and faces using the unit circle as a starting point.
%
% The function currently only plots ellipsoids.
%
% The quadric is of the form,
%   x'Ax + b'x + c
%
% Input:
%   A - 3x3 matrix
%   b - 3x1 vector
%   c - a scalar
%   N - square root of the number of points on the surface
%   MinSigma - prevent eigenvalues of A from going below this number.
%
% Output:
%   f - 3xK matrix of faces
%   v - 3xN matrix of vertices
% 

function [f, v] =  PlotQuadricFromUnitCircle(A, b, c, N=50, MinSigma=1e-8)
  [M,~] = size(A);
  assert(all([M M]==[size(A,2), size(b,1)]),...
         "PlotQuadricFromUnitCircle: Incorrect sizes" );
  
  if( M==3 )
    [v1,f] = GenerateUnitSphere( N );
  elseif( M==2 )
    theta = linspace(-pi, pi, N);
    v1 = [sin(theta); cos(theta)];
    f=[];
  end

  mu = -0.5*(A\b);
  d  =  mu'*A*mu - c;
  [~,o] = chol(A);
  if(o~=0)
    warning("PlotQuadricFromUnitCircle: Matrix is not spd");
    [U1,S1] = eig(A);
    S1(S1<0)=MinSigma;
    A = U1*S1*inv(U1);
  end
  
  [U S] = eig(A);
  if(d<0)
    d = -d;
    U = -U;
  end
  
  C = (U*sqrt(S))';
  v = sqrt(d)*inv(C)*v1 + repmat(mu,1,size(v1,2));
endfunction
