%
% Generate a random quadric surface. Solutions of,
%   (x-mu)' A (x-mu) = d  
%
%  A is drawn from a wishart distribution (unless specified)
%  mu is drawn from a multivariate gaussian
%  d is drawn from a multivariate gaussian
%
% Input:
%  V      - 3x3 matrix a prior for the wishart function
%  n      - The degrees of freedom for the wishart distribution
%  Sigmab - The covariance of the center
%  Sigmac - The variance for the threshold
%  alg    - can be one of "wish", "gauss", "uniform". 
%           In the first case the first two parameters will be used
%
function [A, b, c] = GenerateQuadric(V, n, Sigmab, Sigmac, alg="wish")
  [M, N] = size(V);
  assert( all(M==[N size(Sigmab)]), "GenerateQuadric: Sizes do not match");
  switch alg
    case "wish"
      A = wishrnd(V,n);
    case "gauss"
      a = randn(6,1); 
      A=[a(1), a(2), a(3); ...
         a(2), a(4), a(5); ...
         a(3), a(5), a(6)];
    case "uniform"
      a = rand(6,1); 
      A=[a(1), a(2), a(3); ...
         a(2), a(4), a(5); ...
         a(3), a(5), a(6)];
  end
  mu = mvnrnd( M, Sigmab, 1 );
  b = -2*A*mu';
  d = normrnd( 0, Sigmac );
  c  =  mu*A*mu' -d;
endfunction
