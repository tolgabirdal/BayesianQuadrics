%
% Fit a general quadric to data.
%
% The code is from the paper,
%   "Fitting quadrics with a Bayesian prior" - D Beale, YL Yang, N Campbell, D Cosker, P Hall - Computational Visual Media, 2016
%
% Input:
%   Data   - MxN matrix of datapoints
%   Params - A struct containing
%            sigma - The noise hyper parameter
%            Mu    - The mean matrix (usually the identity)
%   method - Can be one of "MVN", "lietal" or "none"
%
% Output:
%   A - MxM matrix
%   b - Mx1 vector
%   c - 1x1 value
%
function [A b c] = QuadricFit( Data, Params, method="MVN" )
  [M N] = size(Data);
  if( N < 6 )
    warning("QuadricFit: The system is underconstrained");
  end
  
  % Create parameters from struct
  switch method
    case "MVN"
      sigma = Params.sigma;
      Mu = Params.Mu;
    case "lietal"
    case "none"
    otherwise
      error("QuadricFit: The method is not yet implemented");
  end
    
  Mp = M*(M+1)/2;
  % Create the data scatter matrix
  S = DataToScatter(Data, M);
  
  switch method
    case "MVN"
      EI = zeros(Mp+M+1, Mp+M+1);
      vvv=MatrixToVector(2*ones(M,M)-eye(M,M));
      EI(1:Mp,1:Mp) = diag(vvv);
      a=pinv((sigma^2)*S + EI)*[MatrixToVector(Mu + Mu' - diag(diag(Mu)));zeros(M+1,1)];
    case "lietal"
      Di = [Data(1,:)'.^2, ...
            Data(2,:)'.^2, ...
            Data(3,:)'.^2, ...
            2*Data(2,:)'.*Data(3,:)', ...
            2*Data(1,:)'.*Data(3,:)', ...
            2*Data(1,:)'.*Data(2,:)', ...
            Data(1,:)', Data(2,:)', Data(3,:)',...
            ones(N,1)];
      D = Di'*Di;
      C = zeros(10,10);
      C(1:3,1:3) = [-1 1 1; ...
                    1 -1 1; ...
                    1 1 -1];
      C(4:6,4:6) = -4*eye(3);
      
      [U S]= eig( ...
            inv(C(1:6,1:6)) * ...
            (D(1:6,1:6) - D(7:10,1:6)'*(inv(D(7:10,7:10))*D(7:10,1:6))) ...
            );
      
      [~,ii] = max(diag(S));
      aa = U(:,ii);
      aa2 = -(D(7:10,7:10)\D(7:10,1:6))*aa;
      aa3 = [aa;aa2];
    case "none"
      [Ui Si] = eig(S);
      [~,ii] = min(diag(abs(Si)));
      a = Ui(:,ii);
  end
 
  
  if(strcmp("lietal", method))
    A = [aa3([1, 6, 5])'; aa3([6, 2, 4])'; aa3([5, 4, 3])'];
    b = [aa3([7,8,9])];
    c = aa3(10);
  else
    [A, b, c] = GetParameters( a, M, Mp );
  end
  
  if(all(eig(A)<=0))
    %warning("QuadricFit: All eigenvalues are less than zero");
    A=-A;
    b=-b;
    c=-c;
  end
endfunction

function [A, b, c] = GetParameters( a, M, Mp );
  A = VectorToMatrix(a(1:Mp),M);
  b = [a(Mp+1:end-1)];
  c = a(end);
endfunction

function err = FindErrors(Data, A, b, c)
  err = sum(Data.*(A*Data)) + b'*Data + c; 
endfunction

function [S,D] = DataToScatter(Data, M)
  [K,N] = size(Data);
  ii = tril(ones(M));
  [x y] = find(ii);
  
  D = [(Data(x,:)'.*Data(y,:)').*repmat(((x~=y) + 1)', N,1), ...
       Data(x(x==y),:)',...
       ones(N,1)];

  S = D' * D;
endfunction

function theta = MatrixToVector(A)
  [M N] = size(A);
  At = A';
  
  assert(At(:),A(:), 1e-5);
  assert(M==N);
  ii = tril(ones(M));

  theta = zeros(M*(N+1)/2,1);
  theta = A(ii(:)>0);
endfunction

function A = VectorToMatrix(theta, M)
  K = numel(theta);
  assert(K == M*(M+1)/2);
  On = tril(ones(M));
  A = zeros(M);
  
  A(On(:)>0) = theta;
  B =A';
  A(On==0) = B(~On);
endfunction
