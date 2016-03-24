%
% Test the quadric fitting algorithm. 
%
% This test will run the matlb implementation with 10 different values for the 
% prior. It tests the Li ang Griffiths quadric fitting and also the c++ implementation.
%
% The stanford ply files which are generated as output can be opened using meshlab.
%

clear all;
close all;

pkg load statistics;

mkdir Output


% Generate an arbitrary quadric and write to file
[A,b,c] = GenerateQuadric(0.01*eye(3),20, 0.01*eye(3), 0.01);
 
[f, v] = PlotQuadricFromUnitCircle(A, b, c);
WritePly(v',[],uint32([f(1,:); f(3,:); f(2,:)])', 'Output/Original.ply');

% Generate some points on the surface and add gaussian noise
sigma = 0.01;
vtrain = v + sigma*randn(size(v));
Data = vtrain(:,1:10:end,:);
WritePly(Data',80*ones(size(Data')), [], 'Output/Samples.ply');

% For values of sigma between 1 and 100, fit a quadric.
for prior = 1 : 10 : 100
  Params = struct;
  Params.sigma = prior;
  Params.Mu = 20*eye(3);
  
  [Afit bfit cfit] = QuadricFit( Data, Params, "MVN" );
  
  [f, v] =  PlotQuadricFromUnitCircle(Afit, bfit, cfit);
  
  WritePly(v',[],uint32([f(1,:); f(3,:); f(2,:)])', sprintf('Output/Fitted%d.ply',prior));
end

% Fit a quadric to the data without a prior
[Afit bfit cfit] = QuadricFit( Data, struct(), "none" );
mi = min(Data,[],2);
ma = max(Data,[],2);
[f, v, c] = PlotQuadratic( Afit, bfit, cfit, ...
                           linspace(mi(1),ma(1),40), ...
                           linspace(mi(2),ma(2),40), ...
                           linspace(mi(3), ma(3),40));
WritePly(v,[],uint32([f(:,1), f(:,2), f(:,3)])-1, sprintf('Output/NoPrior.ply',prior));

% Fit a quadric using the c++ implementation
[A b c] = Quadric( Data, 10 );
[f, v] =  PlotQuadricFromUnitCircle(A, b, c);
WritePly(v',[],uint32([f(1,:); f(3,:); f(2,:)])', 'Output/cQuad.ply');
