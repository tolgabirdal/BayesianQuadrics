%
% 
%

clear all;
close all;

pkg load statistics;

mkdir Output

[A,b,c] = GenerateQuadric(0.01*eye(3),20, 0.01*eye(3), 0.01);
 
[f, v] = PlotQuadricFromUnitCircle(A, b, c);
WritePly(v',[],uint32([f(1,:); f(3,:); f(2,:)])', 'Output/Original.ply');

sigma = 0.01;
vtrain = v + sigma*randn(size(v));
Data = vtrain(:,1:10:end,:);
WritePly(Data',80*ones(size(Data')), [], 'Output/Samples.ply');

for prior = 1 : 10 : 100
  [Afit bfit cfit] = QuadricFit( Data, struct("sigma",prior, "Mu", 20*eye(3)), "MVN" );
  
  [f, v] =  PlotQuadricFromUnitCircle(Afit, bfit, cfit);
  
  WritePly(v',[],uint32([f(1,:); f(3,:); f(2,:)])', sprintf('Output/Fitted%d.ply',prior));
end

[Afit bfit cfit] = QuadricFit( Data, struct(), "none" );

mi = min(Data,[],2);
ma = max(Data,[],2);
[f, v, c] = PlotQuadratic( Afit, bfit, cfit, linspace(mi(1),ma(1),40), linspace(mi(2),ma(2),40), linspace(mi(3), ma(3),40));

WritePly(v,[],uint32([f(:,1), f(:,2), f(:,3)])-1, sprintf('Output/NoPrior.ply',prior));

