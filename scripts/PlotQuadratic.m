%
% Plots an iso surface on to the current figure.
%

function [f, v, c] = PlotQuadratic( A, b, c, x, y, z, colour='k' )
  [M N C] = size(A);
  assert( all([M N size(b,1)]==3), "PlotQuadric: Sizes do not match" );
  [a(:,:,:,1), a(:,:,:,2), a(:,:,:,3)] = meshgrid (x, y, z);

  ss = permute(a,[4 1 2 3]);
  
  Q = ones(length(x), length(y), length(z));
  for k = 1 : C
    Q = Q.*reshape( ...
              sum( ss(:,:).*(A(:,:,k)*ss(:,:)))+ b(:,k)'*ss(:,:) + c(k), ...
              [length(x), length(y), length(z)]);
  end
  [f, v] = isosurface (a(:,:,:,1), a(:,:,:,2), a(:,:,:,3), Q, 0);
  
  Nx = length(x);
  r = g = b = repmat ([1:Nx] / Nx, [Nx, 1, Nx]);
  p = patch ("Faces", f, "Vertices", v, "EdgeColor", "none");
  
  if(numel(colour)==1)
    switch colour
      case 'k'
      case 'r'
        g = b = zeros(size(r));
      case 'b'
        g = r = zeros(size(b));
      case 'g'
        b = r = zeros(size(g));
      otherwise
        warning('PlotQuadratic: Colour option not implemented.');
    end
  else
    r = r*colour(1);
    g = g*colour(2);
    b = b*colour(3);
  end
  
  cdat = isocolors (a(:,:,:,1), a(:,:,:,2), a(:,:,:,3), r, g, b, v);
  set (p, "FaceVertexCData", cdat);
  
  isofinish(p);
endfunction

function [] = isofinish (p)
  set (gca, "PlotBoxAspectRatioMode", "manual", ...
            "PlotBoxAspectRatio", [1 1 1]);
  set (p, "VertexNormals", -get (p,"VertexNormals")); # Revert normals
  set (p, "FaceColor", "interp");
  ## set (p, "FaceLighting", "phong");
  ## light ("Position", [1 1 5]); # Available with JHandles
endfunction
