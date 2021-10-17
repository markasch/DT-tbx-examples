function G = mynmgrid(R,n)
%MYNMGRID Number the grid points in a two dimensional region.
%	G = MYNMGRID('R',n) numbers the points on an n-by-n grid in
%	the subregion of -1<=x<=1 and -1<=y<=1 determined by 'R'.
%	SPY(MYNMGRID('R',n)) plots the points.
%	MYDELSQ(MYNMGRID('R',n)) generates the 5-point discrete Laplacian.
%	The regions are:
%	   'S'  - entire square.
%	   'L'  - L-shaped domain made from 3/4 of the entire square.
%	   'C'  - like the 'L', but with a quarter circle in the 4-th square.
%	   'D'  - unit disc.
%	   'A'  - annulus.
%	   'H'  - heart-shaped cardioid.
%	   'No' - notch.
%    'Sh' - hollow square.
%    'El' - ellipse.
%    'Db' - dumb-bell.
%		 'Ds' - disc with 2 slits.
%    'Ar' - arch.
% Based on: NUMGRID of
%	C. Moler, 7-16-91, 12-22-93.
%	Copyright (c) 1984-94 by The MathWorks, Inc.
x = ones(n,1)*[-1, (-(n-3):2:(n-3))/(n-1), 1];
y = flipud(x');
if R == 'S'
   G = (x > -1) & (x < 1) & (y > -1) & (y < 1);
elseif R == 'Sh'
   G = ((x > -1) & (x < 1) & (y > -1) & (y < 1)) ... 
   & ( ((x > .2) | (x < -.2)) | ((y > .2) | (y < -.2)));
elseif R == 'L'
   G = (x > -1) & (x < 1) & (y > -1) & (y < 1) & ( (x > 0) | (y > 0));
elseif R == 'C'
   G = (x > -1) & (x < 1) & (y > -1) & (y < 1) & ((x+1).^2+(y+1).^2 > 1);
elseif R == 'No'
   G = (x > -1) & (x < 1) & (y > -1) & (y < 1) & ...
       ( (x<=-.5 | x>=.5) | (y<0));
elseif R == 'D'
   G = x.^2 + y.^2 < (1);%-4/(n-1)); % so that bdy is really at r=1
elseif R == 'Db'
   x0=0.65; r=0.25; y0=0;
   G = (x-x0).^2 + (y-y0).^2 < r^2 | ...
         (x+x0).^2 + (y+y0).^2 < r^2 | ...
         ((x < 0.6) & (x > -0.6) & ...
         (y < 0.1) & (y > -0.1)); 
elseif R == 'Ar'
   a=1;b=.5; h=0.; xa=-.5+h; xb=-.2; ya=0; yb=-1;
   m=(yb-ya)/(xb-xa); c=ya-m*xa;
   G =  (((x/a).^2 + (y/b).^2 < 1) & (y > 0)) ...
         | ( (x> -1 & x < -0.55) & (y> -0.75 &y < 0))...
         | ( (x> 0.55 & x < 1) & (y> -0.75 &y < 0) );
elseif R == 'El'
   a=1;b=.5;
   G =  (x/a).^2 + (y/b).^2 < 1;
elseif R == 'A'
   G = ( x.^2 + y.^2 < 1) & ( x.^2 + y.^2 > 1/3);
elseif R == 'H'
   RHO = .75; SIGMA = .75;
   G = (x.^2+y.^2).*(x.^2+y.^2-SIGMA*y) < RHO*x.^2;
elseif R == 'Ds'
   h = 0.05; r2=x.^2 + y.^2;
   G = ((r2 < 1) & (y>=h | y<=-h)) | ( r2 <= 0.125) ;%& (0.5<x<1 | -1<x<-.5));
else
      error('Invalid region type.');
end
k = find(G);
G = zeros(size(G));      % Convert from logical to double matrix
G(k) = (1:length(k))';