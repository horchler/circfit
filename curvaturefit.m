function  [k,err]=curvaturefit(x,y)
%CURVATUREFIT  Least squares fit of X-Y data to find absolute curvature.
%   K = CURVATUREFIT(X,Y) returns absolute value of thescalar curvature K of
%   data fitted to a circle. X and Y are 1-D arrays of position data in a
%   rectilinear coordinate system. X and Y must be the same length and must
%   contain at least three non-colinear points in order for a valid solution to
%   be found.
%
%   [K,ERR] = CURVATUREFIT(X,Y) additionally returns the scalar root mean
%   squared error of the fitted curvature relative to the position data.
%
%   Examples:
%       % Fit of just five noisy points
%       x1=[1 0 -1 0 1]+0.05*randn(1,5); y1=[0 1 0 -1 0]+0.05*randn(1,5);
%       k1=curvaturefit(x1,y1)
%
%       % CURVATUREFIT can sometimes perfom poorly with less than 180-degree arc
%       t=0:0.1:pi; lt=length(t);
%       x2=cos(t)+0.04*randn(1,lt); y2=sin(t)+0.04*randn(1,lt);
%       k2_90deg=curvaturefit(x2(1:floor(lt/2)),y2(1:floor(lt/2)))
%       k2_180deg=curvaturefit(x2,y2)

%   Andrew D. Horchler, horchler @ gmail . com, Created 5-12-7
%   Revision: 1.3, 4-6-16


% Check inputs
if ~isvector(x) || ~isfloat(x) || ~isreal(x) || ~all(isfinite(x))
    error('circfit:curvaturefit:NonFiniteRealVectorX',...
          'X must be a finite real vector of floating point numbers.');
end
if ~isvector(y) || ~isfloat(y) || ~isreal(y) || ~all(isfinite(y))
    error('circfit:curvaturefit:NonFiniteRealVectorY',...
          'Y must be a finite real vector of floating point numbers.');
end

lx=length(x);
if lx ~= length(y)
    error('circfit:curvaturefit:LengthMismatch',...
          'The vectors X and Y must have the same length.');
end
if lx < 3
    error('circfit:curvaturefit:Min3Points',...
          'The vectors X and Y must contain at least three points.');
end
x=x(:);
y=y(:);

% Check collinearity, assume with sufficient points, some will be non-collinear
if rank(diff([x([1:min(50,lx) 1]) y([1:min(50,lx) 1])])) == 1
    if lx <= 50 || rank(diff([x y;x(1) y(1)])) == 1
        k = 0;
        if nargout > 1
            err = NaN;
        end
        return;
    end
end

xx=x.*x;
yy=y.*y;
xy=x.*y;
xxyy=xx+yy;
sx=sum(x);
sy=sum(y);
sxx=sum(xx);
syy=sum(yy);
sxy=sum(xy);

% Solve linear system without inverting
% a=[sx sy lx;sxy syy sy;sxx sxy sx]\[sxx+syy;sum(xxyy.*y);sum(xxyy.*x)];
[L,U]=lu([sx sy lx;sxy syy sy;sxx sxy sx]);
a=U\(L\[sxx+syy;sum(xxyy.*y);sum(xxyy.*x)]);

xc=0.5*a(1);                % X-position of center of fitted circle
yc=0.5*a(2);                % Y-position of center of fitted circle
k=1/sqrt(xc^2+yc^2+a(3));	% Curvature of fitted circle

% RMSE
if nargout > 1
    err=sqrt(mean((1./sqrt((x-xc).^2+(y-yc).^2)-k).^2));
end