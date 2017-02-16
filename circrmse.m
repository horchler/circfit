function err=circrmse(x,y,r,xc,yc)
%RMSE  Root mean squared error for a circle
%   ERR = CIRCRMSE(X,Y,R,XC,YC) returns the scalar root mean squared error of a
%   circle radius and center relative to position data. X and Y are equal
%   length 1-D arrays of position data in a rectilinear coordinate system. R,
%   XC, and YC are the scalar radius, X and Y center positions of the circle,
%   respectively.
%
%   ERR = CIRCRMSE(X,Y,R) assumes XC and YC are equal to zero.
%
%   See also CIRCFIT, PLOTCIRCFIT

%   Andrew D. Horchler, horchler @ gmail . com, Created 5-12-7
%   Revision: 1.2, 4-6-16


% Check inputs
if ~isvector(x) || ~isfloat(x) || ~isreal(x) || ~all(isfinite(x))
    error('circfit:circrmse:NonFiniteRealVectorX',...
          'X must be a finite real vector of floating point numbers.');
end
if ~isvector(y) || ~isfloat(y) || ~isreal(y) || ~all(isfinite(y))
    error('circfit:circrmse:NonFiniteRealVectorY',...
          'Y must be a finite real vector of floating point numbers.');
end
lx=length(x);
if lx ~= length(y)
    error('circfit:circrmse:LengthMismatch',...
          'The vectors X and Y must have the same length.');
end
if lx < 3
    error('circfit:circrmse:Min3Points',...
          'The vectors X and Y must contain at least three points.');
end

if ~isscalar(r) || ~isfloat(r) || ~isreal(r) || ~isfinite(r)
    error('circfit:circrmse:NonFiniteRealScalarR',...
          'R must be a finite real scalar.');
end
if r < 0
    error('circfit:circrmse:NegativeR','R must be a positive value.');
end

if nargin == 3
    xc=0;
    yc=0;
elseif nargin == 5
    if ~isscalar(xc) || ~isfloat(xc) || ~isreal(xc) || ~isfinite(xc)
        error('circfit:circrmse:NonFiniteRealScalarXC',...
              'XC must be a finite real scalar.');
    end
    if ~isscalar(yc) || ~isfloat(yc) || ~isreal(yc) || ~isfinite(yc)
        error('circfit:circrmse:NonFiniteRealScalarYC',...
              'YC must be a finite real scalar.');
    end
elseif nargin == 4
    error('circfit:circrmse:InvalidInputNumber',...
          'Either both XC and YC must be specified or neither.');
else
    error('circfit:circrmse:TooManyInputs','Too many input arguments.');
end

x=x(:);
y=y(:);
if lx < 20 && rank(diff([x y;x(1) y(1)])) == 1
    error('circfit:circrmse:Collinearity',...
         ['The points in vectors X and Y must not all be collinear, or '...
          'nearly collinear, with each other.']);
end

% Calculate RMSE
err=sqrt(mean((sqrt((x-xc).^2+(y-yc).^2)-r).^2));