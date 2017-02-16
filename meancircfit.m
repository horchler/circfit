function r=meancircfit(x,y,w)
%MEANRADIUS  Fit portions of data to circles and find mean of radii
%   R = MEANRADIUS(X,Y,W) calls CIRCFIT in a loop and returns the scalar mean
%   radius. X and Y are equal length 1-D arrays of position data in a
%   rectilinear coordinate system. W is a scalar that sets the windwow width,
%   the number of points used for each circle fit. W must be greater than or
%   equal to two so that CIRCFIT has at least three points to fit.
%
%   This function can be helpul in cases where a systematic drift is present in
%   the data (e.g., a slowly spiraling trajectory) and one wishes to estimate
%   just the mean local radius of curvature.
%
%   See also CIRCFIT

%   Andrew D. Horchler, horchler @ gmail . com, Created 5-12-7
%   Revision: 1.2, 4-6-16


% Check inputs
if ~isvector(x) || ~isfloat(x) || ~isreal(x) || ~all(isfinite(x))
    error('circfit:meancircfit:NonFiniteRealVectorX',...
          'X must be a finite real vector of floating point numbers.');
end
if ~isvector(y) || ~isfloat(y) || ~isreal(y) || ~all(isfinite(y))
    error('circfit:meancircfit:NonFiniteRealVectorY',...
          'Y must be a finite real vector of floating point numbers.');
end
lx=length(x);
if lx ~= length(y)
    error('circfit:meancircfit:LengthMismatch',...
          'The vectors X and Y must have the same length.');
end
if lx < 3
    error('circfit:meancircfit:Min3Points',...
          'The vectors X and Y must contain at least three points.');
end

if ~isscalar(w) || ~isnumeric(w) || ~isreal(w) || ~isfinite(w)
    error('circfit:meancircfit:NonFiniteRealScalarW',...
          'W must be a finite real integer.');
end
if w < 2 || w-floor(w) ~= 0
    error('circfit:meancircfit:InvalidIntegerW',...
          'W must be a finite real integer greater than or equal to two.');
end

rw=zeros(1,lx-w);
hw=floor(0.5*w);
for i=hw+1:lx-hw
    rw(i)=circfit(x(i-hw:i+hw),y(i-hw:i+hw));
end
r=mean(rw);