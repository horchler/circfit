function bool=iscolinear(varargin)
%ISCOLINEAR  Check colinearity of N-dimesional rectilinear data points.
%   ISCOLINEAR(X,Y) returns logical 1 (true) if the arrays of planar position
%   data, X and Y, are colinear (or nearly colinear) with respect to each other
%   and logical 0 (false) otherwise. The data in X and Y must be finite and
%   real. If X and Y are both vectors, they must be the same length. If X and Y
%   are both arrays, they must have the same dimensions. If X and Y contain two
%   or fewer elements each, logical 1 is always returned.
%
%   ISCOLINEAR(X,Y,Z) returns logical 1 (true) if the arrays of 3-dimensional
%   position data, X, Y, and Z are colinear (or nearly colinear) with respect to
%   each other and logical 0 (false) otherwise. The data in X, Y, and Z must be
%   finite and real. If X, Y, and Z are all vectors, they must be the same
%   length. If X, Y, and Z are all arrays, they must have the same dimensions.
%   If X, Y, and Z contain two or fewer elements each, logical 1 is always
%   returned.
%
%   ISCOLINEAR(V) returns logical 1 (true) if the M-by-N matrix, V, of
%   N-dimensional position data are colinear (or nearly colinear) with respect
%   to each other and logical 0 (false) otherwise. The data in V must be finite
%   and real. The M rows of V correspond to M points in the N-dimensional
%   rectilinear coordinate system. The N columns of V correspond to the
%   coordinates, X1, X2, ... XN. If M <= 2 (two or fewer points) or N <= 1 (one
%   or fewer dimensions), logical 1 is always returned.
%
%   ISCOLINEAR(Z) returns logical 1 (true) if the M element vector or array, Z,
%   of complex coordinates are colinear (or nearly colinear) with respect to
%   each other on the complex plane and logical 0 (false) otherwise. If M <= 2
%   (two or fewer points), logical 1 is always returned.
%
%   Note:
%       ISCOLINEAR relies upon singular-value-decomposition (SVD) via the RANK
%       function to robustly deterimine if data is colinear or not. As a result,
%       floating-point data that is nearly colinear may be labelled as colinear.
%       For example, in the following
%
%           iscolinear(1:3,1:3)
%           iscolinear(1:3,[1 2+2*eps(2) 3])
%           iscolinear(1:3,[1 2+4*eps(2) 3])
%
%       the first two operations return true while the third returns false. In
%       many applications, e.g., those that invlove solving linear systems, this
%       is desirable as nearly colinear data produces the same singularities an
%       numerical inaccuracies as exactly colinear data. The symbolic and
%       variable precision arithmetic capabilities of ISCOLINEAR to more
%       precisely discern colinearity if necessary, e.g., the following now
%       returns false
%
%           iscolinear(1:3,[1 2+sym(2*eps(2)) 3])
%
%   Class support for input X, Y, Z, V:
%       float: double, single
%       int: int8, uint8, int16, uint16, int32, uint32, int64, uint64
%       logical
%       symbolic
%
%   See also RANK, SYM/RANK, SVD, SYM, VPA.

%   Andrew D. Horchler, horchler @ gmail . com, Created 5-12-7
%   Revision: 1.4, 4-6-16


% Handle variable arguments
if nargin == 1
    v=varargin{1};
elseif nargin == 2
    x=varargin{1};
    y=varargin{2};
elseif nargin == 3
    x=varargin{1};
    y=varargin{2};
    z=varargin{3};
else
    if nargin == 0
        error('circfit:iscolinear:TooFewInputs','Too few input arguments.');
    else
        error('circfit:iscolinear:TooManyInputs','Too many input arguments.');
    end
end

% Check inputs
if nargin == 1
    if ~isvector(v) && isreal(v)
        if ndims(v) ~= 2 || ~(isnumeric(v) || islogical(v) || ...
                isa(v,'sym'))   %#ok<ISMAT>
            error('circfit:iscolinear:InvalidDatatypeV',...
                  'V must be a matrix of numeric, logical, or symbolic values.');
        end
        if isnumeric(v) && ~all(isfinite(v(:)))
            error('circfit:iscolinear:InavlidV','V must be finite and real.');
        end
    else
        if ndims(v) ~= 2 || ~(isnumeric(v) || isa(v,'sym'))	%#ok<ISMAT>
            error('circfit:iscolinear:InvalidComplexZ',...
                 ['Z must be a matrix of numeric, logical, or symbolic '...
                  'values.']);
        end
        if isnumeric(v) && ~all(isfinite(v(:)))
            error('circfit:iscolinear:InavlidV','V must be finite.');
        end
        
        v=[real(v(:)) imag(v(:))];
    end
    [m,n]=size(v);
    
    % RANK requires floating-point or symbolic values
    if ~(isfloat(v) || isa(v,'sym'))
        v=double(v);
    end
else
    if ~(isnumeric(x) || islogical(x) || isa(x,'sym'))
        error('circfit:iscolinear:InvalidDatatypeX',...
             ['X must be a vector or array of numeric, logical, or symbolic '...
              'values.']);
    end
    if isnumeric(x) && ~(isreal(x) && all(isfinite(x(:))))
        error('circfit:iscolinear:InavlidX','X must be finite and real.');
    end
    if ~(isnumeric(y) || islogical(y) || isa(y,'sym'))
        error('circfit:iscolinear:InavlidDatatypeY',...
             ['Y must be a vector or array of numeric, logical, or symbolic '...
              'values.']);
    end
    if isnumeric(y) && ~(isreal(y) && all(isfinite(y(:))))
        error('circfit:iscolinear:InavlidY','Y must be finite and real.');
    end
    
    if nargin == 2
        if isvector(x) && ~isvector(y) || ~isvector(x) && isvector(y)
            error('circfit:iscolinear:VectorArrayMismatchXY',...
                 ['X and Y must both be vectors or arrays of equal '...
                  'dimensions.']);
        end
        
        % Check size and dimensions
        if isvector(x)
            m=length(x);
            if m ~= length(y)
                error('circfit:iscolinear:VectorDimensionMismatchXY',...
                      'The vectors X and Y must have the same length.');
            end
        else
            if ndims(x) ~= ndims(y) || ~isequal(size(x),size(y))
                error('circfit:iscolinear:ArrayDimensionMismatchXY',...
                      'The arrays X and Y must have the same dimensions.');
            end
            m=numel(x);
        end
        n=2;
    else
        if isscalar(z) || isempty(z) || ~(isnumeric(z) || islogical(z) ...
                || isa(z,'sym'))
            error('circfit:iscolinear:InavidDatatypeZ',...
                 ['Z must be a non-empty vector or array containing at '...
                  'least two numeric, logical, or symbolic values.']);
        end
        if isnumeric(z) && ~(isreal(z) && all(isfinite(z(:))))
            error('circfit:iscolinear:InavlidZ','Z must be finite and real.');
        end
        if isvector(x) && (~isvector(y) || ~isvector(z)) || ~isvector(x) ...
                && (isvector(y) || isvector(z))
            error('circfit:iscolinear:VectorArrayMismatchXYZ',...
                 ['X, Y, and Z must all be vectors or they must all be '...
                  'arrays of equal dimensions.']);
        end
        
        % Check size and dimensions
        if isvector(x)
            m=length(x);
            if m ~= length(y) || m ~= length(z)
                error('circfit:iscolinear:VectorDimensionMismatchXYZ',...
                      'The vectors X, Y, and Z must have the same length.');
            end
        else
            if ~isequal(ndims(x),ndims(y),ndims(z)) ...
                    || ~isequal(size(x),size(y),size(z))
                error('circfit:iscolinear:ArrayDimensionMismatchXYZ',...
                      'The arrays X, Y, and Z must have the same dimensions.');
            end
            m=numel(x);
        end
        n=3;
    end
    
    % RANK requires floating-point or symbolic values
    if ~(isfloat(x) || isa(x,'sym'))
        x=double(x);
    end
    if ~(isfloat(y) || isa(y,'sym'))
        y=double(y);
    end
    if nargin == 2
        v = [x(:) y(:)];
    else
        if ~(isfloat(z) || isa(z,'sym'))
            z=double(z);
        end
        v = [x(:) y(:) z(:)];
    end
end

% Check collinearity, try a small number of points first
if m <= 2 || n <= 1
    bool=true;
else
    mx=min(m,64);
    if isa(v,'sym')
        bool=(rank(v([2:mx 1],:)-v(1:mx,:)) == 1);
        if bool && m > 64
            bool=(rank(v([2:end 1],:)-v) == 1);
        end
    else
        bool=(rank(diff(v([1:mx 1],:))) == 1);
        if bool && m > 64
            bool=(rank(diff([v;v(1,:)])) == 1);
        end
    end
end