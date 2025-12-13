function test = isdmatrix(x,varargin)

% Check number of parameters
if nargin < 1,
   error('Incorrect number of parameters (type ''help <a href="matlab:help isdmatrix">isdmatrix</a>'' for details).');
 end
 
 % Test: doubles, two dimensions, two or more columns?
 test = isa(x,'double') & length(size(x)) == 2 & size(x,2) >= 2;
 
 % Optional tests
 for i = 1:length(varargin),
     try
         if varargin{i}(1) == '#',
             if size(x,1) ~= str2num(varargin{i}(2:end)), test = false; return; end
         elseif varargin{i}(1) == '@',
             if size(x,2) ~= str2num(varargin{i}(2:end)), test = false; return; end
         elseif ~eval(['all(x(~isnan(x))' varargin{i} ');']), test = false; return; end
     catch err
         error(['Incorrect test ''' varargin{i} ''' (type ''help <a href="matlab:help isdmatrix">isdmatrix</a>'' for details).']);
     end
 end