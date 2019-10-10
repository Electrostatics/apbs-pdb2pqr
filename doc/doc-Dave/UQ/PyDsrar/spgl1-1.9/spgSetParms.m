function options = spgSetParms(varargin)
%SPGSETPARMS  Set options for SPGL1
%
%   options = spgSetParms('param1',val1,'param2',val2,...) creates an
%   options structure in which the named parameters have the specified
%   values.  Unspecified parameters are empty and their default
%   values are used.
%   
%   spgSetParms with no input arguments and no output arguments
%   displays all parameter names and their possible values.
%
%   options = spgSetParms (with no input arguments) creates an options
%   structure where all the fields are empty.
%
%   spgSetParms.m
%   $Id: spgSetParms.m 1407 2009-06-30 20:00:54Z ewout78 $

% Print out possible values of properties.
if nargin == 0 && nargout == 0
   fprintf(' Default parameters for l1Set.m:\n');
   fprintf('        fid : [ positive integer        |     1 ]\n');
   fprintf('  verbosity : [ integer: 1, 2, or 3     |     3 ]\n');
   fprintf(' iterations : [ positive integer        |  10*m ]\n');
   fprintf('  nPrevVals : [ positive integer        |    10 ]\n');
   fprintf('      bpTol : [ positive scalar         | 1e-06 ]\n');
   fprintf('      lsTol : [ positive scalar         | 1e-06 ]\n');
   fprintf('     optTol : [ positive scalar         | 1e-04 ]\n');
   fprintf('     decTol : [ positive scalar         | 1e-04 ]\n');   
   fprintf('    stepMin : [ positive scalar         | 1e-16 ]\n');
   fprintf('    stepMax : [ positive scalar         | 1e+05 ]\n');
   fprintf(' rootMethod : [ 1=linear, 2=quadratic   |     2 ]\n');
   fprintf('activeSetIt : [ positive integer        |   Inf ]\n');
   fprintf('subspaceMin : [ 0=no, 1=yes             |     0 ]\n');
   fprintf('  iscomplex : [ 0=no, 1=yes, NaN=auto   |   NaN ]\n');
   fprintf('  maxMatvec : [ positive integer        |   Inf ]\n');
   fprintf('    weights : [ vector                  |     1 ]\n');
   fprintf('    project : [ projection function     |    @()]\n');
   fprintf('primal_norm : [ primal norm eval fun    |    @()]\n');
   fprintf('  dual_norm : [ dual norm eval fun      |    @()]\n');
   fprintf('\n');
   return;
end

Names = [
    'fid               '
    'verbosity         '
    'iterations        '
    'nPrevVals         '
    'bpTol             '
    'lsTol             '
    'optTol            '
    'decTol            '
    'stepMin           '
    'stepMax           '
    'rootMethod        '
    'activeSetIt       '
    'subspaceMin       '
    'iscomplex         '
    'maxMatvec         '
    'weights           '
    'project           '
    'primal_norm       '
    'dual_norm         '
	];
[m,n] = size(Names);
names = lower(Names);

% Combine all leading options structures o1, o2, ... in l1Set(o1,o2,...).
options = [];
for j = 1:m
   eval(['options.' Names(j,:) '= [];']);
end
i = 1;
while i <= nargin
   arg = varargin{i};
   if ischar(arg), break; end
   if ~isempty(arg)                      % [] is a valid options argument
       if ~isa(arg,'struct')
          error(sprintf(['Expected argument %d to be a string parameter name ' ...
               'or an options structure\ncreated with OPTIMSET.'], i));
      end
      for j = 1:m
          if any(strcmp(fieldnames(arg),deblank(Names(j,:))))
             eval(['val = arg.' Names(j,:) ';']);
          else
             val = [];
          end
          if ~isempty(val)
             eval(['options.' Names(j,:) '= val;']);
         end
      end
   end
   i = i + 1;
end

% A finite state machine to parse name-value pairs.
if rem(nargin-i+1,2) ~= 0
   error('Arguments must occur in name-value pairs.');
end
expectval = 0;                          % start expecting a name, not a value
while i <= nargin
   arg = varargin{i};
   
   if ~expectval
      if ~ischar(arg)
         error(sprintf('Expected argument %d to be a string parameter name.', i));
      end
      
      lowArg = lower(arg);
      j = strmatch(lowArg,names);
      if isempty(j)                       % if no matches
         error(sprintf('Unrecognized parameter name ''%s''.', arg));
      elseif length(j) > 1                % if more than one match
         % Check for any exact matches (in case any names are subsets of others)
         k = strmatch(lowArg,names,'exact');
         if length(k) == 1
            j = k;
         else
            msg = sprintf('Ambiguous parameter name ''%s'' ', arg);
            msg = [msg '(' deblank(Names(j(1),:))];
            for k = j(2:length(j))'
               msg = [msg ', ' deblank(Names(k,:))];
            end
            msg = sprintf('%s).', msg);
            error(msg);
         end
      end
      expectval = 1;                      % we expect a value next
      
   else
      eval(['options.' Names(j,:) '= arg;']);
      expectval = 0;
      
   end
   i = i + 1;
end

if expectval
   error(sprintf('Expected value for parameter ''%s''.', arg));
end

