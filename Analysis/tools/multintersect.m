function varargout= multintersect(varargin)
% insert vectors A,B,C,D...
% returns [common values, idx common in A, idx common in B,...]

if length(varargin)==1 && iscell(varargin{1}) % this way you can pass a cell array
    varargin= varargin{:};
end


varargout= cell(1,length(varargin)+1);
indices= cell(1,length(varargin));
[C{1},indices{1},indices{2}] = intersect(varargin{1},varargin{2});

for i=2:length(varargin)-1
    
    [C{i},ic,indices{i+1}] = intersect(C{i-1},varargin{i+1});
    for j=1:i
        indices{j}= indices{j}(ic);
    end

end
varargout{1}=C{end};
varargout(2:end)= indices;

end