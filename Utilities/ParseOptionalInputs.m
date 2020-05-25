function opts = ParseOptionalInputs(opts,myvarargin)
%Process optional inputs
if mod(length(myvarargin), 2) ~= 0, error('Must pass key/value pairs for options.'); end
for i = 1:2:length(myvarargin)
    try
        opts.(myvarargin{i}) = myvarargin{i+1};
    catch
        error('Couldn''t set option ''%s''.', myvarargin{2*i-1});
    end
end
end %function