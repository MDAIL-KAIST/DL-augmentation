%%  bin %%

%% Wrapper for n-d binning
%%  obj        - object to bin
%%  bin_factor - integer factor to bin by
%%  n          - (optional) binning dimension. n=2 will perform 2D binning and
            %% n=3 will perform 3D. 2D binning across a 3D stack, such as binning of projections
            %% can be accomplished by passing an N x N x num_projections array to obj and 
            %% setting n=2

%%outputs:
%%  out      - binned object

%% Author: Alan (AJ) Pryor, Jr.
%% Jianwei (John) Miao Coherent Imaging Group
%% University of California, Los Angeles
%% Copyright (c) 2015-2016. All Rights Reserved.

function out = bin(obj,bin_factor,n)
if nargin<3
    if ndims(obj)==2
        out = bin2(obj,bin_factor);
    elseif ndims(obj)==3
        out = bin3(obj,bin_factor);
    else
        error('Only 2 and 3 dimensional binning allowed')
    end
else
    switch n
        case 2
            out = bin2(obj,bin_factor);
        case 3
            out = bin3(obj,bin_factor);
    end
end

end

function out = bin2(obj,bin_factor)
obj = convn(obj,ones(bin_factor,bin_factor)./bin_factor^2,'same');
out = obj(1:bin_factor:end,1:bin_factor:end,:);
end


function out = bin3(obj,bin_factor)
obj = convn(obj,ones(bin_factor,bin_factor,bin_factor)./bin_factor^3,'same');
out = obj(1:bin_factor:end,1:bin_factor:end,1:bin_factor:end);
end
