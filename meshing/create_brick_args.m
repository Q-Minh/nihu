function args = create_brick_args(varargin)

switch nargin
    case 1 % Cx (= Cy = Cz) mode
        Cx = sort(varargin{1});
        Cy = Cx;
        Cz = Cx;
        N = length(Cx) * [1 1 1] - 1;
    case 2 % R N mode
        % Process corners
        R = varargin{1};
        switch size(R,1)
            case 1
                if isscalar(R)
                    R = [R R R];
                end
                R = [
                    0    0    0
                    R(1) 0    0
                    R(1) R(2) 0
                    0    R(2) 0
                    0    0    R(3)
                    R(1) 0    R(3)
                    R(1) R(2) R(3)
                    0    R(2) R(3)
                    ];
            case 2
                R = [
                    R(1,1) R(1,2) R(1,3)
                    R(2,1) R(1,2) R(1,3)
                    R(2,1) R(2,2) R(1,3)
                    R(1,1) R(2,2) R(1,3)
                    R(1,1) R(1,2) R(2,3)
                    R(2,1) R(1,2) R(2,3)
                    R(2,1) R(2,2) R(2,3)
                    R(1,1) R(2,2) R(2,3)
                    ];
            case 8
            otherwise
                error('NiHu:create_brick:argFormat',...
                    'Unsupported format of input arguments.');
        end
        N = varargin{2};
        if isscalar(N)
            N = [N N N];
        end
    case 3 % Cx Cy Cz mode
        Cx = sort(varargin{1});
        Cy = sort(varargin{2});
        Cz = sort(varargin{3});
        N = [length(Cx) length(Cy) length(Cz)] - 1;
    otherwise
        error('NiHu:create_brick:argNumber',...
            'Unsupported number of arguments: %d.', nargin);
end

args.N = N;
if exist('Cx', 'var')
    args.Cx = Cx;
    args.Cy = Cy;
    args.Cz = Cz;
else
    args.R = R;
end

end
