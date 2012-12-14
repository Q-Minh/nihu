function args = create_slab_args(varargin)

switch nargin
    case 1 % One argument mode (Cx, Cx mode)
        Cx = sort(varargin{1});
        Cy = Cx;
        N = length(Cx) * [1 1] - 1;
    case 2 % Two argument mode
        % Cx, Cy mode
        if (size(varargin{1},1) > 2 && size(varargin{1},2) == 1) || size(varargin{2},1) > 2
            Cx = sort(varargin{1});
            Cy = sort(varargin{2});
            N = [length(Cx)-1 length(Cy)-1];
        % R, N mode
        else
            % Process corners
            R = varargin{1};
            switch size(R,1)
                case 1
                    if isscalar(R)
                        R = [R R];
                    end
                    R = [
                        0    0    0
                        R(1) 0    0
                        R(1) R(2) 0
                        0    R(2) 0
                        ];
                case 2
                    R = [
                        R(1,1) R(1,2) 0
                        R(2,1) R(1,2) 0
                        R(2,1) R(2,2) 0
                        R(1,1) R(2,2) 0
                        ];
                case 4
                    R = [R zeros(size(R,1), 3-size(R,2))];
                otherwise
                    error('NiHu:create_slab:argFormat',...
                        'Unsupported format of input arguments.');
            end
            N = varargin{2};
            if isscalar(N)
                 N = [N N];
            end
        end
    otherwise
        error('NiHu:create_slab:argNumber',...
            'Unsupported number of arguments: %d.', nargin);
end

args.N = N;
if exist('R', 'var')
    args.R = R;
else
    args.Cx = Cx;
    args.Cy = Cy;
end

end

