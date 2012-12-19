function T = rotation_matrix(alpha, v)
%ROTATION_MATRIX Matrix of rotation around a vector
%   T = ROTATION_MATRIX(ALPHA, V) generates the 3x3 matrix of rotation
%   around vector V by an angle ALPHA. The matrix can be applied as
%   x2 = x * T

%   Copyright 2012 P. Fiala
%   Budapest University of Technology and Economics
%   Dept. of Telecommunications

% Last modifed: 2012.12.14.


% convert vector to spherical coordinates
[phi,theta] = cart2sph(v(1),v(2),v(3)); % matlab notation
theta = pi/2-theta; % conventional math notation

% rotation components
Tz = [ % rotate v into xz plane
    cos(phi) -sin(phi) 0
    sin(phi) cos(phi) 0
    0 0 1
    ];
Ty = [ % rotate v into z axis
    cos(theta) 0 sin(theta)
    0 1 0
    -sin(theta) 0 cos(theta)
    ];

Tmain = [ % rotate around z axis
    cos(alpha) -sin(alpha) 0
    sin(alpha) cos(alpha) 0
    0 0 1
    ];

T = Tz * Ty * Tmain' * Ty' * Tz';

end
