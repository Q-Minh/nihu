function [xx, ww] = gaussquad2(order, nvert)
%GAUSSQUAD2 2D Gaussian quadrature integration
%   [X, W] = GAUSSQUAD2(P, Nvert) returns Gaussian integration base
%   points and weights for numerical integration over stadard 2D elements.
% Parameters:
%   P     : Order of quadrature
%   Nvert : Number of vertices. 4 for standard QUAD, 3 for standard TRIA
%           elements. The standard QUAD's coordinates are
%           [(-1,-1), (+1,-1), (+1,+1), (-1,+1)]
%           The standard TRIA element's coordinates are
%           [(0,0), (+1,0), (0,+1)]
%   X     : Gaussian quadrature base points
%   W     : Gaussian quadrature weights

%   Copyright 2008-2010 P. Fiala
%   Budapest University of Technology and Economics
%   Dept. of Telecommunications

% Last modified 2015.03.05.
%% Argument check
narginchk(1, 2);
if nargin < 2
    nvert = 4;
end

%%
switch nvert
    case 4
        [xi, w] = gaussquad1(order); % Linear Gauss quad over (-1, +1)
        [eta, xi] = meshgrid(xi,xi);  % base points over unit rectangle
        w = (w * w.');                % weights over rectangle
        xx = [xi(:) eta(:)];
        ww = w(:);
    case 3
        switch order
            case {0, 1}
                xx = [1/3 1/3];
                ww = 1/2;
            case 2
                xx = [
                    0.166666666666667   0.666666666666667
                    0.166666666666667   0.166666666666667
                    0.666666666666667   0.166666666666667
                    ];
                ww = [
                    0.166666666666666
                    0.166666666666666
                    0.166666666666666
                    ];
            case 3
                xx = [
                    0.333333333333333   0.333333333333333
                    0.200000000000000   0.600000000000000
                    0.200000000000000   0.200000000000000
                    0.600000000000000   0.200000000000000
                    ];
                ww = [
                    -0.281250000000000
                    0.260416666666667
                    0.260416666666667
                    0.260416666666667
                    ];
            case 4
                xx = [
                    0.445948490915965   0.108103018168070
                    0.445948490915965   0.445948490915965
                    0.108103018168070   0.445948490915965
                    0.091576213509771   0.816847572980459
                    0.091576213509771   0.091576213509771
                    0.816847572980459   0.091576213509771
                    ];
                ww = [
                    0.111690794839006
                    0.111690794839006
                    0.111690794839006
                    0.054975871827661
                    0.054975871827661
                    0.054975871827661
                    ];
            case 5
                xx = [
                    0.333333333333333   0.333333333333333
                    0.470142064105115   0.059715871789770
                    0.470142064105115   0.470142064105115
                    0.059715871789770   0.470142064105115
                    0.101286507323456   0.797426985353087
                    0.101286507323456   0.101286507323456
                    0.797426985353087   0.101286507323456
                    ];
                ww = [
                    0.112500000000000
                    0.066197076394253
                    0.066197076394253
                    0.066197076394253
                    0.062969590272414
                    0.062969590272414
                    0.062969590272414
                    ];
            case 6
                xx = [
                    0.249286745170910   0.501426509658179
                    0.249286745170910   0.249286745170910
                    0.501426509658179   0.249286745170910
                    0.063089014491502   0.873821971016996
                    0.063089014491502   0.063089014491502
                    0.873821971016996   0.063089014491502
                    0.310352451033784   0.053145049844817
                    0.636502499121399   0.310352451033784
                    0.053145049844817   0.636502499121399
                    0.053145049844817   0.310352451033784
                    0.310352451033784   0.636502499121399
                    0.636502499121399   0.053145049844817
                    ];
                ww = [
                    0.058393137863189
                    0.058393137863189
                    0.058393137863189
                    0.025422453185104
                    0.025422453185104
                    0.025422453185104
                    0.041425537809187
                    0.041425537809187
                    0.041425537809187
                    0.041425537809187
                    0.041425537809187
                    0.041425537809187
                    ];
            case 7
                xx = [
                    0.333333333333333   0.333333333333333
                    0.260345966079040   0.479308067841920
                    0.260345966079040   0.260345966079040
                    0.479308067841920   0.260345966079040
                    0.065130102902216   0.869739794195568
                    0.065130102902216   0.065130102902216
                    0.869739794195568   0.065130102902216
                    0.312865496004874   0.048690315425316
                    0.638444188569810   0.312865496004874
                    0.048690315425316   0.638444188569810
                    0.048690315425316   0.312865496004874
                    0.312865496004874   0.638444188569810
                    0.638444188569810   0.048690315425316
                    ];
                ww = [
                    -0.074785022233841
                    0.087807628716604
                    0.087807628716604
                    0.087807628716604
                    0.026673617804419
                    0.026673617804419
                    0.026673617804419
                    0.038556880445128
                    0.038556880445128
                    0.038556880445128
                    0.038556880445128
                    0.038556880445128
                    0.038556880445128
                    ];
            case 8
                xx = [
                    0.333333333333333   0.333333333333333
                    0.459292588292723   0.081414823414554
                    0.459292588292723   0.459292588292723
                    0.081414823414554   0.459292588292723
                    0.170569307751760   0.658861384496480
                    0.170569307751760   0.170569307751760
                    0.658861384496480   0.170569307751760
                    0.050547228317031   0.898905543365938
                    0.050547228317031   0.050547228317031
                    0.898905543365938   0.050547228317031
                    0.263112829634638   0.008394777409958
                    0.728492392955404   0.263112829634638
                    0.008394777409958   0.728492392955404
                    0.008394777409958   0.263112829634638
                    0.263112829634638   0.728492392955404
                    0.728492392955404   0.008394777409958
                    ];
                ww = [
                    0.072157803838894
                    0.047545817133642
                    0.047545817133642
                    0.047545817133642
                    0.051608685267359
                    0.051608685267359
                    0.051608685267359
                    0.016229248811599
                    0.016229248811599
                    0.016229248811599
                    0.013615157087218
                    0.013615157087218
                    0.013615157087218
                    0.013615157087218
                    0.013615157087218
                    0.013615157087218
                    ];
            case 9
                xx = [
                    0.333333333333333   0.333333333333333
                    0.489682519198738   0.020634961602525
                    0.489682519198738   0.489682519198738
                    0.020634961602525   0.489682519198738
                    0.437089591492937   0.125820817014127
                    0.437089591492937   0.437089591492937
                    0.125820817014127   0.437089591492937
                    0.188203535619033   0.623592928761935
                    0.188203535619033   0.188203535619033
                    0.623592928761935   0.188203535619033
                    0.044729513394453   0.910540973211095
                    0.044729513394453   0.044729513394453
                    0.910540973211095   0.044729513394453
                    0.221962989160766   0.036838412054736
                    0.741198598784498   0.221962989160766
                    0.036838412054736   0.741198598784498
                    0.036838412054736   0.221962989160766
                    0.221962989160766   0.741198598784498
                    0.741198598784498   0.036838412054736
                    ];
                ww = [
                    0.048567898141400
                    0.015667350113570
                    0.015667350113570
                    0.015667350113570
                    0.038913770502387
                    0.038913770502387
                    0.038913770502387
                    0.039823869463605
                    0.039823869463605
                    0.039823869463605
                    0.012788837829349
                    0.012788837829349
                    0.012788837829349
                    0.021641769688645
                    0.021641769688645
                    0.021641769688645
                    0.021641769688645
                    0.021641769688645
                    0.021641769688645
                    ];
            otherwise
                [xx, ww] = dunavant_rule(order);
                xx = fliplr(xx.');
                ww = ww(:)/2;
        end
    otherwise
        error('NiHu:gaussquad2:argValue',...
            'Input argument Nvert should be 3 for TRIA and or 4 for QUAD');
end
end