function xi = invmap134b(xvert, x0)

x0 = x0(:);

% Distance toleration constant
dtol = 1.8;
% Check if inside
m = x0.';
midn = sum(xvert(1:4,:),1)/4;   % near midpoint
midf = sum(xvert(5:8,:),1)/4;   % far midpoint
dstm = norm(midn-midf);         % distances of midpoints
vecm = (midf-midn)/dstm;        % unit direction vector of midline
% Distances from midpoints
devn = max(abs([midn;midn;midn;midn]-xvert(1:4,:)),[],1);
devf = max(abs([midf;midf;midf;midf]-xvert(5:8,:)),[],1);
% Distance from midline
dstl = (m-midn)*vecm';
% Check if in the cone
if(norm(m-(midn+dstl*vecm)) > (devn + (devf - devn) * dstl/dstm)*dtol)
    xi = [inf inf inf];
    return
end

maxiter = 50;
converged = false;
tol = 1e-5;

xivec = [
    0 0 0
    -1 -1 -1
    0 -1 -1
    0 1 -1
    -1 1 -1
    -1 -1 1
    0 -1 1
    0 1 1
    -1 1 1
    ];

for p = 1 : 9
    xi = xivec(p,:);
    iIter = 0;
    while true
        iIter = iIter+1;
        [N, dN] = shape(xi);
        %         [N, dN, ddN] = shape(xi);
        x = xvert' * N;
        d = x - x0;
        f = d' * d;
        if f < tol^2
            converged = true;
            break;
        end
        dx = xvert' * dN;
        df = (2 * d' * dx)';
        ddf = dx'* dx;
        %         for k = 1 : 8
        %             s = ddN(:,:,k);
        %             for i = 1 : 3
        %                 ddf = ddf + d(i) * xvert(k,i) * s;
        %             end
        %         end
        ddf = 2*ddf;
        dxi = (ddf \ df);
        if(xi(1) - dxi(1) < 1)
            xi = (xi' - dxi)';
        else
            xi = (xi' - 1/2*(1-abs(xi(1)))*dxi)';
        end
        if iIter > maxiter
            break;
        end
    end
    if converged
        break;
    end
end
 if ~converged
     xi = [inf inf inf];
     %disp('nihu');
 end

end

function [N, dN] = shape(xivec)
s = xivec(1);
t = xivec(2);
u = xivec(3);

N = [
     1/2*(1-u).*(1-t).*s./(s-1);
     1/2*(1-u).*(1+t).*s./(s-1);
     1/2*(1+u).*(1+t).*s./(s-1);
     1/2*(1+u).*(1-t).*s./(s-1);
    -1/2*(1-u).*(1-t).*(1+s)./(2*(s-1));
    -1/2*(1-u).*(1+t).*(1+s)./(2*(s-1));
    -1/2*(1+u).*(1+t).*(1+s)./(2*(s-1));
    -1/2*(1+u).*(1-t).*(1+s)./(2*(s-1));
    ];

dN = [
    -1/2*(1-u).*(1-t)./(s-1).^2, -1/2*(1-u).*s./(s-1), -1/2*(1-t).*s./(s-1);
    -1/2*(1-u).*(1+t)./(s-1).^2,  1/2*(1-u).*s./(s-1), -1/2*(1+t).*s./(s-1);
    -1/2*(1+u).*(1+t)./(s-1).^2,  1/2*(1+u).*s./(s-1),  1/2*(1+t).*s./(s-1);
    -1/2*(1+u).*(1-t)./(s-1).^2, -1/2*(1+u).*s./(s-1),  1/2*(1-t).*s./(s-1);
   
     1/2*(1-u).*(1-t)./(s-1).^2,  1/4*(1-u).*(1+s)./(s-1),  1/4*(1-t).*(1+s)./(s-1);
     1/2*(1-u).*(1+t)./(s-1).^2, -1/4*(1-u).*(1+s)./(s-1),  1/4*(1+t).*(1+s)./(s-1);
     1/2*(1+u).*(1+t)./(s-1).^2, -1/2*(1+u).*(1+s)./(s-1), -1/2*(1+t).*(1+s)./(s-1);
     1/2*(1+u).*(1-t)./(s-1).^2,  1/2*(1+u).*(1+s)./(s-1), -1/2*(1-t).*(1+s)./(s-1);
    ];

end