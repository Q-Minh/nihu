clear;

corners = [
    0 0 0
    1 0 0
    0 1 0
    ];
lset = ShapeSet.LinearTria;

thetavec = pi * [1 .75 .5 .25] + pi/180;
kernelvec = {'SLP', 'DLP', 'DLPt', 'HSP'};
rvec = [5e-1 2e-1 5e-2 2e-2 5e-3];
rvec = [2e-1 2e-2 2e-3];
ordervec = 5 : 4 : 80;

g = nan(length(thetavec), length(rvec), length(kernelvec), length(ordervec));
gs = nan(length(thetavec), length(rvec), length(kernelvec), length(ordervec));
gst = nan(length(thetavec), length(rvec), length(kernelvec), length(ordervec));

for iT = 1 : length(thetavec)
    theta = thetavec(iT);
    disp(iT);
    for iX = 1 : length(rvec)
        x = [rvec(iX)*cos(theta) .3 rvec(iX)*sin(theta)];
        nx = [sin(theta), 0, -cos(theta)];
        for iK = 1 : length(kernelvec)
            kernel_str = kernelvec{iK};
            for iOrd = 1 : length(ordervec)
                order = ordervec(iOrd);
                
                g(iT, iX, iK, iOrd) =...
                    nearly_singular_laplace_3d_planar_constant(...
                    kernel_str, x, nx, corners, order);
                
                if strcmp(kernel_str, 'HSP')
                    gs(iT, iX, iK, iOrd) =...
                        nearly_singular_laplace_3d_constant_stokes(...
                        x, nx, corners, lset, order, false);
                    
                    gst(iT, iX, iK, iOrd) =...
                        nearly_singular_laplace_3d_constant_stokes(...
                        x, nx, corners, lset, order, true);
                end
                
            end % loop over orders
        end % loop over kernels
    end % loop over distances
end
error = log10(abs(bsxfun(@times, g, 1./g(:,:,:,end)) - 1));

%% plot results
for iX = 1 : length(rvec)
    figure;
    for iT = 1 : length(thetavec)
        subplot(2,2,iT);
        plot(ordervec, squeeze(error(iT, iX, :, :)));
        legend(kernelvec);
        title(sprintf('Theta: %g, R = %g', thetavec(iT)/pi*180, rvec(iX)));
        ylim([-15 2]);
    end
end

%% compare theta contour vs stokes for plane tria
for iX = 1 : length(rvec)
    figure;
    for iT = 1 : length(thetavec)
        subplot(2,2,iT);
        plot(ordervec,...
            log10(abs([
            squeeze(g(iT,iX,4,:))/g(iT,iX,4,end), ... 
            squeeze(gs(iT,iX,4,:))/g(iT,iX,4,end), ...
            squeeze(gst(iT,iX,4,:))/g(iT,iX,4,end)
            ]-1)));
        title(sprintf('Theta: %g, R = %g', thetavec(iT)/pi*180, rvec(iX)));
        legend({'Contour \theta', 'Stokes', 'Stokes-Telles'});
        ylim([-15 2]);
    end
end
