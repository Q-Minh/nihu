function [phat] = itrans_inf_invar(ptil, kz)

dk = kz(2)-kz(1);
phat = dk/2/pi*fft(fftshift(ptil,2), [], 2);
phat(:,end+1,:) = phat(:,1,:);