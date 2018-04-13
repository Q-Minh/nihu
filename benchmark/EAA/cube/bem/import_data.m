function [freqs, pf, ps, iters] = import_data(pattern)

directory = 'data';

files = dir(fullfile(directory, sprintf('%s_*pf.res', pattern)));
for f = 1 : length(files)
    progbar(1, length(files), f);
    pfname = files(f).name;
    psname = strrep(pfname, 'pf.res', 'ps.res');
    data = importdata(fullfile(directory, pfname));
    freqs(f) = data(1);
    pf(:,f) = complex(data(2:2:end), data(3:2:end));
    data = importdata(fullfile(directory, psname));
    if (data(1,1) ~= freqs(f))
        aaa = 8;
    end
    iters(f) = data(1,2);
    ps(:,f) = complex(data(2:end,1), data(2:end,2));
end

[freqs, i] = sort(freqs);
pf = pf(:,i);
ps = ps(:,i);
iters = iters(i);

end
