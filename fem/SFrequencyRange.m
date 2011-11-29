function data = SFrequencyRange(data, type, start, stop, steps)

switch lower(type)
    case 'lin'
        data.Excitation.FrequencyVector = linspace(start, stop, steps+1);
    case 'log'
        data.Excitation.FrequencyVector = logspace(log10(start), log10(stop), steps+1);
end
