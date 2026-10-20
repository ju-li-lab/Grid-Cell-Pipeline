function s = clock_of(sec)
if isnan(sec)
    s = 'unknown';
else
    s = sprintf('%02d:%02d:%02d', floor(sec/3600), floor(mod(sec,3600)/60), floor(mod(sec,60)));
end
end
