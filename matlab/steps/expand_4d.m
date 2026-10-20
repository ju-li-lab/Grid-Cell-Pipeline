function vols = expand_4d(niiFile)
V = spm_vol(niiFile);
vols = cell(numel(V),1);
for k=1:numel(V)
    vols{k} = sprintf('%s,%d', niiFile, k);
end
end

%
