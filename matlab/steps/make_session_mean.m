function make_session_mean(meanu_imgs, outFile)
% Average the run-level mean images only (fast)
matlabbatch = {};
matlabbatch{1}.spm.util.imcalc.input = meanu_imgs;
matlabbatch{1}.spm.util.imcalc.output = spm_file(outFile,'filename'); % just filename
matlabbatch{1}.spm.util.imcalc.outdir = {spm_file(outFile,'path')};

expr_terms = arrayfun(@(k)sprintf('i%d',k), 1:numel(meanu_imgs), 'UniformOutput', false);
matlabbatch{1}.spm.util.imcalc.expression = sprintf('(%s)/%d', strjoin(expr_terms,'+'), numel(meanu_imgs));

matlabbatch{1}.spm.util.imcalc.var = struct('name', {}, 'value', {});
matlabbatch{1}.spm.util.imcalc.options.dmtx = 0;
matlabbatch{1}.spm.util.imcalc.options.mask = 0;
matlabbatch{1}.spm.util.imcalc.options.interp = 1;
matlabbatch{1}.spm.util.imcalc.options.dtype = 16;

spm_jobman('run', matlabbatch);
end
