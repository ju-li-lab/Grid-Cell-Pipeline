function vec = pe_dir_to_vector(peDir)
% PE_DIR_TO_VECTOR  Convert FSL PE direction string to acqparams vector.
%   'x'  -> [1 0 0],  'x-' -> [-1 0 0]
%   'y'  -> [0 1 0],  'y-' -> [0 -1 0]
%   'z'  -> [0 0 1],  'z-' -> [0 0 -1]

switch peDir
    case 'x',  vec = [1  0  0];
    case 'x-', vec = [-1 0  0];
    case 'y',  vec = [0  1  0];
    case 'y-', vec = [0 -1  0];
    case 'z',  vec = [0  0  1];
    case 'z-', vec = [0  0 -1];
    otherwise, error('Invalid PE direction: %s (use x, x-, y, y-, z, or z-)', peDir);
end
end
