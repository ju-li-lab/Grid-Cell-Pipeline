function run_shell(cmd)
% RUN_SHELL  Execute a shell command and assert success.
[status, result] = system(cmd);
if status ~= 0
    error('Shell command failed (exit %d):\n  Command: %s\n  Output: %s', status, cmd, result);
end
end

% ===================== BIDS JSON VALIDATION =====================
