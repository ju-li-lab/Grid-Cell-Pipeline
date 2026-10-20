function v = cfg_or(cfg, key, default)
% CFG_OR  Read a setting out of a config.sh struct, with a fallback.
%
%   v = cfg_or(cfg, 'ROI_MODE', 'both')
%
% Returns DEFAULT when the key is absent or set to an empty value, which is
% how config.sh spells "use the default" (e.g. CIRCSTAT_DIR=).
%
% Values are strings, as read_config returns them. Convert at the point of
% use: str2double(cfg_or(cfg,'MAX_WORKERS','1')).

if isfield(cfg, key) && ~isempty(cfg.(key))
    v = cfg.(key);
else
    v = default;
end
end
