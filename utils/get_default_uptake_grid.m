function uptakeBounds = get_default_uptake_grid(mode)
    if nargin < 1
        mode = 'sweep';
    end

    switch lower(mode)
        case 'sweep'
            uptakeBounds = 0.2 + 5 * linspace(0, 1, 100) + 30.5 * (linspace(0, 1, 100) .^ 3);
        case 'fva'
            uptakeBounds = [1000];
        otherwise
            error('Unknown uptake grid mode: %s', mode);
    end
end
