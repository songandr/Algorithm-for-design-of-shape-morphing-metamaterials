function mathematica_coords = makeMathematicaCoords(data, config)
    n = length(data);
    if config == "x" || config == "xi"
        % Check even length
        if mod(n,2) ~= 0
            error('Length of input must be even (2N)');
        end
    
        N = n / 2;
        mathematica_coords = strings(N,1);
    
        for i = 1:N
            x = data(2*i - 1);
            y = data(2*i);
    
            % Format with high precision (Mathematica-friendly)
            if config == "xi"
                str = sprintf('xi%d = {%.15g, %.15g};', i, x, y);
            else
                str = sprintf('x%d = {%.15g, %.15g};', i, x, y);
            end
            % Replace scientific notation: e± → *^±
            str = regexprep(str, 'e([+-]?\d+)', '*^$1');
            mathematica_coords(i) = str;
        end
    elseif config == "y" || config == "yi"
        % Check length
        if mod(n,3) ~= 0
            error('Length of input must be (3N)');
        end
    
        N = n / 3;
        mathematica_coords = strings(N,1);
    
        for i = 1:N
            x = data(3*i - 2);
            y = data(3*i - 1);
            z = data(3*i);
    
            % Format with high precision (Mathematica-friendly)
            if config == "yi"
                str = sprintf('yi%d = {%.15g, %.15g, %.15g};', i, x, y, z);
            else
                str = sprintf('y%d = {%.15g, %.15g, %.15g};', i, x, y, z);
            end
            % Replace scientific notation: e± → *^±
            str = regexprep(str, 'e([+-]?\d+)', '*^$1');
            mathematica_coords(i) = str;
        end
    else
        error('Specify x or y coords ');
    end
end