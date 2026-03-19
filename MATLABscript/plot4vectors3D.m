function plot4vectors3D(vectors, titles, vis, test)

% plots four column vectors of 3-dimensional points (x1, y1, z1, x2, y2,
% z2....) separately in 2x2 subplot

% vis is a boolean that determines whether the visualization of the
% symmetry vectors is off or on


for i = 1:length(titles)

    if strcmp(titles{i}, 'Initial Y') || strcmp(titles{i}, 'Final Y')
        plot1vector3D(vectors(:, i), '', titles{i}, vis, test);
    end

    if strcmp(titles{i}, 'Initial X') || strcmp(titles{i}, 'Final X')
        plot1vector(vectors(:, i), '', titles{i}, vis, test);
    end

end
