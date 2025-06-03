function polar = loadPolarData(filename)
    % Load polar data from XFOIL file
    data = readmatrix(filename, 'FileType', 'text', 'NumHeaderLines', 12);  % skip XFOIL header
    polar.alpha = data(:, 1);
    polar.Cl    = data(:, 2);
    polar.Cd    = data(:, 3);
end