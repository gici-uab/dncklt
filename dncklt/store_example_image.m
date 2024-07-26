function [] = store_example_image()
% STORE_EXAMPLE_IMAGE stores an image to disk.
%   It generates an image that is used in DNCKLT_EXAMPLE

    % image_geometry(1): height of a component of the image
    % image_geometry(2): width of a component of the image
    % image_geometry(3): number of components of the image
    image_geometry = [2, 3, 8];
    
    % input_matrix is a three dimensional matrix, so, it must be
    % reshaped to a two dimensional matrix
    example_matrix = cat(3, [[1 3 5];[2 6 10]], ... % first band
                            [[4 9 14];[7 13 19]], ... % second band
                            [[8 15 22];[11 20 29]], ... % third band
                            [[12 23 34];[16 28 40]], ... % fourth band
                            [[17 30 43];[18 32 46]], ... % fifth band
                            [[21 36 51];[24 40 56]], ... % sixth band
                            [[25 42 59];[26 44 62]], ... % seventh band
                            [[27 46 65];[31 52 73]]); % eigth band

    [fid, msg] = fopen('example_image.raw', 'w');     
    fwrite(fid, example_matrix, 'int16');
    fclose(fid);
end

