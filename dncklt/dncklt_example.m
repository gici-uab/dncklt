% Implementation of seven different divide and conquer strategies in order to provide an approximation to the KLT 
% transformation. Every strategy provides different trade-off among its characteristics. Input images are expected to 
% be two dimensional matrix in which every row is a whole component of a three dimensional image.
%
% The different strategies are describes in the article: I. Blanes, Joan Serra-Sagristà, Michael W. Marcellin and 
% Joan Batrina-Rapesta, "Divide-and-conquer strategies for hyperspectral image processing", Signal Processing Magazine
%
% Coded by: Jose Enrique Sánchez and Estanislau Augé
%
% License: This file is distributed under the terms of the GNU Affero General Public License (AGPL) version 3,  WITH AN 
% ADDITIONAL CLAUSE: if you find it useful, please send an email to any of the authors (so that we can include it in our 
% grant reports). Contact the author for other licensing terms.

% if subsampling must be used when calculating the covariance matrix
% of a cluster. It must be considered that when subsampling is used, the
% inverse transformation of the transformed matrix will not obtaint the 
% same matrix that was transformed.
subsampling = false;

% we must set cluster_mode and cluster_size which will be used to determine in which way the transformation will be 
% performed. Cluster_mode possible values:
% - 'single-level'. Single level transformation
% - 'regular-multi-level'. Multi level transformation
% - 'pot'. Multi level transformation
% - 'variable-size'. Single level transformation
% - 'recursive'. Recursive structure transformation
% - 'static-two-level'. Two level structure using a static strategy
% - 'static-allocation'. Multi level transformation
cluster_mode = 'static-allocation';

% clusters, correct values:
% - 'single-level'. An integer is required
% - 'static-multi-level'. An integer is required
% - 'pot'. No value is required, so the value given will be ignored
% - 'variable-size'. An array in which every element is the size of a cluster. The sum of all elements in this 
%    array must be equal to the number of components of the image that is going to be transformed.
% - 'recursive-structure'. An integer is required.
% - 'static-two-level'.
% - 'static-allocation'.
clusters = [4 2 1];

% image_geometry(1): height of a component of the image
% image_geometry(2): width of a component of the image
% image_geometry(3): number of components of the image
image_geometry = [2, 3, 8];

% we write the image that we would like to transform to disk as an example
store_example_image();

% first, we must read the file 'image.raw'
[fid, msg] = fopen('example_image.raw', 'r');
[data, count] = fread(fid, prod(image_geometry), 'int16');
fclose(fid);

% we reshape the data to a three-dimensional image, which is what the image really is
input_matrix = reshape(data, image_geometry);
band_size = prod(image_geometry(1:2));

% we need to reshape the three-dimensional matrix to a two-dimensional matrix where every line is a component of the 
% original image
M = reshape(input_matrix, [band_size, image_geometry(3)])';
p = [1; 3; 5; 2; 4; 6];
M = M(:, p);

% we may need to center the data in input_matrix before processing
m = mean(M');
M = bsxfun(@minus, M', m)';

% once everything is ready, we use the function to perform the klt transformation using a divide and conquer strategy
[transformed_data, side_information, E] = dncklt(cluster_mode, clusters, M, subsampling);

% now we do the inverse tranformation in order to obtain the original matrix
recovered_matrix = idncklt(transformed_data, side_information);

% at the end, we must undo the center operation
recovered_matrix = bsxfun(@plus, recovered_matrix', m)';
ip = [1; 4; 2; 5; 3; 6];
recovered_matrix = reshape(recovered_matrix(:, ip)', image_geometry);

% at the end we remove the 'image.raw', as it is not needed anymore
delete('example_image.raw');
