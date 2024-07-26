function dst = idncklt(src, side_information)
% IDNCKLT(SRC, SIDE_INFORMATION), performs the inverse KLT on SRC.
%   It applies the inverse KLT transformation on SRC, which is a two 
%   dimensional matrix, considering the  information in SIDE_INFORMATION 
%   which is an array of structs that is generated when DNCKLT is called.
%
%   The resulf of IDNCKLT is a two dimensional matrix which represents the
%   inverse transformation of SRC. If subsampling has not been used, then
%   DST is the same matrix that was transformed using DNCKLT. However, if
%   subsampling has been used, then DST is a matrix that should be similar
%   to the matrix transformed using DNCKLT but not the same.
%
% The different strategies are describes in the article: I. Blanes, 
% Joan Serra-Sagrista, Michael W. Marcellin and Joan Batrina-Rapesta, 
% "Divide-and-conquer strategies for hyperspectral image processing", 
% Signal Processing Magazine
%
% Coded by: Jose Enrique Sanchez and Estanislau Auge
%
% License: This file is distributed under the terms of the GNU Affero General Public License (AGPL) version 3,  WITH 
% AN ADDITIONAL CLAUSE: if you find it useful, please send an email to any of the authors (so that we can include it in 
% our grant reports). Contact the author for other licensing terms.
  
dst = src;
total_clusters = size(side_information, 2);

% To apply the inverse transform, each cluster must be inverted in opposite order as it was applied.
for current_cluster = total_clusters:-1:1
    p = side_information(current_cluster).permutation(1:side_information(current_cluster).cluster_size);
    Q = side_information(current_cluster).eigen_matrix;
    permuted_matrix = dst(p, :);
    dst(p, :) = Q*permuted_matrix;
end

end
