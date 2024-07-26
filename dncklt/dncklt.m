function [dst, side_information, eigenvalues] = dncklt(strategy, arguments, src, subsampling)
% DNCKLT(STRATEGY, ARGUMENTS, SRC), returns a KLT transformation of SRC.
%   It applies the klt transformation to SRC using a divide-and-conquer  
%   STRATEGY with the required ARGUMENTS.
%   DNCKLT(STRATEGY, ARGUMENTS, SRC), returns DST, the transformation of 
%   SRC using a divide and conquer strategy, SIDE_INFORMATION, the 
%   information required in order to perform the inverse transformation of 
%   DST and EIGENVALUES, which is a vector that contains the eigenvalue 
%   relative to every cluster.
%
%   STRATEGY. Is the divide and conquer strategy used to perform the KLT 
%   transformation on the clusters. Valid values are:
%       - 'single-level'. Single level transformation
%       - 'regular-multi-level'. Multi level transformation
%       - 'pot'. Multi level transformation
%       - 'variable-size-cluster'. Single level strategy.
%       - 'recursive'. Recursive multilevel strategy.
%       - 'static-two-level'. Two level strategy.
%       - 'static-allocation'. Multi level strategy.
%
%   ARGUMENTS. A vector that contains the parameters used to configure  the 
%   STRATEGY applied. Valid values are:
%       - 'single-level'. An integer that indicates the number of clusters 
%           considered.
%       - 'regular-multi-level'. An integer that indicates the number of  
%           clusters considered in the first level.
%       - 'pot'. An empty vector.
%       - 'variable-size-cluster'. A vector of integers that indicates the 
%           size of every cluster.
%       - 'recursive'. An integer that indicates the maximum recursion depth
%       - 'static-two-level'. Two integers, indicating the size of a 
%           cluster in the first level and the second level respectively.
%       - 'static-allocation'. A vector with an odd number of elements, 
%           each of one divide the number of components of SRC.
%
%   SRC. A two-dimensional matrix
%
%   SUBSAMPLING. If subsampling can be performed, it is a boolean
%
%   This funcion returns three values. For instance, 
%  
%   [DST, SIDE_INFORMATION, EIGENVALUES] = DNCKLT(STRATEGY, ARGUMENTS, SRC, SUBSAMPLING)
%
%   DST is a matrix which represents the transformed matrix, SIDE_INFORMATION 
%   is needed to perform the inverse transform. It  contains the size of 
%   every cluster, the permutation applied to this  cluster, and its covariance 
%   matrix. Finally, EIGENVALUES can be used to order the clusters in DST by 
%   its energy level
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

sources = size(src, 1);

side_information = [];
dst = src;
eigenvalues = zeros(sources, 1);

switch strategy
  case 'single-level'
    % single-level strategy

    if(size(arguments, 2) ~= 1 || rem(log2(arguments), 1) ~= 0 || rem(sources, arguments) ~= 0)
        error('dncklt:InvalidArguments', 'For this strategy arguments must be an integer power of two that divides the number of sources.');
    end
    cluster_size = floor(sources/arguments);
    % we must apply the klt transformation to every cluster and store the information needed in the decompression
    % process
    for current_cluster = 1:arguments        
        % we create the needed permutation to move the current cluster to the beginning of the matrix
        source = (current_cluster - 1)*cluster_size + 1:current_cluster*cluster_size;
        cluster_permutation = permutation_bring_to(1:sources, 1:cluster_size, source, false);
        [dst, cluster_struct, A] = apply_transformation(dst, cluster_permutation, cluster_size, subsampling);
        eigenvalues(cluster_permutation(1:cluster_size)) = A;
        side_information = [side_information, cluster_struct];
    end
  case 'regular-multi-level'
    % multi-level strategy

    if(size(arguments, 2) ~= 1 || rem(log2(arguments), 1) ~= 0 || rem(sources, arguments) ~= 0)
        error('dncklt:InvalidArguments', 'For this strategy arguments must be an integer power of two that divides the number of sources.');
    end

    cluster_size = floor(sources/arguments);
    number_of_clusters = arguments;
    
    levels = log2(number_of_clusters) + 1;
    
    % this permutation will be the product of the permutations needed in every level. At first level, it is the 
    % identity.
    permutation = 1:sources;
    for level = 0:levels - 1
        % this is the permutation which contains the component exchaning in this level
        level_permutation = 1:sources;
        
        % we must create a different permutation for every level. The first half of a cluster will be moved to the 
        % front, and the other half will be moved to the back
        for current_cluster = 1:number_of_clusters
            
            if(current_cluster == 1 && level > 0)
                for i = 0:number_of_clusters*2-1
                    % First part of the block
                    len = floor((cluster_size + rem(i, 2))/2);
                    
                    source = i*cluster_size + 1;
                    destiny = floor(i/2)*cluster_size + 1;
                    
                    if(rem(i,2) == 1)
                        destiny = destiny + cluster_size - len;
                    end
                    
                    level_permutation(destiny:destiny+len - 1) = source:source+len - 1;
                    
                    % Second part of the block
                    len = cluster_size - len;
                    
                    source = source + cluster_size - len;
                    destiny = destiny + number_of_clusters*cluster_size;
                    
                    if (rem(i,2) == 1)
                        destiny = destiny + cluster_size - 2*len;
                    end
                    
                    level_permutation(destiny:destiny+len - 1) = source:source+len - 1;
                end
                permutation = permutation(level_permutation);
            end
            
            % the permutation over the cluster at this level must be calculated
            source = (current_cluster - 1)*cluster_size + 1:current_cluster*cluster_size;
            cluster_permutation = permutation_bring_to(1:sources, 1:cluster_size, source, false);

            % and we must consider the previous permutations done
            cluster_permutation = permutation(cluster_permutation);
            [dst, cluster_struct, A] = apply_transformation(dst, cluster_permutation, cluster_size, subsampling);
            eigenvalues(cluster_permutation(1:cluster_size)) = A;
            side_information = [side_information, cluster_struct];
        end
        
        % Finally update number_of_clusters
        number_of_clusters = number_of_clusters/2;
    end
  case 'static-allocation'
    % static cluster size allocation

    if(rem(size(arguments, 2), 2) == 0 || sources < arguments(1))
        error('dncklt:InvalidArguments', 'arguments must have an odd length and cluster size of level 0 has to be smaller than the number of components');
    end
    level_number = 0;
    level_cluster_count = arguments(1);
    level_cluster_size = floor(sources / level_cluster_count);

    old_level_cluster_size = -1;
    
    pick = [];
    eigenvalues = zeros(1, sources);
    permutation = 1:sources;
    
    while (true)

        if(level_cluster_count*level_cluster_size > sources)
            error('dncklt:InvalidComp', 'there are too many arguments');
        end
        
        level_permutation = 1:sources;
        
        if (level_number > 0)
            acc = 0;
            for i = 1:size(pick, 2)
                source = old_level_cluster_size*(i - 1) + 1:old_level_cluster_size*(i - 1) + pick(i);
                destiny = acc + 1:acc + pick(i);
                level_permutation = permutation_bring_to(level_permutation, destiny, source, true);
                acc = acc + pick(i);
            end
        end
        permutation = permutation(level_permutation);
        
        for current_cluster = 0:level_cluster_count - 1
            
            cluster_permutation = 1:sources;
            source = level_cluster_size*current_cluster + 1:level_cluster_size*(current_cluster + 1);
            cluster_permutation = permutation(permutation_bring_to(cluster_permutation, source, 1:level_cluster_size,...
                                                                   false));
            [dst, cluster_struct, A] = apply_transformation(dst, cluster_permutation, level_cluster_size, subsampling);
            eigenvalues(cluster_permutation(1:level_cluster_size)) = A;
            side_information = [side_information, cluster_struct];
        end
        
        old_level_cluster_count = level_cluster_count;
        old_level_cluster_size = level_cluster_size;
        if(2*(level_number + 1) >= size(arguments,2))
            return
        end
        level_cluster_count = arguments(2*(level_number + 1) + 1);
        pick(1:old_level_cluster_count) = arguments(2*level_number + 2);
        level_cluster_size = floor(old_level_cluster_count*arguments(2*level_number + 2)/level_cluster_count);
        level_number = level_number + 1;
    end
    
  case 'pot'
    % pot
    if(arguments ~= 2)
        error('dncklt:invalidValue', 'in pot strategy, arguments must be equal to 2');
    end
    cluster_size = 2;
    level_number = 0;
    
    level_cluster_count = floor(sources/cluster_size);
    remaining_component = rem(sources, cluster_size);
    
    % In this permutation we accumulate the composition of the permutations in every level (level_permutation)
    permutation = 1:sources;
    
    while (true)
        
        if (level_cluster_count * cluster_size > sources)
            error('dncklt:InvalidComp', 'there are too many arguments');
        end
        % we must apply the klt transformation to all the clusters considered
        for current_cluster = 1:level_cluster_count
            source = (current_cluster - 1)*cluster_size + 1:current_cluster*cluster_size;
            destiny = 1:cluster_size;
            cluster_permutation = permutation(permutation_bring_to(1:sources, destiny, source, false));
            [dst, cluster_struct, A] = apply_transformation(dst, cluster_permutation, cluster_size, subsampling);
            eigenvalues(cluster_permutation(1:cluster_size)) = A;
            side_information = [side_information, cluster_struct];
        end
        
        %  Are we done yet?
        if (level_cluster_count <= 1 && remaining_component == 0)
            if (level_cluster_count ~= 1)
                error('dncklt:Rec', 'a cluster is expected in the last level');
            end
            return
        end
        level_number = level_number + 1;
        old_level_cluster_count = level_cluster_count;
        
        %  this is the permutation in a concrete level
        level_permutation = 1:sources;
        level_permutation(1:old_level_cluster_count) = 1:cluster_size:1+cluster_size*(old_level_cluster_count-1);
        level_permutation(1+old_level_cluster_count:2*old_level_cluster_count) = ...
                         level_permutation(1:old_level_cluster_count)+1;
        
        %  check for the remaining component
        
        if(remaining_component > 0)
            source = cluster_size*old_level_cluster_count + 1;
            destiny = old_level_cluster_count + 1;
            level_permutation = permutation_bring_to(level_permutation, destiny, source, false);
        elseif(remaining_component < 0)
            source = 1:old_level_cluster_count;
            destiny = 2:old_level_cluster_count + 1;
            level_permutation = permutation_bring_to(level_permutation, destiny, source, true);
            level_permutation = permutation_bring_to(level_permutation, cluster_size*...
                                                     old_level_cluster_count + 1, ...
                                                     old_level_cluster_count + 1, false);
        end
        kept_components = level_cluster_count + abs(remaining_component);
        level_cluster_count = floor(kept_components/cluster_size);
        
        if(remaining_component > 0)
            remaining_component = -rem(kept_components, cluster_size);
        else
            remaining_component = rem(kept_components, cluster_size);
        end
        
        if(remaining_component < 0)
            permutation_size = level_cluster_count*cluster_size;
            source = 2:permutation_size + 1;
            destiny = 1:permutation_size;
            level_permutation = permutation_bring_to(level_permutation,destiny, source, true);
        end
        permutation = permutation(level_permutation);
    end
    
  case 'variable-size'
    % variable size cluster of one level
    
    if(size(arguments, 2) == 1 || sum(arguments) ~= sources)
        error('dncklt:InvalidSize', 'cluster_size must be a vector the sum of the rows must be equal to the sum of the elements');
    end
    
    acc = 0;
    for current_cluster=1:size(arguments, 2)
        % In this case cluster_size is a vector
        current_cluster_size = arguments(current_cluster);
        
        % First, we create the permutation for a cluster, which brings the acc + 1:acc + current_cluster_size to the
        % first positions of the matrix
        cluster_permutation = permutation_bring_to(1:sources, 1:current_cluster_size,... 
                                                   acc + 1:acc + current_cluster_size, true);
       
        [dst, cluster_struct, A] = apply_transformation(dst, cluster_permutation, current_cluster_size, subsampling);
        eigenvalues(cluster_permutation(1:current_cluster_size)) = A;
        side_information = [side_information, cluster_struct];
      
        acc = acc + current_cluster_size;
    end
  case 'recursive'
    if(size(arguments, 2) > 1)
        error('dncklt:InvalidSize', 'cluster_size must be an integer');
    end
    
    [dst, side_information] = recursive_mode(sources, arguments, dst, subsampling);

  case 'static-two-level'
    % static-two-level structure
    
    %  we must check that cluster_size is a vector
    if(size(arguments, 2) ~= 2)
        error('dncklt:InvalidSize', 'cluster_size must be a vector of two components');
    end
    first_level_clusters = arguments(1);
    second_level_clusters = arguments(2);
    first_level_cluster_size = floor(sources / first_level_clusters);
    
    if(rem(sources, first_level_clusters) ~= 0)
        first_level_cluster_size = first_level_cluster_size + 1;
    end
    second_level_sources = [];
    for i = 1:second_level_clusters
        second_level_sources = [second_level_sources, struct('list', [])];
    end
    
    
    % first level
    for current_start=0:first_level_cluster_size:sources - 1
        current_cluster_size = first_level_cluster_size;
        if (current_start + current_cluster_size > sources)
            current_cluster_size = sources - current_start;
        end
        source = current_start + 1:current_start + current_cluster_size;
        destiny = 1:current_cluster_size;
        
        cluster_permutation = permutation_bring_to(1:sources, destiny, source, true);
        [dst, cluster_struct, A] = apply_transformation(dst, cluster_permutation, current_cluster_size, subsampling);
        eigenvalues(cluster_permutation(1:current_cluster_size)) = A;
        side_information = [side_information, cluster_struct];
     
        for i = 0:min(second_level_clusters, current_cluster_size) - 1
            second_level_sources(i + 1).list = cat(2, second_level_sources(i + 1).list, current_start + 1 + i);
        end
    end
    
    % second level
    for i = 1:size(second_level_sources, 2)
        if (size(second_level_sources(i).list, 2) > sources)
            error('dncklt:st2lev', 'too many elements in permutation of second level');
        end
        if (size(second_level_sources(i).list, 2) > 0)
            cluster_permutation = permutation_bring_list_to(1:sources, second_level_sources(i).list, 1);
            [dst, cluster_struct, A] = apply_transformation(dst, cluster_permutation,... 
                                                            size(second_level_sources(i).list, 2), subsampling);
            eigenvalues(cluster_permutation(1:size(second_level_sources(i).list, 2))) = A;
            side_information = [side_information, cluster_struct];
        end
    end
  otherwise
    warning('dncklt:warning', 'invalid strategy, empty klt_matrix');
    dst = [];
end
end


function [dst, side_information, eigenvalues] = recursive_mode(components, recursion_depth, M, subsampling)
% recursive_mode performs the recursive strategy in order to apply the klt to the clusters created.
% Args:
% * components. Number of components of the image
% * recursion_depth. Maxium recursion depth allowed in the strategy
% * M. Original matrix
% * subsampling. If subsampling must be considered or not
% Returns:
% * dst. The klt transformation of M
% * side_information. The structure which contains the needed information in order to apply the inverse transformation 
%   to dst.
% * eigenvalues. Eigenvalues of every component of the image
[dst, side_information, eigenvalues, p] = recursive_structure(0, components, recursion_depth, M, [], 1:components,...
                                                              zeros(components, 1), subsampling);
end

function [dst, side_information, E, p] = recursive_structure(first, last, recursion_depth, original_matrix, structure_array, permutation, eigen_values, subsampling)
% recursive_structure is a recursive function that allows to perform a recursive divide and conquer strategy
% Args:
% * first. First component of the cluster
% * last. Last component of the cluster
% * recursion_depth. Maxium recursion depth allowed in the strategy
% * original_matrix. Original matrix
% * structure_array. Array which contains all the information about the clusters Already processed.
% * permutation. Accumulated permutation
% * eigen_values. Eigenvalue of every component.
% * subsampling. If subsampling must be considered or not
%
% Returns:
% * dst. The klt transformation of M
% * side_information. The structure which contains the needed information in order to apply the inverse transformation 
%   to dst.
% * E. Eigenvalues of every component of the image
% * p. The accumulated permutation

dst = original_matrix;
sources = size(original_matrix, 1);
side_information = structure_array;
s = size(dst);
p = permutation;
E = eigen_values;
count = last - first;
if (count < 0)
    error('dncklt:rec', 'first component can not be greater than last component');
end
% If there are not more than one component to be considered, we have finished the recursion process
if(count < 2)
    return;
end

% if this is the last time we must apply the klt transformation
if (count == 2 || recursion_depth == 0)
    cluster_permutation = 1:sources;
    source = first + 1:first + count;
    destiny = 1:count;
    cluster_permutation = permutation_bring_to(cluster_permutation, destiny, source, true);
    cluster_permutation = permutation(cluster_permutation);
    [dst, cluster_struct, A] = apply_transformation(dst, cluster_permutation, count, subsampling);
    E(cluster_permutation(1:count)) = A;
    side_information = [side_information, cluster_struct];
    return
end
first_half = floor(count/2);
second_half = floor(count/2);

if(rem(recursion_depth, 2) == 0)
    first_half = first_half + rem(count, 2);
else
    second_half = second_half + rem(count, 2);
end
[dst, side_information, E, p] = recursive_structure(first, first + first_half, recursion_depth - 1, dst,...
                                                    side_information, p, E, subsampling);
[dst, side_information, E, p] = recursive_structure(first + first_half, last, recursion_depth - 1, dst,...
                                                    side_information, p, E, subsampling);

first_stage_result_count = floor(count/2) + rem(count,2);
first_half_results = floor(first_stage_result_count/2);
second_half_results = floor(first_stage_result_count/2);

if(rem(recursion_depth, 2) == 0)
    first_half_results = first_half_results + rem(first_stage_result_count, 2);
else
    second_half_results = second_half_results + rem(first_stage_result_count, 2);
end

% Move second half results next to the first half results
level_permutation = 1:sources;
source = first + 1 + first_half:first + first_half + second_half_results;
destiny = first + 1 + first_half_results:first + first_half_results + second_half_results;
level_permutation = permutation_bring_to(level_permutation, destiny, source, true);
level_permutation = permutation_interleave_range(level_permutation, first, first_stage_result_count, ...
                                                 recursion_depth);
p = p(level_permutation);

% Do the second stage
[dst, side_information, E, p] = recursive_structure(first, first + first_stage_result_count, recursion_depth - 1,...
                                                    dst, side_information, p, E, subsampling);
end


function p = permutation_bring_to(permutation, destiny, source, handle_intersection)
% permutation_bring_to performs a permutation between source and destiny
% Args:
% * permutation. Permutation that will be modified
% * destiny. Position where source components will be stored
% * source. Source components that will be moved to destiny
% * handle_intersection. If intersection between destiny and source must be handled
%
% Returns:
% * p. The permutation 'permutation' but in which 'source' and 'destiny' have been permuted

p = permutation;
p(destiny) = permutation(source);
if (handle_intersection)
    s = setdiff(source, destiny);
    destiny = setdiff(destiny, source);
    source = s;
end
p(source) = permutation(destiny);
end

function p = permutation_interleave_range(permutation, from, len, recursion_depth)
% Interleaves positions of the permutation
% Args:
% * permutation. Permutation that will be modified
% * from. First component that will be considered for the operation
% * len. Number of components that will be considered for the operation
% * recursion_depth. Recursion depth
% Returns:
% * p. The permutation 'permutation' modified

v = 1:len;
p = permutation;

if(rem(recursion_depth, 2))
    start_a = from;
    start_b = from + floor(len/2);
else
    start_a = from + ceil(len/2);
    start_b = from;
end
array = 1:2:len;
s = size(array, 2);
v(array) = permutation(start_b + 1:start_b + s);
array = 2:2:len;
s = size(array, 2);
v(array) = permutation(start_a + 1:start_a + s);

p(from + 1:from + len) = v(1:len);
end

function p = permutation_bring_list_to(permutation, list, to)
% permutation_bring_list_to add a permutation to the current one
% Args:
% * permutation. Permutation that will be modified
% * list. List of components that will be moved
% * to. Position where list will be stored
% Returns:
% * p. The permutation 'permutation' modified

p = permutation;
list_outside_destination = union(list(list < to), list(list >= to + size(list, 2)));
set_inside_destination = intersect(list(list >= to), list(list < to + size(list, 2)));
values = p(list);
if(size(list_outside_destination, 2) + size(set_inside_destination, 2) ~= size(list, 2))
    error('dncklt:checkError', 'the number of locations inside and outside destination not sum the total');
end
used_out_side = 1;
for pos = to:to + size(list, 2) - 1
    if(~ismember(set_inside_destination, pos))
        dest = list_outside_destination(used_out_side);
        used_out_side = used_out_side + 1;
        p(dest) = p(pos);
    end
end
p(to:size(list, 2)) = values(1:size(list, 2));
end


function [dst, cluster_struct, A] = apply_transformation(dst, cluster_permutation, cluster_size, subsampling)
% Prepares the cluster for the transformation and applies the KLT
% Args:
% * dst. Original matrix that contains the cluster
% * cluster_permutation. Permutation that must applied to dst in order to obtain the desired cluster
% * cluster_size. Size of the cluster
% * subsampling. If the covariance matrix must be calculated considering a
%   subsampling of the cluster
%
% Returns:
% * dst. The matrix with the cluster transformed
% * cluster_struct. Information about the transformation of the cluster
% * A. Vector with the eigenvalues that correspond to each component of the cluster.

    sources = size(dst, 1);
    samples = size(dst, 2);
    M = dst(cluster_permutation(1:cluster_size), :);
    if (subsampling)
        side_information = sampled_cov(M, ceil((sources*samples)/100));
    else
        side_information = cov(M');
    end
    [Q, A] = eig(side_information);
    v = diag(A);
    [A, p] = sort(v, 'descend');
    Q = Q(:, p);
    R = Q'*M; 
    dst(cluster_permutation(1:cluster_size), :) = R;
    cluster_struct = struct('cluster_size', cluster_size,  'permutation', cluster_permutation, 'eigen_matrix', Q);
end

function C = sampled_cov(cluster, length)
% Return an aproximation of covariance a matrix using a sample subset of length 'length'.
% Args:
% * cluster. Cluster used to calculate the covariance matrix
% * length. Length of the sampling considered
%
% Returns:
% * C. Approximated covariance matrix of cluster

s = size(cluster);
% we must check that length is not bigger than the number of pixels in the
% matrix
rand_cluster_size = ceil(length/s(1));

if rand_cluster_size < 2 || length > s(1)*s(2) || s(1) == 1
    % if there are not enough samples in the
    % cluster
    C = cov(cluster');
else
    rand_matrix = rand_matrix_gen([s(1) rand_cluster_size ], s(2));
    S = cluster';
    S = S(rand_matrix);
    C = S*S'/(s(1) - 1);  
end
end

function rand_matrix = rand_matrix_gen(matrix_size, cluster_size)
% Generates a pseudorandom matrix in which every element of the matrix
% is an index of a sample of a cluster.
% Args:
%  * matrix_size. Size of random_matrix to be generated.
%  * cluster_size. Number of samples in a cluster
% Returns:
%  * rand_matrix. Matrix of pseudorandom integers 

m = floor(rand(1, matrix_size(2))*matrix_size(2)) + 1;
rand_matrix = zeros(matrix_size(1), matrix_size(2));
for j = 1:matrix_size(1)
    rand_matrix(j, :) = m + cluster_size*(j - 1);
end
end
