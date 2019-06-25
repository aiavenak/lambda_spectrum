
function [D,P,B] = f_FastFloyd_und(D)

% works for directed!! :)

%   Input:      D,      weighted/binary undirected matrix, where dij is a
%                       distance measure (not a connection weight).
%
%   Output:     D,      distance (shortest weighted path) matrix
%               P,      number of edges in shortest weighted path matrix
%               B,      elements b_ij indicate the last node
%                       in the shortest path between i and j
%
%   Notes:
%       There may be more than one shortest path between two nodes in the
%       network; the number of edges in the path, given in P, and the
%       elements in B correspond to the first shortest path found by the
%       algorithm.
%
%       Function find_path.m takes as input matrix B to construct the
%       secuence of nodes in the shortest path.
%
%   Algorithm: Fast Floyd
%
%   Andrea Avena-Koenigsberger, IU, 2012
 
P = double(D>0);
D(D == 0) = inf;
n=size(D,2);
B = zeros(n,n);
for i=1:n
  B(:,i) = i;
end

for k=1:n
    
    i2k_k2j = bsxfun(@plus, D(:,k),D(k,:));
    
    path = bsxfun(@gt, D, i2k_k2j);
    [i j] = find(path);  
        
    P(path) = P(i,k) + P(k,j)';
    B(path) = B(i,k);
    
    D = min(D, i2k_k2j);   
    
       
end

D(eye(n)>0)=0;
P(eye(n)>0)=0;



if nargout > 2
    %B = B';
    B(eye(n)>0)=0;
end

    
    
    
    
    
    
    
    
    
    
    