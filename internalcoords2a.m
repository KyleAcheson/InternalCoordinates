clear all

Q = dlmread('bentgeom.xyz'); % ref geometry
[natom, nc] = size(Q);

for i=1:natom
    for j=i+1:natom
        D(i,j) = norm(Q(i,:)-Q(j,:)); % distance matrix - numerical values
        D(j,i) = 0; % upper triangular
    end                                     
end

d = D(:); % flatten distance matrix

Nabla=symcart(natom);          
R(1,:) = Nabla(1:3); % define symbolic cartesian vectors for each of the three atoms
R(2,:) = Nabla(4:6);
R(3,:) = Nabla(7:9);

        
B = sym(zeros(natom^2, 3*natom)); % init symbolic B matrix (derivative of flattened distance mat. wrt cartesian coordinates)
BB = zeros(natom^2, 3*natom); % init gradient matrix for storing numerical values of evaluated analytical elements of symbolic B
[Nsq, tN] = size(B);  % rows of B/BB corrospond to elements of flattened distance matrix (natom^2 rows)
                      % columns of B/BB corrospond to each derivative wrt
                      % nabla - therefore (3*natom cols) - 3/9 are 0
                     
count = 0;
for i=1:natom % Loop over all cartesian coordinates
    for j=1:natom
        if i == j % avoids dividing by 0 for diag elements
            Dij = 0;
        else
            Dij = 1/norm(R(i,:)-R(j,:)); % take 1/Euclidean distance for each off-diag
                                         % element - this is symbolic    
        end
        count = count + 1;
        B(count,:) = gradient(Dij,Nabla); % for each symb distance take gradient wrt nabla and store the 9 derivatives in cols
        ax = 0; % axis counter (xi,yi,zi)
        for k=1:tN
            if mod(ax,3) == 0 % must reset every 3 iterations to loop over axes three times
                ax = 0;
            end
            ax = ax + 1; 
            b = subs(B(count,k),Dij,1/(norm(Q(i,:)-Q(j,:)))); % evaluate denominators of symb matrix B i.e. ||Ri-Rj||^-3
            c = subs(b,(R(i,ax)-R(j,ax)),(Q(i,ax)-Q(j,ax)));  % evaluate numerators i.e. +/-(Ri_ax-Rj_ax) for each i and j and where
            BB(count,k) = c;                                  % ax runs over (x,y,z)
        end
    end
end


G = BB*BB'; 
G(isnan(G))=0;
G(isinf(G))=0;

[eigVec, eigVal] = eig(G);  % Diag matrix & extract eigvals and vecs
diagVals = diag(eigVal);

count = 1;
for i=1:tN  % Extract eigenvectors that have a non-zero eigenvalue 
    if diagVals(i) > 1.d-3 % non-zero defined by some thresh
        U(:, count) = eigVec(:, i); % U contains vecs of non-zero eigvals
        fVal(count) = diagVals(i);  % should be natom (3) no-zero eigvals
        count = count + 1;
    else
        count = count;
    end
end

q = U'*d;  % Transformation of distance matrix into new internal coordinate system