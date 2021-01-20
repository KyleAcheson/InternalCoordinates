clear all

Q = dlmread('bentgeom.xyz'); % ref geometry 
[natom, nc] = size(Q);
        
B = sym(zeros(natom^2, 3*natom)); % init B matrix (derivative of flattened distance mat. wrt cartesian coordinates)
BB = zeros(natom^2, 3*natom);
[Nsq, tN] = size(B);

for i=1:natom
    for j=i+1:natom
        D(i,j) = norm(Q(i,:)-Q(j,:)); % distance matrix 
        D(j,i) = 0; 
    end                                     
end

d = D(:); % flatten distance matrix
        
Nabla = symcart(3);
R(1,:) = Nabla(1:3);
R(2,:) = Nabla(4:6);
R(3,:) = Nabla(7:9);

count = 0;
for i=1:natom
    for j=1:natom
        if i == j
            Dij = 0;
        else
            Dij = 1/norm(R(i,:)-R(j,:));
        end
        count = count + 1;
        B(count,:) = gradient(Dij,Nabla);
        for k=1:tN
            BB(count, k) = subs(B(count,k),Dij,1/(norm(Q(i,:)-Q(j,:))));
        end
    end
end

G = BB*BB';