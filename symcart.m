function [Fext] = symcart(j)
stot=[' '];
for i=1:j
    A = i;
      s1=sprintf(' x%d, ',A);
      s2=sprintf(' y%d, ',A);
      s3=sprintf(' z%d ',A);
      stot=strcat(stot, s1, s2, s3);
end
Fext1=['[' stot ']'];
Fext=str2sym(Fext1);
end

