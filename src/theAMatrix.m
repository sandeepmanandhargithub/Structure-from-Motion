function A = theAMatrix(x1, x2)
x1 = x1';
x2 = x2';
npts = size(x1,2);
  A = [x1(1,:)'.*x2(1,:)'   x1(1,:)'.*x2(2,:)'  x1(1,:)' ...
         x1(2,:)'.*x2(1,:)'   x1(2,:)'.*x2(2,:)'  x1(2,:)' ...
         x2(1,:)'             x2(2,:)'            ones(npts,1) ];    
end