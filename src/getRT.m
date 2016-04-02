function [Rt] = getRT(E)
 
 E = E';
 ea =E(:,1);
 eb =E(:,2);
 ec =E(:,3);
 
[v idx]= max([norm(cross(ea,eb)), norm(cross(ea,ec)), norm(cross(eb,ec))]);
if idx == 1 a=ea; b=eb;
elseif idx == 2 a = ea; b = ec;
elseif idx == 3 a = eb; b = ec;
end

vc = cross(a,b)/norm(cross(a,b));
va = a/norm(a);
vb = cross(vc,va);

ua = E*va/norm(E*va);
ub = E*vb/norm(E*vb);
uc = cross(ua,ub);

D = [0 1 0; -1 0 0; 0 0 1];
V = [va vb vc];
U = [ua ub uc];

tu = [U(1,3) U(2,3) U(3,3)]';
Ra = U*D*V'
Rb = U*D'*V'

Pa = [Ra tu; 0 0 0 1];
Pb = [Ra -tu; 0 0 0 1];
Pc = [Rb tu; 0 0 0 1];
Pd = [Rb -tu; 0 0 0 1];
Rt = [Pa Pb Pc Pd];