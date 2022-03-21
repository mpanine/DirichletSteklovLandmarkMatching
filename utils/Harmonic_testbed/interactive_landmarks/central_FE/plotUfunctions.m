function plotUfunctions(ind1, ind2, l1, a, b, Rmin)
%PLOTUFUNCTIONS plots the u functions given by ind1 aqnd ind2

u1 = U(ind1, l1, a, b, Rmin);
u2 = U(ind2, l1, a, b, Rmin);

rM = maxR(l1, a, b);

tmax = atan(b/a);

tdivs = 100;
rdivs = 50;

pret = linspace(0,tmax,tdivs);

T = repmat(pret,rdivs,1);
R = zeros(size(T));

for i = 1:tdivs
    
    R(:,i) = linspace(Rmin, rM(pret(i)), rdivs);  
    
end

U1 = u1(T,R);
U2 = u2(T,R);


X = R .* cos(T);
Y = R .* sin(T);


surf(X,Y,U1)
hold on
surf(X,Y,U2)

axis equal
xlim([0,max(l1,a)])
ylim([0, b])
zlim([0,1.5])



end