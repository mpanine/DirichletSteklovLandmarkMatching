function u = U(index, l1, a, b, R)
%U: one of the three standard finite element functions
%Vertices at (0,0), (l1,0), (a,b)
%R: internal radius of the central element (sometimes called epsilon)
%index = 0,1,2 the index of the desired function (0: (0,0), 1: (l1,0), 2: (a,b) )

if index == 0
    
    u = @(t,r) 1+((-1).*r+R).*((-1).*R+b.*l1.*(b.*cos(t)+((-1).*a+l1).*sin(t)).^( ...
              -1)).^(-1);

   
elseif index == 1
    
    u = @(t,r) (-1).*(r+(-1).*R).*(b.*cos(t)+(-1).*a.*sin(t)).*((-1).*b.*l1+b.* ...
               R.*cos(t)+(-1).*(a+(-1).*l1).*R.*sin(t)).^(-1);
    
    
elseif index == 2
    
    u = @(t,r) (-1).*l1.*(r+(-1).*R).*sin(t).*((-1).*b.*l1+b.*R.*cos(t)+(-1).*(a+ ...
               (-1).*l1).*R.*sin(t)).^(-1);
    
else
    
    error('The index of the desired function must be equal to 0, 1 or 2!\n')

end

end

