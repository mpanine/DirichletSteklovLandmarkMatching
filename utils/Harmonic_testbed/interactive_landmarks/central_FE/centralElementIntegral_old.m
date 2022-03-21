function I = centralElementIntegral_old(t1, r1, t2, r2, tmax, R, rM)
%CENTRALELEMENTINTEGRAL computes the integral of the product of gradients
%t1, r1: theta and radius derivative of function 1 (function handles)
%t2, r2: theta and radius derivative of function 2 (function handles)
%tmax: max angle
%R: inner radius of finite element (aka epsilon)
%rM: max radius as function of angle (function handle)

f = @(t,r) t1(t,r).*t2(t,r)./r + r1(t,r).*r2(t,r).*r;


figure
[TT,RR] = meshgrid(0:tmax/100:tmax, R:0.001:0.02);
FF = f(TT,RR);
surf(RR.*cos(TT),RR.*sin(TT), FF )




% I = integral2(f, 0, tmax, R, rM, 'AbsTol', 0, 'RelTol', 1e-9);
I = quad2d(f, 0, tmax, R, rM, 'FailurePlot',true,'Singular',true, 'MaxFunEvals', 2000);


end