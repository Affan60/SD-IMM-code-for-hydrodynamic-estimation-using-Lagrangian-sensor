function [x,y]=steady_values(L,yo)
y0 = yo;

yp0 = [0; 0];
options = odeset('RelTol',1e-4,'AbsTol',[1e-6 1e-10], ...
   'Jacobian',{[],[1 y0(1)/y0(2); 0 1]});


%[y0,yp0] = decic(@robertsidae,0,y0,[1 1],yp0,[],options);

xspan = 1:L;
[x,y] = ode15i(@robertsidae,xspan,y0,yp0,options);
% figure(14)
% plot(x,y)
% xlim auto
% ylim auto
% ylabel('y')
% reset(gca)
% reset(gcf)
% xlabel('Distance (m)')
% legend('Water Velcoity','Water Level','Location','northwest')
% legend('boxoff')
% title('Backwater steady states')
end