function res = robertsidae(x,y,yp)
g=9.8;
res = [yp(1) + (y(1)/y(2))*yp(2);
       yp(2) - ((0.001-((0.02^2)*y(1)*abs(y(1)))/(y(2)^(4/3)))/(1-(y(1)^2/(g*y(2)))))];
end
