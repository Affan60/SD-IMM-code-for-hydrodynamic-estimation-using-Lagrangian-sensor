%its a stond alone file to compute boundary conditions using water levels
%and gate equations. 

function [y2,y3,v1,v3]=boundary_cond()
fileID = fopen('0KM_rangedata.txt');
C = textscan(fileID,'%f');
fclose(fileID);
whos C;
y2=[C{1}];
y2=flipud(y2);
fileID = fopen('MBL_rangedata.txt');
C = textscan(fileID,'%f');
fclose(fileID);
whos C;
y1=[C{1}];
y1=flipud(y1);
fileID = fopen('3KM_rangedata.txt');
C = textscan(fileID,'%f');
fclose(fileID);
whos C;
y3=[C{1}];
y3=flipud(y3);

diff=size(y2)-size(y1);
new_y1=zeros(length(y1)*9,1);
j=1;

for i=1:length(y1)
    new_y1(j)=y1(i);
    j=j+9;
end


%Interpolation
[act_data]=interpolation(new_y1);
%plot(act_data)
y1=act_data(1:90);
y2=y2./1000;
y1=y1./1000;
bed=5;
y2=bed-y2;
y1=bed-y1;
g=9.8;

%Calculation of water velocity at  upstream gate
for i=1:1:length(y1)
if i>=20&& i<50
    p1=100;
else
    p1=200;
end
Cc = y2(i)/p1;

Cd = Cc/sqrt((1 + Cc*(p1/y1(i))));

v1(i) = (Cd*sqrt((2*g*y1(i))))*3;

q(i) = v1(i)*p1;
end

[y3]=interpolation(y3);
%making bed same for all sensors. 
y3=y3./1000;
y3=bed-y3-0.2;

T=6;
B=6;
p3=3.87;
for i=1:length(y3)
    v3(i)=((0.6*sqrt(g)*B)/(y3(i)*T))*((y3(i)-p3)^(1.5))*3;     
end
end