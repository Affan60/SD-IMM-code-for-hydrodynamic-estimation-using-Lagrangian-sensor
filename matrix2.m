function [A,B,A_size,N,hbar,vbar]=matrix2(L,t,del_x,Vo,Ho)
yo = [Vo; Ho];
[x,y]=steady_values(L,yo);
N=ceil(L/del_x);
pick=1;
size_y=size(y);
for k=1:1:N
hbar(k)=y(pick,2);
vbar(k)=y(pick,1);
pick=pick+del_x-1;
end
m=0.002;
g=9.8;           %Gravitaional Force
%Length of Each Grid Point
del_t=t;         %Time Step
A_size=(N*2)-4;
B_size=4;
A=zeros(A_size,A_size);
alpha(1)=(hbar(2)-hbar(1))/del_x;
alpha(N)=(hbar(N)-hbar(N-1))/(del_x);
beta(1)=-(vbar(1)/hbar(1))*((hbar(2)-hbar(1))/(del_x));
beta(N)=-(vbar(N)/hbar(N))*((hbar(N)-hbar(N-1))/(del_x));
gema(1)=2*g*(m^2)*(vbar(1)/(hbar(1)^(4/3)))-(vbar(1)/hbar(1))*((hbar(2)-hbar(1))/(del_x));
gema(N)=2*g*(m^2)*(vbar(N)/(hbar(N)^(4/3)))-(vbar(N)/hbar(N))*((hbar(N)-hbar(N-1))/(del_x));
eta(1)=-(4/3)*g*(m^2)*((vbar(1)*abs(vbar(1)))/(hbar(1)^(7/3)));
eta(N)=-(4/3)*g*(m^2)*((vbar(N)*abs(vbar(N)))/(hbar(N)^(7/3)));
w=2;
for i=2:N-1
alpha(i)=(hbar(w+1)-hbar(w-1))/(2*del_x);
beta(i)=(-vbar(w)/hbar(w))*((hbar(w+1)-hbar(w-1))/(2*del_x));
gema(i)=2*g*(m^2)*(vbar(w)/(hbar(w)^(4/3)))-(-vbar(w)/hbar(w))*((hbar(w+1)-hbar(w-1))/(2*del_x));
eta(i)=-(4/3)*g*(m^2)*((vbar(w)*abs(vbar(w)))/(hbar(w)^(7/3)));
w=w+1;
end
vshift=1;
hshift=1;
A(1,:)=[0 (1/2)-(del_t/(4*del_x))*(vbar(3)+vbar(1))-((del_t/2)*gema(3)) zeros(1,N-3) (-g*(del_t/(2*del_x))-(del_t/2)*eta(3)) zeros(1,N-4)];
A(A_size/2,:)=[zeros(1,N-4) (1/2)+(del_t/(4*del_x))*(vbar(N)-vbar(N-2))-(del_t/2)*gema(N-2) zeros(1,N-3) (g*(del_t/(2*del_x))-(del_t/2)*eta(N)) 0];
w=3;
for i=2:1:(A_size/2)-1    
A(i,:)=[(1/2)+(del_t/(4*del_x))*(vbar(w+1)-vbar(w-1))-(del_t/2)*gema(w-1) 0 (1/2)-(del_t/(4*del_x))*(vbar(w+1)+vbar(w-1))-(del_t/2)*gema(w+1) zeros(1,N-5) (g*(del_t/(2*del_x))-(del_t/2)*eta(w-1)) 0 (-g*(del_t/(2*del_x))-(del_t/2)*eta(w+1)) zeros(1,N-5)];
if i>2
A(i,:)=circshift(A(i,:),[i,vshift]);
vshift=vshift+1;
end
w=w+1;
end
A((A_size/2)+1,:)=[0 (-del_t/(4*del_x))*(hbar(3)+hbar(1))-(del_t/2)*alpha(3) zeros(1,N-3) ((1/2)-(del_t/(4*del_x)*(vbar(3)+vbar(1)))-(del_t/2)*beta(3)) zeros(1,N-4)];
hrow=(A_size/2)+2;
for i=3:1:(A_size/2)
A(hrow,:)=[(del_t/(4*del_x))*(hbar(i+1)+hbar(i-1))-(del_t/2)*alpha(i-1) 0 (-del_t/(4*del_x))*(hbar(i+1)+hbar(i-1))-(del_t/2)*alpha(i+1) zeros(1,N-5) ((1/2)+(del_t/(4*del_x)*(vbar(i+1)+vbar(i-1)))-(del_t/2)*beta(i-1)) 0 ((1/2)-(del_t/(4*del_x)*(vbar(i+1)+vbar(i-1)))-(del_t/2)*beta(i+1)) zeros(1,N-5)];
if hrow>(A_size/2)+2
A(hrow,:)=circshift(A(hrow,:),[hrow,hshift]);
hshift=hshift+1;
end
hrow=hrow+1;
end
A(A_size,:)=[zeros(1,N-4) ((del_t/(4*del_x))*(hbar(N)+hbar(N-2)))-((del_t/2)*alpha(N-2)) zeros(1,N-3) (1/2)+((del_t/(4*del_x)*(vbar(N)+vbar(N-2))))-((del_t/2)*beta(N-2)) 0];
B=zeros(A_size,B_size);
 B(1,:)=[(1/2)+(del_t/(4*del_x))*(vbar(3)-vbar(1))-(del_t/2)*gema(1) 0 (g*(del_t/(2*del_x)))-(del_t/2)*eta(1) 0];
 B(A_size/2,:)=[0 (1/2)-((del_t/(4*del_x))*(vbar(N)+vbar(N-2)))-((del_t/2)*gema(N)) 0 (-g*(del_t/(2*del_x))-(del_t/2)*eta(N))];
 B((A_size/2)+1,:)=[((del_t/(4*del_x))*(hbar(3)+hbar(1)))-((del_t/2)*alpha(1)) 0 (1/2)+(del_t/(4*del_x))*(vbar(3)+vbar(1))-(del_t/2)*beta(1) 0];
 B(A_size,:)=[0 ((-del_t/(4*del_x))*(hbar(N)+hbar(N-2)))-((del_t/2)*alpha(N)) 0 (1/2)-((del_t/(4*del_x))*(vbar(N)+vbar(N-2)))-((del_t/2)*beta(N))];


