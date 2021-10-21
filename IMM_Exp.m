%%---- SD-IMM code to estimate water levels and water velocities using
%%upstream and downstream water level and Lagrangian GPS positions as
%%measurements
%%Coded by Affan (email: affan.005@outlook.com)%%%%%%%%%%%%%%%%%%%%%%%

%Reading GPS data from Lagrangian sensor
fileID = fopen('position2.txt');
C = textscan(fileID,'%f');
fclose(fileID);
whos C;
pos=[C{1}];
%Reading Sonar sensor data at 1km range.  
fileID = fopen('1KM_rangedata.txt');
C = textscan(fileID,'%f');
fclose(fileID);
whos C;
h1km=flipud([C{1}]);
h1km=h1km./1000; %Converting to meters from milimeters
%Reading sonar sensor data at 2km  
fileID = fopen('2KM_rangedata.txt');
C = textscan(fileID,'%f');
fclose(fileID);
whos C;
h2km=flipud([C{1}]);
[h2km]=interpolation(h2km);
h2km=h2km./1000;
del_x=270; %Cell size
[y1,yn,v1,vn]=boundary_cond(); %Computing boundary conditions. 

%making water level same at t=0. 
h2km=4.78-h2km;
h1km=4.8-h1km;
% plotting water level sensor data
figure(1)
plot(y1,'c')
hold on
plot(h1km,'g')
plot(h2km,'r')
plot(yn,'b')
title('Water Level Sensor Measurement Data.')
ylabel('m')
xlabel('Time (minutes)')
legend('Upstream', '1 Km','2 Km','Downstream')
legend box off
grid on
%%
%Definig parameters to compute matrix A and B 
time=length(pos);
L=pos(end); %Channel Length in meters
t=1; %time step between two cell for discretization 
%Initial values for steady state back water curve
Vo=v1(1);
Ho=y1(1);
% Defining system
[A,B,A_size,n,hbar,vbar]=matrix2(L,t,del_x,Vo,Ho);
len_A=length(A);
cell_1_km=ceil(1000/del_x); %index in vector
cell_2_km=ceil(2000/del_x);
%Defining System matrix H for upstream water level, downstream water level
%and lagrangian sensor.
H(1,:)=[0 zeros(1,(A_size/2)) 1 zeros(1,(A_size/2)-1)];             
H(2,:)=[1 zeros(1,(A_size))];
H(3,:)=[zeros(1,(A_size)) 1];
%Initial states
X=[pos(1) repmat(Vo,1,10) repmat(Ho,1,10)]';
%Defining state transition matrix 
Pt=n;
T=conv2(eye(Pt),[0.1 0.7 0.2],'same');
T(1,2)=0.5;
T(1,1)=0.5;
T(n,n-1)=0.5;
T(n,n)=0.5;
T=T';
% Definig ui, probability that target is in the state i as computed just
% after the data is received 
ui=[0.75 0.15 0.01 0.01 0.01 0.01 0.01 0.01 0.01 0.01 0.01 0.01]';
%Conditional Probability given that the tragte is in state j that the
%transition occured from state i
uij=zeros(Pt,Pt);
%Defining Kalman filter coveriance matrixes 
%tune these values for better estimation results. 
 q_var=1.0000e-03;
 p_var=0.0032;
 r_var=3.1623e-04;
 Q=(q_var^2)*eye(A_size+1);
 P=(p_var^2)*eye(A_size+1);
P(1,1)=10;% for lagrangian sensor
R=(r_var^2)*eye(3);

%Creating n state vectors 
for i=1:n
    Xi(:,i)=X;
    Pi{i}=P;
end
state=X;

%% Defining IMM 
for k=1:time      %loop on time. 
    if k<40
        R(2,2)=200;
    else
        R(2,2)=100;
    end
    output(:,k)=H*state;
    State(:,k)=state;
    %measurement vector for Kalman Filter
    z=[y1(k) pos(k) yn(k)]; 
    %z=[h2km(k) pos(k)];
    %Input vector for System 
    u=[v1(k) vn(k) y1(k) yn(k)];
    %% IMM Mixing 
    %Computing up
    up=T*ui;
    %Computing uij
    for i=1:n
        for j=1:n
            uij(i,j)=(T(j,i)*ui(i))/up(j);
        end
    end
    for cell=1:n
        %Computing Xj
        Xtemp=[];
        for i=1:n
            Xtemp(:,i)=Xi(:,i)*uij(i,cell);
        end
        Xj(:,cell)=sum(Xtemp,2);
        %Computing Pj
        Dx=0;
        for i=1:n
            Dx=Dx+((Xi(:,i)-Xj(:,cell))*(Xi(:,i)-Xj(:,cell))')*uij(i,cell);
            Ptemp{i}=cell2mat(Pi(i))*uij(i,cell); 
        end
        Pj{cell}=0;
        for i=1:n
            Pj{cell}=cell2mat(Pj(cell))+cell2mat(Ptemp(i));        
        end
        Pj{cell}=cell2mat(Pj(cell))+Dx;
    end
    %% Step-2 Run Individual Kalman Filter for each model
    %Defining Augmented models 
    for cell=1:n
        if cell==1
            A_aug=[1 zeros(1,(len_A));
                  zeros(len_A,1) A];
            B_aug=[60 zeros(1,3); B];
        elseif cell==n
            A_aug=[1 zeros(1,(len_A));
                  zeros(len_A,1) A];
            B_aug=[0 60 zeros(1,2); B] ;
        else
            A_aug=[1 zeros(1,cell-2) 60 zeros(1,(len_A)-(cell)+1);
                  zeros(len_A,1) A]; 
            B_aug=[zeros(1,4); B];
        end
        %Applying Kalman filter      
        if pos(k)~=-1  %Condition to cehck if GPS data is available or not
            [Xj_temp, Pj_temp] = predict(Xj(:,cell), Pj(cell), A_aug, Q, B_aug, u);        %State Prediction
            [nu_temp, S_temp] = innovation(Xj_temp, Pj_temp, z, H, R);               %Computing Innovaton/error
            [Xj_temp, Pj_temp,c_out] = innovation_update(Xj_temp, Pj_temp, nu_temp, S_temp, H);  %Updating states by using kalman gain and innovation
            
            Pj{cell}=Pj_temp; Xj(:,cell)=Xj_temp; nu(:,cell)=nu_temp; S{cell}=S_temp; output(:,cell)=c_out;
            
            %computing statistical distance of an obervation-to-track assignment
            dtemp(cell)=nu(:,cell)'*inv(cell2mat(S(cell)))*nu(:,cell);
            sigma(cell)=exp(-dtemp(cell)/2)/sqrt((2*pi)*det(cell2mat(S(cell))));
        else
            [Xj_temp, Pj_temp] = predict(Xj(:,cell), Pj(cell), A_aug, Q, B_aug, u);        %State Prediction
            Pj{cell}=Pj_temp; Xj(:,cell)=Xj_temp;
        end
    end
    Xi=Xj;
    Pi=Pj;
    %% Step-3 Remix
    %Updating Probabilities
    if pos(k)~=-1
        Ctemp=sigma*up;
        for i=1:n
            ui(i)=(sigma(i)*up(i))/Ctemp; 
        end
        P_temp=[];
        %Combining State vector and Covariance matrix from all models to display output.  
        for i=1:n
            St_temp(:,i)=ui(i)*Xj(:,i);
            %P_temp{i}=up(i)*cell2mat(Pj(i));
        end
        state=sum(St_temp,2);   
        %Computing Pj
        PDx=0;
        for i=1:n
            PDx=PDx+((state-Xj(:,i))*(state-Xj(:,i))')*ui(i);
            PTemp{i}=cell2mat(Pj(i))*ui(i); 
        end
        Pnew=0;
        for i=1:n
            Pnew=cell2mat(Pj(i))+cell2mat(PTemp(i));        
        end
        Pnew=Pnew+PDx;
    end     
end
IMM_state=State;
%% Plotting
%load 'D:\Water Dynamics\MATLAB Files\KF_newest.mat'
%Calculating error
h1km_error=immse(h1km,State((A_size/2)+cell_1_km+1,:)');
h2km_error=immse(h2km,State((A_size/2)+cell_2_km+1,:)');
figure(2)
%subplot(211)
hold on
plot(State((A_size/2)+cell_1_km+1,:),'k')
plot(State((A_size/2)+cell_2_km+1,:),'r')
plot(h1km,'m')
plot(h2km,'g')
title('Water Level')
ylabel('m')
xlabel('Time (s)')
legend('1Km IMM','2Km IMM','1Km Actual','2Km Actual')
legend box off

%subplot(212)
figure(3)
hold on
plot(State(cell_1_km,:),'k')
plot(State(cell_2_km,:),'r')
title('Water Velocity')
ylabel('m/s')
xlabel('Time (s)')
legend('1Km IMM','2Km IMM')
legend box off
for i=1:length(pos)
    if pos(i)==-1
        pos(i)=NaN;
    end
end
figure(4)
hold on
plot(State(1,:),'g')
plot(pos,'c-*')
title('Lagrangian Sensor Trajectory')
ylabel('m')
xlabel('Time (minutes)')
legend('Estimate','Actual')
grid on
legend box off





%%
%Defining Kalman Filter Functions 
function [Xpred, Ppred]=predict(x,P,F,Q,B,u)
P=cell2mat(P);
Xpred=F*x+B*u';
Ppred=F*P*F'+Q;
end
function [nu,S]=innovation(Xpred,Ppred,z,H,R)
nu=z'-(H*Xpred);
S=R'+(H*Ppred*H');
end
function [Xnew,Pnew,c_out]=innovation_update(Xpred,Ppred,nu,s,H)
K=Ppred*(H'*inv(s));
Xnew=Xpred+(K*nu);
c_out=H*Xnew;
Pnew=Ppred-(K*s*K');
end