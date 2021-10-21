function [act_data]=interpolation(raw_data)
% fileID = fopen('data.txt');
% C = textscan(fileID,'%f');
% fclose(fileID);
% whos C;
% out=[C{1}] ;
% raw_data=flipud(out);
act_data=raw_data;
j=0;
for i=2:length(raw_data) % loop for data size
    if j==1
        xo=i-2;          %Getting first point for interpolation
    end
    if raw_data(i)==0 || isnan(raw_data(i))==1 %Condition to check for 0 or Nan Values
        j=j+1;
    %elseif raw_data(i)>7  %Condition to check for unexpected higher values
     %   j=j+1;
    elseif raw_data(i)==300   %Condition to check for unexpected lower values
        j=j+1;      
    else
            if j>0    %Code for lagrange interpolation
            x1=i;
            for k=1:j
                x=xo+k;
                Lo=(x-x1)/(xo-x1);
                L1=(x-xo)/(x1-xo);
                act_data(x)=Lo*raw_data(xo)+L1*raw_data(x1);
            end
            j=0;
            end
    end
end
%Plotting 
% plot(raw_data,'ro')
% hold on
% plot(act_data,'go-')
% hold off
% str = sprintf(' Jamrao Branch Data');
% title(str)
% ylabel('Water Level (ft)')
% xlabel('Time index')
% xlim auto
% ylim auto
% legend('Raw Data', 'Interpolated Data')



