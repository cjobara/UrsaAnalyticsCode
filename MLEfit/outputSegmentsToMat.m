load('ChrisC')
N=25; %make N the same as number of xxx.h5 files (where xxx is an integer without zero padding) in the same folder
for i=1:N;d=h5read([num2str(i) '.h5']);ChrisC(i).cp=d(:,end);end
save('ChrisC_CP.mat','ChrisC')