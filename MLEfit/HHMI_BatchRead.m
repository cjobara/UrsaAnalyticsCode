function [X,Y,z,T,d]=HHMI_BatchRead(folderPathandName,casei);
% simple routine for reading names .mat files from Nov 6 delivery (see notes). HDP-SLDS in matlab
% unpacks the .mat provided and massages
%
% folderPathandName:  string like ../datafolder/usethisfolder   [don't put os file sep as last entry of input...automatically added]


load([folderPathandName  filesep  'ChrisC.mat']) %make more general later...for now, struc is hard coded



d=2;  %just set dimension in this routine to make code agnostic to hard coded info in main routine (but app specific stuff in routine like this)


i=casei;ind=isfinite(ChrisC(i).t);



X=[ChrisC(i).x(ind) ChrisC(i).y(ind)];  %put into microns (units prior is in)
for dimi =1:d
	X(:,dimi)=X(:,dimi)-mean(X(:,dimi));
end

Y=X;
z=ones(1,length(X));
% T=length(X);
T=ChrisC(i).t(ind); %if some have NaN, problem
errorCheck=max(isnan(T))
if(errorCheck>0)
    display('error')
    pause
end
