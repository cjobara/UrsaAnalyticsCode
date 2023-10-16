function [X,Y,z,T,d,csvName]=BatchRead(folderPathandName,casei);
% sample routine for reading names of all csv files in "folder" and returning info needed to run HDP-SLDS routines in matlab
% temporary hacky solution.  ulitmately should create python version 
%
% folderPathandName:  string like ../datafolder/usethisfolder   [don't put os file sep as last entry of input...automatically added]

%pull file structure list
fileStruc=dir([folderPathandName  filesep '*.csv']);

%expected input format of csv
% column 1.  linear observation index
% column 2.  time spacing (some highly nonuniform)
% column 3-5 x,y,z [nm]
% column 6-8 std x,y,z [nm] (optional)

d=3;  %hacky hard coded spatial dimenesion of measurement.  adjust accordingly for 1d,2d or 3d models

csvName = [folderPathandName filesep fileStruc(casei).name];
m=csvread(csvName);

X=m(:,3:3+d-1)/1000;  %put into microns (units prior is in)
for i =1:d
	X(:,i)=X(:,i)-mean(X(:,i));
end

Y=X;
z=ones(1,length(X));
T=length(X);






