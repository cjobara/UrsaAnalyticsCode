
floc='./MaySegmentation/' %Reference the top level HDP-SLDS output dir here


%%%%%%%%%%%%%%%%%%%%%%%%%
%create python like crude python dict objects to cycle through cases of interest (do not include abs or relative path in keys, just filename)
dataFiles = containers.Map



 for i=1:401;
   dataFiles(num2str(i)) = [floc filesep num2str(i) '/changePntLocsMU_trial_1.mat'];
end

%%%%%%%%%%%%%%%%%%%%%%%%%


%nor march through all files above and spit out hdf5 files for further processing
myKeys=keys(dataFiles)


frac=.4;%parameter determining how often CP needs to show up in MCMC sampler output to count to be recorded
istart=50;%determine point for burn in of MCMC sampler (code setup to output cp ever 100 iterations, so setting this to 10 results in waiting 1000 MCMC steps)

for fi = 1:length(myKeys)


	dataFile=[num2str(fi) '.mat'];
    outFile=dataFiles(num2str(fi));
    



	% g=load(dataFile);
  [X,YY,z,T,d]=HHMI_BatchRead('./Obara_May10_2020/',fi); %this block is repeated later;

	
	
	load(outFile);
	meanData=mean(YY)*0;
    %scaleFac=.16; %to have in microns, set to 1
    scaleFac=1.; %to have in microns, set to 1
	gY=(YY-repmat(meanData,length(YY),1))/scaleFac;
	

	
	% y=mean(changePntHist);
	y=mean(changePntHist(istart:end,:));
	l=length(y);x=[1:l];ym=[0 0 y];yp=[y 0 0];Y=(6*y+2*yp(2:end-1)+1*yp(1:end-2)+2*ym(2:end-1)+1*ym(1:end-2))/12;
	
	
	% dz=y.*(y>frac)./y;
    dz=(y>frac);
%time series plot%%%%%%%%%%%%%%%%%%%%%    
% figure;[haxes,hline1,hline2] = plotyy(x,gY(1:l,:),x,dz);
trim=5;
figure;[haxes,hline1,hline2] = plotyy(x(1:end-trim),gY(1:l-trim,:),x(1:end-trim),dz(1:end-trim));
      %title(myKeys{fi})
  	set(haxes(1),'LineWidth',2);set(hline1,'LineWidth',2);set(hline2,'LineWidth',2);xlabel('Obs #','fontsize',20);
    set(gca,'fontsize',20);ylabel(gca,'X/Y [microns]','fontsize',20);%title(myKeys{fi})
    set(gca,'xlim',[0,length(dz)]);
    dhi=floor(max([gY(1:l,1);gY(1:l,2)])); dlo=floor(min([gY(1:l,1);gY(1:l,2)]));dd=ceil((dhi-dlo)/10);yrange=[dlo:dd:dhi];
    set(gca,'ylim',[dlo-2,dhi+2]);
    set(gca,'ytick',yrange)
    set(haxes(2),'xlim',[0,length(dz)])
    set(haxes(1),'xlim',[0,length(dz)])
    %pseVfunc([num2str(fi) '_traj'])
%%%%%%%%%%%%%%%%%%%%% 
dz=[0 dz 0]; %only modify after plot above is done
    
%phase portrait plot%%%%%%%%%%%%%%%%%%%%     
    figure;plot((gY(1:l,1)+meanData(1))/1,(gY(1:l,2)+meanData(2))/1,'LineWidth',2);hold on;plot((gY(dz>0,1)+meanData(1))/1,(gY(dz>0,2)+meanData(2))/1,'ro','Markersize',12,'LineWidth',2); 
    set(gca,'fontsize',20);xlabel(gca,'X [microns]','fontsize',20);ylabel(gca,'Y [microns]','fontsize',20);
     XY=[gY(:,1)+meanData(1) gY(:,2)+meanData(2)];changePoint=dz>0;
     save([floc myKeys{fi}(18:end) '_traj_and_ChangePoints.mat'],'XY','changePoint');
%%%%%%%%%%%%%%%%%%%%%  
	
	dataMat=[T T*1/100 YY];
	dataMat=[dataMat dz'];
	
	h5out(dataFile(1:end-4),dataMat) %only save if needed
    close all

end

