


for i=1:size(Tracks,2)
    A(1:size(Tracks(i).lengths),i)=Tracks(i).lengths;
end


B=A>300; %Set the minimum step threshold here
f=sum(B,1);
[a,b]=find(B);
g=NaN(max(f,[],'all'),size(Tracks,2));
e=NaN(max(f,[],'all'),size(Tracks,2));

for i=1:size(Tracks,2)

    figure
    subplot(1,2,1)
     filebase=Tracks(i).file;
    [~,L]=size(filebase);
    filename=strcat(filebase(1:L-25),'2_TA_BC.tif');
    imS=ChrisPrograms.loadtiff(fullfile('E:\Cris\Cris-003\HighQuality\ER',filename));
    imshow(imS,'Border','tight');
    %hold on
    plot(6.25*Tracks(i).matrix(:,:,2),6.25*Tracks(i).matrix(:,:,3),'Color','b')
    c=find(b==i);
    d=a(c);
    plot(6.25*Tracks(i).matrix(:,d',2),6.25*Tracks(i).matrix(:,d',3),'Color',[ 0.9100 0.4100 0.1700]);
    e(1:size(d),i)=Tracks(i).lengths(d);
    g(1:size(d),i)=d;
    
    subplot(1,2,2)
   
    %imshow(imS,'Border','tight');
    %hold on
    plot(6.25*Tracks(i).matrix(:,d',2),6.25*Tracks(i).matrix(:,d',3));
    
end