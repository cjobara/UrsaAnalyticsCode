


for i=1:size(Tracks,2)
    A(1:size(Tracks(i).lengths),i)=Tracks(i).lengths;
end


B=A>500; %Set the minimum step threshold here
f=sum(B,1);
[a,b,h]=find(B.*A);
g=NaN(max(f,[],'all'),size(Tracks,2)); %records track indexes
e=NaN(max(f,[],'all'),size(Tracks,2)); %records track lengths

for i=1:size(Tracks,2)

    figure
    subplot(1,2,1)
     filebase=Tracks(i).file;
    [~,L]=size(filebase);
    filename=strcat(filebase(1:L-18),'2_maxS2N_max_RGB.tif');
    imS=ChrisPrograms.loadtiff(fullfile('E:\Part2\JNA6-018\HighQuality\Channel2\Endosomes\MaxInt\RGB',filename));
    imshow(imS,'Border','tight');
    hold on
    plot(9.375*Tracks(i).matrix(:,:,2),9.375*Tracks(i).matrix(:,:,3),'Color','b')
    c=find(b==i);
    d=a(c);
    plot(9.375*Tracks(i).matrix(:,d',2),9.375*Tracks(i).matrix(:,d',3),'Color','g');
    e(1:size(d),i)=Tracks(i).lengths(d);
    g(1:size(d),i)=d;
    
    subplot(1,2,2)
   
    imshow(imS,'Border','tight');
    hold on
    plot(9.375*Tracks(i).matrix(:,d',2),9.375*Tracks(i).matrix(:,d',3));
    
end

%Notice this is not the same as the previous 'B' which was just a dummy
%variable so I reused it

B=[b a h];

% ususally you manually save 'B' as a .mat file named for the condition

Long=NaN(max(e,[],'all'),sum(f),3); %sum(f) = size(B,1) just FYI

for i=1:size(B,1)
    
    Long(1:B(i,3),i,:)=Tracks(B(i,1)).matrix(1:B(i,3),B(i,2),:);
    
end

OutFile='TempFile.xlsx';
xlswrite(OutFile,B);