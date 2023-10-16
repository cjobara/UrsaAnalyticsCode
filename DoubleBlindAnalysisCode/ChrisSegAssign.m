
function [ChrisCseg,SegTable]=ChrisSegAssign(ChrisC,SegmentDetails)

    % SegTable is a table of the format of SegmentDetails that contains all of the
    % separated SegIDs as well as the specific track information for each
    % segment.
    
    % ChrisCseg is a structured array of the format of ChrisC that adds the segment
    % details to the input ChrisC so all of the data is accumulated in once
    % place. It also provides linear indexing coordinates to all of the
    % change points, which are easier to work with than the logical arrays
    % in the original structure.
    
    % Note, if you are looking to translate to M and N space, this is
    % performed by another script called after this one that unblinds the
    % ChrisC indexing to true track and cell values.
    
    % A note on definitions:segIDs refer to what number a segment is in the
    % segments result file. segNum refers to what number a segment is in
    % its own track. "TempSegNum" refers to the seg number assigned by
    % ChrisC, it is out of date so don't use it because it doesn't account
    % for data he couldn't fit.
    
    % Define indexing variables for ChrisC struct
    
    for i=1:size(ChrisC,2)
        
        %Define linear indexes and identities of segments in each track
        ChrisC(i).ChPts=[1 transpose(find(ChrisC(i).cp)) size(ChrisC(i).cp,1)];
        ChrisC(i).lengths=diff(ChrisC(i).ChPts);
        ChrisC(i).segIDbySeg=zeros(size(ChrisC(i).lengths));
        %Allocate space to store which IDs each step belongs to for
        %explicit indexing later
        ChrisC(i).segNum=NaN(size(ChrisC(i).cp));
        ChrisC(i).segID=zeros(size(ChrisC(i).cp));
        for j=1:size(ChrisC(i).lengths,2)
           ChrisC(i).segNum(ChrisC(i).ChPts(j):ChrisC(i).ChPts(j+1)-1)= j*ones(ChrisC(i).lengths(j),1);
        end
        
    end
    
    % Run through SegDetails and find the track and segment number for each
    % analyzed track
    
    [A,B,C,D,E,F,G]=ChrisCconvert(SegmentDetails);
    F=permute(F,[3,1,2]);
    G=permute(G,[3,1,2]);
    
    CCindex=NaN(size(A));
    TempSegNum=NaN(size(A));
    SegNumList=zeros(size(A));
    
    for i=1:size(A,1)
        
       SegStrng=char(A{i}); 
       CCindex(i)=str2num(SegStrng(1:strfind(SegStrng,'_')-1));
       TempSegNum(i)=str2num(SegStrng(strfind(SegStrng,'_')+1:end));
       if TempSegNum(i)==1
           DoubleFlag=0;
       end
       
       if numel(find(ChrisC(CCindex(i)).lengths==B(i)))==1
            SegNumList(i)=find(ChrisC(CCindex(i)).lengths==B(i));
       elseif numel(find(ChrisC(CCindex(i)).lengths-1==B(i)))==1
            SegNumList(i)=find(ChrisC(CCindex(i)).lengths-1==B(i));
       elseif numel(find(ChrisC(CCindex(i)).lengths+1==B(i)))==1
            SegNumList(i)=find(ChrisC(CCindex(i)).lengths+1==B(i));
       else
           % Use these for finding doubles
%            disp(strcat('i=',num2str(i)));
%            disp(ChrisC(CCindex(i)).lengths);
%            disp(strcat('B=',num2str(B(i))));
%            SegNumList(i)=NaN;

           % Use this section for matching doubles you know about
           Options=find((ChrisC(CCindex(i)).lengths==B(i))+(ChrisC(CCindex(i)).lengths+1==B(i))+(ChrisC(CCindex(i)).lengths-1==B(i)));
           if DoubleFlag==0
               SegNumList(i)=Options(1);
               DoubleFlag=1;
           elseif size(Options,2)==2
               SegNumList(i)=Options(2);
           else
               error(strcat('There are three segments of same length? i=',num2str(i)));
           end
           
       end
       % Mark which segments in ChrisC have analysis results associated
       ChrisC(CCindex(i)).segIDbySeg(SegNumList(i))=i;
       ChrisC(CCindex(i)).segID(ChrisC(CCindex(i)).ChPts(SegNumList(i)):ChrisC(CCindex(i)).ChPts(SegNumList(i)+1))=i;
        
    end

    % Prepare Outputs    
    ChrisCseg=ChrisC;
    T=table(CCindex,SegNumList,B,C,D,E,F,G);
    T=renamevars(T,{'CCindex','SegNumList','B','C','D','E','F','G'},{'ChrisC index','Segment Number','NumSteps','D to F ratio','F Evals', 'D evals', 'F evecs', 'D evecs'});
    SegTable=T;
    
    % Backup the Data
    writetable(T,'ChrisCSegDetails.xlsx');
    save('ChrisCwithSegInfo.mat','ChrisC','T');
end