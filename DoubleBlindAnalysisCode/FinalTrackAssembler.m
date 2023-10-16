
function TrackFinal=FinalTrackAssembler(TrackStruct,NPBstruct)

    TrackFinal=TrackStruct;
    
    for i=1:size(TrackStruct,2)
        
       TrackFinal(i).CCindex=NaN(size(TrackStruct(i).matrix(:,:,1)));
       TrackFinal(i).segNum=NaN(size(TrackStruct(i).matrix(:,:,1)));
       TrackFinal(i).segID=NaN(size(TrackStruct(i).matrix(:,:,1)));
       TrackFinal(i).cp=zeros(size(TrackStruct(i).matrix(:,:,1)));     
        
    end
        
    for i=1:size(NPBstruct,2)
        TrackLength=size(NPBstruct(i).x,1);
        if TrackLength==sum(isfinite(TrackStruct(NPBstruct(i).CellIndex).matrix(:,NPBstruct(i).track,2)),'all','omitnan')
            TrackFinal(NPBstruct(i).CellIndex).CCindex(1:TrackLength,NPBstruct(i).track)=NPBstruct(i).index*ones(TrackLength,1);
            TrackFinal(NPBstruct(i).CellIndex).segNum(1:TrackLength,NPBstruct(i).track)=NPBstruct(i).segNum;
            TrackFinal(NPBstruct(i).CellIndex).segID(1:TrackLength,NPBstruct(i).track)=NPBstruct(i).segID;
            TrackFinal(NPBstruct(i).CellIndex).cp(1:TrackLength,NPBstruct(i).track)=NPBstruct(i).cp;
           
        else
            error(strcat('NPB track number: ', num2str(i),'does not match'));
        end
        
    end



end