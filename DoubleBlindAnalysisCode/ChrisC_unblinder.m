
function Final=ChrisC_unblinder(ChrisC,IndexMatrix)

    Final=ChrisC(IndexMatrix(:,3)');

    for i=1:size(IndexMatrix,1)
        Final(i).CellIndex=IndexMatrix(i,1);
        Final(i).track=IndexMatrix(i,2);
    end

end