function [TCMPercentileM,Percentile]...
    = MobilityPDFMOB_v1(GSDlength,taucrit,taucritless,PsiMIN,PsiMAX,R2GSeval,VRandom)

%% THIS FUNCTION COMPUTES PROBABILITY DISTRIBUTIONS OF MOBILITY DEFINED BY
% one-half PSI class intervals.  The results are used to compute sediment
% transport.

    % Provide an increment for the Psi class interval
    PsiIncrement = 0.5;
    % Define the length of the Psi class array based on the inpurt GSD
    PsiArrayLength = ((PsiMAX - PsiMIN) / PsiIncrement) + 1;
    % Use the logical indexing to build the Psi class array
    j = 1:PsiArrayLength;
    PsiArray(j) = PsiMIN:PsiIncrement:PsiMAX;
    % Convert the random grains used to compute mobility to Psi units
    R2GSPsi = log2(R2GSeval);
    
    MobilityMatrix = [R2GSPsi,taucrit,taucritless];
    MM = sortrows(MobilityMatrix,1);
    taucritlessarray = MM(:,3);
    taucritlessmatrix = repmat(taucritlessarray,1,GSDlength);
    taucritlessmatrix(isnan(taucritlessmatrix)) = 0;
    
    PsiMatrix(VRandom,PsiArrayLength-1) = 0;
    
    for j = 1:PsiArrayLength - 1
               
       if j == 1
           
           PsiMatrix(:,j) = MM(:,1) < PsiArray(j+1); 
           
       elseif j > 1 && j < PsiArrayLength - 1
           
          PsiMatrix(:,j) = MM(:,1)>= PsiArray(j) & MM(:,1) < PsiArray(j+1); 
        
       else
           
          PsiMatrix(:,j) =  MM(:,1) > PsiArray(j); 
          
       end
       
    end
    
    PsiMatrix(all(PsiMatrix==0,2),:)=[];
    taucritlessmatrix(all(taucritlessmatrix==0,2),:)=[];
    
    % Now build a taucritless matrix that is filtered by Psi class.  This
    % result will be used to build PDFs of the taucritless by Psi class.
    taucritlessmatrix(~PsiMatrix)=0;
    % Transpose the matrix so it lines up the the terms of GSD matrix
    TCM = taucritlessmatrix';
    % Sort each row of the TCM matrix to make the statistics easier
    TCMrowsort = sort(TCM,2,'descend');
    % Replace zeros with NaN so the statistics compute correctly
    TCMrowsort(TCMrowsort == 0) = NaN;
    
    % Now compute percentile statistics of dimensionless critical Shields
    % stress for each Psi class.  This is made easy with the matlab prctile
    % function.  Specify the percentile classes to compute and then compute
    % the percenitle results. The last number in parentheses tells matlab to 
    % operate on the rows, not the columns
    Percentile = [10 20 30 40 50 60 70 80 90];
    TCMPercentileM = prctile(TCMrowsort,Percentile,2);
       
end

