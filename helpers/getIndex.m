function [k,l2,j,i,sameCell] = getIndex(l, numS, numTheta)
    quotient = floor((l-1)/(numTheta * numS));
    k = floor(quotient/numTheta)+1; % Integration segment s
    l2 = rem(quotient,numTheta)+1; % Integration segment theta
    
    remainder = rem((l-1),numTheta * numS);
    j = floor(remainder/numTheta)+1; % Reference segment s
    i = rem(remainder,numTheta)+1; % Reference segment theta

    sameCell = (j == k) && (l2 == i);
end