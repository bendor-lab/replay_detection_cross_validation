function score=fast_pacman(decoded_event)

%https://blog.usejournal.com/dynamic-programming-finding-optimal-paths-2cf89f3d1437

[score1,maxsum1]=maxsum_algorithm(decoded_event); %start from top and go to bottom right
[score2,maxsum2]=maxsum_algorithm(flipud(decoded_event)); %start from bottom and go to top right


if score1>score2;
    score=score1;
    maxsum=maxsum1;
else
    score=score2;
     maxsum=maxsum2;
end  

score=score/size(decoded_event,2);  %normalise by number of time bins

% score=score/(size(decoded_event,2)+sum(sum(decoded_event)));  
% normalise by the sum of all posterior probabilities
 
% score=score/sum(sum(decoded_event));  
% normalise by the sum of all posterior probabilities
 
end

function [score, maxsum]=maxsum_algorithm(matrix)
maxsum=zeros(size(matrix));
 %cumulative sum of first column and row
maxsum(1,:)=cumsum(matrix(1,:)); 
maxsum(:,1)=cumsum(matrix(:,1));


%maxsum matrix is calculated by taking the maximum between the previous
%node (prior column and prior row), and adding it to the value of the
%current matrix node
for i=2:size(maxsum,1)
	for j=2:size(maxsum,2)
		maxsum(i,j)=matrix(i,j)+max([maxsum(i-1,j) maxsum(i,j-1)]);
    end
end

%score of optimal path is the final node (last column and last row)
score=maxsum(size(maxsum,1),size(maxsum,2));

end
