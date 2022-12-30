function [path_score,best_path]=pacman(decoded_event)

psize = size(decoded_event);
[forward_score,forward_path] = find_path(decoded_event);
[reverse_score,reverse_path] = find_path(flipud(decoded_event));

%compare the upper-left lower-right score to the lower-left
%upper-right score, take the higher one
if forward_score > reverse_score
    path_score = forward_score / psize(2);  %divde by number of time bins in path
    best_path = forward_path;
else
    path_score = reverse_score / psize(2);
    reverse_path(:,1) = psize(1) + 1 - reverse_path(:,1);
    best_path = reverse_path;
end

end


function [score,path] = find_path(decoded_event)
%FIND_PATH returns the maximum score and the corresponding path of a frame
%from the upper left corner [1,1] to the lower right corner
    psize = size(decoded_event);
    score_map = zeros(psize);
    score_map(1,1) = decoded_event(1,1);
    path_map = cell(psize);
    path_map{1,1} = [1,1];
    for i = 2 : psize(1)
        score_map(i,1) = score_map(i-1,1) + decoded_event(i,1);
        path_map{i,1} = [path_map{i-1,1}; [i,1]];
    end
    for i = 2 : psize(2)
        score_map(1,i) = score_map(1,i-1) + decoded_event(1,i);
        path_map{1,i} = [path_map{1,i-1}; [1,i]];
    end
    for i = 2 : min(psize)
        for j = i : psize(1)
            if score_map(j-1,i) > score_map(j,i-1)
                score_map(j,i) = score_map(j-1,i) + decoded_event(j,i);
                path_map{j,i} = [path_map{j-1,i}; [j,i]];
            else
                score_map(j,i) = score_map(j,i-1) + decoded_event(j,i);
                path_map{j,i} = [path_map{j,i-1}; [j,i]];
            end
        end
        if i == psize(2)
            break;
        end
        for j = i+1 : psize(2)
            if score_map(i,j-1) > score_map(i-1,j)
                score_map(i,j) = score_map(i,j-1) + decoded_event(i,j);
                path_map{i,j} = [path_map{i,j-1}; [i,j]];
            else
                score_map(i,j) = score_map(i-1,j) + decoded_event(i-1,j);
                path_map{i,j} = [path_map{i-1,j}; [i,j]];
            end
        end
    end
    score = score_map(psize(1),psize(2));
    path = path_map{psize(1),psize(2)};
end


