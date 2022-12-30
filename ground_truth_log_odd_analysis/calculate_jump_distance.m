function jump_distance = calculate_jump_distance(folders,option)

jump_distance = [];
if strcmp(option,'common') || strcmp(option,'original')
    % Load parameters
    for nfolder = 1:length(folders)
        cd(folders{nfolder})
        parameters = list_of_parameters;
        load('decoded_replay_events.mat')
%         load('decoded_replay_events_segments.mat')

        for j = 1:length(decoded_replay_events)
            for i = 1 : length(decoded_replay_events(1).replay_events)
                % find the max probability
                est_pos_tmp= decoded_replay_events(j).replay_events(i).decoded_position;
                est_pos_tmp(est_pos_tmp< 0.02)=0;
                [max_prob,max_prob_idx]= max(est_pos_tmp,[],1);
                max_prob_idx(max_prob==0)= []; % ignore those with no spikes (max still gives an index)
                gaps= abs(diff(max_prob_idx));
                if length(gaps)>3 % ignore if jumps at beginning or end of event
                    gaps(1)=0;
                    gaps(end)=0;
                end
                jump_distance{nfolder}{j}{i} = max(gaps);
            end
        end
        
        cd ..
    end
    
    cd ground_truth_original
    save jump_distance jump_distance
    cd ..
elseif strcmp(option,'global remapped')
    for nfolder = 1:length(folders)
        cd(folders{nfolder})
        cd .\global_remapped_shuffles
        DIR = dir('shuffle_*');
        DataPath = natsortfiles({DIR.name})'
        cd ..

        for f = 1:length(DataPath)
            destination = ['global_remapped_shuffles\shuffle_' num2str(f)]
            cd(destination)

            parameters = list_of_parameters;
            load decoded_replay_events
%             load decoded_replay_events_segments

            

            for j = 1:length(decoded_replay_events)
                for i = 1 : length(decoded_replay_events(1).replay_events)
                    % find the max probability
                    est_pos_tmp= decoded_replay_events(j).replay_events(i).decoded_position;
                    est_pos_tmp(est_pos_tmp< 0.02)=0;
                    [max_prob,max_prob_idx]= max(est_pos_tmp,[],1);
                    max_prob_idx(max_prob==0)= []; % ignore those with no spikes (max still gives an index)
                    gaps= abs(diff(max_prob_idx));
                    if length(gaps)>3 % ignore if jumps at beginning or end of event
                        gaps(1)=0;
                        gaps(end)=0;
                    end
                    jump_distance{f}{nfolder}{j}{i} = max(gaps);
                end
            end
            
            cd ..
            cd ..
        end
        cd ..
    end
    cd ground_truth_original
    save jump_distance_global_remapped jump_distance
    cd ..
end
end