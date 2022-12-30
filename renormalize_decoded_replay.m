function [decoded_replay_events2, bayesian_bias]=renormalize_decoded_replay(decoded_replay_events,tracks_to_compare)

if ~isempty(tracks_to_compare)
    for track_pair = 1:length(tracks_to_compare)
        
        for i=1:length(decoded_replay_events(1).replay_events)
            rescaling_matrix=sum(decoded_replay_events(tracks_to_compare{track_pair}(1)).replay_events(i).decoded_position...
                +decoded_replay_events(tracks_to_compare{track_pair}(2)).replay_events(i).decoded_position,1);
            
            decoded_replay_events2(tracks_to_compare{track_pair}(1)).replay_events(i).decoded_position=decoded_replay_events(tracks_to_compare{track_pair}(1)).replay_events(i).decoded_position./rescaling_matrix;
            decoded_replay_events2(tracks_to_compare{track_pair}(2)).replay_events(i).decoded_position=decoded_replay_events(tracks_to_compare{track_pair}(2)).replay_events(i).decoded_position./rescaling_matrix;
            
            bayesian_bias(tracks_to_compare{track_pair}(1)).replay_events(i)=sum(sum(decoded_replay_events2(tracks_to_compare{track_pair}(1)).replay_events(i).decoded_position));
            bayesian_bias(tracks_to_compare{track_pair}(2)).replay_events(i)=sum(sum(decoded_replay_events2(tracks_to_compare{track_pair}(2)).replay_events(i).decoded_position));
            
            bayesian_sum= bayesian_bias(tracks_to_compare{track_pair}(1)).replay_events(i)+bayesian_bias(tracks_to_compare{track_pair}(2)).replay_events(i);  %add both tracks together, and use this to normalize so that the sum adds to one
            bayesian_bias(tracks_to_compare{track_pair}(1)).replay_events(i)=bayesian_bias(tracks_to_compare{track_pair}(1)).replay_events(i)/bayesian_sum;
            bayesian_bias(tracks_to_compare{track_pair}(2)).replay_events(i)=bayesian_bias(tracks_to_compare{track_pair}(2)).replay_events(i)/bayesian_sum;
            
        end
        
    end
    
else
    for i=1:length(decoded_replay_events(1).replay_events)
        rescaling_matrix=sum(decoded_replay_events(1).replay_events(i).decoded_position+decoded_replay_events(2).replay_events(i).decoded_position,1);
        
        decoded_replay_events2(1).replay_events(i).decoded_position=decoded_replay_events(1).replay_events(i).decoded_position./rescaling_matrix;
        decoded_replay_events2(2).replay_events(i).decoded_position=decoded_replay_events(2).replay_events(i).decoded_position./rescaling_matrix;
        
        bayesian_bias(1).replay_events(i)=sum(sum(decoded_replay_events2(1).replay_events(i).decoded_position));
        bayesian_bias(2).replay_events(i)=sum(sum(decoded_replay_events2(2).replay_events(i).decoded_position));
        
        bayesian_sum= bayesian_bias(1).replay_events(i)+bayesian_bias(2).replay_events(i);  %add both tracks together, and use this to normalize so that the sum adds to one
        bayesian_bias(1).replay_events(i)=bayesian_bias(1).replay_events(i)/bayesian_sum;
        bayesian_bias(2).replay_events(i)=bayesian_bias(2).replay_events(i)/bayesian_sum;
    end
end