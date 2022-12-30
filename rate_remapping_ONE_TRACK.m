function out=rate_remapping_ONE_TRACK(folders)

if isempty(folders)
    a=pwd;
    folders={a};
    cd ..
end
cd(folders{1});

data=compare_replay_across_tracks;
for epoch=1:size(data,1)
    for track=1:size(data,2)
        out(epoch,track).experiment=[];
        out(epoch,track).place_fields_BAYESIAN=[];   
        out(epoch,track).mean_replay_rate=[];
       out(epoch,track).mean_replay_spikes=[];
        out(epoch,track).mean_replay_spikes_nonZero=[];
    end
end
cd ..

for i=1:length(folders)
    cd(folders{i});
    data=compare_replay_across_tracks;
    
    for epoch=1:size(data,1)
        for track=1:size(data,2)
            out(epoch,track).folder=folders{i};
            out(epoch,track).experiment=[out(epoch,track).experiment; (i*ones(size(data(epoch,track).mean_replay_rate)))];
            out(epoch,track).place_fields_BAYESIAN=[out(epoch,track).place_fields_BAYESIAN data(epoch,track).place_fields_BAYESIAN];
            out(epoch,track).mean_replay_rate=[out(epoch,track).mean_replay_rate; data(epoch,track).mean_replay_rate];
            out(epoch,track).mean_replay_spikes=[out(epoch,track).mean_replay_spikes; data(epoch,track).mean_replay_spikes];
            out(epoch,track).mean_replay_spikes_nonZero=[out(epoch,track).mean_replay_spikes_nonZero; data(epoch,track).mean_replay_spikes_nonZero];
        end
    end
    cd ..
end

save rate_remapping_analysis_ONE_TRACK out
end


function out=compare_replay_across_tracks
load extracted_place_fields_BAYESIAN
load significant_replay_events
load sorted_replay

number_of_tracks=length(sorted_replay);
for track=1:number_of_tracks
    for epoch=1:2  %PRE or POST
        if epoch==1
            index1=sorted_replay(track).index.sleepPRE;
            events1=significant_replay_events.track(track);
        elseif epoch==2
            index1=sorted_replay(track).index.sleepPOST;
            events1=significant_replay_events.track(track);
        end
      
        out(epoch,track).mean_replay_rate=[];
        out(epoch,track).mean_replay_spikes=[];
        out(epoch,track).place_fields_BAYESIAN=[];
    
        out(epoch,track).all_cells=place_fields_BAYESIAN.pyramidal_cells;
        
        if length(index1)>0 
            % all pyramidal cells for track 1
            for j=1:length(index1)
                for i=1:length(out(epoch,track).all_cells)
                    out(epoch,track).spikes(i,j)=length(find(events1.spikes{index1(j)}(:,1)==out(epoch,track).all_cells(i)));
                    out(epoch,track).rate(i,j)=length(find(events1.spikes{index1(j)}(:,1)==out(epoch,track).all_cells(i)))/events1.event_duration(index1(j));
                end
            end
                  
            %track 1 only - all pyramidal cells
            
            out(epoch,track).place_fields_BAYESIAN=place_fields_BAYESIAN.track(track).raw_peak(out(epoch,track).all_cells);
            out(epoch,track).mean_replay_rate=mean(out(epoch,track).rate,2);
            out(epoch,track).mean_replay_spikes=mean(out(epoch,track).spikes,2); 
            out(epoch,track).mean_replay_spikes_nonZero=sum(out(epoch,track).spikes,2)./sum(sign(out(epoch,track).spikes),2); %mean of events with 1 or more spikes
        end
    end
end
end