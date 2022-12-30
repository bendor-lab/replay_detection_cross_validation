
function out=batch_folders(folders, option)
rexposure_exception=[];  %if you have four tracks but they are not rexposures, make this not empty...
rexposure=[];  %if not empty, assumes that there is rexposure to track

%order of running batch analysis
%EXTRACT_SPIKES_MARTA or EXTRACT_SPIKES_MARGOT
%EXTRACT_VIDEO
%PRE
%POST

current_directory=pwd;
out=[];

if isempty(folders)
    folders{1}=NaN;
end
for i=1:length(folders)
    if ~isnan(folders{i})
        cd(folders{i})
    end
    if isempty(option)   %default
        
    elseif strcmp(option,'EXTRACT_SPIKES_MARTA') | strcmp(option,'EXTRACT_SPIKES_MARGOT') 
         batch_analysis(option); 
    elseif strcmp(option,'EXTRACT_VIDEO') 
         batch_analysis(option); 
    elseif strcmp(option,'PRE')
        batch_analysis(option); 
    elseif strcmp(option,'POST')
        batch_analysis('POST');
            batch_analysis('CSC');
                batch_analysis('PLACE FIELDS');
                    batch_analysis('REPLAY');
                        batch_analysis('ANALYSIS');
    elseif strcmp(option,'REPLAY')
        batch_analysis('PLACE FIELDS');
            batch_analysis('REPLAY');
                batch_analysis('ANALYSIS');
      elseif strcmp(option,'REPLAY ONLY')
            batch_analysis('REPLAY');
                
    elseif strcmp(option,'SLEEP')
        batch_analysis(option);
    elseif strcmp(option,'ANALYSIS')
         batch_analysis(option);
    elseif strcmp(option,'plot_cleaning_steps') | strcmp(option,'plot_place_fields') | strcmp(option,'plot_significant_events')
         batch_analysis(option);
    elseif strcmp(option,'CHIMERIC')
        chimeric_track_analysis;
    elseif strcmp(option,'RATE_REMAPPED')
        rate_remapped_track_analysis;
        elseif strcmp(option,'GLOBAL_REMAPPED')
        global_remapped_track_analysis;
    end
    if ~isnan(folders{i})
        cd(current_directory);
    end
end

%ANALYSIS OPTIONS
if strcmp(option,'RATE_REMAP_ONE_TRACK')
    out=rate_remapping_ONE_TRACK(folders);
    plot_rate_remapping('ONE_TRACK');
elseif strcmp(option,'RATE_REMAP_TRACK_PAIRS')
    out=rate_remapping_TRACK_PAIRS(folders);
    plot_rate_remapping('TRACK_PAIRS');
elseif strcmp(option,'CUMULATIVE_REPLAY')
    cumulative_replay(folders);
    plot_cumulative_replay
end



