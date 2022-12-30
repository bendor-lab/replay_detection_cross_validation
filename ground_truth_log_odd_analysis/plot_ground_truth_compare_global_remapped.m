function [] = plot_ground_truth_compare_global_remapped(folders,method,option)

total_number{1}(1:3) = 0;
total_number{2}(1:3) = 0;

for nfolder = 1:length(folders) % Get total number of replay events
    for nshuffle = 1:2 % 1 is original and 2 is global remapped
        session_total_number{nshuffle}(nfolder,1:3) = 0;
    end
end

for nfolder = 1:length(folders) % Get total number of replay events
    for nshuffle = 1:2 % 1 is original and 2 is global remapped

        cd(folders{nfolder})

        if nshuffle == 2
            cd global_remapped_shuffles
            number_of_global_remapped_shuffles = length(dir('shuffle_*'));
            cd ..
        elseif nshuffle == 1 % Original data, no shuffle
            number_of_global_remapped_shuffles = 1;
        end

        load decoded_replay_events

        [~,time_range]=sort_replay_events([],'wcorr');

        for event = 1:length(decoded_replay_events(1).replay_events)

            event_time = decoded_replay_events(1).replay_events(event).timebins_edges(1);

           if event_time >= time_range.pre(1) & event_time <= time_range.pre(2) %If PRE

                total_number{nshuffle}(1) = total_number{nshuffle}(1) + 1*number_of_global_remapped_shuffles;
                session_total_number{nshuffle}(nfolder,1) = session_total_number{nshuffle}(nfolder,1) + 1*number_of_global_remapped_shuffles;

            elseif event_time >= time_range.post(1) & event_time <= time_range.post(2) %If POST

                total_number{nshuffle}(3) = total_number{nshuffle}(3) + 1*number_of_global_remapped_shuffles;
                session_total_number{nshuffle}(nfolder,3) = session_total_number{nshuffle}(nfolder,3) + 1*number_of_global_remapped_shuffles;

            elseif event_time >= time_range.track(1).behaviour(1) & event_time <= time_range.track(1).behaviour(2)
                total_number{nshuffle}(2) = total_number{nshuffle}(2) + 1*number_of_global_remapped_shuffles;
                session_total_number{nshuffle}(nfolder,2) = session_total_number{nshuffle}(nfolder,2) + 1*number_of_global_remapped_shuffles;

            elseif event_time >= time_range.track(2).behaviour(1) & event_time <= time_range.track(2).behaviour(2)
                total_number{nshuffle}(2) = total_number{nshuffle}(2) + 1*number_of_global_remapped_shuffles;
                session_total_number{nshuffle}(nfolder,2) = session_total_number{nshuffle}(nfolder,2) + 1*number_of_global_remapped_shuffles;

            end
        end
        cd ..
    end
end


index = [];
states = [-1 0 1 2];
epoch_index = [];
session_index = [];

% method = {'wcorr one shuffle','wcorr three shuffles','wcorr one shuffle jump distance','linear one shuffle','spearman','spearman all spikes'};
jump_threshold = 0.4;

for nmethod = 1:length(method)
    log_odd_compare = [];
    cd ground_truth_original
    if strcmp(method{nmethod},'wcorr one shuffle')
        load log_odd_wcorr_place_bin_circular_shift
        log_odd_compare{1} = log_odd
        load log_odd_wcorr_one_shuffle_global_remapped

    elseif strcmp(method{nmethod},'wcorr three shuffles')
        load log_odd_wcorr
        log_odd_compare{1} = log_odd
        load log_odd_wcorr_global_remapped

    elseif strcmp(method{nmethod},'wcorr PRE Place POST Place')
        load log_odd_wcorr_PRE_place_POST_place
        log_odd_compare{1} = log_odd
        load log_odd_wcorr_PRE_place_POST_place_global_remapped

    elseif strcmp(method{nmethod},'wcorr PRE Place POST Time')
        load log_odd_wcorr_PRE_place_POST_time
        log_odd_compare{1} = log_odd
        load log_odd_wcorr_PRE_place_POST_time_global_remapped

    elseif strcmp(method{nmethod},'wcorr one shuffle jump distance')
        load log_odd_wcorr_place_bin_circular_shift
        log_odd_compare{1} = log_odd
        load log_odd_wcorr_one_shuffle_global_remapped

    elseif strcmp(method{nmethod},'linear PRE Place POST Time')
        load log_odd_linear_PRE_place_POST_time
        log_odd_compare{1} = log_odd
        load log_odd_linear_PRE_place_POST_time_global_remapped

    elseif strcmp(method{nmethod},'spearman')
        load log_odd_spearman
        log_odd_compare{1} = log_odd
        load log_odd_spearman_global_remapped

    elseif strcmp(method{nmethod},'spearman all spikes')
        load log_odd_spearman_all_spikes
        log_odd_compare{1} = log_odd
        load log_odd_spearman_all_spikes_global_remapped

    end

    log_odd_compare{2} = log_odd

    cd ..

    % Get index for PRE,  RUN, POST
    % states = [-1 0 1 2 3 4 5]; % sleepPRE, awakePRE, RUN T1, RUN T2, sleepPOST and awakePOST and multi-track
    % states = [-1 0 1 2 3 4]; % sleepPRE, awakePRE, RUN T1, RUN T2, sleepPOST and awakePOST

    states = [-1 0 1 2];

    for nshuffle = 1:2
        for k=1:length(states) % sort track id and behavioral states (pooled from 5 sessions)
            state_index = find(log_odd_compare{nshuffle}.behavioural_state==states(k));

            if k == 1 % PRE
                epoch_index{nmethod}{nshuffle}{1} = state_index;
            elseif k == 2 % POST
                epoch_index{nmethod}{nshuffle}{3} = state_index;
            elseif k == 3 % RUN Track 1
                epoch_index{nmethod}{nshuffle}{2} = state_index;
            elseif k == 4 % RUN Track 2
                epoch_index{nmethod}{nshuffle}{2} = [epoch_index{nmethod}{nshuffle}{2} state_index];
            end
        end
    end


    % index = epoch_index;

    if strcmp(option,'common')
        data{nmethod}{1} = log_odd_compare{1}.common_zscored.original;
        log_pval{nmethod}{1} = log_odd_compare{1}.pvalue;
        segment_id{nmethod}{1} = log_odd_compare{1}.segment_id;
        if strcmp(method{nmethod},'wcorr one shuffle jump distance')
            max_jump{1} = log_odd_compare{1}.max_jump;
        end

        data{nmethod}{2} = log_odd_compare{2}.common_zscored.global_remapped_original;
        log_pval{nmethod}{2} = log_odd_compare{2}.pvalue;
        segment_id{nmethod}{2} = log_odd_compare{2}.segment_id;
        if strcmp(method{nmethod},'wcorr one shuffle jump distance')
            max_jump{2} = log_odd_compare{2}.max_jump;
        end

%             for session = 1:max(log_odd_compare{nshuffle}.experiment)
%                 session_index{nmethod}{nshuffle}{session} = find(log_odd_compare{nshuffle}.experiment == session);
%             end


    elseif strcmp(option,'original')
        data{nmethod}{1} = log_odd_compare{1}.normal_zscored.original;
        log_pval{nmethod}{1} = log_odd_compare{1}.pvalue;
        segment_id{nmethod}{1} = log_odd_compare{1}.segment_id;

        data{nmethod}{2} = log_odd_compare{2}.normal_zscored.global_remapped_original;
        log_pval{nmethod}{2} = log_odd_compare{2}.pvalue;
        segment_id{nmethod}{2} = log_odd_compare{2}.segment_id;
        %             for session = 1:max(log_odd_compare{nshuffle}.experiment)
        %                 session_index{nmethod}{nshuffle}{session} = log_odd_compare{nshuffle}.index(find(log_odd_compare{nshuffle}.experiment == session));
        %             end

    end

end

% log_odd_compare = [];

p_val_threshold = round(10.^[-3:0.01:log10(0.2)],4);
% p_val_threshold = 0.001:0.0001:0.2;
p_val_threshold = flip(unique(p_val_threshold));
p_val_threshold(1) = 0.2;

log_odd_difference = [];
log_odd_difference_CI = [];
percent_sig_events = [];

cd ground_truth_original
if exist('log_odd_difference_cell_id_shuffles.mat', 'file') ~= 2
    load log_odd_difference_compare_spearman
    % copy spearman all spikes
    log_odd_difference1{1}{1} = log_odd_difference{2};
    log_odd_difference_CI1{1}{1} = log_odd_difference_CI{2};
    percent_sig_events1{1}{1} = percent_sig_events{2};
    percent_sig_events_CI1{1}{1} = percent_sig_events_CI{2};
    percent_multi_events1{1}{1} = percent_multi_events{2};
    percent_multi_events_CI1{1}{1} = percent_multi_events_CI{2};

    % copy spearman median spike
    log_odd_difference1{2}{1} = log_odd_difference{1};
    log_odd_difference_CI1{2}{1} = log_odd_difference_CI{1};
    percent_sig_events1{2}{1} = percent_sig_events{1};
    percent_sig_events_CI1{2}{1} = percent_sig_events_CI{1};
    percent_multi_events1{2}{1} = percent_multi_events{1};
    percent_multi_events_CI1{2}{1} = percent_multi_events_CI{1};
    
    load log_odd_difference_jump_distance
    % copy wcorr 1 shuffle jump distance 0.4 data
    for epoch = 1:3
        log_odd_difference1{3}{1}{epoch} = log_odd_difference{1}{3}{epoch}(40:end,:);
        log_odd_difference_CI1{3}{1}{epoch} = log_odd_difference_CI{1}{3}{epoch}(40:end,:);
        percent_sig_events1{3}{1}{epoch} = percent_sig_events{1}{3}{epoch}(40:end,:);
        percent_sig_events_CI1{3}{1}{epoch} = percent_sig_events_CI{1}{3}{epoch}(40:end,:);
        percent_multi_events1{3}{1}{epoch} = percent_multi_events{1}{3}{epoch}(40:end,:);
        percent_multi_events_CI1{3}{1}{epoch} = percent_multi_events_CI{1}{3}{epoch}(40:end,:);
    end

    load log_odd_difference_multiple_shuffles
    % copy wcorr PRE place POST time data
    log_odd_difference1{4}{1} = log_odd_difference{1}{3};
    log_odd_difference_CI1{4}{1} = log_odd_difference_CI{1}{3};
    percent_sig_events1{4}{1} = percent_sig_events{1}{3};
    percent_sig_events_CI1{4}{1} = percent_sig_events_CI{1}{3};
    percent_multi_events1{4}{1} = percent_multi_events{1}{3};
    percent_multi_events_CI1{4}{1} = percent_multi_events_CI{1}{3};


    % copy wcorr PRE spike PRE place POST place data
    log_odd_difference1{6}{1} = log_odd_difference{1}{6};
    log_odd_difference_CI1{6}{1} = log_odd_difference_CI{1}{6};
    percent_sig_events1{6}{1} = percent_sig_events{1}{6};
    percent_sig_events_CI1{6}{1} = percent_sig_events_CI{1}{6};
    percent_multi_events1{6}{1} = percent_multi_events{1}{6};
    percent_multi_events_CI1{6}{1} = percent_multi_events_CI{1}{6};

    load log_odd_difference_compare_methods
    % copy linear PRE place POST time data
    log_odd_difference1{5}{1} = log_odd_difference{3};
    log_odd_difference_CI1{5}{1} = log_odd_difference_CI{3};
    percent_sig_events1{5}{1} = percent_sig_events{3};
    percent_sig_events_CI1{5}{1} = percent_sig_events_CI{3};
    percent_multi_events1{5}{1} = percent_multi_events{3};
    percent_multi_events_CI1{5}{1} = percent_multi_events_CI{3};


    for nmethod = 1:6
        for epoch = 1:3
            for nshuffle = 2:length(log_pval{nmethod})
                %                 multi_event_index = find(log_odd_compare{nshuffle}.multi_track == 1);
                tic
                parfor threshold = 1:length(p_val_threshold)
                    for nboot = 1:1000 % Bootstrapping 1000 times
                        %     resampled_event = datasample([epoch_index{nmethod}{nshuffle}{epoch}{1} epoch_index{nmethod}{nshuffle}{epoch}{2}],...
                        %         length([epoch_index{nmethod}{nshuffle}{epoch}{1} epoch_index{nmethod}{nshuffle}{epoch}{2}]));
                        %     resampled_event = [datasample(epoch_index{nmethod}{nshuffle}{epoch}{1},...
                        %         length(epoch_index{nmethod}{nshuffle}{epoch}{1})) datasample(epoch_index{nmethod}{nshuffle}{epoch}{2},...
                        %         length(epoch_index{nmethod}{nshuffle}{epoch}{2}))];
                        %                         resampled_event = epoch_index{nmethod}{nshuffle}{epoch};
                        resampled_event = datasample(epoch_index{nmethod}{nshuffle}{epoch},length(epoch_index{nmethod}{nshuffle}{epoch}));

                        %     for track = 1:2
                        this_segment_event = resampled_event(find(segment_id{nmethod}{nshuffle}(1,resampled_event) == 1));
                        track_1_index = this_segment_event(find(log_pval{nmethod}{nshuffle}(1,this_segment_event) <= p_val_threshold(threshold)));
                        this_segment_event = resampled_event(find(segment_id{nmethod}{nshuffle}(1,resampled_event) > 1));
                        track_1_index = [track_1_index this_segment_event(find(log_pval{nmethod}{nshuffle}(1,this_segment_event) <= p_val_threshold(threshold)/2))];
                        if nmethod == 3
                            track_1_index = track_1_index(find(max_jump{nshuffle}(1,track_1_index) <= 20*0.4));
                        end

                        this_segment_event = resampled_event(find(segment_id{nmethod}{nshuffle}(2,resampled_event) == 1));
                        track_2_index = this_segment_event(find(log_pval{nmethod}{nshuffle}(2,this_segment_event) <= p_val_threshold(threshold)));
                        this_segment_event = resampled_event(find(segment_id{nmethod}{nshuffle}(2,resampled_event) > 1));
                        track_2_index = [track_2_index this_segment_event(find(log_pval{nmethod}{nshuffle}(2,this_segment_event) <= p_val_threshold(threshold)/2))];
                        if nmethod == 3
                            track_2_index = track_2_index(find(max_jump{nshuffle}(2,track_2_index) <= 20*0.4));
                        end

                        %     end
                        multi_event_number = sum(ismember(track_1_index,track_2_index));
                        multi_event_percent = multi_event_number/(length(track_1_index)+length(track_2_index)-multi_event_number);

                        boot_log_odd_difference(threshold,nboot) = mean(data{nmethod}{nshuffle}(1,track_1_index)) ...
                            - mean(data{nmethod}{nshuffle}(2,track_2_index));

                        boot_percent_multi_events(threshold,nboot) = multi_event_percent;
                        % Significant event proportion (minus multitrack event number to avoid double counting)
                        boot_percent_sig_events(threshold,nboot) = (length(track_1_index)+length(track_2_index)-multi_event_number)/total_number{nshuffle}(epoch);
                    end
                end

                log_odd_difference1{nmethod}{nshuffle}{epoch} = boot_log_odd_difference;
                log_odd_difference_CI1{nmethod}{nshuffle}{epoch} = prctile(boot_log_odd_difference,[2.5 97.5],2);
                percent_sig_events1{nmethod}{nshuffle}{epoch} = boot_percent_sig_events;
                percent_sig_events_CI1{nmethod}{nshuffle}{epoch} = prctile(boot_percent_sig_events,[2.5 97.5],2);
                percent_multi_events1{nmethod}{nshuffle}{epoch} = boot_percent_multi_events;
                percent_multi_events_CI1{nmethod}{nshuffle}{epoch} = prctile(boot_percent_multi_events,[2.5 97.5],2);
                toc
            end

        end
    end

%     nmethod = 3;
%     for epoch = 1:3
%         for nshuffle = 2:2
%             %                 multi_event_index = find(log_odd_compare{nshuffle}.multi_track == 1);
%             tic
%             parfor threshold = 1:length(p_val_threshold)
%                 for nboot = 1:1000 % Bootstrapping 1000 times
%                     %     resampled_event = datasample([epoch_index{nmethod}{nshuffle}{epoch}{1} epoch_index{nmethod}{nshuffle}{epoch}{2}],...
%                     %         length([epoch_index{nmethod}{nshuffle}{epoch}{1} epoch_index{nmethod}{nshuffle}{epoch}{2}]));
%                     %     resampled_event = [datasample(epoch_index{nmethod}{nshuffle}{epoch}{1},...
%                     %         length(epoch_index{nmethod}{nshuffle}{epoch}{1})) datasample(epoch_index{nmethod}{nshuffle}{epoch}{2},...
%                     %         length(epoch_index{nmethod}{nshuffle}{epoch}{2}))];
%                     %                         resampled_event = epoch_index{nmethod}{nshuffle}{epoch};
%                     resampled_event = datasample(epoch_index{nmethod}{nshuffle}{epoch},length(epoch_index{nmethod}{nshuffle}{epoch}));
% 
%                     %     for track = 1:2
%                     this_segment_event = resampled_event(find(segment_id{nmethod}{nshuffle}(1,resampled_event) == 1));
%                     track_1_index = this_segment_event(find(log_pval{nmethod}{nshuffle}(1,this_segment_event) <= p_val_threshold(threshold)));
%                     this_segment_event = resampled_event(find(segment_id{nmethod}{nshuffle}(1,resampled_event) > 1));
%                     track_1_index = [track_1_index this_segment_event(find(log_pval{nmethod}{nshuffle}(1,this_segment_event) <= p_val_threshold(threshold)/2))];
%                     track_1_index = track_1_index(find(max_jump{nshuffle}(1,track_1_index) <= 20*0.4));
% 
%                     this_segment_event = resampled_event(find(segment_id{nmethod}{nshuffle}(2,resampled_event) == 1));
%                     track_2_index = this_segment_event(find(log_pval{nmethod}{nshuffle}(2,this_segment_event) <= p_val_threshold(threshold)));
%                     this_segment_event = resampled_event(find(segment_id{nmethod}{nshuffle}(2,resampled_event) > 1));
%                     track_2_index = [track_2_index this_segment_event(find(log_pval{nmethod}{nshuffle}(2,this_segment_event) <= p_val_threshold(threshold)/2))];
%                     track_2_index = track_2_index(find(max_jump{nshuffle}(2,track_2_index) <= 20*0.4));
% 
%                     %     end
%                     multi_event_number = sum(ismember(track_1_index,track_2_index));
%                     multi_event_percent = multi_event_number/(length(track_1_index)+length(track_2_index)-multi_event_number);
% 
%                     boot_log_odd_difference(threshold,nboot) = mean(data{nmethod}{nshuffle}(1,track_1_index)) ...
%                         - mean(data{nmethod}{nshuffle}(2,track_2_index));
% 
%                     boot_percent_multi_events(threshold,nboot) = multi_event_percent;
%                     % Significant event proportion (minus multitrack event number to avoid double counting)
%                     boot_percent_sig_events(threshold,nboot) = (length(track_1_index)+length(track_2_index)-multi_event_number)/total_number{nshuffle}(epoch);
%                 end
%             end
% 
%             log_odd_difference1{nmethod}{nshuffle}{epoch} = boot_log_odd_difference;
%             log_odd_difference_CI1{nmethod}{nshuffle}{epoch} = prctile(boot_log_odd_difference,[2.5 97.5],2);
%             percent_sig_events1{nmethod}{nshuffle}{epoch} = boot_percent_sig_events;
%             percent_sig_events_CI1{nmethod}{nshuffle}{epoch} = prctile(boot_percent_sig_events,[2.5 97.5],2);
%             percent_multi_events1{nmethod}{nshuffle}{epoch} = boot_percent_multi_events;
%             percent_multi_events_CI1{nmethod}{nshuffle}{epoch} = prctile(boot_percent_multi_events,[2.5 97.5],2);
%             toc
%         end
% 
%     end

    log_odd_difference = log_odd_difference1;
    log_odd_difference_CI = log_odd_difference_CI1;
    percent_sig_events = percent_sig_events1;
    percent_sig_events_CI = percent_sig_events_CI1;
    percent_multi_events = percent_multi_events1;
    percent_multi_events_CI = percent_multi_events_CI1;

    save log_odd_difference_cell_id_shuffles log_odd_difference log_odd_difference_CI percent_sig_events percent_multi_events percent_sig_events_CI percent_multi_events_CI
    clear boot_log_odd_difference boot_percent_sig_events boot_percent_multi_events
    clear log_odd_difference1 log_odd_difference_CI1 percent_sig_events1 percent_sig_events_CI1 percent_multi_events1 percent_multi_events_CI1
    cd ..

else
    load log_odd_difference_cell_id_shuffles
    cd ..
end


cd ground_truth_original
if exist('log_odd_difference_cell_id_shuffles_original.mat', 'file') ~= 2
    for epoch = 1:3
        for nmethod = 1:6
            for nshuffle = 1:2
                %                 multi_event_index = find(log_odd_compare{nshuffle}.multi_track == 1);
                tic
                for threshold = 1:length(p_val_threshold)

                    %     resampled_event = datasample([epoch_index{nmethod}{nshuffle}{epoch}{1} epoch_index{nmethod}{nshuffle}{epoch}{2}],...
                    %         length([epoch_index{nmethod}{nshuffle}{epoch}{1} epoch_index{nmethod}{nshuffle}{epoch}{2}]));
                    %     resampled_event = [datasample(epoch_index{nmethod}{nshuffle}{epoch}{1},...
                    %         length(epoch_index{nmethod}{nshuffle}{epoch}{1})) datasample(epoch_index{nmethod}{nshuffle}{epoch}{2},...
                    %         length(epoch_index{nmethod}{nshuffle}{epoch}{2}))];
                    resampled_event = epoch_index{nmethod}{nshuffle}{epoch};
                    %                     resampled_event = datasample(epoch_index{nmethod}{nshuffle}{epoch},length(epoch_index{nmethod}{nshuffle}{epoch}));

                    %     for track = 1:2
                    this_segment_event = resampled_event(find(segment_id{nmethod}{nshuffle}(1,resampled_event) == 1));
                    track_1_index = this_segment_event(find(log_pval{nmethod}{nshuffle}(1,this_segment_event) <= p_val_threshold(threshold)));
                    this_segment_event = resampled_event(find(segment_id{nmethod}{nshuffle}(1,resampled_event) > 1));
                    track_1_index = [track_1_index this_segment_event(find(log_pval{nmethod}{nshuffle}(1,this_segment_event) <= p_val_threshold(threshold)/2))];
                    if nmethod == 3
                        track_1_index = track_1_index(find(max_jump{nshuffle}(1,track_1_index) <= 20*0.4));
                    end

                    this_segment_event = resampled_event(find(segment_id{nmethod}{nshuffle}(2,resampled_event) == 1));
                    track_2_index = this_segment_event(find(log_pval{nmethod}{nshuffle}(2,this_segment_event) <= p_val_threshold(threshold)));
                    this_segment_event = resampled_event(find(segment_id{nmethod}{nshuffle}(2,resampled_event) > 1));
                    track_2_index = [track_2_index this_segment_event(find(log_pval{nmethod}{nshuffle}(2,this_segment_event) <= p_val_threshold(threshold)/2))];
                    if nmethod == 3
                        track_2_index = track_2_index(find(max_jump{nshuffle}(2,track_2_index) <= 20*0.4));
                    end
                    %     end
                    multi_event_number = sum(ismember(track_1_index,track_2_index));
                    multi_event_percent = multi_event_number/(length(track_1_index)+length(track_2_index)-multi_event_number);

                    boot_log_odd_difference(threshold) = mean(data{nmethod}{nshuffle}(1,track_1_index)) ...
                        - mean(data{nmethod}{nshuffle}(2,track_2_index));

                    boot_percent_multi_events(threshold) = multi_event_percent;
                    % Significant event proportion (minus multitrack event number to avoid double counting)
                    boot_percent_sig_events(threshold) = (length(track_1_index)+length(track_2_index)-multi_event_number)/total_number{nshuffle}(epoch);

                end

                log_odd_difference_original{nmethod}{nshuffle}{epoch} = boot_log_odd_difference;
                percent_sig_events_original{nmethod}{nshuffle}{epoch} = boot_percent_sig_events;
                percent_multi_events_original{nmethod}{nshuffle}{epoch} = boot_percent_multi_events;
                toc

            end

        end
    end
% 
%     for epoch = 1:3
%         for nmethod = 3:3
%             for nshuffle = 1:2
%                 %                 multi_event_index = find(log_odd_compare{nshuffle}.multi_track == 1);
%                 tic
%                 for threshold = 1:length(p_val_threshold)
% 
%                     %     resampled_event = datasample([epoch_index{nmethod}{nshuffle}{epoch}{1} epoch_index{nmethod}{nshuffle}{epoch}{2}],...
%                     %         length([epoch_index{nmethod}{nshuffle}{epoch}{1} epoch_index{nmethod}{nshuffle}{epoch}{2}]));
%                     %     resampled_event = [datasample(epoch_index{nmethod}{nshuffle}{epoch}{1},...
%                     %         length(epoch_index{nmethod}{nshuffle}{epoch}{1})) datasample(epoch_index{nmethod}{nshuffle}{epoch}{2},...
%                     %         length(epoch_index{nmethod}{nshuffle}{epoch}{2}))];
%                     resampled_event = epoch_index{nmethod}{nshuffle}{epoch};
%                     %                     resampled_event = datasample(epoch_index{nmethod}{nshuffle}{epoch},length(epoch_index{nmethod}{nshuffle}{epoch}));
% 
%                     %     for track = 1:2
%                     track_1_index = this_segment_event(find(log_pval{nmethod}{nshuffle}(1,this_segment_event) <= p_val_threshold(threshold)));
%                     this_segment_event = resampled_event(find(segment_id{nmethod}{nshuffle}(1,resampled_event) > 1));
%                     track_1_index = [track_1_index this_segment_event(find(log_pval{nmethod}{nshuffle}(1,this_segment_event) <= p_val_threshold(threshold)/2))];
%                     track_1_index = track_1_index(find(max_jump{nshuffle}(1,track_1_index) <= 20*0.4));
% 
%                     this_segment_event = resampled_event(find(segment_id{nmethod}{nshuffle}(2,resampled_event) == 1));
%                     track_2_index = this_segment_event(find(log_pval{nmethod}{nshuffle}(2,this_segment_event) <= p_val_threshold(threshold)));
%                     this_segment_event = resampled_event(find(segment_id{nmethod}{nshuffle}(2,resampled_event) > 1));
%                     track_2_index = [track_2_index this_segment_event(find(log_pval{nmethod}{nshuffle}(2,this_segment_event) <= p_val_threshold(threshold)/2))];
%                     track_1_index = track_1_index(find(max_jump{nshuffle}(2,track_2_index) <= 20*0.4));
% 
%                     %     end
%                     multi_event_number = sum(ismember(track_1_index,track_2_index));
%                     multi_event_percent = multi_event_number/(length(track_1_index)+length(track_2_index)-multi_event_number);
% 
%                     boot_log_odd_difference(threshold) = mean(data{nmethod}{nshuffle}(1,track_1_index)) ...
%                         - mean(data{nmethod}{nshuffle}(2,track_2_index));
% 
%                     boot_percent_multi_events(threshold) = multi_event_percent;
%                     % Significant event proportion (minus multitrack event number to avoid double counting)
%                     boot_percent_sig_events(threshold) = (length(track_1_index)+length(track_2_index)-multi_event_number)/total_number{nshuffle}(epoch);
% 
%                 end
% 
%                 log_odd_difference_original{nmethod}{nshuffle}{epoch} = boot_log_odd_difference;
%                 percent_sig_events_original{nmethod}{nshuffle}{epoch} = boot_percent_sig_events;
%                 percent_multi_events_original{nmethod}{nshuffle}{epoch} = boot_percent_multi_events;
%                 toc
% 
%             end
% 
%         end
%     end

    save log_odd_difference_cell_id_shuffles_original log_odd_difference_original percent_sig_events_original percent_multi_events_original
    cd ..
else
    load log_odd_difference_cell_id_shuffles_original
    cd ..
end

% colour_line= {[254,178,76]/255,[198,219,239]/255,[227,26,28]/255,...
%     [158,202,225]/255,[107,174,214]/255,[49,130,189]/255,[8,81,156]/255};
colour_line= {[254,178,76]/255,[208,209,230]/255,[227,26,28]/255,...
    [166,189,219]/255,[116,169,207]/255,[43,140,190]/255,[4,90,141]/255};



% 


colour_line = {[244,165,130]/255,[127,205,187]/255,[29,145,192]/255,...
    [215,48,39]/255,[37,52,148]/255};
    
% colour_line= {[254,178,76]/255,[208,209,230]/255,[227,26,28]/255,...
%     [166,189,219]/255,[116,169,207]/255,[43,140,190]/255,[4,90,141]/255};

% for n = 1:length(colour_line)
%     plot([0 1],[0 n],'Color',colour_line{n})
%     hold on
% end


p_val_threshold = round(10.^[-3:0.01:log10(0.2)],4);
p_val_threshold = flip(unique(p_val_threshold));
p_val_threshold(1) = 0.2;
low_threshold = find(ismembertol(p_val_threshold,[0.0501 0.02 0.01 0.005 0.002 0.001],eps) == 1);
% colour_line1 ={'#ffffd4','#fee391','#fec44f','#fe9929','#d95f0e','#993404'};
Behavioural_epoches = {'PRE','RUN','POST'};
marker_shape = {'o','^','s'};
fig = figure(2)
count = 1;
fig.Position = [834 116 850 885];
for nmethod = [2 3 4 5 6]
    for epoch = 1:3
        subplot(2,2,2)
        m(count) = scatter(mean(percent_sig_events{nmethod}{2}{epoch}(low_threshold(1),:)),mean(percent_sig_events{nmethod}{1}{epoch}(low_threshold(1),:)),10,...
            marker_shape{epoch},'MarkerFaceColor',colour_line{count},'MarkerEdgeColor',colour_line{count},'MarkerFaceAlpha','0.75')
        hold on
        error_x = [percent_sig_events_CI{nmethod}{2}{epoch}(low_threshold(1),1) percent_sig_events_CI{nmethod}{2}{epoch}(low_threshold(1),2)...
            percent_sig_events_CI{nmethod}{2}{epoch}(low_threshold(1),2) percent_sig_events_CI{nmethod}{2}{epoch}(low_threshold(1),1)];
        error_y = [percent_sig_events_CI{nmethod}{1}{epoch}(low_threshold(1),1) percent_sig_events_CI{nmethod}{1}{epoch}(low_threshold(1),1)...
            percent_sig_events_CI{nmethod}{1}{epoch}(low_threshold(1),2) percent_sig_events_CI{nmethod}{1}{epoch}(low_threshold(1),2)];

        patch(error_x,error_y,colour_line{count},'EdgeColor',colour_line{count},'FaceAlpha','0.2','EdgeAlpha','0.2')
        ylabel('Proportion of significant events in original data')
        xlabel('Proportion of significant events in cell id shuffled data')
        title('Significant events')
    end
    plot([0 0.5],[0 0.5],'r')


    for epoch = 1:3
        subplot(2,2,1)
        m(count) = scatter(mean(percent_multi_events{nmethod}{2}{epoch}(low_threshold(1),:)),mean(percent_multi_events{nmethod}{1}{epoch}(low_threshold(1),:)),10,...
            marker_shape{epoch},'MarkerFaceColor',colour_line{count},'MarkerEdgeColor',colour_line{count},'MarkerFaceAlpha','0.75')
        hold on 
        error_x = [percent_multi_events_CI{nmethod}{2}{epoch}(low_threshold(1),1) percent_multi_events_CI{nmethod}{2}{epoch}(low_threshold(1),2)...
            percent_multi_events_CI{nmethod}{2}{epoch}(low_threshold(1),2) percent_multi_events_CI{nmethod}{2}{epoch}(low_threshold(1),1)];
        error_y = [percent_multi_events_CI{nmethod}{1}{epoch}(low_threshold(1),1) percent_multi_events_CI{nmethod}{1}{epoch}(low_threshold(1),1)...
            percent_multi_events_CI{nmethod}{1}{epoch}(low_threshold(1),2) percent_multi_events_CI{nmethod}{1}{epoch}(low_threshold(1),2)];


        patch(error_x,error_y,colour_line{count},'EdgeColor',colour_line{count},'FaceAlpha','0.2','EdgeAlpha','0.2')
        ylabel('Proportion of multitrack events in original data')
        hold on
        xlabel('Proportion of multitrack events in cell id shuffled data')
        title('Multitrack events')
    end
    plot([0 0.2],[0 0.2],'r')



    for epoch = 1:3
        subplot(2,2,3)
        m(count) = scatter(mean(log_odd_difference{nmethod}{2}{epoch}(low_threshold(1),:)),mean(log_odd_difference{nmethod}{1}{epoch}(low_threshold(1),:)),10,...
            marker_shape{epoch},'MarkerFaceColor',colour_line{count},'MarkerEdgeColor',colour_line{count},'MarkerFaceAlpha','0.75')
        hold on
        error_x = [log_odd_difference_CI{nmethod}{2}{epoch}(low_threshold(1),1) log_odd_difference_CI{nmethod}{2}{epoch}(low_threshold(1),2)...
            log_odd_difference_CI{nmethod}{2}{epoch}(low_threshold(1),2) log_odd_difference_CI{nmethod}{2}{epoch}(low_threshold(1),1)];
        error_y = [log_odd_difference_CI{nmethod}{1}{epoch}(low_threshold(1),1) log_odd_difference_CI{nmethod}{1}{epoch}(low_threshold(1),1)...
            log_odd_difference_CI{nmethod}{1}{epoch}(low_threshold(1),2) log_odd_difference_CI{nmethod}{1}{epoch}(low_threshold(1),2)];

        patch(error_x,error_y,colour_line{count},'EdgeColor',colour_line{count},'FaceAlpha','0.2','EdgeAlpha','0.2')
        ylabel('Mean log odd difference in original data')
        hold on
        xlabel('Mean log odd difference in cell id shuffled data')
        title('Mean log odd difference')
    end
    plot([-0.5 3.5],[-0.5 3.5],'r')
    ylim([-0.5 3.5])
    xlim([-0.5 3.5])


    for epoch = 1:3
        [c index(1)] = min(abs(mean(percent_multi_events{nmethod}{1}{epoch},2) - 0.05));
        [c index(2)] = min(abs(mean(percent_multi_events{nmethod}{2}{epoch},2) - 0.05));
        x = log_odd_difference{nmethod}{1}{epoch}(index(1),:) - log_odd_difference{nmethod}{2}{epoch}(index(2),:);
        y = percent_sig_events{nmethod}{1}{epoch}(index(1),:)-percent_sig_events{nmethod}{2}{epoch}(index(2),:);

        subplot(2,2,4)
        m(count) = scatter(mean(x),mean(y),10,marker_shape{epoch},'MarkerFaceColor',colour_line{count},'MarkerEdgeColor',colour_line{count},'MarkerFaceAlpha','0.75');
        
        hold on
        x_CI = prctile(x,[2.5 97.5],2);
        y_CI = prctile(y,[2.5 97.5],2);
        error_x = [x_CI(1) x_CI(2)...
            x_CI(2) x_CI(1)];
        error_y = [y_CI(1) y_CI(1)...
            y_CI(2) y_CI(2)];

        patch(error_x,error_y,colour_line{count},'EdgeColor',colour_line{count},'FaceAlpha','0.2','EdgeAlpha','0.2')
        ylabel('Shuffle-subtracted proportion of significant events')
        hold on
        xlabel('Shuffle-subtracted log odd difference')
        title('Mean log odd difference vs Significant event difference')
    end

    count = count + 1;
end

legend(m(1:5),{method{2},method{3},method{4},method{5},method{6}},'Position',[0.35 0.2 0.05 0.05])
sgtitle('Replay Detection method Comparision Summary')
cd ground_truth_original\Figure
filename = 'Replay Detection Mehtod Comparision Summary.pdf'
saveas(gcf,filename)
clf
cd ..
cd ..
clf




colour_line = {'',[244,165,130]/255,[127,205,187]/255,[29,145,192]/255,...
[215,48,39]/255,[37,52,148]/255}
p_value_to_plot = round(10.^[-3:0.1:log10(0.2)],4);
% p_val_threshold = 0.001:0.0001:0.2;
p_value_to_plot = flip(unique(p_value_to_plot));
p_value_to_plot(1) = 0.2;
alpha_level = linspace(0.1,0.70,length(p_value_to_plot));
alpha_level1 = linspace(0.1,0.3,length(p_value_to_plot));
low_threshold = find(ismembertol(p_val_threshold,p_value_to_plot,eps) == 1);


for nmethod = [2 3 4 5 6]
    for epoch = 1:3
        fig = figure(1)
        fig.Position = [834 116 850 700];
        subplot(2,2,epoch)
        %
        x = log_odd_difference_original{nmethod}{1}{epoch}';
        y = percent_sig_events_original{nmethod}{1}{epoch}';

        for threshold = 1:length(p_value_to_plot)
            if p_value_to_plot(threshold) > 0.0501
                sc(1) = scatter(x(low_threshold(threshold)),y(low_threshold(threshold)),'filled','MarkerFaceColor',colour_line{nmethod},...
                    'MarkerFaceAlpha',alpha_level(threshold))
                hold on
            else
                sc(1) = scatter(x(low_threshold(threshold)),y(low_threshold(threshold)),'filled','MarkerFaceColor',colour_line{nmethod},...
                    'MarkerFaceAlpha',alpha_level(threshold),'MarkerEdgeColor',colour_line{nmethod})
                hold on
            end
        end


        x = log_odd_difference_original{nmethod}{2}{epoch}';
        y = percent_sig_events_original{nmethod}{2}{epoch}';

        for threshold = 1:length(p_value_to_plot)
            if p_value_to_plot(threshold) > 0.0501
                sc(2) = scatter(x(low_threshold(threshold)),y(low_threshold(threshold)),'filled','MarkerFaceColor','k',...
                    'MarkerFaceAlpha',alpha_level1(threshold))
                hold on
            else
                sc(2) = scatter(x(low_threshold(threshold)),y(low_threshold(threshold)),'filled','MarkerFaceColor','k',...
                    'MarkerFaceAlpha',alpha_level1(threshold),'MarkerEdgeColor',colour_line{nmethod})
                hold on
            end
        end


        xlabel('mean log odd difference')
        ylabel('Proportion of all replay events')
        ylim([0 0.7])
        title(Behavioural_epoches{epoch});
    end


    legend([sc(1) sc(2)], {'Original','Cell ID shuffled'},'Position',[0.7 0.3 0.05 0.05])
    sgtitle(method{nmethod})
    cd ground_truth_original\Figure
    filename = sprintf('%s  Orignal vs shuffle original.pdf',method{nmethod})
    saveas(gcf,filename)
    cd ..
    cd ..
    clf
end




fig = figure(2)
fig.Position = [834 116 850 700];
p_val_threshold = round(10.^[-3:0.01:log10(0.2)],4);
p_val_threshold = flip(unique(p_val_threshold));
p_val_threshold(1) = 0.2;
low_threshold = find(ismembertol(p_val_threshold,[0.0501 0.02 0.01 0.005 0.002 0.001],eps) == 1);
alpha_level = linspace(0.2,0.7,length(low_threshold));
alpha_level1 = linspace(0.1,0.3,length(low_threshold));

colour_line = {'',[244,165,130]/255,[127,205,187]/255,[29,145,192]/255,...
    [215,48,39]/255,[37,52,148]/255}
% colour_line = {[69,117,180]/255,[244,165,130]/255,[215,48,39]/255};
%     colour_line = {[37,52,148]/255,[34,94,168]/255,[29,145,192]/255,[65,182,196]/255};
%     colour_line = flip({[12,44,132]/255,'',[29,145,192]/255,'','',[127,205,187]/255});

for nmethod = [2 3 4 5 6]
    for epoch = 1:3
        subplot(2,2,epoch)

        x = mean(log_odd_difference{nmethod}{1}{epoch},2)';
        y = mean(percent_sig_events{nmethod}{1}{epoch},2)';

        for threshold = 1:length(low_threshold)
            scatter(x(low_threshold(threshold)),y(low_threshold(threshold)),'filled','MarkerFaceColor',colour_line{nmethod},...
                'MarkerFaceAlpha',alpha_level(threshold),'MarkerEdgeColor',colour_line{nmethod})
            hold on
        end

        UCI = log_odd_difference_CI{nmethod}{1}{epoch}(:,2)';
        UCI(isnan(x)) = [];
        LCI = log_odd_difference_CI{nmethod}{1}{epoch}(:,1)';
        LCI(isnan(x)) = [];

        y(isnan(x)) = [];
        x(isnan(x)) = [];
        p(1) = patch([LCI fliplr(UCI)], [y fliplr(y)], colour_line{nmethod},'FaceAlpha','0.3','LineStyle','none');
        hold on



        x = mean(log_odd_difference{nmethod}{2}{epoch},2)';
        y = mean(percent_sig_events{nmethod}{2}{epoch},2)';

        for threshold = 1:length(low_threshold)
            scatter(x(low_threshold(threshold)),y(low_threshold(threshold)),'filled','MarkerFaceColor','k',...
                'MarkerFaceAlpha',alpha_level(threshold),'MarkerEdgeColor',colour_line{nmethod})
            hold on
        end

        UCI = log_odd_difference_CI{nmethod}{2}{epoch}(:,2)';
        UCI(isnan(x)) = [];
        LCI = log_odd_difference_CI{nmethod}{2}{epoch}(:,1)';
        LCI(isnan(x)) = [];

        y(isnan(x)) = [];
        x(isnan(x)) = [];
        p(2) = patch([LCI fliplr(UCI)], [y fliplr(y)], 'k','FaceAlpha','0.3','LineStyle','none');
        hold on
        %             plot(x,y,colour_line{nshuffle});

        xlabel('mean log odd difference')
        ylabel('Proportion of all replay events')
        title(Behavioural_epoches{epoch});
        ylim([0 0.7])
    end


    %     legend([p(1) p(3) p(6)], {'Shuffle 3','Shuffle 2+3','Shuffle 1+2+3'})
    % legend([p(1) p(2) p(3)], {method{1},method{2},method{3}},'Position',[0.7 0.3 0.05 0.05])
    legend([p(1) p(2)], {'Original','Cell ID shuffled'},'Position',[0.7 0.3 0.05 0.05])
    %     legend([p(1) p(3) p(6)], {'PRE Place','PRE Place + POST Place','PRE Place + PRE Spike + POST Place'})
    %     legend([p(1),p(2),p(3),p(4),p(5),p(6)], {'Shuffle 3','Shuffle 3+4','Shuffle 2+3','Shuffle 1+3','Shuffle 1+2','Shuffle 1+2+3'})
    sgtitle(method{nmethod})
    cd ground_truth_original\Figure
    filename = sprintf('%s  Orignal vs shuffle CI.pdf',method{nmethod})
    saveas(gcf,filename)
    cd ..
    cd ..
    clf
end
%



end

