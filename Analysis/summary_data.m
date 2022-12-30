function summary_data(folders,option)
% goes through each folder and gets info of interest, compiles into a table
switch option
    case 'wcorr'
        load('rate_remapping_analysis_TRACK_PAIRS_wcorr.mat');
    case 'spearman'
        load('rate_remapping_analysis_TRACK_PAIRS_spearman.mat');
end
load('.\Tables\decoding_error_comparison.mat');

data_summary= table;
counts=1;
parameters= list_of_parameters;
master_folder = pwd;
 for this_folder=1:length(folders)
            cd([master_folder '\' folders{this_folder}]);
            disp(['       ' folders{this_folder}]);
            
            % load files of interest
            load('extracted_place_fields_BAYESIAN.mat');
            load('extracted_sleep_state.mat');
            load('lap_times.mat');
            load('estimated_position_leave_one_out.mat');
            
            for this_track=1:length(place_fields_BAYESIAN.track)
                % Identificators
                if contains(pwd,'Rate remapping')
                    data_summary.session{counts}= folders(this_folder,1);
                else
                    fd= folders{this_folder};
                    fd= strsplit(fd, '\');
                    data_summary.session{counts}= fd(2);
                    data_summary.rat{counts}= fd(1);
                end
                data_summary.track(counts)= this_track;
                
                % Behaviour and Sleep
                data_summary.number_laps(counts)= floor(lap_times(this_track).total_number_of_laps/2);
                if length(lap_times)==4 && this_track<3 % re exposure session
                    data_summary.sleep_PRE(counts)= sleep_state.time_slept.PRE;
                    data_summary.sleep_POST(counts)= sleep_state.time_slept.INTER_post;
                elseif length(lap_times)==4 && this_track>2 % re exposure session
                    data_summary.sleep_PRE(counts)= sleep_state.time_slept.INTER_post;
                    data_summary.sleep_POST(counts)= sleep_state.time_slept.FINAL_post;
                else % other sessions
                    data_summary.sleep_PRE(counts)= sleep_state.time_slept.PRE;
                    if isfield(sleep_state.time_slept,'FINAL_POST')
                        data_summary.sleep_POST(counts)= sleep_state.time_slept.FINAL_POST;
                    elseif isfield(sleep_state.time_slept,'FINAL_post')
                        data_summary.sleep_POST(counts)= sleep_state.time_slept.FINAL_post;
                    end
                end
                
                load time_range.mat;
                data_summary.post(counts)= diff(time_range.post)./60;
                data_summary.pre(counts)= diff(time_range.pre)./60;
                 
                % Cells
                good_cells_track= place_fields_BAYESIAN.track(this_track).good_cells;
                data_summary.total_num_cells(counts)= numel(good_cells_track);
                data_summary.num_common_cells(counts)= numel(remapping(2).common_good_cells(remapping(2).experiment == this_folder));
                data_summary.included_cells(counts)= numel(intersect(remapping(2).ID_active_cells_during_replay(remapping(2).experiment == this_folder),good_cells_track));
                
                % Decoding
                row_idx= find(contains(decoding_error.session, data_summary.session{counts}) & decoding_error.track== this_track & contains(decoding_error.type,'standard'));
                data_summary.classification_accuracy(counts)= decoding_error.accuracy(row_idx);
                data_summary.local_decoding_accuracy(counts)= decoding_error.percent_small(row_idx);
                
                counts= counts+1;
            end
            
 end
        cd(master_folder)
        writetable(data_summary,'.\Tables\summary_data.csv')   
            
end