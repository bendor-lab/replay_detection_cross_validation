function example_neuron_raw_data(varargin)
% cell 1058 looks nice
% cell 3063
% 1037 + 1107

p=inputParser;
addParameter(p,'cell_id',[]);
addParameter(p,'method','wcorr');
parse(p,varargin{:});

load('.\Tables\subsets_of_cells.mat');
load(['rate_remapping_analysis_TRACK_PAIRS_' p.Results.method '.mat']);
load('folders_to_process_remapping.mat');
% epochs of interest
row_idx= strcmp([remapping.epoch],'POST');
all_cell_IDs= [remapping(row_idx).new_ID];
rate_diff= (remapping(row_idx).mean_max_FR_replay_diff);
peak_diff= (remapping(row_idx).place_field_diff);
stable_cells= subset_of_cells.cell_IDs{strcmp(subset_of_cells.subset,'stable cells laps')};
PRE_to_POST_ids= remapping(row_idx).PRE_to_POST_active_cells;

if isempty(varargin) % pseudo- randomly select neuron in upper quadrant
    % need both to be in the same session
    n_attempts= 1;
    sess_id= randi(max(round(all_cell_IDs/1000)),1);
    possible_idx= all_cell_IDs(rate_diff >=0.1*max(rate_diff) & peak_diff >=0.1*max(peak_diff) &...
                    round(all_cell_IDs'/1000)==sess_id & ismember(all_cell_IDs,stable_cells)' &  ismember(all_cell_IDs,PRE_to_POST_ids)');
    possible_idx2= all_cell_IDs(rate_diff <=  0.1*min(rate_diff) & peak_diff <= 0.1*min(peak_diff) &...
                round(all_cell_IDs'/1000)==sess_id  & ismember(all_cell_IDs,stable_cells)' &  ismember(all_cell_IDs,PRE_to_POST_ids)');
    while isempty(possible_idx) && n_attempts<100
        possible_idx= all_cell_IDs(rate_diff >=0.1*max(rate_diff) & peak_diff >=0.1*max(peak_diff) &...
                    round(all_cell_IDs'/1000)==sess_id & ismember(all_cell_IDs,stable_cells)' &  ismember(all_cell_IDs,PRE_to_POST_ids)');
        n_attempts= n_attempts+1;
    end
    while isempty(possible_idx2) && n_attempts<100
        possible_idx2= all_cell_IDs(rate_diff <=  0.1*min(rate_diff) & peak_diff <= 0.1*min(peak_diff) &...
                round(all_cell_IDs'/1000)==sess_id  & ismember(all_cell_IDs,stable_cells)' &  ismember(all_cell_IDs,PRE_to_POST_ids)');
        n_attempts= n_attempts+1;
    end
    cell_of_interest= possible_idx(randi(length(possible_idx),1));
    cell_of_interest2= possible_idx2(randi(length(possible_idx2),1));
else
    cell_of_interest= p.Results.cell_id(1);
    cell_of_interest2= p.Results.cell_id(2);
end
cell_of_interest
% disp('rate diff')
% rate_diff(all_cell_IDs==cell_of_interest)
% disp('peak diff')
% peak_diff(all_cell_IDs==cell_of_interest)
cell_of_interest2

% first create a cdata for that cell in subset_of_cells
start_row= height(subset_of_cells);
if ~sum(contains(subset_of_cells.subset,'example neuron'))
    new_row= start_row+1;
else
    new_row= find(contains(subset_of_cells.subset,'example neuron'));
end
subset_of_cells.subset{new_row}= 'example neuron';
subset_of_cells.cell_IDs{new_row}= stable_cells;
cdata= cell(1,length(remapping))'; 
% for POST
cdata{row_idx}=repmat([0.5 0.5 0.5],length(all_cell_IDs),1);
cdata{row_idx}(all_cell_IDs==cell_of_interest,1)=1; cdata{row_idx}(all_cell_IDs==cell_of_interest,2:3)=0; 
cdata{row_idx}(all_cell_IDs==cell_of_interest2,1:2)=0; cdata{row_idx}(all_cell_IDs==cell_of_interest2,3)=1; 
% for PRE
pre_row_idx= strcmp([remapping.epoch],'PRE');
pre_all_cell_IDs= [remapping(pre_row_idx).new_ID];
cdata{pre_row_idx}=repmat([0.5 0.5 0.5],length(pre_all_cell_IDs),1);
cdata{pre_row_idx}(pre_all_cell_IDs==cell_of_interest,1)=1; cdata{pre_row_idx}(pre_all_cell_IDs==cell_of_interest,2:3)=0; 
cdata{pre_row_idx}(pre_all_cell_IDs==cell_of_interest2,1:2)=0; cdata{pre_row_idx}(pre_all_cell_IDs==cell_of_interest2,3)=1; 
subset_of_cells.cdata{new_row}= cdata;
save('.\Tables\subsets_of_cells.mat','subset_of_cells');

% create corr with example neuron colour coded in red
plot_rate_remapping_NEW('subset','example neuron','epoch',{'PRE','POST'},'cdata',cdata,...
    'x_label','Peak in-field Rate Difference (Hz)','y_label','Peak Replay Rate Difference (Hz)');

% now modify for plotting raw session data
% find session 
sess_idx= floor(cell_of_interest/1000);
sess_cell_id= mod(cell_of_interest,1000);
sess_cell_id2= mod(cell_of_interest2,1000);
cd(folders{sess_idx});

load('extracted_place_fields_BAYESIAN.mat');
figure('Color','w');
subplot(1,2,1)
plot(place_fields_BAYESIAN.track(1).x_bin_centres,place_fields_BAYESIAN.track(1).raw{sess_cell_id},'k','LineWidth',1.5); hold on;
plot(place_fields_BAYESIAN.track(2).x_bin_centres,place_fields_BAYESIAN.track(2).raw{sess_cell_id},'Color',[0.5 0.5 0.5],'LineWidth',1.5);
xlabel('location (cm)'); ylabel('Firing Rate (Hz)'); title(['Cell ID: ' num2str(cell_of_interest)]);
subplot(1,2,2)
plot(place_fields_BAYESIAN.track(1).x_bin_centres,place_fields_BAYESIAN.track(1).raw{sess_cell_id2},'k','LineWidth',1.5); hold on;
plot(place_fields_BAYESIAN.track(2).x_bin_centres,place_fields_BAYESIAN.track(2).raw{sess_cell_id2},'Color',[0.5 0.5 0.5],'LineWidth',1.5);
xlabel('location (cm)'); ylabel('Firing Rate (Hz)'); title(['Cell ID: ' num2str(cell_of_interest2)]);

plot_color_line= [{repmat([0.5 0.5 0.5],length(place_fields_BAYESIAN.track(1).sorted_good_cells),1)},...
    {repmat([0.5 0.5 0.5],length(place_fields_BAYESIAN.track(2).sorted_good_cells),1)}];
plot_color_line{1}(place_fields_BAYESIAN.track(1).sorted_good_cells== mod(cell_of_interest,1000),1)=1;
plot_color_line{1}(place_fields_BAYESIAN.track(1).sorted_good_cells== mod(cell_of_interest,1000),2:3)=0;
plot_color_line{1}(place_fields_BAYESIAN.track(1).sorted_good_cells== mod(cell_of_interest2,1000),1:2)=0;
plot_color_line{1}(place_fields_BAYESIAN.track(1).sorted_good_cells== mod(cell_of_interest2,1000),3)=1;
plot_color_line{2}(place_fields_BAYESIAN.track(2).sorted_good_cells== mod(cell_of_interest,1000),1)=1;
plot_color_line{2}(place_fields_BAYESIAN.track(2).sorted_good_cells== mod(cell_of_interest,1000),2:3)=0;
plot_color_line{2}(place_fields_BAYESIAN.track(2).sorted_good_cells== mod(cell_of_interest2,1000),1:2)=0;
plot_color_line{2}(place_fields_BAYESIAN.track(2).sorted_good_cells== mod(cell_of_interest2,1000),3)=1;
save('example_neuron_cdata.mat','plot_color_line');
data_view;

% cd ..

% get linear positions
load('extracted_position.mat');
figure('Color','w');
plot(position.t,position.linear(1).linear,'k','LineWidth',1.5);
xlim([min(position.t)+[0.2 0.5].*(max(position.t)-min(position.t))]);
ylim([-50 250]);
figure('Color','w');
plot(position.t,position.linear(2).linear,'Color',[0.5 0.5 0.5],'LineWidth',1.5);
xlim([min(position.t)+[0.2 0.5].*(max(position.t)-min(position.t))]);
ylim([-50 250]);

% get mean max FR for each cell
FR_cell1(1)= remapping_raw(row_idx).track1_mean_replay_inst_FR_nonZero(all_cell_IDs==cell_of_interest);
FR_cell1(2)= remapping_raw(row_idx).track2_mean_replay_inst_FR_nonZero(all_cell_IDs==cell_of_interest);
FR_cell2(1)= remapping_raw(row_idx).track1_mean_replay_inst_FR_nonZero(all_cell_IDs==cell_of_interest2);
FR_cell2(2)= remapping_raw(row_idx).track2_mean_replay_inst_FR_nonZero(all_cell_IDs==cell_of_interest2);

end