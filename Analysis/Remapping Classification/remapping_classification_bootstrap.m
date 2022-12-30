function remapping_classification_bootstrap(folders)
parameters= list_of_parameters;
master_folder= pwd;
cell_remapping= table;
% load('cell_remapping.mat');

for this_folder=1:size(folders,1)
    cd([master_folder '\' folders{this_folder,1}]);
    disp([master_folder '\' folders{this_folder,1}])
    tic;
    
    load('lap_times.mat','lap_times');
    load('extracted_place_fields.mat','place_fields');
    load('extracted_place_fields_BAYESIAN.mat','place_fields_BAYESIAN');
%     all_tracks=1:length(place_fields_BAYESIAN.track);

    good_cells= place_fields_BAYESIAN.good_place_cells;
    
    if contains(pwd,'Rate remapping')
        session_id= folders(this_folder,1);
        rat_id= folders(this_folder,1); % repeat because there is no ratID for remapping
    else
        fd= folders{this_folder};
        fd= strsplit(fd, '\');
        session_id= fd(2);
        rat_id= fd(1);
    end
        
    % do the bootstrapping
    % the idea is to take subsets of laps and calculate place fields
    % to get a distribution of peak FR and locations
    % choose which laps will be used
    n_boot= 1000;
    track_laps_ts= cell(length(lap_times),n_boot);
  
    for this_track=1:length(lap_times)
        % we want full laps ! (back and forth!)
        available_laps= floor(lap_times(this_track).total_number_of_laps/2);   
       % for bootstrapping you need to resample the exact number of samples you have each time
       % so if you have 20 laps you need to randomly choose 20 laps with replacement
       full_laps_idx= [(1:2:2*available_laps-1)' (2:2:2*available_laps)'];
 
       for i=1:n_boot
            rand_idx= randi(available_laps,1, available_laps);
            % need to remove duplicates and order
            rand_idx= sort(unique(rand_idx));
            lap_idx= full_laps_idx(rand_idx,:);
            lap_idx= sort(lap_idx(:));
            track_laps_ts{this_track,i}= [lap_times(this_track).start(lap_idx)'  lap_times(this_track).end(lap_idx)' ];
       end
        % OLD VERSION 
%         available_laps= floor(lap_times(this_track).total_number_of_laps/2);   
%         typ= mod(length(lap_times(this_track).lap),2);
%             if typ % odd number of laps
%                 even_laps=4:2:length(lap_times(this_track).lap)-1;
%                 odd_laps=3:2:length(lap_times(this_track).lap);
%             else % even number laps
%                 even_laps=4:2:length(lap_times(this_track).lap);
%                 odd_laps=3:2:length(lap_times(this_track).lap)-1;
%             end
%             % create all combinations of even and odd laps
%             subset_laps1= nchoosek(even_laps,3) ; % had to reduce.... to 3 full laps
%             subset_laps2=  nchoosek(odd_laps,3);
%             n_boot= 100; % given our dataset, max is 200
%             % create all combinations of vectors
%             all_vec_comb= combvec(1:size(subset_laps1,1),1:size(subset_laps2,1));
%             % take a subset
%             sub_idx= randi(size(all_vec_comb,2),1,n_boot);
%             subset_laps= [subset_laps1(all_vec_comb(1,sub_idx),:) subset_laps2(all_vec_comb(2,sub_idx),:)];
%             for i=1:size(subset_laps,1)
%                 subset_laps(i,:)= sort(subset_laps(i,:));
%                 track_laps_ts{this_track,i}= [ lap_times(this_track).start(subset_laps(i,:))'  lap_times(this_track).end(subset_laps(i,:))' ]; 
%             end
    end
    
    
    % keep track of which cells are 'good cells' on which tracks
    % hard coded but there should be a better option with nchoosek(n,k),
    % with k=2:max number tracks
    track_pairs=  nchoosek(1:length(lap_times),2);
    common_cells_tmp= zeros(length(good_cells),length(lap_times));
    for this_track=1:length(lap_times)
          [~,cell_idx,~]= multintersect([{good_cells} arrayfun(@(x) place_fields_BAYESIAN.track(x).good_cells,this_track,'UniformOutput',0)]);
          common_cells_tmp(cell_idx,this_track)= 1;
    end
    
     cell_remapping_tmp= table(repmat(rat_id,length(good_cells),1),repmat(session_id,length(good_cells),1),good_cells',...
                                                     cell(length(good_cells),1),cell(length(good_cells),1),cell(length(good_cells),1),...
                                                     cell(length(good_cells),1),cell(length(good_cells),1),cell(length(good_cells),1),...
                                                     cell(length(good_cells),1),cell(length(good_cells),1),cell(length(good_cells),1),cell(length(good_cells),1),cell(length(good_cells),1),...
                                                     cell(length(good_cells),1),cell(length(good_cells),1),cell(length(good_cells),1),cell(length(good_cells),1),cell(length(good_cells),1),...
                                                     cell(length(good_cells),1),cell(length(good_cells),1),cell(length(good_cells),1),cell(length(good_cells),1),cell(length(good_cells),1),...
                                                     cell(length(good_cells),1),cell(length(good_cells),1),cell(length(good_cells),1),cell(length(good_cells),1),cell(length(good_cells),1),cell(length(good_cells),1),cell(length(good_cells),1),cell(length(good_cells),1),...
                                                      'VariableNames',{'rat','day','cell_id',...
                                                      'spatially_tuned','stability','place_fields'...
                                                      'track_FR','track_Centre','track_COM',...
                                                      'rate_modulation','centre_modulation','COM_modulation','overlap_modulation','between_overlap'...
                                                      'peak_FR','centre','COM','pl_field','pl_field_norm','overlap'...
                                                      'FR_dist','centre_dist','COM_dist','overlap_dist',...
                                                      'FR_median','FR_CI','centre_median','centre_CI','COM_median','COM_CI','overlap_median','overlap_CI'});
                     
      % since n_boot might change with the track, we need to do per
      % track ?
     tic;
     num_tracks= length(lap_times);
     num_cells= length(good_cells);
     parfor i=1:size(track_laps_ts,2)
         
         s{i}= table(cell(num_cells,1),cell(num_cells,1),cell(num_cells,1),cell(num_cells,1),cell(num_cells,1),cell(num_cells,1),...
             'VariableNames',{'spatially_tuned','peak_FR','centre','COM','pl_field','pl_field_norm'});
%          
         
             track_laps_ts_tmp=  track_laps_ts(:,i)';
%              track_idx_this_i= find(~cellfun(@isempty,track_laps_ts_tmp)) ;
             pl_fld_tmp= calculate_place_fields_epochs(parameters.x_bins_width,track_laps_ts_tmp);
                
            % get distr peak,COM,centre, field (for overlap)
            
             for this_cell=1: length(good_cells)
                    peak=  arrayfun(@(x)pl_fld_tmp.track(x).peak(good_cells(this_cell)),1:num_tracks)';
                    centre=  arrayfun(@(x)pl_fld_tmp.track(x).centre(good_cells(this_cell)),1:num_tracks)';
                    COM= arrayfun(@(x) pl_fld_tmp.track(x).centre_of_mass(good_cells(this_cell)),1:num_tracks)';
                    fld= arrayfun(@(x) pl_fld_tmp.track(x).smooth{good_cells(this_cell)},1:num_tracks,'UniformOutput',0)';
                    fld_norm= arrayfun(@(x) pl_fld_tmp.track(x).smooth{good_cells(this_cell)}/peak(x),1:num_tracks,'UniformOutput',0)';
                    nan_idx=  peak ==0;
%                     if sum(nan_idx)>=1
%                         disp('no field found in subset of five laps')
%                         keyboard;
%                     end
                     % if didn't fire at all, do not count - although this
                     % could also say something ?
%                     peak(nan_idx)= NaN;
                    centre(nan_idx)= NaN;
                    COM(nan_idx)= NaN;
                    fld(nan_idx)= {NaN};
                    fld_norm(nan_idx)= {NaN};
                    
                    s{i}.spatially_tuned{this_cell}= common_cells_tmp(this_cell,:);
                    s{i}.peak_FR{this_cell}= [s{i}.peak_FR{this_cell}  peak];
                    s{i}.centre{this_cell}= [s{i}.centre{this_cell}  centre];
                    s{i}.COM{this_cell}= [s{i}.COM{this_cell} COM];
                    s{i}.pl_field{this_cell}= [s{i}.pl_field{this_cell} fld];
                    s{i}.pl_field_norm{this_cell}= [s{i}.pl_field_norm{this_cell} fld_norm];
   
                    
             end
     end
    disp(['time to get ' num2str(n_boot) 'fields took... ']);
    toc
    for i=1:size(track_laps_ts,2)
        cell_remapping_tmp.spatially_tuned= s{1}.spatially_tuned;
        cell_remapping_tmp.peak_FR= arrayfun(@(x) [cell_remapping_tmp.peak_FR{x} s{i}.peak_FR{x}],1:height(cell_remapping_tmp),'UniformOutput',0)';
        cell_remapping_tmp.centre= arrayfun(@(x) [cell_remapping_tmp.centre{x} s{i}.centre{x}],1:height(cell_remapping_tmp),'UniformOutput',0)';
        cell_remapping_tmp.COM= arrayfun(@(x) [cell_remapping_tmp.COM{x} s{i}.COM{x}],1:height(cell_remapping_tmp),'UniformOutput',0)';
        cell_remapping_tmp.pl_field= arrayfun(@(x) [cell_remapping_tmp.pl_field{x} s{i}.pl_field{x}],1:height(cell_remapping_tmp),'UniformOutput',0)';
        cell_remapping_tmp.pl_field_norm= arrayfun(@(x) [cell_remapping_tmp.pl_field_norm{x} s{i}.pl_field_norm{x}],1:height(cell_remapping_tmp),'UniformOutput',0)';
    end
    
    %%% now do the classification

    % create distributions from the N samples
    % find out about stability using overlap
    for this_cell=1: length(good_cells)
            for this_track=1:length(lap_times) % in case tracks are not the same length
                pos_edges= place_fields.track(this_track).x_bin_edges;
                pos_ctrs= place_fields.track(this_track).x_bin_centres;
                rate_edges= [0:50];
                rate_ctrs= rate_edges(1:end-1)+0.05;
                
                % get place fields from whole track experience
                cell_remapping_tmp.place_fields{this_cell}{this_track}= place_fields.track(this_track).smooth{good_cells(this_cell)};
                cell_remapping_tmp.track_FR{this_cell}(this_track)= place_fields.track(this_track).peak(good_cells(this_cell));
                cell_remapping_tmp.track_Centre{this_cell}(this_track)= place_fields.track(this_track).centre(good_cells(this_cell));
                cell_remapping_tmp.track_COM{this_cell}(this_track)= place_fields.track(this_track).centre_of_mass(good_cells(this_cell));
                
                
                % overlap - has to be done now when all the ratemaps are
                % available
                all_comb_overlap= nchoosek(1:length(cell_remapping_tmp.pl_field{this_cell}(this_track,:)),2);
                % we are going to sub-sample a bit to reduce compute time -
                % keep only 10 000 fields
                all_comb_overlap= all_comb_overlap(randsample(size(all_comb_overlap,1),10000,1),:);
                for j=1:size(all_comb_overlap,1) % get all possible overlaps
                    f1= cell_remapping_tmp.pl_field{this_cell}{this_track,all_comb_overlap(j,1)};
                    f2= cell_remapping_tmp.pl_field{this_cell}{this_track,all_comb_overlap(j,2)};
                    cell_remapping_tmp.overlap{this_cell}(this_track,j)= sum(min(f1,f2))/(sum(f1)+sum(f2)-sum(min(f1,f2)));
                end
                % get distr etc as normal
                overlap_edges= [0:0.1:1];
                overlap_ctrs=  overlap_edges(1:end-1)+0.05;
                cell_remapping_tmp.overlap_dist{this_cell}(this_track)= {histcounts(cell_remapping_tmp.overlap{this_cell}(this_track,:),overlap_edges,'Normalization','pdf')};
                cell_remapping_tmp.overlap_median{this_cell}(this_track)= nanmedian(cell_remapping_tmp.overlap{this_cell}(this_track,:));
                [~,cell_remapping_tmp.overlap_CI{this_cell}(this_track,1)]= get_percentiles_pdf(cell_remapping_tmp.overlap_dist{this_cell}{this_track},overlap_ctrs,5);
                [~,cell_remapping_tmp.overlap_CI{this_cell}(this_track,2)]= get_percentiles_pdf(cell_remapping_tmp.overlap_dist{this_cell}{this_track},overlap_ctrs,95);
                
                cell_remapping_tmp.stability{this_cell}(this_track)= cell_remapping_tmp.overlap_median{this_cell}(this_track) > 0.5;
            end
    end
    
    % still calculate even if not "stable"
    for this_cell=1: length(good_cells)
            for this_track=1:length(lap_times) % in case tracks are not the same length
%                 if cell_remapping_tmp.stability{this_cell}(this_track) % if stable
                    % get PDFs
                    cell_remapping_tmp.FR_dist{this_cell}(this_track)= {histcounts(cell_remapping_tmp.peak_FR{this_cell}(this_track,:),rate_edges,'Normalization','pdf')};
                    cell_remapping_tmp.centre_dist{this_cell}(this_track)= {histcounts(cell_remapping_tmp.centre{this_cell}(this_track,:),pos_edges,'Normalization','pdf')};
                    cell_remapping_tmp.COM_dist{this_cell}(this_track)= {histcounts(cell_remapping_tmp.COM{this_cell}(this_track,:),pos_edges,'Normalization','pdf')};

                    % get median and CI
                    cell_remapping_tmp.FR_median{this_cell}(this_track)= nanmedian(cell_remapping_tmp.peak_FR{this_cell}(this_track,:));
                    [~,cell_remapping_tmp.FR_CI{this_cell}(this_track,1)]= get_percentiles_pdf(cell_remapping_tmp.FR_dist{this_cell}{this_track},rate_ctrs,5);
                    [~,cell_remapping_tmp.FR_CI{this_cell}(this_track,2)]= get_percentiles_pdf(cell_remapping_tmp.FR_dist{this_cell}{this_track},rate_ctrs,95);

                    cell_remapping_tmp.centre_median{this_cell}(this_track)= nanmedian(cell_remapping_tmp.centre{this_cell}(this_track,:));
                    [~,cell_remapping_tmp.centre_CI{this_cell}(this_track,1)]= get_percentiles_pdf(cell_remapping_tmp.centre_dist{this_cell}{this_track},pos_ctrs,5);
                    [~,cell_remapping_tmp.centre_CI{this_cell}(this_track,2)]= get_percentiles_pdf(cell_remapping_tmp.centre_dist{this_cell}{this_track},pos_ctrs,95);

                    cell_remapping_tmp.COM_median{this_cell}(this_track)= nanmedian(cell_remapping_tmp.COM{this_cell}(this_track,:));
                    [~,cell_remapping_tmp.COM_CI{this_cell}(this_track,1)]= get_percentiles_pdf(cell_remapping_tmp.COM_dist{this_cell}{this_track},pos_ctrs,5);
                    [~,cell_remapping_tmp.COM_CI{this_cell}(this_track,2)]= get_percentiles_pdf(cell_remapping_tmp.COM_dist{this_cell}{this_track},pos_ctrs,95);
%                 else % not stable
%                     cell_remapping_tmp.FR_dist{this_cell}(this_track)= {NaN(1,length(rate_edges)-1)};
%                     cell_remapping_tmp.centre_dist{this_cell}(this_track)= {NaN(1,length(pos_edges)-1)};
%                     cell_remapping_tmp.COM_dist{this_cell}(this_track)= {NaN(1,length(pos_edges)-1)};
%                     cell_remapping_tmp.FR_median{this_cell}(this_track)= NaN;
%                     cell_remapping_tmp.FR_CI{this_cell}(this_track,:)=  [NaN NaN];
%                     cell_remapping_tmp.centre_median{this_cell}(this_track)= NaN;
%                     cell_remapping_tmp.centre_CI{this_cell}(this_track,:)= [NaN, NaN];
%                     cell_remapping_tmp.COM_median{this_cell}(this_track)= NaN;
%                     cell_remapping_tmp.COM_CI{this_cell}(this_track,:)=   [NaN NaN];
%                 end
            end
    
            % classify
           for this_pair=1:size(track_pairs,1)
                % if either mean is larger/smaller than the 95th/5th
                % percentiles
                
                 % + compare to overlap between tracks
                 f1= place_fields.track(track_pairs(this_pair,1)).smooth{good_cells(this_cell)};
                 f2= place_fields.track(track_pairs(this_pair,2)).smooth{good_cells(this_cell)};
                 cell_remapping_tmp.between_overlap{this_cell}(this_pair)= sum(min(f1,f2))/(sum(f1)+sum(f2)-sum(min(f1,f2)));
                
             % FR - rate modulation
             if sum(cell_remapping_tmp.spatially_tuned{this_cell}(track_pairs(this_pair,:))) > 1
                  if cell_remapping_tmp.FR_median{this_cell}(track_pairs(this_pair,1)) < cell_remapping_tmp.FR_CI{this_cell}(track_pairs(this_pair,2),1) || ...
                     cell_remapping_tmp.FR_median{this_cell}(track_pairs(this_pair,2)) < cell_remapping_tmp.FR_CI{this_cell}(track_pairs(this_pair,1),1) || ...
                     cell_remapping_tmp.FR_median{this_cell}(track_pairs(this_pair,1)) > cell_remapping_tmp.FR_CI{this_cell}(track_pairs(this_pair,2),2) || ...
                     cell_remapping_tmp.FR_median{this_cell}(track_pairs(this_pair,2)) > cell_remapping_tmp.FR_CI{this_cell}(track_pairs(this_pair,1),2)

                            cell_remapping_tmp.rate_modulation{this_cell}(this_pair)=1;
                  else
                            cell_remapping_tmp.rate_modulation{this_cell}(this_pair)=0;
                  end
                              
                 % Centre + COM - location modulation
                 if cell_remapping_tmp.centre_median{this_cell}(track_pairs(this_pair,1)) < cell_remapping_tmp.centre_CI{this_cell}(track_pairs(this_pair,2),1) || ...
                     cell_remapping_tmp.centre_median{this_cell}(track_pairs(this_pair,2)) < cell_remapping_tmp.centre_CI{this_cell}(track_pairs(this_pair,1),1) || ...
                     cell_remapping_tmp.centre_median{this_cell}(track_pairs(this_pair,1)) > cell_remapping_tmp.centre_CI{this_cell}(track_pairs(this_pair,2),2) || ...
                     cell_remapping_tmp.centre_median{this_cell}(track_pairs(this_pair,2)) > cell_remapping_tmp.centre_CI{this_cell}(track_pairs(this_pair,1),2) 

                            cell_remapping_tmp.centre_modulation{this_cell}(this_pair)=1;
                  else
                            cell_remapping_tmp.centre_modulation{this_cell}(this_pair)=0;
                 end
                  if cell_remapping_tmp.COM_median{this_cell}(track_pairs(this_pair,1)) < cell_remapping_tmp.COM_CI{this_cell}(track_pairs(this_pair,2),1) || ...
                     cell_remapping_tmp.COM_median{this_cell}(track_pairs(this_pair,2)) < cell_remapping_tmp.COM_CI{this_cell}(track_pairs(this_pair,1),1) || ...
                     cell_remapping_tmp.COM_median{this_cell}(track_pairs(this_pair,1)) > cell_remapping_tmp.COM_CI{this_cell}(track_pairs(this_pair,2),2) || ...
                     cell_remapping_tmp.COM_median{this_cell}(track_pairs(this_pair,2)) > cell_remapping_tmp.COM_CI{this_cell}(track_pairs(this_pair,1),2) 

                            cell_remapping_tmp.COM_modulation{this_cell}(this_pair)=1;
                  else
                            cell_remapping_tmp.COM_modulation{this_cell}(this_pair)=0;
                  end

                 % overlap - location modulation, alternative method
                 % we consider only if the overlap between is smaller than whithin
                 % (whithin is a measure of stability)
                 if  cell_remapping_tmp.between_overlap{this_cell}(this_pair)< cell_remapping_tmp.overlap_CI{this_cell}(track_pairs(this_pair,2),1) || ...
                      cell_remapping_tmp.between_overlap{this_cell}(this_pair) < cell_remapping_tmp.overlap_CI{this_cell}(track_pairs(this_pair,1),1)
                            cell_remapping_tmp.overlap_modulation{this_cell}(this_pair)=1;
                  else
                            cell_remapping_tmp.overlap_modulation{this_cell}(this_pair)=0;
                 end
             
             else
                  cell_remapping_tmp.rate_modulation{this_cell}(this_pair)=NaN;
                  cell_remapping_tmp.centre_modulation{this_cell}(this_pair)=NaN;
                  cell_remapping_tmp.COM_modulation{this_cell}(this_pair)=NaN;
                  cell_remapping_tmp.overlap_modulation{this_cell}(this_pair)=NaN;
             end
              
        end
    end  
toc
                
        cell_remapping= [cell_remapping; cell_remapping_tmp]; 
        clearvars -except folders this_folder master_folder cell_remapping parameters
  
end

% corrective option - as is will calculate overlap even for non spatially
% tuned cells
for this_cell= 1:height (cell_remapping)
    non_tuned = find(cell_remapping.spatially_tuned{this_cell} == 0);
    if ~isempty(non_tuned)
        cell_remapping.stability{this_cell}(non_tuned)= 0;
        cell_remapping.overlap_median{this_cell}(non_tuned)= NaN;
        cell_remapping.overlap_CI{this_cell}(non_tuned,:)= NaN;
%         cell_remapping.overlap_dist{this_cell}{non_tuned}= NaN;
            for kk= 1:length(non_tuned)
                 pairs_to_remove= find(sum(ismember(nchoosek(1:3,2),non_tuned(kk)),2));
                 cell_remapping.between_overlap{this_cell}(pairs_to_remove)= NaN;
            end
    end
end

cd(master_folder)
if exist('./Tables')~=7
    mkdir('./Tables/');
end
save('./Tables/cell_remapping.mat','cell_remapping','-v7.3');


%% now select cells from the remapping corr analysis
load('rate_remapping_analysis_TRACK_PAIRS_wcorr.mat');
cell_remapping_incl_cells= cell_remapping;
cell_remapping_incl_cells(:,[2 5 13:end])= []; % remove less necessary columns to make variable lighter
clear cell_remapping;
unique_sess= unique(cell_remapping_incl_cells.rat,'stable');
for this_session= 1:length(unique_sess)
%     incl_cells= [];
%     for this_epoch=1:length(remapping)
%         incl_cells= [incl_cells remapping(this_epoch).common_good_cells(remapping(this_epoch).experiment == this_session)'];
%     end
%     incl_cells= unique(incl_cells);
    incl_cells= remapping(1).common_good_cells(remapping(1).experiment_all == this_session); % epoch doesn't matter
    session_idx= find(strcmp(cell_remapping_incl_cells.rat,remapping(1).folder{this_session}));
    cells_in_table= cell_remapping_incl_cells.cell_id(session_idx);
    
    cell_remapping_incl_cells(session_idx(~ismember(cells_in_table,incl_cells)),:)=[];
    
    if strcmp(remapping(1).folder{this_session},'N-BLU_Day7_Ctrl-16x30') % for this session only keep from first two tracks
%         session_idx= find(contains(cell_remapping_incl_cells.rat,remapping(1).folder{this_session}));
%         tuned= vertcat(cell_remapping_incl_cells.spatially_tuned{session_idx});
%         to_exclude= sum(tuned(:,1:2),2) ~= 2;
%         cell_remapping_incl_cells(session_idx(to_exclude),:)=[];
%         % new indices
%         session_idx= find(contains(cell_remapping_incl_cells.rat,remapping(1).folder{this_session}));
        for this_idx=1:length(session_idx)
            cell_remapping_incl_cells.spatially_tuned{session_idx(this_idx)}=  cell_remapping_incl_cells.spatially_tuned{session_idx(this_idx)}(1:2);
            cell_remapping_incl_cells.place_fields{session_idx(this_idx)}=  cell_remapping_incl_cells.place_fields{session_idx(this_idx)}(1:2);
            cell_remapping_incl_cells.track_FR{session_idx(this_idx)}=  cell_remapping_incl_cells.track_FR{session_idx(this_idx)}(1:2);
            cell_remapping_incl_cells.track_Centre{session_idx(this_idx)}=  cell_remapping_incl_cells.track_FR{session_idx(this_idx)}(1:2);
            cell_remapping_incl_cells.track_COM{session_idx(this_idx)}=  cell_remapping_incl_cells.track_FR{session_idx(this_idx)}(1:2);
            cell_remapping_incl_cells.rate_modulation{session_idx(this_idx)}=  cell_remapping_incl_cells.rate_modulation{session_idx(this_idx)}(1);
            cell_remapping_incl_cells.centre_modulation{session_idx(this_idx)}=  cell_remapping_incl_cells.centre_modulation{session_idx(this_idx)}(1);
            cell_remapping_incl_cells.COM_modulation{session_idx(this_idx)}=  cell_remapping_incl_cells.COM_modulation{session_idx(this_idx)}(1);
        end
    end
end

cell_remapping_incl_cells.new_ID = create_new_remapping_IDs('wcorr',cell_remapping_incl_cells.rat,cell_remapping_incl_cells.cell_id)';
cell_remapping_incl_cells = movevars(cell_remapping_incl_cells, 'new_ID', 'After', 'cell_id');
save('./Tables/cell_remapping_incl_cells.mat','cell_remapping_incl_cells','-v7.3');

%% plot corr with colour code
rate_mod_cells= [cell_remapping_incl_cells.rate_modulation{:}] ==1;
non_rate_mod_cells= [cell_remapping_incl_cells.rate_modulation{:}] ==0;
glob_cells=  [cell_remapping_incl_cells.centre_modulation{:}] ==1;
non_glob_cells=  [cell_remapping_incl_cells.centre_modulation{:}] ==0;
pure_rate_cells= find(rate_mod_cells & non_glob_cells);
cdata= zeros(height(cell_remapping_incl_cells),3);
cdata(rate_mod_cells,2)=0.8;
cdata(non_rate_mod_cells,:)= 0.5;
cdata(pure_rate_cells,1)= 1;
cdata(pure_rate_cells,2)= 0;

if exist('.\Tables\subsets_of_cells.mat')
    load('.\Tables\subsets_of_cells.mat');
    start_row= height(subset_of_cells);
else
    start_row=0;
    subset_of_cells= table(cell(1,1),cell(1,1),cell(1,1),'VariableNames',{'subset','cell_IDs','cdata'});
end

subset_of_cells.subset{start_row+1}= 'rate modulated cells';
subset_of_cells.cell_IDs{start_row+1}= cell_remapping_incl_cells.new_ID(rate_mod_cells);
subset_of_cells.subset{start_row+2}= 'non rate modulated cells';
subset_of_cells.cell_IDs{start_row+2}= cell_remapping_incl_cells.new_ID(non_rate_mod_cells);
subset_of_cells.subset{start_row+3}= 'pure rate modulated cells';
subset_of_cells.cell_IDs{start_row+3}= cell_remapping_incl_cells.new_ID(pure_rate_cells);
subset_of_cells.subset{start_row+4}= 'global_mod_cells';
subset_of_cells.cell_IDs{start_row+4}= cell_remapping_incl_cells.new_ID(glob_cells);
save('.\Tables\subsets_of_cells.mat','subset_of_cells');

end