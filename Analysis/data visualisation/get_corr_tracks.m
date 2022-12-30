function get_corr_tracks(folders,varargin)

p= inputParser;
addParameter(p,'epoch','first_second_half',@ischar);
addParameter(p,'method','pop_vector',@ischar);
parse(p,varargin{:});

parameters= list_of_parameters;
master_folder= pwd;
for this_folder= 1:length(folders)
     cd([master_folder '\' folders{this_folder}]);
     disp(['       ' folders{this_folder}])
     
     load('extracted_place_fields_BAYESIAN.mat');
     
     switch p.Results.epoch
         case 'even_odd'
             if exist('place_fields_even_odd_laps.mat') ~=2
                 [place_fields_even,place_fields_odd]= calculate_place_fields_laps('grouping','even_odd','size',parameters.x_bins_width_bayesian,'save_option',0);
                 save('place_fields_even_odd_laps.mat','place_fields_even','place_fields_odd');
             else
                 load('place_fields_even_odd_laps.mat');
             end
         case 'first_second_half'
             load('lap_times.mat');
             for this_track=1:length(lap_times)
                middle_lap= round(lap_times(this_track).total_number_of_laps/2); % roughly 50% of laps
                % find starts & stops 
                h1{this_track,:}= [lap_times(this_track).start(1) lap_times(this_track).end(middle_lap-1) ];
                h2{this_track,:}= [lap_times(this_track).start(middle_lap) lap_times(this_track).end(end)];
             end
            % calculate place fields for each half
            pl_fld_h1= calculate_place_fields_epochs(parameters.x_bins_width_bayesian,h1);
            pl_fld_h2= calculate_place_fields_epochs(parameters.x_bins_width_bayesian,h2);
     end
  
     % between tracks
     track_pairs= nchoosek(1:2,2);
     good_cells= place_fields_BAYESIAN.good_place_cells;
     for this_pair= 1:size(track_pairs,1)
         track1= track_pairs(this_pair,1);
         track2= track_pairs(this_pair,2);
          
            switch p.Results.epoch
                case 'even_odd'
                    switch p.Results.method
                            case 'ratemaps'
                             % correlate even vs odd
                             even_odd_pair_corr(this_folder,this_pair)= corr(cell2mat(place_fields_even.track(track1).raw(good_cells))',...
                                cell2mat(place_fields_odd.track(track2).raw(good_cells))','type','Pearson');
                             % corr odd vs even
                             odd_even_pair_corr(this_folder,this_pair)= corr(cell2mat(place_fields_odd.track(track1).raw(good_cells))',...
                                cell2mat(place_fields_even.track(track2).raw(good_cells))','type','Pearson');
                            case 'pop_vector'
                                rmap_even= cell2mat(place_fields_even.track(track1).raw(good_cells)');
                                rmap_odd= cell2mat(place_fields_odd.track(track1).raw(good_cells)');
                                
                                even_odd_pair_mean_corr(this_folder,this_pair)= mean(arrayfun(@(x) corr(rmap_even(:,x),rmap_odd(:,x),'type','Pearson'),1:size(rmap_even,2)));
                                odd_even_pair_mean_corr(this_folder,this_pair)= mean(arrayfun(@(x) corr(rmap_odd(:,x),rmap_even(:,x),'type','Pearson'),1:size(rmap_odd,2)));
                                
                                even_odd_pair_std_corr(this_folder,this_pair)= std(arrayfun(@(x) corr(rmap_even(:,x),rmap_odd(:,x),'type','Pearson'),1:size(rmap_even,2)));
                                odd_even_pair_std_corr(this_folder,this_pair)= std(arrayfun(@(x) corr(rmap_odd(:,x),rmap_even(:,x),'type','Pearson'),1:size(rmap_odd,2)));
                    end
                case 'first_second_half'
                    switch p.Results.method
                            case 'ratemaps'
                            % first vs second
                            first_second_pair_corr(this_folder,this_pair)= corr(cell2mat(pl_fld_h1.track(track1).raw(good_cells))',...
                                cell2mat(pl_fld_h2.track(track2).raw(good_cells))','type','Pearson');
                            % second vs first
                            second_first_pair_corr(this_folder,this_pair)= corr(cell2mat(pl_fld_h2.track(track1).raw(good_cells))',...
                                cell2mat(pl_fld_h1.track(track2).raw(good_cells))','type','Pearson');
                            case 'pop_vector'
                                rmap1_1= cell2mat(pl_fld_h1.track(track1).raw(good_cells)');
                                rmap1_2= cell2mat(pl_fld_h2.track(track1).raw(good_cells)');
                                rmap2_1= cell2mat(pl_fld_h1.track(track2).raw(good_cells)');
                                rmap2_2= cell2mat(pl_fld_h2.track(track2).raw(good_cells)');
                                
                                first_second_pair_mean_corr(this_folder,this_pair)= mean(arrayfun(@(x) corr(rmap1_1(:,x),rmap2_2(:,x),'type','Pearson'),1:size(rmap1_1,2)));
                                second_first_pair_mean_corr(this_folder,this_pair)= mean(arrayfun(@(x) corr(rmap1_2(:,x),rmap2_1(:,x),'type','Pearson'),1:size(rmap1_1,2)));
                                first_first_pair_mean_corr(this_folder,this_pair)= mean(arrayfun(@(x) corr(rmap1_1(:,x),rmap2_1(:,x),'type','Pearson'),1:size(rmap1_1,2)));
                                second_second_pair_mean_corr(this_folder,this_pair)= mean(arrayfun(@(x) corr(rmap1_2(:,x),rmap2_2(:,x),'type','Pearson'),1:size(rmap1_1,2)));
                                
                                first_second_pair_std_corr(this_folder,this_pair)= std(arrayfun(@(x) corr(rmap1_1(:,x),rmap2_2(:,x),'type','Pearson'),1:size(rmap1_1,2)));
                                second_first_pair_std_corr(this_folder,this_pair)= std(arrayfun(@(x) corr(rmap1_2(:,x),rmap2_1(:,x),'type','Pearson'),1:size(rmap1_1,2)));
                                first_first_pair_std_corr(this_folder,this_pair)= std(arrayfun(@(x) corr(rmap1_1(:,x),rmap2_1(:,x),'type','Pearson'),1:size(rmap1_1,2)));
                                second_second_pair_std_corr(this_folder,this_pair)= std(arrayfun(@(x) corr(rmap1_2(:,x),rmap2_2(:,x),'type','Pearson'),1:size(rmap1_1,2)));
                               
                    end
            end
                    switch p.Results.method
                        case 'ratemaps'
                            % all laps
                            pair_corr(this_folder,this_pair)= corr(cell2mat(place_fields_BAYESIAN.track(track_pairs(this_pair,1)).raw(good_cells))',...
                                cell2mat(place_fields_BAYESIAN.track(track_pairs(this_pair,2)).raw(good_cells))','type','Pearson');
                        case 'pop_vector'
                            rmap_1= cell2mat(place_fields_BAYESIAN.track(track_pairs(this_pair,1)).raw(good_cells)');
                            rmap_2= cell2mat(place_fields_BAYESIAN.track(track_pairs(this_pair,2)).raw(good_cells)');
                         
                            pair_mean_corr(this_folder,this_pair)= mean(arrayfun(@(x) corr(rmap_1(:,x),rmap_2(:,x),'type','Pearson'),1:size(rmap_2,2)));
                            pair_std_corr(this_folder,this_pair)= std(arrayfun(@(x) corr(rmap_1(:,x),rmap_2(:,x),'type','Pearson'),1:size(rmap_2,2)));
                    end
              
     end
     
     % within track
     for this_track=1:2
         switch p.Results.epoch
                case 'even_odd'
                switch p.Results.method
                        case 'ratemaps'
                        % corr for track
                        even_odd_track_corr(this_folder,this_track)= corr(cell2mat(place_fields_even.track(this_track).raw(good_cells))',...
                            cell2mat(place_fields_odd.track(this_track).raw(good_cells))','type','Pearson');
                         % corr odd vs even
%                          odd_even_track_corr(this_folder,this_track)= corr(cell2mat(place_fields_odd.track(this_track).raw(good_cells))',...
%                             cell2mat(place_fields_even.track(this_track).raw(good_cells))','type','Pearson');
                    case 'pop_vector'
                        rmap_even= cell2mat(place_fields_even.track(this_track).raw(good_cells)');
                        rmap_odd= cell2mat(place_fields_odd.track(this_track).raw(good_cells)');
                        
                        even_odd_track_mean_corr(this_folder,this_track)= mean(arrayfun(@(x) corr(rmap_even(:,x),rmap_odd(:,x),'type','Pearson'),1:size(rmap_even,2)));
                        even_odd_track_std_corr(this_folder,this_track)= std(arrayfun(@(x) corr(rmap_even(:,x),rmap_odd(:,x),'type','Pearson'),1:size(rmap_even,2)));
                        
%                         odd_even_track_corr(this_folder,this_track)= corr(sum(cell2mat(place_fields_odd.track(this_track).raw(good_cells)'))',...
%                             sum(cell2mat(place_fields_even.track(this_track).raw(good_cells)'))','type','Pearson');
                end
                case 'first_second_half'
                switch p.Results.method
                    case 'ratemaps'
                        first_second_track_corr(this_folder,this_track)= corr(cell2mat(pl_fld_h1.track(this_track).raw(good_cells))',...
                            cell2mat(pl_fld_h2.track(this_track).raw(good_cells))','type','Pearson');
                    case 'pop_vector'
                        rmap_1= cell2mat(pl_fld_h1.track(this_track).raw(good_cells)');
                        rmap_2= cell2mat(pl_fld_h2.track(this_track).raw(good_cells)');
                        
                        first_second_track_mean_corr(this_folder,this_track)= mean(arrayfun(@(x) corr(rmap_1(:,x),rmap_2(:,x),'type','Pearson'),1:size(rmap_1,2)));
                        first_second_track_std_corr(this_folder,this_track)= std(arrayfun(@(x) corr(rmap_1(:,x),rmap_2(:,x),'type','Pearson'),1:size(rmap_1,2)));
                end          
         end
        end
         
end

cd(master_folder);

 switch p.Results.epoch
        case 'even_odd'
            corr_table= table;
            corr_table.even_odd_pair_mean_corr= even_odd_pair_mean_corr;
            corr_table.even_odd_pair_std_corr= even_odd_pair_std_corr;
            corr_table.odd_even_pair_mean_corr= odd_even_pair_mean_corr;
            corr_table.odd_even_pair_std_corr= odd_even_pair_std_corr;
            corr_table.all_laps_pair_mean_corr= pair_mean_corr;
            corr_table.all_laps_pair_std_corr= pair_std_corr;
            corr_table.even_odd_track_mean_corr= even_odd_track_mean_corr;
            corr_table.even_odd_track_std_corr= even_odd_track_std_corr;
            corr_table.odd_even_track_mean_corr= odd_even_track_mean_corr;
            corr_table.odd_even_track_std_corr= odd_even_track_std_corr;
            corr_table.all_laps_track_mean_corr= track_mean_corr;
            corr_table.all_laps_track_std_corr= track_std_corr;
            save('.\Tables\Pearson_corr_tracks_laps_even_odd.mat','corr_table');
            mat= NaN(2,2);
            mat(1,1)= mean(corr_table.even_odd_track_mean_corr(:,1));
            mat(2,2)= mean(corr_table.even_odd_track_mean_corr(:,2));
            mat(1,2)= mean(corr_table.all_laps_pair_mean_corr(:,1));
            mat(2,1)= NaN;
     case 'first_second_half'
         corr_table= table;
            corr_table.first_second_pair_mean_corr= first_second_pair_mean_corr;
            corr_table.first_second_pair_std_corr= first_second_pair_std_corr;
            corr_table.second_first_pair_mean_corr= second_first_pair_mean_corr;
            corr_table.second_first_pair_std_corr= second_first_pair_std_corr;
            corr_table.all_laps_pair_mean_corr= pair_mean_corr;
            corr_table.all_laps_pair_std_corr= pair_std_corr;
            corr_table.first_second_track_mean_corr= first_second_track_mean_corr;
            corr_table.first_second_track_std_corr= first_second_track_std_corr;
%             corr_table.all_laps_track_mean_corr= track_mean_corr;
%             corr_table.all_laps_track_std_corr= track_std_corr;
            save('.\Tables\Pearson_corr_tracks_laps_first_second_half.mat','corr_table');
            mat= NaN(2,2,1);
            mat(1,1,1)= mean(corr_table.first_second_track_mean_corr(:,1));
            mat(2,2,1)= mean(corr_table.first_second_track_mean_corr(:,2));
            mat(1,2,1)= mean(corr_table.all_laps_pair_mean_corr(:,1));
            mat(2,1,1)= NaN;
            mat(1,1,2)= std(corr_table.first_second_track_mean_corr(:,1));
            mat(2,2,2)= std(corr_table.first_second_track_mean_corr(:,2));
            mat(1,2,2)= std(corr_table.all_laps_pair_mean_corr(:,1));
            mat(2,1,2)= NaN;
 end
 
 
% nan_C= mat; nan_C(~isnan(nan_C))=1; nan_C(isnan(nan_C))=0;
figure('Color','w');
imagesc(mat);
colormap(gray); caxis([0 1]);
textStrings= num2str(mat(:),2);
textStrings =strtrim(cellstr(textStrings));  % Remove any space padding
nan_idx= cellfun(@(x) strcmp(x,'NaN'),textStrings);
textStrings(nan_idx)={ ''};
textStrings(strcmp(textStrings,'0'))= {'1'};
[x, y] = meshgrid(1:size(mat,2));  % Create x and y coordinates for the p_val_textstrings
hStrings = text(x(:), y(:), textStrings(:), ...  % Plot the p_val_testrings
                'HorizontalAlignment', 'center');
midValue = mean(get(gca, 'CLim'));  % Get the middle value of the color range
textColors = repmat(mat(:) < midValue, 1, 3);  % Choose white or black for the
                                               %   text color of the p_val_textstrings so
                                               %   they can be easily seen over
                                               %   the background color
set(hStrings, {'Color'}, num2cell(textColors, 2));
xticks([1:2]);yticks([1:2]); set(gca,'TickDir','out','Box','off','YColor','none','XColor','none');
 switch p.Results.epoch
        case 'even_odd'
            title('Pearson Correlation Even/Odd Laps')
        case 'first_second_half'
            title('Pearson Correlation first 50% / last 50% Laps')
 end

end