function place_fields_BAYESIAN_combined = generate_cross_experiment_cell_id_shuffles(folders)

workingDir = pwd;
cross_experiment_good_cell_id = [];
place_fields_BAYESIAN_combined = [];
place_fields_BAYESIAN_combined.good_place_cells = [];

for track = 1:2
    place_fields_BAYESIAN_combined.track(track).raw = [];
    place_fields_BAYESIAN_combined.track(track).centre = [];
    place_fields_BAYESIAN_combined.track(track).sorted_good_cells = [];
    place_fields_BAYESIAN_combined.track(track).session_id = [];
end

for f = 1:length(folders)
    cd(folders{f})
    load('extracted_place_fields_BAYESIAN.mat')
    place_fields_BAYESIAN_combined.good_place_cells = [place_fields_BAYESIAN_combined.good_place_cells (1000*f + place_fields_BAYESIAN.good_place_cells)];

    for track = 1:2
        place_fields_BAYESIAN_combined.track(track).raw = [place_fields_BAYESIAN_combined.track(track).raw ...
            {place_fields_BAYESIAN.track(track).raw{place_fields_BAYESIAN.track(track).sorted_good_cells}}];
        place_fields_BAYESIAN_combined.track(track).centre = [place_fields_BAYESIAN_combined.track(track).centre ...
            place_fields_BAYESIAN.track(track).centre(place_fields_BAYESIAN.track(track).sorted_good_cells)];

        place_fields_BAYESIAN_combined.track(track).sorted_good_cells = [place_fields_BAYESIAN_combined.track(track).sorted_good_cells ...
            (1000*f + place_fields_BAYESIAN.track(track).sorted_good_cells)];
        place_fields_BAYESIAN_combined.track(track).session_id = [place_fields_BAYESIAN_combined.track(track).session_id ...
            f*ones(1,length(place_fields_BAYESIAN.track(track).sorted_good_cells))];
    end

    %     place_fields_BAYESIAN_combined.track(1).raw = [place_fields_BAYESIAN_combined.track(1).raw {place_fields_BAYESIAN.track(1).raw{place_fields_BAYESIAN.good_place_cells}}];
    %     place_fields_BAYESIAN_combined.track(2).raw = [place_fields_BAYESIAN_combined.track(2).raw {place_fields_BAYESIAN.track(2).raw{place_fields_BAYESIAN.good_place_cells}}];
    cd(workingDir)
end

for track = 1:2
    % Randomise the cell id first to avoid cells from the same session to
    % cluster together
    index = randperm(length(place_fields_BAYESIAN_combined.track(track).centre));
    place_fields_BAYESIAN_combined.track(track).sorted_good_cells = place_fields_BAYESIAN_combined.track(track).sorted_good_cells(index);
    place_fields_BAYESIAN_combined.track(track).centre = place_fields_BAYESIAN_combined.track(track).centre(index);
    place_fields_BAYESIAN_combined.track(track).session_id = place_fields_BAYESIAN_combined.track(track).session_id(index);
    place_fields_BAYESIAN_combined.track(track).raw = {place_fields_BAYESIAN_combined.track(track).raw{index}};

    % sort according to place field peak centre
    [~,index] = sort(place_fields_BAYESIAN_combined.track(track).centre);
    place_fields_BAYESIAN_combined.track(track).sorted_good_cells = place_fields_BAYESIAN_combined.track(track).sorted_good_cells(index);
    place_fields_BAYESIAN_combined.track(track).centre = place_fields_BAYESIAN_combined.track(track).centre(index);
    place_fields_BAYESIAN_combined.track(track).session_id = place_fields_BAYESIAN_combined.track(track).session_id(index);
    place_fields_BAYESIAN_combined.track(track).raw = {place_fields_BAYESIAN_combined.track(track).raw{index}};
end

save place_fields_BAYESIAN_combined place_fields_BAYESIAN_combined
end



 