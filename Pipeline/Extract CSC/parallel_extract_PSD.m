function all_PSD=parallel_extract_PSD(mode)
%performs extract_PSD over multiple channels in parallel
%total number of channels per core is equal to 
%ceil(available channels/number of cores)
%channel numbers greater than the total number of channels (created by using the ceiling function)
%are given a value of NaN and not processed

available_channels=get_available_channels;
p = gcp; % Starting new parallel pool
if isempty(p)
    num_cores = 0;
    disp('parallel processing not possible');
else
    num_cores = p.NumWorkers; 
end
number_of_loops=num_cores
added_empty_channels=number_of_loops*ceil(length(available_channels)/number_of_loops)-length(available_channels);
available_channels=[available_channels NaN(1,added_empty_channels)];
channel_matrix=reshape(available_channels,length(available_channels)/number_of_loops,number_of_loops);

%remove cores processing only NaN channels (caused by an uneven split of channels between existing cores)
while length(channel_matrix(:,number_of_loops))==length(find(isnan(channel_matrix(:,number_of_loops))))  %check if last column is all NaNs
       channel_matrix(:,number_of_loops)=[];
       number_of_loops=number_of_loops-1;
end   
parfor i=1:number_of_loops
   all_PSD_subset{i}= extract_PSD(channel_matrix(:,i),mode);
end
[states,state_name]=get_awake_and_sleep_times;
PSD_states.states = states;
PSD_states.states_name = state_name;

all_PSD=all_PSD_subset{1};

for i=2:number_of_loops
    all_PSD=[all_PSD all_PSD_subset{i}];
end
save('PSD_data','all_PSD','PSD_states');

end


function available_channels = get_available_channels
if exist('raw')==7
    cd raw
    LFP_files=dir('*.ncs');
    dupl_LFP_files= dir('CSC*_*.ncs');
    cd ..
else
    LFP_files=dir('*.ncs');
    dupl_LFP_files= dir('CSC*_*.ncs');
end

% check if there are multiple .ncs files  
if ~isempty(dupl_LFP_files)
    CSC1_i= find(contains({LFP_files.name},'CSC1.ncs'));
    CSC_dupl_i= find(contains({dupl_LFP_files.name},'CSC1_'));
    multiples_CSC= [LFP_files(CSC1_i), dupl_LFP_files(CSC_dupl_i)];
    [~,largest_file]= max([multiples_CSC.bytes]);
    extension= strsplit(multiples_CSC(largest_file).name,'_');
    disp(['selecting _'  extension{2} ' files...']);
    LFP_files= LFP_files([LFP_files.bytes] == max(multiples_CSC.bytes));
end

for c = 1:length(LFP_files)
    LFP_filename = LFP_files(c).name;
    if isempty(dupl_LFP_files)
        name_split = strsplit(LFP_filename,{'.','CSC'},'CollapseDelimiters',true);
        CSC_number = str2num(name_split{2});
        ncs_name= [{'CSC'} {'.ncs'}];
    elseif ~isempty(dupl_LFP_files)
        name_split = strsplit(LFP_filename,{'_','CSC'});
        CSC_number = str2num(name_split{2});
        ncs_name= [{'CSC'} {['_' name_split{3}]}];
    else
        disp('something went wrong in the filenames');
        keyboard;
    end

    available_channels(c) = CSC_number;
end

available_channels= unique(available_channels);
end

