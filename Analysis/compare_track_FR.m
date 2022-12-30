function compare_track_FR
% Quick comparison between peak FR distribution of place fields in the
% tracks. Compares within session, merging all sessions for both real and
% shuffled data, respectively. Also compares distributions between real and shuffle data.
% Runs Two-sample Kolmogorov-Smirnov test.


load('X:\BendorLab\Drobo\Neural and Behavioural Data\Rate remapping\Data\rate_remapping_analysis_TRACK_PAIRS_wcorr.mat')
rate_shuffle = load('X:\BendorLab\Drobo\Neural and Behavioural Data\Rate remapping\Data\CONTROLS\rate_remapped\rate_remapping_analysis_TRACK_PAIRS_wcorr_INTRINSIC_RATE.mat');

exp = unique(remapping_raw(6).experiment); %Experiments ID
 
f1= figure('units','normalized','Color','w');
f1.Name = 'Tracks peak FR distributions';
subplot(floor(length(exp)/2)+1,2,1)
histogram(remapping_raw(6).raw_peak_BAYESIAN_plfield_1)
hold on
histogram(remapping_raw(6).raw_peak_BAYESIAN_plfield_2)
[~,p] = kstest2(remapping_raw(5).raw_peak_BAYESIAN_plfield_1,remapping_raw(5).raw_peak_BAYESIAN_plfield_2);
xlabel('Peak FR')
title(['ALL SESSIONS - KS test p= ' num2str(p)])
legend({'Track 1', 'Track 2'})

for ses = 1 : length(exp)
    
    T1_peakFR = remapping_raw(6).raw_peak_BAYESIAN_plfield_1(remapping_raw(6).experiment == ses);
    T2_peakFR = remapping_raw(6).raw_peak_BAYESIAN_plfield_2(remapping_raw(6).experiment == ses);
    
    figure(f1)
    subplot(floor(length(exp)/2)+1,2,ses+1)
    histogram(T1_peakFR) 
    hold on
    histogram(T2_peakFR) 
    [~,p] = kstest2(T1_peakFR,T2_peakFR);
    xlabel('Peak FR')
    title([remapping_raw(6).folder{1,ses} '- KS test p= ' num2str(p)])
    legend({'Track 1', 'Track 2'})
    
end
    
% Compares distribution from real and rate shuffle data
f2= figure('units','normalized','Color','w');
f2.Name = 'Peak FR distributions - real vs rate shuffle';

real_dist = [remapping_raw(6).raw_peak_BAYESIAN_plfield_1 remapping_raw(6).raw_peak_BAYESIAN_plfield_2];
shuffle_dist = [rate_shuffle.remapping_raw(6).raw_peak_BAYESIAN_plfield_1 rate_shuffle.remapping_raw(6).raw_peak_BAYESIAN_plfield_2];
histogram(real_dist)
hold on
histogram(shuffle_dist)
[~,p] = kstest2(real_dist,shuffle_dist);
xlabel('Peak FR')
title(['Real vs Rate Shuffle - KS test p= ' num2str(p)])
legend({'Real', 'Rate shuffle'})



% Compares peak FR in track from shuffled data
f3= figure('units','normalized','Color','w');
f3.Name = 'Rate Shuffle - Tracks peak FR distributions';
subplot(floor(length(exp)/2)+1,2,1)
histogram(rate_shuffle.remapping_raw(6).raw_peak_BAYESIAN_plfield_1)
hold on
histogram(rate_shuffle.remapping_raw(6).raw_peak_BAYESIAN_plfield_2)
[~,p] = kstest2(rate_shuffle.remapping_raw(5).raw_peak_BAYESIAN_plfield_1,rate_shuffle.remapping_raw(5).raw_peak_BAYESIAN_plfield_2);
xlabel('Peak FR')
title(['ALL SESSIONS - KS test p= ' num2str(p)])
legend({'Track 1', 'Track 2'})

for ses = 1 : length(exp)
    
    T1_peakFR = rate_shuffle.remapping_raw(6).raw_peak_BAYESIAN_plfield_1(rate_shuffle.remapping_raw(6).experiment == ses);
    T2_peakFR = rate_shuffle.remapping_raw(6).raw_peak_BAYESIAN_plfield_2(rate_shuffle.remapping_raw(6).experiment == ses);
    
    figure(f3)
    subplot(floor(length(exp)/2)+1,2,ses+1)
    histogram(T1_peakFR) 
    hold on
    histogram(T2_peakFR) 
    [~,p] = kstest2(T1_peakFR,T2_peakFR);
    xlabel('Peak FR')
    title([remapping_raw(6).folder{1,ses} '- KS test p= ' num2str(p)])
    legend({'Track 1', 'Track 2'})
    
end








end