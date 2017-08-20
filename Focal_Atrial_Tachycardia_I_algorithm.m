I = [];
I_abs = abs(I);
expected_p_wave_amplitude = 0.25;
%Determine this value
fraction_of_p_wave_amplitude = 0.3;
%Determine this value. The electrical activity of the heart results from
%the depolarization of the atria, depolarization of the ventricles (which
%masks the repolarization of the atria, and the repolarization of the
%ventricles. A patient who experiences atrial fibrillation will experience
%a quiver and there will be no discernable p-wave. we must determine the
%fraction of the atria that is unlikely to contract indpendently and use
%this as the threshold to determine if the signal is too noisy.
threshold_fraction_of_signal = 0.2;
%Determine this value
counter = 0;
%Determine this value
for i = 1:length(I)
    if I_abs(i) > fraction_of_p_wave_amplitude*expected_p_wave_amplitude;
        counter = counter + 1;
    end
end
if counter/length(I) > threshold_fraction_of_signal
    isClean = 0;
else
    isClean = 1;
end
if isClean
%Sorry about the lack of indentation below. The is a lot of code.
R_wave_threshold = 0.6;
R_wave_values = zeros(length(I),1);
R_wave_indices = zeros(length(I),1);
T_wave_threshold_abs_value = 0.4;
%Confirm this value 
T_wave_values = zeros(length(I),1);
T_wave_indices = zeros(length(I),1);
I_flat_ekg_array = [];
%Construct this array with clean signals.
qrs_complex_duration = 0.1;
T_wave_duration = 0.25;
p_wave_duration = 0.12;
sampling_frequency = 1;
%Determine these values above
qrs_complex_index_length = sampling_frequency*qrs_complex_duration;
T_wave_index_length = sampling_frequency*T_wave_duration;
p_wave_index_length = sampling_frequency*p_wave_duration;
p_wave_boundaries = zeros(length(I),1);
p_wave_boundaries_clone = zeros(length(I),1);
I_diff = zeros(length(I),1);
I_diff2 = zeros(length(I),1);
I_der = diff(I);
I_der2 = diff(I,2);
zero_threshold = 0.005;
discernable_zeros_threshold_time = 0.004;
discernable_zeros_threshold_index = floor(discernable_zeros_threshold_time*sampling_frequency);
for i = 1:length(I)
    if I(i) > R_wave_threshold
        R_wave_values(i) = I(i);
    end
end
R_wave_values_clone = R_wave_values;
for i = 1:length(I)
    if i < length(I) && i > 1
        if ~((R_wave_values(i+1) < R_wave_values(i)) && (R_wave_values(i-1) < R_wave_values(i)))
            R_wave_values_clone(i) = 0;
        end
    else
        R_wave_values_clone(i) = 0;
    end
end
R_wave_values = R_wave_values_clone;
for i = 1:length(I)
    if R_wave_values(i) > 0
        R_wave_indices(i) = i;
        T_wave_values(i) = 100;
    end
end
T_wave_values_clone = T_wave_values;
for i = 1:length(T_wave_values)
    if T_wave_values(i) == 100
        for j = i-(floor(qrs_index_length/2)):i
            T_wave_values_clone(j) = 100;
        end
        %This for loop above makes the values of the T-wave values array
        %before the R-wave peak and within the qrs-complex equal to 100
        for j = i:i+(floor(qrs_index_length/2))
            T_wave_values_clone(j) = 100;
        end
        %This for loop above makes the values of the T-wave values array
        %after the R-wave peak and within the qrs-complex equal to 100
    end
end
T_wave_values = T_wave_values_clone;
for i = 1:length(T_wave_values)
    if (abs(I(i)) > T_wave_threshold_abs_value) && ~(T_wave_values(i) == 100)
        T_wave_values_clone(i) = I(i);
    end
end
T_wave_values = T_wave_values_clone;
for i = 1:length(T_wave_values)
    if T_wave_values(i) == 100
        T_wave_values_clone(i) = 0;
    end
end
T_wave_values = T_wave_values_clone;
for i = 1:length(T_wave_values)
    if i > 1 && i < length(T_wave_values)
        if ~(T_wave_values(i+1) < T_wave_values(i)) && (T_wave_values(i-1) < T_wave_values(i)) || ~((T_wave_values(i+1) > T_wave_values(i)) && (T_wave_values(i-1) > T_wave_values(i))) 
            T_wave_values_clone(i) = 0;
        end
    else
        T_wave_values_clone(i) = 0;
    end
end
T_wave_values = T_wave_values_clone;
for i = 1:length(T_wave_values)
    if abs(T_wave_values(i)) > 0 
        T_wave_indices(i) = i;
    end
end
R_wave_values = R_wave_values(R_wave_values~=0);
T_wave_values = T_wave_values(T_wave_values~=0);
R_wave_indices = R_wave_indices(R_wave_indices~=0);
T_wave_indices = T_wave_indices(T_wave_indices~=0);
R_wave_index_differences = diff(R_wave_indices);
RR_intervals = R_wave_index_differences/sampling_frequency;
Average_RR_interval = sum(RR_intervals)/length(RR_intervals);
heart_rate = 60/Average_RR_interval;
RR_dev_squared = zeros(length(RR_intervals),1);
for i = 1:length(RR_intervals)
    RR_dev_squared(i) = (RR_intervals(i)-Average_RR_interval)^2;
end
RR_interval_SD = sqrt(sum(RR_dev_squared)/(length(RR_dev_squared)-1));
theshold_RR_deviation_avg_ratio = 0.3;
%Confirm this above value
if RR_interval_SD/Average_RR_interval > theshold_RR_deviation_avg_ratio
    disp('Irregular Ventricular Contraction')
else
    disp('Regular Ventricular Contraction')
end
if heart_rate > 140 && heart_rate < 220
    disp('Heart rate is elevated and in the range that is consistent with Atrial Tachycardia.')
else
    disp('Heart is not in the range that is consistent with Atrial Tachycardia.')
end
if heart_rate > 140 && heart_rate < 220 && RR_interval_SD/Average_RR_interval > theshold_RR_deviation_avg_ratio
    disp('Patient is likely afflicted with Atrial Tachycardia')
else
    disp('Heart is likely not afflicted with Atrial Tachycardia')
end
for i = 1:length(I)
    if i == length(I)
        I_diff(i) = I(i);
    else
        I_diff(i) = I_der(i);
    end
end
for i = 1:length(I)
    if i == length(I) || i == length(I) - 1
        I_diff2(i) = I_diff(i);
    else 
        I_diff2(i) = I_der2(i);
    end
end
I_flat_ekg_array_avg = sum(I_flat_ekg_array)/length(I_flat_ekg_array);
I_flat_ekg_array_sq_dev = zeros(length(I_flat_ekg_array),1);
for i = 1:length(I_flat_ekg_array_sq_dev)
    I_flat_ekg_array_sq_dev(i) = (I_flat_ekg_array(i) - I_flat_ekg_array_avg)^2;
end
I_flat_ekg_array_SD = sqrt(sum(I_flat_ekg_array_sq_dev)/(length(I_flat_ekg_array_sq_dev)-1));
I_flat_ekg_array_sorted = sort(I_flat_ekg_array);
I_flat_ekg_array_dev_frm_avg_min_1SD = zeros(length(I_flat_ekg_array),1);
I_flat_ekg_array_dev_frm_avg_pls_1SD = zeros(length(I_flat_ekg_array),1);
I_flat_ekg_array_dev_frm_avg_min_2SD = zeros(length(I_flat_ekg_array),1);
I_flat_ekg_array_dev_frm_avg_pls_2SD = zeros(length(I_flat_ekg_array),1);
I_flat_ekg_array_dev_frm_avg_min_3SD = zeros(length(I_flat_ekg_array),1);
I_flat_ekg_array_dev_frm_avg_pls_3SD = zeros(length(I_flat_ekg_array),1);
for i = 1:length(I_flat_ekg_array)
    I_flat_ekg_array_dev_frm_avg_min_1SD(i) = abs(I_flat_ekg_array_sorted(i) - (I_flat_ekg_array_avg - I_flat_ekg_array_SD));
    I_flat_ekg_array_dev_frm_avg_pls_1SD(i) = abs(I_flat_ekg_array_sorted(i) - (I_flat_ekg_array_avg + I_flat_ekg_array_SD));
    I_flat_ekg_array_dev_frm_avg_min_2SD(i) = abs(I_flat_ekg_array_sorted(i) - (I_flat_ekg_array_avg - 2*I_flat_ekg_array_SD));
    I_flat_ekg_array_dev_frm_avg_pls_2SD(i) = abs(I_flat_ekg_array_sorted(i) - (I_flat_ekg_array_avg + 2*I_flat_ekg_array_SD));
    I_flat_ekg_array_dev_frm_avg_min_3SD(i) = abs(I_flat_ekg_array_sorted(i) -(I_flat_ekg_array_avg - 3*I_flat_ekg_array_SD));
    I_flat_ekg_array_dev_frm_avg_pls_3SD(i) = abs(I_flat_ekg_array_sorted(i) -(I_flat_ekg_array_avg + 3*I_flat_ekg_array_SD));
end
for i = 1:length(I_flat_ekg_array)
    if I_flat_ekg_array_dev_frm_avg_min_1SD(i) == min(I_flat_ekg_array_dev_frm_avg_min_1SD)
        avg_min_1SD_index = i;
    end
    if I_flat_ekg_array_dev_frm_avg_pls_1SD(i) == min(I_flat_ekg_array_dev_frm_avg_pls_1SD)
        avg_pls_1SD_index = i;
    end
    if I_flat_ekg_array_dev_frm_avg_min_2SD(i) == min(I_flat_ekg_array_dev_frm_avg_min_2SD)
        avg_min_2SD_index = i;
    end
    if I_flat_ekg_array_dev_frm_avg_pls_2SD(i) == min(I_flat_ekg_array_dev_frm_avg_pls_2SD)
        avg_pls_2SD_index = i;
    end
    if I_flat_ekg_array_dev_frm_avg_min_3SD(i) == min(I_flat_ekg_array_dev_frm_avg_min_3SD)
        avg_min_3SD_index = i;
    end
    if I_flat_ekg_array_dev_frm_avg_pls_3SD(i) == min(I_flat_ekg_array_dev_frm_avg_pls_3SD)
        avg_pls_3SD_index = i;
    end
end
quotient_1SD = (avg_pls_1SD_index - avg_min_1SD_index)/length(I_flat_ekg_array_sorted);
quotient_2SD = (avg_pls_2SD_index - avg_min_2SD_index)/length(I_flat_ekg_array_sorted);
quotient_3SD = (avg_pls_3SD_index - avg_min_3SD_index)/length(I_flat_ekg_array_sorted);
if (quotient_1SD > 0.66 && quotient_1SD < 0.70) && (quotient_2SD > 0.93 && quotient_2SD < 0.97) &&(quotient_3SD > 0.97 && quotient_3SD < 1)
    isNormal_flatline = 1;
else
    isNormal_flatline = 0;
end
%Here we have a way to test the normality of the flat line data point
%distribution. 1 is taken as true in Matlab, 0 is taken as false 
%First part of strategy here needs to be to find the indices that are off
%limits. These include the ones that correspond to the T-waves and
%QRS-waves. Set the values for p-wave boundaries equal to 1000 or some
%other absurd integer and then remove those specific values later
if isNormal_flatline
    %Here isNormal_flatline is true
    for i = 1:length(I)
        for j = 1:length(R_wave_indices)
            if (i > R_wave_indices(j) - qrs_complex_index_length/2) && (i < T_wave_indices(j) + T_wave_index_length/2) 
                p_wave_boundaries(i) = -1000;
                %It is very important to make sure that we account for the
                %fact that the p-wave could get 'nicked' by this filtering
                %mechanism. If the p-wave is nicked if a patient experiences
                %heart block, then only its starting
                %index will be recorded. This starting index should be
                %filtered out by seeing how far the index pairs are from
                %each other. If one is isolated, then it should be far
                %greater than the one previous and far greater than the one
                %after. Account for this and the qrs-complex and T_wave indices 
                %should be filtered adequately.
            end
        end
    end
    for i = 1:length(I)
        if ~(p_wave_boundaries(i) == -1000) && abs(I(i)) > I_flat_ekg_array_avg + 2*I_flat_ekg_array_SD
            p_wave_boundaries_clone(i) = i;   
        end
    end
    p_wave_boundaries = p_wave_boundaries_clone;
    for i = 1:length(p_wave_boundaries)
        if p_wave_boundaries(i) == -1000
            p_wave_boundaries_clone(i) = 0;
        end
    end
    p_wave_boundaries_clone = p_wave_boundaries;
    recursion_fraction = 0.1;
    test_period = floor(recursion_fraction*p_wave_duration*sampling_frequency);
    for i = 1:length(p_wave_boundaries)
        if (i < length(p_wave_boundaries) - (test_period + 1)) && (i > test_period + 1)
            counter_before_p_wave_boundary = 0;
            counter_after_p_wave_boundary = 0;
            for j = i-test_period:i
                %Filter out all but the end boundaries
                if p_wave_boundaries(j) > 0
                    counter_before_p_wave_boundary = counter_before_p_wave_boundary + 1;
                end
            end
            for j = i:i+test_period
                %Filter out all but the end boundaries
                if p_wave_boundaries(j) > 0
                    counter_after_p_wave_boundary = counter_after_p_wave_boundary + 1;
                end
            end
            if ~((counter_before_p_wave_boundary == 0) || (counter_after_p_wave_boundary == 0))
                p_wave_boundaries_clone(i) = 0;
            end
        end
    end
    p_wave_boundaries = p_wave_boundaries_clone;
    p_wave_boundaries = p_wave_boundaries(p_wave_boundaries~=0);
    p_wave_boundaries_clone = p_wave_boundaries;
    for i = 1:length(p_wave_boundaries)
        if (i > 1) && (i < length(p_wave_boundaries))
            if (p_wave_boundaries(i) > p_wave_boundaries(i-1) + 3*p_wave_index_length) || (p_wave_boundaries(i) < p_wave_boundaries(i+1) - 3*p_wave_index_length)
                p_wave_boundaries_clone(i) = 0;
                %This filters out the boundaries of 'nicked' p-waves
            end
        end
    end
    p_wave_boundaries = p_wave_boundaries_clone;
    p_wave_boundaries = p_wave_boundaries(p_wave_boundaries~=0);
    if ~(mod(length(p_wave_boundaries),2) == 0)
        disp('p_wave boundaries array has a lone p-wave boundary. Array length should not be odd')
    end
    p_wave_start_indices = zeros(length(p_wave_boundaries));
    p_wave_end_indices = zeros(length(p_wave_boundaries));
    for i = 1:length(p_wave_boundaries)
        if mod(i,2) == 1
            p_wave_start_indices(i) = p_wave_boundaries(i);
        else
            p_wave_end_indices(i) = p_wave_boundaries(i);
        end
    end
    p_wave_start_indices = p_wave_start_indices(p_wave_start_indices~=0);
    p_wave_end_indices = p_wave_end_indices(p_wave_end_indices~=0);
else
    %Here isNormal_flatline is false
    disp('The distribution for flat line data points is not sufficiently normal. More data points are likely needed.')
end
I_maxima = zeros(p_wave_start_indices,1);
I_minima = zeros(p_wave_start_indices,1);
I_mdslpes = zeros(p_wave_start_indices,1);
I_concavs = zeros(p_wave_start_indices,1);
I_sym_inds = zeros(p_wave_start_indices,1);
I_number_of_zeros_array = zeros(p_wave_start_indices,1);
I_p_wave_areas_1 = zeros(p_wave_start_indices,1);
I_p_wave_areas_2 = zeros(p_wave_start_indices,1);
I_p_wave_areas_3 = zeros(p_wave_start_indices,1);
I_p_wave_areas_4 = zeros(p_wave_start_indices,1);
I_p_wave_areas_5 = zeros(p_wave_start_indices,1);
for i = 1:length(p_wave_start_indices)
    I_maxima(i) = max(I(p_wave_start_indices(i):p_wave_end_indices(i)));
    I_minima(i) = min(I(p_wave_start_indices(i):p_wave_end_indices(i)));
    I_mdslpes(i) = I_diff(floor((p_wave_start_indices(i)+p_wave_end_indices(i))/2));
    I_concavs(i) = I_diff2(floor((p_wave_start_indices(i)+p_wave_end_indices(i))/2));
    %The code below creates the symmetry index array
    p_wave_pt1 = I(p_wave_start_indices(i):floor((p_wave_start_indices(i)+p_wave_end_indices(i))/2));
    p_wave_pt2 = I(floor((p_wave_start_indices(i)+p_wave_end_indices(i)/2)):p_wave_end_indices(i));
    I_sym_inds(i) = (sum(p_wave_pt2) - sum(p_wave_pt1))/sum(I(p_wave_start_indices(i):p_wave_end_indices(i)));
    %The code above creates the symmetry index array
    %The code below creates the number of zeros array
    I_zeros = zeros((p_wave_end_indices(i)-p_wave_start_indices(i)),1);
    for j = p_wave_start_indices(i):p_wave_end_indices(i)
        if abs(I(j)) < zero_threshold
            I_zeros(j) = j;
        end
    end
    %Take a look at the below section again later
    I_zeros_clone = I_zeros;
    for j = 1:length(I_zeros)
        if I_zeros(j) > 0
            for k = j:j+discernable_zeros_threshold_index
                I_zeros(k) = 0;
            end
        end
    end
    %Take a look at the above section again later
    I_zeros = I_zeros(I_zeros~=0);
    I_number_of_zeros_array(i) = length(I_zeros);
    %The code above creates the number of zeros array
    if ~(isEmpty(I_zeros))
        I_p_wave_areas_1(i) = trapz(I(p_wave_start_indices(i):I_zeros(1)));
    else
        I_p_wave_areas_1(i) = trapz(I(p_wave_start_indices(i):p_wave_end_indices(i)));
    end
    if length(I_zeros) > 1
        I_p_wave_areas_2(i) = trapz(I(I_zeros(1):I_zeros(2)));
    else
        if ~(isEmpty(I_zeros))
            I_p_wave_areas_2(i) = trapz(I(I_zeros(1):p_wave_end_indices(i)));
        end
    end
    if length(I_zeros) > 2
        I_p_wave_areas_3(i) = trapz(I(I_zeros(2):I_zeros(3)));
    else
        if length(I_zeros) > 1
            I_p_wave_areas_3(i) = trapz(I(I_zeros(2):p_wave_end_indices(i)));
        end
    end
    if length(I_zeros) > 3
        I_p_wave_areas_4(i) = trapz(I(I_zeros(3):I_zeros(4)));
    else
        if length(I_zeros) > 2
            I_p_wave_areas_4(i) = trapz(I(I_zeros(3):p_wave_end_indices(i)));
        end
    end
    if length(I_zeros) > 4
        I_p_wave_areas_5(i) = trapz(I(I_zeros(4):I_zeros(5)));
    else
        if length(I_zeros) > 3
            I_p_wave_areas_5(i) = trapz(I(I_zeros(4):p_wave_end_indices(i)));
        end
    end
end
%The code above creates the arrays that correspond to a given patient's
%data
%The code below analyzes the patient's data with respect to a number of
%distributions.
z_score_lower_threshold = 2;
I_num_zeros_SA_node = [];
I_maxima_SA_node = [];
I_minima_SA_node = [];
I_sym_inds_SA_node = [];
I_concavs_SA_node = [];
I_mdslpes_SA_node = [];
I_p_wave_areas_1_SA_node = [];
I_p_wave_areas_2_SA_node = [];
I_p_wave_areas_3_SA_node = [];
I_p_wave_areas_4_SA_node = [];
I_p_wave_areas_5_SA_node = [];
I_num_zeros_crista_terminalis = [];
I_maxima_crista_terminalis = [];
I_minima_crista_terminalis = [];
I_sym_inds_crista_terminalis = [];
I_concavs_crista_terminalis = [];
I_mdslpes_crista_terminalis = [];
I_p_wave_areas_1_crista_terminalis = [];
I_p_wave_areas_2_crista_terminalis = [];
I_p_wave_areas_3_crista_terminalis = [];
I_p_wave_areas_4_crista_terminalis = [];
I_p_wave_areas_5_crista_terminalis = [];
I_num_zeros_tricuspid_annulus = [];
I_maxima_tricuspid_annulus = [];
I_minima_tricuspid_annulus = [];
I_sym_inds_tricuspid_annulus = [];
I_concavs_tricuspid_annulus = [];
I_mdslpes_tricuspid_annulus = [];
I_p_wave_areas_1_tricuspid_annulus = [];
I_p_wave_areas_2_tricuspid_annulus = [];
I_p_wave_areas_3_tricuspid_annulus = [];
I_p_wave_areas_4_tricuspid_annulus = [];
I_p_wave_areas_5_tricuspid_annulus = [];
I_num_zeros_coronary_sinus = [];
I_maxima_coronary_sinus = [];
I_minima_coronary_sinus = [];
I_sym_inds_coronary_sinus = [];
I_concavs_coronary_sinus = [];
I_mdslpes_coronary_sinus = [];
I_p_wave_areas_1_coronary_sinus = [];
I_p_wave_areas_2_coronary_sinus = [];
I_p_wave_areas_3_coronary_sinus = [];
I_p_wave_areas_4_coronary_sinus = [];
I_p_wave_areas_5_coronary_sinus = [];
I_num_zeros_ostium = [];
I_maxima_ostium = [];
I_minima_ostium = [];
I_sym_inds_ostium = [];
I_concavs_ostium = [];
I_mdslpes_ostium = [];
I_p_wave_areas_1_ostium = [];
I_p_wave_areas_2_ostium = [];
I_p_wave_areas_3_ostium = [];
I_p_wave_areas_4_ostium = [];
I_p_wave_areas_5_ostium = [];
I_num_zeros_perinodal = [];
I_maxima_perinodal = [];
I_minima_perinodal = [];
I_sym_inds_perinodal = [];
I_concavs_perinodal = [];
I_mdslpes_perinodal = [];
I_p_wave_areas_1_perinodal = [];
I_p_wave_areas_2_perinodal = [];
I_p_wave_areas_3_perinodal = [];
I_p_wave_areas_4_perinodal = [];
I_p_wave_areas_5_perinodal = [];
I_num_zeros_right_septum = [];
I_maxima_right_septum = [];
I_minima_right_septum = [];
I_sym_inds_right_septum = [];
I_concavs_right_septum = [];
I_mdslpes_right_septum = [];
I_p_wave_areas_1_right_septum = [];
I_p_wave_areas_2_right_septum = [];
I_p_wave_areas_3_right_septum = [];
I_p_wave_areas_4_right_septum = [];
I_p_wave_areas_5_right_septum = [];
I_num_zeros_right_atrial_appendage = [];
I_maxima_right_atrial_appendage = [];
I_minima_right_atrial_appendage = [];
I_sym_inds_right_atrial_appendage = [];
I_concavs_right_atrial_appendage = [];
I_mdslpes_right_atrial_appendage = [];
I_p_wave_areas_1_right_atrial_appendage = [];
I_p_wave_areas_2_right_atrial_appendage = [];
I_p_wave_areas_3_right_atrial_appendage = [];
I_p_wave_areas_4_right_atrial_appendage = [];
I_p_wave_areas_5_right_atrial_appendage = [];
I_num_zeros_pulmonary_veins = [];
I_maxima_pulmonary_veins = [];
I_minima_pulmonary_veins = [];
I_sym_inds_pulmonary_veins = [];
I_concavs_pulmonary_veins = [];
I_mdslpes_pulmonary_veins = [];
I_p_wave_areas_1_pulmonary_veins = [];
I_p_wave_areas_2_pulmonary_veins = [];
I_p_wave_areas_3_pulmonary_veins = [];
I_p_wave_areas_4_pulmonary_veins = [];
I_p_wave_areas_5_pulmonary_veins = [];
I_num_zeros_mitral_annulus = [];
I_maxima_mitral_annulus = [];
I_minima_mitral_annulus = [];
I_sym_inds_mitral_annulus = [];
I_concavs_mitral_annulus = [];
I_mdslpes_mitral_annulus = [];
I_p_wave_areas_1_mitral_annulus = [];
I_p_wave_areas_2_mitral_annulus = [];
I_p_wave_areas_3_mitral_annulus = [];
I_p_wave_areas_4_mitral_annulus = [];
I_p_wave_areas_5_mitral_annulus = [];
I_num_zeros_CS_body = [];
I_maxima_CS_body = [];
I_minima_CS_body = [];
I_sym_inds_CS_body = [];
I_concavs_CS_body = [];
I_mdslpes_CS_body = [];
I_p_wave_areas_1_CS_body = [];
I_p_wave_areas_2_CS_body = [];
I_p_wave_areas_3_CS_body = [];
I_p_wave_areas_4_CS_body = [];
I_p_wave_areas_5_CS_body = [];
I_num_zeros_left_septum = [];
I_maxima_left_septum = [];
I_minima_left_septum = [];
I_sym_inds_left_septum = [];
I_concavs_left_septum = [];
I_mdslpes_left_septum = [];
I_p_wave_areas_1_left_septum = [];
I_p_wave_areas_2_left_septum = [];
I_p_wave_areas_3_left_septum = [];
I_p_wave_areas_4_left_septum = [];
I_p_wave_areas_5_left_septum = [];
I_num_zeros_left_atrial_appendage = [];
I_maxima_left_atrial_appendage = [];
I_minima_left_atrial_appendage = [];
I_sym_inds_left_atrial_appendage = [];
I_concavs_left_atrial_appendage = [];
I_mdslpes_left_atrial_appendage = [];
I_p_wave_areas_1_left_atrial_appendage = [];
I_p_wave_areas_2_left_atrial_appendage = [];
I_p_wave_areas_3_left_atrial_appendage = [];
I_p_wave_areas_4_left_atrial_appendage = [];
I_p_wave_areas_5_left_atrial_appendage = [];
I_num_zeros_z_scores_SA_node = zeros(length(I_number_of_zeros_array),1);
I_num_zeros_SA_node_average = sum(I_num_zeros_SA_node)/length(I_num_zeros_SA_node);
I_num_zeros_SA_node_sq_dev = zeros(length(I_num_zeros_SA_node),1);
for i = 1:length(I_num_zeros_SA_node_sq_dev)
    I_num_zeros_SA_node_sq_dev(i) = (I_num_zeros_SA_node(i)-I_num_zeros_SA_node_average)^2;
end
I_num_zeros_SA_node_SD = sqrt(sum(I_num_zeros_SA_node_sq_dev)/(length(I_num_zeros_SA_node)-1));
for i = 1:length(I_num_zeros_z_scores_SA_node)
    I_num_zeros_z_scores_SA_node(i) = (I_number_of_zeros_array(i) - I_num_zeros_SA_node_average)/I_num_zeros_SA_node_SD;
end
I_maxima_z_scores_SA_node = zeros(length(I_maxima),1);
I_maxima_SA_node_average = sum(I_maxima_SA_node)/length(I_maxima_SA_node);
I_maxima_SA_node_sq_dev = zeros(length(I_maxima_SA_node),1);
for i = 1:length(I_maxima_SA_node_sq_dev)
    I_maxima_SA_node_sq_dev(i) = (I_maxima_SA_node(i)-I_maxima_SA_node_average)^2;
end
I_maxima_SA_node_SD = sqrt(sum(I_maxima_SA_node_sq_dev)/(length(I_maxima_SA_node)-1));
for i = 1:length(I_maxima_z_scores_SA_node)
    I_maxima_z_scores_SA_node(i) = (I_maxima(i) - I_maxima_SA_node_average)/I_maxima_SA_node_SD;
end
I_minima_z_scores_SA_node = zeros(length(I_minima),1);
I_minima_SA_node_average = sum(I_minima_SA_node)/length(I_minima_SA_node);
I_minima_SA_node_sq_dev = zeros(length(I_minima_SA_node),1);
for i = 1:length(I_minima_SA_node_sq_dev)
    I_minima_SA_node_sq_dev(i) = (I_minima_SA_node(i)-I_minima_SA_node_average)^2;
end
I_minima_SA_node_SD = sqrt(sum(I_minima_SA_node_sq_dev)/(length(I_minima_SA_node)-1));
for i = 1:length(I_minima_z_scores_SA_node)
    I_minima_z_scores_SA_node(i) = (I_minima(i) - I_minima_SA_node_average)/I_minima_SA_node_SD;
end
I_mdslpes_z_scores_SA_node = zeros(length(I_mdslpes),1);
I_mdslpes_SA_node_average = sum(I_mdslpes_SA_node)/length(I_mdslpes_SA_node);
I_mdslpes_SA_node_sq_dev = zeros(length(I_mdslpes_SA_node),1);
for i = 1:length(I_mdslpes_SA_node_sq_dev)
    I_mdslpes_SA_node_sq_dev(i) = (I_mdslpes_SA_node(i)-I_mdslpes_SA_node_average)^2;
end
I_mdslpes_SA_node_SD = sqrt(sum(I_mdslpes_SA_node_sq_dev)/(length(I_mdslpes_SA_node)-1));
for i = 1:length(I_mdslpes_z_scores_SA_node)
    I_mdslpes_z_scores_SA_node(i) = (I_mdslpes(i) - I_mdslpes_SA_node_average)/I_mdslpes_SA_node_SD;
end
I_concavs_z_scores_SA_node = zeros(length(I_concavs),1);
I_concavs_SA_node_average = sum(I_concavs_SA_node)/length(I_concavs_SA_node);
I_concavs_SA_node_sq_dev = zeros(length(I_concavs_SA_node),1);
for i = 1:length(I_concavs_SA_node_sq_dev)
    I_concavs_SA_node_sq_dev(i) = (I_concavs_SA_node(i)-I_concavs_SA_node_average)^2;
end
I_concavs_SA_node_SD = sqrt(sum(I_concavs_SA_node_sq_dev)/(length(I_concavs_SA_node)-1));
for i = 1:length(I_concavs_z_scores_SA_node)
    I_concavs_z_scores_SA_node(i) = (I_concavs(i) - I_concavs_SA_node_average)/I_concavs_SA_node_SD;
end
I_sym_inds_z_scores_SA_node = zeros(length(I_sym_inds),1);
I_sym_inds_SA_node_average = sum(I_sym_inds_SA_node)/length(I_sym_inds_SA_node);
I_sym_inds_SA_node_sq_dev = zeros(length(I_sym_inds_SA_node),1);
for i = 1:length(I_sym_inds_SA_node_sq_dev)
    I_sym_inds_SA_node_sq_dev(i) = (I_sym_inds_SA_node(i)-I_sym_inds_SA_node_average)^2;
end
I_sym_inds_SA_node_SD = sqrt(sum(I_sym_inds_SA_node_sq_dev)/(length(I_sym_inds_SA_node)-1));
for i = 1:length(I_sym_inds_z_scores_SA_node)
    I_sym_inds_z_scores_SA_node(i) = (I_sym_inds(i) - I_sym_inds_SA_node_average)/I_sym_inds_SA_node_SD;
end
I_p_wave_areas_1_z_scores_SA_node = zeros(length(I_p_wave_areas_1),1);
I_p_wave_areas_1_SA_node_average = sum(I_p_wave_areas_1_SA_node)/length(I_p_wave_areas_1_SA_node);
I_p_wave_areas_1_SA_node_sq_dev = zeros(length(I_p_wave_areas_1_SA_node),1);
for i = 1:length(I_p_wave_areas_1_SA_node_sq_dev)
    I_p_wave_areas_1_SA_node_sq_dev(i) = (I_p_wave_areas_1_SA_node(i)-I_p_wave_areas_1_SA_node_average)^2;
end
I_p_wave_areas_1_SA_node_SD = sqrt(sum(I_p_wave_areas_1_SA_node_sq_dev)/(length(I_p_wave_areas_1_SA_node)-1));
for i = 1:length(I_p_wave_areas_1_z_scores_SA_node)
    I_p_wave_areas_1_z_scores_SA_node(i) = (I_p_wave_areas_1(i) - I_p_wave_areas_1_SA_node_average)/I_p_wave_areas_1_SA_node_SD;
end
I_p_wave_areas_2_z_scores_SA_node = zeros(length(I_p_wave_areas_2),1);
I_p_wave_areas_2_SA_node_average = sum(I_p_wave_areas_2_SA_node)/length(I_p_wave_areas_2_SA_node);
I_p_wave_areas_2_SA_node_sq_dev = zeros(length(I_p_wave_areas_2_SA_node),1);
for i = 1:length(I_p_wave_areas_2_SA_node_sq_dev)
    I_p_wave_areas_2_SA_node_sq_dev(i) = (I_p_wave_areas_2_SA_node(i)-I_p_wave_areas_2_SA_node_average)^2;
end
I_p_wave_areas_2_SA_node_SD = sqrt(sum(I_p_wave_areas_2_SA_node_sq_dev)/(length(I_p_wave_areas_2_SA_node)-1));
for i = 1:length(I_p_wave_areas_2_z_scores_SA_node)
    I_p_wave_areas_2_z_scores_SA_node(i) = (I_p_wave_areas_2(i) - I_p_wave_areas_2_SA_node_average)/I_p_wave_areas_2_SA_node_SD;
end
I_p_wave_areas_3_z_scores_SA_node = zeros(length(I_p_wave_areas_3),1);
I_p_wave_areas_3_SA_node_average = sum(I_p_wave_areas_3_SA_node)/length(I_p_wave_areas_3_SA_node);
I_p_wave_areas_3_SA_node_sq_dev = zeros(length(I_p_wave_areas_3_SA_node),1);
for i = 1:length(I_p_wave_areas_3_SA_node_sq_dev)
    I_p_wave_areas_3_SA_node_sq_dev(i) = (I_p_wave_areas_3_SA_node(i)-I_p_wave_areas_3_SA_node_average)^3;
end
I_p_wave_areas_3_SA_node_SD = sqrt(sum(I_p_wave_areas_3_SA_node_sq_dev)/(length(I_p_wave_areas_3_SA_node)-1));
for i = 1:length(I_p_wave_areas_3_z_scores_SA_node)
    I_p_wave_areas_3_z_scores_SA_node(i) = (I_p_wave_areas_3(i) - I_p_wave_areas_3_SA_node_average)/I_p_wave_areas_3_SA_node_SD;
end
I_p_wave_areas_4_z_scores_SA_node = zeros(length(I_p_wave_areas_4),1);
I_p_wave_areas_4_SA_node_average = sum(I_p_wave_areas_4_SA_node)/length(I_p_wave_areas_4_SA_node);
I_p_wave_areas_4_SA_node_sq_dev = zeros(length(I_p_wave_areas_4_SA_node),1);
for i = 1:length(I_p_wave_areas_4_SA_node_sq_dev)
    I_p_wave_areas_4_SA_node_sq_dev(i) = (I_p_wave_areas_4_SA_node(i)-I_p_wave_areas_4_SA_node_average)^4;
end
I_p_wave_areas_4_SA_node_SD = sqrt(sum(I_p_wave_areas_4_SA_node_sq_dev)/(length(I_p_wave_areas_4_SA_node)-1));
for i = 1:length(I_p_wave_areas_4_z_scores_SA_node)
    I_p_wave_areas_4_z_scores_SA_node(i) = (I_p_wave_areas_4(i) - I_p_wave_areas_4_SA_node_average)/I_p_wave_areas_4_SA_node_SD;
end
I_p_wave_areas_5_z_scores_SA_node = zeros(length(I_p_wave_areas_5),1);
I_p_wave_areas_5_SA_node_average = sum(I_p_wave_areas_5_SA_node)/length(I_p_wave_areas_5_SA_node);
I_p_wave_areas_5_SA_node_sq_dev = zeros(length(I_p_wave_areas_5_SA_node),1);
for i = 1:length(I_p_wave_areas_5_SA_node_sq_dev)
    I_p_wave_areas_5_SA_node_sq_dev(i) = (I_p_wave_areas_5_SA_node(i)-I_p_wave_areas_5_SA_node_average)^5;
end
I_p_wave_areas_5_SA_node_SD = sqrt(sum(I_p_wave_areas_5_SA_node_sq_dev)/(length(I_p_wave_areas_5_SA_node)-1));
for i = 1:length(I_p_wave_areas_5_z_scores_SA_node)
    I_p_wave_areas_5_z_scores_SA_node(i) = (I_p_wave_areas_5(i) - I_p_wave_areas_5_SA_node_average)/I_p_wave_areas_5_SA_node_SD;
end
I_num_zeros_z_scores_crista_terminalis = zeros(length(I_number_of_zeros_array),1);
I_num_zeros_crista_terminalis_average = sum(I_num_zeros_crista_terminalis)/length(I_num_zeros_crista_terminalis);
I_num_zeros_crista_terminalis_sq_dev = zeros(length(I_num_zeros_crista_terminalis),1);
for i = 1:length(I_num_zeros_crista_terminalis_sq_dev)
    I_num_zeros_crista_terminalis_sq_dev(i) = (I_num_zeros_crista_terminalis(i)-I_num_zeros_crista_terminalis_average)^2;
end
I_num_zeros_crista_terminalis_SD = sqrt(sum(I_num_zeros_crista_terminalis_sq_dev)/(length(I_num_zeros_crista_terminalis)-1));
for i = 1:length(I_num_zeros_z_scores_crista_terminalis)
    I_num_zeros_z_scores_crista_terminalis(i) = (I_number_of_zeros_array(i) - I_num_zeros_crista_terminalis_average)/I_num_zeros_crista_terminalis_SD;
end
I_maxima_z_scores_crista_terminalis = zeros(length(I_maxima),1);
I_maxima_crista_terminalis_average = sum(I_maxima_crista_terminalis)/length(I_maxima_crista_terminalis);
I_maxima_crista_terminalis_sq_dev = zeros(length(I_maxima_crista_terminalis),1);
for i = 1:length(I_maxima_crista_terminalis_sq_dev)
    I_maxima_crista_terminalis_sq_dev(i) = (I_maxima_crista_terminalis(i)-I_maxima_crista_terminalis_average)^2;
end
I_maxima_crista_terminalis_SD = sqrt(sum(I_maxima_crista_terminalis_sq_dev)/(length(I_maxima_crista_terminalis)-1));
for i = 1:length(I_maxima_z_scores_crista_terminalis)
    I_maxima_z_scores_crista_terminalis(i) = (I_maxima(i) - I_maxima_crista_terminalis_average)/I_maxima_crista_terminalis_SD;
end
I_minima_z_scores_crista_terminalis = zeros(length(I_minima),1);
I_minima_crista_terminalis_average = sum(I_minima_crista_terminalis)/length(I_minima_crista_terminalis);
I_minima_crista_terminalis_sq_dev = zeros(length(I_minima_crista_terminalis),1);
for i = 1:length(I_minima_crista_terminalis_sq_dev)
    I_minima_crista_terminalis_sq_dev(i) = (I_minima_crista_terminalis(i)-I_minima_crista_terminalis_average)^2;
end
I_minima_crista_terminalis_SD = sqrt(sum(I_minima_crista_terminalis_sq_dev)/(length(I_minima_crista_terminalis)-1));
for i = 1:length(I_minima_z_scores_crista_terminalis)
    I_minima_z_scores_crista_terminalis(i) = (I_minima(i) - I_minima_crista_terminalis_average)/I_minima_crista_terminalis_SD;
end
I_mdslpes_z_scores_crista_terminalis = zeros(length(I_mdslpes),1);
I_mdslpes_crista_terminalis_average = sum(I_mdslpes_crista_terminalis)/length(I_mdslpes_crista_terminalis);
I_mdslpes_crista_terminalis_sq_dev = zeros(length(I_mdslpes_crista_terminalis),1);
for i = 1:length(I_mdslpes_crista_terminalis_sq_dev)
    I_mdslpes_crista_terminalis_sq_dev(i) = (I_mdslpes_crista_terminalis(i)-I_mdslpes_crista_terminalis_average)^2;
end
I_mdslpes_crista_terminalis_SD = sqrt(sum(I_mdslpes_crista_terminalis_sq_dev)/(length(I_mdslpes_crista_terminalis)-1));
for i = 1:length(I_mdslpes_z_scores_crista_terminalis)
    I_mdslpes_z_scores_crista_terminalis(i) = (I_mdslpes(i) - I_mdslpes_crista_terminalis_average)/I_mdslpes_crista_terminalis_SD;
end
I_concavs_z_scores_crista_terminalis = zeros(length(I_concavs),1);
I_concavs_crista_terminalis_average = sum(I_concavs_crista_terminalis)/length(I_concavs_crista_terminalis);
I_concavs_crista_terminalis_sq_dev = zeros(length(I_concavs_crista_terminalis),1);
for i = 1:length(I_concavs_crista_terminalis_sq_dev)
    I_concavs_crista_terminalis_sq_dev(i) = (I_concavs_crista_terminalis(i)-I_concavs_crista_terminalis_average)^2;
end
I_concavs_crista_terminalis_SD = sqrt(sum(I_concavs_crista_terminalis_sq_dev)/(length(I_concavs_crista_terminalis)-1));
for i = 1:length(I_concavs_z_scores_crista_terminalis)
    I_concavs_z_scores_crista_terminalis(i) = (I_concavs(i) - I_concavs_crista_terminalis_average)/I_concavs_crista_terminalis_SD;
end
I_sym_inds_z_scores_crista_terminalis = zeros(length(I_sym_inds),1);
I_sym_inds_crista_terminalis_average = sum(I_sym_inds_crista_terminalis)/length(I_sym_inds_crista_terminalis);
I_sym_inds_crista_terminalis_sq_dev = zeros(length(I_sym_inds_crista_terminalis),1);
for i = 1:length(I_sym_inds_crista_terminalis_sq_dev)
    I_sym_inds_crista_terminalis_sq_dev(i) = (I_sym_inds_crista_terminalis(i)-I_sym_inds_crista_terminalis_average)^2;
end
I_sym_inds_crista_terminalis_SD = sqrt(sum(I_sym_inds_crista_terminalis_sq_dev)/(length(I_sym_inds_crista_terminalis)-1));
for i = 1:length(I_sym_inds_z_scores_crista_terminalis)
    I_sym_inds_z_scores_crista_terminalis(i) = (I_sym_inds(i) - I_sym_inds_crista_terminalis_average)/I_sym_inds_crista_terminalis_SD;
end
I_p_wave_areas_1_z_scores_crista_terminalis = zeros(length(I_p_wave_areas_1),1);
I_p_wave_areas_1_crista_terminalis_average = sum(I_p_wave_areas_1_crista_terminalis)/length(I_p_wave_areas_1_crista_terminalis);
I_p_wave_areas_1_crista_terminalis_sq_dev = zeros(length(I_p_wave_areas_1_crista_terminalis),1);
for i = 1:length(I_p_wave_areas_1_crista_terminalis_sq_dev)
    I_p_wave_areas_1_crista_terminalis_sq_dev(i) = (I_p_wave_areas_1_crista_terminalis(i)-I_p_wave_areas_1_crista_terminalis_average)^2;
end
I_p_wave_areas_1_crista_terminalis_SD = sqrt(sum(I_p_wave_areas_1_crista_terminalis_sq_dev)/(length(I_p_wave_areas_1_crista_terminalis)-1));
for i = 1:length(I_p_wave_areas_1_z_scores_crista_terminalis)
    I_p_wave_areas_1_z_scores_crista_terminalis(i) = (I_p_wave_areas_1(i) - I_p_wave_areas_1_crista_terminalis_average)/I_p_wave_areas_1_crista_terminalis_SD;
end
I_p_wave_areas_2_z_scores_crista_terminalis = zeros(length(I_p_wave_areas_2),1);
I_p_wave_areas_2_crista_terminalis_average = sum(I_p_wave_areas_2_crista_terminalis)/length(I_p_wave_areas_2_crista_terminalis);
I_p_wave_areas_2_crista_terminalis_sq_dev = zeros(length(I_p_wave_areas_2_crista_terminalis),1);
for i = 1:length(I_p_wave_areas_2_crista_terminalis_sq_dev)
    I_p_wave_areas_2_crista_terminalis_sq_dev(i) = (I_p_wave_areas_2_crista_terminalis(i)-I_p_wave_areas_2_crista_terminalis_average)^2;
end
I_p_wave_areas_2_crista_terminalis_SD = sqrt(sum(I_p_wave_areas_2_crista_terminalis_sq_dev)/(length(I_p_wave_areas_2_crista_terminalis)-1));
for i = 1:length(I_p_wave_areas_2_z_scores_crista_terminalis)
    I_p_wave_areas_2_z_scores_crista_terminalis(i) = (I_p_wave_areas_2(i) - I_p_wave_areas_2_crista_terminalis_average)/I_p_wave_areas_2_crista_terminalis_SD;
end
I_p_wave_areas_3_z_scores_crista_terminalis = zeros(length(I_p_wave_areas_3),1);
I_p_wave_areas_3_crista_terminalis_average = sum(I_p_wave_areas_3_crista_terminalis)/length(I_p_wave_areas_3_crista_terminalis);
I_p_wave_areas_3_crista_terminalis_sq_dev = zeros(length(I_p_wave_areas_3_crista_terminalis),1);
for i = 1:length(I_p_wave_areas_3_crista_terminalis_sq_dev)
    I_p_wave_areas_3_crista_terminalis_sq_dev(i) = (I_p_wave_areas_3_crista_terminalis(i)-I_p_wave_areas_3_crista_terminalis_average)^3;
end
I_p_wave_areas_3_crista_terminalis_SD = sqrt(sum(I_p_wave_areas_3_crista_terminalis_sq_dev)/(length(I_p_wave_areas_3_crista_terminalis)-1));
for i = 1:length(I_p_wave_areas_3_z_scores_crista_terminalis)
    I_p_wave_areas_3_z_scores_crista_terminalis(i) = (I_p_wave_areas_3(i) - I_p_wave_areas_3_crista_terminalis_average)/I_p_wave_areas_3_crista_terminalis_SD;
end
I_p_wave_areas_4_z_scores_crista_terminalis = zeros(length(I_p_wave_areas_4),1);
I_p_wave_areas_4_crista_terminalis_average = sum(I_p_wave_areas_4_crista_terminalis)/length(I_p_wave_areas_4_crista_terminalis);
I_p_wave_areas_4_crista_terminalis_sq_dev = zeros(length(I_p_wave_areas_4_crista_terminalis),1);
for i = 1:length(I_p_wave_areas_4_crista_terminalis_sq_dev)
    I_p_wave_areas_4_crista_terminalis_sq_dev(i) = (I_p_wave_areas_4_crista_terminalis(i)-I_p_wave_areas_4_crista_terminalis_average)^4;
end
I_p_wave_areas_4_crista_terminalis_SD = sqrt(sum(I_p_wave_areas_4_crista_terminalis_sq_dev)/(length(I_p_wave_areas_4_crista_terminalis)-1));
for i = 1:length(I_p_wave_areas_4_z_scores_crista_terminalis)
    I_p_wave_areas_4_z_scores_crista_terminalis(i) = (I_p_wave_areas_4(i) - I_p_wave_areas_4_crista_terminalis_average)/I_p_wave_areas_4_crista_terminalis_SD;
end
I_p_wave_areas_5_z_scores_crista_terminalis = zeros(length(I_p_wave_areas_5),1);
I_p_wave_areas_5_crista_terminalis_average = sum(I_p_wave_areas_5_crista_terminalis)/length(I_p_wave_areas_5_crista_terminalis);
I_p_wave_areas_5_crista_terminalis_sq_dev = zeros(length(I_p_wave_areas_5_crista_terminalis),1);
for i = 1:length(I_p_wave_areas_5_crista_terminalis_sq_dev)
    I_p_wave_areas_5_crista_terminalis_sq_dev(i) = (I_p_wave_areas_5_crista_terminalis(i)-I_p_wave_areas_5_crista_terminalis_average)^5;
end
I_p_wave_areas_5_crista_terminalis_SD = sqrt(sum(I_p_wave_areas_5_crista_terminalis_sq_dev)/(length(I_p_wave_areas_5_crista_terminalis)-1));
for i = 1:length(I_p_wave_areas_5_z_scores_crista_terminalis)
    I_p_wave_areas_5_z_scores_crista_terminalis(i) = (I_p_wave_areas_5(i) - I_p_wave_areas_5_crista_terminalis_average)/I_p_wave_areas_5_crista_terminalis_SD;
end
I_num_zeros_z_scores_tricuspid_annulus = zeros(length(I_number_of_zeros_array),1);
I_num_zeros_tricuspid_annulus_average = sum(I_num_zeros_tricuspid_annulus)/length(I_num_zeros_tricuspid_annulus);
I_num_zeros_tricuspid_annulus_sq_dev = zeros(length(I_num_zeros_tricuspid_annulus),1);
for i = 1:length(I_num_zeros_tricuspid_annulus_sq_dev)
    I_num_zeros_tricuspid_annulus_sq_dev(i) = (I_num_zeros_tricuspid_annulus(i)-I_num_zeros_tricuspid_annulus_average)^2;
end
I_num_zeros_tricuspid_annulus_SD = sqrt(sum(I_num_zeros_tricuspid_annulus_sq_dev)/(length(I_num_zeros_tricuspid_annulus)-1));
for i = 1:length(I_num_zeros_z_scores_tricuspid_annulus)
    I_num_zeros_z_scores_tricuspid_annulus(i) = (I_number_of_zeros_array(i) - I_num_zeros_tricuspid_annulus_average)/I_num_zeros_tricuspid_annulus_SD;
end
I_maxima_z_scores_tricuspid_annulus = zeros(length(I_maxima),1);
I_maxima_tricuspid_annulus_average = sum(I_maxima_tricuspid_annulus)/length(I_maxima_tricuspid_annulus);
I_maxima_tricuspid_annulus_sq_dev = zeros(length(I_maxima_tricuspid_annulus),1);
for i = 1:length(I_maxima_tricuspid_annulus_sq_dev)
    I_maxima_tricuspid_annulus_sq_dev(i) = (I_maxima_tricuspid_annulus(i)-I_maxima_tricuspid_annulus_average)^2;
end
I_maxima_tricuspid_annulus_SD = sqrt(sum(I_maxima_tricuspid_annulus_sq_dev)/(length(I_maxima_tricuspid_annulus)-1));
for i = 1:length(I_maxima_z_scores_tricuspid_annulus)
    I_maxima_z_scores_tricuspid_annulus(i) = (I_maxima(i) - I_maxima_tricuspid_annulus_average)/I_maxima_tricuspid_annulus_SD;
end
I_minima_z_scores_tricuspid_annulus = zeros(length(I_minima),1);
I_minima_tricuspid_annulus_average = sum(I_minima_tricuspid_annulus)/length(I_minima_tricuspid_annulus);
I_minima_tricuspid_annulus_sq_dev = zeros(length(I_minima_tricuspid_annulus),1);
for i = 1:length(I_minima_tricuspid_annulus_sq_dev)
    I_minima_tricuspid_annulus_sq_dev(i) = (I_minima_tricuspid_annulus(i)-I_minima_tricuspid_annulus_average)^2;
end
I_minima_tricuspid_annulus_SD = sqrt(sum(I_minima_tricuspid_annulus_sq_dev)/(length(I_minima_tricuspid_annulus)-1));
for i = 1:length(I_minima_z_scores_tricuspid_annulus)
    I_minima_z_scores_tricuspid_annulus(i) = (I_minima(i) - I_minima_tricuspid_annulus_average)/I_minima_tricuspid_annulus_SD;
end
I_mdslpes_z_scores_tricuspid_annulus = zeros(length(I_mdslpes),1);
I_mdslpes_tricuspid_annulus_average = sum(I_mdslpes_tricuspid_annulus)/length(I_mdslpes_tricuspid_annulus);
I_mdslpes_tricuspid_annulus_sq_dev = zeros(length(I_mdslpes_tricuspid_annulus),1);
for i = 1:length(I_mdslpes_tricuspid_annulus_sq_dev)
    I_mdslpes_tricuspid_annulus_sq_dev(i) = (I_mdslpes_tricuspid_annulus(i)-I_mdslpes_tricuspid_annulus_average)^2;
end
I_mdslpes_tricuspid_annulus_SD = sqrt(sum(I_mdslpes_tricuspid_annulus_sq_dev)/(length(I_mdslpes_tricuspid_annulus)-1));
for i = 1:length(I_mdslpes_z_scores_tricuspid_annulus)
    I_mdslpes_z_scores_tricuspid_annulus(i) = (I_mdslpes(i) - I_mdslpes_tricuspid_annulus_average)/I_mdslpes_tricuspid_annulus_SD;
end
I_concavs_z_scores_tricuspid_annulus = zeros(length(I_concavs),1);
I_concavs_tricuspid_annulus_average = sum(I_concavs_tricuspid_annulus)/length(I_concavs_tricuspid_annulus);
I_concavs_tricuspid_annulus_sq_dev = zeros(length(I_concavs_tricuspid_annulus),1);
for i = 1:length(I_concavs_tricuspid_annulus_sq_dev)
    I_concavs_tricuspid_annulus_sq_dev(i) = (I_concavs_tricuspid_annulus(i)-I_concavs_tricuspid_annulus_average)^2;
end
I_concavs_tricuspid_annulus_SD = sqrt(sum(I_concavs_tricuspid_annulus_sq_dev)/(length(I_concavs_tricuspid_annulus)-1));
for i = 1:length(I_concavs_z_scores_tricuspid_annulus)
    I_concavs_z_scores_tricuspid_annulus(i) = (I_concavs(i) - I_concavs_tricuspid_annulus_average)/I_concavs_tricuspid_annulus_SD;
end
I_sym_inds_z_scores_tricuspid_annulus = zeros(length(I_sym_inds),1);
I_sym_inds_tricuspid_annulus_average = sum(I_sym_inds_tricuspid_annulus)/length(I_sym_inds_tricuspid_annulus);
I_sym_inds_tricuspid_annulus_sq_dev = zeros(length(I_sym_inds_tricuspid_annulus),1);
for i = 1:length(I_sym_inds_tricuspid_annulus_sq_dev)
    I_sym_inds_tricuspid_annulus_sq_dev(i) = (I_sym_inds_tricuspid_annulus(i)-I_sym_inds_tricuspid_annulus_average)^2;
end
I_sym_inds_tricuspid_annulus_SD = sqrt(sum(I_sym_inds_tricuspid_annulus_sq_dev)/(length(I_sym_inds_tricuspid_annulus)-1));
for i = 1:length(I_sym_inds_z_scores_tricuspid_annulus)
    I_sym_inds_z_scores_tricuspid_annulus(i) = (I_sym_inds(i) - I_sym_inds_tricuspid_annulus_average)/I_sym_inds_tricuspid_annulus_SD;
end
I_p_wave_areas_1_z_scores_tricuspid_annulus = zeros(length(I_p_wave_areas_1),1);
I_p_wave_areas_1_tricuspid_annulus_average = sum(I_p_wave_areas_1_tricuspid_annulus)/length(I_p_wave_areas_1_tricuspid_annulus);
I_p_wave_areas_1_tricuspid_annulus_sq_dev = zeros(length(I_p_wave_areas_1_tricuspid_annulus),1);
for i = 1:length(I_p_wave_areas_1_tricuspid_annulus_sq_dev)
    I_p_wave_areas_1_tricuspid_annulus_sq_dev(i) = (I_p_wave_areas_1_tricuspid_annulus(i)-I_p_wave_areas_1_tricuspid_annulus_average)^2;
end
I_p_wave_areas_1_tricuspid_annulus_SD = sqrt(sum(I_p_wave_areas_1_tricuspid_annulus_sq_dev)/(length(I_p_wave_areas_1_tricuspid_annulus)-1));
for i = 1:length(I_p_wave_areas_1_z_scores_tricuspid_annulus)
    I_p_wave_areas_1_z_scores_tricuspid_annulus(i) = (I_p_wave_areas_1(i) - I_p_wave_areas_1_tricuspid_annulus_average)/I_p_wave_areas_1_tricuspid_annulus_SD;
end
I_p_wave_areas_2_z_scores_tricuspid_annulus = zeros(length(I_p_wave_areas_2),1);
I_p_wave_areas_2_tricuspid_annulus_average = sum(I_p_wave_areas_2_tricuspid_annulus)/length(I_p_wave_areas_2_tricuspid_annulus);
I_p_wave_areas_2_tricuspid_annulus_sq_dev = zeros(length(I_p_wave_areas_2_tricuspid_annulus),1);
for i = 1:length(I_p_wave_areas_2_tricuspid_annulus_sq_dev)
    I_p_wave_areas_2_tricuspid_annulus_sq_dev(i) = (I_p_wave_areas_2_tricuspid_annulus(i)-I_p_wave_areas_2_tricuspid_annulus_average)^2;
end
I_p_wave_areas_2_tricuspid_annulus_SD = sqrt(sum(I_p_wave_areas_2_tricuspid_annulus_sq_dev)/(length(I_p_wave_areas_2_tricuspid_annulus)-1));
for i = 1:length(I_p_wave_areas_2_z_scores_tricuspid_annulus)
    I_p_wave_areas_2_z_scores_tricuspid_annulus(i) = (I_p_wave_areas_2(i) - I_p_wave_areas_2_tricuspid_annulus_average)/I_p_wave_areas_2_tricuspid_annulus_SD;
end
I_p_wave_areas_3_z_scores_tricuspid_annulus = zeros(length(I_p_wave_areas_3),1);
I_p_wave_areas_3_tricuspid_annulus_average = sum(I_p_wave_areas_3_tricuspid_annulus)/length(I_p_wave_areas_3_tricuspid_annulus);
I_p_wave_areas_3_tricuspid_annulus_sq_dev = zeros(length(I_p_wave_areas_3_tricuspid_annulus),1);
for i = 1:length(I_p_wave_areas_3_tricuspid_annulus_sq_dev)
    I_p_wave_areas_3_tricuspid_annulus_sq_dev(i) = (I_p_wave_areas_3_tricuspid_annulus(i)-I_p_wave_areas_3_tricuspid_annulus_average)^3;
end
I_p_wave_areas_3_tricuspid_annulus_SD = sqrt(sum(I_p_wave_areas_3_tricuspid_annulus_sq_dev)/(length(I_p_wave_areas_3_tricuspid_annulus)-1));
for i = 1:length(I_p_wave_areas_3_z_scores_tricuspid_annulus)
    I_p_wave_areas_3_z_scores_tricuspid_annulus(i) = (I_p_wave_areas_3(i) - I_p_wave_areas_3_tricuspid_annulus_average)/I_p_wave_areas_3_tricuspid_annulus_SD;
end
I_p_wave_areas_4_z_scores_tricuspid_annulus = zeros(length(I_p_wave_areas_4),1);
I_p_wave_areas_4_tricuspid_annulus_average = sum(I_p_wave_areas_4_tricuspid_annulus)/length(I_p_wave_areas_4_tricuspid_annulus);
I_p_wave_areas_4_tricuspid_annulus_sq_dev = zeros(length(I_p_wave_areas_4_tricuspid_annulus),1);
for i = 1:length(I_p_wave_areas_4_tricuspid_annulus_sq_dev)
    I_p_wave_areas_4_tricuspid_annulus_sq_dev(i) = (I_p_wave_areas_4_tricuspid_annulus(i)-I_p_wave_areas_4_tricuspid_annulus_average)^4;
end
I_p_wave_areas_4_tricuspid_annulus_SD = sqrt(sum(I_p_wave_areas_4_tricuspid_annulus_sq_dev)/(length(I_p_wave_areas_4_tricuspid_annulus)-1));
for i = 1:length(I_p_wave_areas_4_z_scores_tricuspid_annulus)
    I_p_wave_areas_4_z_scores_tricuspid_annulus(i) = (I_p_wave_areas_4(i) - I_p_wave_areas_4_tricuspid_annulus_average)/I_p_wave_areas_4_tricuspid_annulus_SD;
end
I_p_wave_areas_5_z_scores_tricuspid_annulus = zeros(length(I_p_wave_areas_5),1);
I_p_wave_areas_5_tricuspid_annulus_average = sum(I_p_wave_areas_5_tricuspid_annulus)/length(I_p_wave_areas_5_tricuspid_annulus);
I_p_wave_areas_5_tricuspid_annulus_sq_dev = zeros(length(I_p_wave_areas_5_tricuspid_annulus),1);
for i = 1:length(I_p_wave_areas_5_tricuspid_annulus_sq_dev)
    I_p_wave_areas_5_tricuspid_annulus_sq_dev(i) = (I_p_wave_areas_5_tricuspid_annulus(i)-I_p_wave_areas_5_tricuspid_annulus_average)^5;
end
I_p_wave_areas_5_tricuspid_annulus_SD = sqrt(sum(I_p_wave_areas_5_tricuspid_annulus_sq_dev)/(length(I_p_wave_areas_5_tricuspid_annulus)-1));
for i = 1:length(I_p_wave_areas_5_z_scores_tricuspid_annulus)
    I_p_wave_areas_5_z_scores_tricuspid_annulus(i) = (I_p_wave_areas_5(i) - I_p_wave_areas_5_tricuspid_annulus_average)/I_p_wave_areas_5_tricuspid_annulus_SD;
end
I_num_zeros_z_scores_coronary_sinus = zeros(length(I_number_of_zeros_array),1);
I_num_zeros_coronary_sinus_average = sum(I_num_zeros_coronary_sinus)/length(I_num_zeros_coronary_sinus);
I_num_zeros_coronary_sinus_sq_dev = zeros(length(I_num_zeros_coronary_sinus),1);
for i = 1:length(I_num_zeros_coronary_sinus_sq_dev)
    I_num_zeros_coronary_sinus_sq_dev(i) = (I_num_zeros_coronary_sinus(i)-I_num_zeros_coronary_sinus_average)^2;
end
I_num_zeros_coronary_sinus_SD = sqrt(sum(I_num_zeros_coronary_sinus_sq_dev)/(length(I_num_zeros_coronary_sinus)-1));
for i = 1:length(I_num_zeros_z_scores_coronary_sinus)
    I_num_zeros_z_scores_coronary_sinus(i) = (I_number_of_zeros_array(i) - I_num_zeros_coronary_sinus_average)/I_num_zeros_coronary_sinus_SD;
end
I_maxima_z_scores_coronary_sinus = zeros(length(I_maxima),1);
I_maxima_coronary_sinus_average = sum(I_maxima_coronary_sinus)/length(I_maxima_coronary_sinus);
I_maxima_coronary_sinus_sq_dev = zeros(length(I_maxima_coronary_sinus),1);
for i = 1:length(I_maxima_coronary_sinus_sq_dev)
    I_maxima_coronary_sinus_sq_dev(i) = (I_maxima_coronary_sinus(i)-I_maxima_coronary_sinus_average)^2;
end
I_maxima_coronary_sinus_SD = sqrt(sum(I_maxima_coronary_sinus_sq_dev)/(length(I_maxima_coronary_sinus)-1));
for i = 1:length(I_maxima_z_scores_coronary_sinus)
    I_maxima_z_scores_coronary_sinus(i) = (I_maxima(i) - I_maxima_coronary_sinus_average)/I_maxima_coronary_sinus_SD;
end
I_minima_z_scores_coronary_sinus = zeros(length(I_minima),1);
I_minima_coronary_sinus_average = sum(I_minima_coronary_sinus)/length(I_minima_coronary_sinus);
I_minima_coronary_sinus_sq_dev = zeros(length(I_minima_coronary_sinus),1);
for i = 1:length(I_minima_coronary_sinus_sq_dev)
    I_minima_coronary_sinus_sq_dev(i) = (I_minima_coronary_sinus(i)-I_minima_coronary_sinus_average)^2;
end
I_minima_coronary_sinus_SD = sqrt(sum(I_minima_coronary_sinus_sq_dev)/(length(I_minima_coronary_sinus)-1));
for i = 1:length(I_minima_z_scores_coronary_sinus)
    I_minima_z_scores_coronary_sinus(i) = (I_minima(i) - I_minima_coronary_sinus_average)/I_minima_coronary_sinus_SD;
end
I_mdslpes_z_scores_coronary_sinus = zeros(length(I_mdslpes),1);
I_mdslpes_coronary_sinus_average = sum(I_mdslpes_coronary_sinus)/length(I_mdslpes_coronary_sinus);
I_mdslpes_coronary_sinus_sq_dev = zeros(length(I_mdslpes_coronary_sinus),1);
for i = 1:length(I_mdslpes_coronary_sinus_sq_dev)
    I_mdslpes_coronary_sinus_sq_dev(i) = (I_mdslpes_coronary_sinus(i)-I_mdslpes_coronary_sinus_average)^2;
end
I_mdslpes_coronary_sinus_SD = sqrt(sum(I_mdslpes_coronary_sinus_sq_dev)/(length(I_mdslpes_coronary_sinus)-1));
for i = 1:length(I_mdslpes_z_scores_coronary_sinus)
    I_mdslpes_z_scores_coronary_sinus(i) = (I_mdslpes(i) - I_mdslpes_coronary_sinus_average)/I_mdslpes_coronary_sinus_SD;
end
I_concavs_z_scores_coronary_sinus = zeros(length(I_concavs),1);
I_concavs_coronary_sinus_average = sum(I_concavs_coronary_sinus)/length(I_concavs_coronary_sinus);
I_concavs_coronary_sinus_sq_dev = zeros(length(I_concavs_coronary_sinus),1);
for i = 1:length(I_concavs_coronary_sinus_sq_dev)
    I_concavs_coronary_sinus_sq_dev(i) = (I_concavs_coronary_sinus(i)-I_concavs_coronary_sinus_average)^2;
end
I_concavs_coronary_sinus_SD = sqrt(sum(I_concavs_coronary_sinus_sq_dev)/(length(I_concavs_coronary_sinus)-1));
for i = 1:length(I_concavs_z_scores_coronary_sinus)
    I_concavs_z_scores_coronary_sinus(i) = (I_concavs(i) - I_concavs_coronary_sinus_average)/I_concavs_coronary_sinus_SD;
end
I_sym_inds_z_scores_coronary_sinus = zeros(length(I_sym_inds),1);
I_sym_inds_coronary_sinus_average = sum(I_sym_inds_coronary_sinus)/length(I_sym_inds_coronary_sinus);
I_sym_inds_coronary_sinus_sq_dev = zeros(length(I_sym_inds_coronary_sinus),1);
for i = 1:length(I_sym_inds_coronary_sinus_sq_dev)
    I_sym_inds_coronary_sinus_sq_dev(i) = (I_sym_inds_coronary_sinus(i)-I_sym_inds_coronary_sinus_average)^2;
end
I_sym_inds_coronary_sinus_SD = sqrt(sum(I_sym_inds_coronary_sinus_sq_dev)/(length(I_sym_inds_coronary_sinus)-1));
for i = 1:length(I_sym_inds_z_scores_coronary_sinus)
    I_sym_inds_z_scores_coronary_sinus(i) = (I_sym_inds(i) - I_sym_inds_coronary_sinus_average)/I_sym_inds_coronary_sinus_SD;
end
I_p_wave_areas_1_z_scores_coronary_sinus = zeros(length(I_p_wave_areas_1),1);
I_p_wave_areas_1_coronary_sinus_average = sum(I_p_wave_areas_1_coronary_sinus)/length(I_p_wave_areas_1_coronary_sinus);
I_p_wave_areas_1_coronary_sinus_sq_dev = zeros(length(I_p_wave_areas_1_coronary_sinus),1);
for i = 1:length(I_p_wave_areas_1_coronary_sinus_sq_dev)
    I_p_wave_areas_1_coronary_sinus_sq_dev(i) = (I_p_wave_areas_1_coronary_sinus(i)-I_p_wave_areas_1_coronary_sinus_average)^2;
end
I_p_wave_areas_1_coronary_sinus_SD = sqrt(sum(I_p_wave_areas_1_coronary_sinus_sq_dev)/(length(I_p_wave_areas_1_coronary_sinus)-1));
for i = 1:length(I_p_wave_areas_1_z_scores_coronary_sinus)
    I_p_wave_areas_1_z_scores_coronary_sinus(i) = (I_p_wave_areas_1(i) - I_p_wave_areas_1_coronary_sinus_average)/I_p_wave_areas_1_coronary_sinus_SD;
end
I_p_wave_areas_2_z_scores_coronary_sinus = zeros(length(I_p_wave_areas_2),1);
I_p_wave_areas_2_coronary_sinus_average = sum(I_p_wave_areas_2_coronary_sinus)/length(I_p_wave_areas_2_coronary_sinus);
I_p_wave_areas_2_coronary_sinus_sq_dev = zeros(length(I_p_wave_areas_2_coronary_sinus),1);
for i = 1:length(I_p_wave_areas_2_coronary_sinus_sq_dev)
    I_p_wave_areas_2_coronary_sinus_sq_dev(i) = (I_p_wave_areas_2_coronary_sinus(i)-I_p_wave_areas_2_coronary_sinus_average)^2;
end
I_p_wave_areas_2_coronary_sinus_SD = sqrt(sum(I_p_wave_areas_2_coronary_sinus_sq_dev)/(length(I_p_wave_areas_2_coronary_sinus)-1));
for i = 1:length(I_p_wave_areas_2_z_scores_coronary_sinus)
    I_p_wave_areas_2_z_scores_coronary_sinus(i) = (I_p_wave_areas_2(i) - I_p_wave_areas_2_coronary_sinus_average)/I_p_wave_areas_2_coronary_sinus_SD;
end
I_p_wave_areas_3_z_scores_coronary_sinus = zeros(length(I_p_wave_areas_3),1);
I_p_wave_areas_3_coronary_sinus_average = sum(I_p_wave_areas_3_coronary_sinus)/length(I_p_wave_areas_3_coronary_sinus);
I_p_wave_areas_3_coronary_sinus_sq_dev = zeros(length(I_p_wave_areas_3_coronary_sinus),1);
for i = 1:length(I_p_wave_areas_3_coronary_sinus_sq_dev)
    I_p_wave_areas_3_coronary_sinus_sq_dev(i) = (I_p_wave_areas_3_coronary_sinus(i)-I_p_wave_areas_3_coronary_sinus_average)^3;
end
I_p_wave_areas_3_coronary_sinus_SD = sqrt(sum(I_p_wave_areas_3_coronary_sinus_sq_dev)/(length(I_p_wave_areas_3_coronary_sinus)-1));
for i = 1:length(I_p_wave_areas_3_z_scores_coronary_sinus)
    I_p_wave_areas_3_z_scores_coronary_sinus(i) = (I_p_wave_areas_3(i) - I_p_wave_areas_3_coronary_sinus_average)/I_p_wave_areas_3_coronary_sinus_SD;
end
I_p_wave_areas_4_z_scores_coronary_sinus = zeros(length(I_p_wave_areas_4),1);
I_p_wave_areas_4_coronary_sinus_average = sum(I_p_wave_areas_4_coronary_sinus)/length(I_p_wave_areas_4_coronary_sinus);
I_p_wave_areas_4_coronary_sinus_sq_dev = zeros(length(I_p_wave_areas_4_coronary_sinus),1);
for i = 1:length(I_p_wave_areas_4_coronary_sinus_sq_dev)
    I_p_wave_areas_4_coronary_sinus_sq_dev(i) = (I_p_wave_areas_4_coronary_sinus(i)-I_p_wave_areas_4_coronary_sinus_average)^4;
end
I_p_wave_areas_4_coronary_sinus_SD = sqrt(sum(I_p_wave_areas_4_coronary_sinus_sq_dev)/(length(I_p_wave_areas_4_coronary_sinus)-1));
for i = 1:length(I_p_wave_areas_4_z_scores_coronary_sinus)
    I_p_wave_areas_4_z_scores_coronary_sinus(i) = (I_p_wave_areas_4(i) - I_p_wave_areas_4_coronary_sinus_average)/I_p_wave_areas_4_coronary_sinus_SD;
end
I_p_wave_areas_5_z_scores_coronary_sinus = zeros(length(I_p_wave_areas_5),1);
I_p_wave_areas_5_coronary_sinus_average = sum(I_p_wave_areas_5_coronary_sinus)/length(I_p_wave_areas_5_coronary_sinus);
I_p_wave_areas_5_coronary_sinus_sq_dev = zeros(length(I_p_wave_areas_5_coronary_sinus),1);
for i = 1:length(I_p_wave_areas_5_coronary_sinus_sq_dev)
    I_p_wave_areas_5_coronary_sinus_sq_dev(i) = (I_p_wave_areas_5_coronary_sinus(i)-I_p_wave_areas_5_coronary_sinus_average)^5;
end
I_p_wave_areas_5_coronary_sinus_SD = sqrt(sum(I_p_wave_areas_5_coronary_sinus_sq_dev)/(length(I_p_wave_areas_5_coronary_sinus)-1));
for i = 1:length(I_p_wave_areas_5_z_scores_coronary_sinus)
    I_p_wave_areas_5_z_scores_coronary_sinus(i) = (I_p_wave_areas_5(i) - I_p_wave_areas_5_coronary_sinus_average)/I_p_wave_areas_5_coronary_sinus_SD;
end
I_num_zeros_z_scores_ostium = zeros(length(I_number_of_zeros_array),1);
I_num_zeros_ostium_average = sum(I_num_zeros_ostium)/length(I_num_zeros_ostium);
I_num_zeros_ostium_sq_dev = zeros(length(I_num_zeros_ostium),1);
for i = 1:length(I_num_zeros_ostium_sq_dev)
    I_num_zeros_ostium_sq_dev(i) = (I_num_zeros_ostium(i)-I_num_zeros_ostium_average)^2;
end
I_num_zeros_ostium_SD = sqrt(sum(I_num_zeros_ostium_sq_dev)/(length(I_num_zeros_ostium)-1));
for i = 1:length(I_num_zeros_z_scores_ostium)
    I_num_zeros_z_scores_ostium(i) = (I_number_of_zeros_array(i) - I_num_zeros_ostium_average)/I_num_zeros_ostium_SD;
end
I_maxima_z_scores_ostium = zeros(length(I_maxima),1);
I_maxima_ostium_average = sum(I_maxima_ostium)/length(I_maxima_ostium);
I_maxima_ostium_sq_dev = zeros(length(I_maxima_ostium),1);
for i = 1:length(I_maxima_ostium_sq_dev)
    I_maxima_ostium_sq_dev(i) = (I_maxima_ostium(i)-I_maxima_ostium_average)^2;
end
I_maxima_ostium_SD = sqrt(sum(I_maxima_ostium_sq_dev)/(length(I_maxima_ostium)-1));
for i = 1:length(I_maxima_z_scores_ostium)
    I_maxima_z_scores_ostium(i) = (I_maxima(i) - I_maxima_ostium_average)/I_maxima_ostium_SD;
end
I_minima_z_scores_ostium = zeros(length(I_minima),1);
I_minima_ostium_average = sum(I_minima_ostium)/length(I_minima_ostium);
I_minima_ostium_sq_dev = zeros(length(I_minima_ostium),1);
for i = 1:length(I_minima_ostium_sq_dev)
    I_minima_ostium_sq_dev(i) = (I_minima_ostium(i)-I_minima_ostium_average)^2;
end
I_minima_ostium_SD = sqrt(sum(I_minima_ostium_sq_dev)/(length(I_minima_ostium)-1));
for i = 1:length(I_minima_z_scores_ostium)
    I_minima_z_scores_ostium(i) = (I_minima(i) - I_minima_ostium_average)/I_minima_ostium_SD;
end
I_mdslpes_z_scores_ostium = zeros(length(I_mdslpes),1);
I_mdslpes_ostium_average = sum(I_mdslpes_ostium)/length(I_mdslpes_ostium);
I_mdslpes_ostium_sq_dev = zeros(length(I_mdslpes_ostium),1);
for i = 1:length(I_mdslpes_ostium_sq_dev)
    I_mdslpes_ostium_sq_dev(i) = (I_mdslpes_ostium(i)-I_mdslpes_ostium_average)^2;
end
I_mdslpes_ostium_SD = sqrt(sum(I_mdslpes_ostium_sq_dev)/(length(I_mdslpes_ostium)-1));
for i = 1:length(I_mdslpes_z_scores_ostium)
    I_mdslpes_z_scores_ostium(i) = (I_mdslpes(i) - I_mdslpes_ostium_average)/I_mdslpes_ostium_SD;
end
I_concavs_z_scores_ostium = zeros(length(I_concavs),1);
I_concavs_ostium_average = sum(I_concavs_ostium)/length(I_concavs_ostium);
I_concavs_ostium_sq_dev = zeros(length(I_concavs_ostium),1);
for i = 1:length(I_concavs_ostium_sq_dev)
    I_concavs_ostium_sq_dev(i) = (I_concavs_ostium(i)-I_concavs_ostium_average)^2;
end
I_concavs_ostium_SD = sqrt(sum(I_concavs_ostium_sq_dev)/(length(I_concavs_ostium)-1));
for i = 1:length(I_concavs_z_scores_ostium)
    I_concavs_z_scores_ostium(i) = (I_concavs(i) - I_concavs_ostium_average)/I_concavs_ostium_SD;
end
I_sym_inds_z_scores_ostium = zeros(length(I_sym_inds),1);
I_sym_inds_ostium_average = sum(I_sym_inds_ostium)/length(I_sym_inds_ostium);
I_sym_inds_ostium_sq_dev = zeros(length(I_sym_inds_ostium),1);
for i = 1:length(I_sym_inds_ostium_sq_dev)
    I_sym_inds_ostium_sq_dev(i) = (I_sym_inds_ostium(i)-I_sym_inds_ostium_average)^2;
end
I_sym_inds_ostium_SD = sqrt(sum(I_sym_inds_ostium_sq_dev)/(length(I_sym_inds_ostium)-1));
for i = 1:length(I_sym_inds_z_scores_ostium)
    I_sym_inds_z_scores_ostium(i) = (I_sym_inds(i) - I_sym_inds_ostium_average)/I_sym_inds_ostium_SD;
end
I_p_wave_areas_1_z_scores_ostium = zeros(length(I_p_wave_areas_1),1);
I_p_wave_areas_1_ostium_average = sum(I_p_wave_areas_1_ostium)/length(I_p_wave_areas_1_ostium);
I_p_wave_areas_1_ostium_sq_dev = zeros(length(I_p_wave_areas_1_ostium),1);
for i = 1:length(I_p_wave_areas_1_ostium_sq_dev)
    I_p_wave_areas_1_ostium_sq_dev(i) = (I_p_wave_areas_1_ostium(i)-I_p_wave_areas_1_ostium_average)^2;
end
I_p_wave_areas_1_ostium_SD = sqrt(sum(I_p_wave_areas_1_ostium_sq_dev)/(length(I_p_wave_areas_1_ostium)-1));
for i = 1:length(I_p_wave_areas_1_z_scores_ostium)
    I_p_wave_areas_1_z_scores_ostium(i) = (I_p_wave_areas_1(i) - I_p_wave_areas_1_ostium_average)/I_p_wave_areas_1_ostium_SD;
end
I_p_wave_areas_2_z_scores_ostium = zeros(length(I_p_wave_areas_2),1);
I_p_wave_areas_2_ostium_average = sum(I_p_wave_areas_2_ostium)/length(I_p_wave_areas_2_ostium);
I_p_wave_areas_2_ostium_sq_dev = zeros(length(I_p_wave_areas_2_ostium),1);
for i = 1:length(I_p_wave_areas_2_ostium_sq_dev)
    I_p_wave_areas_2_ostium_sq_dev(i) = (I_p_wave_areas_2_ostium(i)-I_p_wave_areas_2_ostium_average)^2;
end
I_p_wave_areas_2_ostium_SD = sqrt(sum(I_p_wave_areas_2_ostium_sq_dev)/(length(I_p_wave_areas_2_ostium)-1));
for i = 1:length(I_p_wave_areas_2_z_scores_ostium)
    I_p_wave_areas_2_z_scores_ostium(i) = (I_p_wave_areas_2(i) - I_p_wave_areas_2_ostium_average)/I_p_wave_areas_2_ostium_SD;
end
I_p_wave_areas_3_z_scores_ostium = zeros(length(I_p_wave_areas_3),1);
I_p_wave_areas_3_ostium_average = sum(I_p_wave_areas_3_ostium)/length(I_p_wave_areas_3_ostium);
I_p_wave_areas_3_ostium_sq_dev = zeros(length(I_p_wave_areas_3_ostium),1);
for i = 1:length(I_p_wave_areas_3_ostium_sq_dev)
    I_p_wave_areas_3_ostium_sq_dev(i) = (I_p_wave_areas_3_ostium(i)-I_p_wave_areas_3_ostium_average)^3;
end
I_p_wave_areas_3_ostium_SD = sqrt(sum(I_p_wave_areas_3_ostium_sq_dev)/(length(I_p_wave_areas_3_ostium)-1));
for i = 1:length(I_p_wave_areas_3_z_scores_ostium)
    I_p_wave_areas_3_z_scores_ostium(i) = (I_p_wave_areas_3(i) - I_p_wave_areas_3_ostium_average)/I_p_wave_areas_3_ostium_SD;
end
I_p_wave_areas_4_z_scores_ostium = zeros(length(I_p_wave_areas_4),1);
I_p_wave_areas_4_ostium_average = sum(I_p_wave_areas_4_ostium)/length(I_p_wave_areas_4_ostium);
I_p_wave_areas_4_ostium_sq_dev = zeros(length(I_p_wave_areas_4_ostium),1);
for i = 1:length(I_p_wave_areas_4_ostium_sq_dev)
    I_p_wave_areas_4_ostium_sq_dev(i) = (I_p_wave_areas_4_ostium(i)-I_p_wave_areas_4_ostium_average)^4;
end
I_p_wave_areas_4_ostium_SD = sqrt(sum(I_p_wave_areas_4_ostium_sq_dev)/(length(I_p_wave_areas_4_ostium)-1));
for i = 1:length(I_p_wave_areas_4_z_scores_ostium)
    I_p_wave_areas_4_z_scores_ostium(i) = (I_p_wave_areas_4(i) - I_p_wave_areas_4_ostium_average)/I_p_wave_areas_4_ostium_SD;
end
I_p_wave_areas_5_z_scores_ostium = zeros(length(I_p_wave_areas_5),1);
I_p_wave_areas_5_ostium_average = sum(I_p_wave_areas_5_ostium)/length(I_p_wave_areas_5_ostium);
I_p_wave_areas_5_ostium_sq_dev = zeros(length(I_p_wave_areas_5_ostium),1);
for i = 1:length(I_p_wave_areas_5_ostium_sq_dev)
    I_p_wave_areas_5_ostium_sq_dev(i) = (I_p_wave_areas_5_ostium(i)-I_p_wave_areas_5_ostium_average)^5;
end
I_p_wave_areas_5_ostium_SD = sqrt(sum(I_p_wave_areas_5_ostium_sq_dev)/(length(I_p_wave_areas_5_ostium)-1));
for i = 1:length(I_p_wave_areas_5_z_scores_ostium)
    I_p_wave_areas_5_z_scores_ostium(i) = (I_p_wave_areas_5(i) - I_p_wave_areas_5_ostium_average)/I_p_wave_areas_5_ostium_SD;
end
I_num_zeros_z_scores_perinodal = zeros(length(I_number_of_zeros_array),1);
I_num_zeros_perinodal_average = sum(I_num_zeros_perinodal)/length(I_num_zeros_perinodal);
I_num_zeros_perinodal_sq_dev = zeros(length(I_num_zeros_perinodal),1);
for i = 1:length(I_num_zeros_perinodal_sq_dev)
    I_num_zeros_perinodal_sq_dev(i) = (I_num_zeros_perinodal(i)-I_num_zeros_perinodal_average)^2;
end
I_num_zeros_perinodal_SD = sqrt(sum(I_num_zeros_perinodal_sq_dev)/(length(I_num_zeros_perinodal)-1));
for i = 1:length(I_num_zeros_z_scores_perinodal)
    I_num_zeros_z_scores_perinodal(i) = (I_number_of_zeros_array(i) - I_num_zeros_perinodal_average)/I_num_zeros_perinodal_SD;
end
I_maxima_z_scores_perinodal = zeros(length(I_maxima),1);
I_maxima_perinodal_average = sum(I_maxima_perinodal)/length(I_maxima_perinodal);
I_maxima_perinodal_sq_dev = zeros(length(I_maxima_perinodal),1);
for i = 1:length(I_maxima_perinodal_sq_dev)
    I_maxima_perinodal_sq_dev(i) = (I_maxima_perinodal(i)-I_maxima_perinodal_average)^2;
end
I_maxima_perinodal_SD = sqrt(sum(I_maxima_perinodal_sq_dev)/(length(I_maxima_perinodal)-1));
for i = 1:length(I_maxima_z_scores_perinodal)
    I_maxima_z_scores_perinodal(i) = (I_maxima(i) - I_maxima_perinodal_average)/I_maxima_perinodal_SD;
end
I_minima_z_scores_perinodal = zeros(length(I_minima),1);
I_minima_perinodal_average = sum(I_minima_perinodal)/length(I_minima_perinodal);
I_minima_perinodal_sq_dev = zeros(length(I_minima_perinodal),1);
for i = 1:length(I_minima_perinodal_sq_dev)
    I_minima_perinodal_sq_dev(i) = (I_minima_perinodal(i)-I_minima_perinodal_average)^2;
end
I_minima_perinodal_SD = sqrt(sum(I_minima_perinodal_sq_dev)/(length(I_minima_perinodal)-1));
for i = 1:length(I_minima_z_scores_perinodal)
    I_minima_z_scores_perinodal(i) = (I_minima(i) - I_minima_perinodal_average)/I_minima_perinodal_SD;
end
I_mdslpes_z_scores_perinodal = zeros(length(I_mdslpes),1);
I_mdslpes_perinodal_average = sum(I_mdslpes_perinodal)/length(I_mdslpes_perinodal);
I_mdslpes_perinodal_sq_dev = zeros(length(I_mdslpes_perinodal),1);
for i = 1:length(I_mdslpes_perinodal_sq_dev)
    I_mdslpes_perinodal_sq_dev(i) = (I_mdslpes_perinodal(i)-I_mdslpes_perinodal_average)^2;
end
I_mdslpes_perinodal_SD = sqrt(sum(I_mdslpes_perinodal_sq_dev)/(length(I_mdslpes_perinodal)-1));
for i = 1:length(I_mdslpes_z_scores_perinodal)
    I_mdslpes_z_scores_perinodal(i) = (I_mdslpes(i) - I_mdslpes_perinodal_average)/I_mdslpes_perinodal_SD;
end
I_concavs_z_scores_perinodal = zeros(length(I_concavs),1);
I_concavs_perinodal_average = sum(I_concavs_perinodal)/length(I_concavs_perinodal);
I_concavs_perinodal_sq_dev = zeros(length(I_concavs_perinodal),1);
for i = 1:length(I_concavs_perinodal_sq_dev)
    I_concavs_perinodal_sq_dev(i) = (I_concavs_perinodal(i)-I_concavs_perinodal_average)^2;
end
I_concavs_perinodal_SD = sqrt(sum(I_concavs_perinodal_sq_dev)/(length(I_concavs_perinodal)-1));
for i = 1:length(I_concavs_z_scores_perinodal)
    I_concavs_z_scores_perinodal(i) = (I_concavs(i) - I_concavs_perinodal_average)/I_concavs_perinodal_SD;
end
I_sym_inds_z_scores_perinodal = zeros(length(I_sym_inds),1);
I_sym_inds_perinodal_average = sum(I_sym_inds_perinodal)/length(I_sym_inds_perinodal);
I_sym_inds_perinodal_sq_dev = zeros(length(I_sym_inds_perinodal),1);
for i = 1:length(I_sym_inds_perinodal_sq_dev)
    I_sym_inds_perinodal_sq_dev(i) = (I_sym_inds_perinodal(i)-I_sym_inds_perinodal_average)^2;
end
I_sym_inds_perinodal_SD = sqrt(sum(I_sym_inds_perinodal_sq_dev)/(length(I_sym_inds_perinodal)-1));
for i = 1:length(I_sym_inds_z_scores_perinodal)
    I_sym_inds_z_scores_perinodal(i) = (I_sym_inds(i) - I_sym_inds_perinodal_average)/I_sym_inds_perinodal_SD;
end
I_p_wave_areas_1_z_scores_perinodal = zeros(length(I_p_wave_areas_1),1);
I_p_wave_areas_1_perinodal_average = sum(I_p_wave_areas_1_perinodal)/length(I_p_wave_areas_1_perinodal);
I_p_wave_areas_1_perinodal_sq_dev = zeros(length(I_p_wave_areas_1_perinodal),1);
for i = 1:length(I_p_wave_areas_1_perinodal_sq_dev)
    I_p_wave_areas_1_perinodal_sq_dev(i) = (I_p_wave_areas_1_perinodal(i)-I_p_wave_areas_1_perinodal_average)^2;
end
I_p_wave_areas_1_perinodal_SD = sqrt(sum(I_p_wave_areas_1_perinodal_sq_dev)/(length(I_p_wave_areas_1_perinodal)-1));
for i = 1:length(I_p_wave_areas_1_z_scores_perinodal)
    I_p_wave_areas_1_z_scores_perinodal(i) = (I_p_wave_areas_1(i) - I_p_wave_areas_1_perinodal_average)/I_p_wave_areas_1_perinodal_SD;
end
I_p_wave_areas_2_z_scores_perinodal = zeros(length(I_p_wave_areas_2),1);
I_p_wave_areas_2_perinodal_average = sum(I_p_wave_areas_2_perinodal)/length(I_p_wave_areas_2_perinodal);
I_p_wave_areas_2_perinodal_sq_dev = zeros(length(I_p_wave_areas_2_perinodal),1);
for i = 1:length(I_p_wave_areas_2_perinodal_sq_dev)
    I_p_wave_areas_2_perinodal_sq_dev(i) = (I_p_wave_areas_2_perinodal(i)-I_p_wave_areas_2_perinodal_average)^2;
end
I_p_wave_areas_2_perinodal_SD = sqrt(sum(I_p_wave_areas_2_perinodal_sq_dev)/(length(I_p_wave_areas_2_perinodal)-1));
for i = 1:length(I_p_wave_areas_2_z_scores_perinodal)
    I_p_wave_areas_2_z_scores_perinodal(i) = (I_p_wave_areas_2(i) - I_p_wave_areas_2_perinodal_average)/I_p_wave_areas_2_perinodal_SD;
end
I_p_wave_areas_3_z_scores_perinodal = zeros(length(I_p_wave_areas_3),1);
I_p_wave_areas_3_perinodal_average = sum(I_p_wave_areas_3_perinodal)/length(I_p_wave_areas_3_perinodal);
I_p_wave_areas_3_perinodal_sq_dev = zeros(length(I_p_wave_areas_3_perinodal),1);
for i = 1:length(I_p_wave_areas_3_perinodal_sq_dev)
    I_p_wave_areas_3_perinodal_sq_dev(i) = (I_p_wave_areas_3_perinodal(i)-I_p_wave_areas_3_perinodal_average)^3;
end
I_p_wave_areas_3_perinodal_SD = sqrt(sum(I_p_wave_areas_3_perinodal_sq_dev)/(length(I_p_wave_areas_3_perinodal)-1));
for i = 1:length(I_p_wave_areas_3_z_scores_perinodal)
    I_p_wave_areas_3_z_scores_perinodal(i) = (I_p_wave_areas_3(i) - I_p_wave_areas_3_perinodal_average)/I_p_wave_areas_3_perinodal_SD;
end
I_p_wave_areas_4_z_scores_perinodal = zeros(length(I_p_wave_areas_4),1);
I_p_wave_areas_4_perinodal_average = sum(I_p_wave_areas_4_perinodal)/length(I_p_wave_areas_4_perinodal);
I_p_wave_areas_4_perinodal_sq_dev = zeros(length(I_p_wave_areas_4_perinodal),1);
for i = 1:length(I_p_wave_areas_4_perinodal_sq_dev)
    I_p_wave_areas_4_perinodal_sq_dev(i) = (I_p_wave_areas_4_perinodal(i)-I_p_wave_areas_4_perinodal_average)^4;
end
I_p_wave_areas_4_perinodal_SD = sqrt(sum(I_p_wave_areas_4_perinodal_sq_dev)/(length(I_p_wave_areas_4_perinodal)-1));
for i = 1:length(I_p_wave_areas_4_z_scores_perinodal)
    I_p_wave_areas_4_z_scores_perinodal(i) = (I_p_wave_areas_4(i) - I_p_wave_areas_4_perinodal_average)/I_p_wave_areas_4_perinodal_SD;
end
I_p_wave_areas_5_z_scores_perinodal = zeros(length(I_p_wave_areas_5),1);
I_p_wave_areas_5_perinodal_average = sum(I_p_wave_areas_5_perinodal)/length(I_p_wave_areas_5_perinodal);
I_p_wave_areas_5_perinodal_sq_dev = zeros(length(I_p_wave_areas_5_perinodal),1);
for i = 1:length(I_p_wave_areas_5_perinodal_sq_dev)
    I_p_wave_areas_5_perinodal_sq_dev(i) = (I_p_wave_areas_5_perinodal(i)-I_p_wave_areas_5_perinodal_average)^5;
end
I_p_wave_areas_5_perinodal_SD = sqrt(sum(I_p_wave_areas_5_perinodal_sq_dev)/(length(I_p_wave_areas_5_perinodal)-1));
for i = 1:length(I_p_wave_areas_5_z_scores_perinodal)
    I_p_wave_areas_5_z_scores_perinodal(i) = (I_p_wave_areas_5(i) - I_p_wave_areas_5_perinodal_average)/I_p_wave_areas_5_perinodal_SD;
end
I_num_zeros_z_scores_right_septum = zeros(length(I_number_of_zeros_array),1);
I_num_zeros_right_septum_average = sum(I_num_zeros_right_septum)/length(I_num_zeros_right_septum);
I_num_zeros_right_septum_sq_dev = zeros(length(I_num_zeros_right_septum),1);
for i = 1:length(I_num_zeros_right_septum_sq_dev)
    I_num_zeros_right_septum_sq_dev(i) = (I_num_zeros_right_septum(i)-I_num_zeros_right_septum_average)^2;
end
I_num_zeros_right_septum_SD = sqrt(sum(I_num_zeros_right_septum_sq_dev)/(length(I_num_zeros_right_septum)-1));
for i = 1:length(I_num_zeros_z_scores_right_septum)
    I_num_zeros_z_scores_right_septum(i) = (I_number_of_zeros_array(i) - I_num_zeros_right_septum_average)/I_num_zeros_right_septum_SD;
end
I_maxima_z_scores_right_septum = zeros(length(I_maxima),1);
I_maxima_right_septum_average = sum(I_maxima_right_septum)/length(I_maxima_right_septum);
I_maxima_right_septum_sq_dev = zeros(length(I_maxima_right_septum),1);
for i = 1:length(I_maxima_right_septum_sq_dev)
    I_maxima_right_septum_sq_dev(i) = (I_maxima_right_septum(i)-I_maxima_right_septum_average)^2;
end
I_maxima_right_septum_SD = sqrt(sum(I_maxima_right_septum_sq_dev)/(length(I_maxima_right_septum)-1));
for i = 1:length(I_maxima_z_scores_right_septum)
    I_maxima_z_scores_right_septum(i) = (I_maxima(i) - I_maxima_right_septum_average)/I_maxima_right_septum_SD;
end
I_minima_z_scores_right_septum = zeros(length(I_minima),1);
I_minima_right_septum_average = sum(I_minima_right_septum)/length(I_minima_right_septum);
I_minima_right_septum_sq_dev = zeros(length(I_minima_right_septum),1);
for i = 1:length(I_minima_right_septum_sq_dev)
    I_minima_right_septum_sq_dev(i) = (I_minima_right_septum(i)-I_minima_right_septum_average)^2;
end
I_minima_right_septum_SD = sqrt(sum(I_minima_right_septum_sq_dev)/(length(I_minima_right_septum)-1));
for i = 1:length(I_minima_z_scores_right_septum)
    I_minima_z_scores_right_septum(i) = (I_minima(i) - I_minima_right_septum_average)/I_minima_right_septum_SD;
end
I_mdslpes_z_scores_right_septum = zeros(length(I_mdslpes),1);
I_mdslpes_right_septum_average = sum(I_mdslpes_right_septum)/length(I_mdslpes_right_septum);
I_mdslpes_right_septum_sq_dev = zeros(length(I_mdslpes_right_septum),1);
for i = 1:length(I_mdslpes_right_septum_sq_dev)
    I_mdslpes_right_septum_sq_dev(i) = (I_mdslpes_right_septum(i)-I_mdslpes_right_septum_average)^2;
end
I_mdslpes_right_septum_SD = sqrt(sum(I_mdslpes_right_septum_sq_dev)/(length(I_mdslpes_right_septum)-1));
for i = 1:length(I_mdslpes_z_scores_right_septum)
    I_mdslpes_z_scores_right_septum(i) = (I_mdslpes(i) - I_mdslpes_right_septum_average)/I_mdslpes_right_septum_SD;
end
I_concavs_z_scores_right_septum = zeros(length(I_concavs),1);
I_concavs_right_septum_average = sum(I_concavs_right_septum)/length(I_concavs_right_septum);
I_concavs_right_septum_sq_dev = zeros(length(I_concavs_right_septum),1);
for i = 1:length(I_concavs_right_septum_sq_dev)
    I_concavs_right_septum_sq_dev(i) = (I_concavs_right_septum(i)-I_concavs_right_septum_average)^2;
end
I_concavs_right_septum_SD = sqrt(sum(I_concavs_right_septum_sq_dev)/(length(I_concavs_right_septum)-1));
for i = 1:length(I_concavs_z_scores_right_septum)
    I_concavs_z_scores_right_septum(i) = (I_concavs(i) - I_concavs_right_septum_average)/I_concavs_right_septum_SD;
end
I_sym_inds_z_scores_right_septum = zeros(length(I_sym_inds),1);
I_sym_inds_right_septum_average = sum(I_sym_inds_right_septum)/length(I_sym_inds_right_septum);
I_sym_inds_right_septum_sq_dev = zeros(length(I_sym_inds_right_septum),1);
for i = 1:length(I_sym_inds_right_septum_sq_dev)
    I_sym_inds_right_septum_sq_dev(i) = (I_sym_inds_right_septum(i)-I_sym_inds_right_septum_average)^2;
end
I_sym_inds_right_septum_SD = sqrt(sum(I_sym_inds_right_septum_sq_dev)/(length(I_sym_inds_right_septum)-1));
for i = 1:length(I_sym_inds_z_scores_right_septum)
    I_sym_inds_z_scores_right_septum(i) = (I_sym_inds(i) - I_sym_inds_right_septum_average)/I_sym_inds_right_septum_SD;
end
I_p_wave_areas_1_z_scores_right_septum = zeros(length(I_p_wave_areas_1),1);
I_p_wave_areas_1_right_septum_average = sum(I_p_wave_areas_1_right_septum)/length(I_p_wave_areas_1_right_septum);
I_p_wave_areas_1_right_septum_sq_dev = zeros(length(I_p_wave_areas_1_right_septum),1);
for i = 1:length(I_p_wave_areas_1_right_septum_sq_dev)
    I_p_wave_areas_1_right_septum_sq_dev(i) = (I_p_wave_areas_1_right_septum(i)-I_p_wave_areas_1_right_septum_average)^2;
end
I_p_wave_areas_1_right_septum_SD = sqrt(sum(I_p_wave_areas_1_right_septum_sq_dev)/(length(I_p_wave_areas_1_right_septum)-1));
for i = 1:length(I_p_wave_areas_1_z_scores_right_septum)
    I_p_wave_areas_1_z_scores_right_septum(i) = (I_p_wave_areas_1(i) - I_p_wave_areas_1_right_septum_average)/I_p_wave_areas_1_right_septum_SD;
end
I_p_wave_areas_2_z_scores_right_septum = zeros(length(I_p_wave_areas_2),1);
I_p_wave_areas_2_right_septum_average = sum(I_p_wave_areas_2_right_septum)/length(I_p_wave_areas_2_right_septum);
I_p_wave_areas_2_right_septum_sq_dev = zeros(length(I_p_wave_areas_2_right_septum),1);
for i = 1:length(I_p_wave_areas_2_right_septum_sq_dev)
    I_p_wave_areas_2_right_septum_sq_dev(i) = (I_p_wave_areas_2_right_septum(i)-I_p_wave_areas_2_right_septum_average)^2;
end
I_p_wave_areas_2_right_septum_SD = sqrt(sum(I_p_wave_areas_2_right_septum_sq_dev)/(length(I_p_wave_areas_2_right_septum)-1));
for i = 1:length(I_p_wave_areas_2_z_scores_right_septum)
    I_p_wave_areas_2_z_scores_right_septum(i) = (I_p_wave_areas_2(i) - I_p_wave_areas_2_right_septum_average)/I_p_wave_areas_2_right_septum_SD;
end
I_p_wave_areas_3_z_scores_right_septum = zeros(length(I_p_wave_areas_3),1);
I_p_wave_areas_3_right_septum_average = sum(I_p_wave_areas_3_right_septum)/length(I_p_wave_areas_3_right_septum);
I_p_wave_areas_3_right_septum_sq_dev = zeros(length(I_p_wave_areas_3_right_septum),1);
for i = 1:length(I_p_wave_areas_3_right_septum_sq_dev)
    I_p_wave_areas_3_right_septum_sq_dev(i) = (I_p_wave_areas_3_right_septum(i)-I_p_wave_areas_3_right_septum_average)^3;
end
I_p_wave_areas_3_right_septum_SD = sqrt(sum(I_p_wave_areas_3_right_septum_sq_dev)/(length(I_p_wave_areas_3_right_septum)-1));
for i = 1:length(I_p_wave_areas_3_z_scores_right_septum)
    I_p_wave_areas_3_z_scores_right_septum(i) = (I_p_wave_areas_3(i) - I_p_wave_areas_3_right_septum_average)/I_p_wave_areas_3_right_septum_SD;
end
I_p_wave_areas_4_z_scores_right_septum = zeros(length(I_p_wave_areas_4),1);
I_p_wave_areas_4_right_septum_average = sum(I_p_wave_areas_4_right_septum)/length(I_p_wave_areas_4_right_septum);
I_p_wave_areas_4_right_septum_sq_dev = zeros(length(I_p_wave_areas_4_right_septum),1);
for i = 1:length(I_p_wave_areas_4_right_septum_sq_dev)
    I_p_wave_areas_4_right_septum_sq_dev(i) = (I_p_wave_areas_4_right_septum(i)-I_p_wave_areas_4_right_septum_average)^4;
end
I_p_wave_areas_4_right_septum_SD = sqrt(sum(I_p_wave_areas_4_right_septum_sq_dev)/(length(I_p_wave_areas_4_right_septum)-1));
for i = 1:length(I_p_wave_areas_4_z_scores_right_septum)
    I_p_wave_areas_4_z_scores_right_septum(i) = (I_p_wave_areas_4(i) - I_p_wave_areas_4_right_septum_average)/I_p_wave_areas_4_right_septum_SD;
end
I_p_wave_areas_5_z_scores_right_septum = zeros(length(I_p_wave_areas_5),1);
I_p_wave_areas_5_right_septum_average = sum(I_p_wave_areas_5_right_septum)/length(I_p_wave_areas_5_right_septum);
I_p_wave_areas_5_right_septum_sq_dev = zeros(length(I_p_wave_areas_5_right_septum),1);
for i = 1:length(I_p_wave_areas_5_right_septum_sq_dev)
    I_p_wave_areas_5_right_septum_sq_dev(i) = (I_p_wave_areas_5_right_septum(i)-I_p_wave_areas_5_right_septum_average)^5;
end
I_p_wave_areas_5_right_septum_SD = sqrt(sum(I_p_wave_areas_5_right_septum_sq_dev)/(length(I_p_wave_areas_5_right_septum)-1));
for i = 1:length(I_p_wave_areas_5_z_scores_right_septum)
    I_p_wave_areas_5_z_scores_right_septum(i) = (I_p_wave_areas_5(i) - I_p_wave_areas_5_right_septum_average)/I_p_wave_areas_5_right_septum_SD;
end
I_num_zeros_z_scores_right_atrial_appendage = zeros(length(I_number_of_zeros_array),1);
I_num_zeros_right_atrial_appendage_average = sum(I_num_zeros_right_atrial_appendage)/length(I_num_zeros_right_atrial_appendage);
I_num_zeros_right_atrial_appendage_sq_dev = zeros(length(I_num_zeros_right_atrial_appendage),1);
for i = 1:length(I_num_zeros_right_atrial_appendage_sq_dev)
    I_num_zeros_right_atrial_appendage_sq_dev(i) = (I_num_zeros_right_atrial_appendage(i)-I_num_zeros_right_atrial_appendage_average)^2;
end
I_num_zeros_right_atrial_appendage_SD = sqrt(sum(I_num_zeros_right_atrial_appendage_sq_dev)/(length(I_num_zeros_right_atrial_appendage)-1));
for i = 1:length(I_num_zeros_z_scores_right_atrial_appendage)
    I_num_zeros_z_scores_right_atrial_appendage(i) = (I_number_of_zeros_array(i) - I_num_zeros_right_atrial_appendage_average)/I_num_zeros_right_atrial_appendage_SD;
end
I_maxima_z_scores_right_atrial_appendage = zeros(length(I_maxima),1);
I_maxima_right_atrial_appendage_average = sum(I_maxima_right_atrial_appendage)/length(I_maxima_right_atrial_appendage);
I_maxima_right_atrial_appendage_sq_dev = zeros(length(I_maxima_right_atrial_appendage),1);
for i = 1:length(I_maxima_right_atrial_appendage_sq_dev)
    I_maxima_right_atrial_appendage_sq_dev(i) = (I_maxima_right_atrial_appendage(i)-I_maxima_right_atrial_appendage_average)^2;
end
I_maxima_right_atrial_appendage_SD = sqrt(sum(I_maxima_right_atrial_appendage_sq_dev)/(length(I_maxima_right_atrial_appendage)-1));
for i = 1:length(I_maxima_z_scores_right_atrial_appendage)
    I_maxima_z_scores_right_atrial_appendage(i) = (I_maxima(i) - I_maxima_right_atrial_appendage_average)/I_maxima_right_atrial_appendage_SD;
end
I_minima_z_scores_right_atrial_appendage = zeros(length(I_minima),1);
I_minima_right_atrial_appendage_average = sum(I_minima_right_atrial_appendage)/length(I_minima_right_atrial_appendage);
I_minima_right_atrial_appendage_sq_dev = zeros(length(I_minima_right_atrial_appendage),1);
for i = 1:length(I_minima_right_atrial_appendage_sq_dev)
    I_minima_right_atrial_appendage_sq_dev(i) = (I_minima_right_atrial_appendage(i)-I_minima_right_atrial_appendage_average)^2;
end
I_minima_right_atrial_appendage_SD = sqrt(sum(I_minima_right_atrial_appendage_sq_dev)/(length(I_minima_right_atrial_appendage)-1));
for i = 1:length(I_minima_z_scores_right_atrial_appendage)
    I_minima_z_scores_right_atrial_appendage(i) = (I_minima(i) - I_minima_right_atrial_appendage_average)/I_minima_right_atrial_appendage_SD;
end
I_mdslpes_z_scores_right_atrial_appendage = zeros(length(I_mdslpes),1);
I_mdslpes_right_atrial_appendage_average = sum(I_mdslpes_right_atrial_appendage)/length(I_mdslpes_right_atrial_appendage);
I_mdslpes_right_atrial_appendage_sq_dev = zeros(length(I_mdslpes_right_atrial_appendage),1);
for i = 1:length(I_mdslpes_right_atrial_appendage_sq_dev)
    I_mdslpes_right_atrial_appendage_sq_dev(i) = (I_mdslpes_right_atrial_appendage(i)-I_mdslpes_right_atrial_appendage_average)^2;
end
I_mdslpes_right_atrial_appendage_SD = sqrt(sum(I_mdslpes_right_atrial_appendage_sq_dev)/(length(I_mdslpes_right_atrial_appendage)-1));
for i = 1:length(I_mdslpes_z_scores_right_atrial_appendage)
    I_mdslpes_z_scores_right_atrial_appendage(i) = (I_mdslpes(i) - I_mdslpes_right_atrial_appendage_average)/I_mdslpes_right_atrial_appendage_SD;
end
I_concavs_z_scores_right_atrial_appendage = zeros(length(I_concavs),1);
I_concavs_right_atrial_appendage_average = sum(I_concavs_right_atrial_appendage)/length(I_concavs_right_atrial_appendage);
I_concavs_right_atrial_appendage_sq_dev = zeros(length(I_concavs_right_atrial_appendage),1);
for i = 1:length(I_concavs_right_atrial_appendage_sq_dev)
    I_concavs_right_atrial_appendage_sq_dev(i) = (I_concavs_right_atrial_appendage(i)-I_concavs_right_atrial_appendage_average)^2;
end
I_concavs_right_atrial_appendage_SD = sqrt(sum(I_concavs_right_atrial_appendage_sq_dev)/(length(I_concavs_right_atrial_appendage)-1));
for i = 1:length(I_concavs_z_scores_right_atrial_appendage)
    I_concavs_z_scores_right_atrial_appendage(i) = (I_concavs(i) - I_concavs_right_atrial_appendage_average)/I_concavs_right_atrial_appendage_SD;
end
I_sym_inds_z_scores_right_atrial_appendage = zeros(length(I_sym_inds),1);
I_sym_inds_right_atrial_appendage_average = sum(I_sym_inds_right_atrial_appendage)/length(I_sym_inds_right_atrial_appendage);
I_sym_inds_right_atrial_appendage_sq_dev = zeros(length(I_sym_inds_right_atrial_appendage),1);
for i = 1:length(I_sym_inds_right_atrial_appendage_sq_dev)
    I_sym_inds_right_atrial_appendage_sq_dev(i) = (I_sym_inds_right_atrial_appendage(i)-I_sym_inds_right_atrial_appendage_average)^2;
end
I_sym_inds_right_atrial_appendage_SD = sqrt(sum(I_sym_inds_right_atrial_appendage_sq_dev)/(length(I_sym_inds_right_atrial_appendage)-1));
for i = 1:length(I_sym_inds_z_scores_right_atrial_appendage)
    I_sym_inds_z_scores_right_atrial_appendage(i) = (I_sym_inds(i) - I_sym_inds_right_atrial_appendage_average)/I_sym_inds_right_atrial_appendage_SD;
end
I_p_wave_areas_1_z_scores_right_atrial_appendage = zeros(length(I_p_wave_areas_1),1);
I_p_wave_areas_1_right_atrial_appendage_average = sum(I_p_wave_areas_1_right_atrial_appendage)/length(I_p_wave_areas_1_right_atrial_appendage);
I_p_wave_areas_1_right_atrial_appendage_sq_dev = zeros(length(I_p_wave_areas_1_right_atrial_appendage),1);
for i = 1:length(I_p_wave_areas_1_right_atrial_appendage_sq_dev)
    I_p_wave_areas_1_right_atrial_appendage_sq_dev(i) = (I_p_wave_areas_1_right_atrial_appendage(i)-I_p_wave_areas_1_right_atrial_appendage_average)^2;
end
I_p_wave_areas_1_right_atrial_appendage_SD = sqrt(sum(I_p_wave_areas_1_right_atrial_appendage_sq_dev)/(length(I_p_wave_areas_1_right_atrial_appendage)-1));
for i = 1:length(I_p_wave_areas_1_z_scores_right_atrial_appendage)
    I_p_wave_areas_1_z_scores_right_atrial_appendage(i) = (I_p_wave_areas_1(i) - I_p_wave_areas_1_right_atrial_appendage_average)/I_p_wave_areas_1_right_atrial_appendage_SD;
end
I_p_wave_areas_2_z_scores_right_atrial_appendage = zeros(length(I_p_wave_areas_2),1);
I_p_wave_areas_2_right_atrial_appendage_average = sum(I_p_wave_areas_2_right_atrial_appendage)/length(I_p_wave_areas_2_right_atrial_appendage);
I_p_wave_areas_2_right_atrial_appendage_sq_dev = zeros(length(I_p_wave_areas_2_right_atrial_appendage),1);
for i = 1:length(I_p_wave_areas_2_right_atrial_appendage_sq_dev)
    I_p_wave_areas_2_right_atrial_appendage_sq_dev(i) = (I_p_wave_areas_2_right_atrial_appendage(i)-I_p_wave_areas_2_right_atrial_appendage_average)^2;
end
I_p_wave_areas_2_right_atrial_appendage_SD = sqrt(sum(I_p_wave_areas_2_right_atrial_appendage_sq_dev)/(length(I_p_wave_areas_2_right_atrial_appendage)-1));
for i = 1:length(I_p_wave_areas_2_z_scores_right_atrial_appendage)
    I_p_wave_areas_2_z_scores_right_atrial_appendage(i) = (I_p_wave_areas_2(i) - I_p_wave_areas_2_right_atrial_appendage_average)/I_p_wave_areas_2_right_atrial_appendage_SD;
end
I_p_wave_areas_3_z_scores_right_atrial_appendage = zeros(length(I_p_wave_areas_3),1);
I_p_wave_areas_3_right_atrial_appendage_average = sum(I_p_wave_areas_3_right_atrial_appendage)/length(I_p_wave_areas_3_right_atrial_appendage);
I_p_wave_areas_3_right_atrial_appendage_sq_dev = zeros(length(I_p_wave_areas_3_right_atrial_appendage),1);
for i = 1:length(I_p_wave_areas_3_right_atrial_appendage_sq_dev)
    I_p_wave_areas_3_right_atrial_appendage_sq_dev(i) = (I_p_wave_areas_3_right_atrial_appendage(i)-I_p_wave_areas_3_right_atrial_appendage_average)^3;
end
I_p_wave_areas_3_right_atrial_appendage_SD = sqrt(sum(I_p_wave_areas_3_right_atrial_appendage_sq_dev)/(length(I_p_wave_areas_3_right_atrial_appendage)-1));
for i = 1:length(I_p_wave_areas_3_z_scores_right_atrial_appendage)
    I_p_wave_areas_3_z_scores_right_atrial_appendage(i) = (I_p_wave_areas_3(i) - I_p_wave_areas_3_right_atrial_appendage_average)/I_p_wave_areas_3_right_atrial_appendage_SD;
end
I_p_wave_areas_4_z_scores_right_atrial_appendage = zeros(length(I_p_wave_areas_4),1);
I_p_wave_areas_4_right_atrial_appendage_average = sum(I_p_wave_areas_4_right_atrial_appendage)/length(I_p_wave_areas_4_right_atrial_appendage);
I_p_wave_areas_4_right_atrial_appendage_sq_dev = zeros(length(I_p_wave_areas_4_right_atrial_appendage),1);
for i = 1:length(I_p_wave_areas_4_right_atrial_appendage_sq_dev)
    I_p_wave_areas_4_right_atrial_appendage_sq_dev(i) = (I_p_wave_areas_4_right_atrial_appendage(i)-I_p_wave_areas_4_right_atrial_appendage_average)^4;
end
I_p_wave_areas_4_right_atrial_appendage_SD = sqrt(sum(I_p_wave_areas_4_right_atrial_appendage_sq_dev)/(length(I_p_wave_areas_4_right_atrial_appendage)-1));
for i = 1:length(I_p_wave_areas_4_z_scores_right_atrial_appendage)
    I_p_wave_areas_4_z_scores_right_atrial_appendage(i) = (I_p_wave_areas_4(i) - I_p_wave_areas_4_right_atrial_appendage_average)/I_p_wave_areas_4_right_atrial_appendage_SD;
end
I_p_wave_areas_5_z_scores_right_atrial_appendage = zeros(length(I_p_wave_areas_5),1);
I_p_wave_areas_5_right_atrial_appendage_average = sum(I_p_wave_areas_5_right_atrial_appendage)/length(I_p_wave_areas_5_right_atrial_appendage);
I_p_wave_areas_5_right_atrial_appendage_sq_dev = zeros(length(I_p_wave_areas_5_right_atrial_appendage),1);
for i = 1:length(I_p_wave_areas_5_right_atrial_appendage_sq_dev)
    I_p_wave_areas_5_right_atrial_appendage_sq_dev(i) = (I_p_wave_areas_5_right_atrial_appendage(i)-I_p_wave_areas_5_right_atrial_appendage_average)^5;
end
I_p_wave_areas_5_right_atrial_appendage_SD = sqrt(sum(I_p_wave_areas_5_right_atrial_appendage_sq_dev)/(length(I_p_wave_areas_5_right_atrial_appendage)-1));
for i = 1:length(I_p_wave_areas_5_z_scores_right_atrial_appendage)
    I_p_wave_areas_5_z_scores_right_atrial_appendage(i) = (I_p_wave_areas_5(i) - I_p_wave_areas_5_right_atrial_appendage_average)/I_p_wave_areas_5_right_atrial_appendage_SD;
end
I_num_zeros_z_scores_pulmonary_veins = zeros(length(I_number_of_zeros_array),1);
I_num_zeros_pulmonary_veins_average = sum(I_num_zeros_pulmonary_veins)/length(I_num_zeros_pulmonary_veins);
I_num_zeros_pulmonary_veins_sq_dev = zeros(length(I_num_zeros_pulmonary_veins),1);
for i = 1:length(I_num_zeros_pulmonary_veins_sq_dev)
    I_num_zeros_pulmonary_veins_sq_dev(i) = (I_num_zeros_pulmonary_veins(i)-I_num_zeros_pulmonary_veins_average)^2;
end
I_num_zeros_pulmonary_veins_SD = sqrt(sum(I_num_zeros_pulmonary_veins_sq_dev)/(length(I_num_zeros_pulmonary_veins)-1));
for i = 1:length(I_num_zeros_z_scores_pulmonary_veins)
    I_num_zeros_z_scores_pulmonary_veins(i) = (I_number_of_zeros_array(i) - I_num_zeros_pulmonary_veins_average)/I_num_zeros_pulmonary_veins_SD;
end
I_maxima_z_scores_pulmonary_veins = zeros(length(I_maxima),1);
I_maxima_pulmonary_veins_average = sum(I_maxima_pulmonary_veins)/length(I_maxima_pulmonary_veins);
I_maxima_pulmonary_veins_sq_dev = zeros(length(I_maxima_pulmonary_veins),1);
for i = 1:length(I_maxima_pulmonary_veins_sq_dev)
    I_maxima_pulmonary_veins_sq_dev(i) = (I_maxima_pulmonary_veins(i)-I_maxima_pulmonary_veins_average)^2;
end
I_maxima_pulmonary_veins_SD = sqrt(sum(I_maxima_pulmonary_veins_sq_dev)/(length(I_maxima_pulmonary_veins)-1));
for i = 1:length(I_maxima_z_scores_pulmonary_veins)
    I_maxima_z_scores_pulmonary_veins(i) = (I_maxima(i) - I_maxima_pulmonary_veins_average)/I_maxima_pulmonary_veins_SD;
end
I_minima_z_scores_pulmonary_veins = zeros(length(I_minima),1);
I_minima_pulmonary_veins_average = sum(I_minima_pulmonary_veins)/length(I_minima_pulmonary_veins);
I_minima_pulmonary_veins_sq_dev = zeros(length(I_minima_pulmonary_veins),1);
for i = 1:length(I_minima_pulmonary_veins_sq_dev)
    I_minima_pulmonary_veins_sq_dev(i) = (I_minima_pulmonary_veins(i)-I_minima_pulmonary_veins_average)^2;
end
I_minima_pulmonary_veins_SD = sqrt(sum(I_minima_pulmonary_veins_sq_dev)/(length(I_minima_pulmonary_veins)-1));
for i = 1:length(I_minima_z_scores_pulmonary_veins)
    I_minima_z_scores_pulmonary_veins(i) = (I_minima(i) - I_minima_pulmonary_veins_average)/I_minima_pulmonary_veins_SD;
end
I_mdslpes_z_scores_pulmonary_veins = zeros(length(I_mdslpes),1);
I_mdslpes_pulmonary_veins_average = sum(I_mdslpes_pulmonary_veins)/length(I_mdslpes_pulmonary_veins);
I_mdslpes_pulmonary_veins_sq_dev = zeros(length(I_mdslpes_pulmonary_veins),1);
for i = 1:length(I_mdslpes_pulmonary_veins_sq_dev)
    I_mdslpes_pulmonary_veins_sq_dev(i) = (I_mdslpes_pulmonary_veins(i)-I_mdslpes_pulmonary_veins_average)^2;
end
I_mdslpes_pulmonary_veins_SD = sqrt(sum(I_mdslpes_pulmonary_veins_sq_dev)/(length(I_mdslpes_pulmonary_veins)-1));
for i = 1:length(I_mdslpes_z_scores_pulmonary_veins)
    I_mdslpes_z_scores_pulmonary_veins(i) = (I_mdslpes(i) - I_mdslpes_pulmonary_veins_average)/I_mdslpes_pulmonary_veins_SD;
end
I_concavs_z_scores_pulmonary_veins = zeros(length(I_concavs),1);
I_concavs_pulmonary_veins_average = sum(I_concavs_pulmonary_veins)/length(I_concavs_pulmonary_veins);
I_concavs_pulmonary_veins_sq_dev = zeros(length(I_concavs_pulmonary_veins),1);
for i = 1:length(I_concavs_pulmonary_veins_sq_dev)
    I_concavs_pulmonary_veins_sq_dev(i) = (I_concavs_pulmonary_veins(i)-I_concavs_pulmonary_veins_average)^2;
end
I_concavs_pulmonary_veins_SD = sqrt(sum(I_concavs_pulmonary_veins_sq_dev)/(length(I_concavs_pulmonary_veins)-1));
for i = 1:length(I_concavs_z_scores_pulmonary_veins)
    I_concavs_z_scores_pulmonary_veins(i) = (I_concavs(i) - I_concavs_pulmonary_veins_average)/I_concavs_pulmonary_veins_SD;
end
I_sym_inds_z_scores_pulmonary_veins = zeros(length(I_sym_inds),1);
I_sym_inds_pulmonary_veins_average = sum(I_sym_inds_pulmonary_veins)/length(I_sym_inds_pulmonary_veins);
I_sym_inds_pulmonary_veins_sq_dev = zeros(length(I_sym_inds_pulmonary_veins),1);
for i = 1:length(I_sym_inds_pulmonary_veins_sq_dev)
    I_sym_inds_pulmonary_veins_sq_dev(i) = (I_sym_inds_pulmonary_veins(i)-I_sym_inds_pulmonary_veins_average)^2;
end
I_sym_inds_pulmonary_veins_SD = sqrt(sum(I_sym_inds_pulmonary_veins_sq_dev)/(length(I_sym_inds_pulmonary_veins)-1));
for i = 1:length(I_sym_inds_z_scores_pulmonary_veins)
    I_sym_inds_z_scores_pulmonary_veins(i) = (I_sym_inds(i) - I_sym_inds_pulmonary_veins_average)/I_sym_inds_pulmonary_veins_SD;
end
I_p_wave_areas_1_z_scores_pulmonary_veins = zeros(length(I_p_wave_areas_1),1);
I_p_wave_areas_1_pulmonary_veins_average = sum(I_p_wave_areas_1_pulmonary_veins)/length(I_p_wave_areas_1_pulmonary_veins);
I_p_wave_areas_1_pulmonary_veins_sq_dev = zeros(length(I_p_wave_areas_1_pulmonary_veins),1);
for i = 1:length(I_p_wave_areas_1_pulmonary_veins_sq_dev)
    I_p_wave_areas_1_pulmonary_veins_sq_dev(i) = (I_p_wave_areas_1_pulmonary_veins(i)-I_p_wave_areas_1_pulmonary_veins_average)^2;
end
I_p_wave_areas_1_pulmonary_veins_SD = sqrt(sum(I_p_wave_areas_1_pulmonary_veins_sq_dev)/(length(I_p_wave_areas_1_pulmonary_veins)-1));
for i = 1:length(I_p_wave_areas_1_z_scores_pulmonary_veins)
    I_p_wave_areas_1_z_scores_pulmonary_veins(i) = (I_p_wave_areas_1(i) - I_p_wave_areas_1_pulmonary_veins_average)/I_p_wave_areas_1_pulmonary_veins_SD;
end
I_p_wave_areas_2_z_scores_pulmonary_veins = zeros(length(I_p_wave_areas_2),1);
I_p_wave_areas_2_pulmonary_veins_average = sum(I_p_wave_areas_2_pulmonary_veins)/length(I_p_wave_areas_2_pulmonary_veins);
I_p_wave_areas_2_pulmonary_veins_sq_dev = zeros(length(I_p_wave_areas_2_pulmonary_veins),1);
for i = 1:length(I_p_wave_areas_2_pulmonary_veins_sq_dev)
    I_p_wave_areas_2_pulmonary_veins_sq_dev(i) = (I_p_wave_areas_2_pulmonary_veins(i)-I_p_wave_areas_2_pulmonary_veins_average)^2;
end
I_p_wave_areas_2_pulmonary_veins_SD = sqrt(sum(I_p_wave_areas_2_pulmonary_veins_sq_dev)/(length(I_p_wave_areas_2_pulmonary_veins)-1));
for i = 1:length(I_p_wave_areas_2_z_scores_pulmonary_veins)
    I_p_wave_areas_2_z_scores_pulmonary_veins(i) = (I_p_wave_areas_2(i) - I_p_wave_areas_2_pulmonary_veins_average)/I_p_wave_areas_2_pulmonary_veins_SD;
end
I_p_wave_areas_3_z_scores_pulmonary_veins = zeros(length(I_p_wave_areas_3),1);
I_p_wave_areas_3_pulmonary_veins_average = sum(I_p_wave_areas_3_pulmonary_veins)/length(I_p_wave_areas_3_pulmonary_veins);
I_p_wave_areas_3_pulmonary_veins_sq_dev = zeros(length(I_p_wave_areas_3_pulmonary_veins),1);
for i = 1:length(I_p_wave_areas_3_pulmonary_veins_sq_dev)
    I_p_wave_areas_3_pulmonary_veins_sq_dev(i) = (I_p_wave_areas_3_pulmonary_veins(i)-I_p_wave_areas_3_pulmonary_veins_average)^3;
end
I_p_wave_areas_3_pulmonary_veins_SD = sqrt(sum(I_p_wave_areas_3_pulmonary_veins_sq_dev)/(length(I_p_wave_areas_3_pulmonary_veins)-1));
for i = 1:length(I_p_wave_areas_3_z_scores_pulmonary_veins)
    I_p_wave_areas_3_z_scores_pulmonary_veins(i) = (I_p_wave_areas_3(i) - I_p_wave_areas_3_pulmonary_veins_average)/I_p_wave_areas_3_pulmonary_veins_SD;
end
I_p_wave_areas_4_z_scores_pulmonary_veins = zeros(length(I_p_wave_areas_4),1);
I_p_wave_areas_4_pulmonary_veins_average = sum(I_p_wave_areas_4_pulmonary_veins)/length(I_p_wave_areas_4_pulmonary_veins);
I_p_wave_areas_4_pulmonary_veins_sq_dev = zeros(length(I_p_wave_areas_4_pulmonary_veins),1);
for i = 1:length(I_p_wave_areas_4_pulmonary_veins_sq_dev)
    I_p_wave_areas_4_pulmonary_veins_sq_dev(i) = (I_p_wave_areas_4_pulmonary_veins(i)-I_p_wave_areas_4_pulmonary_veins_average)^4;
end
I_p_wave_areas_4_pulmonary_veins_SD = sqrt(sum(I_p_wave_areas_4_pulmonary_veins_sq_dev)/(length(I_p_wave_areas_4_pulmonary_veins)-1));
for i = 1:length(I_p_wave_areas_4_z_scores_pulmonary_veins)
    I_p_wave_areas_4_z_scores_pulmonary_veins(i) = (I_p_wave_areas_4(i) - I_p_wave_areas_4_pulmonary_veins_average)/I_p_wave_areas_4_pulmonary_veins_SD;
end
I_p_wave_areas_5_z_scores_pulmonary_veins = zeros(length(I_p_wave_areas_5),1);
I_p_wave_areas_5_pulmonary_veins_average = sum(I_p_wave_areas_5_pulmonary_veins)/length(I_p_wave_areas_5_pulmonary_veins);
I_p_wave_areas_5_pulmonary_veins_sq_dev = zeros(length(I_p_wave_areas_5_pulmonary_veins),1);
for i = 1:length(I_p_wave_areas_5_pulmonary_veins_sq_dev)
    I_p_wave_areas_5_pulmonary_veins_sq_dev(i) = (I_p_wave_areas_5_pulmonary_veins(i)-I_p_wave_areas_5_pulmonary_veins_average)^5;
end
I_p_wave_areas_5_pulmonary_veins_SD = sqrt(sum(I_p_wave_areas_5_pulmonary_veins_sq_dev)/(length(I_p_wave_areas_5_pulmonary_veins)-1));
for i = 1:length(I_p_wave_areas_5_z_scores_pulmonary_veins)
    I_p_wave_areas_5_z_scores_pulmonary_veins(i) = (I_p_wave_areas_5(i) - I_p_wave_areas_5_pulmonary_veins_average)/I_p_wave_areas_5_pulmonary_veins_SD;
end
I_num_zeros_z_scores_mitral_annulus = zeros(length(I_number_of_zeros_array),1);
I_num_zeros_mitral_annulus_average = sum(I_num_zeros_mitral_annulus)/length(I_num_zeros_mitral_annulus);
I_num_zeros_mitral_annulus_sq_dev = zeros(length(I_num_zeros_mitral_annulus),1);
for i = 1:length(I_num_zeros_mitral_annulus_sq_dev)
    I_num_zeros_mitral_annulus_sq_dev(i) = (I_num_zeros_mitral_annulus(i)-I_num_zeros_mitral_annulus_average)^2;
end
I_num_zeros_mitral_annulus_SD = sqrt(sum(I_num_zeros_mitral_annulus_sq_dev)/(length(I_num_zeros_mitral_annulus)-1));
for i = 1:length(I_num_zeros_z_scores_mitral_annulus)
    I_num_zeros_z_scores_mitral_annulus(i) = (I_number_of_zeros_array(i) - I_num_zeros_mitral_annulus_average)/I_num_zeros_mitral_annulus_SD;
end
I_maxima_z_scores_mitral_annulus = zeros(length(I_maxima),1);
I_maxima_mitral_annulus_average = sum(I_maxima_mitral_annulus)/length(I_maxima_mitral_annulus);
I_maxima_mitral_annulus_sq_dev = zeros(length(I_maxima_mitral_annulus),1);
for i = 1:length(I_maxima_mitral_annulus_sq_dev)
    I_maxima_mitral_annulus_sq_dev(i) = (I_maxima_mitral_annulus(i)-I_maxima_mitral_annulus_average)^2;
end
I_maxima_mitral_annulus_SD = sqrt(sum(I_maxima_mitral_annulus_sq_dev)/(length(I_maxima_mitral_annulus)-1));
for i = 1:length(I_maxima_z_scores_mitral_annulus)
    I_maxima_z_scores_mitral_annulus(i) = (I_maxima(i) - I_maxima_mitral_annulus_average)/I_maxima_mitral_annulus_SD;
end
I_minima_z_scores_mitral_annulus = zeros(length(I_minima),1);
I_minima_mitral_annulus_average = sum(I_minima_mitral_annulus)/length(I_minima_mitral_annulus);
I_minima_mitral_annulus_sq_dev = zeros(length(I_minima_mitral_annulus),1);
for i = 1:length(I_minima_mitral_annulus_sq_dev)
    I_minima_mitral_annulus_sq_dev(i) = (I_minima_mitral_annulus(i)-I_minima_mitral_annulus_average)^2;
end
I_minima_mitral_annulus_SD = sqrt(sum(I_minima_mitral_annulus_sq_dev)/(length(I_minima_mitral_annulus)-1));
for i = 1:length(I_minima_z_scores_mitral_annulus)
    I_minima_z_scores_mitral_annulus(i) = (I_minima(i) - I_minima_mitral_annulus_average)/I_minima_mitral_annulus_SD;
end
I_mdslpes_z_scores_mitral_annulus = zeros(length(I_mdslpes),1);
I_mdslpes_mitral_annulus_average = sum(I_mdslpes_mitral_annulus)/length(I_mdslpes_mitral_annulus);
I_mdslpes_mitral_annulus_sq_dev = zeros(length(I_mdslpes_mitral_annulus),1);
for i = 1:length(I_mdslpes_mitral_annulus_sq_dev)
    I_mdslpes_mitral_annulus_sq_dev(i) = (I_mdslpes_mitral_annulus(i)-I_mdslpes_mitral_annulus_average)^2;
end
I_mdslpes_mitral_annulus_SD = sqrt(sum(I_mdslpes_mitral_annulus_sq_dev)/(length(I_mdslpes_mitral_annulus)-1));
for i = 1:length(I_mdslpes_z_scores_mitral_annulus)
    I_mdslpes_z_scores_mitral_annulus(i) = (I_mdslpes(i) - I_mdslpes_mitral_annulus_average)/I_mdslpes_mitral_annulus_SD;
end
I_concavs_z_scores_mitral_annulus = zeros(length(I_concavs),1);
I_concavs_mitral_annulus_average = sum(I_concavs_mitral_annulus)/length(I_concavs_mitral_annulus);
I_concavs_mitral_annulus_sq_dev = zeros(length(I_concavs_mitral_annulus),1);
for i = 1:length(I_concavs_mitral_annulus_sq_dev)
    I_concavs_mitral_annulus_sq_dev(i) = (I_concavs_mitral_annulus(i)-I_concavs_mitral_annulus_average)^2;
end
I_concavs_mitral_annulus_SD = sqrt(sum(I_concavs_mitral_annulus_sq_dev)/(length(I_concavs_mitral_annulus)-1));
for i = 1:length(I_concavs_z_scores_mitral_annulus)
    I_concavs_z_scores_mitral_annulus(i) = (I_concavs(i) - I_concavs_mitral_annulus_average)/I_concavs_mitral_annulus_SD;
end
I_sym_inds_z_scores_mitral_annulus = zeros(length(I_sym_inds),1);
I_sym_inds_mitral_annulus_average = sum(I_sym_inds_mitral_annulus)/length(I_sym_inds_mitral_annulus);
I_sym_inds_mitral_annulus_sq_dev = zeros(length(I_sym_inds_mitral_annulus),1);
for i = 1:length(I_sym_inds_mitral_annulus_sq_dev)
    I_sym_inds_mitral_annulus_sq_dev(i) = (I_sym_inds_mitral_annulus(i)-I_sym_inds_mitral_annulus_average)^2;
end
I_sym_inds_mitral_annulus_SD = sqrt(sum(I_sym_inds_mitral_annulus_sq_dev)/(length(I_sym_inds_mitral_annulus)-1));
for i = 1:length(I_sym_inds_z_scores_mitral_annulus)
    I_sym_inds_z_scores_mitral_annulus(i) = (I_sym_inds(i) - I_sym_inds_mitral_annulus_average)/I_sym_inds_mitral_annulus_SD;
end
I_p_wave_areas_1_z_scores_mitral_annulus = zeros(length(I_p_wave_areas_1),1);
I_p_wave_areas_1_mitral_annulus_average = sum(I_p_wave_areas_1_mitral_annulus)/length(I_p_wave_areas_1_mitral_annulus);
I_p_wave_areas_1_mitral_annulus_sq_dev = zeros(length(I_p_wave_areas_1_mitral_annulus),1);
for i = 1:length(I_p_wave_areas_1_mitral_annulus_sq_dev)
    I_p_wave_areas_1_mitral_annulus_sq_dev(i) = (I_p_wave_areas_1_mitral_annulus(i)-I_p_wave_areas_1_mitral_annulus_average)^2;
end
I_p_wave_areas_1_mitral_annulus_SD = sqrt(sum(I_p_wave_areas_1_mitral_annulus_sq_dev)/(length(I_p_wave_areas_1_mitral_annulus)-1));
for i = 1:length(I_p_wave_areas_1_z_scores_mitral_annulus)
    I_p_wave_areas_1_z_scores_mitral_annulus(i) = (I_p_wave_areas_1(i) - I_p_wave_areas_1_mitral_annulus_average)/I_p_wave_areas_1_mitral_annulus_SD;
end
I_p_wave_areas_2_z_scores_mitral_annulus = zeros(length(I_p_wave_areas_2),1);
I_p_wave_areas_2_mitral_annulus_average = sum(I_p_wave_areas_2_mitral_annulus)/length(I_p_wave_areas_2_mitral_annulus);
I_p_wave_areas_2_mitral_annulus_sq_dev = zeros(length(I_p_wave_areas_2_mitral_annulus),1);
for i = 1:length(I_p_wave_areas_2_mitral_annulus_sq_dev)
    I_p_wave_areas_2_mitral_annulus_sq_dev(i) = (I_p_wave_areas_2_mitral_annulus(i)-I_p_wave_areas_2_mitral_annulus_average)^2;
end
I_p_wave_areas_2_mitral_annulus_SD = sqrt(sum(I_p_wave_areas_2_mitral_annulus_sq_dev)/(length(I_p_wave_areas_2_mitral_annulus)-1));
for i = 1:length(I_p_wave_areas_2_z_scores_mitral_annulus)
    I_p_wave_areas_2_z_scores_mitral_annulus(i) = (I_p_wave_areas_2(i) - I_p_wave_areas_2_mitral_annulus_average)/I_p_wave_areas_2_mitral_annulus_SD;
end
I_p_wave_areas_3_z_scores_mitral_annulus = zeros(length(I_p_wave_areas_3),1);
I_p_wave_areas_3_mitral_annulus_average = sum(I_p_wave_areas_3_mitral_annulus)/length(I_p_wave_areas_3_mitral_annulus);
I_p_wave_areas_3_mitral_annulus_sq_dev = zeros(length(I_p_wave_areas_3_mitral_annulus),1);
for i = 1:length(I_p_wave_areas_3_mitral_annulus_sq_dev)
    I_p_wave_areas_3_mitral_annulus_sq_dev(i) = (I_p_wave_areas_3_mitral_annulus(i)-I_p_wave_areas_3_mitral_annulus_average)^3;
end
I_p_wave_areas_3_mitral_annulus_SD = sqrt(sum(I_p_wave_areas_3_mitral_annulus_sq_dev)/(length(I_p_wave_areas_3_mitral_annulus)-1));
for i = 1:length(I_p_wave_areas_3_z_scores_mitral_annulus)
    I_p_wave_areas_3_z_scores_mitral_annulus(i) = (I_p_wave_areas_3(i) - I_p_wave_areas_3_mitral_annulus_average)/I_p_wave_areas_3_mitral_annulus_SD;
end
I_p_wave_areas_4_z_scores_mitral_annulus = zeros(length(I_p_wave_areas_4),1);
I_p_wave_areas_4_mitral_annulus_average = sum(I_p_wave_areas_4_mitral_annulus)/length(I_p_wave_areas_4_mitral_annulus);
I_p_wave_areas_4_mitral_annulus_sq_dev = zeros(length(I_p_wave_areas_4_mitral_annulus),1);
for i = 1:length(I_p_wave_areas_4_mitral_annulus_sq_dev)
    I_p_wave_areas_4_mitral_annulus_sq_dev(i) = (I_p_wave_areas_4_mitral_annulus(i)-I_p_wave_areas_4_mitral_annulus_average)^4;
end
I_p_wave_areas_4_mitral_annulus_SD = sqrt(sum(I_p_wave_areas_4_mitral_annulus_sq_dev)/(length(I_p_wave_areas_4_mitral_annulus)-1));
for i = 1:length(I_p_wave_areas_4_z_scores_mitral_annulus)
    I_p_wave_areas_4_z_scores_mitral_annulus(i) = (I_p_wave_areas_4(i) - I_p_wave_areas_4_mitral_annulus_average)/I_p_wave_areas_4_mitral_annulus_SD;
end
I_p_wave_areas_5_z_scores_mitral_annulus = zeros(length(I_p_wave_areas_5),1);
I_p_wave_areas_5_mitral_annulus_average = sum(I_p_wave_areas_5_mitral_annulus)/length(I_p_wave_areas_5_mitral_annulus);
I_p_wave_areas_5_mitral_annulus_sq_dev = zeros(length(I_p_wave_areas_5_mitral_annulus),1);
for i = 1:length(I_p_wave_areas_5_mitral_annulus_sq_dev)
    I_p_wave_areas_5_mitral_annulus_sq_dev(i) = (I_p_wave_areas_5_mitral_annulus(i)-I_p_wave_areas_5_mitral_annulus_average)^5;
end
I_p_wave_areas_5_mitral_annulus_SD = sqrt(sum(I_p_wave_areas_5_mitral_annulus_sq_dev)/(length(I_p_wave_areas_5_mitral_annulus)-1));
for i = 1:length(I_p_wave_areas_5_z_scores_mitral_annulus)
    I_p_wave_areas_5_z_scores_mitral_annulus(i) = (I_p_wave_areas_5(i) - I_p_wave_areas_5_mitral_annulus_average)/I_p_wave_areas_5_mitral_annulus_SD;
end
I_num_zeros_z_scores_CS_body = zeros(length(I_number_of_zeros_array),1);
I_num_zeros_CS_body_average = sum(I_num_zeros_CS_body)/length(I_num_zeros_CS_body);
I_num_zeros_CS_body_sq_dev = zeros(length(I_num_zeros_CS_body),1);
for i = 1:length(I_num_zeros_CS_body_sq_dev)
    I_num_zeros_CS_body_sq_dev(i) = (I_num_zeros_CS_body(i)-I_num_zeros_CS_body_average)^2;
end
I_num_zeros_CS_body_SD = sqrt(sum(I_num_zeros_CS_body_sq_dev)/(length(I_num_zeros_CS_body)-1));
for i = 1:length(I_num_zeros_z_scores_CS_body)
    I_num_zeros_z_scores_CS_body(i) = (I_number_of_zeros_array(i) - I_num_zeros_CS_body_average)/I_num_zeros_CS_body_SD;
end
I_maxima_z_scores_CS_body = zeros(length(I_maxima),1);
I_maxima_CS_body_average = sum(I_maxima_CS_body)/length(I_maxima_CS_body);
I_maxima_CS_body_sq_dev = zeros(length(I_maxima_CS_body),1);
for i = 1:length(I_maxima_CS_body_sq_dev)
    I_maxima_CS_body_sq_dev(i) = (I_maxima_CS_body(i)-I_maxima_CS_body_average)^2;
end
I_maxima_CS_body_SD = sqrt(sum(I_maxima_CS_body_sq_dev)/(length(I_maxima_CS_body)-1));
for i = 1:length(I_maxima_z_scores_CS_body)
    I_maxima_z_scores_CS_body(i) = (I_maxima(i) - I_maxima_CS_body_average)/I_maxima_CS_body_SD;
end
I_minima_z_scores_CS_body = zeros(length(I_minima),1);
I_minima_CS_body_average = sum(I_minima_CS_body)/length(I_minima_CS_body);
I_minima_CS_body_sq_dev = zeros(length(I_minima_CS_body),1);
for i = 1:length(I_minima_CS_body_sq_dev)
    I_minima_CS_body_sq_dev(i) = (I_minima_CS_body(i)-I_minima_CS_body_average)^2;
end
I_minima_CS_body_SD = sqrt(sum(I_minima_CS_body_sq_dev)/(length(I_minima_CS_body)-1));
for i = 1:length(I_minima_z_scores_CS_body)
    I_minima_z_scores_CS_body(i) = (I_minima(i) - I_minima_CS_body_average)/I_minima_CS_body_SD;
end
I_mdslpes_z_scores_CS_body = zeros(length(I_mdslpes),1);
I_mdslpes_CS_body_average = sum(I_mdslpes_CS_body)/length(I_mdslpes_CS_body);
I_mdslpes_CS_body_sq_dev = zeros(length(I_mdslpes_CS_body),1);
for i = 1:length(I_mdslpes_CS_body_sq_dev)
    I_mdslpes_CS_body_sq_dev(i) = (I_mdslpes_CS_body(i)-I_mdslpes_CS_body_average)^2;
end
I_mdslpes_CS_body_SD = sqrt(sum(I_mdslpes_CS_body_sq_dev)/(length(I_mdslpes_CS_body)-1));
for i = 1:length(I_mdslpes_z_scores_CS_body)
    I_mdslpes_z_scores_CS_body(i) = (I_mdslpes(i) - I_mdslpes_CS_body_average)/I_mdslpes_CS_body_SD;
end
I_concavs_z_scores_CS_body = zeros(length(I_concavs),1);
I_concavs_CS_body_average = sum(I_concavs_CS_body)/length(I_concavs_CS_body);
I_concavs_CS_body_sq_dev = zeros(length(I_concavs_CS_body),1);
for i = 1:length(I_concavs_CS_body_sq_dev)
    I_concavs_CS_body_sq_dev(i) = (I_concavs_CS_body(i)-I_concavs_CS_body_average)^2;
end
I_concavs_CS_body_SD = sqrt(sum(I_concavs_CS_body_sq_dev)/(length(I_concavs_CS_body)-1));
for i = 1:length(I_concavs_z_scores_CS_body)
    I_concavs_z_scores_CS_body(i) = (I_concavs(i) - I_concavs_CS_body_average)/I_concavs_CS_body_SD;
end
I_sym_inds_z_scores_CS_body = zeros(length(I_sym_inds),1);
I_sym_inds_CS_body_average = sum(I_sym_inds_CS_body)/length(I_sym_inds_CS_body);
I_sym_inds_CS_body_sq_dev = zeros(length(I_sym_inds_CS_body),1);
for i = 1:length(I_sym_inds_CS_body_sq_dev)
    I_sym_inds_CS_body_sq_dev(i) = (I_sym_inds_CS_body(i)-I_sym_inds_CS_body_average)^2;
end
I_sym_inds_CS_body_SD = sqrt(sum(I_sym_inds_CS_body_sq_dev)/(length(I_sym_inds_CS_body)-1));
for i = 1:length(I_sym_inds_z_scores_CS_body)
    I_sym_inds_z_scores_CS_body(i) = (I_sym_inds(i) - I_sym_inds_CS_body_average)/I_sym_inds_CS_body_SD;
end
I_p_wave_areas_1_z_scores_CS_body = zeros(length(I_p_wave_areas_1),1);
I_p_wave_areas_1_CS_body_average = sum(I_p_wave_areas_1_CS_body)/length(I_p_wave_areas_1_CS_body);
I_p_wave_areas_1_CS_body_sq_dev = zeros(length(I_p_wave_areas_1_CS_body),1);
for i = 1:length(I_p_wave_areas_1_CS_body_sq_dev)
    I_p_wave_areas_1_CS_body_sq_dev(i) = (I_p_wave_areas_1_CS_body(i)-I_p_wave_areas_1_CS_body_average)^2;
end
I_p_wave_areas_1_CS_body_SD = sqrt(sum(I_p_wave_areas_1_CS_body_sq_dev)/(length(I_p_wave_areas_1_CS_body)-1));
for i = 1:length(I_p_wave_areas_1_z_scores_CS_body)
    I_p_wave_areas_1_z_scores_CS_body(i) = (I_p_wave_areas_1(i) - I_p_wave_areas_1_CS_body_average)/I_p_wave_areas_1_CS_body_SD;
end
I_p_wave_areas_2_z_scores_CS_body = zeros(length(I_p_wave_areas_2),1);
I_p_wave_areas_2_CS_body_average = sum(I_p_wave_areas_2_CS_body)/length(I_p_wave_areas_2_CS_body);
I_p_wave_areas_2_CS_body_sq_dev = zeros(length(I_p_wave_areas_2_CS_body),1);
for i = 1:length(I_p_wave_areas_2_CS_body_sq_dev)
    I_p_wave_areas_2_CS_body_sq_dev(i) = (I_p_wave_areas_2_CS_body(i)-I_p_wave_areas_2_CS_body_average)^2;
end
I_p_wave_areas_2_CS_body_SD = sqrt(sum(I_p_wave_areas_2_CS_body_sq_dev)/(length(I_p_wave_areas_2_CS_body)-1));
for i = 1:length(I_p_wave_areas_2_z_scores_CS_body)
    I_p_wave_areas_2_z_scores_CS_body(i) = (I_p_wave_areas_2(i) - I_p_wave_areas_2_CS_body_average)/I_p_wave_areas_2_CS_body_SD;
end
I_p_wave_areas_3_z_scores_CS_body = zeros(length(I_p_wave_areas_3),1);
I_p_wave_areas_3_CS_body_average = sum(I_p_wave_areas_3_CS_body)/length(I_p_wave_areas_3_CS_body);
I_p_wave_areas_3_CS_body_sq_dev = zeros(length(I_p_wave_areas_3_CS_body),1);
for i = 1:length(I_p_wave_areas_3_CS_body_sq_dev)
    I_p_wave_areas_3_CS_body_sq_dev(i) = (I_p_wave_areas_3_CS_body(i)-I_p_wave_areas_3_CS_body_average)^3;
end
I_p_wave_areas_3_CS_body_SD = sqrt(sum(I_p_wave_areas_3_CS_body_sq_dev)/(length(I_p_wave_areas_3_CS_body)-1));
for i = 1:length(I_p_wave_areas_3_z_scores_CS_body)
    I_p_wave_areas_3_z_scores_CS_body(i) = (I_p_wave_areas_3(i) - I_p_wave_areas_3_CS_body_average)/I_p_wave_areas_3_CS_body_SD;
end
I_p_wave_areas_4_z_scores_CS_body = zeros(length(I_p_wave_areas_4),1);
I_p_wave_areas_4_CS_body_average = sum(I_p_wave_areas_4_CS_body)/length(I_p_wave_areas_4_CS_body);
I_p_wave_areas_4_CS_body_sq_dev = zeros(length(I_p_wave_areas_4_CS_body),1);
for i = 1:length(I_p_wave_areas_4_CS_body_sq_dev)
    I_p_wave_areas_4_CS_body_sq_dev(i) = (I_p_wave_areas_4_CS_body(i)-I_p_wave_areas_4_CS_body_average)^4;
end
I_p_wave_areas_4_CS_body_SD = sqrt(sum(I_p_wave_areas_4_CS_body_sq_dev)/(length(I_p_wave_areas_4_CS_body)-1));
for i = 1:length(I_p_wave_areas_4_z_scores_CS_body)
    I_p_wave_areas_4_z_scores_CS_body(i) = (I_p_wave_areas_4(i) - I_p_wave_areas_4_CS_body_average)/I_p_wave_areas_4_CS_body_SD;
end
I_p_wave_areas_5_z_scores_CS_body = zeros(length(I_p_wave_areas_5),1);
I_p_wave_areas_5_CS_body_average = sum(I_p_wave_areas_5_CS_body)/length(I_p_wave_areas_5_CS_body);
I_p_wave_areas_5_CS_body_sq_dev = zeros(length(I_p_wave_areas_5_CS_body),1);
for i = 1:length(I_p_wave_areas_5_CS_body_sq_dev)
    I_p_wave_areas_5_CS_body_sq_dev(i) = (I_p_wave_areas_5_CS_body(i)-I_p_wave_areas_5_CS_body_average)^5;
end
I_p_wave_areas_5_CS_body_SD = sqrt(sum(I_p_wave_areas_5_CS_body_sq_dev)/(length(I_p_wave_areas_5_CS_body)-1));
for i = 1:length(I_p_wave_areas_5_z_scores_CS_body)
    I_p_wave_areas_5_z_scores_CS_body(i) = (I_p_wave_areas_5(i) - I_p_wave_areas_5_CS_body_average)/I_p_wave_areas_5_CS_body_SD;
end
I_num_zeros_z_scores_left_septum = zeros(length(I_number_of_zeros_array),1);
I_num_zeros_left_septum_average = sum(I_num_zeros_left_septum)/length(I_num_zeros_left_septum);
I_num_zeros_left_septum_sq_dev = zeros(length(I_num_zeros_left_septum),1);
for i = 1:length(I_num_zeros_left_septum_sq_dev)
    I_num_zeros_left_septum_sq_dev(i) = (I_num_zeros_left_septum(i)-I_num_zeros_left_septum_average)^2;
end
I_num_zeros_left_septum_SD = sqrt(sum(I_num_zeros_left_septum_sq_dev)/(length(I_num_zeros_left_septum)-1));
for i = 1:length(I_num_zeros_z_scores_left_septum)
    I_num_zeros_z_scores_left_septum(i) = (I_number_of_zeros_array(i) - I_num_zeros_left_septum_average)/I_num_zeros_left_septum_SD;
end
I_maxima_z_scores_left_septum = zeros(length(I_maxima),1);
I_maxima_left_septum_average = sum(I_maxima_left_septum)/length(I_maxima_left_septum);
I_maxima_left_septum_sq_dev = zeros(length(I_maxima_left_septum),1);
for i = 1:length(I_maxima_left_septum_sq_dev)
    I_maxima_left_septum_sq_dev(i) = (I_maxima_left_septum(i)-I_maxima_left_septum_average)^2;
end
I_maxima_left_septum_SD = sqrt(sum(I_maxima_left_septum_sq_dev)/(length(I_maxima_left_septum)-1));
for i = 1:length(I_maxima_z_scores_left_septum)
    I_maxima_z_scores_left_septum(i) = (I_maxima(i) - I_maxima_left_septum_average)/I_maxima_left_septum_SD;
end
I_minima_z_scores_left_septum = zeros(length(I_minima),1);
I_minima_left_septum_average = sum(I_minima_left_septum)/length(I_minima_left_septum);
I_minima_left_septum_sq_dev = zeros(length(I_minima_left_septum),1);
for i = 1:length(I_minima_left_septum_sq_dev)
    I_minima_left_septum_sq_dev(i) = (I_minima_left_septum(i)-I_minima_left_septum_average)^2;
end
I_minima_left_septum_SD = sqrt(sum(I_minima_left_septum_sq_dev)/(length(I_minima_left_septum)-1));
for i = 1:length(I_minima_z_scores_left_septum)
    I_minima_z_scores_left_septum(i) = (I_minima(i) - I_minima_left_septum_average)/I_minima_left_septum_SD;
end
I_mdslpes_z_scores_left_septum = zeros(length(I_mdslpes),1);
I_mdslpes_left_septum_average = sum(I_mdslpes_left_septum)/length(I_mdslpes_left_septum);
I_mdslpes_left_septum_sq_dev = zeros(length(I_mdslpes_left_septum),1);
for i = 1:length(I_mdslpes_left_septum_sq_dev)
    I_mdslpes_left_septum_sq_dev(i) = (I_mdslpes_left_septum(i)-I_mdslpes_left_septum_average)^2;
end
I_mdslpes_left_septum_SD = sqrt(sum(I_mdslpes_left_septum_sq_dev)/(length(I_mdslpes_left_septum)-1));
for i = 1:length(I_mdslpes_z_scores_left_septum)
    I_mdslpes_z_scores_left_septum(i) = (I_mdslpes(i) - I_mdslpes_left_septum_average)/I_mdslpes_left_septum_SD;
end
I_concavs_z_scores_left_septum = zeros(length(I_concavs),1);
I_concavs_left_septum_average = sum(I_concavs_left_septum)/length(I_concavs_left_septum);
I_concavs_left_septum_sq_dev = zeros(length(I_concavs_left_septum),1);
for i = 1:length(I_concavs_left_septum_sq_dev)
    I_concavs_left_septum_sq_dev(i) = (I_concavs_left_septum(i)-I_concavs_left_septum_average)^2;
end
I_concavs_left_septum_SD = sqrt(sum(I_concavs_left_septum_sq_dev)/(length(I_concavs_left_septum)-1));
for i = 1:length(I_concavs_z_scores_left_septum)
    I_concavs_z_scores_left_septum(i) = (I_concavs(i) - I_concavs_left_septum_average)/I_concavs_left_septum_SD;
end
I_sym_inds_z_scores_left_septum = zeros(length(I_sym_inds),1);
I_sym_inds_left_septum_average = sum(I_sym_inds_left_septum)/length(I_sym_inds_left_septum);
I_sym_inds_left_septum_sq_dev = zeros(length(I_sym_inds_left_septum),1);
for i = 1:length(I_sym_inds_left_septum_sq_dev)
    I_sym_inds_left_septum_sq_dev(i) = (I_sym_inds_left_septum(i)-I_sym_inds_left_septum_average)^2;
end
I_sym_inds_left_septum_SD = sqrt(sum(I_sym_inds_left_septum_sq_dev)/(length(I_sym_inds_left_septum)-1));
for i = 1:length(I_sym_inds_z_scores_left_septum)
    I_sym_inds_z_scores_left_septum(i) = (I_sym_inds(i) - I_sym_inds_left_septum_average)/I_sym_inds_left_septum_SD;
end
I_p_wave_areas_1_z_scores_left_septum = zeros(length(I_p_wave_areas_1),1);
I_p_wave_areas_1_left_septum_average = sum(I_p_wave_areas_1_left_septum)/length(I_p_wave_areas_1_left_septum);
I_p_wave_areas_1_left_septum_sq_dev = zeros(length(I_p_wave_areas_1_left_septum),1);
for i = 1:length(I_p_wave_areas_1_left_septum_sq_dev)
    I_p_wave_areas_1_left_septum_sq_dev(i) = (I_p_wave_areas_1_left_septum(i)-I_p_wave_areas_1_left_septum_average)^2;
end
I_p_wave_areas_1_left_septum_SD = sqrt(sum(I_p_wave_areas_1_left_septum_sq_dev)/(length(I_p_wave_areas_1_left_septum)-1));
for i = 1:length(I_p_wave_areas_1_z_scores_left_septum)
    I_p_wave_areas_1_z_scores_left_septum(i) = (I_p_wave_areas_1(i) - I_p_wave_areas_1_left_septum_average)/I_p_wave_areas_1_left_septum_SD;
end
I_p_wave_areas_2_z_scores_left_septum = zeros(length(I_p_wave_areas_2),1);
I_p_wave_areas_2_left_septum_average = sum(I_p_wave_areas_2_left_septum)/length(I_p_wave_areas_2_left_septum);
I_p_wave_areas_2_left_septum_sq_dev = zeros(length(I_p_wave_areas_2_left_septum),1);
for i = 1:length(I_p_wave_areas_2_left_septum_sq_dev)
    I_p_wave_areas_2_left_septum_sq_dev(i) = (I_p_wave_areas_2_left_septum(i)-I_p_wave_areas_2_left_septum_average)^2;
end
I_p_wave_areas_2_left_septum_SD = sqrt(sum(I_p_wave_areas_2_left_septum_sq_dev)/(length(I_p_wave_areas_2_left_septum)-1));
for i = 1:length(I_p_wave_areas_2_z_scores_left_septum)
    I_p_wave_areas_2_z_scores_left_septum(i) = (I_p_wave_areas_2(i) - I_p_wave_areas_2_left_septum_average)/I_p_wave_areas_2_left_septum_SD;
end
I_p_wave_areas_3_z_scores_left_septum = zeros(length(I_p_wave_areas_3),1);
I_p_wave_areas_3_left_septum_average = sum(I_p_wave_areas_3_left_septum)/length(I_p_wave_areas_3_left_septum);
I_p_wave_areas_3_left_septum_sq_dev = zeros(length(I_p_wave_areas_3_left_septum),1);
for i = 1:length(I_p_wave_areas_3_left_septum_sq_dev)
    I_p_wave_areas_3_left_septum_sq_dev(i) = (I_p_wave_areas_3_left_septum(i)-I_p_wave_areas_3_left_septum_average)^3;
end
I_p_wave_areas_3_left_septum_SD = sqrt(sum(I_p_wave_areas_3_left_septum_sq_dev)/(length(I_p_wave_areas_3_left_septum)-1));
for i = 1:length(I_p_wave_areas_3_z_scores_left_septum)
    I_p_wave_areas_3_z_scores_left_septum(i) = (I_p_wave_areas_3(i) - I_p_wave_areas_3_left_septum_average)/I_p_wave_areas_3_left_septum_SD;
end
I_p_wave_areas_4_z_scores_left_septum = zeros(length(I_p_wave_areas_4),1);
I_p_wave_areas_4_left_septum_average = sum(I_p_wave_areas_4_left_septum)/length(I_p_wave_areas_4_left_septum);
I_p_wave_areas_4_left_septum_sq_dev = zeros(length(I_p_wave_areas_4_left_septum),1);
for i = 1:length(I_p_wave_areas_4_left_septum_sq_dev)
    I_p_wave_areas_4_left_septum_sq_dev(i) = (I_p_wave_areas_4_left_septum(i)-I_p_wave_areas_4_left_septum_average)^4;
end
I_p_wave_areas_4_left_septum_SD = sqrt(sum(I_p_wave_areas_4_left_septum_sq_dev)/(length(I_p_wave_areas_4_left_septum)-1));
for i = 1:length(I_p_wave_areas_4_z_scores_left_septum)
    I_p_wave_areas_4_z_scores_left_septum(i) = (I_p_wave_areas_4(i) - I_p_wave_areas_4_left_septum_average)/I_p_wave_areas_4_left_septum_SD;
end
I_p_wave_areas_5_z_scores_left_septum = zeros(length(I_p_wave_areas_5),1);
I_p_wave_areas_5_left_septum_average = sum(I_p_wave_areas_5_left_septum)/length(I_p_wave_areas_5_left_septum);
I_p_wave_areas_5_left_septum_sq_dev = zeros(length(I_p_wave_areas_5_left_septum),1);
for i = 1:length(I_p_wave_areas_5_left_septum_sq_dev)
    I_p_wave_areas_5_left_septum_sq_dev(i) = (I_p_wave_areas_5_left_septum(i)-I_p_wave_areas_5_left_septum_average)^5;
end
I_p_wave_areas_5_left_septum_SD = sqrt(sum(I_p_wave_areas_5_left_septum_sq_dev)/(length(I_p_wave_areas_5_left_septum)-1));
for i = 1:length(I_p_wave_areas_5_z_scores_left_septum)
    I_p_wave_areas_5_z_scores_left_septum(i) = (I_p_wave_areas_5(i) - I_p_wave_areas_5_left_septum_average)/I_p_wave_areas_5_left_septum_SD;
end
I_num_zeros_z_scores_left_atrial_appendage = zeros(length(I_number_of_zeros_array),1);
I_num_zeros_left_atrial_appendage_average = sum(I_num_zeros_left_atrial_appendage)/length(I_num_zeros_left_atrial_appendage);
I_num_zeros_left_atrial_appendage_sq_dev = zeros(length(I_num_zeros_left_atrial_appendage),1);
for i = 1:length(I_num_zeros_left_atrial_appendage_sq_dev)
    I_num_zeros_left_atrial_appendage_sq_dev(i) = (I_num_zeros_left_atrial_appendage(i)-I_num_zeros_left_atrial_appendage_average)^2;
end
I_num_zeros_left_atrial_appendage_SD = sqrt(sum(I_num_zeros_left_atrial_appendage_sq_dev)/(length(I_num_zeros_left_atrial_appendage)-1));
for i = 1:length(I_num_zeros_z_scores_left_atrial_appendage)
    I_num_zeros_z_scores_left_atrial_appendage(i) = (I_number_of_zeros_array(i) - I_num_zeros_left_atrial_appendage_average)/I_num_zeros_left_atrial_appendage_SD;
end
I_maxima_z_scores_left_atrial_appendage = zeros(length(I_maxima),1);
I_maxima_left_atrial_appendage_average = sum(I_maxima_left_atrial_appendage)/length(I_maxima_left_atrial_appendage);
I_maxima_left_atrial_appendage_sq_dev = zeros(length(I_maxima_left_atrial_appendage),1);
for i = 1:length(I_maxima_left_atrial_appendage_sq_dev)
    I_maxima_left_atrial_appendage_sq_dev(i) = (I_maxima_left_atrial_appendage(i)-I_maxima_left_atrial_appendage_average)^2;
end
I_maxima_left_atrial_appendage_SD = sqrt(sum(I_maxima_left_atrial_appendage_sq_dev)/(length(I_maxima_left_atrial_appendage)-1));
for i = 1:length(I_maxima_z_scores_left_atrial_appendage)
    I_maxima_z_scores_left_atrial_appendage(i) = (I_maxima(i) - I_maxima_left_atrial_appendage_average)/I_maxima_left_atrial_appendage_SD;
end
I_minima_z_scores_left_atrial_appendage = zeros(length(I_minima),1);
I_minima_left_atrial_appendage_average = sum(I_minima_left_atrial_appendage)/length(I_minima_left_atrial_appendage);
I_minima_left_atrial_appendage_sq_dev = zeros(length(I_minima_left_atrial_appendage),1);
for i = 1:length(I_minima_left_atrial_appendage_sq_dev)
    I_minima_left_atrial_appendage_sq_dev(i) = (I_minima_left_atrial_appendage(i)-I_minima_left_atrial_appendage_average)^2;
end
I_minima_left_atrial_appendage_SD = sqrt(sum(I_minima_left_atrial_appendage_sq_dev)/(length(I_minima_left_atrial_appendage)-1));
for i = 1:length(I_minima_z_scores_left_atrial_appendage)
    I_minima_z_scores_left_atrial_appendage(i) = (I_minima(i) - I_minima_left_atrial_appendage_average)/I_minima_left_atrial_appendage_SD;
end
I_mdslpes_z_scores_left_atrial_appendage = zeros(length(I_mdslpes),1);
I_mdslpes_left_atrial_appendage_average = sum(I_mdslpes_left_atrial_appendage)/length(I_mdslpes_left_atrial_appendage);
I_mdslpes_left_atrial_appendage_sq_dev = zeros(length(I_mdslpes_left_atrial_appendage),1);
for i = 1:length(I_mdslpes_left_atrial_appendage_sq_dev)
    I_mdslpes_left_atrial_appendage_sq_dev(i) = (I_mdslpes_left_atrial_appendage(i)-I_mdslpes_left_atrial_appendage_average)^2;
end
I_mdslpes_left_atrial_appendage_SD = sqrt(sum(I_mdslpes_left_atrial_appendage_sq_dev)/(length(I_mdslpes_left_atrial_appendage)-1));
for i = 1:length(I_mdslpes_z_scores_left_atrial_appendage)
    I_mdslpes_z_scores_left_atrial_appendage(i) = (I_mdslpes(i) - I_mdslpes_left_atrial_appendage_average)/I_mdslpes_left_atrial_appendage_SD;
end
I_concavs_z_scores_left_atrial_appendage = zeros(length(I_concavs),1);
I_concavs_left_atrial_appendage_average = sum(I_concavs_left_atrial_appendage)/length(I_concavs_left_atrial_appendage);
I_concavs_left_atrial_appendage_sq_dev = zeros(length(I_concavs_left_atrial_appendage),1);
for i = 1:length(I_concavs_left_atrial_appendage_sq_dev)
    I_concavs_left_atrial_appendage_sq_dev(i) = (I_concavs_left_atrial_appendage(i)-I_concavs_left_atrial_appendage_average)^2;
end
I_concavs_left_atrial_appendage_SD = sqrt(sum(I_concavs_left_atrial_appendage_sq_dev)/(length(I_concavs_left_atrial_appendage)-1));
for i = 1:length(I_concavs_z_scores_left_atrial_appendage)
    I_concavs_z_scores_left_atrial_appendage(i) = (I_concavs(i) - I_concavs_left_atrial_appendage_average)/I_concavs_left_atrial_appendage_SD;
end
I_sym_inds_z_scores_left_atrial_appendage = zeros(length(I_sym_inds),1);
I_sym_inds_left_atrial_appendage_average = sum(I_sym_inds_left_atrial_appendage)/length(I_sym_inds_left_atrial_appendage);
I_sym_inds_left_atrial_appendage_sq_dev = zeros(length(I_sym_inds_left_atrial_appendage),1);
for i = 1:length(I_sym_inds_left_atrial_appendage_sq_dev)
    I_sym_inds_left_atrial_appendage_sq_dev(i) = (I_sym_inds_left_atrial_appendage(i)-I_sym_inds_left_atrial_appendage_average)^2;
end
I_sym_inds_left_atrial_appendage_SD = sqrt(sum(I_sym_inds_left_atrial_appendage_sq_dev)/(length(I_sym_inds_left_atrial_appendage)-1));
for i = 1:length(I_sym_inds_z_scores_left_atrial_appendage)
    I_sym_inds_z_scores_left_atrial_appendage(i) = (I_sym_inds(i) - I_sym_inds_left_atrial_appendage_average)/I_sym_inds_left_atrial_appendage_SD;
end
I_p_wave_areas_1_z_scores_left_atrial_appendage = zeros(length(I_p_wave_areas_1),1);
I_p_wave_areas_1_left_atrial_appendage_average = sum(I_p_wave_areas_1_left_atrial_appendage)/length(I_p_wave_areas_1_left_atrial_appendage);
I_p_wave_areas_1_left_atrial_appendage_sq_dev = zeros(length(I_p_wave_areas_1_left_atrial_appendage),1);
for i = 1:length(I_p_wave_areas_1_left_atrial_appendage_sq_dev)
    I_p_wave_areas_1_left_atrial_appendage_sq_dev(i) = (I_p_wave_areas_1_left_atrial_appendage(i)-I_p_wave_areas_1_left_atrial_appendage_average)^2;
end
I_p_wave_areas_1_left_atrial_appendage_SD = sqrt(sum(I_p_wave_areas_1_left_atrial_appendage_sq_dev)/(length(I_p_wave_areas_1_left_atrial_appendage)-1));
for i = 1:length(I_p_wave_areas_1_z_scores_left_atrial_appendage)
    I_p_wave_areas_1_z_scores_left_atrial_appendage(i) = (I_p_wave_areas_1(i) - I_p_wave_areas_1_left_atrial_appendage_average)/I_p_wave_areas_1_left_atrial_appendage_SD;
end
I_p_wave_areas_2_z_scores_left_atrial_appendage = zeros(length(I_p_wave_areas_2),1);
I_p_wave_areas_2_left_atrial_appendage_average = sum(I_p_wave_areas_2_left_atrial_appendage)/length(I_p_wave_areas_2_left_atrial_appendage);
I_p_wave_areas_2_left_atrial_appendage_sq_dev = zeros(length(I_p_wave_areas_2_left_atrial_appendage),1);
for i = 1:length(I_p_wave_areas_2_left_atrial_appendage_sq_dev)
    I_p_wave_areas_2_left_atrial_appendage_sq_dev(i) = (I_p_wave_areas_2_left_atrial_appendage(i)-I_p_wave_areas_2_left_atrial_appendage_average)^2;
end
I_p_wave_areas_2_left_atrial_appendage_SD = sqrt(sum(I_p_wave_areas_2_left_atrial_appendage_sq_dev)/(length(I_p_wave_areas_2_left_atrial_appendage)-1));
for i = 1:length(I_p_wave_areas_2_z_scores_left_atrial_appendage)
    I_p_wave_areas_2_z_scores_left_atrial_appendage(i) = (I_p_wave_areas_2(i) - I_p_wave_areas_2_left_atrial_appendage_average)/I_p_wave_areas_2_left_atrial_appendage_SD;
end
I_p_wave_areas_3_z_scores_left_atrial_appendage = zeros(length(I_p_wave_areas_3),1);
I_p_wave_areas_3_left_atrial_appendage_average = sum(I_p_wave_areas_3_left_atrial_appendage)/length(I_p_wave_areas_3_left_atrial_appendage);
I_p_wave_areas_3_left_atrial_appendage_sq_dev = zeros(length(I_p_wave_areas_3_left_atrial_appendage),1);
for i = 1:length(I_p_wave_areas_3_left_atrial_appendage_sq_dev)
    I_p_wave_areas_3_left_atrial_appendage_sq_dev(i) = (I_p_wave_areas_3_left_atrial_appendage(i)-I_p_wave_areas_3_left_atrial_appendage_average)^3;
end
I_p_wave_areas_3_left_atrial_appendage_SD = sqrt(sum(I_p_wave_areas_3_left_atrial_appendage_sq_dev)/(length(I_p_wave_areas_3_left_atrial_appendage)-1));
for i = 1:length(I_p_wave_areas_3_z_scores_left_atrial_appendage)
    I_p_wave_areas_3_z_scores_left_atrial_appendage(i) = (I_p_wave_areas_3(i) - I_p_wave_areas_3_left_atrial_appendage_average)/I_p_wave_areas_3_left_atrial_appendage_SD;
end
I_p_wave_areas_4_z_scores_left_atrial_appendage = zeros(length(I_p_wave_areas_4),1);
I_p_wave_areas_4_left_atrial_appendage_average = sum(I_p_wave_areas_4_left_atrial_appendage)/length(I_p_wave_areas_4_left_atrial_appendage);
I_p_wave_areas_4_left_atrial_appendage_sq_dev = zeros(length(I_p_wave_areas_4_left_atrial_appendage),1);
for i = 1:length(I_p_wave_areas_4_left_atrial_appendage_sq_dev)
    I_p_wave_areas_4_left_atrial_appendage_sq_dev(i) = (I_p_wave_areas_4_left_atrial_appendage(i)-I_p_wave_areas_4_left_atrial_appendage_average)^4;
end
I_p_wave_areas_4_left_atrial_appendage_SD = sqrt(sum(I_p_wave_areas_4_left_atrial_appendage_sq_dev)/(length(I_p_wave_areas_4_left_atrial_appendage)-1));
for i = 1:length(I_p_wave_areas_4_z_scores_left_atrial_appendage)
    I_p_wave_areas_4_z_scores_left_atrial_appendage(i) = (I_p_wave_areas_4(i) - I_p_wave_areas_4_left_atrial_appendage_average)/I_p_wave_areas_4_left_atrial_appendage_SD;
end
I_p_wave_areas_5_z_scores_left_atrial_appendage = zeros(length(I_p_wave_areas_5),1);
I_p_wave_areas_5_left_atrial_appendage_average = sum(I_p_wave_areas_5_left_atrial_appendage)/length(I_p_wave_areas_5_left_atrial_appendage);
I_p_wave_areas_5_left_atrial_appendage_sq_dev = zeros(length(I_p_wave_areas_5_left_atrial_appendage),1);
for i = 1:length(I_p_wave_areas_5_left_atrial_appendage_sq_dev)
    I_p_wave_areas_5_left_atrial_appendage_sq_dev(i) = (I_p_wave_areas_5_left_atrial_appendage(i)-I_p_wave_areas_5_left_atrial_appendage_average)^5;
end
I_p_wave_areas_5_left_atrial_appendage_SD = sqrt(sum(I_p_wave_areas_5_left_atrial_appendage_sq_dev)/(length(I_p_wave_areas_5_left_atrial_appendage)-1));
for i = 1:length(I_p_wave_areas_5_z_scores_left_atrial_appendage)
    I_p_wave_areas_5_z_scores_left_atrial_appendage(i) = (I_p_wave_areas_5(i) - I_p_wave_areas_5_left_atrial_appendage_average)/I_p_wave_areas_5_left_atrial_appendage_SD;
end
number_of_parameters = 11;
number_of_regions = 13;
isNormal = zeros(number_of_parameters, number_of_regions);
isSpike = zeros(number_of_parameters, number_of_regions);
%isNormal/isSpike has 11 rows (1 for each parameter) and 13 columns (1 for each
%possible region of origin)
I_num_zeros_SA_node_sorted = sort(I_num_zeros_SA_node);
I_num_zeros_SA_node_dev_avg_min_1SD = zeros(length(I_num_zeros_SA_node_sorted),1);
for i = 1:length(I_num_zeros_SA_node_sorted)
    I_num_zeros_SA_node_dev_avg_min_1SD(i) = I_num_zeros_SA_node_sorted(i) - (I_num_zeros_SA_node_average-I_num_zeros_SA_node_SD);
end
I_num_zeros_SA_node_dev_avg_min_1SD = abs(I_num_zeros_SA_node_dev_avg_min_1SD);
for i = 1:length(I_num_zeros_SA_node_sorted)
    if I_num_zeros_SA_node_dev_avg_min_1SD(i) == min(I_num_zeros_SA_node_dev_avg_min_1SD)
        I_num_zeros_SA_node_dev_avg_min_1SD_index = i;
    end
end
I_num_zeros_SA_node_dev_avg_pls_1SD = zeros(length(I_num_zeros_SA_node_sorted),1);
for i = 1:length(I_num_zeros_SA_node_sorted)
    I_num_zeros_SA_node_dev_avg_pls_1SD(i) = I_num_zeros_SA_node_sorted(i) - (I_num_zeros_SA_node_average+I_num_zeros_SA_node_SD);
end
I_num_zeros_SA_node_dev_avg_pls_1SD = abs(I_num_zeros_SA_node_dev_avg_pls_1SD);
for i = 1:length(I_num_zeros_SA_node_sorted)
    if I_num_zeros_SA_node_dev_avg_pls_1SD(i) == min(I_num_zeros_SA_node_dev_avg_pls_1SD)
        I_num_zeros_SA_node_dev_avg_pls_1SD_index = i;
    end
end
I_num_zeros_SA_node_dev_avg_min_2SD = zeros(length(I_num_zeros_SA_node_sorted),1);
for i = 1:length(I_num_zeros_SA_node_sorted)
    I_num_zeros_SA_node_dev_avg_min_2SD(i) = I_num_zeros_SA_node_sorted(i) - (I_num_zeros_SA_node_average-2*I_num_zeros_SA_node_SD);
end
I_num_zeros_SA_node_dev_avg_min_2SD = abs(I_num_zeros_SA_node_dev_avg_min_2SD);
for i = 1:length(I_num_zeros_SA_node_sorted)
    if I_num_zeros_SA_node_dev_avg_min_2SD(i) == min(I_num_zeros_SA_node_dev_avg_min_2SD)
        I_num_zeros_dev_avg_min_2SD_index = i;
    end
end
I_num_zeros_SA_node_dev_avg_pls_2SD = zeros(length(I_num_zeros_SA_node_sorted),1);
for i = 1:length(I_num_zeros_SA_node_sorted)
    I_num_zeros_SA_node_dev_avg_pls_2SD(i) = I_num_zeros_SA_node_sorted(i) - (I_num_zeros_SA_node_average+2*I_num_zeros_SA_node_SD);
end
I_num_zeros_SA_node_dev_avg_pls_2SD = abs(I_num_zeros_SA_node_dev_avg_pls_2SD);
for i = 1:length(I_num_zeros_SA_node_sorted)
    if I_num_zeros_SA_node_dev_avg_pls_2SD(i) == min(I_num_zeros_SA_node_dev_avg_pls_2SD)
        I_num_zeros_SA_node_dev_avg_pls_2SD_index = i;
    end
end
I_num_zeros_SA_node_dev_avg_min_3SD = zeros(length(I_num_zeros_SA_node_sorted),1);
for i = 1:length(I_num_zeros_SA_node_sorted)
    I_num_zeros_SA_node_dev_avg_min_3SD(i) = I_num_zeros_SA_node_sorted(i) - (I_num_zeros_SA_node_average-3*I_num_zeros_SA_node_SD);
end
I_num_zeros_SA_node_dev_avg_min_3SD = abs(I_num_zeros_SA_node_dev_avg_min_3SD);
for i = 1:length(I_num_zeros_SA_node_sorted)
    if I_num_zeros_SA_node_dev_avg_min_3SD(i) == min(I_num_zeros_SA_node_dev_avg_min_3SD)
        I_num_zeros_SA_node_dev_avg_min_3SD_index = i;
    end
end
I_num_zeros_SA_node_dev_avg_pls_3SD = zeros(length(I_num_zeros_SA_node_sorted),1);
for i = 1:length(I_num_zeros_SA_node_sorted)
    I_num_zeros_SA_node_dev_avg_pls_3SD(i) = I_num_zeros_SA_node_sorted(i) - (I_num_zeros_SA_node_average+3*I_num_zeros_SA_node_SD);
end
I_num_zeros_SA_node_dev_avg_pls_3SD = abs(I_num_zeros_SA_node_dev_avg_pls_3SD);
for i = 1:length(I_num_zeros_SA_node_sorted)
    if I_num_zeros_SA_node_dev_avg_pls_3SD(i) == min(I_num_zeros_SA_node_dev_avg_pls_3SD)
        I_num_zeros_SA_node_dev_avg_pls_3SD_index = i;
    end
end
I_num_zeros_SA_node_tot_number_of_indices = length(I_num_zeros_SA_node_sorted);
%This is the approximation of the tot "area under curve." 
I_num_zeros_SA_node_quotient_within_1SD = (I_num_zeros_SA_node_dev_avg_pls_1SD_index - I_num_zeros_SA_nodedev_avg_min_1SD_index)/I_num_zeros_SA_node_tot_number_of_indices;
I_num_zeros_SA_node_quotient_within_2SD = (I_num_zeros_SA_node_dev_avg_pls_2SD_index - I_num_zeros_dev_avg_min_2SD_index)/I_num_zeros_SA_node_tot_number_of_indices;
I_num_zeros_SA_node_quotient_within_3SD = (I_num_zeros_SA_node_dev_avg_pls_3SD_index - I_num_zeros_SA_node_dev_avg_min_3SD_index)/I_num_zeros_SA_node_tot_number_of_indices;
if ((I_num_zeros_SA_node_quotient_within_1SD > 0.66 && I_num_zeros_SA_node_quotient_within_1SD < 0.70) && (I_num_zeros_SA_node_quotient_within_2SD > 0.93 && I_num_zeros_SA_node_quotient_within_2SD < 0.97) && (I_num_zeros_SA_node_quotient_within_3SD > 0.98 && I_num_zeros_SA_node_quotient_within_3SD < 1)) 
    isNormal(1,1) = 1;
end
if I_num_zeros_SA_node_quotient_within_1SD == 0
    isSpike(1,1) = 1;
end
I_num_zeros_crista_terminalis_sorted = sort(I_num_zeros_crista_terminalis);
I_num_zeros_crista_terminalis_dev_avg_min_1SD = zeros(length(I_num_zeros_crista_terminalis_sorted),1);
for i = 1:length(I_num_zeros_crista_terminalis_sorted)
    I_num_zeros_crista_terminalis_dev_avg_min_1SD(i) = I_num_zeros_crista_terminalis_sorted(i) - (I_num_zeros_crista_terminalis_average-I_num_zeros_crista_terminalis_SD);
end
I_num_zeros_crista_terminalis_dev_avg_min_1SD = abs(I_num_zeros_crista_terminalis_dev_avg_min_1SD);
for i = 1:length(I_num_zeros_crista_terminalis_sorted)
    if I_num_zeros_crista_terminalis_dev_avg_min_1SD(i) == min(I_num_zeros_crista_terminalis_dev_avg_min_1SD)
        I_num_zeros_crista_terminalis_dev_avg_min_1SD_index = i;
    end
end
I_num_zeros_crista_terminalis_dev_avg_pls_1SD = zeros(length(I_num_zeros_crista_terminalis_sorted),1);
for i = 1:length(I_num_zeros_crista_terminalis_sorted)
    I_num_zeros_crista_terminalis_dev_avg_pls_1SD(i) = I_num_zeros_crista_terminalis_sorted(i) - (I_num_zeros_crista_terminalis_average+I_num_zeros_crista_terminalis_SD);
end
I_num_zeros_crista_terminalis_dev_avg_pls_1SD = abs(I_num_zeros_crista_terminalis_dev_avg_pls_1SD);
for i = 1:length(I_num_zeros_crista_terminalis_sorted)
    if I_num_zeros_crista_terminalis_dev_avg_pls_1SD(i) == min(I_num_zeros_crista_terminalis_dev_avg_pls_1SD)
        I_num_zeros_crista_terminalis_dev_avg_pls_1SD_index = i;
    end
end
I_num_zeros_crista_terminalis_dev_avg_min_2SD = zeros(length(I_num_zeros_crista_terminalis_sorted),1);
for i = 1:length(I_num_zeros_crista_terminalis_sorted)
    I_num_zeros_crista_terminalis_dev_avg_min_2SD(i) = I_num_zeros_crista_terminalis_sorted(i) - (I_num_zeros_crista_terminalis_average-2*I_num_zeros_crista_terminalis_SD);
end
I_num_zeros_crista_terminalis_dev_avg_min_2SD = abs(I_num_zeros_crista_terminalis_dev_avg_min_2SD);
for i = 1:length(I_num_zeros_crista_terminalis_sorted)
    if I_num_zeros_crista_terminalis_dev_avg_min_2SD(i) == min(I_num_zeros_crista_terminalis_dev_avg_min_2SD)
        I_num_zeros_dev_avg_min_2SD_index = i;
    end
end
I_num_zeros_crista_terminalis_dev_avg_pls_2SD = zeros(length(I_num_zeros_crista_terminalis_sorted),1);
for i = 1:length(I_num_zeros_crista_terminalis_sorted)
    I_num_zeros_crista_terminalis_dev_avg_pls_2SD(i) = I_num_zeros_crista_terminalis_sorted(i) - (I_num_zeros_crista_terminalis_average+2*I_num_zeros_crista_terminalis_SD);
end
I_num_zeros_crista_terminalis_dev_avg_pls_2SD = abs(I_num_zeros_crista_terminalis_dev_avg_pls_2SD);
for i = 1:length(I_num_zeros_crista_terminalis_sorted)
    if I_num_zeros_crista_terminalis_dev_avg_pls_2SD(i) == min(I_num_zeros_crista_terminalis_dev_avg_pls_2SD)
        I_num_zeros_crista_terminalis_dev_avg_pls_2SD_index = i;
    end
end
I_num_zeros_crista_terminalis_dev_avg_min_3SD = zeros(length(I_num_zeros_crista_terminalis_sorted),1);
for i = 1:length(I_num_zeros_crista_terminalis_sorted)
    I_num_zeros_crista_terminalis_dev_avg_min_3SD(i) = I_num_zeros_crista_terminalis_sorted(i) - (I_num_zeros_crista_terminalis_average-3*I_num_zeros_crista_terminalis_SD);
end
I_num_zeros_crista_terminalis_dev_avg_min_3SD = abs(I_num_zeros_crista_terminalis_dev_avg_min_3SD);
for i = 1:length(I_num_zeros_crista_terminalis_sorted)
    if I_num_zeros_crista_terminalis_dev_avg_min_3SD(i) == min(I_num_zeros_crista_terminalis_dev_avg_min_3SD)
        I_num_zeros_crista_terminalis_dev_avg_min_3SD_index = i;
    end
end
I_num_zeros_crista_terminalis_dev_avg_pls_3SD = zeros(length(I_num_zeros_crista_terminalis_sorted),1);
for i = 1:length(I_num_zeros_crista_terminalis_sorted)
    I_num_zeros_crista_terminalis_dev_avg_pls_3SD(i) = I_num_zeros_crista_terminalis_sorted(i) - (I_num_zeros_crista_terminalis_average+3*I_num_zeros_crista_terminalis_SD);
end
I_num_zeros_crista_terminalis_dev_avg_pls_3SD = abs(I_num_zeros_crista_terminalis_dev_avg_pls_3SD);
for i = 1:length(I_num_zeros_crista_terminalis_sorted)
    if I_num_zeros_crista_terminalis_dev_avg_pls_3SD(i) == min(I_num_zeros_crista_terminalis_dev_avg_pls_3SD)
        I_num_zeros_crista_terminalis_dev_avg_pls_3SD_index = i;
    end
end
I_num_zeros_crista_terminalis_tot_number_of_indices = length(I_num_zeros_crista_terminalis_sorted);
%This is the approximation of the tot "area under curve." 
I_num_zeros_crista_terminalis_quotient_within_1SD = (I_num_zeros_crista_terminalis_dev_avg_pls_1SD_index - I_num_zeros_crista_terminalisdev_avg_min_1SD_index)/I_num_zeros_crista_terminalis_tot_number_of_indices;
I_num_zeros_crista_terminalis_quotient_within_2SD = (I_num_zeros_crista_terminalis_dev_avg_pls_2SD_index - I_num_zeros_dev_avg_min_2SD_index)/I_num_zeros_crista_terminalis_tot_number_of_indices;
I_num_zeros_crista_terminalis_quotient_within_3SD = (I_num_zeros_crista_terminalis_dev_avg_pls_3SD_index - I_num_zeros_crista_terminalis_dev_avg_min_3SD_index)/I_num_zeros_crista_terminalis_tot_number_of_indices;
if (I_num_zeros_crista_terminalis_quotient_within_1SD > 0.66 && I_num_zeros_crista_terminalis_quotient_within_1SD < 0.70) && (I_num_zeros_crista_terminalis_quotient_within_2SD > 0.93 && I_num_zeros_crista_terminalis_quotient_within_2SD < 0.97) && (I_num_zeros_crista_terminalis_quotient_within_3SD > 0.98 && I_num_zeros_crista_terminalis_quotient_within_3SD < 1)
    isNormal(1,2) = 1;
end
if I_num_zeros_crista_terminalis_quotient_within_1SD == 0
    isSpike(1,2) = 1;
end
I_num_zeros_tricuspid_annulus_sorted = sort(I_num_zeros_tricuspid_annulus);
I_num_zeros_tricuspid_annulus_dev_avg_min_1SD = zeros(length(I_num_zeros_tricuspid_annulus_sorted),1);
for i = 1:length(I_num_zeros_tricuspid_annulus_sorted)
    I_num_zeros_tricuspid_annulus_dev_avg_min_1SD(i) = I_num_zeros_tricuspid_annulus_sorted(i) - (I_num_zeros_tricuspid_annulus_average-I_num_zeros_tricuspid_annulus_SD);
end
I_num_zeros_tricuspid_annulus_dev_avg_min_1SD = abs(I_num_zeros_tricuspid_annulus_dev_avg_min_1SD);
for i = 1:length(I_num_zeros_tricuspid_annulus_sorted)
    if I_num_zeros_tricuspid_annulus_dev_avg_min_1SD(i) == min(I_num_zeros_tricuspid_annulus_dev_avg_min_1SD)
        I_num_zeros_tricuspid_annulus_dev_avg_min_1SD_index = i;
    end
end
I_num_zeros_tricuspid_annulus_dev_avg_pls_1SD = zeros(length(I_num_zeros_tricuspid_annulus_sorted),1);
for i = 1:length(I_num_zeros_tricuspid_annulus_sorted)
    I_num_zeros_tricuspid_annulus_dev_avg_pls_1SD(i) = I_num_zeros_tricuspid_annulus_sorted(i) - (I_num_zeros_tricuspid_annulus_average+I_num_zeros_tricuspid_annulus_SD);
end
I_num_zeros_tricuspid_annulus_dev_avg_pls_1SD = abs(I_num_zeros_tricuspid_annulus_dev_avg_pls_1SD);
for i = 1:length(I_num_zeros_tricuspid_annulus_sorted)
    if I_num_zeros_tricuspid_annulus_dev_avg_pls_1SD(i) == min(I_num_zeros_tricuspid_annulus_dev_avg_pls_1SD)
        I_num_zeros_tricuspid_annulus_dev_avg_pls_1SD_index = i;
    end
end
I_num_zeros_tricuspid_annulus_dev_avg_min_2SD = zeros(length(I_num_zeros_tricuspid_annulus_sorted),1);
for i = 1:length(I_num_zeros_tricuspid_annulus_sorted)
    I_num_zeros_tricuspid_annulus_dev_avg_min_2SD(i) = I_num_zeros_tricuspid_annulus_sorted(i) - (I_num_zeros_tricuspid_annulus_average-2*I_num_zeros_tricuspid_annulus_SD);
end
I_num_zeros_tricuspid_annulus_dev_avg_min_2SD = abs(I_num_zeros_tricuspid_annulus_dev_avg_min_2SD);
for i = 1:length(I_num_zeros_tricuspid_annulus_sorted)
    if I_num_zeros_tricuspid_annulus_dev_avg_min_2SD(i) == min(I_num_zeros_tricuspid_annulus_dev_avg_min_2SD)
        I_num_zeros_dev_avg_min_2SD_index = i;
    end
end
I_num_zeros_tricuspid_annulus_dev_avg_pls_2SD = zeros(length(I_num_zeros_tricuspid_annulus_sorted),1);
for i = 1:length(I_num_zeros_tricuspid_annulus_sorted)
    I_num_zeros_tricuspid_annulus_dev_avg_pls_2SD(i) = I_num_zeros_tricuspid_annulus_sorted(i) - (I_num_zeros_tricuspid_annulus_average+2*I_num_zeros_tricuspid_annulus_SD);
end
I_num_zeros_tricuspid_annulus_dev_avg_pls_2SD = abs(I_num_zeros_tricuspid_annulus_dev_avg_pls_2SD);
for i = 1:length(I_num_zeros_tricuspid_annulus_sorted)
    if I_num_zeros_tricuspid_annulus_dev_avg_pls_2SD(i) == min(I_num_zeros_tricuspid_annulus_dev_avg_pls_2SD)
        I_num_zeros_tricuspid_annulus_dev_avg_pls_2SD_index = i;
    end
end
I_num_zeros_tricuspid_annulus_dev_avg_min_3SD = zeros(length(I_num_zeros_tricuspid_annulus_sorted),1);
for i = 1:length(I_num_zeros_tricuspid_annulus_sorted)
    I_num_zeros_tricuspid_annulus_dev_avg_min_3SD(i) = I_num_zeros_tricuspid_annulus_sorted(i) - (I_num_zeros_tricuspid_annulus_average-3*I_num_zeros_tricuspid_annulus_SD);
end
I_num_zeros_tricuspid_annulus_dev_avg_min_3SD = abs(I_num_zeros_tricuspid_annulus_dev_avg_min_3SD);
for i = 1:length(I_num_zeros_tricuspid_annulus_sorted)
    if I_num_zeros_tricuspid_annulus_dev_avg_min_3SD(i) == min(I_num_zeros_tricuspid_annulus_dev_avg_min_3SD)
        I_num_zeros_tricuspid_annulus_dev_avg_min_3SD_index = i;
    end
end
I_num_zeros_tricuspid_annulus_dev_avg_pls_3SD = zeros(length(I_num_zeros_tricuspid_annulus_sorted),1);
for i = 1:length(I_num_zeros_tricuspid_annulus_sorted)
    I_num_zeros_tricuspid_annulus_dev_avg_pls_3SD(i) = I_num_zeros_tricuspid_annulus_sorted(i) - (I_num_zeros_tricuspid_annulus_average+3*I_num_zeros_tricuspid_annulus_SD);
end
I_num_zeros_tricuspid_annulus_dev_avg_pls_3SD = abs(I_num_zeros_tricuspid_annulus_dev_avg_pls_3SD);
for i = 1:length(I_num_zeros_tricuspid_annulus_sorted)
    if I_num_zeros_tricuspid_annulus_dev_avg_pls_3SD(i) == min(I_num_zeros_tricuspid_annulus_dev_avg_pls_3SD)
        I_num_zeros_tricuspid_annulus_dev_avg_pls_3SD_index = i;
    end
end
I_num_zeros_tricuspid_annulus_tot_number_of_indices = length(I_num_zeros_tricuspid_annulus_sorted);
%This is the approximation of the tot "area under curve." 
I_num_zeros_tricuspid_annulus_quotient_within_1SD = (I_num_zeros_tricuspid_annulus_dev_avg_pls_1SD_index - I_num_zeros_tricuspid_annulusdev_avg_min_1SD_index)/I_num_zeros_tricuspid_annulus_tot_number_of_indices;
I_num_zeros_tricuspid_annulus_quotient_within_2SD = (I_num_zeros_tricuspid_annulus_dev_avg_pls_2SD_index - I_num_zeros_dev_avg_min_2SD_index)/I_num_zeros_tricuspid_annulus_tot_number_of_indices;
I_num_zeros_tricuspid_annulus_quotient_within_3SD = (I_num_zeros_tricuspid_annulus_dev_avg_pls_3SD_index - I_num_zeros_tricuspid_annulus_dev_avg_min_3SD_index)/I_num_zeros_tricuspid_annulus_tot_number_of_indices;
if (I_num_zeros_tricuspid_annulus_quotient_within_1SD > 0.66 && I_num_zeros_tricuspid_annulus_quotient_within_1SD < 0.70) && (I_num_zeros_tricuspid_annulus_quotient_within_2SD > 0.93 && I_num_zeros_tricuspid_annulus_quotient_within_2SD < 0.97) && (I_num_zeros_tricuspid_annulus_quotient_within_3SD > 0.98 && I_num_zeros_tricuspid_annulus_quotient_within_3SD < 1)
    isNormal(1,3) = 1;
end
if I_num_zeros_tricuspid_annulus_quotient_within_1SD == 0
    isSpike(1,3) = 1;
end
I_num_zeros_coronary_sinus_sorted = sort(I_num_zeros_coronary_sinus);
I_num_zeros_coronary_sinus_dev_avg_min_1SD = zeros(length(I_num_zeros_coronary_sinus_sorted),1);
for i = 1:length(I_num_zeros_coronary_sinus_sorted)
    I_num_zeros_coronary_sinus_dev_avg_min_1SD(i) = I_num_zeros_coronary_sinus_sorted(i) - (I_num_zeros_coronary_sinus_average-I_num_zeros_coronary_sinus_SD);
end
I_num_zeros_coronary_sinus_dev_avg_min_1SD = abs(I_num_zeros_coronary_sinus_dev_avg_min_1SD);
for i = 1:length(I_num_zeros_coronary_sinus_sorted)
    if I_num_zeros_coronary_sinus_dev_avg_min_1SD(i) == min(I_num_zeros_coronary_sinus_dev_avg_min_1SD)
        I_num_zeros_coronary_sinus_dev_avg_min_1SD_index = i;
    end
end
I_num_zeros_coronary_sinus_dev_avg_pls_1SD = zeros(length(I_num_zeros_coronary_sinus_sorted),1);
for i = 1:length(I_num_zeros_coronary_sinus_sorted)
    I_num_zeros_coronary_sinus_dev_avg_pls_1SD(i) = I_num_zeros_coronary_sinus_sorted(i) - (I_num_zeros_coronary_sinus_average+I_num_zeros_coronary_sinus_SD);
end
I_num_zeros_coronary_sinus_dev_avg_pls_1SD = abs(I_num_zeros_coronary_sinus_dev_avg_pls_1SD);
for i = 1:length(I_num_zeros_coronary_sinus_sorted)
    if I_num_zeros_coronary_sinus_dev_avg_pls_1SD(i) == min(I_num_zeros_coronary_sinus_dev_avg_pls_1SD)
        I_num_zeros_coronary_sinus_dev_avg_pls_1SD_index = i;
    end
end
I_num_zeros_coronary_sinus_dev_avg_min_2SD = zeros(length(I_num_zeros_coronary_sinus_sorted),1);
for i = 1:length(I_num_zeros_coronary_sinus_sorted)
    I_num_zeros_coronary_sinus_dev_avg_min_2SD(i) = I_num_zeros_coronary_sinus_sorted(i) - (I_num_zeros_coronary_sinus_average-2*I_num_zeros_coronary_sinus_SD);
end
I_num_zeros_coronary_sinus_dev_avg_min_2SD = abs(I_num_zeros_coronary_sinus_dev_avg_min_2SD);
for i = 1:length(I_num_zeros_coronary_sinus_sorted)
    if I_num_zeros_coronary_sinus_dev_avg_min_2SD(i) == min(I_num_zeros_coronary_sinus_dev_avg_min_2SD)
        I_num_zeros_dev_avg_min_2SD_index = i;
    end
end
I_num_zeros_coronary_sinus_dev_avg_pls_2SD = zeros(length(I_num_zeros_coronary_sinus_sorted),1);
for i = 1:length(I_num_zeros_coronary_sinus_sorted)
    I_num_zeros_coronary_sinus_dev_avg_pls_2SD(i) = I_num_zeros_coronary_sinus_sorted(i) - (I_num_zeros_coronary_sinus_average+2*I_num_zeros_coronary_sinus_SD);
end
I_num_zeros_coronary_sinus_dev_avg_pls_2SD = abs(I_num_zeros_coronary_sinus_dev_avg_pls_2SD);
for i = 1:length(I_num_zeros_coronary_sinus_sorted)
    if I_num_zeros_coronary_sinus_dev_avg_pls_2SD(i) == min(I_num_zeros_coronary_sinus_dev_avg_pls_2SD)
        I_num_zeros_coronary_sinus_dev_avg_pls_2SD_index = i;
    end
end
I_num_zeros_coronary_sinus_dev_avg_min_3SD = zeros(length(I_num_zeros_coronary_sinus_sorted),1);
for i = 1:length(I_num_zeros_coronary_sinus_sorted)
    I_num_zeros_coronary_sinus_dev_avg_min_3SD(i) = I_num_zeros_coronary_sinus_sorted(i) - (I_num_zeros_coronary_sinus_average-3*I_num_zeros_coronary_sinus_SD);
end
I_num_zeros_coronary_sinus_dev_avg_min_3SD = abs(I_num_zeros_coronary_sinus_dev_avg_min_3SD);
for i = 1:length(I_num_zeros_coronary_sinus_sorted)
    if I_num_zeros_coronary_sinus_dev_avg_min_3SD(i) == min(I_num_zeros_coronary_sinus_dev_avg_min_3SD)
        I_num_zeros_coronary_sinus_dev_avg_min_3SD_index = i;
    end
end
I_num_zeros_coronary_sinus_dev_avg_pls_3SD = zeros(length(I_num_zeros_coronary_sinus_sorted),1);
for i = 1:length(I_num_zeros_coronary_sinus_sorted)
    I_num_zeros_coronary_sinus_dev_avg_pls_3SD(i) = I_num_zeros_coronary_sinus_sorted(i) - (I_num_zeros_coronary_sinus_average+3*I_num_zeros_coronary_sinus_SD);
end
I_num_zeros_coronary_sinus_dev_avg_pls_3SD = abs(I_num_zeros_coronary_sinus_dev_avg_pls_3SD);
for i = 1:length(I_num_zeros_coronary_sinus_sorted)
    if I_num_zeros_coronary_sinus_dev_avg_pls_3SD(i) == min(I_num_zeros_coronary_sinus_dev_avg_pls_3SD)
        I_num_zeros_coronary_sinus_dev_avg_pls_3SD_index = i;
    end
end
I_num_zeros_coronary_sinus_tot_number_of_indices = length(I_num_zeros_coronary_sinus_sorted);
%This is the approximation of the tot "area under curve." 
I_num_zeros_coronary_sinus_quotient_within_1SD = (I_num_zeros_coronary_sinus_dev_avg_pls_1SD_index - I_num_zeros_coronary_sinusdev_avg_min_1SD_index)/I_num_zeros_coronary_sinus_tot_number_of_indices;
I_num_zeros_coronary_sinus_quotient_within_2SD = (I_num_zeros_coronary_sinus_dev_avg_pls_2SD_index - I_num_zeros_dev_avg_min_2SD_index)/I_num_zeros_coronary_sinus_tot_number_of_indices;
I_num_zeros_coronary_sinus_quotient_within_3SD = (I_num_zeros_coronary_sinus_dev_avg_pls_3SD_index - I_num_zeros_coronary_sinus_dev_avg_min_3SD_index)/I_num_zeros_coronary_sinus_tot_number_of_indices;
if (I_num_zeros_coronary_sinus_quotient_within_1SD > 0.66 && I_num_zeros_coronary_sinus_quotient_within_1SD < 0.70) && (I_num_zeros_coronary_sinus_quotient_within_2SD > 0.93 && I_num_zeros_coronary_sinus_quotient_within_2SD < 0.97) && (I_num_zeros_coronary_sinus_quotient_within_3SD > 0.98 && I_num_zeros_coronary_sinus_quotient_within_3SD < 1)
    isNormal(1,4) = 1;
end
if I_num_zeros_coronary_sinus_quotient_within_1SD == 0
    isSpike(1,4) = 1;
end
I_num_zeros_ostium_sorted = sort(I_num_zeros_ostium);
I_num_zeros_ostium_dev_avg_min_1SD = zeros(length(I_num_zeros_ostium_sorted),1);
for i = 1:length(I_num_zeros_ostium_sorted)
    I_num_zeros_ostium_dev_avg_min_1SD(i) = I_num_zeros_ostium_sorted(i) - (I_num_zeros_ostium_average-I_num_zeros_ostium_SD);
end
I_num_zeros_ostium_dev_avg_min_1SD = abs(I_num_zeros_ostium_dev_avg_min_1SD);
for i = 1:length(I_num_zeros_ostium_sorted)
    if I_num_zeros_ostium_dev_avg_min_1SD(i) == min(I_num_zeros_ostium_dev_avg_min_1SD)
        I_num_zeros_ostium_dev_avg_min_1SD_index = i;
    end
end
I_num_zeros_ostium_dev_avg_pls_1SD = zeros(length(I_num_zeros_ostium_sorted),1);
for i = 1:length(I_num_zeros_ostium_sorted)
    I_num_zeros_ostium_dev_avg_pls_1SD(i) = I_num_zeros_ostium_sorted(i) - (I_num_zeros_ostium_average+I_num_zeros_ostium_SD);
end
I_num_zeros_ostium_dev_avg_pls_1SD = abs(I_num_zeros_ostium_dev_avg_pls_1SD);
for i = 1:length(I_num_zeros_ostium_sorted)
    if I_num_zeros_ostium_dev_avg_pls_1SD(i) == min(I_num_zeros_ostium_dev_avg_pls_1SD)
        I_num_zeros_ostium_dev_avg_pls_1SD_index = i;
    end
end
I_num_zeros_ostium_dev_avg_min_2SD = zeros(length(I_num_zeros_ostium_sorted),1);
for i = 1:length(I_num_zeros_ostium_sorted)
    I_num_zeros_ostium_dev_avg_min_2SD(i) = I_num_zeros_ostium_sorted(i) - (I_num_zeros_ostium_average-2*I_num_zeros_ostium_SD);
end
I_num_zeros_ostium_dev_avg_min_2SD = abs(I_num_zeros_ostium_dev_avg_min_2SD);
for i = 1:length(I_num_zeros_ostium_sorted)
    if I_num_zeros_ostium_dev_avg_min_2SD(i) == min(I_num_zeros_ostium_dev_avg_min_2SD)
        I_num_zeros_dev_avg_min_2SD_index = i;
    end
end
I_num_zeros_ostium_dev_avg_pls_2SD = zeros(length(I_num_zeros_ostium_sorted),1);
for i = 1:length(I_num_zeros_ostium_sorted)
    I_num_zeros_ostium_dev_avg_pls_2SD(i) = I_num_zeros_ostium_sorted(i) - (I_num_zeros_ostium_average+2*I_num_zeros_ostium_SD);
end
I_num_zeros_ostium_dev_avg_pls_2SD = abs(I_num_zeros_ostium_dev_avg_pls_2SD);
for i = 1:length(I_num_zeros_ostium_sorted)
    if I_num_zeros_ostium_dev_avg_pls_2SD(i) == min(I_num_zeros_ostium_dev_avg_pls_2SD)
        I_num_zeros_ostium_dev_avg_pls_2SD_index = i;
    end
end
I_num_zeros_ostium_dev_avg_min_3SD = zeros(length(I_num_zeros_ostium_sorted),1);
for i = 1:length(I_num_zeros_ostium_sorted)
    I_num_zeros_ostium_dev_avg_min_3SD(i) = I_num_zeros_ostium_sorted(i) - (I_num_zeros_ostium_average-3*I_num_zeros_ostium_SD);
end
I_num_zeros_ostium_dev_avg_min_3SD = abs(I_num_zeros_ostium_dev_avg_min_3SD);
for i = 1:length(I_num_zeros_ostium_sorted)
    if I_num_zeros_ostium_dev_avg_min_3SD(i) == min(I_num_zeros_ostium_dev_avg_min_3SD)
        I_num_zeros_ostium_dev_avg_min_3SD_index = i;
    end
end
I_num_zeros_ostium_dev_avg_pls_3SD = zeros(length(I_num_zeros_ostium_sorted),1);
for i = 1:length(I_num_zeros_ostium_sorted)
    I_num_zeros_ostium_dev_avg_pls_3SD(i) = I_num_zeros_ostium_sorted(i) - (I_num_zeros_ostium_average+3*I_num_zeros_ostium_SD);
end
I_num_zeros_ostium_dev_avg_pls_3SD = abs(I_num_zeros_ostium_dev_avg_pls_3SD);
for i = 1:length(I_num_zeros_ostium_sorted)
    if I_num_zeros_ostium_dev_avg_pls_3SD(i) == min(I_num_zeros_ostium_dev_avg_pls_3SD)
        I_num_zeros_ostium_dev_avg_pls_3SD_index = i;
    end
end
I_num_zeros_ostium_tot_number_of_indices = length(I_num_zeros_ostium_sorted);
%This is the approximation of the tot "area under curve." 
I_num_zeros_ostium_quotient_within_1SD = (I_num_zeros_ostium_dev_avg_pls_1SD_index - I_num_zeros_ostiumdev_avg_min_1SD_index)/I_num_zeros_ostium_tot_number_of_indices;
I_num_zeros_ostium_quotient_within_2SD = (I_num_zeros_ostium_dev_avg_pls_2SD_index - I_num_zeros_dev_avg_min_2SD_index)/I_num_zeros_ostium_tot_number_of_indices;
I_num_zeros_ostium_quotient_within_3SD = (I_num_zeros_ostium_dev_avg_pls_3SD_index - I_num_zeros_ostium_dev_avg_min_3SD_index)/I_num_zeros_ostium_tot_number_of_indices;
if (I_num_zeros_ostium_quotient_within_1SD > 0.66 && I_num_zeros_ostium_quotient_within_1SD < 0.70) && (I_num_zeros_ostium_quotient_within_2SD > 0.93 && I_num_zeros_ostium_quotient_within_2SD < 0.97) && (I_num_zeros_ostium_quotient_within_3SD > 0.98 && I_num_zeros_ostium_quotient_within_3SD < 1)
    isNormal(1,5) = 1;
end
if I_num_zeros_ostium_quotient_within_1SD == 0
    isSpike(1,5) = 1;
end
I_num_zeros_perinodal_sorted = sort(I_num_zeros_perinodal);
I_num_zeros_perinodal_dev_avg_min_1SD = zeros(length(I_num_zeros_perinodal_sorted),1);
for i = 1:length(I_num_zeros_perinodal_sorted)
    I_num_zeros_perinodal_dev_avg_min_1SD(i) = I_num_zeros_perinodal_sorted(i) - (I_num_zeros_perinodal_average-I_num_zeros_perinodal_SD);
end
I_num_zeros_perinodal_dev_avg_min_1SD = abs(I_num_zeros_perinodal_dev_avg_min_1SD);
for i = 1:length(I_num_zeros_perinodal_sorted)
    if I_num_zeros_perinodal_dev_avg_min_1SD(i) == min(I_num_zeros_perinodal_dev_avg_min_1SD)
        I_num_zeros_perinodal_dev_avg_min_1SD_index = i;
    end
end
I_num_zeros_perinodal_dev_avg_pls_1SD = zeros(length(I_num_zeros_perinodal_sorted),1);
for i = 1:length(I_num_zeros_perinodal_sorted)
    I_num_zeros_perinodal_dev_avg_pls_1SD(i) = I_num_zeros_perinodal_sorted(i) - (I_num_zeros_perinodal_average+I_num_zeros_perinodal_SD);
end
I_num_zeros_perinodal_dev_avg_pls_1SD = abs(I_num_zeros_perinodal_dev_avg_pls_1SD);
for i = 1:length(I_num_zeros_perinodal_sorted)
    if I_num_zeros_perinodal_dev_avg_pls_1SD(i) == min(I_num_zeros_perinodal_dev_avg_pls_1SD)
        I_num_zeros_perinodal_dev_avg_pls_1SD_index = i;
    end
end
I_num_zeros_perinodal_dev_avg_min_2SD = zeros(length(I_num_zeros_perinodal_sorted),1);
for i = 1:length(I_num_zeros_perinodal_sorted)
    I_num_zeros_perinodal_dev_avg_min_2SD(i) = I_num_zeros_perinodal_sorted(i) - (I_num_zeros_perinodal_average-2*I_num_zeros_perinodal_SD);
end
I_num_zeros_perinodal_dev_avg_min_2SD = abs(I_num_zeros_perinodal_dev_avg_min_2SD);
for i = 1:length(I_num_zeros_perinodal_sorted)
    if I_num_zeros_perinodal_dev_avg_min_2SD(i) == min(I_num_zeros_perinodal_dev_avg_min_2SD)
        I_num_zeros_dev_avg_min_2SD_index = i;
    end
end
I_num_zeros_perinodal_dev_avg_pls_2SD = zeros(length(I_num_zeros_perinodal_sorted),1);
for i = 1:length(I_num_zeros_perinodal_sorted)
    I_num_zeros_perinodal_dev_avg_pls_2SD(i) = I_num_zeros_perinodal_sorted(i) - (I_num_zeros_perinodal_average+2*I_num_zeros_perinodal_SD);
end
I_num_zeros_perinodal_dev_avg_pls_2SD = abs(I_num_zeros_perinodal_dev_avg_pls_2SD);
for i = 1:length(I_num_zeros_perinodal_sorted)
    if I_num_zeros_perinodal_dev_avg_pls_2SD(i) == min(I_num_zeros_perinodal_dev_avg_pls_2SD)
        I_num_zeros_perinodal_dev_avg_pls_2SD_index = i;
    end
end
I_num_zeros_perinodal_dev_avg_min_3SD = zeros(length(I_num_zeros_perinodal_sorted),1);
for i = 1:length(I_num_zeros_perinodal_sorted)
    I_num_zeros_perinodal_dev_avg_min_3SD(i) = I_num_zeros_perinodal_sorted(i) - (I_num_zeros_perinodal_average-3*I_num_zeros_perinodal_SD);
end
I_num_zeros_perinodal_dev_avg_min_3SD = abs(I_num_zeros_perinodal_dev_avg_min_3SD);
for i = 1:length(I_num_zeros_perinodal_sorted)
    if I_num_zeros_perinodal_dev_avg_min_3SD(i) == min(I_num_zeros_perinodal_dev_avg_min_3SD)
        I_num_zeros_perinodal_dev_avg_min_3SD_index = i;
    end
end
I_num_zeros_perinodal_dev_avg_pls_3SD = zeros(length(I_num_zeros_perinodal_sorted),1);
for i = 1:length(I_num_zeros_perinodal_sorted)
    I_num_zeros_perinodal_dev_avg_pls_3SD(i) = I_num_zeros_perinodal_sorted(i) - (I_num_zeros_perinodal_average+3*I_num_zeros_perinodal_SD);
end
I_num_zeros_perinodal_dev_avg_pls_3SD = abs(I_num_zeros_perinodal_dev_avg_pls_3SD);
for i = 1:length(I_num_zeros_perinodal_sorted)
    if I_num_zeros_perinodal_dev_avg_pls_3SD(i) == min(I_num_zeros_perinodal_dev_avg_pls_3SD)
        I_num_zeros_perinodal_dev_avg_pls_3SD_index = i;
    end
end
I_num_zeros_perinodal_tot_number_of_indices = length(I_num_zeros_perinodal_sorted);
%This is the approximation of the tot "area under curve." 
I_num_zeros_perinodal_quotient_within_1SD = (I_num_zeros_perinodal_dev_avg_pls_1SD_index - I_num_zeros_perinodaldev_avg_min_1SD_index)/I_num_zeros_perinodal_tot_number_of_indices;
I_num_zeros_perinodal_quotient_within_2SD = (I_num_zeros_perinodal_dev_avg_pls_2SD_index - I_num_zeros_dev_avg_min_2SD_index)/I_num_zeros_perinodal_tot_number_of_indices;
I_num_zeros_perinodal_quotient_within_3SD = (I_num_zeros_perinodal_dev_avg_pls_3SD_index - I_num_zeros_perinodal_dev_avg_min_3SD_index)/I_num_zeros_perinodal_tot_number_of_indices;
if (I_num_zeros_perinodal_quotient_within_1SD > 0.66 && I_num_zeros_perinodal_quotient_within_1SD < 0.70) && (I_num_zeros_perinodal_quotient_within_2SD > 0.93 && I_num_zeros_perinodal_quotient_within_2SD < 0.97) && (I_num_zeros_perinodal_quotient_within_3SD > 0.98 && I_num_zeros_perinodal_quotient_within_3SD < 1)
    isNormal(1,6) = 1;
end
if I_num_zeros_perinodal_quotient_within_1SD == 0
    isSpike(1,6) = 1;
end
I_num_zeros_right_septum_sorted = sort(I_num_zeros_right_septum);
I_num_zeros_right_septum_dev_avg_min_1SD = zeros(length(I_num_zeros_right_septum_sorted),1);
for i = 1:length(I_num_zeros_right_septum_sorted)
    I_num_zeros_right_septum_dev_avg_min_1SD(i) = I_num_zeros_right_septum_sorted(i) - (I_num_zeros_right_septum_average-I_num_zeros_right_septum_SD);
end
I_num_zeros_right_septum_dev_avg_min_1SD = abs(I_num_zeros_right_septum_dev_avg_min_1SD);
for i = 1:length(I_num_zeros_right_septum_sorted)
    if I_num_zeros_right_septum_dev_avg_min_1SD(i) == min(I_num_zeros_right_septum_dev_avg_min_1SD)
        I_num_zeros_right_septum_dev_avg_min_1SD_index = i;
    end
end
I_num_zeros_right_septum_dev_avg_pls_1SD = zeros(length(I_num_zeros_right_septum_sorted),1);
for i = 1:length(I_num_zeros_right_septum_sorted)
    I_num_zeros_right_septum_dev_avg_pls_1SD(i) = I_num_zeros_right_septum_sorted(i) - (I_num_zeros_right_septum_average+I_num_zeros_right_septum_SD);
end
I_num_zeros_right_septum_dev_avg_pls_1SD = abs(I_num_zeros_right_septum_dev_avg_pls_1SD);
for i = 1:length(I_num_zeros_right_septum_sorted)
    if I_num_zeros_right_septum_dev_avg_pls_1SD(i) == min(I_num_zeros_right_septum_dev_avg_pls_1SD)
        I_num_zeros_right_septum_dev_avg_pls_1SD_index = i;
    end
end
I_num_zeros_right_septum_dev_avg_min_2SD = zeros(length(I_num_zeros_right_septum_sorted),1);
for i = 1:length(I_num_zeros_right_septum_sorted)
    I_num_zeros_right_septum_dev_avg_min_2SD(i) = I_num_zeros_right_septum_sorted(i) - (I_num_zeros_right_septum_average-2*I_num_zeros_right_septum_SD);
end
I_num_zeros_right_septum_dev_avg_min_2SD = abs(I_num_zeros_right_septum_dev_avg_min_2SD);
for i = 1:length(I_num_zeros_right_septum_sorted)
    if I_num_zeros_right_septum_dev_avg_min_2SD(i) == min(I_num_zeros_right_septum_dev_avg_min_2SD)
        I_num_zeros_dev_avg_min_2SD_index = i;
    end
end
I_num_zeros_right_septum_dev_avg_pls_2SD = zeros(length(I_num_zeros_right_septum_sorted),1);
for i = 1:length(I_num_zeros_right_septum_sorted)
    I_num_zeros_right_septum_dev_avg_pls_2SD(i) = I_num_zeros_right_septum_sorted(i) - (I_num_zeros_right_septum_average+2*I_num_zeros_right_septum_SD);
end
I_num_zeros_right_septum_dev_avg_pls_2SD = abs(I_num_zeros_right_septum_dev_avg_pls_2SD);
for i = 1:length(I_num_zeros_right_septum_sorted)
    if I_num_zeros_right_septum_dev_avg_pls_2SD(i) == min(I_num_zeros_right_septum_dev_avg_pls_2SD)
        I_num_zeros_right_septum_dev_avg_pls_2SD_index = i;
    end
end
I_num_zeros_right_septum_dev_avg_min_3SD = zeros(length(I_num_zeros_right_septum_sorted),1);
for i = 1:length(I_num_zeros_right_septum_sorted)
    I_num_zeros_right_septum_dev_avg_min_3SD(i) = I_num_zeros_right_septum_sorted(i) - (I_num_zeros_right_septum_average-3*I_num_zeros_right_septum_SD);
end
I_num_zeros_right_septum_dev_avg_min_3SD = abs(I_num_zeros_right_septum_dev_avg_min_3SD);
for i = 1:length(I_num_zeros_right_septum_sorted)
    if I_num_zeros_right_septum_dev_avg_min_3SD(i) == min(I_num_zeros_right_septum_dev_avg_min_3SD)
        I_num_zeros_right_septum_dev_avg_min_3SD_index = i;
    end
end
I_num_zeros_right_septum_dev_avg_pls_3SD = zeros(length(I_num_zeros_right_septum_sorted),1);
for i = 1:length(I_num_zeros_right_septum_sorted)
    I_num_zeros_right_septum_dev_avg_pls_3SD(i) = I_num_zeros_right_septum_sorted(i) - (I_num_zeros_right_septum_average+3*I_num_zeros_right_septum_SD);
end
I_num_zeros_right_septum_dev_avg_pls_3SD = abs(I_num_zeros_right_septum_dev_avg_pls_3SD);
for i = 1:length(I_num_zeros_right_septum_sorted)
    if I_num_zeros_right_septum_dev_avg_pls_3SD(i) == min(I_num_zeros_right_septum_dev_avg_pls_3SD)
        I_num_zeros_right_septum_dev_avg_pls_3SD_index = i;
    end
end
I_num_zeros_right_septum_tot_number_of_indices = length(I_num_zeros_right_septum_sorted);
%This is the approximation of the tot "area under curve." 
I_num_zeros_right_septum_quotient_within_1SD = (I_num_zeros_right_septum_dev_avg_pls_1SD_index - I_num_zeros_right_septumdev_avg_min_1SD_index)/I_num_zeros_right_septum_tot_number_of_indices;
I_num_zeros_right_septum_quotient_within_2SD = (I_num_zeros_right_septum_dev_avg_pls_2SD_index - I_num_zeros_dev_avg_min_2SD_index)/I_num_zeros_right_septum_tot_number_of_indices;
I_num_zeros_right_septum_quotient_within_3SD = (I_num_zeros_right_septum_dev_avg_pls_3SD_index - I_num_zeros_right_septum_dev_avg_min_3SD_index)/I_num_zeros_right_septum_tot_number_of_indices;
if (I_num_zeros_right_septum_quotient_within_1SD > 0.66 && I_num_zeros_right_septum_quotient_within_1SD < 0.70) && (I_num_zeros_right_septum_quotient_within_2SD > 0.93 && I_num_zeros_right_septum_quotient_within_2SD < 0.97) && (I_num_zeros_right_septum_quotient_within_3SD > 0.98 && I_num_zeros_right_septum_quotient_within_3SD < 1)
    isNormal(1,7) = 1;
end
if I_num_zeros_right_septum_quotient_within_1SD == 0
    isSpike(1,7) = 1;
end
I_num_zeros_right_atrial_appendage_sorted = sort(I_num_zeros_right_atrial_appendage);
I_num_zeros_right_atrial_appendage_dev_avg_min_1SD = zeros(length(I_num_zeros_right_atrial_appendage_sorted),1);
for i = 1:length(I_num_zeros_right_atrial_appendage_sorted)
    I_num_zeros_right_atrial_appendage_dev_avg_min_1SD(i) = I_num_zeros_right_atrial_appendage_sorted(i) - (I_num_zeros_right_atrial_appendage_average-I_num_zeros_right_atrial_appendage_SD);
end
I_num_zeros_right_atrial_appendage_dev_avg_min_1SD = abs(I_num_zeros_right_atrial_appendage_dev_avg_min_1SD);
for i = 1:length(I_num_zeros_right_atrial_appendage_sorted)
    if I_num_zeros_right_atrial_appendage_dev_avg_min_1SD(i) == min(I_num_zeros_right_atrial_appendage_dev_avg_min_1SD)
        I_num_zeros_right_atrial_appendage_dev_avg_min_1SD_index = i;
    end
end
I_num_zeros_right_atrial_appendage_dev_avg_pls_1SD = zeros(length(I_num_zeros_right_atrial_appendage_sorted),1);
for i = 1:length(I_num_zeros_right_atrial_appendage_sorted)
    I_num_zeros_right_atrial_appendage_dev_avg_pls_1SD(i) = I_num_zeros_right_atrial_appendage_sorted(i) - (I_num_zeros_right_atrial_appendage_average+I_num_zeros_right_atrial_appendage_SD);
end
I_num_zeros_right_atrial_appendage_dev_avg_pls_1SD = abs(I_num_zeros_right_atrial_appendage_dev_avg_pls_1SD);
for i = 1:length(I_num_zeros_right_atrial_appendage_sorted)
    if I_num_zeros_right_atrial_appendage_dev_avg_pls_1SD(i) == min(I_num_zeros_right_atrial_appendage_dev_avg_pls_1SD)
        I_num_zeros_right_atrial_appendage_dev_avg_pls_1SD_index = i;
    end
end
I_num_zeros_right_atrial_appendage_dev_avg_min_2SD = zeros(length(I_num_zeros_right_atrial_appendage_sorted),1);
for i = 1:length(I_num_zeros_right_atrial_appendage_sorted)
    I_num_zeros_right_atrial_appendage_dev_avg_min_2SD(i) = I_num_zeros_right_atrial_appendage_sorted(i) - (I_num_zeros_right_atrial_appendage_average-2*I_num_zeros_right_atrial_appendage_SD);
end
I_num_zeros_right_atrial_appendage_dev_avg_min_2SD = abs(I_num_zeros_right_atrial_appendage_dev_avg_min_2SD);
for i = 1:length(I_num_zeros_right_atrial_appendage_sorted)
    if I_num_zeros_right_atrial_appendage_dev_avg_min_2SD(i) == min(I_num_zeros_right_atrial_appendage_dev_avg_min_2SD)
        I_num_zeros_dev_avg_min_2SD_index = i;
    end
end
I_num_zeros_right_atrial_appendage_dev_avg_pls_2SD = zeros(length(I_num_zeros_right_atrial_appendage_sorted),1);
for i = 1:length(I_num_zeros_right_atrial_appendage_sorted)
    I_num_zeros_right_atrial_appendage_dev_avg_pls_2SD(i) = I_num_zeros_right_atrial_appendage_sorted(i) - (I_num_zeros_right_atrial_appendage_average+2*I_num_zeros_right_atrial_appendage_SD);
end
I_num_zeros_right_atrial_appendage_dev_avg_pls_2SD = abs(I_num_zeros_right_atrial_appendage_dev_avg_pls_2SD);
for i = 1:length(I_num_zeros_right_atrial_appendage_sorted)
    if I_num_zeros_right_atrial_appendage_dev_avg_pls_2SD(i) == min(I_num_zeros_right_atrial_appendage_dev_avg_pls_2SD)
        I_num_zeros_right_atrial_appendage_dev_avg_pls_2SD_index = i;
    end
end
I_num_zeros_right_atrial_appendage_dev_avg_min_3SD = zeros(length(I_num_zeros_right_atrial_appendage_sorted),1);
for i = 1:length(I_num_zeros_right_atrial_appendage_sorted)
    I_num_zeros_right_atrial_appendage_dev_avg_min_3SD(i) = I_num_zeros_right_atrial_appendage_sorted(i) - (I_num_zeros_right_atrial_appendage_average-3*I_num_zeros_right_atrial_appendage_SD);
end
I_num_zeros_right_atrial_appendage_dev_avg_min_3SD = abs(I_num_zeros_right_atrial_appendage_dev_avg_min_3SD);
for i = 1:length(I_num_zeros_right_atrial_appendage_sorted)
    if I_num_zeros_right_atrial_appendage_dev_avg_min_3SD(i) == min(I_num_zeros_right_atrial_appendage_dev_avg_min_3SD)
        I_num_zeros_right_atrial_appendage_dev_avg_min_3SD_index = i;
    end
end
I_num_zeros_right_atrial_appendage_dev_avg_pls_3SD = zeros(length(I_num_zeros_right_atrial_appendage_sorted),1);
for i = 1:length(I_num_zeros_right_atrial_appendage_sorted)
    I_num_zeros_right_atrial_appendage_dev_avg_pls_3SD(i) = I_num_zeros_right_atrial_appendage_sorted(i) - (I_num_zeros_right_atrial_appendage_average+3*I_num_zeros_right_atrial_appendage_SD);
end
I_num_zeros_right_atrial_appendage_dev_avg_pls_3SD = abs(I_num_zeros_right_atrial_appendage_dev_avg_pls_3SD);
for i = 1:length(I_num_zeros_right_atrial_appendage_sorted)
    if I_num_zeros_right_atrial_appendage_dev_avg_pls_3SD(i) == min(I_num_zeros_right_atrial_appendage_dev_avg_pls_3SD)
        I_num_zeros_right_atrial_appendage_dev_avg_pls_3SD_index = i;
    end
end
I_num_zeros_right_atrial_appendage_tot_number_of_indices = length(I_num_zeros_right_atrial_appendage_sorted);
%This is the approximation of the tot "area under curve." 
I_num_zeros_right_atrial_appendage_quotient_within_1SD = (I_num_zeros_right_atrial_appendage_dev_avg_pls_1SD_index - I_num_zeros_right_atrial_appendagedev_avg_min_1SD_index)/I_num_zeros_right_atrial_appendage_tot_number_of_indices;
I_num_zeros_right_atrial_appendage_quotient_within_2SD = (I_num_zeros_right_atrial_appendage_dev_avg_pls_2SD_index - I_num_zeros_dev_avg_min_2SD_index)/I_num_zeros_right_atrial_appendage_tot_number_of_indices;
I_num_zeros_right_atrial_appendage_quotient_within_3SD = (I_num_zeros_right_atrial_appendage_dev_avg_pls_3SD_index - I_num_zeros_right_atrial_appendage_dev_avg_min_3SD_index)/I_num_zeros_right_atrial_appendage_tot_number_of_indices;
if (I_num_zeros_right_atrial_appendage_quotient_within_1SD > 0.66 && I_num_zeros_right_atrial_appendage_quotient_within_1SD < 0.70) && (I_num_zeros_right_atrial_appendage_quotient_within_2SD > 0.93 && I_num_zeros_right_atrial_appendage_quotient_within_2SD < 0.97) && (I_num_zeros_right_atrial_appendage_quotient_within_3SD > 0.98 && I_num_zeros_right_atrial_appendage_quotient_within_3SD < 1)
    isNormal(1,8) = 1;
end
if I_num_zeros_right_atrial_appendage_quotient_within_1SD == 0
    isSpike(1,8) = 1;
end
I_num_zeros_pulmonary_veins_sorted = sort(I_num_zeros_pulmonary_veins);
I_num_zeros_pulmonary_veins_dev_avg_min_1SD = zeros(length(I_num_zeros_pulmonary_veins_sorted),1);
for i = 1:length(I_num_zeros_pulmonary_veins_sorted)
    I_num_zeros_pulmonary_veins_dev_avg_min_1SD(i) = I_num_zeros_pulmonary_veins_sorted(i) - (I_num_zeros_pulmonary_veins_average-I_num_zeros_pulmonary_veins_SD);
end
I_num_zeros_pulmonary_veins_dev_avg_min_1SD = abs(I_num_zeros_pulmonary_veins_dev_avg_min_1SD);
for i = 1:length(I_num_zeros_pulmonary_veins_sorted)
    if I_num_zeros_pulmonary_veins_dev_avg_min_1SD(i) == min(I_num_zeros_pulmonary_veins_dev_avg_min_1SD)
        I_num_zeros_pulmonary_veins_dev_avg_min_1SD_index = i;
    end
end
I_num_zeros_pulmonary_veins_dev_avg_pls_1SD = zeros(length(I_num_zeros_pulmonary_veins_sorted),1);
for i = 1:length(I_num_zeros_pulmonary_veins_sorted)
    I_num_zeros_pulmonary_veins_dev_avg_pls_1SD(i) = I_num_zeros_pulmonary_veins_sorted(i) - (I_num_zeros_pulmonary_veins_average+I_num_zeros_pulmonary_veins_SD);
end
I_num_zeros_pulmonary_veins_dev_avg_pls_1SD = abs(I_num_zeros_pulmonary_veins_dev_avg_pls_1SD);
for i = 1:length(I_num_zeros_pulmonary_veins_sorted)
    if I_num_zeros_pulmonary_veins_dev_avg_pls_1SD(i) == min(I_num_zeros_pulmonary_veins_dev_avg_pls_1SD)
        I_num_zeros_pulmonary_veins_dev_avg_pls_1SD_index = i;
    end
end
I_num_zeros_pulmonary_veins_dev_avg_min_2SD = zeros(length(I_num_zeros_pulmonary_veins_sorted),1);
for i = 1:length(I_num_zeros_pulmonary_veins_sorted)
    I_num_zeros_pulmonary_veins_dev_avg_min_2SD(i) = I_num_zeros_pulmonary_veins_sorted(i) - (I_num_zeros_pulmonary_veins_average-2*I_num_zeros_pulmonary_veins_SD);
end
I_num_zeros_pulmonary_veins_dev_avg_min_2SD = abs(I_num_zeros_pulmonary_veins_dev_avg_min_2SD);
for i = 1:length(I_num_zeros_pulmonary_veins_sorted)
    if I_num_zeros_pulmonary_veins_dev_avg_min_2SD(i) == min(I_num_zeros_pulmonary_veins_dev_avg_min_2SD)
        I_num_zeros_dev_avg_min_2SD_index = i;
    end
end
I_num_zeros_pulmonary_veins_dev_avg_pls_2SD = zeros(length(I_num_zeros_pulmonary_veins_sorted),1);
for i = 1:length(I_num_zeros_pulmonary_veins_sorted)
    I_num_zeros_pulmonary_veins_dev_avg_pls_2SD(i) = I_num_zeros_pulmonary_veins_sorted(i) - (I_num_zeros_pulmonary_veins_average+2*I_num_zeros_pulmonary_veins_SD);
end
I_num_zeros_pulmonary_veins_dev_avg_pls_2SD = abs(I_num_zeros_pulmonary_veins_dev_avg_pls_2SD);
for i = 1:length(I_num_zeros_pulmonary_veins_sorted)
    if I_num_zeros_pulmonary_veins_dev_avg_pls_2SD(i) == min(I_num_zeros_pulmonary_veins_dev_avg_pls_2SD)
        I_num_zeros_pulmonary_veins_dev_avg_pls_2SD_index = i;
    end
end
I_num_zeros_pulmonary_veins_dev_avg_min_3SD = zeros(length(I_num_zeros_pulmonary_veins_sorted),1);
for i = 1:length(I_num_zeros_pulmonary_veins_sorted)
    I_num_zeros_pulmonary_veins_dev_avg_min_3SD(i) = I_num_zeros_pulmonary_veins_sorted(i) - (I_num_zeros_pulmonary_veins_average-3*I_num_zeros_pulmonary_veins_SD);
end
I_num_zeros_pulmonary_veins_dev_avg_min_3SD = abs(I_num_zeros_pulmonary_veins_dev_avg_min_3SD);
for i = 1:length(I_num_zeros_pulmonary_veins_sorted)
    if I_num_zeros_pulmonary_veins_dev_avg_min_3SD(i) == min(I_num_zeros_pulmonary_veins_dev_avg_min_3SD)
        I_num_zeros_pulmonary_veins_dev_avg_min_3SD_index = i;
    end
end
I_num_zeros_pulmonary_veins_dev_avg_pls_3SD = zeros(length(I_num_zeros_pulmonary_veins_sorted),1);
for i = 1:length(I_num_zeros_pulmonary_veins_sorted)
    I_num_zeros_pulmonary_veins_dev_avg_pls_3SD(i) = I_num_zeros_pulmonary_veins_sorted(i) - (I_num_zeros_pulmonary_veins_average+3*I_num_zeros_pulmonary_veins_SD);
end
I_num_zeros_pulmonary_veins_dev_avg_pls_3SD = abs(I_num_zeros_pulmonary_veins_dev_avg_pls_3SD);
for i = 1:length(I_num_zeros_pulmonary_veins_sorted)
    if I_num_zeros_pulmonary_veins_dev_avg_pls_3SD(i) == min(I_num_zeros_pulmonary_veins_dev_avg_pls_3SD)
        I_num_zeros_pulmonary_veins_dev_avg_pls_3SD_index = i;
    end
end
I_num_zeros_pulmonary_veins_tot_number_of_indices = length(I_num_zeros_pulmonary_veins_sorted);
%This is the approximation of the tot "area under curve." 
I_num_zeros_pulmonary_veins_quotient_within_1SD = (I_num_zeros_pulmonary_veins_dev_avg_pls_1SD_index - I_num_zeros_pulmonary_veinsdev_avg_min_1SD_index)/I_num_zeros_pulmonary_veins_tot_number_of_indices;
I_num_zeros_pulmonary_veins_quotient_within_2SD = (I_num_zeros_pulmonary_veins_dev_avg_pls_2SD_index - I_num_zeros_dev_avg_min_2SD_index)/I_num_zeros_pulmonary_veins_tot_number_of_indices;
I_num_zeros_pulmonary_veins_quotient_within_3SD = (I_num_zeros_pulmonary_veins_dev_avg_pls_3SD_index - I_num_zeros_pulmonary_veins_dev_avg_min_3SD_index)/I_num_zeros_pulmonary_veins_tot_number_of_indices;
if (I_num_zeros_pulmonary_veins_quotient_within_1SD > 0.66 && I_num_zeros_pulmonary_veins_quotient_within_1SD < 0.70) && (I_num_zeros_pulmonary_veins_quotient_within_2SD > 0.93 && I_num_zeros_pulmonary_veins_quotient_within_2SD < 0.97) && (I_num_zeros_pulmonary_veins_quotient_within_3SD > 0.98 && I_num_zeros_pulmonary_veins_quotient_within_3SD < 1)
    isNormal(1,9) = 1;
end
if I_num_zeros_pulmonary_veins_quotient_within_1SD == 0
    isSpike(1,9) = 1;
end
I_num_zeros_mitral_annulus_sorted = sort(I_num_zeros_mitral_annulus);
I_num_zeros_mitral_annulus_dev_avg_min_1SD = zeros(length(I_num_zeros_mitral_annulus_sorted),1);
for i = 1:length(I_num_zeros_mitral_annulus_sorted)
    I_num_zeros_mitral_annulus_dev_avg_min_1SD(i) = I_num_zeros_mitral_annulus_sorted(i) - (I_num_zeros_mitral_annulus_average-I_num_zeros_mitral_annulus_SD);
end
I_num_zeros_mitral_annulus_dev_avg_min_1SD = abs(I_num_zeros_mitral_annulus_dev_avg_min_1SD);
for i = 1:length(I_num_zeros_mitral_annulus_sorted)
    if I_num_zeros_mitral_annulus_dev_avg_min_1SD(i) == min(I_num_zeros_mitral_annulus_dev_avg_min_1SD)
        I_num_zeros_mitral_annulus_dev_avg_min_1SD_index = i;
    end
end
I_num_zeros_mitral_annulus_dev_avg_pls_1SD = zeros(length(I_num_zeros_mitral_annulus_sorted),1);
for i = 1:length(I_num_zeros_mitral_annulus_sorted)
    I_num_zeros_mitral_annulus_dev_avg_pls_1SD(i) = I_num_zeros_mitral_annulus_sorted(i) - (I_num_zeros_mitral_annulus_average+I_num_zeros_mitral_annulus_SD);
end
I_num_zeros_mitral_annulus_dev_avg_pls_1SD = abs(I_num_zeros_mitral_annulus_dev_avg_pls_1SD);
for i = 1:length(I_num_zeros_mitral_annulus_sorted)
    if I_num_zeros_mitral_annulus_dev_avg_pls_1SD(i) == min(I_num_zeros_mitral_annulus_dev_avg_pls_1SD)
        I_num_zeros_mitral_annulus_dev_avg_pls_1SD_index = i;
    end
end
I_num_zeros_mitral_annulus_dev_avg_min_2SD = zeros(length(I_num_zeros_mitral_annulus_sorted),1);
for i = 1:length(I_num_zeros_mitral_annulus_sorted)
    I_num_zeros_mitral_annulus_dev_avg_min_2SD(i) = I_num_zeros_mitral_annulus_sorted(i) - (I_num_zeros_mitral_annulus_average-2*I_num_zeros_mitral_annulus_SD);
end
I_num_zeros_mitral_annulus_dev_avg_min_2SD = abs(I_num_zeros_mitral_annulus_dev_avg_min_2SD);
for i = 1:length(I_num_zeros_mitral_annulus_sorted)
    if I_num_zeros_mitral_annulus_dev_avg_min_2SD(i) == min(I_num_zeros_mitral_annulus_dev_avg_min_2SD)
        I_num_zeros_dev_avg_min_2SD_index = i;
    end
end
I_num_zeros_mitral_annulus_dev_avg_pls_2SD = zeros(length(I_num_zeros_mitral_annulus_sorted),1);
for i = 1:length(I_num_zeros_mitral_annulus_sorted)
    I_num_zeros_mitral_annulus_dev_avg_pls_2SD(i) = I_num_zeros_mitral_annulus_sorted(i) - (I_num_zeros_mitral_annulus_average+2*I_num_zeros_mitral_annulus_SD);
end
I_num_zeros_mitral_annulus_dev_avg_pls_2SD = abs(I_num_zeros_mitral_annulus_dev_avg_pls_2SD);
for i = 1:length(I_num_zeros_mitral_annulus_sorted)
    if I_num_zeros_mitral_annulus_dev_avg_pls_2SD(i) == min(I_num_zeros_mitral_annulus_dev_avg_pls_2SD)
        I_num_zeros_mitral_annulus_dev_avg_pls_2SD_index = i;
    end
end
I_num_zeros_mitral_annulus_dev_avg_min_3SD = zeros(length(I_num_zeros_mitral_annulus_sorted),1);
for i = 1:length(I_num_zeros_mitral_annulus_sorted)
    I_num_zeros_mitral_annulus_dev_avg_min_3SD(i) = I_num_zeros_mitral_annulus_sorted(i) - (I_num_zeros_mitral_annulus_average-3*I_num_zeros_mitral_annulus_SD);
end
I_num_zeros_mitral_annulus_dev_avg_min_3SD = abs(I_num_zeros_mitral_annulus_dev_avg_min_3SD);
for i = 1:length(I_num_zeros_mitral_annulus_sorted)
    if I_num_zeros_mitral_annulus_dev_avg_min_3SD(i) == min(I_num_zeros_mitral_annulus_dev_avg_min_3SD)
        I_num_zeros_mitral_annulus_dev_avg_min_3SD_index = i;
    end
end
I_num_zeros_mitral_annulus_dev_avg_pls_3SD = zeros(length(I_num_zeros_mitral_annulus_sorted),1);
for i = 1:length(I_num_zeros_mitral_annulus_sorted)
    I_num_zeros_mitral_annulus_dev_avg_pls_3SD(i) = I_num_zeros_mitral_annulus_sorted(i) - (I_num_zeros_mitral_annulus_average+3*I_num_zeros_mitral_annulus_SD);
end
I_num_zeros_mitral_annulus_dev_avg_pls_3SD = abs(I_num_zeros_mitral_annulus_dev_avg_pls_3SD);
for i = 1:length(I_num_zeros_mitral_annulus_sorted)
    if I_num_zeros_mitral_annulus_dev_avg_pls_3SD(i) == min(I_num_zeros_mitral_annulus_dev_avg_pls_3SD)
        I_num_zeros_mitral_annulus_dev_avg_pls_3SD_index = i;
    end
end
I_num_zeros_mitral_annulus_tot_number_of_indices = length(I_num_zeros_mitral_annulus_sorted);
%This is the approximation of the tot "area under curve." 
I_num_zeros_mitral_annulus_quotient_within_1SD = (I_num_zeros_mitral_annulus_dev_avg_pls_1SD_index - I_num_zeros_mitral_annulusdev_avg_min_1SD_index)/I_num_zeros_mitral_annulus_tot_number_of_indices;
I_num_zeros_mitral_annulus_quotient_within_2SD = (I_num_zeros_mitral_annulus_dev_avg_pls_2SD_index - I_num_zeros_dev_avg_min_2SD_index)/I_num_zeros_mitral_annulus_tot_number_of_indices;
I_num_zeros_mitral_annulus_quotient_within_3SD = (I_num_zeros_mitral_annulus_dev_avg_pls_3SD_index - I_num_zeros_mitral_annulus_dev_avg_min_3SD_index)/I_num_zeros_mitral_annulus_tot_number_of_indices;
if (I_num_zeros_mitral_annulus_quotient_within_1SD > 0.66 && I_num_zeros_mitral_annulus_quotient_within_1SD < 0.70) && (I_num_zeros_mitral_annulus_quotient_within_2SD > 0.93 && I_num_zeros_mitral_annulus_quotient_within_2SD < 0.97) && (I_num_zeros_mitral_annulus_quotient_within_3SD > 0.98 && I_num_zeros_mitral_annulus_quotient_within_3SD < 1)
    isNormal(1,10) = 1;
end
if I_num_zeros_mitral_annulus_quotient_within_1SD == 0
    isSpike(1,10) = 1;
end
I_num_zeros_CS_body_sorted = sort(I_num_zeros_CS_body);
I_num_zeros_CS_body_dev_avg_min_1SD = zeros(length(I_num_zeros_CS_body_sorted),1);
for i = 1:length(I_num_zeros_CS_body_sorted)
    I_num_zeros_CS_body_dev_avg_min_1SD(i) = I_num_zeros_CS_body_sorted(i) - (I_num_zeros_CS_body_average-I_num_zeros_CS_body_SD);
end
I_num_zeros_CS_body_dev_avg_min_1SD = abs(I_num_zeros_CS_body_dev_avg_min_1SD);
for i = 1:length(I_num_zeros_CS_body_sorted)
    if I_num_zeros_CS_body_dev_avg_min_1SD(i) == min(I_num_zeros_CS_body_dev_avg_min_1SD)
        I_num_zeros_CS_body_dev_avg_min_1SD_index = i;
    end
end
I_num_zeros_CS_body_dev_avg_pls_1SD = zeros(length(I_num_zeros_CS_body_sorted),1);
for i = 1:length(I_num_zeros_CS_body_sorted)
    I_num_zeros_CS_body_dev_avg_pls_1SD(i) = I_num_zeros_CS_body_sorted(i) - (I_num_zeros_CS_body_average+I_num_zeros_CS_body_SD);
end
I_num_zeros_CS_body_dev_avg_pls_1SD = abs(I_num_zeros_CS_body_dev_avg_pls_1SD);
for i = 1:length(I_num_zeros_CS_body_sorted)
    if I_num_zeros_CS_body_dev_avg_pls_1SD(i) == min(I_num_zeros_CS_body_dev_avg_pls_1SD)
        I_num_zeros_CS_body_dev_avg_pls_1SD_index = i;
    end
end
I_num_zeros_CS_body_dev_avg_min_2SD = zeros(length(I_num_zeros_CS_body_sorted),1);
for i = 1:length(I_num_zeros_CS_body_sorted)
    I_num_zeros_CS_body_dev_avg_min_2SD(i) = I_num_zeros_CS_body_sorted(i) - (I_num_zeros_CS_body_average-2*I_num_zeros_CS_body_SD);
end
I_num_zeros_CS_body_dev_avg_min_2SD = abs(I_num_zeros_CS_body_dev_avg_min_2SD);
for i = 1:length(I_num_zeros_CS_body_sorted)
    if I_num_zeros_CS_body_dev_avg_min_2SD(i) == min(I_num_zeros_CS_body_dev_avg_min_2SD)
        I_num_zeros_dev_avg_min_2SD_index = i;
    end
end
I_num_zeros_CS_body_dev_avg_pls_2SD = zeros(length(I_num_zeros_CS_body_sorted),1);
for i = 1:length(I_num_zeros_CS_body_sorted)
    I_num_zeros_CS_body_dev_avg_pls_2SD(i) = I_num_zeros_CS_body_sorted(i) - (I_num_zeros_CS_body_average+2*I_num_zeros_CS_body_SD);
end
I_num_zeros_CS_body_dev_avg_pls_2SD = abs(I_num_zeros_CS_body_dev_avg_pls_2SD);
for i = 1:length(I_num_zeros_CS_body_sorted)
    if I_num_zeros_CS_body_dev_avg_pls_2SD(i) == min(I_num_zeros_CS_body_dev_avg_pls_2SD)
        I_num_zeros_CS_body_dev_avg_pls_2SD_index = i;
    end
end
I_num_zeros_CS_body_dev_avg_min_3SD = zeros(length(I_num_zeros_CS_body_sorted),1);
for i = 1:length(I_num_zeros_CS_body_sorted)
    I_num_zeros_CS_body_dev_avg_min_3SD(i) = I_num_zeros_CS_body_sorted(i) - (I_num_zeros_CS_body_average-3*I_num_zeros_CS_body_SD);
end
I_num_zeros_CS_body_dev_avg_min_3SD = abs(I_num_zeros_CS_body_dev_avg_min_3SD);
for i = 1:length(I_num_zeros_CS_body_sorted)
    if I_num_zeros_CS_body_dev_avg_min_3SD(i) == min(I_num_zeros_CS_body_dev_avg_min_3SD)
        I_num_zeros_CS_body_dev_avg_min_3SD_index = i;
    end
end
I_num_zeros_CS_body_dev_avg_pls_3SD = zeros(length(I_num_zeros_CS_body_sorted),1);
for i = 1:length(I_num_zeros_CS_body_sorted)
    I_num_zeros_CS_body_dev_avg_pls_3SD(i) = I_num_zeros_CS_body_sorted(i) - (I_num_zeros_CS_body_average+3*I_num_zeros_CS_body_SD);
end
I_num_zeros_CS_body_dev_avg_pls_3SD = abs(I_num_zeros_CS_body_dev_avg_pls_3SD);
for i = 1:length(I_num_zeros_CS_body_sorted)
    if I_num_zeros_CS_body_dev_avg_pls_3SD(i) == min(I_num_zeros_CS_body_dev_avg_pls_3SD)
        I_num_zeros_CS_body_dev_avg_pls_3SD_index = i;
    end
end
I_num_zeros_CS_body_tot_number_of_indices = length(I_num_zeros_CS_body_sorted);
%This is the approximation of the tot "area under curve." 
I_num_zeros_CS_body_quotient_within_1SD = (I_num_zeros_CS_body_dev_avg_pls_1SD_index - I_num_zeros_CS_bodydev_avg_min_1SD_index)/I_num_zeros_CS_body_tot_number_of_indices;
I_num_zeros_CS_body_quotient_within_2SD = (I_num_zeros_CS_body_dev_avg_pls_2SD_index - I_num_zeros_dev_avg_min_2SD_index)/I_num_zeros_CS_body_tot_number_of_indices;
I_num_zeros_CS_body_quotient_within_3SD = (I_num_zeros_CS_body_dev_avg_pls_3SD_index - I_num_zeros_CS_body_dev_avg_min_3SD_index)/I_num_zeros_CS_body_tot_number_of_indices;
if (I_num_zeros_CS_body_quotient_within_1SD > 0.66 && I_num_zeros_CS_body_quotient_within_1SD < 0.70) && (I_num_zeros_CS_body_quotient_within_2SD > 0.93 && I_num_zeros_CS_body_quotient_within_2SD < 0.97) && (I_num_zeros_CS_body_quotient_within_3SD > 0.98 && I_num_zeros_CS_body_quotient_within_3SD < 1)
    isNormal(1,11) = 1;
end
if I_num_zeros_CS_body_quotient_within_1SD == 0
    isSpike(1,11) = 1;
end
I_num_zeros_left_septum_sorted = sort(I_num_zeros_left_septum);
I_num_zeros_left_septum_dev_avg_min_1SD = zeros(length(I_num_zeros_left_septum_sorted),1);
for i = 1:length(I_num_zeros_left_septum_sorted)
    I_num_zeros_left_septum_dev_avg_min_1SD(i) = I_num_zeros_left_septum_sorted(i) - (I_num_zeros_left_septum_average-I_num_zeros_left_septum_SD);
end
I_num_zeros_left_septum_dev_avg_min_1SD = abs(I_num_zeros_left_septum_dev_avg_min_1SD);
for i = 1:length(I_num_zeros_left_septum_sorted)
    if I_num_zeros_left_septum_dev_avg_min_1SD(i) == min(I_num_zeros_left_septum_dev_avg_min_1SD)
        I_num_zeros_left_septum_dev_avg_min_1SD_index = i;
    end
end
I_num_zeros_left_septum_dev_avg_pls_1SD = zeros(length(I_num_zeros_left_septum_sorted),1);
for i = 1:length(I_num_zeros_left_septum_sorted)
    I_num_zeros_left_septum_dev_avg_pls_1SD(i) = I_num_zeros_left_septum_sorted(i) - (I_num_zeros_left_septum_average+I_num_zeros_left_septum_SD);
end
I_num_zeros_left_septum_dev_avg_pls_1SD = abs(I_num_zeros_left_septum_dev_avg_pls_1SD);
for i = 1:length(I_num_zeros_left_septum_sorted)
    if I_num_zeros_left_septum_dev_avg_pls_1SD(i) == min(I_num_zeros_left_septum_dev_avg_pls_1SD)
        I_num_zeros_left_septum_dev_avg_pls_1SD_index = i;
    end
end
I_num_zeros_left_septum_dev_avg_min_2SD = zeros(length(I_num_zeros_left_septum_sorted),1);
for i = 1:length(I_num_zeros_left_septum_sorted)
    I_num_zeros_left_septum_dev_avg_min_2SD(i) = I_num_zeros_left_septum_sorted(i) - (I_num_zeros_left_septum_average-2*I_num_zeros_left_septum_SD);
end
I_num_zeros_left_septum_dev_avg_min_2SD = abs(I_num_zeros_left_septum_dev_avg_min_2SD);
for i = 1:length(I_num_zeros_left_septum_sorted)
    if I_num_zeros_left_septum_dev_avg_min_2SD(i) == min(I_num_zeros_left_septum_dev_avg_min_2SD)
        I_num_zeros_dev_avg_min_2SD_index = i;
    end
end
I_num_zeros_left_septum_dev_avg_pls_2SD = zeros(length(I_num_zeros_left_septum_sorted),1);
for i = 1:length(I_num_zeros_left_septum_sorted)
    I_num_zeros_left_septum_dev_avg_pls_2SD(i) = I_num_zeros_left_septum_sorted(i) - (I_num_zeros_left_septum_average+2*I_num_zeros_left_septum_SD);
end
I_num_zeros_left_septum_dev_avg_pls_2SD = abs(I_num_zeros_left_septum_dev_avg_pls_2SD);
for i = 1:length(I_num_zeros_left_septum_sorted)
    if I_num_zeros_left_septum_dev_avg_pls_2SD(i) == min(I_num_zeros_left_septum_dev_avg_pls_2SD)
        I_num_zeros_left_septum_dev_avg_pls_2SD_index = i;
    end
end
I_num_zeros_left_septum_dev_avg_min_3SD = zeros(length(I_num_zeros_left_septum_sorted),1);
for i = 1:length(I_num_zeros_left_septum_sorted)
    I_num_zeros_left_septum_dev_avg_min_3SD(i) = I_num_zeros_left_septum_sorted(i) - (I_num_zeros_left_septum_average-3*I_num_zeros_left_septum_SD);
end
I_num_zeros_left_septum_dev_avg_min_3SD = abs(I_num_zeros_left_septum_dev_avg_min_3SD);
for i = 1:length(I_num_zeros_left_septum_sorted)
    if I_num_zeros_left_septum_dev_avg_min_3SD(i) == min(I_num_zeros_left_septum_dev_avg_min_3SD)
        I_num_zeros_left_septum_dev_avg_min_3SD_index = i;
    end
end
I_num_zeros_left_septum_dev_avg_pls_3SD = zeros(length(I_num_zeros_left_septum_sorted),1);
for i = 1:length(I_num_zeros_left_septum_sorted)
    I_num_zeros_left_septum_dev_avg_pls_3SD(i) = I_num_zeros_left_septum_sorted(i) - (I_num_zeros_left_septum_average+3*I_num_zeros_left_septum_SD);
end
I_num_zeros_left_septum_dev_avg_pls_3SD = abs(I_num_zeros_left_septum_dev_avg_pls_3SD);
for i = 1:length(I_num_zeros_left_septum_sorted)
    if I_num_zeros_left_septum_dev_avg_pls_3SD(i) == min(I_num_zeros_left_septum_dev_avg_pls_3SD)
        I_num_zeros_left_septum_dev_avg_pls_3SD_index = i;
    end
end
I_num_zeros_left_septum_tot_number_of_indices = length(I_num_zeros_left_septum_sorted);
%This is the approximation of the tot "area under curve." 
I_num_zeros_left_septum_quotient_within_1SD = (I_num_zeros_left_septum_dev_avg_pls_1SD_index - I_num_zeros_left_septumdev_avg_min_1SD_index)/I_num_zeros_left_septum_tot_number_of_indices;
I_num_zeros_left_septum_quotient_within_2SD = (I_num_zeros_left_septum_dev_avg_pls_2SD_index - I_num_zeros_dev_avg_min_2SD_index)/I_num_zeros_left_septum_tot_number_of_indices;
I_num_zeros_left_septum_quotient_within_3SD = (I_num_zeros_left_septum_dev_avg_pls_3SD_index - I_num_zeros_left_septum_dev_avg_min_3SD_index)/I_num_zeros_left_septum_tot_number_of_indices;
if (I_num_zeros_left_septum_quotient_within_1SD > 0.66 && I_num_zeros_left_septum_quotient_within_1SD < 0.70) && (I_num_zeros_left_septum_quotient_within_2SD > 0.93 && I_num_zeros_left_septum_quotient_within_2SD < 0.97) && (I_num_zeros_left_septum_quotient_within_3SD > 0.98 && I_num_zeros_left_septum_quotient_within_3SD < 1)
    isNormal(1,12) = 1;
end
if I_num_zeros_left_septum_quotient_within_1SD == 0
    isSpike(1,12) = 1;
end
I_num_zeros_left_atrial_appendage_sorted = sort(I_num_zeros_left_atrial_appendage);
I_num_zeros_left_atrial_appendage_dev_avg_min_1SD = zeros(length(I_num_zeros_left_atrial_appendage_sorted),1);
for i = 1:length(I_num_zeros_left_atrial_appendage_sorted)
    I_num_zeros_left_atrial_appendage_dev_avg_min_1SD(i) = I_num_zeros_left_atrial_appendage_sorted(i) - (I_num_zeros_left_atrial_appendage_average-I_num_zeros_left_atrial_appendage_SD);
end
I_num_zeros_left_atrial_appendage_dev_avg_min_1SD = abs(I_num_zeros_left_atrial_appendage_dev_avg_min_1SD);
for i = 1:length(I_num_zeros_left_atrial_appendage_sorted)
    if I_num_zeros_left_atrial_appendage_dev_avg_min_1SD(i) == min(I_num_zeros_left_atrial_appendage_dev_avg_min_1SD)
        I_num_zeros_left_atrial_appendage_dev_avg_min_1SD_index = i;
    end
end
I_num_zeros_left_atrial_appendage_dev_avg_pls_1SD = zeros(length(I_num_zeros_left_atrial_appendage_sorted),1);
for i = 1:length(I_num_zeros_left_atrial_appendage_sorted)
    I_num_zeros_left_atrial_appendage_dev_avg_pls_1SD(i) = I_num_zeros_left_atrial_appendage_sorted(i) - (I_num_zeros_left_atrial_appendage_average+I_num_zeros_left_atrial_appendage_SD);
end
I_num_zeros_left_atrial_appendage_dev_avg_pls_1SD = abs(I_num_zeros_left_atrial_appendage_dev_avg_pls_1SD);
for i = 1:length(I_num_zeros_left_atrial_appendage_sorted)
    if I_num_zeros_left_atrial_appendage_dev_avg_pls_1SD(i) == min(I_num_zeros_left_atrial_appendage_dev_avg_pls_1SD)
        I_num_zeros_left_atrial_appendage_dev_avg_pls_1SD_index = i;
    end
end
I_num_zeros_left_atrial_appendage_dev_avg_min_2SD = zeros(length(I_num_zeros_left_atrial_appendage_sorted),1);
for i = 1:length(I_num_zeros_left_atrial_appendage_sorted)
    I_num_zeros_left_atrial_appendage_dev_avg_min_2SD(i) = I_num_zeros_left_atrial_appendage_sorted(i) - (I_num_zeros_left_atrial_appendage_average-2*I_num_zeros_left_atrial_appendage_SD);
end
I_num_zeros_left_atrial_appendage_dev_avg_min_2SD = abs(I_num_zeros_left_atrial_appendage_dev_avg_min_2SD);
for i = 1:length(I_num_zeros_left_atrial_appendage_sorted)
    if I_num_zeros_left_atrial_appendage_dev_avg_min_2SD(i) == min(I_num_zeros_left_atrial_appendage_dev_avg_min_2SD)
        I_num_zeros_dev_avg_min_2SD_index = i;
    end
end
I_num_zeros_left_atrial_appendage_dev_avg_pls_2SD = zeros(length(I_num_zeros_left_atrial_appendage_sorted),1);
for i = 1:length(I_num_zeros_left_atrial_appendage_sorted)
    I_num_zeros_left_atrial_appendage_dev_avg_pls_2SD(i) = I_num_zeros_left_atrial_appendage_sorted(i) - (I_num_zeros_left_atrial_appendage_average+2*I_num_zeros_left_atrial_appendage_SD);
end
I_num_zeros_left_atrial_appendage_dev_avg_pls_2SD = abs(I_num_zeros_left_atrial_appendage_dev_avg_pls_2SD);
for i = 1:length(I_num_zeros_left_atrial_appendage_sorted)
    if I_num_zeros_left_atrial_appendage_dev_avg_pls_2SD(i) == min(I_num_zeros_left_atrial_appendage_dev_avg_pls_2SD)
        I_num_zeros_left_atrial_appendage_dev_avg_pls_2SD_index = i;
    end
end
I_num_zeros_left_atrial_appendage_dev_avg_min_3SD = zeros(length(I_num_zeros_left_atrial_appendage_sorted),1);
for i = 1:length(I_num_zeros_left_atrial_appendage_sorted)
    I_num_zeros_left_atrial_appendage_dev_avg_min_3SD(i) = I_num_zeros_left_atrial_appendage_sorted(i) - (I_num_zeros_left_atrial_appendage_average-3*I_num_zeros_left_atrial_appendage_SD);
end
I_num_zeros_left_atrial_appendage_dev_avg_min_3SD = abs(I_num_zeros_left_atrial_appendage_dev_avg_min_3SD);
for i = 1:length(I_num_zeros_left_atrial_appendage_sorted)
    if I_num_zeros_left_atrial_appendage_dev_avg_min_3SD(i) == min(I_num_zeros_left_atrial_appendage_dev_avg_min_3SD)
        I_num_zeros_left_atrial_appendage_dev_avg_min_3SD_index = i;
    end
end
I_num_zeros_left_atrial_appendage_dev_avg_pls_3SD = zeros(length(I_num_zeros_left_atrial_appendage_sorted),1);
for i = 1:length(I_num_zeros_left_atrial_appendage_sorted)
    I_num_zeros_left_atrial_appendage_dev_avg_pls_3SD(i) = I_num_zeros_left_atrial_appendage_sorted(i) - (I_num_zeros_left_atrial_appendage_average+3*I_num_zeros_left_atrial_appendage_SD);
end
I_num_zeros_left_atrial_appendage_dev_avg_pls_3SD = abs(I_num_zeros_left_atrial_appendage_dev_avg_pls_3SD);
for i = 1:length(I_num_zeros_left_atrial_appendage_sorted)
    if I_num_zeros_left_atrial_appendage_dev_avg_pls_3SD(i) == min(I_num_zeros_left_atrial_appendage_dev_avg_pls_3SD)
        I_num_zeros_left_atrial_appendage_dev_avg_pls_3SD_index = i;
    end
end
I_num_zeros_left_atrial_appendage_tot_number_of_indices = length(I_num_zeros_left_atrial_appendage_sorted);
%This is the approximation of the tot "area under curve." 
I_num_zeros_left_atrial_appendage_quotient_within_1SD = (I_num_zeros_left_atrial_appendage_dev_avg_pls_1SD_index - I_num_zeros_left_atrial_appendagedev_avg_min_1SD_index)/I_num_zeros_left_atrial_appendage_tot_number_of_indices;
I_num_zeros_left_atrial_appendage_quotient_within_2SD = (I_num_zeros_left_atrial_appendage_dev_avg_pls_2SD_index - I_num_zeros_dev_avg_min_2SD_index)/I_num_zeros_left_atrial_appendage_tot_number_of_indices;
I_num_zeros_left_atrial_appendage_quotient_within_3SD = (I_num_zeros_left_atrial_appendage_dev_avg_pls_3SD_index - I_num_zeros_left_atrial_appendage_dev_avg_min_3SD_index)/I_num_zeros_left_atrial_appendage_tot_number_of_indices;
if (I_num_zeros_left_atrial_appendage_quotient_within_1SD > 0.66 && I_num_zeros_left_atrial_appendage_quotient_within_1SD < 0.70) && (I_num_zeros_left_atrial_appendage_quotient_within_2SD > 0.93 && I_num_zeros_left_atrial_appendage_quotient_within_2SD < 0.97) && (I_num_zeros_left_atrial_appendage_quotient_within_3SD > 0.98 && I_num_zeros_left_atrial_appendage_quotient_within_3SD < 1)
    isNormal(1,13) = 1;
end
if I_num_zeros_left_atrial_appendage_quotient_within_1SD == 0
    isSpike(1,13) = 1;
end
I_maxima_SA_node_sorted = sort(I_maxima_SA_node);
I_maxima_SA_node_dev_avg_min_1SD = maxima(length(I_maxima_SA_node_sorted),1);
for i = 1:length(I_maxima_SA_node_sorted)
    I_maxima_SA_node_dev_avg_min_1SD(i) = I_maxima_SA_node_sorted(i) - (I_maxima_SA_node_average-I_maxima_SA_node_SD);
end
I_maxima_SA_node_dev_avg_min_1SD = abs(I_maxima_SA_node_dev_avg_min_1SD);
for i = 1:length(I_maxima_SA_node_sorted)
    if I_maxima_SA_node_dev_avg_min_1SD(i) == min(I_maxima_SA_node_dev_avg_min_1SD)
        I_maxima_SA_node_dev_avg_min_1SD_index = i;
    end
end
I_maxima_SA_node_dev_avg_pls_1SD = maxima(length(I_maxima_SA_node_sorted),1);
for i = 1:length(I_maxima_SA_node_sorted)
    I_maxima_SA_node_dev_avg_pls_1SD(i) = I_maxima_SA_node_sorted(i) - (I_maxima_SA_node_average+I_maxima_SA_node_SD);
end
I_maxima_SA_node_dev_avg_pls_1SD = abs(I_maxima_SA_node_dev_avg_pls_1SD);
for i = 1:length(I_maxima_SA_node_sorted)
    if I_maxima_SA_node_dev_avg_pls_1SD(i) == min(I_maxima_SA_node_dev_avg_pls_1SD)
        I_maxima_SA_node_dev_avg_pls_1SD_index = i;
    end
end
I_maxima_SA_node_dev_avg_min_2SD = maxima(length(I_maxima_SA_node_sorted),1);
for i = 1:length(I_maxima_SA_node_sorted)
    I_maxima_SA_node_dev_avg_min_2SD(i) = I_maxima_SA_node_sorted(i) - (I_maxima_SA_node_average-2*I_maxima_SA_node_SD);
end
I_maxima_SA_node_dev_avg_min_2SD = abs(I_maxima_SA_node_dev_avg_min_2SD);
for i = 1:length(I_maxima_SA_node_sorted)
    if I_maxima_SA_node_dev_avg_min_2SD(i) == min(I_maxima_SA_node_dev_avg_min_2SD)
        I_maxima_dev_avg_min_2SD_index = i;
    end
end
I_maxima_SA_node_dev_avg_pls_2SD = maxima(length(I_maxima_SA_node_sorted),1);
for i = 1:length(I_maxima_SA_node_sorted)
    I_maxima_SA_node_dev_avg_pls_2SD(i) = I_maxima_SA_node_sorted(i) - (I_maxima_SA_node_average+2*I_maxima_SA_node_SD);
end
I_maxima_SA_node_dev_avg_pls_2SD = abs(I_maxima_SA_node_dev_avg_pls_2SD);
for i = 1:length(I_maxima_SA_node_sorted)
    if I_maxima_SA_node_dev_avg_pls_2SD(i) == min(I_maxima_SA_node_dev_avg_pls_2SD)
        I_maxima_SA_node_dev_avg_pls_2SD_index = i;
    end
end
I_maxima_SA_node_dev_avg_min_3SD = maxima(length(I_maxima_SA_node_sorted),1);
for i = 1:length(I_maxima_SA_node_sorted)
    I_maxima_SA_node_dev_avg_min_3SD(i) = I_maxima_SA_node_sorted(i) - (I_maxima_SA_node_average-3*I_maxima_SA_node_SD);
end
I_maxima_SA_node_dev_avg_min_3SD = abs(I_maxima_SA_node_dev_avg_min_3SD);
for i = 1:length(I_maxima_SA_node_sorted)
    if I_maxima_SA_node_dev_avg_min_3SD(i) == min(I_maxima_SA_node_dev_avg_min_3SD)
        I_maxima_SA_node_dev_avg_min_3SD_index = i;
    end
end
I_maxima_SA_node_dev_avg_pls_3SD = maxima(length(I_maxima_SA_node_sorted),1);
for i = 1:length(I_maxima_SA_node_sorted)
    I_maxima_SA_node_dev_avg_pls_3SD(i) = I_maxima_SA_node_sorted(i) - (I_maxima_SA_node_average+3*I_maxima_SA_node_SD);
end
I_maxima_SA_node_dev_avg_pls_3SD = abs(I_maxima_SA_node_dev_avg_pls_3SD);
for i = 1:length(I_maxima_SA_node_sorted)
    if I_maxima_SA_node_dev_avg_pls_3SD(i) == min(I_maxima_SA_node_dev_avg_pls_3SD)
        I_maxima_SA_node_dev_avg_pls_3SD_index = i;
    end
end
I_maxima_SA_node_tot_number_of_indices = length(I_maxima_SA_node_sorted);
%This is the approximation of the tot "area under curve." 
I_maxima_SA_node_quotient_within_1SD = (I_maxima_SA_node_dev_avg_pls_1SD_index - I_maxima_SA_nodedev_avg_min_1SD_index)/I_maxima_SA_node_tot_number_of_indices;
I_maxima_SA_node_quotient_within_2SD = (I_maxima_SA_node_dev_avg_pls_2SD_index - I_maxima_dev_avg_min_2SD_index)/I_maxima_SA_node_tot_number_of_indices;
I_maxima_SA_node_quotient_within_3SD = (I_maxima_SA_node_dev_avg_pls_3SD_index - I_maxima_SA_node_dev_avg_min_3SD_index)/I_maxima_SA_node_tot_number_of_indices;
if ((I_maxima_SA_node_quotient_within_1SD > 0.66 && I_maxima_SA_node_quotient_within_1SD < 0.70) && (I_maxima_SA_node_quotient_within_2SD > 0.93 && I_maxima_SA_node_quotient_within_2SD < 0.97) && (I_maxima_SA_node_quotient_within_3SD > 0.98 && I_maxima_SA_node_quotient_within_3SD < 1)) 
    isNormal(2,1) = 1;
end
if I_maxima_SA_node_quotient_within_1SD == 0
    isSpike(2,1) = 1;
end
I_maxima_crista_terminalis_sorted = sort(I_maxima_crista_terminalis);
I_maxima_crista_terminalis_dev_avg_min_1SD = maxima(length(I_maxima_crista_terminalis_sorted),1);
for i = 1:length(I_maxima_crista_terminalis_sorted)
    I_maxima_crista_terminalis_dev_avg_min_1SD(i) = I_maxima_crista_terminalis_sorted(i) - (I_maxima_crista_terminalis_average-I_maxima_crista_terminalis_SD);
end
I_maxima_crista_terminalis_dev_avg_min_1SD = abs(I_maxima_crista_terminalis_dev_avg_min_1SD);
for i = 1:length(I_maxima_crista_terminalis_sorted)
    if I_maxima_crista_terminalis_dev_avg_min_1SD(i) == min(I_maxima_crista_terminalis_dev_avg_min_1SD)
        I_maxima_crista_terminalis_dev_avg_min_1SD_index = i;
    end
end
I_maxima_crista_terminalis_dev_avg_pls_1SD = maxima(length(I_maxima_crista_terminalis_sorted),1);
for i = 1:length(I_maxima_crista_terminalis_sorted)
    I_maxima_crista_terminalis_dev_avg_pls_1SD(i) = I_maxima_crista_terminalis_sorted(i) - (I_maxima_crista_terminalis_average+I_maxima_crista_terminalis_SD);
end
I_maxima_crista_terminalis_dev_avg_pls_1SD = abs(I_maxima_crista_terminalis_dev_avg_pls_1SD);
for i = 1:length(I_maxima_crista_terminalis_sorted)
    if I_maxima_crista_terminalis_dev_avg_pls_1SD(i) == min(I_maxima_crista_terminalis_dev_avg_pls_1SD)
        I_maxima_crista_terminalis_dev_avg_pls_1SD_index = i;
    end
end
I_maxima_crista_terminalis_dev_avg_min_2SD = maxima(length(I_maxima_crista_terminalis_sorted),1);
for i = 1:length(I_maxima_crista_terminalis_sorted)
    I_maxima_crista_terminalis_dev_avg_min_2SD(i) = I_maxima_crista_terminalis_sorted(i) - (I_maxima_crista_terminalis_average-2*I_maxima_crista_terminalis_SD);
end
I_maxima_crista_terminalis_dev_avg_min_2SD = abs(I_maxima_crista_terminalis_dev_avg_min_2SD);
for i = 1:length(I_maxima_crista_terminalis_sorted)
    if I_maxima_crista_terminalis_dev_avg_min_2SD(i) == min(I_maxima_crista_terminalis_dev_avg_min_2SD)
        I_maxima_dev_avg_min_2SD_index = i;
    end
end
I_maxima_crista_terminalis_dev_avg_pls_2SD = maxima(length(I_maxima_crista_terminalis_sorted),1);
for i = 1:length(I_maxima_crista_terminalis_sorted)
    I_maxima_crista_terminalis_dev_avg_pls_2SD(i) = I_maxima_crista_terminalis_sorted(i) - (I_maxima_crista_terminalis_average+2*I_maxima_crista_terminalis_SD);
end
I_maxima_crista_terminalis_dev_avg_pls_2SD = abs(I_maxima_crista_terminalis_dev_avg_pls_2SD);
for i = 1:length(I_maxima_crista_terminalis_sorted)
    if I_maxima_crista_terminalis_dev_avg_pls_2SD(i) == min(I_maxima_crista_terminalis_dev_avg_pls_2SD)
        I_maxima_crista_terminalis_dev_avg_pls_2SD_index = i;
    end
end
I_maxima_crista_terminalis_dev_avg_min_3SD = maxima(length(I_maxima_crista_terminalis_sorted),1);
for i = 1:length(I_maxima_crista_terminalis_sorted)
    I_maxima_crista_terminalis_dev_avg_min_3SD(i) = I_maxima_crista_terminalis_sorted(i) - (I_maxima_crista_terminalis_average-3*I_maxima_crista_terminalis_SD);
end
I_maxima_crista_terminalis_dev_avg_min_3SD = abs(I_maxima_crista_terminalis_dev_avg_min_3SD);
for i = 1:length(I_maxima_crista_terminalis_sorted)
    if I_maxima_crista_terminalis_dev_avg_min_3SD(i) == min(I_maxima_crista_terminalis_dev_avg_min_3SD)
        I_maxima_crista_terminalis_dev_avg_min_3SD_index = i;
    end
end
I_maxima_crista_terminalis_dev_avg_pls_3SD = maxima(length(I_maxima_crista_terminalis_sorted),1);
for i = 1:length(I_maxima_crista_terminalis_sorted)
    I_maxima_crista_terminalis_dev_avg_pls_3SD(i) = I_maxima_crista_terminalis_sorted(i) - (I_maxima_crista_terminalis_average+3*I_maxima_crista_terminalis_SD);
end
I_maxima_crista_terminalis_dev_avg_pls_3SD = abs(I_maxima_crista_terminalis_dev_avg_pls_3SD);
for i = 1:length(I_maxima_crista_terminalis_sorted)
    if I_maxima_crista_terminalis_dev_avg_pls_3SD(i) == min(I_maxima_crista_terminalis_dev_avg_pls_3SD)
        I_maxima_crista_terminalis_dev_avg_pls_3SD_index = i;
    end
end
I_maxima_crista_terminalis_tot_number_of_indices = length(I_maxima_crista_terminalis_sorted);
%This is the approximation of the tot "area under curve." 
I_maxima_crista_terminalis_quotient_within_1SD = (I_maxima_crista_terminalis_dev_avg_pls_1SD_index - I_maxima_crista_terminalisdev_avg_min_1SD_index)/I_maxima_crista_terminalis_tot_number_of_indices;
I_maxima_crista_terminalis_quotient_within_2SD = (I_maxima_crista_terminalis_dev_avg_pls_2SD_index - I_maxima_dev_avg_min_2SD_index)/I_maxima_crista_terminalis_tot_number_of_indices;
I_maxima_crista_terminalis_quotient_within_3SD = (I_maxima_crista_terminalis_dev_avg_pls_3SD_index - I_maxima_crista_terminalis_dev_avg_min_3SD_index)/I_maxima_crista_terminalis_tot_number_of_indices;
if (I_maxima_crista_terminalis_quotient_within_1SD > 0.66 && I_maxima_crista_terminalis_quotient_within_1SD < 0.70) && (I_maxima_crista_terminalis_quotient_within_2SD > 0.93 && I_maxima_crista_terminalis_quotient_within_2SD < 0.97) && (I_maxima_crista_terminalis_quotient_within_3SD > 0.98 && I_maxima_crista_terminalis_quotient_within_3SD < 1)
    isNormal(2,2) = 1;
end
if I_maxima_crista_terminalis_quotient_within_1SD == 0
    isSpike(2,2) = 1;
end
I_maxima_tricuspid_annulus_sorted = sort(I_maxima_tricuspid_annulus);
I_maxima_tricuspid_annulus_dev_avg_min_1SD = maxima(length(I_maxima_tricuspid_annulus_sorted),1);
for i = 1:length(I_maxima_tricuspid_annulus_sorted)
    I_maxima_tricuspid_annulus_dev_avg_min_1SD(i) = I_maxima_tricuspid_annulus_sorted(i) - (I_maxima_tricuspid_annulus_average-I_maxima_tricuspid_annulus_SD);
end
I_maxima_tricuspid_annulus_dev_avg_min_1SD = abs(I_maxima_tricuspid_annulus_dev_avg_min_1SD);
for i = 1:length(I_maxima_tricuspid_annulus_sorted)
    if I_maxima_tricuspid_annulus_dev_avg_min_1SD(i) == min(I_maxima_tricuspid_annulus_dev_avg_min_1SD)
        I_maxima_tricuspid_annulus_dev_avg_min_1SD_index = i;
    end
end
I_maxima_tricuspid_annulus_dev_avg_pls_1SD = maxima(length(I_maxima_tricuspid_annulus_sorted),1);
for i = 1:length(I_maxima_tricuspid_annulus_sorted)
    I_maxima_tricuspid_annulus_dev_avg_pls_1SD(i) = I_maxima_tricuspid_annulus_sorted(i) - (I_maxima_tricuspid_annulus_average+I_maxima_tricuspid_annulus_SD);
end
I_maxima_tricuspid_annulus_dev_avg_pls_1SD = abs(I_maxima_tricuspid_annulus_dev_avg_pls_1SD);
for i = 1:length(I_maxima_tricuspid_annulus_sorted)
    if I_maxima_tricuspid_annulus_dev_avg_pls_1SD(i) == min(I_maxima_tricuspid_annulus_dev_avg_pls_1SD)
        I_maxima_tricuspid_annulus_dev_avg_pls_1SD_index = i;
    end
end
I_maxima_tricuspid_annulus_dev_avg_min_2SD = maxima(length(I_maxima_tricuspid_annulus_sorted),1);
for i = 1:length(I_maxima_tricuspid_annulus_sorted)
    I_maxima_tricuspid_annulus_dev_avg_min_2SD(i) = I_maxima_tricuspid_annulus_sorted(i) - (I_maxima_tricuspid_annulus_average-2*I_maxima_tricuspid_annulus_SD);
end
I_maxima_tricuspid_annulus_dev_avg_min_2SD = abs(I_maxima_tricuspid_annulus_dev_avg_min_2SD);
for i = 1:length(I_maxima_tricuspid_annulus_sorted)
    if I_maxima_tricuspid_annulus_dev_avg_min_2SD(i) == min(I_maxima_tricuspid_annulus_dev_avg_min_2SD)
        I_maxima_dev_avg_min_2SD_index = i;
    end
end
I_maxima_tricuspid_annulus_dev_avg_pls_2SD = maxima(length(I_maxima_tricuspid_annulus_sorted),1);
for i = 1:length(I_maxima_tricuspid_annulus_sorted)
    I_maxima_tricuspid_annulus_dev_avg_pls_2SD(i) = I_maxima_tricuspid_annulus_sorted(i) - (I_maxima_tricuspid_annulus_average+2*I_maxima_tricuspid_annulus_SD);
end
I_maxima_tricuspid_annulus_dev_avg_pls_2SD = abs(I_maxima_tricuspid_annulus_dev_avg_pls_2SD);
for i = 1:length(I_maxima_tricuspid_annulus_sorted)
    if I_maxima_tricuspid_annulus_dev_avg_pls_2SD(i) == min(I_maxima_tricuspid_annulus_dev_avg_pls_2SD)
        I_maxima_tricuspid_annulus_dev_avg_pls_2SD_index = i;
    end
end
I_maxima_tricuspid_annulus_dev_avg_min_3SD = maxima(length(I_maxima_tricuspid_annulus_sorted),1);
for i = 1:length(I_maxima_tricuspid_annulus_sorted)
    I_maxima_tricuspid_annulus_dev_avg_min_3SD(i) = I_maxima_tricuspid_annulus_sorted(i) - (I_maxima_tricuspid_annulus_average-3*I_maxima_tricuspid_annulus_SD);
end
I_maxima_tricuspid_annulus_dev_avg_min_3SD = abs(I_maxima_tricuspid_annulus_dev_avg_min_3SD);
for i = 1:length(I_maxima_tricuspid_annulus_sorted)
    if I_maxima_tricuspid_annulus_dev_avg_min_3SD(i) == min(I_maxima_tricuspid_annulus_dev_avg_min_3SD)
        I_maxima_tricuspid_annulus_dev_avg_min_3SD_index = i;
    end
end
I_maxima_tricuspid_annulus_dev_avg_pls_3SD = maxima(length(I_maxima_tricuspid_annulus_sorted),1);
for i = 1:length(I_maxima_tricuspid_annulus_sorted)
    I_maxima_tricuspid_annulus_dev_avg_pls_3SD(i) = I_maxima_tricuspid_annulus_sorted(i) - (I_maxima_tricuspid_annulus_average+3*I_maxima_tricuspid_annulus_SD);
end
I_maxima_tricuspid_annulus_dev_avg_pls_3SD = abs(I_maxima_tricuspid_annulus_dev_avg_pls_3SD);
for i = 1:length(I_maxima_tricuspid_annulus_sorted)
    if I_maxima_tricuspid_annulus_dev_avg_pls_3SD(i) == min(I_maxima_tricuspid_annulus_dev_avg_pls_3SD)
        I_maxima_tricuspid_annulus_dev_avg_pls_3SD_index = i;
    end
end
I_maxima_tricuspid_annulus_tot_number_of_indices = length(I_maxima_tricuspid_annulus_sorted);
%This is the approximation of the tot "area under curve." 
I_maxima_tricuspid_annulus_quotient_within_1SD = (I_maxima_tricuspid_annulus_dev_avg_pls_1SD_index - I_maxima_tricuspid_annulusdev_avg_min_1SD_index)/I_maxima_tricuspid_annulus_tot_number_of_indices;
I_maxima_tricuspid_annulus_quotient_within_2SD = (I_maxima_tricuspid_annulus_dev_avg_pls_2SD_index - I_maxima_dev_avg_min_2SD_index)/I_maxima_tricuspid_annulus_tot_number_of_indices;
I_maxima_tricuspid_annulus_quotient_within_3SD = (I_maxima_tricuspid_annulus_dev_avg_pls_3SD_index - I_maxima_tricuspid_annulus_dev_avg_min_3SD_index)/I_maxima_tricuspid_annulus_tot_number_of_indices;
if (I_maxima_tricuspid_annulus_quotient_within_1SD > 0.66 && I_maxima_tricuspid_annulus_quotient_within_1SD < 0.70) && (I_maxima_tricuspid_annulus_quotient_within_2SD > 0.93 && I_maxima_tricuspid_annulus_quotient_within_2SD < 0.97) && (I_maxima_tricuspid_annulus_quotient_within_3SD > 0.98 && I_maxima_tricuspid_annulus_quotient_within_3SD < 1)
    isNormal(2,3) = 1;
end
if I_maxima_tricuspid_annulus_quotient_within_1SD == 0
    isSpike(2,3) = 1;
end
I_maxima_coronary_sinus_sorted = sort(I_maxima_coronary_sinus);
I_maxima_coronary_sinus_dev_avg_min_1SD = maxima(length(I_maxima_coronary_sinus_sorted),1);
for i = 1:length(I_maxima_coronary_sinus_sorted)
    I_maxima_coronary_sinus_dev_avg_min_1SD(i) = I_maxima_coronary_sinus_sorted(i) - (I_maxima_coronary_sinus_average-I_maxima_coronary_sinus_SD);
end
I_maxima_coronary_sinus_dev_avg_min_1SD = abs(I_maxima_coronary_sinus_dev_avg_min_1SD);
for i = 1:length(I_maxima_coronary_sinus_sorted)
    if I_maxima_coronary_sinus_dev_avg_min_1SD(i) == min(I_maxima_coronary_sinus_dev_avg_min_1SD)
        I_maxima_coronary_sinus_dev_avg_min_1SD_index = i;
    end
end
I_maxima_coronary_sinus_dev_avg_pls_1SD = maxima(length(I_maxima_coronary_sinus_sorted),1);
for i = 1:length(I_maxima_coronary_sinus_sorted)
    I_maxima_coronary_sinus_dev_avg_pls_1SD(i) = I_maxima_coronary_sinus_sorted(i) - (I_maxima_coronary_sinus_average+I_maxima_coronary_sinus_SD);
end
I_maxima_coronary_sinus_dev_avg_pls_1SD = abs(I_maxima_coronary_sinus_dev_avg_pls_1SD);
for i = 1:length(I_maxima_coronary_sinus_sorted)
    if I_maxima_coronary_sinus_dev_avg_pls_1SD(i) == min(I_maxima_coronary_sinus_dev_avg_pls_1SD)
        I_maxima_coronary_sinus_dev_avg_pls_1SD_index = i;
    end
end
I_maxima_coronary_sinus_dev_avg_min_2SD = maxima(length(I_maxima_coronary_sinus_sorted),1);
for i = 1:length(I_maxima_coronary_sinus_sorted)
    I_maxima_coronary_sinus_dev_avg_min_2SD(i) = I_maxima_coronary_sinus_sorted(i) - (I_maxima_coronary_sinus_average-2*I_maxima_coronary_sinus_SD);
end
I_maxima_coronary_sinus_dev_avg_min_2SD = abs(I_maxima_coronary_sinus_dev_avg_min_2SD);
for i = 1:length(I_maxima_coronary_sinus_sorted)
    if I_maxima_coronary_sinus_dev_avg_min_2SD(i) == min(I_maxima_coronary_sinus_dev_avg_min_2SD)
        I_maxima_dev_avg_min_2SD_index = i;
    end
end
I_maxima_coronary_sinus_dev_avg_pls_2SD = maxima(length(I_maxima_coronary_sinus_sorted),1);
for i = 1:length(I_maxima_coronary_sinus_sorted)
    I_maxima_coronary_sinus_dev_avg_pls_2SD(i) = I_maxima_coronary_sinus_sorted(i) - (I_maxima_coronary_sinus_average+2*I_maxima_coronary_sinus_SD);
end
I_maxima_coronary_sinus_dev_avg_pls_2SD = abs(I_maxima_coronary_sinus_dev_avg_pls_2SD);
for i = 1:length(I_maxima_coronary_sinus_sorted)
    if I_maxima_coronary_sinus_dev_avg_pls_2SD(i) == min(I_maxima_coronary_sinus_dev_avg_pls_2SD)
        I_maxima_coronary_sinus_dev_avg_pls_2SD_index = i;
    end
end
I_maxima_coronary_sinus_dev_avg_min_3SD = maxima(length(I_maxima_coronary_sinus_sorted),1);
for i = 1:length(I_maxima_coronary_sinus_sorted)
    I_maxima_coronary_sinus_dev_avg_min_3SD(i) = I_maxima_coronary_sinus_sorted(i) - (I_maxima_coronary_sinus_average-3*I_maxima_coronary_sinus_SD);
end
I_maxima_coronary_sinus_dev_avg_min_3SD = abs(I_maxima_coronary_sinus_dev_avg_min_3SD);
for i = 1:length(I_maxima_coronary_sinus_sorted)
    if I_maxima_coronary_sinus_dev_avg_min_3SD(i) == min(I_maxima_coronary_sinus_dev_avg_min_3SD)
        I_maxima_coronary_sinus_dev_avg_min_3SD_index = i;
    end
end
I_maxima_coronary_sinus_dev_avg_pls_3SD = maxima(length(I_maxima_coronary_sinus_sorted),1);
for i = 1:length(I_maxima_coronary_sinus_sorted)
    I_maxima_coronary_sinus_dev_avg_pls_3SD(i) = I_maxima_coronary_sinus_sorted(i) - (I_maxima_coronary_sinus_average+3*I_maxima_coronary_sinus_SD);
end
I_maxima_coronary_sinus_dev_avg_pls_3SD = abs(I_maxima_coronary_sinus_dev_avg_pls_3SD);
for i = 1:length(I_maxima_coronary_sinus_sorted)
    if I_maxima_coronary_sinus_dev_avg_pls_3SD(i) == min(I_maxima_coronary_sinus_dev_avg_pls_3SD)
        I_maxima_coronary_sinus_dev_avg_pls_3SD_index = i;
    end
end
I_maxima_coronary_sinus_tot_number_of_indices = length(I_maxima_coronary_sinus_sorted);
%This is the approximation of the tot "area under curve." 
I_maxima_coronary_sinus_quotient_within_1SD = (I_maxima_coronary_sinus_dev_avg_pls_1SD_index - I_maxima_coronary_sinusdev_avg_min_1SD_index)/I_maxima_coronary_sinus_tot_number_of_indices;
I_maxima_coronary_sinus_quotient_within_2SD = (I_maxima_coronary_sinus_dev_avg_pls_2SD_index - I_maxima_dev_avg_min_2SD_index)/I_maxima_coronary_sinus_tot_number_of_indices;
I_maxima_coronary_sinus_quotient_within_3SD = (I_maxima_coronary_sinus_dev_avg_pls_3SD_index - I_maxima_coronary_sinus_dev_avg_min_3SD_index)/I_maxima_coronary_sinus_tot_number_of_indices;
if (I_maxima_coronary_sinus_quotient_within_1SD > 0.66 && I_maxima_coronary_sinus_quotient_within_1SD < 0.70) && (I_maxima_coronary_sinus_quotient_within_2SD > 0.93 && I_maxima_coronary_sinus_quotient_within_2SD < 0.97) && (I_maxima_coronary_sinus_quotient_within_3SD > 0.98 && I_maxima_coronary_sinus_quotient_within_3SD < 1)
    isNormal(2,4) = 1;
end
if I_maxima_coronary_sinus_quotient_within_1SD == 0
    isSpike(2,4) = 1;
end
I_maxima_ostium_sorted = sort(I_maxima_ostium);
I_maxima_ostium_dev_avg_min_1SD = maxima(length(I_maxima_ostium_sorted),1);
for i = 1:length(I_maxima_ostium_sorted)
    I_maxima_ostium_dev_avg_min_1SD(i) = I_maxima_ostium_sorted(i) - (I_maxima_ostium_average-I_maxima_ostium_SD);
end
I_maxima_ostium_dev_avg_min_1SD = abs(I_maxima_ostium_dev_avg_min_1SD);
for i = 1:length(I_maxima_ostium_sorted)
    if I_maxima_ostium_dev_avg_min_1SD(i) == min(I_maxima_ostium_dev_avg_min_1SD)
        I_maxima_ostium_dev_avg_min_1SD_index = i;
    end
end
I_maxima_ostium_dev_avg_pls_1SD = maxima(length(I_maxima_ostium_sorted),1);
for i = 1:length(I_maxima_ostium_sorted)
    I_maxima_ostium_dev_avg_pls_1SD(i) = I_maxima_ostium_sorted(i) - (I_maxima_ostium_average+I_maxima_ostium_SD);
end
I_maxima_ostium_dev_avg_pls_1SD = abs(I_maxima_ostium_dev_avg_pls_1SD);
for i = 1:length(I_maxima_ostium_sorted)
    if I_maxima_ostium_dev_avg_pls_1SD(i) == min(I_maxima_ostium_dev_avg_pls_1SD)
        I_maxima_ostium_dev_avg_pls_1SD_index = i;
    end
end
I_maxima_ostium_dev_avg_min_2SD = maxima(length(I_maxima_ostium_sorted),1);
for i = 1:length(I_maxima_ostium_sorted)
    I_maxima_ostium_dev_avg_min_2SD(i) = I_maxima_ostium_sorted(i) - (I_maxima_ostium_average-2*I_maxima_ostium_SD);
end
I_maxima_ostium_dev_avg_min_2SD = abs(I_maxima_ostium_dev_avg_min_2SD);
for i = 1:length(I_maxima_ostium_sorted)
    if I_maxima_ostium_dev_avg_min_2SD(i) == min(I_maxima_ostium_dev_avg_min_2SD)
        I_maxima_dev_avg_min_2SD_index = i;
    end
end
I_maxima_ostium_dev_avg_pls_2SD = maxima(length(I_maxima_ostium_sorted),1);
for i = 1:length(I_maxima_ostium_sorted)
    I_maxima_ostium_dev_avg_pls_2SD(i) = I_maxima_ostium_sorted(i) - (I_maxima_ostium_average+2*I_maxima_ostium_SD);
end
I_maxima_ostium_dev_avg_pls_2SD = abs(I_maxima_ostium_dev_avg_pls_2SD);
for i = 1:length(I_maxima_ostium_sorted)
    if I_maxima_ostium_dev_avg_pls_2SD(i) == min(I_maxima_ostium_dev_avg_pls_2SD)
        I_maxima_ostium_dev_avg_pls_2SD_index = i;
    end
end
I_maxima_ostium_dev_avg_min_3SD = maxima(length(I_maxima_ostium_sorted),1);
for i = 1:length(I_maxima_ostium_sorted)
    I_maxima_ostium_dev_avg_min_3SD(i) = I_maxima_ostium_sorted(i) - (I_maxima_ostium_average-3*I_maxima_ostium_SD);
end
I_maxima_ostium_dev_avg_min_3SD = abs(I_maxima_ostium_dev_avg_min_3SD);
for i = 1:length(I_maxima_ostium_sorted)
    if I_maxima_ostium_dev_avg_min_3SD(i) == min(I_maxima_ostium_dev_avg_min_3SD)
        I_maxima_ostium_dev_avg_min_3SD_index = i;
    end
end
I_maxima_ostium_dev_avg_pls_3SD = maxima(length(I_maxima_ostium_sorted),1);
for i = 1:length(I_maxima_ostium_sorted)
    I_maxima_ostium_dev_avg_pls_3SD(i) = I_maxima_ostium_sorted(i) - (I_maxima_ostium_average+3*I_maxima_ostium_SD);
end
I_maxima_ostium_dev_avg_pls_3SD = abs(I_maxima_ostium_dev_avg_pls_3SD);
for i = 1:length(I_maxima_ostium_sorted)
    if I_maxima_ostium_dev_avg_pls_3SD(i) == min(I_maxima_ostium_dev_avg_pls_3SD)
        I_maxima_ostium_dev_avg_pls_3SD_index = i;
    end
end
I_maxima_ostium_tot_number_of_indices = length(I_maxima_ostium_sorted);
%This is the approximation of the tot "area under curve." 
I_maxima_ostium_quotient_within_1SD = (I_maxima_ostium_dev_avg_pls_1SD_index - I_maxima_ostiumdev_avg_min_1SD_index)/I_maxima_ostium_tot_number_of_indices;
I_maxima_ostium_quotient_within_2SD = (I_maxima_ostium_dev_avg_pls_2SD_index - I_maxima_dev_avg_min_2SD_index)/I_maxima_ostium_tot_number_of_indices;
I_maxima_ostium_quotient_within_3SD = (I_maxima_ostium_dev_avg_pls_3SD_index - I_maxima_ostium_dev_avg_min_3SD_index)/I_maxima_ostium_tot_number_of_indices;
if (I_maxima_ostium_quotient_within_1SD > 0.66 && I_maxima_ostium_quotient_within_1SD < 0.70) && (I_maxima_ostium_quotient_within_2SD > 0.93 && I_maxima_ostium_quotient_within_2SD < 0.97) && (I_maxima_ostium_quotient_within_3SD > 0.98 && I_maxima_ostium_quotient_within_3SD < 1)
    isNormal(2,5) = 1;
end
if I_maxima_ostium_quotient_within_1SD == 0
    isSpike(2,5) = 1;
end
I_maxima_perinodal_sorted = sort(I_maxima_perinodal);
I_maxima_perinodal_dev_avg_min_1SD = maxima(length(I_maxima_perinodal_sorted),1);
for i = 1:length(I_maxima_perinodal_sorted)
    I_maxima_perinodal_dev_avg_min_1SD(i) = I_maxima_perinodal_sorted(i) - (I_maxima_perinodal_average-I_maxima_perinodal_SD);
end
I_maxima_perinodal_dev_avg_min_1SD = abs(I_maxima_perinodal_dev_avg_min_1SD);
for i = 1:length(I_maxima_perinodal_sorted)
    if I_maxima_perinodal_dev_avg_min_1SD(i) == min(I_maxima_perinodal_dev_avg_min_1SD)
        I_maxima_perinodal_dev_avg_min_1SD_index = i;
    end
end
I_maxima_perinodal_dev_avg_pls_1SD = maxima(length(I_maxima_perinodal_sorted),1);
for i = 1:length(I_maxima_perinodal_sorted)
    I_maxima_perinodal_dev_avg_pls_1SD(i) = I_maxima_perinodal_sorted(i) - (I_maxima_perinodal_average+I_maxima_perinodal_SD);
end
I_maxima_perinodal_dev_avg_pls_1SD = abs(I_maxima_perinodal_dev_avg_pls_1SD);
for i = 1:length(I_maxima_perinodal_sorted)
    if I_maxima_perinodal_dev_avg_pls_1SD(i) == min(I_maxima_perinodal_dev_avg_pls_1SD)
        I_maxima_perinodal_dev_avg_pls_1SD_index = i;
    end
end
I_maxima_perinodal_dev_avg_min_2SD = maxima(length(I_maxima_perinodal_sorted),1);
for i = 1:length(I_maxima_perinodal_sorted)
    I_maxima_perinodal_dev_avg_min_2SD(i) = I_maxima_perinodal_sorted(i) - (I_maxima_perinodal_average-2*I_maxima_perinodal_SD);
end
I_maxima_perinodal_dev_avg_min_2SD = abs(I_maxima_perinodal_dev_avg_min_2SD);
for i = 1:length(I_maxima_perinodal_sorted)
    if I_maxima_perinodal_dev_avg_min_2SD(i) == min(I_maxima_perinodal_dev_avg_min_2SD)
        I_maxima_dev_avg_min_2SD_index = i;
    end
end
I_maxima_perinodal_dev_avg_pls_2SD = maxima(length(I_maxima_perinodal_sorted),1);
for i = 1:length(I_maxima_perinodal_sorted)
    I_maxima_perinodal_dev_avg_pls_2SD(i) = I_maxima_perinodal_sorted(i) - (I_maxima_perinodal_average+2*I_maxima_perinodal_SD);
end
I_maxima_perinodal_dev_avg_pls_2SD = abs(I_maxima_perinodal_dev_avg_pls_2SD);
for i = 1:length(I_maxima_perinodal_sorted)
    if I_maxima_perinodal_dev_avg_pls_2SD(i) == min(I_maxima_perinodal_dev_avg_pls_2SD)
        I_maxima_perinodal_dev_avg_pls_2SD_index = i;
    end
end
I_maxima_perinodal_dev_avg_min_3SD = maxima(length(I_maxima_perinodal_sorted),1);
for i = 1:length(I_maxima_perinodal_sorted)
    I_maxima_perinodal_dev_avg_min_3SD(i) = I_maxima_perinodal_sorted(i) - (I_maxima_perinodal_average-3*I_maxima_perinodal_SD);
end
I_maxima_perinodal_dev_avg_min_3SD = abs(I_maxima_perinodal_dev_avg_min_3SD);
for i = 1:length(I_maxima_perinodal_sorted)
    if I_maxima_perinodal_dev_avg_min_3SD(i) == min(I_maxima_perinodal_dev_avg_min_3SD)
        I_maxima_perinodal_dev_avg_min_3SD_index = i;
    end
end
I_maxima_perinodal_dev_avg_pls_3SD = maxima(length(I_maxima_perinodal_sorted),1);
for i = 1:length(I_maxima_perinodal_sorted)
    I_maxima_perinodal_dev_avg_pls_3SD(i) = I_maxima_perinodal_sorted(i) - (I_maxima_perinodal_average+3*I_maxima_perinodal_SD);
end
I_maxima_perinodal_dev_avg_pls_3SD = abs(I_maxima_perinodal_dev_avg_pls_3SD);
for i = 1:length(I_maxima_perinodal_sorted)
    if I_maxima_perinodal_dev_avg_pls_3SD(i) == min(I_maxima_perinodal_dev_avg_pls_3SD)
        I_maxima_perinodal_dev_avg_pls_3SD_index = i;
    end
end
I_maxima_perinodal_tot_number_of_indices = length(I_maxima_perinodal_sorted);
%This is the approximation of the tot "area under curve." 
I_maxima_perinodal_quotient_within_1SD = (I_maxima_perinodal_dev_avg_pls_1SD_index - I_maxima_perinodaldev_avg_min_1SD_index)/I_maxima_perinodal_tot_number_of_indices;
I_maxima_perinodal_quotient_within_2SD = (I_maxima_perinodal_dev_avg_pls_2SD_index - I_maxima_dev_avg_min_2SD_index)/I_maxima_perinodal_tot_number_of_indices;
I_maxima_perinodal_quotient_within_3SD = (I_maxima_perinodal_dev_avg_pls_3SD_index - I_maxima_perinodal_dev_avg_min_3SD_index)/I_maxima_perinodal_tot_number_of_indices;
if (I_maxima_perinodal_quotient_within_1SD > 0.66 && I_maxima_perinodal_quotient_within_1SD < 0.70) && (I_maxima_perinodal_quotient_within_2SD > 0.93 && I_maxima_perinodal_quotient_within_2SD < 0.97) && (I_maxima_perinodal_quotient_within_3SD > 0.98 && I_maxima_perinodal_quotient_within_3SD < 1)
    isNormal(2,6) = 1;
end
if I_maxima_perinodal_quotient_within_1SD == 0
    isSpike(2,6) = 1;
end
I_maxima_right_septum_sorted = sort(I_maxima_right_septum);
I_maxima_right_septum_dev_avg_min_1SD = maxima(length(I_maxima_right_septum_sorted),1);
for i = 1:length(I_maxima_right_septum_sorted)
    I_maxima_right_septum_dev_avg_min_1SD(i) = I_maxima_right_septum_sorted(i) - (I_maxima_right_septum_average-I_maxima_right_septum_SD);
end
I_maxima_right_septum_dev_avg_min_1SD = abs(I_maxima_right_septum_dev_avg_min_1SD);
for i = 1:length(I_maxima_right_septum_sorted)
    if I_maxima_right_septum_dev_avg_min_1SD(i) == min(I_maxima_right_septum_dev_avg_min_1SD)
        I_maxima_right_septum_dev_avg_min_1SD_index = i;
    end
end
I_maxima_right_septum_dev_avg_pls_1SD = maxima(length(I_maxima_right_septum_sorted),1);
for i = 1:length(I_maxima_right_septum_sorted)
    I_maxima_right_septum_dev_avg_pls_1SD(i) = I_maxima_right_septum_sorted(i) - (I_maxima_right_septum_average+I_maxima_right_septum_SD);
end
I_maxima_right_septum_dev_avg_pls_1SD = abs(I_maxima_right_septum_dev_avg_pls_1SD);
for i = 1:length(I_maxima_right_septum_sorted)
    if I_maxima_right_septum_dev_avg_pls_1SD(i) == min(I_maxima_right_septum_dev_avg_pls_1SD)
        I_maxima_right_septum_dev_avg_pls_1SD_index = i;
    end
end
I_maxima_right_septum_dev_avg_min_2SD = maxima(length(I_maxima_right_septum_sorted),1);
for i = 1:length(I_maxima_right_septum_sorted)
    I_maxima_right_septum_dev_avg_min_2SD(i) = I_maxima_right_septum_sorted(i) - (I_maxima_right_septum_average-2*I_maxima_right_septum_SD);
end
I_maxima_right_septum_dev_avg_min_2SD = abs(I_maxima_right_septum_dev_avg_min_2SD);
for i = 1:length(I_maxima_right_septum_sorted)
    if I_maxima_right_septum_dev_avg_min_2SD(i) == min(I_maxima_right_septum_dev_avg_min_2SD)
        I_maxima_dev_avg_min_2SD_index = i;
    end
end
I_maxima_right_septum_dev_avg_pls_2SD = maxima(length(I_maxima_right_septum_sorted),1);
for i = 1:length(I_maxima_right_septum_sorted)
    I_maxima_right_septum_dev_avg_pls_2SD(i) = I_maxima_right_septum_sorted(i) - (I_maxima_right_septum_average+2*I_maxima_right_septum_SD);
end
I_maxima_right_septum_dev_avg_pls_2SD = abs(I_maxima_right_septum_dev_avg_pls_2SD);
for i = 1:length(I_maxima_right_septum_sorted)
    if I_maxima_right_septum_dev_avg_pls_2SD(i) == min(I_maxima_right_septum_dev_avg_pls_2SD)
        I_maxima_right_septum_dev_avg_pls_2SD_index = i;
    end
end
I_maxima_right_septum_dev_avg_min_3SD = maxima(length(I_maxima_right_septum_sorted),1);
for i = 1:length(I_maxima_right_septum_sorted)
    I_maxima_right_septum_dev_avg_min_3SD(i) = I_maxima_right_septum_sorted(i) - (I_maxima_right_septum_average-3*I_maxima_right_septum_SD);
end
I_maxima_right_septum_dev_avg_min_3SD = abs(I_maxima_right_septum_dev_avg_min_3SD);
for i = 1:length(I_maxima_right_septum_sorted)
    if I_maxima_right_septum_dev_avg_min_3SD(i) == min(I_maxima_right_septum_dev_avg_min_3SD)
        I_maxima_right_septum_dev_avg_min_3SD_index = i;
    end
end
I_maxima_right_septum_dev_avg_pls_3SD = maxima(length(I_maxima_right_septum_sorted),1);
for i = 1:length(I_maxima_right_septum_sorted)
    I_maxima_right_septum_dev_avg_pls_3SD(i) = I_maxima_right_septum_sorted(i) - (I_maxima_right_septum_average+3*I_maxima_right_septum_SD);
end
I_maxima_right_septum_dev_avg_pls_3SD = abs(I_maxima_right_septum_dev_avg_pls_3SD);
for i = 1:length(I_maxima_right_septum_sorted)
    if I_maxima_right_septum_dev_avg_pls_3SD(i) == min(I_maxima_right_septum_dev_avg_pls_3SD)
        I_maxima_right_septum_dev_avg_pls_3SD_index = i;
    end
end
I_maxima_right_septum_tot_number_of_indices = length(I_maxima_right_septum_sorted);
%This is the approximation of the tot "area under curve." 
I_maxima_right_septum_quotient_within_1SD = (I_maxima_right_septum_dev_avg_pls_1SD_index - I_maxima_right_septumdev_avg_min_1SD_index)/I_maxima_right_septum_tot_number_of_indices;
I_maxima_right_septum_quotient_within_2SD = (I_maxima_right_septum_dev_avg_pls_2SD_index - I_maxima_dev_avg_min_2SD_index)/I_maxima_right_septum_tot_number_of_indices;
I_maxima_right_septum_quotient_within_3SD = (I_maxima_right_septum_dev_avg_pls_3SD_index - I_maxima_right_septum_dev_avg_min_3SD_index)/I_maxima_right_septum_tot_number_of_indices;
if (I_maxima_right_septum_quotient_within_1SD > 0.66 && I_maxima_right_septum_quotient_within_1SD < 0.70) && (I_maxima_right_septum_quotient_within_2SD > 0.93 && I_maxima_right_septum_quotient_within_2SD < 0.97) && (I_maxima_right_septum_quotient_within_3SD > 0.98 && I_maxima_right_septum_quotient_within_3SD < 1)
    isNormal(2,7) = 1;
end
if I_maxima_right_septum_quotient_within_1SD == 0
    isSpike(2,7) = 1;
end
I_maxima_right_atrial_appendage_sorted = sort(I_maxima_right_atrial_appendage);
I_maxima_right_atrial_appendage_dev_avg_min_1SD = maxima(length(I_maxima_right_atrial_appendage_sorted),1);
for i = 1:length(I_maxima_right_atrial_appendage_sorted)
    I_maxima_right_atrial_appendage_dev_avg_min_1SD(i) = I_maxima_right_atrial_appendage_sorted(i) - (I_maxima_right_atrial_appendage_average-I_maxima_right_atrial_appendage_SD);
end
I_maxima_right_atrial_appendage_dev_avg_min_1SD = abs(I_maxima_right_atrial_appendage_dev_avg_min_1SD);
for i = 1:length(I_maxima_right_atrial_appendage_sorted)
    if I_maxima_right_atrial_appendage_dev_avg_min_1SD(i) == min(I_maxima_right_atrial_appendage_dev_avg_min_1SD)
        I_maxima_right_atrial_appendage_dev_avg_min_1SD_index = i;
    end
end
I_maxima_right_atrial_appendage_dev_avg_pls_1SD = maxima(length(I_maxima_right_atrial_appendage_sorted),1);
for i = 1:length(I_maxima_right_atrial_appendage_sorted)
    I_maxima_right_atrial_appendage_dev_avg_pls_1SD(i) = I_maxima_right_atrial_appendage_sorted(i) - (I_maxima_right_atrial_appendage_average+I_maxima_right_atrial_appendage_SD);
end
I_maxima_right_atrial_appendage_dev_avg_pls_1SD = abs(I_maxima_right_atrial_appendage_dev_avg_pls_1SD);
for i = 1:length(I_maxima_right_atrial_appendage_sorted)
    if I_maxima_right_atrial_appendage_dev_avg_pls_1SD(i) == min(I_maxima_right_atrial_appendage_dev_avg_pls_1SD)
        I_maxima_right_atrial_appendage_dev_avg_pls_1SD_index = i;
    end
end
I_maxima_right_atrial_appendage_dev_avg_min_2SD = maxima(length(I_maxima_right_atrial_appendage_sorted),1);
for i = 1:length(I_maxima_right_atrial_appendage_sorted)
    I_maxima_right_atrial_appendage_dev_avg_min_2SD(i) = I_maxima_right_atrial_appendage_sorted(i) - (I_maxima_right_atrial_appendage_average-2*I_maxima_right_atrial_appendage_SD);
end
I_maxima_right_atrial_appendage_dev_avg_min_2SD = abs(I_maxima_right_atrial_appendage_dev_avg_min_2SD);
for i = 1:length(I_maxima_right_atrial_appendage_sorted)
    if I_maxima_right_atrial_appendage_dev_avg_min_2SD(i) == min(I_maxima_right_atrial_appendage_dev_avg_min_2SD)
        I_maxima_dev_avg_min_2SD_index = i;
    end
end
I_maxima_right_atrial_appendage_dev_avg_pls_2SD = maxima(length(I_maxima_right_atrial_appendage_sorted),1);
for i = 1:length(I_maxima_right_atrial_appendage_sorted)
    I_maxima_right_atrial_appendage_dev_avg_pls_2SD(i) = I_maxima_right_atrial_appendage_sorted(i) - (I_maxima_right_atrial_appendage_average+2*I_maxima_right_atrial_appendage_SD);
end
I_maxima_right_atrial_appendage_dev_avg_pls_2SD = abs(I_maxima_right_atrial_appendage_dev_avg_pls_2SD);
for i = 1:length(I_maxima_right_atrial_appendage_sorted)
    if I_maxima_right_atrial_appendage_dev_avg_pls_2SD(i) == min(I_maxima_right_atrial_appendage_dev_avg_pls_2SD)
        I_maxima_right_atrial_appendage_dev_avg_pls_2SD_index = i;
    end
end
I_maxima_right_atrial_appendage_dev_avg_min_3SD = maxima(length(I_maxima_right_atrial_appendage_sorted),1);
for i = 1:length(I_maxima_right_atrial_appendage_sorted)
    I_maxima_right_atrial_appendage_dev_avg_min_3SD(i) = I_maxima_right_atrial_appendage_sorted(i) - (I_maxima_right_atrial_appendage_average-3*I_maxima_right_atrial_appendage_SD);
end
I_maxima_right_atrial_appendage_dev_avg_min_3SD = abs(I_maxima_right_atrial_appendage_dev_avg_min_3SD);
for i = 1:length(I_maxima_right_atrial_appendage_sorted)
    if I_maxima_right_atrial_appendage_dev_avg_min_3SD(i) == min(I_maxima_right_atrial_appendage_dev_avg_min_3SD)
        I_maxima_right_atrial_appendage_dev_avg_min_3SD_index = i;
    end
end
I_maxima_right_atrial_appendage_dev_avg_pls_3SD = maxima(length(I_maxima_right_atrial_appendage_sorted),1);
for i = 1:length(I_maxima_right_atrial_appendage_sorted)
    I_maxima_right_atrial_appendage_dev_avg_pls_3SD(i) = I_maxima_right_atrial_appendage_sorted(i) - (I_maxima_right_atrial_appendage_average+3*I_maxima_right_atrial_appendage_SD);
end
I_maxima_right_atrial_appendage_dev_avg_pls_3SD = abs(I_maxima_right_atrial_appendage_dev_avg_pls_3SD);
for i = 1:length(I_maxima_right_atrial_appendage_sorted)
    if I_maxima_right_atrial_appendage_dev_avg_pls_3SD(i) == min(I_maxima_right_atrial_appendage_dev_avg_pls_3SD)
        I_maxima_right_atrial_appendage_dev_avg_pls_3SD_index = i;
    end
end
I_maxima_right_atrial_appendage_tot_number_of_indices = length(I_maxima_right_atrial_appendage_sorted);
%This is the approximation of the tot "area under curve." 
I_maxima_right_atrial_appendage_quotient_within_1SD = (I_maxima_right_atrial_appendage_dev_avg_pls_1SD_index - I_maxima_right_atrial_appendagedev_avg_min_1SD_index)/I_maxima_right_atrial_appendage_tot_number_of_indices;
I_maxima_right_atrial_appendage_quotient_within_2SD = (I_maxima_right_atrial_appendage_dev_avg_pls_2SD_index - I_maxima_dev_avg_min_2SD_index)/I_maxima_right_atrial_appendage_tot_number_of_indices;
I_maxima_right_atrial_appendage_quotient_within_3SD = (I_maxima_right_atrial_appendage_dev_avg_pls_3SD_index - I_maxima_right_atrial_appendage_dev_avg_min_3SD_index)/I_maxima_right_atrial_appendage_tot_number_of_indices;
if (I_maxima_right_atrial_appendage_quotient_within_1SD > 0.66 && I_maxima_right_atrial_appendage_quotient_within_1SD < 0.70) && (I_maxima_right_atrial_appendage_quotient_within_2SD > 0.93 && I_maxima_right_atrial_appendage_quotient_within_2SD < 0.97) && (I_maxima_right_atrial_appendage_quotient_within_3SD > 0.98 && I_maxima_right_atrial_appendage_quotient_within_3SD < 1)
    isNormal(2,8) = 1;
end
if I_maxima_right_atrial_appendage_quotient_within_1SD == 0
    isSpike(2,8) = 1;
end
I_maxima_pulmonary_veins_sorted = sort(I_maxima_pulmonary_veins);
I_maxima_pulmonary_veins_dev_avg_min_1SD = maxima(length(I_maxima_pulmonary_veins_sorted),1);
for i = 1:length(I_maxima_pulmonary_veins_sorted)
    I_maxima_pulmonary_veins_dev_avg_min_1SD(i) = I_maxima_pulmonary_veins_sorted(i) - (I_maxima_pulmonary_veins_average-I_maxima_pulmonary_veins_SD);
end
I_maxima_pulmonary_veins_dev_avg_min_1SD = abs(I_maxima_pulmonary_veins_dev_avg_min_1SD);
for i = 1:length(I_maxima_pulmonary_veins_sorted)
    if I_maxima_pulmonary_veins_dev_avg_min_1SD(i) == min(I_maxima_pulmonary_veins_dev_avg_min_1SD)
        I_maxima_pulmonary_veins_dev_avg_min_1SD_index = i;
    end
end
I_maxima_pulmonary_veins_dev_avg_pls_1SD = maxima(length(I_maxima_pulmonary_veins_sorted),1);
for i = 1:length(I_maxima_pulmonary_veins_sorted)
    I_maxima_pulmonary_veins_dev_avg_pls_1SD(i) = I_maxima_pulmonary_veins_sorted(i) - (I_maxima_pulmonary_veins_average+I_maxima_pulmonary_veins_SD);
end
I_maxima_pulmonary_veins_dev_avg_pls_1SD = abs(I_maxima_pulmonary_veins_dev_avg_pls_1SD);
for i = 1:length(I_maxima_pulmonary_veins_sorted)
    if I_maxima_pulmonary_veins_dev_avg_pls_1SD(i) == min(I_maxima_pulmonary_veins_dev_avg_pls_1SD)
        I_maxima_pulmonary_veins_dev_avg_pls_1SD_index = i;
    end
end
I_maxima_pulmonary_veins_dev_avg_min_2SD = maxima(length(I_maxima_pulmonary_veins_sorted),1);
for i = 1:length(I_maxima_pulmonary_veins_sorted)
    I_maxima_pulmonary_veins_dev_avg_min_2SD(i) = I_maxima_pulmonary_veins_sorted(i) - (I_maxima_pulmonary_veins_average-2*I_maxima_pulmonary_veins_SD);
end
I_maxima_pulmonary_veins_dev_avg_min_2SD = abs(I_maxima_pulmonary_veins_dev_avg_min_2SD);
for i = 1:length(I_maxima_pulmonary_veins_sorted)
    if I_maxima_pulmonary_veins_dev_avg_min_2SD(i) == min(I_maxima_pulmonary_veins_dev_avg_min_2SD)
        I_maxima_dev_avg_min_2SD_index = i;
    end
end
I_maxima_pulmonary_veins_dev_avg_pls_2SD = maxima(length(I_maxima_pulmonary_veins_sorted),1);
for i = 1:length(I_maxima_pulmonary_veins_sorted)
    I_maxima_pulmonary_veins_dev_avg_pls_2SD(i) = I_maxima_pulmonary_veins_sorted(i) - (I_maxima_pulmonary_veins_average+2*I_maxima_pulmonary_veins_SD);
end
I_maxima_pulmonary_veins_dev_avg_pls_2SD = abs(I_maxima_pulmonary_veins_dev_avg_pls_2SD);
for i = 1:length(I_maxima_pulmonary_veins_sorted)
    if I_maxima_pulmonary_veins_dev_avg_pls_2SD(i) == min(I_maxima_pulmonary_veins_dev_avg_pls_2SD)
        I_maxima_pulmonary_veins_dev_avg_pls_2SD_index = i;
    end
end
I_maxima_pulmonary_veins_dev_avg_min_3SD = maxima(length(I_maxima_pulmonary_veins_sorted),1);
for i = 1:length(I_maxima_pulmonary_veins_sorted)
    I_maxima_pulmonary_veins_dev_avg_min_3SD(i) = I_maxima_pulmonary_veins_sorted(i) - (I_maxima_pulmonary_veins_average-3*I_maxima_pulmonary_veins_SD);
end
I_maxima_pulmonary_veins_dev_avg_min_3SD = abs(I_maxima_pulmonary_veins_dev_avg_min_3SD);
for i = 1:length(I_maxima_pulmonary_veins_sorted)
    if I_maxima_pulmonary_veins_dev_avg_min_3SD(i) == min(I_maxima_pulmonary_veins_dev_avg_min_3SD)
        I_maxima_pulmonary_veins_dev_avg_min_3SD_index = i;
    end
end
I_maxima_pulmonary_veins_dev_avg_pls_3SD = maxima(length(I_maxima_pulmonary_veins_sorted),1);
for i = 1:length(I_maxima_pulmonary_veins_sorted)
    I_maxima_pulmonary_veins_dev_avg_pls_3SD(i) = I_maxima_pulmonary_veins_sorted(i) - (I_maxima_pulmonary_veins_average+3*I_maxima_pulmonary_veins_SD);
end
I_maxima_pulmonary_veins_dev_avg_pls_3SD = abs(I_maxima_pulmonary_veins_dev_avg_pls_3SD);
for i = 1:length(I_maxima_pulmonary_veins_sorted)
    if I_maxima_pulmonary_veins_dev_avg_pls_3SD(i) == min(I_maxima_pulmonary_veins_dev_avg_pls_3SD)
        I_maxima_pulmonary_veins_dev_avg_pls_3SD_index = i;
    end
end
I_maxima_pulmonary_veins_tot_number_of_indices = length(I_maxima_pulmonary_veins_sorted);
%This is the approximation of the tot "area under curve." 
I_maxima_pulmonary_veins_quotient_within_1SD = (I_maxima_pulmonary_veins_dev_avg_pls_1SD_index - I_maxima_pulmonary_veinsdev_avg_min_1SD_index)/I_maxima_pulmonary_veins_tot_number_of_indices;
I_maxima_pulmonary_veins_quotient_within_2SD = (I_maxima_pulmonary_veins_dev_avg_pls_2SD_index - I_maxima_dev_avg_min_2SD_index)/I_maxima_pulmonary_veins_tot_number_of_indices;
I_maxima_pulmonary_veins_quotient_within_3SD = (I_maxima_pulmonary_veins_dev_avg_pls_3SD_index - I_maxima_pulmonary_veins_dev_avg_min_3SD_index)/I_maxima_pulmonary_veins_tot_number_of_indices;
if (I_maxima_pulmonary_veins_quotient_within_1SD > 0.66 && I_maxima_pulmonary_veins_quotient_within_1SD < 0.70) && (I_maxima_pulmonary_veins_quotient_within_2SD > 0.93 && I_maxima_pulmonary_veins_quotient_within_2SD < 0.97) && (I_maxima_pulmonary_veins_quotient_within_3SD > 0.98 && I_maxima_pulmonary_veins_quotient_within_3SD < 1)
    isNormal(2,9) = 1;
end
if I_maxima_pulmonary_veins_quotient_within_1SD == 0
    isSpike(2,9) = 1;
end
I_maxima_mitral_annulus_sorted = sort(I_maxima_mitral_annulus);
I_maxima_mitral_annulus_dev_avg_min_1SD = maxima(length(I_maxima_mitral_annulus_sorted),1);
for i = 1:length(I_maxima_mitral_annulus_sorted)
    I_maxima_mitral_annulus_dev_avg_min_1SD(i) = I_maxima_mitral_annulus_sorted(i) - (I_maxima_mitral_annulus_average-I_maxima_mitral_annulus_SD);
end
I_maxima_mitral_annulus_dev_avg_min_1SD = abs(I_maxima_mitral_annulus_dev_avg_min_1SD);
for i = 1:length(I_maxima_mitral_annulus_sorted)
    if I_maxima_mitral_annulus_dev_avg_min_1SD(i) == min(I_maxima_mitral_annulus_dev_avg_min_1SD)
        I_maxima_mitral_annulus_dev_avg_min_1SD_index = i;
    end
end
I_maxima_mitral_annulus_dev_avg_pls_1SD = maxima(length(I_maxima_mitral_annulus_sorted),1);
for i = 1:length(I_maxima_mitral_annulus_sorted)
    I_maxima_mitral_annulus_dev_avg_pls_1SD(i) = I_maxima_mitral_annulus_sorted(i) - (I_maxima_mitral_annulus_average+I_maxima_mitral_annulus_SD);
end
I_maxima_mitral_annulus_dev_avg_pls_1SD = abs(I_maxima_mitral_annulus_dev_avg_pls_1SD);
for i = 1:length(I_maxima_mitral_annulus_sorted)
    if I_maxima_mitral_annulus_dev_avg_pls_1SD(i) == min(I_maxima_mitral_annulus_dev_avg_pls_1SD)
        I_maxima_mitral_annulus_dev_avg_pls_1SD_index = i;
    end
end
I_maxima_mitral_annulus_dev_avg_min_2SD = maxima(length(I_maxima_mitral_annulus_sorted),1);
for i = 1:length(I_maxima_mitral_annulus_sorted)
    I_maxima_mitral_annulus_dev_avg_min_2SD(i) = I_maxima_mitral_annulus_sorted(i) - (I_maxima_mitral_annulus_average-2*I_maxima_mitral_annulus_SD);
end
I_maxima_mitral_annulus_dev_avg_min_2SD = abs(I_maxima_mitral_annulus_dev_avg_min_2SD);
for i = 1:length(I_maxima_mitral_annulus_sorted)
    if I_maxima_mitral_annulus_dev_avg_min_2SD(i) == min(I_maxima_mitral_annulus_dev_avg_min_2SD)
        I_maxima_dev_avg_min_2SD_index = i;
    end
end
I_maxima_mitral_annulus_dev_avg_pls_2SD = maxima(length(I_maxima_mitral_annulus_sorted),1);
for i = 1:length(I_maxima_mitral_annulus_sorted)
    I_maxima_mitral_annulus_dev_avg_pls_2SD(i) = I_maxima_mitral_annulus_sorted(i) - (I_maxima_mitral_annulus_average+2*I_maxima_mitral_annulus_SD);
end
I_maxima_mitral_annulus_dev_avg_pls_2SD = abs(I_maxima_mitral_annulus_dev_avg_pls_2SD);
for i = 1:length(I_maxima_mitral_annulus_sorted)
    if I_maxima_mitral_annulus_dev_avg_pls_2SD(i) == min(I_maxima_mitral_annulus_dev_avg_pls_2SD)
        I_maxima_mitral_annulus_dev_avg_pls_2SD_index = i;
    end
end
I_maxima_mitral_annulus_dev_avg_min_3SD = maxima(length(I_maxima_mitral_annulus_sorted),1);
for i = 1:length(I_maxima_mitral_annulus_sorted)
    I_maxima_mitral_annulus_dev_avg_min_3SD(i) = I_maxima_mitral_annulus_sorted(i) - (I_maxima_mitral_annulus_average-3*I_maxima_mitral_annulus_SD);
end
I_maxima_mitral_annulus_dev_avg_min_3SD = abs(I_maxima_mitral_annulus_dev_avg_min_3SD);
for i = 1:length(I_maxima_mitral_annulus_sorted)
    if I_maxima_mitral_annulus_dev_avg_min_3SD(i) == min(I_maxima_mitral_annulus_dev_avg_min_3SD)
        I_maxima_mitral_annulus_dev_avg_min_3SD_index = i;
    end
end
I_maxima_mitral_annulus_dev_avg_pls_3SD = maxima(length(I_maxima_mitral_annulus_sorted),1);
for i = 1:length(I_maxima_mitral_annulus_sorted)
    I_maxima_mitral_annulus_dev_avg_pls_3SD(i) = I_maxima_mitral_annulus_sorted(i) - (I_maxima_mitral_annulus_average+3*I_maxima_mitral_annulus_SD);
end
I_maxima_mitral_annulus_dev_avg_pls_3SD = abs(I_maxima_mitral_annulus_dev_avg_pls_3SD);
for i = 1:length(I_maxima_mitral_annulus_sorted)
    if I_maxima_mitral_annulus_dev_avg_pls_3SD(i) == min(I_maxima_mitral_annulus_dev_avg_pls_3SD)
        I_maxima_mitral_annulus_dev_avg_pls_3SD_index = i;
    end
end
I_maxima_mitral_annulus_tot_number_of_indices = length(I_maxima_mitral_annulus_sorted);
%This is the approximation of the tot "area under curve." 
I_maxima_mitral_annulus_quotient_within_1SD = (I_maxima_mitral_annulus_dev_avg_pls_1SD_index - I_maxima_mitral_annulusdev_avg_min_1SD_index)/I_maxima_mitral_annulus_tot_number_of_indices;
I_maxima_mitral_annulus_quotient_within_2SD = (I_maxima_mitral_annulus_dev_avg_pls_2SD_index - I_maxima_dev_avg_min_2SD_index)/I_maxima_mitral_annulus_tot_number_of_indices;
I_maxima_mitral_annulus_quotient_within_3SD = (I_maxima_mitral_annulus_dev_avg_pls_3SD_index - I_maxima_mitral_annulus_dev_avg_min_3SD_index)/I_maxima_mitral_annulus_tot_number_of_indices;
if (I_maxima_mitral_annulus_quotient_within_1SD > 0.66 && I_maxima_mitral_annulus_quotient_within_1SD < 0.70) && (I_maxima_mitral_annulus_quotient_within_2SD > 0.93 && I_maxima_mitral_annulus_quotient_within_2SD < 0.97) && (I_maxima_mitral_annulus_quotient_within_3SD > 0.98 && I_maxima_mitral_annulus_quotient_within_3SD < 1)
    isNormal(2,10) = 1;
end
if I_maxima_mitral_annulus_quotient_within_1SD == 0
    isSpike(2,10) = 1;
end
I_maxima_CS_body_sorted = sort(I_maxima_CS_body);
I_maxima_CS_body_dev_avg_min_1SD = maxima(length(I_maxima_CS_body_sorted),1);
for i = 1:length(I_maxima_CS_body_sorted)
    I_maxima_CS_body_dev_avg_min_1SD(i) = I_maxima_CS_body_sorted(i) - (I_maxima_CS_body_average-I_maxima_CS_body_SD);
end
I_maxima_CS_body_dev_avg_min_1SD = abs(I_maxima_CS_body_dev_avg_min_1SD);
for i = 1:length(I_maxima_CS_body_sorted)
    if I_maxima_CS_body_dev_avg_min_1SD(i) == min(I_maxima_CS_body_dev_avg_min_1SD)
        I_maxima_CS_body_dev_avg_min_1SD_index = i;
    end
end
I_maxima_CS_body_dev_avg_pls_1SD = maxima(length(I_maxima_CS_body_sorted),1);
for i = 1:length(I_maxima_CS_body_sorted)
    I_maxima_CS_body_dev_avg_pls_1SD(i) = I_maxima_CS_body_sorted(i) - (I_maxima_CS_body_average+I_maxima_CS_body_SD);
end
I_maxima_CS_body_dev_avg_pls_1SD = abs(I_maxima_CS_body_dev_avg_pls_1SD);
for i = 1:length(I_maxima_CS_body_sorted)
    if I_maxima_CS_body_dev_avg_pls_1SD(i) == min(I_maxima_CS_body_dev_avg_pls_1SD)
        I_maxima_CS_body_dev_avg_pls_1SD_index = i;
    end
end
I_maxima_CS_body_dev_avg_min_2SD = maxima(length(I_maxima_CS_body_sorted),1);
for i = 1:length(I_maxima_CS_body_sorted)
    I_maxima_CS_body_dev_avg_min_2SD(i) = I_maxima_CS_body_sorted(i) - (I_maxima_CS_body_average-2*I_maxima_CS_body_SD);
end
I_maxima_CS_body_dev_avg_min_2SD = abs(I_maxima_CS_body_dev_avg_min_2SD);
for i = 1:length(I_maxima_CS_body_sorted)
    if I_maxima_CS_body_dev_avg_min_2SD(i) == min(I_maxima_CS_body_dev_avg_min_2SD)
        I_maxima_dev_avg_min_2SD_index = i;
    end
end
I_maxima_CS_body_dev_avg_pls_2SD = maxima(length(I_maxima_CS_body_sorted),1);
for i = 1:length(I_maxima_CS_body_sorted)
    I_maxima_CS_body_dev_avg_pls_2SD(i) = I_maxima_CS_body_sorted(i) - (I_maxima_CS_body_average+2*I_maxima_CS_body_SD);
end
I_maxima_CS_body_dev_avg_pls_2SD = abs(I_maxima_CS_body_dev_avg_pls_2SD);
for i = 1:length(I_maxima_CS_body_sorted)
    if I_maxima_CS_body_dev_avg_pls_2SD(i) == min(I_maxima_CS_body_dev_avg_pls_2SD)
        I_maxima_CS_body_dev_avg_pls_2SD_index = i;
    end
end
I_maxima_CS_body_dev_avg_min_3SD = maxima(length(I_maxima_CS_body_sorted),1);
for i = 1:length(I_maxima_CS_body_sorted)
    I_maxima_CS_body_dev_avg_min_3SD(i) = I_maxima_CS_body_sorted(i) - (I_maxima_CS_body_average-3*I_maxima_CS_body_SD);
end
I_maxima_CS_body_dev_avg_min_3SD = abs(I_maxima_CS_body_dev_avg_min_3SD);
for i = 1:length(I_maxima_CS_body_sorted)
    if I_maxima_CS_body_dev_avg_min_3SD(i) == min(I_maxima_CS_body_dev_avg_min_3SD)
        I_maxima_CS_body_dev_avg_min_3SD_index = i;
    end
end
I_maxima_CS_body_dev_avg_pls_3SD = maxima(length(I_maxima_CS_body_sorted),1);
for i = 1:length(I_maxima_CS_body_sorted)
    I_maxima_CS_body_dev_avg_pls_3SD(i) = I_maxima_CS_body_sorted(i) - (I_maxima_CS_body_average+3*I_maxima_CS_body_SD);
end
I_maxima_CS_body_dev_avg_pls_3SD = abs(I_maxima_CS_body_dev_avg_pls_3SD);
for i = 1:length(I_maxima_CS_body_sorted)
    if I_maxima_CS_body_dev_avg_pls_3SD(i) == min(I_maxima_CS_body_dev_avg_pls_3SD)
        I_maxima_CS_body_dev_avg_pls_3SD_index = i;
    end
end
I_maxima_CS_body_tot_number_of_indices = length(I_maxima_CS_body_sorted);
%This is the approximation of the tot "area under curve." 
I_maxima_CS_body_quotient_within_1SD = (I_maxima_CS_body_dev_avg_pls_1SD_index - I_maxima_CS_bodydev_avg_min_1SD_index)/I_maxima_CS_body_tot_number_of_indices;
I_maxima_CS_body_quotient_within_2SD = (I_maxima_CS_body_dev_avg_pls_2SD_index - I_maxima_dev_avg_min_2SD_index)/I_maxima_CS_body_tot_number_of_indices;
I_maxima_CS_body_quotient_within_3SD = (I_maxima_CS_body_dev_avg_pls_3SD_index - I_maxima_CS_body_dev_avg_min_3SD_index)/I_maxima_CS_body_tot_number_of_indices;
if (I_maxima_CS_body_quotient_within_1SD > 0.66 && I_maxima_CS_body_quotient_within_1SD < 0.70) && (I_maxima_CS_body_quotient_within_2SD > 0.93 && I_maxima_CS_body_quotient_within_2SD < 0.97) && (I_maxima_CS_body_quotient_within_3SD > 0.98 && I_maxima_CS_body_quotient_within_3SD < 1)
    isNormal(2,11) = 1;
end
if I_maxima_CS_body_quotient_within_1SD == 0
    isSpike(2,11) = 1;
end
I_maxima_left_septum_sorted = sort(I_maxima_left_septum);
I_maxima_left_septum_dev_avg_min_1SD = maxima(length(I_maxima_left_septum_sorted),1);
for i = 1:length(I_maxima_left_septum_sorted)
    I_maxima_left_septum_dev_avg_min_1SD(i) = I_maxima_left_septum_sorted(i) - (I_maxima_left_septum_average-I_maxima_left_septum_SD);
end
I_maxima_left_septum_dev_avg_min_1SD = abs(I_maxima_left_septum_dev_avg_min_1SD);
for i = 1:length(I_maxima_left_septum_sorted)
    if I_maxima_left_septum_dev_avg_min_1SD(i) == min(I_maxima_left_septum_dev_avg_min_1SD)
        I_maxima_left_septum_dev_avg_min_1SD_index = i;
    end
end
I_maxima_left_septum_dev_avg_pls_1SD = maxima(length(I_maxima_left_septum_sorted),1);
for i = 1:length(I_maxima_left_septum_sorted)
    I_maxima_left_septum_dev_avg_pls_1SD(i) = I_maxima_left_septum_sorted(i) - (I_maxima_left_septum_average+I_maxima_left_septum_SD);
end
I_maxima_left_septum_dev_avg_pls_1SD = abs(I_maxima_left_septum_dev_avg_pls_1SD);
for i = 1:length(I_maxima_left_septum_sorted)
    if I_maxima_left_septum_dev_avg_pls_1SD(i) == min(I_maxima_left_septum_dev_avg_pls_1SD)
        I_maxima_left_septum_dev_avg_pls_1SD_index = i;
    end
end
I_maxima_left_septum_dev_avg_min_2SD = maxima(length(I_maxima_left_septum_sorted),1);
for i = 1:length(I_maxima_left_septum_sorted)
    I_maxima_left_septum_dev_avg_min_2SD(i) = I_maxima_left_septum_sorted(i) - (I_maxima_left_septum_average-2*I_maxima_left_septum_SD);
end
I_maxima_left_septum_dev_avg_min_2SD = abs(I_maxima_left_septum_dev_avg_min_2SD);
for i = 1:length(I_maxima_left_septum_sorted)
    if I_maxima_left_septum_dev_avg_min_2SD(i) == min(I_maxima_left_septum_dev_avg_min_2SD)
        I_maxima_dev_avg_min_2SD_index = i;
    end
end
I_maxima_left_septum_dev_avg_pls_2SD = maxima(length(I_maxima_left_septum_sorted),1);
for i = 1:length(I_maxima_left_septum_sorted)
    I_maxima_left_septum_dev_avg_pls_2SD(i) = I_maxima_left_septum_sorted(i) - (I_maxima_left_septum_average+2*I_maxima_left_septum_SD);
end
I_maxima_left_septum_dev_avg_pls_2SD = abs(I_maxima_left_septum_dev_avg_pls_2SD);
for i = 1:length(I_maxima_left_septum_sorted)
    if I_maxima_left_septum_dev_avg_pls_2SD(i) == min(I_maxima_left_septum_dev_avg_pls_2SD)
        I_maxima_left_septum_dev_avg_pls_2SD_index = i;
    end
end
I_maxima_left_septum_dev_avg_min_3SD = maxima(length(I_maxima_left_septum_sorted),1);
for i = 1:length(I_maxima_left_septum_sorted)
    I_maxima_left_septum_dev_avg_min_3SD(i) = I_maxima_left_septum_sorted(i) - (I_maxima_left_septum_average-3*I_maxima_left_septum_SD);
end
I_maxima_left_septum_dev_avg_min_3SD = abs(I_maxima_left_septum_dev_avg_min_3SD);
for i = 1:length(I_maxima_left_septum_sorted)
    if I_maxima_left_septum_dev_avg_min_3SD(i) == min(I_maxima_left_septum_dev_avg_min_3SD)
        I_maxima_left_septum_dev_avg_min_3SD_index = i;
    end
end
I_maxima_left_septum_dev_avg_pls_3SD = maxima(length(I_maxima_left_septum_sorted),1);
for i = 1:length(I_maxima_left_septum_sorted)
    I_maxima_left_septum_dev_avg_pls_3SD(i) = I_maxima_left_septum_sorted(i) - (I_maxima_left_septum_average+3*I_maxima_left_septum_SD);
end
I_maxima_left_septum_dev_avg_pls_3SD = abs(I_maxima_left_septum_dev_avg_pls_3SD);
for i = 1:length(I_maxima_left_septum_sorted)
    if I_maxima_left_septum_dev_avg_pls_3SD(i) == min(I_maxima_left_septum_dev_avg_pls_3SD)
        I_maxima_left_septum_dev_avg_pls_3SD_index = i;
    end
end
I_maxima_left_septum_tot_number_of_indices = length(I_maxima_left_septum_sorted);
%This is the approximation of the tot "area under curve." 
I_maxima_left_septum_quotient_within_1SD = (I_maxima_left_septum_dev_avg_pls_1SD_index - I_maxima_left_septumdev_avg_min_1SD_index)/I_maxima_left_septum_tot_number_of_indices;
I_maxima_left_septum_quotient_within_2SD = (I_maxima_left_septum_dev_avg_pls_2SD_index - I_maxima_dev_avg_min_2SD_index)/I_maxima_left_septum_tot_number_of_indices;
I_maxima_left_septum_quotient_within_3SD = (I_maxima_left_septum_dev_avg_pls_3SD_index - I_maxima_left_septum_dev_avg_min_3SD_index)/I_maxima_left_septum_tot_number_of_indices;
if (I_maxima_left_septum_quotient_within_1SD > 0.66 && I_maxima_left_septum_quotient_within_1SD < 0.70) && (I_maxima_left_septum_quotient_within_2SD > 0.93 && I_maxima_left_septum_quotient_within_2SD < 0.97) && (I_maxima_left_septum_quotient_within_3SD > 0.98 && I_maxima_left_septum_quotient_within_3SD < 1)
    isNormal(2,12) = 1;
end
if I_maxima_left_septum_quotient_within_1SD == 0
    isSpike(2,12) = 1;
end
I_maxima_left_atrial_appendage_sorted = sort(I_maxima_left_atrial_appendage);
I_maxima_left_atrial_appendage_dev_avg_min_1SD = maxima(length(I_maxima_left_atrial_appendage_sorted),1);
for i = 1:length(I_maxima_left_atrial_appendage_sorted)
    I_maxima_left_atrial_appendage_dev_avg_min_1SD(i) = I_maxima_left_atrial_appendage_sorted(i) - (I_maxima_left_atrial_appendage_average-I_maxima_left_atrial_appendage_SD);
end
I_maxima_left_atrial_appendage_dev_avg_min_1SD = abs(I_maxima_left_atrial_appendage_dev_avg_min_1SD);
for i = 1:length(I_maxima_left_atrial_appendage_sorted)
    if I_maxima_left_atrial_appendage_dev_avg_min_1SD(i) == min(I_maxima_left_atrial_appendage_dev_avg_min_1SD)
        I_maxima_left_atrial_appendage_dev_avg_min_1SD_index = i;
    end
end
I_maxima_left_atrial_appendage_dev_avg_pls_1SD = maxima(length(I_maxima_left_atrial_appendage_sorted),1);
for i = 1:length(I_maxima_left_atrial_appendage_sorted)
    I_maxima_left_atrial_appendage_dev_avg_pls_1SD(i) = I_maxima_left_atrial_appendage_sorted(i) - (I_maxima_left_atrial_appendage_average+I_maxima_left_atrial_appendage_SD);
end
I_maxima_left_atrial_appendage_dev_avg_pls_1SD = abs(I_maxima_left_atrial_appendage_dev_avg_pls_1SD);
for i = 1:length(I_maxima_left_atrial_appendage_sorted)
    if I_maxima_left_atrial_appendage_dev_avg_pls_1SD(i) == min(I_maxima_left_atrial_appendage_dev_avg_pls_1SD)
        I_maxima_left_atrial_appendage_dev_avg_pls_1SD_index = i;
    end
end
I_maxima_left_atrial_appendage_dev_avg_min_2SD = maxima(length(I_maxima_left_atrial_appendage_sorted),1);
for i = 1:length(I_maxima_left_atrial_appendage_sorted)
    I_maxima_left_atrial_appendage_dev_avg_min_2SD(i) = I_maxima_left_atrial_appendage_sorted(i) - (I_maxima_left_atrial_appendage_average-2*I_maxima_left_atrial_appendage_SD);
end
I_maxima_left_atrial_appendage_dev_avg_min_2SD = abs(I_maxima_left_atrial_appendage_dev_avg_min_2SD);
for i = 1:length(I_maxima_left_atrial_appendage_sorted)
    if I_maxima_left_atrial_appendage_dev_avg_min_2SD(i) == min(I_maxima_left_atrial_appendage_dev_avg_min_2SD)
        I_maxima_dev_avg_min_2SD_index = i;
    end
end
I_maxima_left_atrial_appendage_dev_avg_pls_2SD = maxima(length(I_maxima_left_atrial_appendage_sorted),1);
for i = 1:length(I_maxima_left_atrial_appendage_sorted)
    I_maxima_left_atrial_appendage_dev_avg_pls_2SD(i) = I_maxima_left_atrial_appendage_sorted(i) - (I_maxima_left_atrial_appendage_average+2*I_maxima_left_atrial_appendage_SD);
end
I_maxima_left_atrial_appendage_dev_avg_pls_2SD = abs(I_maxima_left_atrial_appendage_dev_avg_pls_2SD);
for i = 1:length(I_maxima_left_atrial_appendage_sorted)
    if I_maxima_left_atrial_appendage_dev_avg_pls_2SD(i) == min(I_maxima_left_atrial_appendage_dev_avg_pls_2SD)
        I_maxima_left_atrial_appendage_dev_avg_pls_2SD_index = i;
    end
end
I_maxima_left_atrial_appendage_dev_avg_min_3SD = maxima(length(I_maxima_left_atrial_appendage_sorted),1);
for i = 1:length(I_maxima_left_atrial_appendage_sorted)
    I_maxima_left_atrial_appendage_dev_avg_min_3SD(i) = I_maxima_left_atrial_appendage_sorted(i) - (I_maxima_left_atrial_appendage_average-3*I_maxima_left_atrial_appendage_SD);
end
I_maxima_left_atrial_appendage_dev_avg_min_3SD = abs(I_maxima_left_atrial_appendage_dev_avg_min_3SD);
for i = 1:length(I_maxima_left_atrial_appendage_sorted)
    if I_maxima_left_atrial_appendage_dev_avg_min_3SD(i) == min(I_maxima_left_atrial_appendage_dev_avg_min_3SD)
        I_maxima_left_atrial_appendage_dev_avg_min_3SD_index = i;
    end
end
I_maxima_left_atrial_appendage_dev_avg_pls_3SD = maxima(length(I_maxima_left_atrial_appendage_sorted),1);
for i = 1:length(I_maxima_left_atrial_appendage_sorted)
    I_maxima_left_atrial_appendage_dev_avg_pls_3SD(i) = I_maxima_left_atrial_appendage_sorted(i) - (I_maxima_left_atrial_appendage_average+3*I_maxima_left_atrial_appendage_SD);
end
I_maxima_left_atrial_appendage_dev_avg_pls_3SD = abs(I_maxima_left_atrial_appendage_dev_avg_pls_3SD);
for i = 1:length(I_maxima_left_atrial_appendage_sorted)
    if I_maxima_left_atrial_appendage_dev_avg_pls_3SD(i) == min(I_maxima_left_atrial_appendage_dev_avg_pls_3SD)
        I_maxima_left_atrial_appendage_dev_avg_pls_3SD_index = i;
    end
end
I_maxima_left_atrial_appendage_tot_number_of_indices = length(I_maxima_left_atrial_appendage_sorted);
%This is the approximation of the tot "area under curve." 
I_maxima_left_atrial_appendage_quotient_within_1SD = (I_maxima_left_atrial_appendage_dev_avg_pls_1SD_index - I_maxima_left_atrial_appendagedev_avg_min_1SD_index)/I_maxima_left_atrial_appendage_tot_number_of_indices;
I_maxima_left_atrial_appendage_quotient_within_2SD = (I_maxima_left_atrial_appendage_dev_avg_pls_2SD_index - I_maxima_dev_avg_min_2SD_index)/I_maxima_left_atrial_appendage_tot_number_of_indices;
I_maxima_left_atrial_appendage_quotient_within_3SD = (I_maxima_left_atrial_appendage_dev_avg_pls_3SD_index - I_maxima_left_atrial_appendage_dev_avg_min_3SD_index)/I_maxima_left_atrial_appendage_tot_number_of_indices;
if (I_maxima_left_atrial_appendage_quotient_within_1SD > 0.66 && I_maxima_left_atrial_appendage_quotient_within_1SD < 0.70) && (I_maxima_left_atrial_appendage_quotient_within_2SD > 0.93 && I_maxima_left_atrial_appendage_quotient_within_2SD < 0.97) && (I_maxima_left_atrial_appendage_quotient_within_3SD > 0.98 && I_maxima_left_atrial_appendage_quotient_within_3SD < 1)
    isNormal(2,13) = 1;
end
if I_maxima_left_atrial_appendage_quotient_within_1SD == 0
    isSpike(2,13) = 1;
end
I_minima_SA_node_sorted = sort(I_minima_SA_node);
I_minima_SA_node_dev_avg_min_1SD = minima(length(I_minima_SA_node_sorted),1);
for i = 1:length(I_minima_SA_node_sorted)
    I_minima_SA_node_dev_avg_min_1SD(i) = I_minima_SA_node_sorted(i) - (I_minima_SA_node_average-I_minima_SA_node_SD);
end
I_minima_SA_node_dev_avg_min_1SD = abs(I_minima_SA_node_dev_avg_min_1SD);
for i = 1:length(I_minima_SA_node_sorted)
    if I_minima_SA_node_dev_avg_min_1SD(i) == min(I_minima_SA_node_dev_avg_min_1SD)
        I_minima_SA_node_dev_avg_min_1SD_index = i;
    end
end
I_minima_SA_node_dev_avg_pls_1SD = minima(length(I_minima_SA_node_sorted),1);
for i = 1:length(I_minima_SA_node_sorted)
    I_minima_SA_node_dev_avg_pls_1SD(i) = I_minima_SA_node_sorted(i) - (I_minima_SA_node_average+I_minima_SA_node_SD);
end
I_minima_SA_node_dev_avg_pls_1SD = abs(I_minima_SA_node_dev_avg_pls_1SD);
for i = 1:length(I_minima_SA_node_sorted)
    if I_minima_SA_node_dev_avg_pls_1SD(i) == min(I_minima_SA_node_dev_avg_pls_1SD)
        I_minima_SA_node_dev_avg_pls_1SD_index = i;
    end
end
I_minima_SA_node_dev_avg_min_2SD = minima(length(I_minima_SA_node_sorted),1);
for i = 1:length(I_minima_SA_node_sorted)
    I_minima_SA_node_dev_avg_min_2SD(i) = I_minima_SA_node_sorted(i) - (I_minima_SA_node_average-2*I_minima_SA_node_SD);
end
I_minima_SA_node_dev_avg_min_2SD = abs(I_minima_SA_node_dev_avg_min_2SD);
for i = 1:length(I_minima_SA_node_sorted)
    if I_minima_SA_node_dev_avg_min_2SD(i) == min(I_minima_SA_node_dev_avg_min_2SD)
        I_minima_dev_avg_min_2SD_index = i;
    end
end
I_minima_SA_node_dev_avg_pls_2SD = minima(length(I_minima_SA_node_sorted),1);
for i = 1:length(I_minima_SA_node_sorted)
    I_minima_SA_node_dev_avg_pls_2SD(i) = I_minima_SA_node_sorted(i) - (I_minima_SA_node_average+2*I_minima_SA_node_SD);
end
I_minima_SA_node_dev_avg_pls_2SD = abs(I_minima_SA_node_dev_avg_pls_2SD);
for i = 1:length(I_minima_SA_node_sorted)
    if I_minima_SA_node_dev_avg_pls_2SD(i) == min(I_minima_SA_node_dev_avg_pls_2SD)
        I_minima_SA_node_dev_avg_pls_2SD_index = i;
    end
end
I_minima_SA_node_dev_avg_min_3SD = minima(length(I_minima_SA_node_sorted),1);
for i = 1:length(I_minima_SA_node_sorted)
    I_minima_SA_node_dev_avg_min_3SD(i) = I_minima_SA_node_sorted(i) - (I_minima_SA_node_average-3*I_minima_SA_node_SD);
end
I_minima_SA_node_dev_avg_min_3SD = abs(I_minima_SA_node_dev_avg_min_3SD);
for i = 1:length(I_minima_SA_node_sorted)
    if I_minima_SA_node_dev_avg_min_3SD(i) == min(I_minima_SA_node_dev_avg_min_3SD)
        I_minima_SA_node_dev_avg_min_3SD_index = i;
    end
end
I_minima_SA_node_dev_avg_pls_3SD = minima(length(I_minima_SA_node_sorted),1);
for i = 1:length(I_minima_SA_node_sorted)
    I_minima_SA_node_dev_avg_pls_3SD(i) = I_minima_SA_node_sorted(i) - (I_minima_SA_node_average+3*I_minima_SA_node_SD);
end
I_minima_SA_node_dev_avg_pls_3SD = abs(I_minima_SA_node_dev_avg_pls_3SD);
for i = 1:length(I_minima_SA_node_sorted)
    if I_minima_SA_node_dev_avg_pls_3SD(i) == min(I_minima_SA_node_dev_avg_pls_3SD)
        I_minima_SA_node_dev_avg_pls_3SD_index = i;
    end
end
I_minima_SA_node_tot_number_of_indices = length(I_minima_SA_node_sorted);
%This is the approximation of the tot "area under curve." 
I_minima_SA_node_quotient_within_1SD = (I_minima_SA_node_dev_avg_pls_1SD_index - I_minima_SA_nodedev_avg_min_1SD_index)/I_minima_SA_node_tot_number_of_indices;
I_minima_SA_node_quotient_within_2SD = (I_minima_SA_node_dev_avg_pls_2SD_index - I_minima_dev_avg_min_2SD_index)/I_minima_SA_node_tot_number_of_indices;
I_minima_SA_node_quotient_within_3SD = (I_minima_SA_node_dev_avg_pls_3SD_index - I_minima_SA_node_dev_avg_min_3SD_index)/I_minima_SA_node_tot_number_of_indices;
if ((I_minima_SA_node_quotient_within_1SD > 0.66 && I_minima_SA_node_quotient_within_1SD < 0.70) && (I_minima_SA_node_quotient_within_2SD > 0.93 && I_minima_SA_node_quotient_within_2SD < 0.97) && (I_minima_SA_node_quotient_within_3SD > 0.98 && I_minima_SA_node_quotient_within_3SD < 1)) 
    isNormal(3,1) = 1;
end
if I_minima_SA_node_quotient_within_1SD == 0
    isSpike(3,1) = 1;
end
I_minima_crista_terminalis_sorted = sort(I_minima_crista_terminalis);
I_minima_crista_terminalis_dev_avg_min_1SD = minima(length(I_minima_crista_terminalis_sorted),1);
for i = 1:length(I_minima_crista_terminalis_sorted)
    I_minima_crista_terminalis_dev_avg_min_1SD(i) = I_minima_crista_terminalis_sorted(i) - (I_minima_crista_terminalis_average-I_minima_crista_terminalis_SD);
end
I_minima_crista_terminalis_dev_avg_min_1SD = abs(I_minima_crista_terminalis_dev_avg_min_1SD);
for i = 1:length(I_minima_crista_terminalis_sorted)
    if I_minima_crista_terminalis_dev_avg_min_1SD(i) == min(I_minima_crista_terminalis_dev_avg_min_1SD)
        I_minima_crista_terminalis_dev_avg_min_1SD_index = i;
    end
end
I_minima_crista_terminalis_dev_avg_pls_1SD = minima(length(I_minima_crista_terminalis_sorted),1);
for i = 1:length(I_minima_crista_terminalis_sorted)
    I_minima_crista_terminalis_dev_avg_pls_1SD(i) = I_minima_crista_terminalis_sorted(i) - (I_minima_crista_terminalis_average+I_minima_crista_terminalis_SD);
end
I_minima_crista_terminalis_dev_avg_pls_1SD = abs(I_minima_crista_terminalis_dev_avg_pls_1SD);
for i = 1:length(I_minima_crista_terminalis_sorted)
    if I_minima_crista_terminalis_dev_avg_pls_1SD(i) == min(I_minima_crista_terminalis_dev_avg_pls_1SD)
        I_minima_crista_terminalis_dev_avg_pls_1SD_index = i;
    end
end
I_minima_crista_terminalis_dev_avg_min_2SD = minima(length(I_minima_crista_terminalis_sorted),1);
for i = 1:length(I_minima_crista_terminalis_sorted)
    I_minima_crista_terminalis_dev_avg_min_2SD(i) = I_minima_crista_terminalis_sorted(i) - (I_minima_crista_terminalis_average-2*I_minima_crista_terminalis_SD);
end
I_minima_crista_terminalis_dev_avg_min_2SD = abs(I_minima_crista_terminalis_dev_avg_min_2SD);
for i = 1:length(I_minima_crista_terminalis_sorted)
    if I_minima_crista_terminalis_dev_avg_min_2SD(i) == min(I_minima_crista_terminalis_dev_avg_min_2SD)
        I_minima_dev_avg_min_2SD_index = i;
    end
end
I_minima_crista_terminalis_dev_avg_pls_2SD = minima(length(I_minima_crista_terminalis_sorted),1);
for i = 1:length(I_minima_crista_terminalis_sorted)
    I_minima_crista_terminalis_dev_avg_pls_2SD(i) = I_minima_crista_terminalis_sorted(i) - (I_minima_crista_terminalis_average+2*I_minima_crista_terminalis_SD);
end
I_minima_crista_terminalis_dev_avg_pls_2SD = abs(I_minima_crista_terminalis_dev_avg_pls_2SD);
for i = 1:length(I_minima_crista_terminalis_sorted)
    if I_minima_crista_terminalis_dev_avg_pls_2SD(i) == min(I_minima_crista_terminalis_dev_avg_pls_2SD)
        I_minima_crista_terminalis_dev_avg_pls_2SD_index = i;
    end
end
I_minima_crista_terminalis_dev_avg_min_3SD = minima(length(I_minima_crista_terminalis_sorted),1);
for i = 1:length(I_minima_crista_terminalis_sorted)
    I_minima_crista_terminalis_dev_avg_min_3SD(i) = I_minima_crista_terminalis_sorted(i) - (I_minima_crista_terminalis_average-3*I_minima_crista_terminalis_SD);
end
I_minima_crista_terminalis_dev_avg_min_3SD = abs(I_minima_crista_terminalis_dev_avg_min_3SD);
for i = 1:length(I_minima_crista_terminalis_sorted)
    if I_minima_crista_terminalis_dev_avg_min_3SD(i) == min(I_minima_crista_terminalis_dev_avg_min_3SD)
        I_minima_crista_terminalis_dev_avg_min_3SD_index = i;
    end
end
I_minima_crista_terminalis_dev_avg_pls_3SD = minima(length(I_minima_crista_terminalis_sorted),1);
for i = 1:length(I_minima_crista_terminalis_sorted)
    I_minima_crista_terminalis_dev_avg_pls_3SD(i) = I_minima_crista_terminalis_sorted(i) - (I_minima_crista_terminalis_average+3*I_minima_crista_terminalis_SD);
end
I_minima_crista_terminalis_dev_avg_pls_3SD = abs(I_minima_crista_terminalis_dev_avg_pls_3SD);
for i = 1:length(I_minima_crista_terminalis_sorted)
    if I_minima_crista_terminalis_dev_avg_pls_3SD(i) == min(I_minima_crista_terminalis_dev_avg_pls_3SD)
        I_minima_crista_terminalis_dev_avg_pls_3SD_index = i;
    end
end
I_minima_crista_terminalis_tot_number_of_indices = length(I_minima_crista_terminalis_sorted);
%This is the approximation of the tot "area under curve." 
I_minima_crista_terminalis_quotient_within_1SD = (I_minima_crista_terminalis_dev_avg_pls_1SD_index - I_minima_crista_terminalisdev_avg_min_1SD_index)/I_minima_crista_terminalis_tot_number_of_indices;
I_minima_crista_terminalis_quotient_within_2SD = (I_minima_crista_terminalis_dev_avg_pls_2SD_index - I_minima_dev_avg_min_2SD_index)/I_minima_crista_terminalis_tot_number_of_indices;
I_minima_crista_terminalis_quotient_within_3SD = (I_minima_crista_terminalis_dev_avg_pls_3SD_index - I_minima_crista_terminalis_dev_avg_min_3SD_index)/I_minima_crista_terminalis_tot_number_of_indices;
if (I_minima_crista_terminalis_quotient_within_1SD > 0.66 && I_minima_crista_terminalis_quotient_within_1SD < 0.70) && (I_minima_crista_terminalis_quotient_within_2SD > 0.93 && I_minima_crista_terminalis_quotient_within_2SD < 0.97) && (I_minima_crista_terminalis_quotient_within_3SD > 0.98 && I_minima_crista_terminalis_quotient_within_3SD < 1)
    isNormal(3,2) = 1;
end
if I_minima_crista_terminalis_quotient_within_1SD == 0
    isSpike(3,2) = 1;
end
I_minima_tricuspid_annulus_sorted = sort(I_minima_tricuspid_annulus);
I_minima_tricuspid_annulus_dev_avg_min_1SD = minima(length(I_minima_tricuspid_annulus_sorted),1);
for i = 1:length(I_minima_tricuspid_annulus_sorted)
    I_minima_tricuspid_annulus_dev_avg_min_1SD(i) = I_minima_tricuspid_annulus_sorted(i) - (I_minima_tricuspid_annulus_average-I_minima_tricuspid_annulus_SD);
end
I_minima_tricuspid_annulus_dev_avg_min_1SD = abs(I_minima_tricuspid_annulus_dev_avg_min_1SD);
for i = 1:length(I_minima_tricuspid_annulus_sorted)
    if I_minima_tricuspid_annulus_dev_avg_min_1SD(i) == min(I_minima_tricuspid_annulus_dev_avg_min_1SD)
        I_minima_tricuspid_annulus_dev_avg_min_1SD_index = i;
    end
end
I_minima_tricuspid_annulus_dev_avg_pls_1SD = minima(length(I_minima_tricuspid_annulus_sorted),1);
for i = 1:length(I_minima_tricuspid_annulus_sorted)
    I_minima_tricuspid_annulus_dev_avg_pls_1SD(i) = I_minima_tricuspid_annulus_sorted(i) - (I_minima_tricuspid_annulus_average+I_minima_tricuspid_annulus_SD);
end
I_minima_tricuspid_annulus_dev_avg_pls_1SD = abs(I_minima_tricuspid_annulus_dev_avg_pls_1SD);
for i = 1:length(I_minima_tricuspid_annulus_sorted)
    if I_minima_tricuspid_annulus_dev_avg_pls_1SD(i) == min(I_minima_tricuspid_annulus_dev_avg_pls_1SD)
        I_minima_tricuspid_annulus_dev_avg_pls_1SD_index = i;
    end
end
I_minima_tricuspid_annulus_dev_avg_min_2SD = minima(length(I_minima_tricuspid_annulus_sorted),1);
for i = 1:length(I_minima_tricuspid_annulus_sorted)
    I_minima_tricuspid_annulus_dev_avg_min_2SD(i) = I_minima_tricuspid_annulus_sorted(i) - (I_minima_tricuspid_annulus_average-2*I_minima_tricuspid_annulus_SD);
end
I_minima_tricuspid_annulus_dev_avg_min_2SD = abs(I_minima_tricuspid_annulus_dev_avg_min_2SD);
for i = 1:length(I_minima_tricuspid_annulus_sorted)
    if I_minima_tricuspid_annulus_dev_avg_min_2SD(i) == min(I_minima_tricuspid_annulus_dev_avg_min_2SD)
        I_minima_dev_avg_min_2SD_index = i;
    end
end
I_minima_tricuspid_annulus_dev_avg_pls_2SD = minima(length(I_minima_tricuspid_annulus_sorted),1);
for i = 1:length(I_minima_tricuspid_annulus_sorted)
    I_minima_tricuspid_annulus_dev_avg_pls_2SD(i) = I_minima_tricuspid_annulus_sorted(i) - (I_minima_tricuspid_annulus_average+2*I_minima_tricuspid_annulus_SD);
end
I_minima_tricuspid_annulus_dev_avg_pls_2SD = abs(I_minima_tricuspid_annulus_dev_avg_pls_2SD);
for i = 1:length(I_minima_tricuspid_annulus_sorted)
    if I_minima_tricuspid_annulus_dev_avg_pls_2SD(i) == min(I_minima_tricuspid_annulus_dev_avg_pls_2SD)
        I_minima_tricuspid_annulus_dev_avg_pls_2SD_index = i;
    end
end
I_minima_tricuspid_annulus_dev_avg_min_3SD = minima(length(I_minima_tricuspid_annulus_sorted),1);
for i = 1:length(I_minima_tricuspid_annulus_sorted)
    I_minima_tricuspid_annulus_dev_avg_min_3SD(i) = I_minima_tricuspid_annulus_sorted(i) - (I_minima_tricuspid_annulus_average-3*I_minima_tricuspid_annulus_SD);
end
I_minima_tricuspid_annulus_dev_avg_min_3SD = abs(I_minima_tricuspid_annulus_dev_avg_min_3SD);
for i = 1:length(I_minima_tricuspid_annulus_sorted)
    if I_minima_tricuspid_annulus_dev_avg_min_3SD(i) == min(I_minima_tricuspid_annulus_dev_avg_min_3SD)
        I_minima_tricuspid_annulus_dev_avg_min_3SD_index = i;
    end
end
I_minima_tricuspid_annulus_dev_avg_pls_3SD = minima(length(I_minima_tricuspid_annulus_sorted),1);
for i = 1:length(I_minima_tricuspid_annulus_sorted)
    I_minima_tricuspid_annulus_dev_avg_pls_3SD(i) = I_minima_tricuspid_annulus_sorted(i) - (I_minima_tricuspid_annulus_average+3*I_minima_tricuspid_annulus_SD);
end
I_minima_tricuspid_annulus_dev_avg_pls_3SD = abs(I_minima_tricuspid_annulus_dev_avg_pls_3SD);
for i = 1:length(I_minima_tricuspid_annulus_sorted)
    if I_minima_tricuspid_annulus_dev_avg_pls_3SD(i) == min(I_minima_tricuspid_annulus_dev_avg_pls_3SD)
        I_minima_tricuspid_annulus_dev_avg_pls_3SD_index = i;
    end
end
I_minima_tricuspid_annulus_tot_number_of_indices = length(I_minima_tricuspid_annulus_sorted);
%This is the approximation of the tot "area under curve." 
I_minima_tricuspid_annulus_quotient_within_1SD = (I_minima_tricuspid_annulus_dev_avg_pls_1SD_index - I_minima_tricuspid_annulusdev_avg_min_1SD_index)/I_minima_tricuspid_annulus_tot_number_of_indices;
I_minima_tricuspid_annulus_quotient_within_2SD = (I_minima_tricuspid_annulus_dev_avg_pls_2SD_index - I_minima_dev_avg_min_2SD_index)/I_minima_tricuspid_annulus_tot_number_of_indices;
I_minima_tricuspid_annulus_quotient_within_3SD = (I_minima_tricuspid_annulus_dev_avg_pls_3SD_index - I_minima_tricuspid_annulus_dev_avg_min_3SD_index)/I_minima_tricuspid_annulus_tot_number_of_indices;
if (I_minima_tricuspid_annulus_quotient_within_1SD > 0.66 && I_minima_tricuspid_annulus_quotient_within_1SD < 0.70) && (I_minima_tricuspid_annulus_quotient_within_2SD > 0.93 && I_minima_tricuspid_annulus_quotient_within_2SD < 0.97) && (I_minima_tricuspid_annulus_quotient_within_3SD > 0.98 && I_minima_tricuspid_annulus_quotient_within_3SD < 1)
    isNormal(3,3) = 1;
end
if I_minima_tricuspid_annulus_quotient_within_1SD == 0
    isSpike(3,3) = 1;
end
I_minima_coronary_sinus_sorted = sort(I_minima_coronary_sinus);
I_minima_coronary_sinus_dev_avg_min_1SD = minima(length(I_minima_coronary_sinus_sorted),1);
for i = 1:length(I_minima_coronary_sinus_sorted)
    I_minima_coronary_sinus_dev_avg_min_1SD(i) = I_minima_coronary_sinus_sorted(i) - (I_minima_coronary_sinus_average-I_minima_coronary_sinus_SD);
end
I_minima_coronary_sinus_dev_avg_min_1SD = abs(I_minima_coronary_sinus_dev_avg_min_1SD);
for i = 1:length(I_minima_coronary_sinus_sorted)
    if I_minima_coronary_sinus_dev_avg_min_1SD(i) == min(I_minima_coronary_sinus_dev_avg_min_1SD)
        I_minima_coronary_sinus_dev_avg_min_1SD_index = i;
    end
end
I_minima_coronary_sinus_dev_avg_pls_1SD = minima(length(I_minima_coronary_sinus_sorted),1);
for i = 1:length(I_minima_coronary_sinus_sorted)
    I_minima_coronary_sinus_dev_avg_pls_1SD(i) = I_minima_coronary_sinus_sorted(i) - (I_minima_coronary_sinus_average+I_minima_coronary_sinus_SD);
end
I_minima_coronary_sinus_dev_avg_pls_1SD = abs(I_minima_coronary_sinus_dev_avg_pls_1SD);
for i = 1:length(I_minima_coronary_sinus_sorted)
    if I_minima_coronary_sinus_dev_avg_pls_1SD(i) == min(I_minima_coronary_sinus_dev_avg_pls_1SD)
        I_minima_coronary_sinus_dev_avg_pls_1SD_index = i;
    end
end
I_minima_coronary_sinus_dev_avg_min_2SD = minima(length(I_minima_coronary_sinus_sorted),1);
for i = 1:length(I_minima_coronary_sinus_sorted)
    I_minima_coronary_sinus_dev_avg_min_2SD(i) = I_minima_coronary_sinus_sorted(i) - (I_minima_coronary_sinus_average-2*I_minima_coronary_sinus_SD);
end
I_minima_coronary_sinus_dev_avg_min_2SD = abs(I_minima_coronary_sinus_dev_avg_min_2SD);
for i = 1:length(I_minima_coronary_sinus_sorted)
    if I_minima_coronary_sinus_dev_avg_min_2SD(i) == min(I_minima_coronary_sinus_dev_avg_min_2SD)
        I_minima_dev_avg_min_2SD_index = i;
    end
end
I_minima_coronary_sinus_dev_avg_pls_2SD = minima(length(I_minima_coronary_sinus_sorted),1);
for i = 1:length(I_minima_coronary_sinus_sorted)
    I_minima_coronary_sinus_dev_avg_pls_2SD(i) = I_minima_coronary_sinus_sorted(i) - (I_minima_coronary_sinus_average+2*I_minima_coronary_sinus_SD);
end
I_minima_coronary_sinus_dev_avg_pls_2SD = abs(I_minima_coronary_sinus_dev_avg_pls_2SD);
for i = 1:length(I_minima_coronary_sinus_sorted)
    if I_minima_coronary_sinus_dev_avg_pls_2SD(i) == min(I_minima_coronary_sinus_dev_avg_pls_2SD)
        I_minima_coronary_sinus_dev_avg_pls_2SD_index = i;
    end
end
I_minima_coronary_sinus_dev_avg_min_3SD = minima(length(I_minima_coronary_sinus_sorted),1);
for i = 1:length(I_minima_coronary_sinus_sorted)
    I_minima_coronary_sinus_dev_avg_min_3SD(i) = I_minima_coronary_sinus_sorted(i) - (I_minima_coronary_sinus_average-3*I_minima_coronary_sinus_SD);
end
I_minima_coronary_sinus_dev_avg_min_3SD = abs(I_minima_coronary_sinus_dev_avg_min_3SD);
for i = 1:length(I_minima_coronary_sinus_sorted)
    if I_minima_coronary_sinus_dev_avg_min_3SD(i) == min(I_minima_coronary_sinus_dev_avg_min_3SD)
        I_minima_coronary_sinus_dev_avg_min_3SD_index = i;
    end
end
I_minima_coronary_sinus_dev_avg_pls_3SD = minima(length(I_minima_coronary_sinus_sorted),1);
for i = 1:length(I_minima_coronary_sinus_sorted)
    I_minima_coronary_sinus_dev_avg_pls_3SD(i) = I_minima_coronary_sinus_sorted(i) - (I_minima_coronary_sinus_average+3*I_minima_coronary_sinus_SD);
end
I_minima_coronary_sinus_dev_avg_pls_3SD = abs(I_minima_coronary_sinus_dev_avg_pls_3SD);
for i = 1:length(I_minima_coronary_sinus_sorted)
    if I_minima_coronary_sinus_dev_avg_pls_3SD(i) == min(I_minima_coronary_sinus_dev_avg_pls_3SD)
        I_minima_coronary_sinus_dev_avg_pls_3SD_index = i;
    end
end
I_minima_coronary_sinus_tot_number_of_indices = length(I_minima_coronary_sinus_sorted);
%This is the approximation of the tot "area under curve." 
I_minima_coronary_sinus_quotient_within_1SD = (I_minima_coronary_sinus_dev_avg_pls_1SD_index - I_minima_coronary_sinusdev_avg_min_1SD_index)/I_minima_coronary_sinus_tot_number_of_indices;
I_minima_coronary_sinus_quotient_within_2SD = (I_minima_coronary_sinus_dev_avg_pls_2SD_index - I_minima_dev_avg_min_2SD_index)/I_minima_coronary_sinus_tot_number_of_indices;
I_minima_coronary_sinus_quotient_within_3SD = (I_minima_coronary_sinus_dev_avg_pls_3SD_index - I_minima_coronary_sinus_dev_avg_min_3SD_index)/I_minima_coronary_sinus_tot_number_of_indices;
if (I_minima_coronary_sinus_quotient_within_1SD > 0.66 && I_minima_coronary_sinus_quotient_within_1SD < 0.70) && (I_minima_coronary_sinus_quotient_within_2SD > 0.93 && I_minima_coronary_sinus_quotient_within_2SD < 0.97) && (I_minima_coronary_sinus_quotient_within_3SD > 0.98 && I_minima_coronary_sinus_quotient_within_3SD < 1)
    isNormal(3,4) = 1;
end
if I_minima_coronary_sinus_quotient_within_1SD == 0
    isSpike(3,4) = 1;
end
I_minima_ostium_sorted = sort(I_minima_ostium);
I_minima_ostium_dev_avg_min_1SD = minima(length(I_minima_ostium_sorted),1);
for i = 1:length(I_minima_ostium_sorted)
    I_minima_ostium_dev_avg_min_1SD(i) = I_minima_ostium_sorted(i) - (I_minima_ostium_average-I_minima_ostium_SD);
end
I_minima_ostium_dev_avg_min_1SD = abs(I_minima_ostium_dev_avg_min_1SD);
for i = 1:length(I_minima_ostium_sorted)
    if I_minima_ostium_dev_avg_min_1SD(i) == min(I_minima_ostium_dev_avg_min_1SD)
        I_minima_ostium_dev_avg_min_1SD_index = i;
    end
end
I_minima_ostium_dev_avg_pls_1SD = minima(length(I_minima_ostium_sorted),1);
for i = 1:length(I_minima_ostium_sorted)
    I_minima_ostium_dev_avg_pls_1SD(i) = I_minima_ostium_sorted(i) - (I_minima_ostium_average+I_minima_ostium_SD);
end
I_minima_ostium_dev_avg_pls_1SD = abs(I_minima_ostium_dev_avg_pls_1SD);
for i = 1:length(I_minima_ostium_sorted)
    if I_minima_ostium_dev_avg_pls_1SD(i) == min(I_minima_ostium_dev_avg_pls_1SD)
        I_minima_ostium_dev_avg_pls_1SD_index = i;
    end
end
I_minima_ostium_dev_avg_min_2SD = minima(length(I_minima_ostium_sorted),1);
for i = 1:length(I_minima_ostium_sorted)
    I_minima_ostium_dev_avg_min_2SD(i) = I_minima_ostium_sorted(i) - (I_minima_ostium_average-2*I_minima_ostium_SD);
end
I_minima_ostium_dev_avg_min_2SD = abs(I_minima_ostium_dev_avg_min_2SD);
for i = 1:length(I_minima_ostium_sorted)
    if I_minima_ostium_dev_avg_min_2SD(i) == min(I_minima_ostium_dev_avg_min_2SD)
        I_minima_dev_avg_min_2SD_index = i;
    end
end
I_minima_ostium_dev_avg_pls_2SD = minima(length(I_minima_ostium_sorted),1);
for i = 1:length(I_minima_ostium_sorted)
    I_minima_ostium_dev_avg_pls_2SD(i) = I_minima_ostium_sorted(i) - (I_minima_ostium_average+2*I_minima_ostium_SD);
end
I_minima_ostium_dev_avg_pls_2SD = abs(I_minima_ostium_dev_avg_pls_2SD);
for i = 1:length(I_minima_ostium_sorted)
    if I_minima_ostium_dev_avg_pls_2SD(i) == min(I_minima_ostium_dev_avg_pls_2SD)
        I_minima_ostium_dev_avg_pls_2SD_index = i;
    end
end
I_minima_ostium_dev_avg_min_3SD = minima(length(I_minima_ostium_sorted),1);
for i = 1:length(I_minima_ostium_sorted)
    I_minima_ostium_dev_avg_min_3SD(i) = I_minima_ostium_sorted(i) - (I_minima_ostium_average-3*I_minima_ostium_SD);
end
I_minima_ostium_dev_avg_min_3SD = abs(I_minima_ostium_dev_avg_min_3SD);
for i = 1:length(I_minima_ostium_sorted)
    if I_minima_ostium_dev_avg_min_3SD(i) == min(I_minima_ostium_dev_avg_min_3SD)
        I_minima_ostium_dev_avg_min_3SD_index = i;
    end
end
I_minima_ostium_dev_avg_pls_3SD = minima(length(I_minima_ostium_sorted),1);
for i = 1:length(I_minima_ostium_sorted)
    I_minima_ostium_dev_avg_pls_3SD(i) = I_minima_ostium_sorted(i) - (I_minima_ostium_average+3*I_minima_ostium_SD);
end
I_minima_ostium_dev_avg_pls_3SD = abs(I_minima_ostium_dev_avg_pls_3SD);
for i = 1:length(I_minima_ostium_sorted)
    if I_minima_ostium_dev_avg_pls_3SD(i) == min(I_minima_ostium_dev_avg_pls_3SD)
        I_minima_ostium_dev_avg_pls_3SD_index = i;
    end
end
I_minima_ostium_tot_number_of_indices = length(I_minima_ostium_sorted);
%This is the approximation of the tot "area under curve." 
I_minima_ostium_quotient_within_1SD = (I_minima_ostium_dev_avg_pls_1SD_index - I_minima_ostiumdev_avg_min_1SD_index)/I_minima_ostium_tot_number_of_indices;
I_minima_ostium_quotient_within_2SD = (I_minima_ostium_dev_avg_pls_2SD_index - I_minima_dev_avg_min_2SD_index)/I_minima_ostium_tot_number_of_indices;
I_minima_ostium_quotient_within_3SD = (I_minima_ostium_dev_avg_pls_3SD_index - I_minima_ostium_dev_avg_min_3SD_index)/I_minima_ostium_tot_number_of_indices;
if (I_minima_ostium_quotient_within_1SD > 0.66 && I_minima_ostium_quotient_within_1SD < 0.70) && (I_minima_ostium_quotient_within_2SD > 0.93 && I_minima_ostium_quotient_within_2SD < 0.97) && (I_minima_ostium_quotient_within_3SD > 0.98 && I_minima_ostium_quotient_within_3SD < 1)
    isNormal(3,5) = 1;
end
if I_minima_ostium_quotient_within_1SD == 0
    isSpike(3,5) = 1;
end
I_minima_perinodal_sorted = sort(I_minima_perinodal);
I_minima_perinodal_dev_avg_min_1SD = minima(length(I_minima_perinodal_sorted),1);
for i = 1:length(I_minima_perinodal_sorted)
    I_minima_perinodal_dev_avg_min_1SD(i) = I_minima_perinodal_sorted(i) - (I_minima_perinodal_average-I_minima_perinodal_SD);
end
I_minima_perinodal_dev_avg_min_1SD = abs(I_minima_perinodal_dev_avg_min_1SD);
for i = 1:length(I_minima_perinodal_sorted)
    if I_minima_perinodal_dev_avg_min_1SD(i) == min(I_minima_perinodal_dev_avg_min_1SD)
        I_minima_perinodal_dev_avg_min_1SD_index = i;
    end
end
I_minima_perinodal_dev_avg_pls_1SD = minima(length(I_minima_perinodal_sorted),1);
for i = 1:length(I_minima_perinodal_sorted)
    I_minima_perinodal_dev_avg_pls_1SD(i) = I_minima_perinodal_sorted(i) - (I_minima_perinodal_average+I_minima_perinodal_SD);
end
I_minima_perinodal_dev_avg_pls_1SD = abs(I_minima_perinodal_dev_avg_pls_1SD);
for i = 1:length(I_minima_perinodal_sorted)
    if I_minima_perinodal_dev_avg_pls_1SD(i) == min(I_minima_perinodal_dev_avg_pls_1SD)
        I_minima_perinodal_dev_avg_pls_1SD_index = i;
    end
end
I_minima_perinodal_dev_avg_min_2SD = minima(length(I_minima_perinodal_sorted),1);
for i = 1:length(I_minima_perinodal_sorted)
    I_minima_perinodal_dev_avg_min_2SD(i) = I_minima_perinodal_sorted(i) - (I_minima_perinodal_average-2*I_minima_perinodal_SD);
end
I_minima_perinodal_dev_avg_min_2SD = abs(I_minima_perinodal_dev_avg_min_2SD);
for i = 1:length(I_minima_perinodal_sorted)
    if I_minima_perinodal_dev_avg_min_2SD(i) == min(I_minima_perinodal_dev_avg_min_2SD)
        I_minima_dev_avg_min_2SD_index = i;
    end
end
I_minima_perinodal_dev_avg_pls_2SD = minima(length(I_minima_perinodal_sorted),1);
for i = 1:length(I_minima_perinodal_sorted)
    I_minima_perinodal_dev_avg_pls_2SD(i) = I_minima_perinodal_sorted(i) - (I_minima_perinodal_average+2*I_minima_perinodal_SD);
end
I_minima_perinodal_dev_avg_pls_2SD = abs(I_minima_perinodal_dev_avg_pls_2SD);
for i = 1:length(I_minima_perinodal_sorted)
    if I_minima_perinodal_dev_avg_pls_2SD(i) == min(I_minima_perinodal_dev_avg_pls_2SD)
        I_minima_perinodal_dev_avg_pls_2SD_index = i;
    end
end
I_minima_perinodal_dev_avg_min_3SD = minima(length(I_minima_perinodal_sorted),1);
for i = 1:length(I_minima_perinodal_sorted)
    I_minima_perinodal_dev_avg_min_3SD(i) = I_minima_perinodal_sorted(i) - (I_minima_perinodal_average-3*I_minima_perinodal_SD);
end
I_minima_perinodal_dev_avg_min_3SD = abs(I_minima_perinodal_dev_avg_min_3SD);
for i = 1:length(I_minima_perinodal_sorted)
    if I_minima_perinodal_dev_avg_min_3SD(i) == min(I_minima_perinodal_dev_avg_min_3SD)
        I_minima_perinodal_dev_avg_min_3SD_index = i;
    end
end
I_minima_perinodal_dev_avg_pls_3SD = minima(length(I_minima_perinodal_sorted),1);
for i = 1:length(I_minima_perinodal_sorted)
    I_minima_perinodal_dev_avg_pls_3SD(i) = I_minima_perinodal_sorted(i) - (I_minima_perinodal_average+3*I_minima_perinodal_SD);
end
I_minima_perinodal_dev_avg_pls_3SD = abs(I_minima_perinodal_dev_avg_pls_3SD);
for i = 1:length(I_minima_perinodal_sorted)
    if I_minima_perinodal_dev_avg_pls_3SD(i) == min(I_minima_perinodal_dev_avg_pls_3SD)
        I_minima_perinodal_dev_avg_pls_3SD_index = i;
    end
end
I_minima_perinodal_tot_number_of_indices = length(I_minima_perinodal_sorted);
%This is the approximation of the tot "area under curve." 
I_minima_perinodal_quotient_within_1SD = (I_minima_perinodal_dev_avg_pls_1SD_index - I_minima_perinodaldev_avg_min_1SD_index)/I_minima_perinodal_tot_number_of_indices;
I_minima_perinodal_quotient_within_2SD = (I_minima_perinodal_dev_avg_pls_2SD_index - I_minima_dev_avg_min_2SD_index)/I_minima_perinodal_tot_number_of_indices;
I_minima_perinodal_quotient_within_3SD = (I_minima_perinodal_dev_avg_pls_3SD_index - I_minima_perinodal_dev_avg_min_3SD_index)/I_minima_perinodal_tot_number_of_indices;
if (I_minima_perinodal_quotient_within_1SD > 0.66 && I_minima_perinodal_quotient_within_1SD < 0.70) && (I_minima_perinodal_quotient_within_2SD > 0.93 && I_minima_perinodal_quotient_within_2SD < 0.97) && (I_minima_perinodal_quotient_within_3SD > 0.98 && I_minima_perinodal_quotient_within_3SD < 1)
    isNormal(3,6) = 1;
end
if I_minima_perinodal_quotient_within_1SD == 0
    isSpike(3,6) = 1;
end
I_minima_right_septum_sorted = sort(I_minima_right_septum);
I_minima_right_septum_dev_avg_min_1SD = minima(length(I_minima_right_septum_sorted),1);
for i = 1:length(I_minima_right_septum_sorted)
    I_minima_right_septum_dev_avg_min_1SD(i) = I_minima_right_septum_sorted(i) - (I_minima_right_septum_average-I_minima_right_septum_SD);
end
I_minima_right_septum_dev_avg_min_1SD = abs(I_minima_right_septum_dev_avg_min_1SD);
for i = 1:length(I_minima_right_septum_sorted)
    if I_minima_right_septum_dev_avg_min_1SD(i) == min(I_minima_right_septum_dev_avg_min_1SD)
        I_minima_right_septum_dev_avg_min_1SD_index = i;
    end
end
I_minima_right_septum_dev_avg_pls_1SD = minima(length(I_minima_right_septum_sorted),1);
for i = 1:length(I_minima_right_septum_sorted)
    I_minima_right_septum_dev_avg_pls_1SD(i) = I_minima_right_septum_sorted(i) - (I_minima_right_septum_average+I_minima_right_septum_SD);
end
I_minima_right_septum_dev_avg_pls_1SD = abs(I_minima_right_septum_dev_avg_pls_1SD);
for i = 1:length(I_minima_right_septum_sorted)
    if I_minima_right_septum_dev_avg_pls_1SD(i) == min(I_minima_right_septum_dev_avg_pls_1SD)
        I_minima_right_septum_dev_avg_pls_1SD_index = i;
    end
end
I_minima_right_septum_dev_avg_min_2SD = minima(length(I_minima_right_septum_sorted),1);
for i = 1:length(I_minima_right_septum_sorted)
    I_minima_right_septum_dev_avg_min_2SD(i) = I_minima_right_septum_sorted(i) - (I_minima_right_septum_average-2*I_minima_right_septum_SD);
end
I_minima_right_septum_dev_avg_min_2SD = abs(I_minima_right_septum_dev_avg_min_2SD);
for i = 1:length(I_minima_right_septum_sorted)
    if I_minima_right_septum_dev_avg_min_2SD(i) == min(I_minima_right_septum_dev_avg_min_2SD)
        I_minima_dev_avg_min_2SD_index = i;
    end
end
I_minima_right_septum_dev_avg_pls_2SD = minima(length(I_minima_right_septum_sorted),1);
for i = 1:length(I_minima_right_septum_sorted)
    I_minima_right_septum_dev_avg_pls_2SD(i) = I_minima_right_septum_sorted(i) - (I_minima_right_septum_average+2*I_minima_right_septum_SD);
end
I_minima_right_septum_dev_avg_pls_2SD = abs(I_minima_right_septum_dev_avg_pls_2SD);
for i = 1:length(I_minima_right_septum_sorted)
    if I_minima_right_septum_dev_avg_pls_2SD(i) == min(I_minima_right_septum_dev_avg_pls_2SD)
        I_minima_right_septum_dev_avg_pls_2SD_index = i;
    end
end
I_minima_right_septum_dev_avg_min_3SD = minima(length(I_minima_right_septum_sorted),1);
for i = 1:length(I_minima_right_septum_sorted)
    I_minima_right_septum_dev_avg_min_3SD(i) = I_minima_right_septum_sorted(i) - (I_minima_right_septum_average-3*I_minima_right_septum_SD);
end
I_minima_right_septum_dev_avg_min_3SD = abs(I_minima_right_septum_dev_avg_min_3SD);
for i = 1:length(I_minima_right_septum_sorted)
    if I_minima_right_septum_dev_avg_min_3SD(i) == min(I_minima_right_septum_dev_avg_min_3SD)
        I_minima_right_septum_dev_avg_min_3SD_index = i;
    end
end
I_minima_right_septum_dev_avg_pls_3SD = minima(length(I_minima_right_septum_sorted),1);
for i = 1:length(I_minima_right_septum_sorted)
    I_minima_right_septum_dev_avg_pls_3SD(i) = I_minima_right_septum_sorted(i) - (I_minima_right_septum_average+3*I_minima_right_septum_SD);
end
I_minima_right_septum_dev_avg_pls_3SD = abs(I_minima_right_septum_dev_avg_pls_3SD);
for i = 1:length(I_minima_right_septum_sorted)
    if I_minima_right_septum_dev_avg_pls_3SD(i) == min(I_minima_right_septum_dev_avg_pls_3SD)
        I_minima_right_septum_dev_avg_pls_3SD_index = i;
    end
end
I_minima_right_septum_tot_number_of_indices = length(I_minima_right_septum_sorted);
%This is the approximation of the tot "area under curve." 
I_minima_right_septum_quotient_within_1SD = (I_minima_right_septum_dev_avg_pls_1SD_index - I_minima_right_septumdev_avg_min_1SD_index)/I_minima_right_septum_tot_number_of_indices;
I_minima_right_septum_quotient_within_2SD = (I_minima_right_septum_dev_avg_pls_2SD_index - I_minima_dev_avg_min_2SD_index)/I_minima_right_septum_tot_number_of_indices;
I_minima_right_septum_quotient_within_3SD = (I_minima_right_septum_dev_avg_pls_3SD_index - I_minima_right_septum_dev_avg_min_3SD_index)/I_minima_right_septum_tot_number_of_indices;
if (I_minima_right_septum_quotient_within_1SD > 0.66 && I_minima_right_septum_quotient_within_1SD < 0.70) && (I_minima_right_septum_quotient_within_2SD > 0.93 && I_minima_right_septum_quotient_within_2SD < 0.97) && (I_minima_right_septum_quotient_within_3SD > 0.98 && I_minima_right_septum_quotient_within_3SD < 1)
    isNormal(3,7) = 1;
end
if I_minima_right_septum_quotient_within_1SD == 0
    isSpike(3,7) = 1;
end
I_minima_right_atrial_appendage_sorted = sort(I_minima_right_atrial_appendage);
I_minima_right_atrial_appendage_dev_avg_min_1SD = minima(length(I_minima_right_atrial_appendage_sorted),1);
for i = 1:length(I_minima_right_atrial_appendage_sorted)
    I_minima_right_atrial_appendage_dev_avg_min_1SD(i) = I_minima_right_atrial_appendage_sorted(i) - (I_minima_right_atrial_appendage_average-I_minima_right_atrial_appendage_SD);
end
I_minima_right_atrial_appendage_dev_avg_min_1SD = abs(I_minima_right_atrial_appendage_dev_avg_min_1SD);
for i = 1:length(I_minima_right_atrial_appendage_sorted)
    if I_minima_right_atrial_appendage_dev_avg_min_1SD(i) == min(I_minima_right_atrial_appendage_dev_avg_min_1SD)
        I_minima_right_atrial_appendage_dev_avg_min_1SD_index = i;
    end
end
I_minima_right_atrial_appendage_dev_avg_pls_1SD = minima(length(I_minima_right_atrial_appendage_sorted),1);
for i = 1:length(I_minima_right_atrial_appendage_sorted)
    I_minima_right_atrial_appendage_dev_avg_pls_1SD(i) = I_minima_right_atrial_appendage_sorted(i) - (I_minima_right_atrial_appendage_average+I_minima_right_atrial_appendage_SD);
end
I_minima_right_atrial_appendage_dev_avg_pls_1SD = abs(I_minima_right_atrial_appendage_dev_avg_pls_1SD);
for i = 1:length(I_minima_right_atrial_appendage_sorted)
    if I_minima_right_atrial_appendage_dev_avg_pls_1SD(i) == min(I_minima_right_atrial_appendage_dev_avg_pls_1SD)
        I_minima_right_atrial_appendage_dev_avg_pls_1SD_index = i;
    end
end
I_minima_right_atrial_appendage_dev_avg_min_2SD = minima(length(I_minima_right_atrial_appendage_sorted),1);
for i = 1:length(I_minima_right_atrial_appendage_sorted)
    I_minima_right_atrial_appendage_dev_avg_min_2SD(i) = I_minima_right_atrial_appendage_sorted(i) - (I_minima_right_atrial_appendage_average-2*I_minima_right_atrial_appendage_SD);
end
I_minima_right_atrial_appendage_dev_avg_min_2SD = abs(I_minima_right_atrial_appendage_dev_avg_min_2SD);
for i = 1:length(I_minima_right_atrial_appendage_sorted)
    if I_minima_right_atrial_appendage_dev_avg_min_2SD(i) == min(I_minima_right_atrial_appendage_dev_avg_min_2SD)
        I_minima_dev_avg_min_2SD_index = i;
    end
end
I_minima_right_atrial_appendage_dev_avg_pls_2SD = minima(length(I_minima_right_atrial_appendage_sorted),1);
for i = 1:length(I_minima_right_atrial_appendage_sorted)
    I_minima_right_atrial_appendage_dev_avg_pls_2SD(i) = I_minima_right_atrial_appendage_sorted(i) - (I_minima_right_atrial_appendage_average+2*I_minima_right_atrial_appendage_SD);
end
I_minima_right_atrial_appendage_dev_avg_pls_2SD = abs(I_minima_right_atrial_appendage_dev_avg_pls_2SD);
for i = 1:length(I_minima_right_atrial_appendage_sorted)
    if I_minima_right_atrial_appendage_dev_avg_pls_2SD(i) == min(I_minima_right_atrial_appendage_dev_avg_pls_2SD)
        I_minima_right_atrial_appendage_dev_avg_pls_2SD_index = i;
    end
end
I_minima_right_atrial_appendage_dev_avg_min_3SD = minima(length(I_minima_right_atrial_appendage_sorted),1);
for i = 1:length(I_minima_right_atrial_appendage_sorted)
    I_minima_right_atrial_appendage_dev_avg_min_3SD(i) = I_minima_right_atrial_appendage_sorted(i) - (I_minima_right_atrial_appendage_average-3*I_minima_right_atrial_appendage_SD);
end
I_minima_right_atrial_appendage_dev_avg_min_3SD = abs(I_minima_right_atrial_appendage_dev_avg_min_3SD);
for i = 1:length(I_minima_right_atrial_appendage_sorted)
    if I_minima_right_atrial_appendage_dev_avg_min_3SD(i) == min(I_minima_right_atrial_appendage_dev_avg_min_3SD)
        I_minima_right_atrial_appendage_dev_avg_min_3SD_index = i;
    end
end
I_minima_right_atrial_appendage_dev_avg_pls_3SD = minima(length(I_minima_right_atrial_appendage_sorted),1);
for i = 1:length(I_minima_right_atrial_appendage_sorted)
    I_minima_right_atrial_appendage_dev_avg_pls_3SD(i) = I_minima_right_atrial_appendage_sorted(i) - (I_minima_right_atrial_appendage_average+3*I_minima_right_atrial_appendage_SD);
end
I_minima_right_atrial_appendage_dev_avg_pls_3SD = abs(I_minima_right_atrial_appendage_dev_avg_pls_3SD);
for i = 1:length(I_minima_right_atrial_appendage_sorted)
    if I_minima_right_atrial_appendage_dev_avg_pls_3SD(i) == min(I_minima_right_atrial_appendage_dev_avg_pls_3SD)
        I_minima_right_atrial_appendage_dev_avg_pls_3SD_index = i;
    end
end
I_minima_right_atrial_appendage_tot_number_of_indices = length(I_minima_right_atrial_appendage_sorted);
%This is the approximation of the tot "area under curve." 
I_minima_right_atrial_appendage_quotient_within_1SD = (I_minima_right_atrial_appendage_dev_avg_pls_1SD_index - I_minima_right_atrial_appendagedev_avg_min_1SD_index)/I_minima_right_atrial_appendage_tot_number_of_indices;
I_minima_right_atrial_appendage_quotient_within_2SD = (I_minima_right_atrial_appendage_dev_avg_pls_2SD_index - I_minima_dev_avg_min_2SD_index)/I_minima_right_atrial_appendage_tot_number_of_indices;
I_minima_right_atrial_appendage_quotient_within_3SD = (I_minima_right_atrial_appendage_dev_avg_pls_3SD_index - I_minima_right_atrial_appendage_dev_avg_min_3SD_index)/I_minima_right_atrial_appendage_tot_number_of_indices;
if (I_minima_right_atrial_appendage_quotient_within_1SD > 0.66 && I_minima_right_atrial_appendage_quotient_within_1SD < 0.70) && (I_minima_right_atrial_appendage_quotient_within_2SD > 0.93 && I_minima_right_atrial_appendage_quotient_within_2SD < 0.97) && (I_minima_right_atrial_appendage_quotient_within_3SD > 0.98 && I_minima_right_atrial_appendage_quotient_within_3SD < 1)
    isNormal(3,8) = 1;
end
if I_minima_right_atrial_appendage_quotient_within_1SD == 0
    isSpike(3,8) = 1;
end
I_minima_pulmonary_veins_sorted = sort(I_minima_pulmonary_veins);
I_minima_pulmonary_veins_dev_avg_min_1SD = minima(length(I_minima_pulmonary_veins_sorted),1);
for i = 1:length(I_minima_pulmonary_veins_sorted)
    I_minima_pulmonary_veins_dev_avg_min_1SD(i) = I_minima_pulmonary_veins_sorted(i) - (I_minima_pulmonary_veins_average-I_minima_pulmonary_veins_SD);
end
I_minima_pulmonary_veins_dev_avg_min_1SD = abs(I_minima_pulmonary_veins_dev_avg_min_1SD);
for i = 1:length(I_minima_pulmonary_veins_sorted)
    if I_minima_pulmonary_veins_dev_avg_min_1SD(i) == min(I_minima_pulmonary_veins_dev_avg_min_1SD)
        I_minima_pulmonary_veins_dev_avg_min_1SD_index = i;
    end
end
I_minima_pulmonary_veins_dev_avg_pls_1SD = minima(length(I_minima_pulmonary_veins_sorted),1);
for i = 1:length(I_minima_pulmonary_veins_sorted)
    I_minima_pulmonary_veins_dev_avg_pls_1SD(i) = I_minima_pulmonary_veins_sorted(i) - (I_minima_pulmonary_veins_average+I_minima_pulmonary_veins_SD);
end
I_minima_pulmonary_veins_dev_avg_pls_1SD = abs(I_minima_pulmonary_veins_dev_avg_pls_1SD);
for i = 1:length(I_minima_pulmonary_veins_sorted)
    if I_minima_pulmonary_veins_dev_avg_pls_1SD(i) == min(I_minima_pulmonary_veins_dev_avg_pls_1SD)
        I_minima_pulmonary_veins_dev_avg_pls_1SD_index = i;
    end
end
I_minima_pulmonary_veins_dev_avg_min_2SD = minima(length(I_minima_pulmonary_veins_sorted),1);
for i = 1:length(I_minima_pulmonary_veins_sorted)
    I_minima_pulmonary_veins_dev_avg_min_2SD(i) = I_minima_pulmonary_veins_sorted(i) - (I_minima_pulmonary_veins_average-2*I_minima_pulmonary_veins_SD);
end
I_minima_pulmonary_veins_dev_avg_min_2SD = abs(I_minima_pulmonary_veins_dev_avg_min_2SD);
for i = 1:length(I_minima_pulmonary_veins_sorted)
    if I_minima_pulmonary_veins_dev_avg_min_2SD(i) == min(I_minima_pulmonary_veins_dev_avg_min_2SD)
        I_minima_dev_avg_min_2SD_index = i;
    end
end
I_minima_pulmonary_veins_dev_avg_pls_2SD = minima(length(I_minima_pulmonary_veins_sorted),1);
for i = 1:length(I_minima_pulmonary_veins_sorted)
    I_minima_pulmonary_veins_dev_avg_pls_2SD(i) = I_minima_pulmonary_veins_sorted(i) - (I_minima_pulmonary_veins_average+2*I_minima_pulmonary_veins_SD);
end
I_minima_pulmonary_veins_dev_avg_pls_2SD = abs(I_minima_pulmonary_veins_dev_avg_pls_2SD);
for i = 1:length(I_minima_pulmonary_veins_sorted)
    if I_minima_pulmonary_veins_dev_avg_pls_2SD(i) == min(I_minima_pulmonary_veins_dev_avg_pls_2SD)
        I_minima_pulmonary_veins_dev_avg_pls_2SD_index = i;
    end
end
I_minima_pulmonary_veins_dev_avg_min_3SD = minima(length(I_minima_pulmonary_veins_sorted),1);
for i = 1:length(I_minima_pulmonary_veins_sorted)
    I_minima_pulmonary_veins_dev_avg_min_3SD(i) = I_minima_pulmonary_veins_sorted(i) - (I_minima_pulmonary_veins_average-3*I_minima_pulmonary_veins_SD);
end
I_minima_pulmonary_veins_dev_avg_min_3SD = abs(I_minima_pulmonary_veins_dev_avg_min_3SD);
for i = 1:length(I_minima_pulmonary_veins_sorted)
    if I_minima_pulmonary_veins_dev_avg_min_3SD(i) == min(I_minima_pulmonary_veins_dev_avg_min_3SD)
        I_minima_pulmonary_veins_dev_avg_min_3SD_index = i;
    end
end
I_minima_pulmonary_veins_dev_avg_pls_3SD = minima(length(I_minima_pulmonary_veins_sorted),1);
for i = 1:length(I_minima_pulmonary_veins_sorted)
    I_minima_pulmonary_veins_dev_avg_pls_3SD(i) = I_minima_pulmonary_veins_sorted(i) - (I_minima_pulmonary_veins_average+3*I_minima_pulmonary_veins_SD);
end
I_minima_pulmonary_veins_dev_avg_pls_3SD = abs(I_minima_pulmonary_veins_dev_avg_pls_3SD);
for i = 1:length(I_minima_pulmonary_veins_sorted)
    if I_minima_pulmonary_veins_dev_avg_pls_3SD(i) == min(I_minima_pulmonary_veins_dev_avg_pls_3SD)
        I_minima_pulmonary_veins_dev_avg_pls_3SD_index = i;
    end
end
I_minima_pulmonary_veins_tot_number_of_indices = length(I_minima_pulmonary_veins_sorted);
%This is the approximation of the tot "area under curve." 
I_minima_pulmonary_veins_quotient_within_1SD = (I_minima_pulmonary_veins_dev_avg_pls_1SD_index - I_minima_pulmonary_veinsdev_avg_min_1SD_index)/I_minima_pulmonary_veins_tot_number_of_indices;
I_minima_pulmonary_veins_quotient_within_2SD = (I_minima_pulmonary_veins_dev_avg_pls_2SD_index - I_minima_dev_avg_min_2SD_index)/I_minima_pulmonary_veins_tot_number_of_indices;
I_minima_pulmonary_veins_quotient_within_3SD = (I_minima_pulmonary_veins_dev_avg_pls_3SD_index - I_minima_pulmonary_veins_dev_avg_min_3SD_index)/I_minima_pulmonary_veins_tot_number_of_indices;
if (I_minima_pulmonary_veins_quotient_within_1SD > 0.66 && I_minima_pulmonary_veins_quotient_within_1SD < 0.70) && (I_minima_pulmonary_veins_quotient_within_2SD > 0.93 && I_minima_pulmonary_veins_quotient_within_2SD < 0.97) && (I_minima_pulmonary_veins_quotient_within_3SD > 0.98 && I_minima_pulmonary_veins_quotient_within_3SD < 1)
    isNormal(3,9) = 1;
end
if I_minima_pulmonary_veins_quotient_within_1SD == 0
    isSpike(3,9) = 1;
end
I_minima_mitral_annulus_sorted = sort(I_minima_mitral_annulus);
I_minima_mitral_annulus_dev_avg_min_1SD = minima(length(I_minima_mitral_annulus_sorted),1);
for i = 1:length(I_minima_mitral_annulus_sorted)
    I_minima_mitral_annulus_dev_avg_min_1SD(i) = I_minima_mitral_annulus_sorted(i) - (I_minima_mitral_annulus_average-I_minima_mitral_annulus_SD);
end
I_minima_mitral_annulus_dev_avg_min_1SD = abs(I_minima_mitral_annulus_dev_avg_min_1SD);
for i = 1:length(I_minima_mitral_annulus_sorted)
    if I_minima_mitral_annulus_dev_avg_min_1SD(i) == min(I_minima_mitral_annulus_dev_avg_min_1SD)
        I_minima_mitral_annulus_dev_avg_min_1SD_index = i;
    end
end
I_minima_mitral_annulus_dev_avg_pls_1SD = minima(length(I_minima_mitral_annulus_sorted),1);
for i = 1:length(I_minima_mitral_annulus_sorted)
    I_minima_mitral_annulus_dev_avg_pls_1SD(i) = I_minima_mitral_annulus_sorted(i) - (I_minima_mitral_annulus_average+I_minima_mitral_annulus_SD);
end
I_minima_mitral_annulus_dev_avg_pls_1SD = abs(I_minima_mitral_annulus_dev_avg_pls_1SD);
for i = 1:length(I_minima_mitral_annulus_sorted)
    if I_minima_mitral_annulus_dev_avg_pls_1SD(i) == min(I_minima_mitral_annulus_dev_avg_pls_1SD)
        I_minima_mitral_annulus_dev_avg_pls_1SD_index = i;
    end
end
I_minima_mitral_annulus_dev_avg_min_2SD = minima(length(I_minima_mitral_annulus_sorted),1);
for i = 1:length(I_minima_mitral_annulus_sorted)
    I_minima_mitral_annulus_dev_avg_min_2SD(i) = I_minima_mitral_annulus_sorted(i) - (I_minima_mitral_annulus_average-2*I_minima_mitral_annulus_SD);
end
I_minima_mitral_annulus_dev_avg_min_2SD = abs(I_minima_mitral_annulus_dev_avg_min_2SD);
for i = 1:length(I_minima_mitral_annulus_sorted)
    if I_minima_mitral_annulus_dev_avg_min_2SD(i) == min(I_minima_mitral_annulus_dev_avg_min_2SD)
        I_minima_dev_avg_min_2SD_index = i;
    end
end
I_minima_mitral_annulus_dev_avg_pls_2SD = minima(length(I_minima_mitral_annulus_sorted),1);
for i = 1:length(I_minima_mitral_annulus_sorted)
    I_minima_mitral_annulus_dev_avg_pls_2SD(i) = I_minima_mitral_annulus_sorted(i) - (I_minima_mitral_annulus_average+2*I_minima_mitral_annulus_SD);
end
I_minima_mitral_annulus_dev_avg_pls_2SD = abs(I_minima_mitral_annulus_dev_avg_pls_2SD);
for i = 1:length(I_minima_mitral_annulus_sorted)
    if I_minima_mitral_annulus_dev_avg_pls_2SD(i) == min(I_minima_mitral_annulus_dev_avg_pls_2SD)
        I_minima_mitral_annulus_dev_avg_pls_2SD_index = i;
    end
end
I_minima_mitral_annulus_dev_avg_min_3SD = minima(length(I_minima_mitral_annulus_sorted),1);
for i = 1:length(I_minima_mitral_annulus_sorted)
    I_minima_mitral_annulus_dev_avg_min_3SD(i) = I_minima_mitral_annulus_sorted(i) - (I_minima_mitral_annulus_average-3*I_minima_mitral_annulus_SD);
end
I_minima_mitral_annulus_dev_avg_min_3SD = abs(I_minima_mitral_annulus_dev_avg_min_3SD);
for i = 1:length(I_minima_mitral_annulus_sorted)
    if I_minima_mitral_annulus_dev_avg_min_3SD(i) == min(I_minima_mitral_annulus_dev_avg_min_3SD)
        I_minima_mitral_annulus_dev_avg_min_3SD_index = i;
    end
end
I_minima_mitral_annulus_dev_avg_pls_3SD = minima(length(I_minima_mitral_annulus_sorted),1);
for i = 1:length(I_minima_mitral_annulus_sorted)
    I_minima_mitral_annulus_dev_avg_pls_3SD(i) = I_minima_mitral_annulus_sorted(i) - (I_minima_mitral_annulus_average+3*I_minima_mitral_annulus_SD);
end
I_minima_mitral_annulus_dev_avg_pls_3SD = abs(I_minima_mitral_annulus_dev_avg_pls_3SD);
for i = 1:length(I_minima_mitral_annulus_sorted)
    if I_minima_mitral_annulus_dev_avg_pls_3SD(i) == min(I_minima_mitral_annulus_dev_avg_pls_3SD)
        I_minima_mitral_annulus_dev_avg_pls_3SD_index = i;
    end
end
I_minima_mitral_annulus_tot_number_of_indices = length(I_minima_mitral_annulus_sorted);
%This is the approximation of the tot "area under curve." 
I_minima_mitral_annulus_quotient_within_1SD = (I_minima_mitral_annulus_dev_avg_pls_1SD_index - I_minima_mitral_annulusdev_avg_min_1SD_index)/I_minima_mitral_annulus_tot_number_of_indices;
I_minima_mitral_annulus_quotient_within_2SD = (I_minima_mitral_annulus_dev_avg_pls_2SD_index - I_minima_dev_avg_min_2SD_index)/I_minima_mitral_annulus_tot_number_of_indices;
I_minima_mitral_annulus_quotient_within_3SD = (I_minima_mitral_annulus_dev_avg_pls_3SD_index - I_minima_mitral_annulus_dev_avg_min_3SD_index)/I_minima_mitral_annulus_tot_number_of_indices;
if (I_minima_mitral_annulus_quotient_within_1SD > 0.66 && I_minima_mitral_annulus_quotient_within_1SD < 0.70) && (I_minima_mitral_annulus_quotient_within_2SD > 0.93 && I_minima_mitral_annulus_quotient_within_2SD < 0.97) && (I_minima_mitral_annulus_quotient_within_3SD > 0.98 && I_minima_mitral_annulus_quotient_within_3SD < 1)
    isNormal(3,10) = 1;
end
if I_minima_mitral_annulus_quotient_within_1SD == 0
    isSpike(3,10) = 1;
end
I_minima_CS_body_sorted = sort(I_minima_CS_body);
I_minima_CS_body_dev_avg_min_1SD = minima(length(I_minima_CS_body_sorted),1);
for i = 1:length(I_minima_CS_body_sorted)
    I_minima_CS_body_dev_avg_min_1SD(i) = I_minima_CS_body_sorted(i) - (I_minima_CS_body_average-I_minima_CS_body_SD);
end
I_minima_CS_body_dev_avg_min_1SD = abs(I_minima_CS_body_dev_avg_min_1SD);
for i = 1:length(I_minima_CS_body_sorted)
    if I_minima_CS_body_dev_avg_min_1SD(i) == min(I_minima_CS_body_dev_avg_min_1SD)
        I_minima_CS_body_dev_avg_min_1SD_index = i;
    end
end
I_minima_CS_body_dev_avg_pls_1SD = minima(length(I_minima_CS_body_sorted),1);
for i = 1:length(I_minima_CS_body_sorted)
    I_minima_CS_body_dev_avg_pls_1SD(i) = I_minima_CS_body_sorted(i) - (I_minima_CS_body_average+I_minima_CS_body_SD);
end
I_minima_CS_body_dev_avg_pls_1SD = abs(I_minima_CS_body_dev_avg_pls_1SD);
for i = 1:length(I_minima_CS_body_sorted)
    if I_minima_CS_body_dev_avg_pls_1SD(i) == min(I_minima_CS_body_dev_avg_pls_1SD)
        I_minima_CS_body_dev_avg_pls_1SD_index = i;
    end
end
I_minima_CS_body_dev_avg_min_2SD = minima(length(I_minima_CS_body_sorted),1);
for i = 1:length(I_minima_CS_body_sorted)
    I_minima_CS_body_dev_avg_min_2SD(i) = I_minima_CS_body_sorted(i) - (I_minima_CS_body_average-2*I_minima_CS_body_SD);
end
I_minima_CS_body_dev_avg_min_2SD = abs(I_minima_CS_body_dev_avg_min_2SD);
for i = 1:length(I_minima_CS_body_sorted)
    if I_minima_CS_body_dev_avg_min_2SD(i) == min(I_minima_CS_body_dev_avg_min_2SD)
        I_minima_dev_avg_min_2SD_index = i;
    end
end
I_minima_CS_body_dev_avg_pls_2SD = minima(length(I_minima_CS_body_sorted),1);
for i = 1:length(I_minima_CS_body_sorted)
    I_minima_CS_body_dev_avg_pls_2SD(i) = I_minima_CS_body_sorted(i) - (I_minima_CS_body_average+2*I_minima_CS_body_SD);
end
I_minima_CS_body_dev_avg_pls_2SD = abs(I_minima_CS_body_dev_avg_pls_2SD);
for i = 1:length(I_minima_CS_body_sorted)
    if I_minima_CS_body_dev_avg_pls_2SD(i) == min(I_minima_CS_body_dev_avg_pls_2SD)
        I_minima_CS_body_dev_avg_pls_2SD_index = i;
    end
end
I_minima_CS_body_dev_avg_min_3SD = minima(length(I_minima_CS_body_sorted),1);
for i = 1:length(I_minima_CS_body_sorted)
    I_minima_CS_body_dev_avg_min_3SD(i) = I_minima_CS_body_sorted(i) - (I_minima_CS_body_average-3*I_minima_CS_body_SD);
end
I_minima_CS_body_dev_avg_min_3SD = abs(I_minima_CS_body_dev_avg_min_3SD);
for i = 1:length(I_minima_CS_body_sorted)
    if I_minima_CS_body_dev_avg_min_3SD(i) == min(I_minima_CS_body_dev_avg_min_3SD)
        I_minima_CS_body_dev_avg_min_3SD_index = i;
    end
end
I_minima_CS_body_dev_avg_pls_3SD = minima(length(I_minima_CS_body_sorted),1);
for i = 1:length(I_minima_CS_body_sorted)
    I_minima_CS_body_dev_avg_pls_3SD(i) = I_minima_CS_body_sorted(i) - (I_minima_CS_body_average+3*I_minima_CS_body_SD);
end
I_minima_CS_body_dev_avg_pls_3SD = abs(I_minima_CS_body_dev_avg_pls_3SD);
for i = 1:length(I_minima_CS_body_sorted)
    if I_minima_CS_body_dev_avg_pls_3SD(i) == min(I_minima_CS_body_dev_avg_pls_3SD)
        I_minima_CS_body_dev_avg_pls_3SD_index = i;
    end
end
I_minima_CS_body_tot_number_of_indices = length(I_minima_CS_body_sorted);
%This is the approximation of the tot "area under curve." 
I_minima_CS_body_quotient_within_1SD = (I_minima_CS_body_dev_avg_pls_1SD_index - I_minima_CS_bodydev_avg_min_1SD_index)/I_minima_CS_body_tot_number_of_indices;
I_minima_CS_body_quotient_within_2SD = (I_minima_CS_body_dev_avg_pls_2SD_index - I_minima_dev_avg_min_2SD_index)/I_minima_CS_body_tot_number_of_indices;
I_minima_CS_body_quotient_within_3SD = (I_minima_CS_body_dev_avg_pls_3SD_index - I_minima_CS_body_dev_avg_min_3SD_index)/I_minima_CS_body_tot_number_of_indices;
if (I_minima_CS_body_quotient_within_1SD > 0.66 && I_minima_CS_body_quotient_within_1SD < 0.70) && (I_minima_CS_body_quotient_within_2SD > 0.93 && I_minima_CS_body_quotient_within_2SD < 0.97) && (I_minima_CS_body_quotient_within_3SD > 0.98 && I_minima_CS_body_quotient_within_3SD < 1)
    isNormal(3,11) = 1;
end
if I_minima_CS_body_quotient_within_1SD == 0
    isSpike(3,11) = 1;
end
I_minima_left_septum_sorted = sort(I_minima_left_septum);
I_minima_left_septum_dev_avg_min_1SD = minima(length(I_minima_left_septum_sorted),1);
for i = 1:length(I_minima_left_septum_sorted)
    I_minima_left_septum_dev_avg_min_1SD(i) = I_minima_left_septum_sorted(i) - (I_minima_left_septum_average-I_minima_left_septum_SD);
end
I_minima_left_septum_dev_avg_min_1SD = abs(I_minima_left_septum_dev_avg_min_1SD);
for i = 1:length(I_minima_left_septum_sorted)
    if I_minima_left_septum_dev_avg_min_1SD(i) == min(I_minima_left_septum_dev_avg_min_1SD)
        I_minima_left_septum_dev_avg_min_1SD_index = i;
    end
end
I_minima_left_septum_dev_avg_pls_1SD = minima(length(I_minima_left_septum_sorted),1);
for i = 1:length(I_minima_left_septum_sorted)
    I_minima_left_septum_dev_avg_pls_1SD(i) = I_minima_left_septum_sorted(i) - (I_minima_left_septum_average+I_minima_left_septum_SD);
end
I_minima_left_septum_dev_avg_pls_1SD = abs(I_minima_left_septum_dev_avg_pls_1SD);
for i = 1:length(I_minima_left_septum_sorted)
    if I_minima_left_septum_dev_avg_pls_1SD(i) == min(I_minima_left_septum_dev_avg_pls_1SD)
        I_minima_left_septum_dev_avg_pls_1SD_index = i;
    end
end
I_minima_left_septum_dev_avg_min_2SD = minima(length(I_minima_left_septum_sorted),1);
for i = 1:length(I_minima_left_septum_sorted)
    I_minima_left_septum_dev_avg_min_2SD(i) = I_minima_left_septum_sorted(i) - (I_minima_left_septum_average-2*I_minima_left_septum_SD);
end
I_minima_left_septum_dev_avg_min_2SD = abs(I_minima_left_septum_dev_avg_min_2SD);
for i = 1:length(I_minima_left_septum_sorted)
    if I_minima_left_septum_dev_avg_min_2SD(i) == min(I_minima_left_septum_dev_avg_min_2SD)
        I_minima_dev_avg_min_2SD_index = i;
    end
end
I_minima_left_septum_dev_avg_pls_2SD = minima(length(I_minima_left_septum_sorted),1);
for i = 1:length(I_minima_left_septum_sorted)
    I_minima_left_septum_dev_avg_pls_2SD(i) = I_minima_left_septum_sorted(i) - (I_minima_left_septum_average+2*I_minima_left_septum_SD);
end
I_minima_left_septum_dev_avg_pls_2SD = abs(I_minima_left_septum_dev_avg_pls_2SD);
for i = 1:length(I_minima_left_septum_sorted)
    if I_minima_left_septum_dev_avg_pls_2SD(i) == min(I_minima_left_septum_dev_avg_pls_2SD)
        I_minima_left_septum_dev_avg_pls_2SD_index = i;
    end
end
I_minima_left_septum_dev_avg_min_3SD = minima(length(I_minima_left_septum_sorted),1);
for i = 1:length(I_minima_left_septum_sorted)
    I_minima_left_septum_dev_avg_min_3SD(i) = I_minima_left_septum_sorted(i) - (I_minima_left_septum_average-3*I_minima_left_septum_SD);
end
I_minima_left_septum_dev_avg_min_3SD = abs(I_minima_left_septum_dev_avg_min_3SD);
for i = 1:length(I_minima_left_septum_sorted)
    if I_minima_left_septum_dev_avg_min_3SD(i) == min(I_minima_left_septum_dev_avg_min_3SD)
        I_minima_left_septum_dev_avg_min_3SD_index = i;
    end
end
I_minima_left_septum_dev_avg_pls_3SD = minima(length(I_minima_left_septum_sorted),1);
for i = 1:length(I_minima_left_septum_sorted)
    I_minima_left_septum_dev_avg_pls_3SD(i) = I_minima_left_septum_sorted(i) - (I_minima_left_septum_average+3*I_minima_left_septum_SD);
end
I_minima_left_septum_dev_avg_pls_3SD = abs(I_minima_left_septum_dev_avg_pls_3SD);
for i = 1:length(I_minima_left_septum_sorted)
    if I_minima_left_septum_dev_avg_pls_3SD(i) == min(I_minima_left_septum_dev_avg_pls_3SD)
        I_minima_left_septum_dev_avg_pls_3SD_index = i;
    end
end
I_minima_left_septum_tot_number_of_indices = length(I_minima_left_septum_sorted);
%This is the approximation of the tot "area under curve." 
I_minima_left_septum_quotient_within_1SD = (I_minima_left_septum_dev_avg_pls_1SD_index - I_minima_left_septumdev_avg_min_1SD_index)/I_minima_left_septum_tot_number_of_indices;
I_minima_left_septum_quotient_within_2SD = (I_minima_left_septum_dev_avg_pls_2SD_index - I_minima_dev_avg_min_2SD_index)/I_minima_left_septum_tot_number_of_indices;
I_minima_left_septum_quotient_within_3SD = (I_minima_left_septum_dev_avg_pls_3SD_index - I_minima_left_septum_dev_avg_min_3SD_index)/I_minima_left_septum_tot_number_of_indices;
if (I_minima_left_septum_quotient_within_1SD > 0.66 && I_minima_left_septum_quotient_within_1SD < 0.70) && (I_minima_left_septum_quotient_within_2SD > 0.93 && I_minima_left_septum_quotient_within_2SD < 0.97) && (I_minima_left_septum_quotient_within_3SD > 0.98 && I_minima_left_septum_quotient_within_3SD < 1)
    isNormal(3,12) = 1;
end
if I_minima_left_septum_quotient_within_1SD == 0
    isSpike(3,12) = 1;
end
I_minima_left_atrial_appendage_sorted = sort(I_minima_left_atrial_appendage);
I_minima_left_atrial_appendage_dev_avg_min_1SD = minima(length(I_minima_left_atrial_appendage_sorted),1);
for i = 1:length(I_minima_left_atrial_appendage_sorted)
    I_minima_left_atrial_appendage_dev_avg_min_1SD(i) = I_minima_left_atrial_appendage_sorted(i) - (I_minima_left_atrial_appendage_average-I_minima_left_atrial_appendage_SD);
end
I_minima_left_atrial_appendage_dev_avg_min_1SD = abs(I_minima_left_atrial_appendage_dev_avg_min_1SD);
for i = 1:length(I_minima_left_atrial_appendage_sorted)
    if I_minima_left_atrial_appendage_dev_avg_min_1SD(i) == min(I_minima_left_atrial_appendage_dev_avg_min_1SD)
        I_minima_left_atrial_appendage_dev_avg_min_1SD_index = i;
    end
end
I_minima_left_atrial_appendage_dev_avg_pls_1SD = minima(length(I_minima_left_atrial_appendage_sorted),1);
for i = 1:length(I_minima_left_atrial_appendage_sorted)
    I_minima_left_atrial_appendage_dev_avg_pls_1SD(i) = I_minima_left_atrial_appendage_sorted(i) - (I_minima_left_atrial_appendage_average+I_minima_left_atrial_appendage_SD);
end
I_minima_left_atrial_appendage_dev_avg_pls_1SD = abs(I_minima_left_atrial_appendage_dev_avg_pls_1SD);
for i = 1:length(I_minima_left_atrial_appendage_sorted)
    if I_minima_left_atrial_appendage_dev_avg_pls_1SD(i) == min(I_minima_left_atrial_appendage_dev_avg_pls_1SD)
        I_minima_left_atrial_appendage_dev_avg_pls_1SD_index = i;
    end
end
I_minima_left_atrial_appendage_dev_avg_min_2SD = minima(length(I_minima_left_atrial_appendage_sorted),1);
for i = 1:length(I_minima_left_atrial_appendage_sorted)
    I_minima_left_atrial_appendage_dev_avg_min_2SD(i) = I_minima_left_atrial_appendage_sorted(i) - (I_minima_left_atrial_appendage_average-2*I_minima_left_atrial_appendage_SD);
end
I_minima_left_atrial_appendage_dev_avg_min_2SD = abs(I_minima_left_atrial_appendage_dev_avg_min_2SD);
for i = 1:length(I_minima_left_atrial_appendage_sorted)
    if I_minima_left_atrial_appendage_dev_avg_min_2SD(i) == min(I_minima_left_atrial_appendage_dev_avg_min_2SD)
        I_minima_dev_avg_min_2SD_index = i;
    end
end
I_minima_left_atrial_appendage_dev_avg_pls_2SD = minima(length(I_minima_left_atrial_appendage_sorted),1);
for i = 1:length(I_minima_left_atrial_appendage_sorted)
    I_minima_left_atrial_appendage_dev_avg_pls_2SD(i) = I_minima_left_atrial_appendage_sorted(i) - (I_minima_left_atrial_appendage_average+2*I_minima_left_atrial_appendage_SD);
end
I_minima_left_atrial_appendage_dev_avg_pls_2SD = abs(I_minima_left_atrial_appendage_dev_avg_pls_2SD);
for i = 1:length(I_minima_left_atrial_appendage_sorted)
    if I_minima_left_atrial_appendage_dev_avg_pls_2SD(i) == min(I_minima_left_atrial_appendage_dev_avg_pls_2SD)
        I_minima_left_atrial_appendage_dev_avg_pls_2SD_index = i;
    end
end
I_minima_left_atrial_appendage_dev_avg_min_3SD = minima(length(I_minima_left_atrial_appendage_sorted),1);
for i = 1:length(I_minima_left_atrial_appendage_sorted)
    I_minima_left_atrial_appendage_dev_avg_min_3SD(i) = I_minima_left_atrial_appendage_sorted(i) - (I_minima_left_atrial_appendage_average-3*I_minima_left_atrial_appendage_SD);
end
I_minima_left_atrial_appendage_dev_avg_min_3SD = abs(I_minima_left_atrial_appendage_dev_avg_min_3SD);
for i = 1:length(I_minima_left_atrial_appendage_sorted)
    if I_minima_left_atrial_appendage_dev_avg_min_3SD(i) == min(I_minima_left_atrial_appendage_dev_avg_min_3SD)
        I_minima_left_atrial_appendage_dev_avg_min_3SD_index = i;
    end
end
I_minima_left_atrial_appendage_dev_avg_pls_3SD = minima(length(I_minima_left_atrial_appendage_sorted),1);
for i = 1:length(I_minima_left_atrial_appendage_sorted)
    I_minima_left_atrial_appendage_dev_avg_pls_3SD(i) = I_minima_left_atrial_appendage_sorted(i) - (I_minima_left_atrial_appendage_average+3*I_minima_left_atrial_appendage_SD);
end
I_minima_left_atrial_appendage_dev_avg_pls_3SD = abs(I_minima_left_atrial_appendage_dev_avg_pls_3SD);
for i = 1:length(I_minima_left_atrial_appendage_sorted)
    if I_minima_left_atrial_appendage_dev_avg_pls_3SD(i) == min(I_minima_left_atrial_appendage_dev_avg_pls_3SD)
        I_minima_left_atrial_appendage_dev_avg_pls_3SD_index = i;
    end
end
I_minima_left_atrial_appendage_tot_number_of_indices = length(I_minima_left_atrial_appendage_sorted);
%This is the approximation of the tot "area under curve." 
I_minima_left_atrial_appendage_quotient_within_1SD = (I_minima_left_atrial_appendage_dev_avg_pls_1SD_index - I_minima_left_atrial_appendagedev_avg_min_1SD_index)/I_minima_left_atrial_appendage_tot_number_of_indices;
I_minima_left_atrial_appendage_quotient_within_2SD = (I_minima_left_atrial_appendage_dev_avg_pls_2SD_index - I_minima_dev_avg_min_2SD_index)/I_minima_left_atrial_appendage_tot_number_of_indices;
I_minima_left_atrial_appendage_quotient_within_3SD = (I_minima_left_atrial_appendage_dev_avg_pls_3SD_index - I_minima_left_atrial_appendage_dev_avg_min_3SD_index)/I_minima_left_atrial_appendage_tot_number_of_indices;
if (I_minima_left_atrial_appendage_quotient_within_1SD > 0.66 && I_minima_left_atrial_appendage_quotient_within_1SD < 0.70) && (I_minima_left_atrial_appendage_quotient_within_2SD > 0.93 && I_minima_left_atrial_appendage_quotient_within_2SD < 0.97) && (I_minima_left_atrial_appendage_quotient_within_3SD > 0.98 && I_minima_left_atrial_appendage_quotient_within_3SD < 1)
    isNormal(3,13) = 1;
end
if I_minima_left_atrial_appendage_quotient_within_1SD == 0
    isSpike(3,13) = 1;
end
I_sym_inds_SA_node_sorted = sort(I_sym_inds_SA_node);
I_sym_inds_SA_node_dev_avg_min_1SD = sym_ind(length(I_sym_inds_SA_node_sorted),1);
for i = 1:length(I_sym_inds_SA_node_sorted)
    I_sym_inds_SA_node_dev_avg_min_1SD(i) = I_sym_inds_SA_node_sorted(i) - (I_sym_inds_SA_node_average-I_sym_inds_SA_node_SD);
end
I_sym_inds_SA_node_dev_avg_min_1SD = abs(I_sym_inds_SA_node_dev_avg_min_1SD);
for i = 1:length(I_sym_inds_SA_node_sorted)
    if I_sym_inds_SA_node_dev_avg_min_1SD(i) == min(I_sym_inds_SA_node_dev_avg_min_1SD)
        I_sym_inds_SA_node_dev_avg_min_1SD_index = i;
    end
end
I_sym_inds_SA_node_dev_avg_pls_1SD = sym_ind(length(I_sym_inds_SA_node_sorted),1);
for i = 1:length(I_sym_inds_SA_node_sorted)
    I_sym_inds_SA_node_dev_avg_pls_1SD(i) = I_sym_inds_SA_node_sorted(i) - (I_sym_inds_SA_node_average+I_sym_inds_SA_node_SD);
end
I_sym_inds_SA_node_dev_avg_pls_1SD = abs(I_sym_inds_SA_node_dev_avg_pls_1SD);
for i = 1:length(I_sym_inds_SA_node_sorted)
    if I_sym_inds_SA_node_dev_avg_pls_1SD(i) == min(I_sym_inds_SA_node_dev_avg_pls_1SD)
        I_sym_inds_SA_node_dev_avg_pls_1SD_index = i;
    end
end
I_sym_inds_SA_node_dev_avg_min_2SD = sym_ind(length(I_sym_inds_SA_node_sorted),1);
for i = 1:length(I_sym_inds_SA_node_sorted)
    I_sym_inds_SA_node_dev_avg_min_2SD(i) = I_sym_inds_SA_node_sorted(i) - (I_sym_inds_SA_node_average-2*I_sym_inds_SA_node_SD);
end
I_sym_inds_SA_node_dev_avg_min_2SD = abs(I_sym_inds_SA_node_dev_avg_min_2SD);
for i = 1:length(I_sym_inds_SA_node_sorted)
    if I_sym_inds_SA_node_dev_avg_min_2SD(i) == min(I_sym_inds_SA_node_dev_avg_min_2SD)
        I_sym_inds_dev_avg_min_2SD_index = i;
    end
end
I_sym_inds_SA_node_dev_avg_pls_2SD = sym_ind(length(I_sym_inds_SA_node_sorted),1);
for i = 1:length(I_sym_inds_SA_node_sorted)
    I_sym_inds_SA_node_dev_avg_pls_2SD(i) = I_sym_inds_SA_node_sorted(i) - (I_sym_inds_SA_node_average+2*I_sym_inds_SA_node_SD);
end
I_sym_inds_SA_node_dev_avg_pls_2SD = abs(I_sym_inds_SA_node_dev_avg_pls_2SD);
for i = 1:length(I_sym_inds_SA_node_sorted)
    if I_sym_inds_SA_node_dev_avg_pls_2SD(i) == min(I_sym_inds_SA_node_dev_avg_pls_2SD)
        I_sym_inds_SA_node_dev_avg_pls_2SD_index = i;
    end
end
I_sym_inds_SA_node_dev_avg_min_3SD = sym_ind(length(I_sym_inds_SA_node_sorted),1);
for i = 1:length(I_sym_inds_SA_node_sorted)
    I_sym_inds_SA_node_dev_avg_min_3SD(i) = I_sym_inds_SA_node_sorted(i) - (I_sym_inds_SA_node_average-3*I_sym_inds_SA_node_SD);
end
I_sym_inds_SA_node_dev_avg_min_3SD = abs(I_sym_inds_SA_node_dev_avg_min_3SD);
for i = 1:length(I_sym_inds_SA_node_sorted)
    if I_sym_inds_SA_node_dev_avg_min_3SD(i) == min(I_sym_inds_SA_node_dev_avg_min_3SD)
        I_sym_inds_SA_node_dev_avg_min_3SD_index = i;
    end
end
I_sym_inds_SA_node_dev_avg_pls_3SD = sym_ind(length(I_sym_inds_SA_node_sorted),1);
for i = 1:length(I_sym_inds_SA_node_sorted)
    I_sym_inds_SA_node_dev_avg_pls_3SD(i) = I_sym_inds_SA_node_sorted(i) - (I_sym_inds_SA_node_average+3*I_sym_inds_SA_node_SD);
end
I_sym_inds_SA_node_dev_avg_pls_3SD = abs(I_sym_inds_SA_node_dev_avg_pls_3SD);
for i = 1:length(I_sym_inds_SA_node_sorted)
    if I_sym_inds_SA_node_dev_avg_pls_3SD(i) == min(I_sym_inds_SA_node_dev_avg_pls_3SD)
        I_sym_inds_SA_node_dev_avg_pls_3SD_index = i;
    end
end
I_sym_inds_SA_node_tot_number_of_indices = length(I_sym_inds_SA_node_sorted);
%This is the approximation of the tot "area under curve." 
I_sym_inds_SA_node_quotient_within_1SD = (I_sym_inds_SA_node_dev_avg_pls_1SD_index - I_sym_inds_SA_nodedev_avg_min_1SD_index)/I_sym_inds_SA_node_tot_number_of_indices;
I_sym_inds_SA_node_quotient_within_2SD = (I_sym_inds_SA_node_dev_avg_pls_2SD_index - I_sym_inds_dev_avg_min_2SD_index)/I_sym_inds_SA_node_tot_number_of_indices;
I_sym_inds_SA_node_quotient_within_3SD = (I_sym_inds_SA_node_dev_avg_pls_3SD_index - I_sym_inds_SA_node_dev_avg_min_3SD_index)/I_sym_inds_SA_node_tot_number_of_indices;
if ((I_sym_inds_SA_node_quotient_within_1SD > 0.66 && I_sym_inds_SA_node_quotient_within_1SD < 0.70) && (I_sym_inds_SA_node_quotient_within_2SD > 0.93 && I_sym_inds_SA_node_quotient_within_2SD < 0.97) && (I_sym_inds_SA_node_quotient_within_3SD > 0.98 && I_sym_inds_SA_node_quotient_within_3SD < 1)) 
    isNormal(4,1) = 1;
end
if I_sym_inds_SA_node_quotient_within_1SD == 0
    isSpike(4,1) = 1;
end
I_sym_inds_crista_terminalis_sorted = sort(I_sym_inds_crista_terminalis);
I_sym_inds_crista_terminalis_dev_avg_min_1SD = sym_ind(length(I_sym_inds_crista_terminalis_sorted),1);
for i = 1:length(I_sym_inds_crista_terminalis_sorted)
    I_sym_inds_crista_terminalis_dev_avg_min_1SD(i) = I_sym_inds_crista_terminalis_sorted(i) - (I_sym_inds_crista_terminalis_average-I_sym_inds_crista_terminalis_SD);
end
I_sym_inds_crista_terminalis_dev_avg_min_1SD = abs(I_sym_inds_crista_terminalis_dev_avg_min_1SD);
for i = 1:length(I_sym_inds_crista_terminalis_sorted)
    if I_sym_inds_crista_terminalis_dev_avg_min_1SD(i) == min(I_sym_inds_crista_terminalis_dev_avg_min_1SD)
        I_sym_inds_crista_terminalis_dev_avg_min_1SD_index = i;
    end
end
I_sym_inds_crista_terminalis_dev_avg_pls_1SD = sym_ind(length(I_sym_inds_crista_terminalis_sorted),1);
for i = 1:length(I_sym_inds_crista_terminalis_sorted)
    I_sym_inds_crista_terminalis_dev_avg_pls_1SD(i) = I_sym_inds_crista_terminalis_sorted(i) - (I_sym_inds_crista_terminalis_average+I_sym_inds_crista_terminalis_SD);
end
I_sym_inds_crista_terminalis_dev_avg_pls_1SD = abs(I_sym_inds_crista_terminalis_dev_avg_pls_1SD);
for i = 1:length(I_sym_inds_crista_terminalis_sorted)
    if I_sym_inds_crista_terminalis_dev_avg_pls_1SD(i) == min(I_sym_inds_crista_terminalis_dev_avg_pls_1SD)
        I_sym_inds_crista_terminalis_dev_avg_pls_1SD_index = i;
    end
end
I_sym_inds_crista_terminalis_dev_avg_min_2SD = sym_ind(length(I_sym_inds_crista_terminalis_sorted),1);
for i = 1:length(I_sym_inds_crista_terminalis_sorted)
    I_sym_inds_crista_terminalis_dev_avg_min_2SD(i) = I_sym_inds_crista_terminalis_sorted(i) - (I_sym_inds_crista_terminalis_average-2*I_sym_inds_crista_terminalis_SD);
end
I_sym_inds_crista_terminalis_dev_avg_min_2SD = abs(I_sym_inds_crista_terminalis_dev_avg_min_2SD);
for i = 1:length(I_sym_inds_crista_terminalis_sorted)
    if I_sym_inds_crista_terminalis_dev_avg_min_2SD(i) == min(I_sym_inds_crista_terminalis_dev_avg_min_2SD)
        I_sym_inds_dev_avg_min_2SD_index = i;
    end
end
I_sym_inds_crista_terminalis_dev_avg_pls_2SD = sym_ind(length(I_sym_inds_crista_terminalis_sorted),1);
for i = 1:length(I_sym_inds_crista_terminalis_sorted)
    I_sym_inds_crista_terminalis_dev_avg_pls_2SD(i) = I_sym_inds_crista_terminalis_sorted(i) - (I_sym_inds_crista_terminalis_average+2*I_sym_inds_crista_terminalis_SD);
end
I_sym_inds_crista_terminalis_dev_avg_pls_2SD = abs(I_sym_inds_crista_terminalis_dev_avg_pls_2SD);
for i = 1:length(I_sym_inds_crista_terminalis_sorted)
    if I_sym_inds_crista_terminalis_dev_avg_pls_2SD(i) == min(I_sym_inds_crista_terminalis_dev_avg_pls_2SD)
        I_sym_inds_crista_terminalis_dev_avg_pls_2SD_index = i;
    end
end
I_sym_inds_crista_terminalis_dev_avg_min_3SD = sym_ind(length(I_sym_inds_crista_terminalis_sorted),1);
for i = 1:length(I_sym_inds_crista_terminalis_sorted)
    I_sym_inds_crista_terminalis_dev_avg_min_3SD(i) = I_sym_inds_crista_terminalis_sorted(i) - (I_sym_inds_crista_terminalis_average-3*I_sym_inds_crista_terminalis_SD);
end
I_sym_inds_crista_terminalis_dev_avg_min_3SD = abs(I_sym_inds_crista_terminalis_dev_avg_min_3SD);
for i = 1:length(I_sym_inds_crista_terminalis_sorted)
    if I_sym_inds_crista_terminalis_dev_avg_min_3SD(i) == min(I_sym_inds_crista_terminalis_dev_avg_min_3SD)
        I_sym_inds_crista_terminalis_dev_avg_min_3SD_index = i;
    end
end
I_sym_inds_crista_terminalis_dev_avg_pls_3SD = sym_ind(length(I_sym_inds_crista_terminalis_sorted),1);
for i = 1:length(I_sym_inds_crista_terminalis_sorted)
    I_sym_inds_crista_terminalis_dev_avg_pls_3SD(i) = I_sym_inds_crista_terminalis_sorted(i) - (I_sym_inds_crista_terminalis_average+3*I_sym_inds_crista_terminalis_SD);
end
I_sym_inds_crista_terminalis_dev_avg_pls_3SD = abs(I_sym_inds_crista_terminalis_dev_avg_pls_3SD);
for i = 1:length(I_sym_inds_crista_terminalis_sorted)
    if I_sym_inds_crista_terminalis_dev_avg_pls_3SD(i) == min(I_sym_inds_crista_terminalis_dev_avg_pls_3SD)
        I_sym_inds_crista_terminalis_dev_avg_pls_3SD_index = i;
    end
end
I_sym_inds_crista_terminalis_tot_number_of_indices = length(I_sym_inds_crista_terminalis_sorted);
%This is the approximation of the tot "area under curve." 
I_sym_inds_crista_terminalis_quotient_within_1SD = (I_sym_inds_crista_terminalis_dev_avg_pls_1SD_index - I_sym_inds_crista_terminalisdev_avg_min_1SD_index)/I_sym_inds_crista_terminalis_tot_number_of_indices;
I_sym_inds_crista_terminalis_quotient_within_2SD = (I_sym_inds_crista_terminalis_dev_avg_pls_2SD_index - I_sym_inds_dev_avg_min_2SD_index)/I_sym_inds_crista_terminalis_tot_number_of_indices;
I_sym_inds_crista_terminalis_quotient_within_3SD = (I_sym_inds_crista_terminalis_dev_avg_pls_3SD_index - I_sym_inds_crista_terminalis_dev_avg_min_3SD_index)/I_sym_inds_crista_terminalis_tot_number_of_indices;
if (I_sym_inds_crista_terminalis_quotient_within_1SD > 0.66 && I_sym_inds_crista_terminalis_quotient_within_1SD < 0.70) && (I_sym_inds_crista_terminalis_quotient_within_2SD > 0.93 && I_sym_inds_crista_terminalis_quotient_within_2SD < 0.97) && (I_sym_inds_crista_terminalis_quotient_within_3SD > 0.98 && I_sym_inds_crista_terminalis_quotient_within_3SD < 1)
    isNormal(4,2) = 1;
end
if I_sym_inds_crista_terminalis_quotient_within_1SD == 0
    isSpike(4,2) = 1;
end
I_sym_inds_tricuspid_annulus_sorted = sort(I_sym_inds_tricuspid_annulus);
I_sym_inds_tricuspid_annulus_dev_avg_min_1SD = sym_ind(length(I_sym_inds_tricuspid_annulus_sorted),1);
for i = 1:length(I_sym_inds_tricuspid_annulus_sorted)
    I_sym_inds_tricuspid_annulus_dev_avg_min_1SD(i) = I_sym_inds_tricuspid_annulus_sorted(i) - (I_sym_inds_tricuspid_annulus_average-I_sym_inds_tricuspid_annulus_SD);
end
I_sym_inds_tricuspid_annulus_dev_avg_min_1SD = abs(I_sym_inds_tricuspid_annulus_dev_avg_min_1SD);
for i = 1:length(I_sym_inds_tricuspid_annulus_sorted)
    if I_sym_inds_tricuspid_annulus_dev_avg_min_1SD(i) == min(I_sym_inds_tricuspid_annulus_dev_avg_min_1SD)
        I_sym_inds_tricuspid_annulus_dev_avg_min_1SD_index = i;
    end
end
I_sym_inds_tricuspid_annulus_dev_avg_pls_1SD = sym_ind(length(I_sym_inds_tricuspid_annulus_sorted),1);
for i = 1:length(I_sym_inds_tricuspid_annulus_sorted)
    I_sym_inds_tricuspid_annulus_dev_avg_pls_1SD(i) = I_sym_inds_tricuspid_annulus_sorted(i) - (I_sym_inds_tricuspid_annulus_average+I_sym_inds_tricuspid_annulus_SD);
end
I_sym_inds_tricuspid_annulus_dev_avg_pls_1SD = abs(I_sym_inds_tricuspid_annulus_dev_avg_pls_1SD);
for i = 1:length(I_sym_inds_tricuspid_annulus_sorted)
    if I_sym_inds_tricuspid_annulus_dev_avg_pls_1SD(i) == min(I_sym_inds_tricuspid_annulus_dev_avg_pls_1SD)
        I_sym_inds_tricuspid_annulus_dev_avg_pls_1SD_index = i;
    end
end
I_sym_inds_tricuspid_annulus_dev_avg_min_2SD = sym_ind(length(I_sym_inds_tricuspid_annulus_sorted),1);
for i = 1:length(I_sym_inds_tricuspid_annulus_sorted)
    I_sym_inds_tricuspid_annulus_dev_avg_min_2SD(i) = I_sym_inds_tricuspid_annulus_sorted(i) - (I_sym_inds_tricuspid_annulus_average-2*I_sym_inds_tricuspid_annulus_SD);
end
I_sym_inds_tricuspid_annulus_dev_avg_min_2SD = abs(I_sym_inds_tricuspid_annulus_dev_avg_min_2SD);
for i = 1:length(I_sym_inds_tricuspid_annulus_sorted)
    if I_sym_inds_tricuspid_annulus_dev_avg_min_2SD(i) == min(I_sym_inds_tricuspid_annulus_dev_avg_min_2SD)
        I_sym_inds_dev_avg_min_2SD_index = i;
    end
end
I_sym_inds_tricuspid_annulus_dev_avg_pls_2SD = sym_ind(length(I_sym_inds_tricuspid_annulus_sorted),1);
for i = 1:length(I_sym_inds_tricuspid_annulus_sorted)
    I_sym_inds_tricuspid_annulus_dev_avg_pls_2SD(i) = I_sym_inds_tricuspid_annulus_sorted(i) - (I_sym_inds_tricuspid_annulus_average+2*I_sym_inds_tricuspid_annulus_SD);
end
I_sym_inds_tricuspid_annulus_dev_avg_pls_2SD = abs(I_sym_inds_tricuspid_annulus_dev_avg_pls_2SD);
for i = 1:length(I_sym_inds_tricuspid_annulus_sorted)
    if I_sym_inds_tricuspid_annulus_dev_avg_pls_2SD(i) == min(I_sym_inds_tricuspid_annulus_dev_avg_pls_2SD)
        I_sym_inds_tricuspid_annulus_dev_avg_pls_2SD_index = i;
    end
end
I_sym_inds_tricuspid_annulus_dev_avg_min_3SD = sym_ind(length(I_sym_inds_tricuspid_annulus_sorted),1);
for i = 1:length(I_sym_inds_tricuspid_annulus_sorted)
    I_sym_inds_tricuspid_annulus_dev_avg_min_3SD(i) = I_sym_inds_tricuspid_annulus_sorted(i) - (I_sym_inds_tricuspid_annulus_average-3*I_sym_inds_tricuspid_annulus_SD);
end
I_sym_inds_tricuspid_annulus_dev_avg_min_3SD = abs(I_sym_inds_tricuspid_annulus_dev_avg_min_3SD);
for i = 1:length(I_sym_inds_tricuspid_annulus_sorted)
    if I_sym_inds_tricuspid_annulus_dev_avg_min_3SD(i) == min(I_sym_inds_tricuspid_annulus_dev_avg_min_3SD)
        I_sym_inds_tricuspid_annulus_dev_avg_min_3SD_index = i;
    end
end
I_sym_inds_tricuspid_annulus_dev_avg_pls_3SD = sym_ind(length(I_sym_inds_tricuspid_annulus_sorted),1);
for i = 1:length(I_sym_inds_tricuspid_annulus_sorted)
    I_sym_inds_tricuspid_annulus_dev_avg_pls_3SD(i) = I_sym_inds_tricuspid_annulus_sorted(i) - (I_sym_inds_tricuspid_annulus_average+3*I_sym_inds_tricuspid_annulus_SD);
end
I_sym_inds_tricuspid_annulus_dev_avg_pls_3SD = abs(I_sym_inds_tricuspid_annulus_dev_avg_pls_3SD);
for i = 1:length(I_sym_inds_tricuspid_annulus_sorted)
    if I_sym_inds_tricuspid_annulus_dev_avg_pls_3SD(i) == min(I_sym_inds_tricuspid_annulus_dev_avg_pls_3SD)
        I_sym_inds_tricuspid_annulus_dev_avg_pls_3SD_index = i;
    end
end
I_sym_inds_tricuspid_annulus_tot_number_of_indices = length(I_sym_inds_tricuspid_annulus_sorted);
%This is the approximation of the tot "area under curve." 
I_sym_inds_tricuspid_annulus_quotient_within_1SD = (I_sym_inds_tricuspid_annulus_dev_avg_pls_1SD_index - I_sym_inds_tricuspid_annulusdev_avg_min_1SD_index)/I_sym_inds_tricuspid_annulus_tot_number_of_indices;
I_sym_inds_tricuspid_annulus_quotient_within_2SD = (I_sym_inds_tricuspid_annulus_dev_avg_pls_2SD_index - I_sym_inds_dev_avg_min_2SD_index)/I_sym_inds_tricuspid_annulus_tot_number_of_indices;
I_sym_inds_tricuspid_annulus_quotient_within_3SD = (I_sym_inds_tricuspid_annulus_dev_avg_pls_3SD_index - I_sym_inds_tricuspid_annulus_dev_avg_min_3SD_index)/I_sym_inds_tricuspid_annulus_tot_number_of_indices;
if (I_sym_inds_tricuspid_annulus_quotient_within_1SD > 0.66 && I_sym_inds_tricuspid_annulus_quotient_within_1SD < 0.70) && (I_sym_inds_tricuspid_annulus_quotient_within_2SD > 0.93 && I_sym_inds_tricuspid_annulus_quotient_within_2SD < 0.97) && (I_sym_inds_tricuspid_annulus_quotient_within_3SD > 0.98 && I_sym_inds_tricuspid_annulus_quotient_within_3SD < 1)
    isNormal(4,3) = 1;
end
if I_sym_inds_tricuspid_annulus_quotient_within_1SD == 0
    isSpike(4,3) = 1;
end
I_sym_inds_coronary_sinus_sorted = sort(I_sym_inds_coronary_sinus);
I_sym_inds_coronary_sinus_dev_avg_min_1SD = sym_ind(length(I_sym_inds_coronary_sinus_sorted),1);
for i = 1:length(I_sym_inds_coronary_sinus_sorted)
    I_sym_inds_coronary_sinus_dev_avg_min_1SD(i) = I_sym_inds_coronary_sinus_sorted(i) - (I_sym_inds_coronary_sinus_average-I_sym_inds_coronary_sinus_SD);
end
I_sym_inds_coronary_sinus_dev_avg_min_1SD = abs(I_sym_inds_coronary_sinus_dev_avg_min_1SD);
for i = 1:length(I_sym_inds_coronary_sinus_sorted)
    if I_sym_inds_coronary_sinus_dev_avg_min_1SD(i) == min(I_sym_inds_coronary_sinus_dev_avg_min_1SD)
        I_sym_inds_coronary_sinus_dev_avg_min_1SD_index = i;
    end
end
I_sym_inds_coronary_sinus_dev_avg_pls_1SD = sym_ind(length(I_sym_inds_coronary_sinus_sorted),1);
for i = 1:length(I_sym_inds_coronary_sinus_sorted)
    I_sym_inds_coronary_sinus_dev_avg_pls_1SD(i) = I_sym_inds_coronary_sinus_sorted(i) - (I_sym_inds_coronary_sinus_average+I_sym_inds_coronary_sinus_SD);
end
I_sym_inds_coronary_sinus_dev_avg_pls_1SD = abs(I_sym_inds_coronary_sinus_dev_avg_pls_1SD);
for i = 1:length(I_sym_inds_coronary_sinus_sorted)
    if I_sym_inds_coronary_sinus_dev_avg_pls_1SD(i) == min(I_sym_inds_coronary_sinus_dev_avg_pls_1SD)
        I_sym_inds_coronary_sinus_dev_avg_pls_1SD_index = i;
    end
end
I_sym_inds_coronary_sinus_dev_avg_min_2SD = sym_ind(length(I_sym_inds_coronary_sinus_sorted),1);
for i = 1:length(I_sym_inds_coronary_sinus_sorted)
    I_sym_inds_coronary_sinus_dev_avg_min_2SD(i) = I_sym_inds_coronary_sinus_sorted(i) - (I_sym_inds_coronary_sinus_average-2*I_sym_inds_coronary_sinus_SD);
end
I_sym_inds_coronary_sinus_dev_avg_min_2SD = abs(I_sym_inds_coronary_sinus_dev_avg_min_2SD);
for i = 1:length(I_sym_inds_coronary_sinus_sorted)
    if I_sym_inds_coronary_sinus_dev_avg_min_2SD(i) == min(I_sym_inds_coronary_sinus_dev_avg_min_2SD)
        I_sym_inds_dev_avg_min_2SD_index = i;
    end
end
I_sym_inds_coronary_sinus_dev_avg_pls_2SD = sym_ind(length(I_sym_inds_coronary_sinus_sorted),1);
for i = 1:length(I_sym_inds_coronary_sinus_sorted)
    I_sym_inds_coronary_sinus_dev_avg_pls_2SD(i) = I_sym_inds_coronary_sinus_sorted(i) - (I_sym_inds_coronary_sinus_average+2*I_sym_inds_coronary_sinus_SD);
end
I_sym_inds_coronary_sinus_dev_avg_pls_2SD = abs(I_sym_inds_coronary_sinus_dev_avg_pls_2SD);
for i = 1:length(I_sym_inds_coronary_sinus_sorted)
    if I_sym_inds_coronary_sinus_dev_avg_pls_2SD(i) == min(I_sym_inds_coronary_sinus_dev_avg_pls_2SD)
        I_sym_inds_coronary_sinus_dev_avg_pls_2SD_index = i;
    end
end
I_sym_inds_coronary_sinus_dev_avg_min_3SD = sym_ind(length(I_sym_inds_coronary_sinus_sorted),1);
for i = 1:length(I_sym_inds_coronary_sinus_sorted)
    I_sym_inds_coronary_sinus_dev_avg_min_3SD(i) = I_sym_inds_coronary_sinus_sorted(i) - (I_sym_inds_coronary_sinus_average-3*I_sym_inds_coronary_sinus_SD);
end
I_sym_inds_coronary_sinus_dev_avg_min_3SD = abs(I_sym_inds_coronary_sinus_dev_avg_min_3SD);
for i = 1:length(I_sym_inds_coronary_sinus_sorted)
    if I_sym_inds_coronary_sinus_dev_avg_min_3SD(i) == min(I_sym_inds_coronary_sinus_dev_avg_min_3SD)
        I_sym_inds_coronary_sinus_dev_avg_min_3SD_index = i;
    end
end
I_sym_inds_coronary_sinus_dev_avg_pls_3SD = sym_ind(length(I_sym_inds_coronary_sinus_sorted),1);
for i = 1:length(I_sym_inds_coronary_sinus_sorted)
    I_sym_inds_coronary_sinus_dev_avg_pls_3SD(i) = I_sym_inds_coronary_sinus_sorted(i) - (I_sym_inds_coronary_sinus_average+3*I_sym_inds_coronary_sinus_SD);
end
I_sym_inds_coronary_sinus_dev_avg_pls_3SD = abs(I_sym_inds_coronary_sinus_dev_avg_pls_3SD);
for i = 1:length(I_sym_inds_coronary_sinus_sorted)
    if I_sym_inds_coronary_sinus_dev_avg_pls_3SD(i) == min(I_sym_inds_coronary_sinus_dev_avg_pls_3SD)
        I_sym_inds_coronary_sinus_dev_avg_pls_3SD_index = i;
    end
end
I_sym_inds_coronary_sinus_tot_number_of_indices = length(I_sym_inds_coronary_sinus_sorted);
%This is the approximation of the tot "area under curve." 
I_sym_inds_coronary_sinus_quotient_within_1SD = (I_sym_inds_coronary_sinus_dev_avg_pls_1SD_index - I_sym_inds_coronary_sinusdev_avg_min_1SD_index)/I_sym_inds_coronary_sinus_tot_number_of_indices;
I_sym_inds_coronary_sinus_quotient_within_2SD = (I_sym_inds_coronary_sinus_dev_avg_pls_2SD_index - I_sym_inds_dev_avg_min_2SD_index)/I_sym_inds_coronary_sinus_tot_number_of_indices;
I_sym_inds_coronary_sinus_quotient_within_3SD = (I_sym_inds_coronary_sinus_dev_avg_pls_3SD_index - I_sym_inds_coronary_sinus_dev_avg_min_3SD_index)/I_sym_inds_coronary_sinus_tot_number_of_indices;
if (I_sym_inds_coronary_sinus_quotient_within_1SD > 0.66 && I_sym_inds_coronary_sinus_quotient_within_1SD < 0.70) && (I_sym_inds_coronary_sinus_quotient_within_2SD > 0.93 && I_sym_inds_coronary_sinus_quotient_within_2SD < 0.97) && (I_sym_inds_coronary_sinus_quotient_within_3SD > 0.98 && I_sym_inds_coronary_sinus_quotient_within_3SD < 1)
    isNormal(4,4) = 1;
end
if I_sym_inds_coronary_sinus_quotient_within_1SD == 0
    isSpike(4,4) = 1;
end
I_sym_inds_ostium_sorted = sort(I_sym_inds_ostium);
I_sym_inds_ostium_dev_avg_min_1SD = sym_ind(length(I_sym_inds_ostium_sorted),1);
for i = 1:length(I_sym_inds_ostium_sorted)
    I_sym_inds_ostium_dev_avg_min_1SD(i) = I_sym_inds_ostium_sorted(i) - (I_sym_inds_ostium_average-I_sym_inds_ostium_SD);
end
I_sym_inds_ostium_dev_avg_min_1SD = abs(I_sym_inds_ostium_dev_avg_min_1SD);
for i = 1:length(I_sym_inds_ostium_sorted)
    if I_sym_inds_ostium_dev_avg_min_1SD(i) == min(I_sym_inds_ostium_dev_avg_min_1SD)
        I_sym_inds_ostium_dev_avg_min_1SD_index = i;
    end
end
I_sym_inds_ostium_dev_avg_pls_1SD = sym_ind(length(I_sym_inds_ostium_sorted),1);
for i = 1:length(I_sym_inds_ostium_sorted)
    I_sym_inds_ostium_dev_avg_pls_1SD(i) = I_sym_inds_ostium_sorted(i) - (I_sym_inds_ostium_average+I_sym_inds_ostium_SD);
end
I_sym_inds_ostium_dev_avg_pls_1SD = abs(I_sym_inds_ostium_dev_avg_pls_1SD);
for i = 1:length(I_sym_inds_ostium_sorted)
    if I_sym_inds_ostium_dev_avg_pls_1SD(i) == min(I_sym_inds_ostium_dev_avg_pls_1SD)
        I_sym_inds_ostium_dev_avg_pls_1SD_index = i;
    end
end
I_sym_inds_ostium_dev_avg_min_2SD = sym_ind(length(I_sym_inds_ostium_sorted),1);
for i = 1:length(I_sym_inds_ostium_sorted)
    I_sym_inds_ostium_dev_avg_min_2SD(i) = I_sym_inds_ostium_sorted(i) - (I_sym_inds_ostium_average-2*I_sym_inds_ostium_SD);
end
I_sym_inds_ostium_dev_avg_min_2SD = abs(I_sym_inds_ostium_dev_avg_min_2SD);
for i = 1:length(I_sym_inds_ostium_sorted)
    if I_sym_inds_ostium_dev_avg_min_2SD(i) == min(I_sym_inds_ostium_dev_avg_min_2SD)
        I_sym_inds_dev_avg_min_2SD_index = i;
    end
end
I_sym_inds_ostium_dev_avg_pls_2SD = sym_ind(length(I_sym_inds_ostium_sorted),1);
for i = 1:length(I_sym_inds_ostium_sorted)
    I_sym_inds_ostium_dev_avg_pls_2SD(i) = I_sym_inds_ostium_sorted(i) - (I_sym_inds_ostium_average+2*I_sym_inds_ostium_SD);
end
I_sym_inds_ostium_dev_avg_pls_2SD = abs(I_sym_inds_ostium_dev_avg_pls_2SD);
for i = 1:length(I_sym_inds_ostium_sorted)
    if I_sym_inds_ostium_dev_avg_pls_2SD(i) == min(I_sym_inds_ostium_dev_avg_pls_2SD)
        I_sym_inds_ostium_dev_avg_pls_2SD_index = i;
    end
end
I_sym_inds_ostium_dev_avg_min_3SD = sym_ind(length(I_sym_inds_ostium_sorted),1);
for i = 1:length(I_sym_inds_ostium_sorted)
    I_sym_inds_ostium_dev_avg_min_3SD(i) = I_sym_inds_ostium_sorted(i) - (I_sym_inds_ostium_average-3*I_sym_inds_ostium_SD);
end
I_sym_inds_ostium_dev_avg_min_3SD = abs(I_sym_inds_ostium_dev_avg_min_3SD);
for i = 1:length(I_sym_inds_ostium_sorted)
    if I_sym_inds_ostium_dev_avg_min_3SD(i) == min(I_sym_inds_ostium_dev_avg_min_3SD)
        I_sym_inds_ostium_dev_avg_min_3SD_index = i;
    end
end
I_sym_inds_ostium_dev_avg_pls_3SD = sym_ind(length(I_sym_inds_ostium_sorted),1);
for i = 1:length(I_sym_inds_ostium_sorted)
    I_sym_inds_ostium_dev_avg_pls_3SD(i) = I_sym_inds_ostium_sorted(i) - (I_sym_inds_ostium_average+3*I_sym_inds_ostium_SD);
end
I_sym_inds_ostium_dev_avg_pls_3SD = abs(I_sym_inds_ostium_dev_avg_pls_3SD);
for i = 1:length(I_sym_inds_ostium_sorted)
    if I_sym_inds_ostium_dev_avg_pls_3SD(i) == min(I_sym_inds_ostium_dev_avg_pls_3SD)
        I_sym_inds_ostium_dev_avg_pls_3SD_index = i;
    end
end
I_sym_inds_ostium_tot_number_of_indices = length(I_sym_inds_ostium_sorted);
%This is the approximation of the tot "area under curve." 
I_sym_inds_ostium_quotient_within_1SD = (I_sym_inds_ostium_dev_avg_pls_1SD_index - I_sym_inds_ostiumdev_avg_min_1SD_index)/I_sym_inds_ostium_tot_number_of_indices;
I_sym_inds_ostium_quotient_within_2SD = (I_sym_inds_ostium_dev_avg_pls_2SD_index - I_sym_inds_dev_avg_min_2SD_index)/I_sym_inds_ostium_tot_number_of_indices;
I_sym_inds_ostium_quotient_within_3SD = (I_sym_inds_ostium_dev_avg_pls_3SD_index - I_sym_inds_ostium_dev_avg_min_3SD_index)/I_sym_inds_ostium_tot_number_of_indices;
if (I_sym_inds_ostium_quotient_within_1SD > 0.66 && I_sym_inds_ostium_quotient_within_1SD < 0.70) && (I_sym_inds_ostium_quotient_within_2SD > 0.93 && I_sym_inds_ostium_quotient_within_2SD < 0.97) && (I_sym_inds_ostium_quotient_within_3SD > 0.98 && I_sym_inds_ostium_quotient_within_3SD < 1)
    isNormal(4,5) = 1;
end
if I_sym_inds_ostium_quotient_within_1SD == 0
    isSpike(4,5) = 1;
end
I_sym_inds_perinodal_sorted = sort(I_sym_inds_perinodal);
I_sym_inds_perinodal_dev_avg_min_1SD = sym_ind(length(I_sym_inds_perinodal_sorted),1);
for i = 1:length(I_sym_inds_perinodal_sorted)
    I_sym_inds_perinodal_dev_avg_min_1SD(i) = I_sym_inds_perinodal_sorted(i) - (I_sym_inds_perinodal_average-I_sym_inds_perinodal_SD);
end
I_sym_inds_perinodal_dev_avg_min_1SD = abs(I_sym_inds_perinodal_dev_avg_min_1SD);
for i = 1:length(I_sym_inds_perinodal_sorted)
    if I_sym_inds_perinodal_dev_avg_min_1SD(i) == min(I_sym_inds_perinodal_dev_avg_min_1SD)
        I_sym_inds_perinodal_dev_avg_min_1SD_index = i;
    end
end
I_sym_inds_perinodal_dev_avg_pls_1SD = sym_ind(length(I_sym_inds_perinodal_sorted),1);
for i = 1:length(I_sym_inds_perinodal_sorted)
    I_sym_inds_perinodal_dev_avg_pls_1SD(i) = I_sym_inds_perinodal_sorted(i) - (I_sym_inds_perinodal_average+I_sym_inds_perinodal_SD);
end
I_sym_inds_perinodal_dev_avg_pls_1SD = abs(I_sym_inds_perinodal_dev_avg_pls_1SD);
for i = 1:length(I_sym_inds_perinodal_sorted)
    if I_sym_inds_perinodal_dev_avg_pls_1SD(i) == min(I_sym_inds_perinodal_dev_avg_pls_1SD)
        I_sym_inds_perinodal_dev_avg_pls_1SD_index = i;
    end
end
I_sym_inds_perinodal_dev_avg_min_2SD = sym_ind(length(I_sym_inds_perinodal_sorted),1);
for i = 1:length(I_sym_inds_perinodal_sorted)
    I_sym_inds_perinodal_dev_avg_min_2SD(i) = I_sym_inds_perinodal_sorted(i) - (I_sym_inds_perinodal_average-2*I_sym_inds_perinodal_SD);
end
I_sym_inds_perinodal_dev_avg_min_2SD = abs(I_sym_inds_perinodal_dev_avg_min_2SD);
for i = 1:length(I_sym_inds_perinodal_sorted)
    if I_sym_inds_perinodal_dev_avg_min_2SD(i) == min(I_sym_inds_perinodal_dev_avg_min_2SD)
        I_sym_inds_dev_avg_min_2SD_index = i;
    end
end
I_sym_inds_perinodal_dev_avg_pls_2SD = sym_ind(length(I_sym_inds_perinodal_sorted),1);
for i = 1:length(I_sym_inds_perinodal_sorted)
    I_sym_inds_perinodal_dev_avg_pls_2SD(i) = I_sym_inds_perinodal_sorted(i) - (I_sym_inds_perinodal_average+2*I_sym_inds_perinodal_SD);
end
I_sym_inds_perinodal_dev_avg_pls_2SD = abs(I_sym_inds_perinodal_dev_avg_pls_2SD);
for i = 1:length(I_sym_inds_perinodal_sorted)
    if I_sym_inds_perinodal_dev_avg_pls_2SD(i) == min(I_sym_inds_perinodal_dev_avg_pls_2SD)
        I_sym_inds_perinodal_dev_avg_pls_2SD_index = i;
    end
end
I_sym_inds_perinodal_dev_avg_min_3SD = sym_ind(length(I_sym_inds_perinodal_sorted),1);
for i = 1:length(I_sym_inds_perinodal_sorted)
    I_sym_inds_perinodal_dev_avg_min_3SD(i) = I_sym_inds_perinodal_sorted(i) - (I_sym_inds_perinodal_average-3*I_sym_inds_perinodal_SD);
end
I_sym_inds_perinodal_dev_avg_min_3SD = abs(I_sym_inds_perinodal_dev_avg_min_3SD);
for i = 1:length(I_sym_inds_perinodal_sorted)
    if I_sym_inds_perinodal_dev_avg_min_3SD(i) == min(I_sym_inds_perinodal_dev_avg_min_3SD)
        I_sym_inds_perinodal_dev_avg_min_3SD_index = i;
    end
end
I_sym_inds_perinodal_dev_avg_pls_3SD = sym_ind(length(I_sym_inds_perinodal_sorted),1);
for i = 1:length(I_sym_inds_perinodal_sorted)
    I_sym_inds_perinodal_dev_avg_pls_3SD(i) = I_sym_inds_perinodal_sorted(i) - (I_sym_inds_perinodal_average+3*I_sym_inds_perinodal_SD);
end
I_sym_inds_perinodal_dev_avg_pls_3SD = abs(I_sym_inds_perinodal_dev_avg_pls_3SD);
for i = 1:length(I_sym_inds_perinodal_sorted)
    if I_sym_inds_perinodal_dev_avg_pls_3SD(i) == min(I_sym_inds_perinodal_dev_avg_pls_3SD)
        I_sym_inds_perinodal_dev_avg_pls_3SD_index = i;
    end
end
I_sym_inds_perinodal_tot_number_of_indices = length(I_sym_inds_perinodal_sorted);
%This is the approximation of the tot "area under curve." 
I_sym_inds_perinodal_quotient_within_1SD = (I_sym_inds_perinodal_dev_avg_pls_1SD_index - I_sym_inds_perinodaldev_avg_min_1SD_index)/I_sym_inds_perinodal_tot_number_of_indices;
I_sym_inds_perinodal_quotient_within_2SD = (I_sym_inds_perinodal_dev_avg_pls_2SD_index - I_sym_inds_dev_avg_min_2SD_index)/I_sym_inds_perinodal_tot_number_of_indices;
I_sym_inds_perinodal_quotient_within_3SD = (I_sym_inds_perinodal_dev_avg_pls_3SD_index - I_sym_inds_perinodal_dev_avg_min_3SD_index)/I_sym_inds_perinodal_tot_number_of_indices;
if (I_sym_inds_perinodal_quotient_within_1SD > 0.66 && I_sym_inds_perinodal_quotient_within_1SD < 0.70) && (I_sym_inds_perinodal_quotient_within_2SD > 0.93 && I_sym_inds_perinodal_quotient_within_2SD < 0.97) && (I_sym_inds_perinodal_quotient_within_3SD > 0.98 && I_sym_inds_perinodal_quotient_within_3SD < 1)
    isNormal(4,6) = 1;
end
if I_sym_inds_perinodal_quotient_within_1SD == 0
    isSpike(4,6) = 1;
end
I_sym_inds_right_septum_sorted = sort(I_sym_inds_right_septum);
I_sym_inds_right_septum_dev_avg_min_1SD = sym_ind(length(I_sym_inds_right_septum_sorted),1);
for i = 1:length(I_sym_inds_right_septum_sorted)
    I_sym_inds_right_septum_dev_avg_min_1SD(i) = I_sym_inds_right_septum_sorted(i) - (I_sym_inds_right_septum_average-I_sym_inds_right_septum_SD);
end
I_sym_inds_right_septum_dev_avg_min_1SD = abs(I_sym_inds_right_septum_dev_avg_min_1SD);
for i = 1:length(I_sym_inds_right_septum_sorted)
    if I_sym_inds_right_septum_dev_avg_min_1SD(i) == min(I_sym_inds_right_septum_dev_avg_min_1SD)
        I_sym_inds_right_septum_dev_avg_min_1SD_index = i;
    end
end
I_sym_inds_right_septum_dev_avg_pls_1SD = sym_ind(length(I_sym_inds_right_septum_sorted),1);
for i = 1:length(I_sym_inds_right_septum_sorted)
    I_sym_inds_right_septum_dev_avg_pls_1SD(i) = I_sym_inds_right_septum_sorted(i) - (I_sym_inds_right_septum_average+I_sym_inds_right_septum_SD);
end
I_sym_inds_right_septum_dev_avg_pls_1SD = abs(I_sym_inds_right_septum_dev_avg_pls_1SD);
for i = 1:length(I_sym_inds_right_septum_sorted)
    if I_sym_inds_right_septum_dev_avg_pls_1SD(i) == min(I_sym_inds_right_septum_dev_avg_pls_1SD)
        I_sym_inds_right_septum_dev_avg_pls_1SD_index = i;
    end
end
I_sym_inds_right_septum_dev_avg_min_2SD = sym_ind(length(I_sym_inds_right_septum_sorted),1);
for i = 1:length(I_sym_inds_right_septum_sorted)
    I_sym_inds_right_septum_dev_avg_min_2SD(i) = I_sym_inds_right_septum_sorted(i) - (I_sym_inds_right_septum_average-2*I_sym_inds_right_septum_SD);
end
I_sym_inds_right_septum_dev_avg_min_2SD = abs(I_sym_inds_right_septum_dev_avg_min_2SD);
for i = 1:length(I_sym_inds_right_septum_sorted)
    if I_sym_inds_right_septum_dev_avg_min_2SD(i) == min(I_sym_inds_right_septum_dev_avg_min_2SD)
        I_sym_inds_dev_avg_min_2SD_index = i;
    end
end
I_sym_inds_right_septum_dev_avg_pls_2SD = sym_ind(length(I_sym_inds_right_septum_sorted),1);
for i = 1:length(I_sym_inds_right_septum_sorted)
    I_sym_inds_right_septum_dev_avg_pls_2SD(i) = I_sym_inds_right_septum_sorted(i) - (I_sym_inds_right_septum_average+2*I_sym_inds_right_septum_SD);
end
I_sym_inds_right_septum_dev_avg_pls_2SD = abs(I_sym_inds_right_septum_dev_avg_pls_2SD);
for i = 1:length(I_sym_inds_right_septum_sorted)
    if I_sym_inds_right_septum_dev_avg_pls_2SD(i) == min(I_sym_inds_right_septum_dev_avg_pls_2SD)
        I_sym_inds_right_septum_dev_avg_pls_2SD_index = i;
    end
end
I_sym_inds_right_septum_dev_avg_min_3SD = sym_ind(length(I_sym_inds_right_septum_sorted),1);
for i = 1:length(I_sym_inds_right_septum_sorted)
    I_sym_inds_right_septum_dev_avg_min_3SD(i) = I_sym_inds_right_septum_sorted(i) - (I_sym_inds_right_septum_average-3*I_sym_inds_right_septum_SD);
end
I_sym_inds_right_septum_dev_avg_min_3SD = abs(I_sym_inds_right_septum_dev_avg_min_3SD);
for i = 1:length(I_sym_inds_right_septum_sorted)
    if I_sym_inds_right_septum_dev_avg_min_3SD(i) == min(I_sym_inds_right_septum_dev_avg_min_3SD)
        I_sym_inds_right_septum_dev_avg_min_3SD_index = i;
    end
end
I_sym_inds_right_septum_dev_avg_pls_3SD = sym_ind(length(I_sym_inds_right_septum_sorted),1);
for i = 1:length(I_sym_inds_right_septum_sorted)
    I_sym_inds_right_septum_dev_avg_pls_3SD(i) = I_sym_inds_right_septum_sorted(i) - (I_sym_inds_right_septum_average+3*I_sym_inds_right_septum_SD);
end
I_sym_inds_right_septum_dev_avg_pls_3SD = abs(I_sym_inds_right_septum_dev_avg_pls_3SD);
for i = 1:length(I_sym_inds_right_septum_sorted)
    if I_sym_inds_right_septum_dev_avg_pls_3SD(i) == min(I_sym_inds_right_septum_dev_avg_pls_3SD)
        I_sym_inds_right_septum_dev_avg_pls_3SD_index = i;
    end
end
I_sym_inds_right_septum_tot_number_of_indices = length(I_sym_inds_right_septum_sorted);
%This is the approximation of the tot "area under curve." 
I_sym_inds_right_septum_quotient_within_1SD = (I_sym_inds_right_septum_dev_avg_pls_1SD_index - I_sym_inds_right_septumdev_avg_min_1SD_index)/I_sym_inds_right_septum_tot_number_of_indices;
I_sym_inds_right_septum_quotient_within_2SD = (I_sym_inds_right_septum_dev_avg_pls_2SD_index - I_sym_inds_dev_avg_min_2SD_index)/I_sym_inds_right_septum_tot_number_of_indices;
I_sym_inds_right_septum_quotient_within_3SD = (I_sym_inds_right_septum_dev_avg_pls_3SD_index - I_sym_inds_right_septum_dev_avg_min_3SD_index)/I_sym_inds_right_septum_tot_number_of_indices;
if (I_sym_inds_right_septum_quotient_within_1SD > 0.66 && I_sym_inds_right_septum_quotient_within_1SD < 0.70) && (I_sym_inds_right_septum_quotient_within_2SD > 0.93 && I_sym_inds_right_septum_quotient_within_2SD < 0.97) && (I_sym_inds_right_septum_quotient_within_3SD > 0.98 && I_sym_inds_right_septum_quotient_within_3SD < 1)
    isNormal(4,7) = 1;
end
if I_sym_inds_right_septum_quotient_within_1SD == 0
    isSpike(4,7) = 1;
end
I_sym_inds_right_atrial_appendage_sorted = sort(I_sym_inds_right_atrial_appendage);
I_sym_inds_right_atrial_appendage_dev_avg_min_1SD = sym_ind(length(I_sym_inds_right_atrial_appendage_sorted),1);
for i = 1:length(I_sym_inds_right_atrial_appendage_sorted)
    I_sym_inds_right_atrial_appendage_dev_avg_min_1SD(i) = I_sym_inds_right_atrial_appendage_sorted(i) - (I_sym_inds_right_atrial_appendage_average-I_sym_inds_right_atrial_appendage_SD);
end
I_sym_inds_right_atrial_appendage_dev_avg_min_1SD = abs(I_sym_inds_right_atrial_appendage_dev_avg_min_1SD);
for i = 1:length(I_sym_inds_right_atrial_appendage_sorted)
    if I_sym_inds_right_atrial_appendage_dev_avg_min_1SD(i) == min(I_sym_inds_right_atrial_appendage_dev_avg_min_1SD)
        I_sym_inds_right_atrial_appendage_dev_avg_min_1SD_index = i;
    end
end
I_sym_inds_right_atrial_appendage_dev_avg_pls_1SD = sym_ind(length(I_sym_inds_right_atrial_appendage_sorted),1);
for i = 1:length(I_sym_inds_right_atrial_appendage_sorted)
    I_sym_inds_right_atrial_appendage_dev_avg_pls_1SD(i) = I_sym_inds_right_atrial_appendage_sorted(i) - (I_sym_inds_right_atrial_appendage_average+I_sym_inds_right_atrial_appendage_SD);
end
I_sym_inds_right_atrial_appendage_dev_avg_pls_1SD = abs(I_sym_inds_right_atrial_appendage_dev_avg_pls_1SD);
for i = 1:length(I_sym_inds_right_atrial_appendage_sorted)
    if I_sym_inds_right_atrial_appendage_dev_avg_pls_1SD(i) == min(I_sym_inds_right_atrial_appendage_dev_avg_pls_1SD)
        I_sym_inds_right_atrial_appendage_dev_avg_pls_1SD_index = i;
    end
end
I_sym_inds_right_atrial_appendage_dev_avg_min_2SD = sym_ind(length(I_sym_inds_right_atrial_appendage_sorted),1);
for i = 1:length(I_sym_inds_right_atrial_appendage_sorted)
    I_sym_inds_right_atrial_appendage_dev_avg_min_2SD(i) = I_sym_inds_right_atrial_appendage_sorted(i) - (I_sym_inds_right_atrial_appendage_average-2*I_sym_inds_right_atrial_appendage_SD);
end
I_sym_inds_right_atrial_appendage_dev_avg_min_2SD = abs(I_sym_inds_right_atrial_appendage_dev_avg_min_2SD);
for i = 1:length(I_sym_inds_right_atrial_appendage_sorted)
    if I_sym_inds_right_atrial_appendage_dev_avg_min_2SD(i) == min(I_sym_inds_right_atrial_appendage_dev_avg_min_2SD)
        I_sym_inds_dev_avg_min_2SD_index = i;
    end
end
I_sym_inds_right_atrial_appendage_dev_avg_pls_2SD = sym_ind(length(I_sym_inds_right_atrial_appendage_sorted),1);
for i = 1:length(I_sym_inds_right_atrial_appendage_sorted)
    I_sym_inds_right_atrial_appendage_dev_avg_pls_2SD(i) = I_sym_inds_right_atrial_appendage_sorted(i) - (I_sym_inds_right_atrial_appendage_average+2*I_sym_inds_right_atrial_appendage_SD);
end
I_sym_inds_right_atrial_appendage_dev_avg_pls_2SD = abs(I_sym_inds_right_atrial_appendage_dev_avg_pls_2SD);
for i = 1:length(I_sym_inds_right_atrial_appendage_sorted)
    if I_sym_inds_right_atrial_appendage_dev_avg_pls_2SD(i) == min(I_sym_inds_right_atrial_appendage_dev_avg_pls_2SD)
        I_sym_inds_right_atrial_appendage_dev_avg_pls_2SD_index = i;
    end
end
I_sym_inds_right_atrial_appendage_dev_avg_min_3SD = sym_ind(length(I_sym_inds_right_atrial_appendage_sorted),1);
for i = 1:length(I_sym_inds_right_atrial_appendage_sorted)
    I_sym_inds_right_atrial_appendage_dev_avg_min_3SD(i) = I_sym_inds_right_atrial_appendage_sorted(i) - (I_sym_inds_right_atrial_appendage_average-3*I_sym_inds_right_atrial_appendage_SD);
end
I_sym_inds_right_atrial_appendage_dev_avg_min_3SD = abs(I_sym_inds_right_atrial_appendage_dev_avg_min_3SD);
for i = 1:length(I_sym_inds_right_atrial_appendage_sorted)
    if I_sym_inds_right_atrial_appendage_dev_avg_min_3SD(i) == min(I_sym_inds_right_atrial_appendage_dev_avg_min_3SD)
        I_sym_inds_right_atrial_appendage_dev_avg_min_3SD_index = i;
    end
end
I_sym_inds_right_atrial_appendage_dev_avg_pls_3SD = sym_ind(length(I_sym_inds_right_atrial_appendage_sorted),1);
for i = 1:length(I_sym_inds_right_atrial_appendage_sorted)
    I_sym_inds_right_atrial_appendage_dev_avg_pls_3SD(i) = I_sym_inds_right_atrial_appendage_sorted(i) - (I_sym_inds_right_atrial_appendage_average+3*I_sym_inds_right_atrial_appendage_SD);
end
I_sym_inds_right_atrial_appendage_dev_avg_pls_3SD = abs(I_sym_inds_right_atrial_appendage_dev_avg_pls_3SD);
for i = 1:length(I_sym_inds_right_atrial_appendage_sorted)
    if I_sym_inds_right_atrial_appendage_dev_avg_pls_3SD(i) == min(I_sym_inds_right_atrial_appendage_dev_avg_pls_3SD)
        I_sym_inds_right_atrial_appendage_dev_avg_pls_3SD_index = i;
    end
end
I_sym_inds_right_atrial_appendage_tot_number_of_indices = length(I_sym_inds_right_atrial_appendage_sorted);
%This is the approximation of the tot "area under curve." 
I_sym_inds_right_atrial_appendage_quotient_within_1SD = (I_sym_inds_right_atrial_appendage_dev_avg_pls_1SD_index - I_sym_inds_right_atrial_appendagedev_avg_min_1SD_index)/I_sym_inds_right_atrial_appendage_tot_number_of_indices;
I_sym_inds_right_atrial_appendage_quotient_within_2SD = (I_sym_inds_right_atrial_appendage_dev_avg_pls_2SD_index - I_sym_inds_dev_avg_min_2SD_index)/I_sym_inds_right_atrial_appendage_tot_number_of_indices;
I_sym_inds_right_atrial_appendage_quotient_within_3SD = (I_sym_inds_right_atrial_appendage_dev_avg_pls_3SD_index - I_sym_inds_right_atrial_appendage_dev_avg_min_3SD_index)/I_sym_inds_right_atrial_appendage_tot_number_of_indices;
if (I_sym_inds_right_atrial_appendage_quotient_within_1SD > 0.66 && I_sym_inds_right_atrial_appendage_quotient_within_1SD < 0.70) && (I_sym_inds_right_atrial_appendage_quotient_within_2SD > 0.93 && I_sym_inds_right_atrial_appendage_quotient_within_2SD < 0.97) && (I_sym_inds_right_atrial_appendage_quotient_within_3SD > 0.98 && I_sym_inds_right_atrial_appendage_quotient_within_3SD < 1)
    isNormal(4,8) = 1;
end
if I_sym_inds_right_atrial_appendage_quotient_within_1SD == 0
    isSpike(4,8) = 1;
end
I_sym_inds_pulmonary_veins_sorted = sort(I_sym_inds_pulmonary_veins);
I_sym_inds_pulmonary_veins_dev_avg_min_1SD = sym_ind(length(I_sym_inds_pulmonary_veins_sorted),1);
for i = 1:length(I_sym_inds_pulmonary_veins_sorted)
    I_sym_inds_pulmonary_veins_dev_avg_min_1SD(i) = I_sym_inds_pulmonary_veins_sorted(i) - (I_sym_inds_pulmonary_veins_average-I_sym_inds_pulmonary_veins_SD);
end
I_sym_inds_pulmonary_veins_dev_avg_min_1SD = abs(I_sym_inds_pulmonary_veins_dev_avg_min_1SD);
for i = 1:length(I_sym_inds_pulmonary_veins_sorted)
    if I_sym_inds_pulmonary_veins_dev_avg_min_1SD(i) == min(I_sym_inds_pulmonary_veins_dev_avg_min_1SD)
        I_sym_inds_pulmonary_veins_dev_avg_min_1SD_index = i;
    end
end
I_sym_inds_pulmonary_veins_dev_avg_pls_1SD = sym_ind(length(I_sym_inds_pulmonary_veins_sorted),1);
for i = 1:length(I_sym_inds_pulmonary_veins_sorted)
    I_sym_inds_pulmonary_veins_dev_avg_pls_1SD(i) = I_sym_inds_pulmonary_veins_sorted(i) - (I_sym_inds_pulmonary_veins_average+I_sym_inds_pulmonary_veins_SD);
end
I_sym_inds_pulmonary_veins_dev_avg_pls_1SD = abs(I_sym_inds_pulmonary_veins_dev_avg_pls_1SD);
for i = 1:length(I_sym_inds_pulmonary_veins_sorted)
    if I_sym_inds_pulmonary_veins_dev_avg_pls_1SD(i) == min(I_sym_inds_pulmonary_veins_dev_avg_pls_1SD)
        I_sym_inds_pulmonary_veins_dev_avg_pls_1SD_index = i;
    end
end
I_sym_inds_pulmonary_veins_dev_avg_min_2SD = sym_ind(length(I_sym_inds_pulmonary_veins_sorted),1);
for i = 1:length(I_sym_inds_pulmonary_veins_sorted)
    I_sym_inds_pulmonary_veins_dev_avg_min_2SD(i) = I_sym_inds_pulmonary_veins_sorted(i) - (I_sym_inds_pulmonary_veins_average-2*I_sym_inds_pulmonary_veins_SD);
end
I_sym_inds_pulmonary_veins_dev_avg_min_2SD = abs(I_sym_inds_pulmonary_veins_dev_avg_min_2SD);
for i = 1:length(I_sym_inds_pulmonary_veins_sorted)
    if I_sym_inds_pulmonary_veins_dev_avg_min_2SD(i) == min(I_sym_inds_pulmonary_veins_dev_avg_min_2SD)
        I_sym_inds_dev_avg_min_2SD_index = i;
    end
end
I_sym_inds_pulmonary_veins_dev_avg_pls_2SD = sym_ind(length(I_sym_inds_pulmonary_veins_sorted),1);
for i = 1:length(I_sym_inds_pulmonary_veins_sorted)
    I_sym_inds_pulmonary_veins_dev_avg_pls_2SD(i) = I_sym_inds_pulmonary_veins_sorted(i) - (I_sym_inds_pulmonary_veins_average+2*I_sym_inds_pulmonary_veins_SD);
end
I_sym_inds_pulmonary_veins_dev_avg_pls_2SD = abs(I_sym_inds_pulmonary_veins_dev_avg_pls_2SD);
for i = 1:length(I_sym_inds_pulmonary_veins_sorted)
    if I_sym_inds_pulmonary_veins_dev_avg_pls_2SD(i) == min(I_sym_inds_pulmonary_veins_dev_avg_pls_2SD)
        I_sym_inds_pulmonary_veins_dev_avg_pls_2SD_index = i;
    end
end
I_sym_inds_pulmonary_veins_dev_avg_min_3SD = sym_ind(length(I_sym_inds_pulmonary_veins_sorted),1);
for i = 1:length(I_sym_inds_pulmonary_veins_sorted)
    I_sym_inds_pulmonary_veins_dev_avg_min_3SD(i) = I_sym_inds_pulmonary_veins_sorted(i) - (I_sym_inds_pulmonary_veins_average-3*I_sym_inds_pulmonary_veins_SD);
end
I_sym_inds_pulmonary_veins_dev_avg_min_3SD = abs(I_sym_inds_pulmonary_veins_dev_avg_min_3SD);
for i = 1:length(I_sym_inds_pulmonary_veins_sorted)
    if I_sym_inds_pulmonary_veins_dev_avg_min_3SD(i) == min(I_sym_inds_pulmonary_veins_dev_avg_min_3SD)
        I_sym_inds_pulmonary_veins_dev_avg_min_3SD_index = i;
    end
end
I_sym_inds_pulmonary_veins_dev_avg_pls_3SD = sym_ind(length(I_sym_inds_pulmonary_veins_sorted),1);
for i = 1:length(I_sym_inds_pulmonary_veins_sorted)
    I_sym_inds_pulmonary_veins_dev_avg_pls_3SD(i) = I_sym_inds_pulmonary_veins_sorted(i) - (I_sym_inds_pulmonary_veins_average+3*I_sym_inds_pulmonary_veins_SD);
end
I_sym_inds_pulmonary_veins_dev_avg_pls_3SD = abs(I_sym_inds_pulmonary_veins_dev_avg_pls_3SD);
for i = 1:length(I_sym_inds_pulmonary_veins_sorted)
    if I_sym_inds_pulmonary_veins_dev_avg_pls_3SD(i) == min(I_sym_inds_pulmonary_veins_dev_avg_pls_3SD)
        I_sym_inds_pulmonary_veins_dev_avg_pls_3SD_index = i;
    end
end
I_sym_inds_pulmonary_veins_tot_number_of_indices = length(I_sym_inds_pulmonary_veins_sorted);
%This is the approximation of the tot "area under curve." 
I_sym_inds_pulmonary_veins_quotient_within_1SD = (I_sym_inds_pulmonary_veins_dev_avg_pls_1SD_index - I_sym_inds_pulmonary_veinsdev_avg_min_1SD_index)/I_sym_inds_pulmonary_veins_tot_number_of_indices;
I_sym_inds_pulmonary_veins_quotient_within_2SD = (I_sym_inds_pulmonary_veins_dev_avg_pls_2SD_index - I_sym_inds_dev_avg_min_2SD_index)/I_sym_inds_pulmonary_veins_tot_number_of_indices;
I_sym_inds_pulmonary_veins_quotient_within_3SD = (I_sym_inds_pulmonary_veins_dev_avg_pls_3SD_index - I_sym_inds_pulmonary_veins_dev_avg_min_3SD_index)/I_sym_inds_pulmonary_veins_tot_number_of_indices;
if (I_sym_inds_pulmonary_veins_quotient_within_1SD > 0.66 && I_sym_inds_pulmonary_veins_quotient_within_1SD < 0.70) && (I_sym_inds_pulmonary_veins_quotient_within_2SD > 0.93 && I_sym_inds_pulmonary_veins_quotient_within_2SD < 0.97) && (I_sym_inds_pulmonary_veins_quotient_within_3SD > 0.98 && I_sym_inds_pulmonary_veins_quotient_within_3SD < 1)
    isNormal(4,9) = 1;
end
if I_sym_inds_pulmonary_veins_quotient_within_1SD == 0
    isSpike(4,9) = 1;
end
I_sym_inds_mitral_annulus_sorted = sort(I_sym_inds_mitral_annulus);
I_sym_inds_mitral_annulus_dev_avg_min_1SD = sym_ind(length(I_sym_inds_mitral_annulus_sorted),1);
for i = 1:length(I_sym_inds_mitral_annulus_sorted)
    I_sym_inds_mitral_annulus_dev_avg_min_1SD(i) = I_sym_inds_mitral_annulus_sorted(i) - (I_sym_inds_mitral_annulus_average-I_sym_inds_mitral_annulus_SD);
end
I_sym_inds_mitral_annulus_dev_avg_min_1SD = abs(I_sym_inds_mitral_annulus_dev_avg_min_1SD);
for i = 1:length(I_sym_inds_mitral_annulus_sorted)
    if I_sym_inds_mitral_annulus_dev_avg_min_1SD(i) == min(I_sym_inds_mitral_annulus_dev_avg_min_1SD)
        I_sym_inds_mitral_annulus_dev_avg_min_1SD_index = i;
    end
end
I_sym_inds_mitral_annulus_dev_avg_pls_1SD = sym_ind(length(I_sym_inds_mitral_annulus_sorted),1);
for i = 1:length(I_sym_inds_mitral_annulus_sorted)
    I_sym_inds_mitral_annulus_dev_avg_pls_1SD(i) = I_sym_inds_mitral_annulus_sorted(i) - (I_sym_inds_mitral_annulus_average+I_sym_inds_mitral_annulus_SD);
end
I_sym_inds_mitral_annulus_dev_avg_pls_1SD = abs(I_sym_inds_mitral_annulus_dev_avg_pls_1SD);
for i = 1:length(I_sym_inds_mitral_annulus_sorted)
    if I_sym_inds_mitral_annulus_dev_avg_pls_1SD(i) == min(I_sym_inds_mitral_annulus_dev_avg_pls_1SD)
        I_sym_inds_mitral_annulus_dev_avg_pls_1SD_index = i;
    end
end
I_sym_inds_mitral_annulus_dev_avg_min_2SD = sym_ind(length(I_sym_inds_mitral_annulus_sorted),1);
for i = 1:length(I_sym_inds_mitral_annulus_sorted)
    I_sym_inds_mitral_annulus_dev_avg_min_2SD(i) = I_sym_inds_mitral_annulus_sorted(i) - (I_sym_inds_mitral_annulus_average-2*I_sym_inds_mitral_annulus_SD);
end
I_sym_inds_mitral_annulus_dev_avg_min_2SD = abs(I_sym_inds_mitral_annulus_dev_avg_min_2SD);
for i = 1:length(I_sym_inds_mitral_annulus_sorted)
    if I_sym_inds_mitral_annulus_dev_avg_min_2SD(i) == min(I_sym_inds_mitral_annulus_dev_avg_min_2SD)
        I_sym_inds_dev_avg_min_2SD_index = i;
    end
end
I_sym_inds_mitral_annulus_dev_avg_pls_2SD = sym_ind(length(I_sym_inds_mitral_annulus_sorted),1);
for i = 1:length(I_sym_inds_mitral_annulus_sorted)
    I_sym_inds_mitral_annulus_dev_avg_pls_2SD(i) = I_sym_inds_mitral_annulus_sorted(i) - (I_sym_inds_mitral_annulus_average+2*I_sym_inds_mitral_annulus_SD);
end
I_sym_inds_mitral_annulus_dev_avg_pls_2SD = abs(I_sym_inds_mitral_annulus_dev_avg_pls_2SD);
for i = 1:length(I_sym_inds_mitral_annulus_sorted)
    if I_sym_inds_mitral_annulus_dev_avg_pls_2SD(i) == min(I_sym_inds_mitral_annulus_dev_avg_pls_2SD)
        I_sym_inds_mitral_annulus_dev_avg_pls_2SD_index = i;
    end
end
I_sym_inds_mitral_annulus_dev_avg_min_3SD = sym_ind(length(I_sym_inds_mitral_annulus_sorted),1);
for i = 1:length(I_sym_inds_mitral_annulus_sorted)
    I_sym_inds_mitral_annulus_dev_avg_min_3SD(i) = I_sym_inds_mitral_annulus_sorted(i) - (I_sym_inds_mitral_annulus_average-3*I_sym_inds_mitral_annulus_SD);
end
I_sym_inds_mitral_annulus_dev_avg_min_3SD = abs(I_sym_inds_mitral_annulus_dev_avg_min_3SD);
for i = 1:length(I_sym_inds_mitral_annulus_sorted)
    if I_sym_inds_mitral_annulus_dev_avg_min_3SD(i) == min(I_sym_inds_mitral_annulus_dev_avg_min_3SD)
        I_sym_inds_mitral_annulus_dev_avg_min_3SD_index = i;
    end
end
I_sym_inds_mitral_annulus_dev_avg_pls_3SD = sym_ind(length(I_sym_inds_mitral_annulus_sorted),1);
for i = 1:length(I_sym_inds_mitral_annulus_sorted)
    I_sym_inds_mitral_annulus_dev_avg_pls_3SD(i) = I_sym_inds_mitral_annulus_sorted(i) - (I_sym_inds_mitral_annulus_average+3*I_sym_inds_mitral_annulus_SD);
end
I_sym_inds_mitral_annulus_dev_avg_pls_3SD = abs(I_sym_inds_mitral_annulus_dev_avg_pls_3SD);
for i = 1:length(I_sym_inds_mitral_annulus_sorted)
    if I_sym_inds_mitral_annulus_dev_avg_pls_3SD(i) == min(I_sym_inds_mitral_annulus_dev_avg_pls_3SD)
        I_sym_inds_mitral_annulus_dev_avg_pls_3SD_index = i;
    end
end
I_sym_inds_mitral_annulus_tot_number_of_indices = length(I_sym_inds_mitral_annulus_sorted);
%This is the approximation of the tot "area under curve." 
I_sym_inds_mitral_annulus_quotient_within_1SD = (I_sym_inds_mitral_annulus_dev_avg_pls_1SD_index - I_sym_inds_mitral_annulusdev_avg_min_1SD_index)/I_sym_inds_mitral_annulus_tot_number_of_indices;
I_sym_inds_mitral_annulus_quotient_within_2SD = (I_sym_inds_mitral_annulus_dev_avg_pls_2SD_index - I_sym_inds_dev_avg_min_2SD_index)/I_sym_inds_mitral_annulus_tot_number_of_indices;
I_sym_inds_mitral_annulus_quotient_within_3SD = (I_sym_inds_mitral_annulus_dev_avg_pls_3SD_index - I_sym_inds_mitral_annulus_dev_avg_min_3SD_index)/I_sym_inds_mitral_annulus_tot_number_of_indices;
if (I_sym_inds_mitral_annulus_quotient_within_1SD > 0.66 && I_sym_inds_mitral_annulus_quotient_within_1SD < 0.70) && (I_sym_inds_mitral_annulus_quotient_within_2SD > 0.93 && I_sym_inds_mitral_annulus_quotient_within_2SD < 0.97) && (I_sym_inds_mitral_annulus_quotient_within_3SD > 0.98 && I_sym_inds_mitral_annulus_quotient_within_3SD < 1)
    isNormal(4,10) = 1;
end
if I_sym_inds_mitral_annulus_quotient_within_1SD == 0
    isSpike(4,10) = 1;
end
I_sym_inds_CS_body_sorted = sort(I_sym_inds_CS_body);
I_sym_inds_CS_body_dev_avg_min_1SD = sym_ind(length(I_sym_inds_CS_body_sorted),1);
for i = 1:length(I_sym_inds_CS_body_sorted)
    I_sym_inds_CS_body_dev_avg_min_1SD(i) = I_sym_inds_CS_body_sorted(i) - (I_sym_inds_CS_body_average-I_sym_inds_CS_body_SD);
end
I_sym_inds_CS_body_dev_avg_min_1SD = abs(I_sym_inds_CS_body_dev_avg_min_1SD);
for i = 1:length(I_sym_inds_CS_body_sorted)
    if I_sym_inds_CS_body_dev_avg_min_1SD(i) == min(I_sym_inds_CS_body_dev_avg_min_1SD)
        I_sym_inds_CS_body_dev_avg_min_1SD_index = i;
    end
end
I_sym_inds_CS_body_dev_avg_pls_1SD = sym_ind(length(I_sym_inds_CS_body_sorted),1);
for i = 1:length(I_sym_inds_CS_body_sorted)
    I_sym_inds_CS_body_dev_avg_pls_1SD(i) = I_sym_inds_CS_body_sorted(i) - (I_sym_inds_CS_body_average+I_sym_inds_CS_body_SD);
end
I_sym_inds_CS_body_dev_avg_pls_1SD = abs(I_sym_inds_CS_body_dev_avg_pls_1SD);
for i = 1:length(I_sym_inds_CS_body_sorted)
    if I_sym_inds_CS_body_dev_avg_pls_1SD(i) == min(I_sym_inds_CS_body_dev_avg_pls_1SD)
        I_sym_inds_CS_body_dev_avg_pls_1SD_index = i;
    end
end
I_sym_inds_CS_body_dev_avg_min_2SD = sym_ind(length(I_sym_inds_CS_body_sorted),1);
for i = 1:length(I_sym_inds_CS_body_sorted)
    I_sym_inds_CS_body_dev_avg_min_2SD(i) = I_sym_inds_CS_body_sorted(i) - (I_sym_inds_CS_body_average-2*I_sym_inds_CS_body_SD);
end
I_sym_inds_CS_body_dev_avg_min_2SD = abs(I_sym_inds_CS_body_dev_avg_min_2SD);
for i = 1:length(I_sym_inds_CS_body_sorted)
    if I_sym_inds_CS_body_dev_avg_min_2SD(i) == min(I_sym_inds_CS_body_dev_avg_min_2SD)
        I_sym_inds_dev_avg_min_2SD_index = i;
    end
end
I_sym_inds_CS_body_dev_avg_pls_2SD = sym_ind(length(I_sym_inds_CS_body_sorted),1);
for i = 1:length(I_sym_inds_CS_body_sorted)
    I_sym_inds_CS_body_dev_avg_pls_2SD(i) = I_sym_inds_CS_body_sorted(i) - (I_sym_inds_CS_body_average+2*I_sym_inds_CS_body_SD);
end
I_sym_inds_CS_body_dev_avg_pls_2SD = abs(I_sym_inds_CS_body_dev_avg_pls_2SD);
for i = 1:length(I_sym_inds_CS_body_sorted)
    if I_sym_inds_CS_body_dev_avg_pls_2SD(i) == min(I_sym_inds_CS_body_dev_avg_pls_2SD)
        I_sym_inds_CS_body_dev_avg_pls_2SD_index = i;
    end
end
I_sym_inds_CS_body_dev_avg_min_3SD = sym_ind(length(I_sym_inds_CS_body_sorted),1);
for i = 1:length(I_sym_inds_CS_body_sorted)
    I_sym_inds_CS_body_dev_avg_min_3SD(i) = I_sym_inds_CS_body_sorted(i) - (I_sym_inds_CS_body_average-3*I_sym_inds_CS_body_SD);
end
I_sym_inds_CS_body_dev_avg_min_3SD = abs(I_sym_inds_CS_body_dev_avg_min_3SD);
for i = 1:length(I_sym_inds_CS_body_sorted)
    if I_sym_inds_CS_body_dev_avg_min_3SD(i) == min(I_sym_inds_CS_body_dev_avg_min_3SD)
        I_sym_inds_CS_body_dev_avg_min_3SD_index = i;
    end
end
I_sym_inds_CS_body_dev_avg_pls_3SD = sym_ind(length(I_sym_inds_CS_body_sorted),1);
for i = 1:length(I_sym_inds_CS_body_sorted)
    I_sym_inds_CS_body_dev_avg_pls_3SD(i) = I_sym_inds_CS_body_sorted(i) - (I_sym_inds_CS_body_average+3*I_sym_inds_CS_body_SD);
end
I_sym_inds_CS_body_dev_avg_pls_3SD = abs(I_sym_inds_CS_body_dev_avg_pls_3SD);
for i = 1:length(I_sym_inds_CS_body_sorted)
    if I_sym_inds_CS_body_dev_avg_pls_3SD(i) == min(I_sym_inds_CS_body_dev_avg_pls_3SD)
        I_sym_inds_CS_body_dev_avg_pls_3SD_index = i;
    end
end
I_sym_inds_CS_body_tot_number_of_indices = length(I_sym_inds_CS_body_sorted);
%This is the approximation of the tot "area under curve." 
I_sym_inds_CS_body_quotient_within_1SD = (I_sym_inds_CS_body_dev_avg_pls_1SD_index - I_sym_inds_CS_bodydev_avg_min_1SD_index)/I_sym_inds_CS_body_tot_number_of_indices;
I_sym_inds_CS_body_quotient_within_2SD = (I_sym_inds_CS_body_dev_avg_pls_2SD_index - I_sym_inds_dev_avg_min_2SD_index)/I_sym_inds_CS_body_tot_number_of_indices;
I_sym_inds_CS_body_quotient_within_3SD = (I_sym_inds_CS_body_dev_avg_pls_3SD_index - I_sym_inds_CS_body_dev_avg_min_3SD_index)/I_sym_inds_CS_body_tot_number_of_indices;
if (I_sym_inds_CS_body_quotient_within_1SD > 0.66 && I_sym_inds_CS_body_quotient_within_1SD < 0.70) && (I_sym_inds_CS_body_quotient_within_2SD > 0.93 && I_sym_inds_CS_body_quotient_within_2SD < 0.97) && (I_sym_inds_CS_body_quotient_within_3SD > 0.98 && I_sym_inds_CS_body_quotient_within_3SD < 1)
    isNormal(4,11) = 1;
end
if I_sym_inds_CS_body_quotient_within_1SD == 0
    isSpike(4,11) = 1;
end
I_sym_inds_left_septum_sorted = sort(I_sym_inds_left_septum);
I_sym_inds_left_septum_dev_avg_min_1SD = sym_ind(length(I_sym_inds_left_septum_sorted),1);
for i = 1:length(I_sym_inds_left_septum_sorted)
    I_sym_inds_left_septum_dev_avg_min_1SD(i) = I_sym_inds_left_septum_sorted(i) - (I_sym_inds_left_septum_average-I_sym_inds_left_septum_SD);
end
I_sym_inds_left_septum_dev_avg_min_1SD = abs(I_sym_inds_left_septum_dev_avg_min_1SD);
for i = 1:length(I_sym_inds_left_septum_sorted)
    if I_sym_inds_left_septum_dev_avg_min_1SD(i) == min(I_sym_inds_left_septum_dev_avg_min_1SD)
        I_sym_inds_left_septum_dev_avg_min_1SD_index = i;
    end
end
I_sym_inds_left_septum_dev_avg_pls_1SD = sym_ind(length(I_sym_inds_left_septum_sorted),1);
for i = 1:length(I_sym_inds_left_septum_sorted)
    I_sym_inds_left_septum_dev_avg_pls_1SD(i) = I_sym_inds_left_septum_sorted(i) - (I_sym_inds_left_septum_average+I_sym_inds_left_septum_SD);
end
I_sym_inds_left_septum_dev_avg_pls_1SD = abs(I_sym_inds_left_septum_dev_avg_pls_1SD);
for i = 1:length(I_sym_inds_left_septum_sorted)
    if I_sym_inds_left_septum_dev_avg_pls_1SD(i) == min(I_sym_inds_left_septum_dev_avg_pls_1SD)
        I_sym_inds_left_septum_dev_avg_pls_1SD_index = i;
    end
end
I_sym_inds_left_septum_dev_avg_min_2SD = sym_ind(length(I_sym_inds_left_septum_sorted),1);
for i = 1:length(I_sym_inds_left_septum_sorted)
    I_sym_inds_left_septum_dev_avg_min_2SD(i) = I_sym_inds_left_septum_sorted(i) - (I_sym_inds_left_septum_average-2*I_sym_inds_left_septum_SD);
end
I_sym_inds_left_septum_dev_avg_min_2SD = abs(I_sym_inds_left_septum_dev_avg_min_2SD);
for i = 1:length(I_sym_inds_left_septum_sorted)
    if I_sym_inds_left_septum_dev_avg_min_2SD(i) == min(I_sym_inds_left_septum_dev_avg_min_2SD)
        I_sym_inds_dev_avg_min_2SD_index = i;
    end
end
I_sym_inds_left_septum_dev_avg_pls_2SD = sym_ind(length(I_sym_inds_left_septum_sorted),1);
for i = 1:length(I_sym_inds_left_septum_sorted)
    I_sym_inds_left_septum_dev_avg_pls_2SD(i) = I_sym_inds_left_septum_sorted(i) - (I_sym_inds_left_septum_average+2*I_sym_inds_left_septum_SD);
end
I_sym_inds_left_septum_dev_avg_pls_2SD = abs(I_sym_inds_left_septum_dev_avg_pls_2SD);
for i = 1:length(I_sym_inds_left_septum_sorted)
    if I_sym_inds_left_septum_dev_avg_pls_2SD(i) == min(I_sym_inds_left_septum_dev_avg_pls_2SD)
        I_sym_inds_left_septum_dev_avg_pls_2SD_index = i;
    end
end
I_sym_inds_left_septum_dev_avg_min_3SD = sym_ind(length(I_sym_inds_left_septum_sorted),1);
for i = 1:length(I_sym_inds_left_septum_sorted)
    I_sym_inds_left_septum_dev_avg_min_3SD(i) = I_sym_inds_left_septum_sorted(i) - (I_sym_inds_left_septum_average-3*I_sym_inds_left_septum_SD);
end
I_sym_inds_left_septum_dev_avg_min_3SD = abs(I_sym_inds_left_septum_dev_avg_min_3SD);
for i = 1:length(I_sym_inds_left_septum_sorted)
    if I_sym_inds_left_septum_dev_avg_min_3SD(i) == min(I_sym_inds_left_septum_dev_avg_min_3SD)
        I_sym_inds_left_septum_dev_avg_min_3SD_index = i;
    end
end
I_sym_inds_left_septum_dev_avg_pls_3SD = sym_ind(length(I_sym_inds_left_septum_sorted),1);
for i = 1:length(I_sym_inds_left_septum_sorted)
    I_sym_inds_left_septum_dev_avg_pls_3SD(i) = I_sym_inds_left_septum_sorted(i) - (I_sym_inds_left_septum_average+3*I_sym_inds_left_septum_SD);
end
I_sym_inds_left_septum_dev_avg_pls_3SD = abs(I_sym_inds_left_septum_dev_avg_pls_3SD);
for i = 1:length(I_sym_inds_left_septum_sorted)
    if I_sym_inds_left_septum_dev_avg_pls_3SD(i) == min(I_sym_inds_left_septum_dev_avg_pls_3SD)
        I_sym_inds_left_septum_dev_avg_pls_3SD_index = i;
    end
end
I_sym_inds_left_septum_tot_number_of_indices = length(I_sym_inds_left_septum_sorted);
%This is the approximation of the tot "area under curve." 
I_sym_inds_left_septum_quotient_within_1SD = (I_sym_inds_left_septum_dev_avg_pls_1SD_index - I_sym_inds_left_septumdev_avg_min_1SD_index)/I_sym_inds_left_septum_tot_number_of_indices;
I_sym_inds_left_septum_quotient_within_2SD = (I_sym_inds_left_septum_dev_avg_pls_2SD_index - I_sym_inds_dev_avg_min_2SD_index)/I_sym_inds_left_septum_tot_number_of_indices;
I_sym_inds_left_septum_quotient_within_3SD = (I_sym_inds_left_septum_dev_avg_pls_3SD_index - I_sym_inds_left_septum_dev_avg_min_3SD_index)/I_sym_inds_left_septum_tot_number_of_indices;
if (I_sym_inds_left_septum_quotient_within_1SD > 0.66 && I_sym_inds_left_septum_quotient_within_1SD < 0.70) && (I_sym_inds_left_septum_quotient_within_2SD > 0.93 && I_sym_inds_left_septum_quotient_within_2SD < 0.97) && (I_sym_inds_left_septum_quotient_within_3SD > 0.98 && I_sym_inds_left_septum_quotient_within_3SD < 1)
    isNormal(4,12) = 1;
end
if I_sym_inds_left_septum_quotient_within_1SD == 0
    isSpike(4,12) = 1;
end
I_sym_inds_left_atrial_appendage_sorted = sort(I_sym_inds_left_atrial_appendage);
I_sym_inds_left_atrial_appendage_dev_avg_min_1SD = sym_ind(length(I_sym_inds_left_atrial_appendage_sorted),1);
for i = 1:length(I_sym_inds_left_atrial_appendage_sorted)
    I_sym_inds_left_atrial_appendage_dev_avg_min_1SD(i) = I_sym_inds_left_atrial_appendage_sorted(i) - (I_sym_inds_left_atrial_appendage_average-I_sym_inds_left_atrial_appendage_SD);
end
I_sym_inds_left_atrial_appendage_dev_avg_min_1SD = abs(I_sym_inds_left_atrial_appendage_dev_avg_min_1SD);
for i = 1:length(I_sym_inds_left_atrial_appendage_sorted)
    if I_sym_inds_left_atrial_appendage_dev_avg_min_1SD(i) == min(I_sym_inds_left_atrial_appendage_dev_avg_min_1SD)
        I_sym_inds_left_atrial_appendage_dev_avg_min_1SD_index = i;
    end
end
I_sym_inds_left_atrial_appendage_dev_avg_pls_1SD = sym_ind(length(I_sym_inds_left_atrial_appendage_sorted),1);
for i = 1:length(I_sym_inds_left_atrial_appendage_sorted)
    I_sym_inds_left_atrial_appendage_dev_avg_pls_1SD(i) = I_sym_inds_left_atrial_appendage_sorted(i) - (I_sym_inds_left_atrial_appendage_average+I_sym_inds_left_atrial_appendage_SD);
end
I_sym_inds_left_atrial_appendage_dev_avg_pls_1SD = abs(I_sym_inds_left_atrial_appendage_dev_avg_pls_1SD);
for i = 1:length(I_sym_inds_left_atrial_appendage_sorted)
    if I_sym_inds_left_atrial_appendage_dev_avg_pls_1SD(i) == min(I_sym_inds_left_atrial_appendage_dev_avg_pls_1SD)
        I_sym_inds_left_atrial_appendage_dev_avg_pls_1SD_index = i;
    end
end
I_sym_inds_left_atrial_appendage_dev_avg_min_2SD = sym_ind(length(I_sym_inds_left_atrial_appendage_sorted),1);
for i = 1:length(I_sym_inds_left_atrial_appendage_sorted)
    I_sym_inds_left_atrial_appendage_dev_avg_min_2SD(i) = I_sym_inds_left_atrial_appendage_sorted(i) - (I_sym_inds_left_atrial_appendage_average-2*I_sym_inds_left_atrial_appendage_SD);
end
I_sym_inds_left_atrial_appendage_dev_avg_min_2SD = abs(I_sym_inds_left_atrial_appendage_dev_avg_min_2SD);
for i = 1:length(I_sym_inds_left_atrial_appendage_sorted)
    if I_sym_inds_left_atrial_appendage_dev_avg_min_2SD(i) == min(I_sym_inds_left_atrial_appendage_dev_avg_min_2SD)
        I_sym_inds_dev_avg_min_2SD_index = i;
    end
end
I_sym_inds_left_atrial_appendage_dev_avg_pls_2SD = sym_ind(length(I_sym_inds_left_atrial_appendage_sorted),1);
for i = 1:length(I_sym_inds_left_atrial_appendage_sorted)
    I_sym_inds_left_atrial_appendage_dev_avg_pls_2SD(i) = I_sym_inds_left_atrial_appendage_sorted(i) - (I_sym_inds_left_atrial_appendage_average+2*I_sym_inds_left_atrial_appendage_SD);
end
I_sym_inds_left_atrial_appendage_dev_avg_pls_2SD = abs(I_sym_inds_left_atrial_appendage_dev_avg_pls_2SD);
for i = 1:length(I_sym_inds_left_atrial_appendage_sorted)
    if I_sym_inds_left_atrial_appendage_dev_avg_pls_2SD(i) == min(I_sym_inds_left_atrial_appendage_dev_avg_pls_2SD)
        I_sym_inds_left_atrial_appendage_dev_avg_pls_2SD_index = i;
    end
end
I_sym_inds_left_atrial_appendage_dev_avg_min_3SD = sym_ind(length(I_sym_inds_left_atrial_appendage_sorted),1);
for i = 1:length(I_sym_inds_left_atrial_appendage_sorted)
    I_sym_inds_left_atrial_appendage_dev_avg_min_3SD(i) = I_sym_inds_left_atrial_appendage_sorted(i) - (I_sym_inds_left_atrial_appendage_average-3*I_sym_inds_left_atrial_appendage_SD);
end
I_sym_inds_left_atrial_appendage_dev_avg_min_3SD = abs(I_sym_inds_left_atrial_appendage_dev_avg_min_3SD);
for i = 1:length(I_sym_inds_left_atrial_appendage_sorted)
    if I_sym_inds_left_atrial_appendage_dev_avg_min_3SD(i) == min(I_sym_inds_left_atrial_appendage_dev_avg_min_3SD)
        I_sym_inds_left_atrial_appendage_dev_avg_min_3SD_index = i;
    end
end
I_sym_inds_left_atrial_appendage_dev_avg_pls_3SD = sym_ind(length(I_sym_inds_left_atrial_appendage_sorted),1);
for i = 1:length(I_sym_inds_left_atrial_appendage_sorted)
    I_sym_inds_left_atrial_appendage_dev_avg_pls_3SD(i) = I_sym_inds_left_atrial_appendage_sorted(i) - (I_sym_inds_left_atrial_appendage_average+3*I_sym_inds_left_atrial_appendage_SD);
end
I_sym_inds_left_atrial_appendage_dev_avg_pls_3SD = abs(I_sym_inds_left_atrial_appendage_dev_avg_pls_3SD);
for i = 1:length(I_sym_inds_left_atrial_appendage_sorted)
    if I_sym_inds_left_atrial_appendage_dev_avg_pls_3SD(i) == min(I_sym_inds_left_atrial_appendage_dev_avg_pls_3SD)
        I_sym_inds_left_atrial_appendage_dev_avg_pls_3SD_index = i;
    end
end
I_sym_inds_left_atrial_appendage_tot_number_of_indices = length(I_sym_inds_left_atrial_appendage_sorted);
%This is the approximation of the tot "area under curve." 
I_sym_inds_left_atrial_appendage_quotient_within_1SD = (I_sym_inds_left_atrial_appendage_dev_avg_pls_1SD_index - I_sym_inds_left_atrial_appendagedev_avg_min_1SD_index)/I_sym_inds_left_atrial_appendage_tot_number_of_indices;
I_sym_inds_left_atrial_appendage_quotient_within_2SD = (I_sym_inds_left_atrial_appendage_dev_avg_pls_2SD_index - I_sym_inds_dev_avg_min_2SD_index)/I_sym_inds_left_atrial_appendage_tot_number_of_indices;
I_sym_inds_left_atrial_appendage_quotient_within_3SD = (I_sym_inds_left_atrial_appendage_dev_avg_pls_3SD_index - I_sym_inds_left_atrial_appendage_dev_avg_min_3SD_index)/I_sym_inds_left_atrial_appendage_tot_number_of_indices;
if (I_sym_inds_left_atrial_appendage_quotient_within_1SD > 0.66 && I_sym_inds_left_atrial_appendage_quotient_within_1SD < 0.70) && (I_sym_inds_left_atrial_appendage_quotient_within_2SD > 0.93 && I_sym_inds_left_atrial_appendage_quotient_within_2SD < 0.97) && (I_sym_inds_left_atrial_appendage_quotient_within_3SD > 0.98 && I_sym_inds_left_atrial_appendage_quotient_within_3SD < 1)
    isNormal(4,13) = 1;
end
if I_sym_inds_left_atrial_appendage_quotient_within_1SD == 0
    isSpike(4,13) = 1;
end
I_concavs_SA_node_sorted = sort(I_concavs_SA_node);
I_concavs_SA_node_dev_avg_min_1SD = concav(length(I_concavs_SA_node_sorted),1);
for i = 1:length(I_concavs_SA_node_sorted)
    I_concavs_SA_node_dev_avg_min_1SD(i) = I_concavs_SA_node_sorted(i) - (I_concavs_SA_node_average-I_concavs_SA_node_SD);
end
I_concavs_SA_node_dev_avg_min_1SD = abs(I_concavs_SA_node_dev_avg_min_1SD);
for i = 1:length(I_concavs_SA_node_sorted)
    if I_concavs_SA_node_dev_avg_min_1SD(i) == min(I_concavs_SA_node_dev_avg_min_1SD)
        I_concavs_SA_node_dev_avg_min_1SD_index = i;
    end
end
I_concavs_SA_node_dev_avg_pls_1SD = concav(length(I_concavs_SA_node_sorted),1);
for i = 1:length(I_concavs_SA_node_sorted)
    I_concavs_SA_node_dev_avg_pls_1SD(i) = I_concavs_SA_node_sorted(i) - (I_concavs_SA_node_average+I_concavs_SA_node_SD);
end
I_concavs_SA_node_dev_avg_pls_1SD = abs(I_concavs_SA_node_dev_avg_pls_1SD);
for i = 1:length(I_concavs_SA_node_sorted)
    if I_concavs_SA_node_dev_avg_pls_1SD(i) == min(I_concavs_SA_node_dev_avg_pls_1SD)
        I_concavs_SA_node_dev_avg_pls_1SD_index = i;
    end
end
I_concavs_SA_node_dev_avg_min_2SD = concav(length(I_concavs_SA_node_sorted),1);
for i = 1:length(I_concavs_SA_node_sorted)
    I_concavs_SA_node_dev_avg_min_2SD(i) = I_concavs_SA_node_sorted(i) - (I_concavs_SA_node_average-2*I_concavs_SA_node_SD);
end
I_concavs_SA_node_dev_avg_min_2SD = abs(I_concavs_SA_node_dev_avg_min_2SD);
for i = 1:length(I_concavs_SA_node_sorted)
    if I_concavs_SA_node_dev_avg_min_2SD(i) == min(I_concavs_SA_node_dev_avg_min_2SD)
        I_concavs_dev_avg_min_2SD_index = i;
    end
end
I_concavs_SA_node_dev_avg_pls_2SD = concav(length(I_concavs_SA_node_sorted),1);
for i = 1:length(I_concavs_SA_node_sorted)
    I_concavs_SA_node_dev_avg_pls_2SD(i) = I_concavs_SA_node_sorted(i) - (I_concavs_SA_node_average+2*I_concavs_SA_node_SD);
end
I_concavs_SA_node_dev_avg_pls_2SD = abs(I_concavs_SA_node_dev_avg_pls_2SD);
for i = 1:length(I_concavs_SA_node_sorted)
    if I_concavs_SA_node_dev_avg_pls_2SD(i) == min(I_concavs_SA_node_dev_avg_pls_2SD)
        I_concavs_SA_node_dev_avg_pls_2SD_index = i;
    end
end
I_concavs_SA_node_dev_avg_min_3SD = concav(length(I_concavs_SA_node_sorted),1);
for i = 1:length(I_concavs_SA_node_sorted)
    I_concavs_SA_node_dev_avg_min_3SD(i) = I_concavs_SA_node_sorted(i) - (I_concavs_SA_node_average-3*I_concavs_SA_node_SD);
end
I_concavs_SA_node_dev_avg_min_3SD = abs(I_concavs_SA_node_dev_avg_min_3SD);
for i = 1:length(I_concavs_SA_node_sorted)
    if I_concavs_SA_node_dev_avg_min_3SD(i) == min(I_concavs_SA_node_dev_avg_min_3SD)
        I_concavs_SA_node_dev_avg_min_3SD_index = i;
    end
end
I_concavs_SA_node_dev_avg_pls_3SD = concav(length(I_concavs_SA_node_sorted),1);
for i = 1:length(I_concavs_SA_node_sorted)
    I_concavs_SA_node_dev_avg_pls_3SD(i) = I_concavs_SA_node_sorted(i) - (I_concavs_SA_node_average+3*I_concavs_SA_node_SD);
end
I_concavs_SA_node_dev_avg_pls_3SD = abs(I_concavs_SA_node_dev_avg_pls_3SD);
for i = 1:length(I_concavs_SA_node_sorted)
    if I_concavs_SA_node_dev_avg_pls_3SD(i) == min(I_concavs_SA_node_dev_avg_pls_3SD)
        I_concavs_SA_node_dev_avg_pls_3SD_index = i;
    end
end
I_concavs_SA_node_tot_number_of_indices = length(I_concavs_SA_node_sorted);
%This is the approximation of the tot "area under curve." 
I_concavs_SA_node_quotient_within_1SD = (I_concavs_SA_node_dev_avg_pls_1SD_index - I_concavs_SA_nodedev_avg_min_1SD_index)/I_concavs_SA_node_tot_number_of_indices;
I_concavs_SA_node_quotient_within_2SD = (I_concavs_SA_node_dev_avg_pls_2SD_index - I_concavs_dev_avg_min_2SD_index)/I_concavs_SA_node_tot_number_of_indices;
I_concavs_SA_node_quotient_within_3SD = (I_concavs_SA_node_dev_avg_pls_3SD_index - I_concavs_SA_node_dev_avg_min_3SD_index)/I_concavs_SA_node_tot_number_of_indices;
if ((I_concavs_SA_node_quotient_within_1SD > 0.66 && I_concavs_SA_node_quotient_within_1SD < 0.70) && (I_concavs_SA_node_quotient_within_2SD > 0.93 && I_concavs_SA_node_quotient_within_2SD < 0.97) && (I_concavs_SA_node_quotient_within_3SD > 0.98 && I_concavs_SA_node_quotient_within_3SD < 1)) 
    isNormal(5,1) = 1;
end
if I_concavs_SA_node_quotient_within_1SD == 0
    isSpike(5,1) = 1;
end
I_concavs_crista_terminalis_sorted = sort(I_concavs_crista_terminalis);
I_concavs_crista_terminalis_dev_avg_min_1SD = concav(length(I_concavs_crista_terminalis_sorted),1);
for i = 1:length(I_concavs_crista_terminalis_sorted)
    I_concavs_crista_terminalis_dev_avg_min_1SD(i) = I_concavs_crista_terminalis_sorted(i) - (I_concavs_crista_terminalis_average-I_concavs_crista_terminalis_SD);
end
I_concavs_crista_terminalis_dev_avg_min_1SD = abs(I_concavs_crista_terminalis_dev_avg_min_1SD);
for i = 1:length(I_concavs_crista_terminalis_sorted)
    if I_concavs_crista_terminalis_dev_avg_min_1SD(i) == min(I_concavs_crista_terminalis_dev_avg_min_1SD)
        I_concavs_crista_terminalis_dev_avg_min_1SD_index = i;
    end
end
I_concavs_crista_terminalis_dev_avg_pls_1SD = concav(length(I_concavs_crista_terminalis_sorted),1);
for i = 1:length(I_concavs_crista_terminalis_sorted)
    I_concavs_crista_terminalis_dev_avg_pls_1SD(i) = I_concavs_crista_terminalis_sorted(i) - (I_concavs_crista_terminalis_average+I_concavs_crista_terminalis_SD);
end
I_concavs_crista_terminalis_dev_avg_pls_1SD = abs(I_concavs_crista_terminalis_dev_avg_pls_1SD);
for i = 1:length(I_concavs_crista_terminalis_sorted)
    if I_concavs_crista_terminalis_dev_avg_pls_1SD(i) == min(I_concavs_crista_terminalis_dev_avg_pls_1SD)
        I_concavs_crista_terminalis_dev_avg_pls_1SD_index = i;
    end
end
I_concavs_crista_terminalis_dev_avg_min_2SD = concav(length(I_concavs_crista_terminalis_sorted),1);
for i = 1:length(I_concavs_crista_terminalis_sorted)
    I_concavs_crista_terminalis_dev_avg_min_2SD(i) = I_concavs_crista_terminalis_sorted(i) - (I_concavs_crista_terminalis_average-2*I_concavs_crista_terminalis_SD);
end
I_concavs_crista_terminalis_dev_avg_min_2SD = abs(I_concavs_crista_terminalis_dev_avg_min_2SD);
for i = 1:length(I_concavs_crista_terminalis_sorted)
    if I_concavs_crista_terminalis_dev_avg_min_2SD(i) == min(I_concavs_crista_terminalis_dev_avg_min_2SD)
        I_concavs_dev_avg_min_2SD_index = i;
    end
end
I_concavs_crista_terminalis_dev_avg_pls_2SD = concav(length(I_concavs_crista_terminalis_sorted),1);
for i = 1:length(I_concavs_crista_terminalis_sorted)
    I_concavs_crista_terminalis_dev_avg_pls_2SD(i) = I_concavs_crista_terminalis_sorted(i) - (I_concavs_crista_terminalis_average+2*I_concavs_crista_terminalis_SD);
end
I_concavs_crista_terminalis_dev_avg_pls_2SD = abs(I_concavs_crista_terminalis_dev_avg_pls_2SD);
for i = 1:length(I_concavs_crista_terminalis_sorted)
    if I_concavs_crista_terminalis_dev_avg_pls_2SD(i) == min(I_concavs_crista_terminalis_dev_avg_pls_2SD)
        I_concavs_crista_terminalis_dev_avg_pls_2SD_index = i;
    end
end
I_concavs_crista_terminalis_dev_avg_min_3SD = concav(length(I_concavs_crista_terminalis_sorted),1);
for i = 1:length(I_concavs_crista_terminalis_sorted)
    I_concavs_crista_terminalis_dev_avg_min_3SD(i) = I_concavs_crista_terminalis_sorted(i) - (I_concavs_crista_terminalis_average-3*I_concavs_crista_terminalis_SD);
end
I_concavs_crista_terminalis_dev_avg_min_3SD = abs(I_concavs_crista_terminalis_dev_avg_min_3SD);
for i = 1:length(I_concavs_crista_terminalis_sorted)
    if I_concavs_crista_terminalis_dev_avg_min_3SD(i) == min(I_concavs_crista_terminalis_dev_avg_min_3SD)
        I_concavs_crista_terminalis_dev_avg_min_3SD_index = i;
    end
end
I_concavs_crista_terminalis_dev_avg_pls_3SD = concav(length(I_concavs_crista_terminalis_sorted),1);
for i = 1:length(I_concavs_crista_terminalis_sorted)
    I_concavs_crista_terminalis_dev_avg_pls_3SD(i) = I_concavs_crista_terminalis_sorted(i) - (I_concavs_crista_terminalis_average+3*I_concavs_crista_terminalis_SD);
end
I_concavs_crista_terminalis_dev_avg_pls_3SD = abs(I_concavs_crista_terminalis_dev_avg_pls_3SD);
for i = 1:length(I_concavs_crista_terminalis_sorted)
    if I_concavs_crista_terminalis_dev_avg_pls_3SD(i) == min(I_concavs_crista_terminalis_dev_avg_pls_3SD)
        I_concavs_crista_terminalis_dev_avg_pls_3SD_index = i;
    end
end
I_concavs_crista_terminalis_tot_number_of_indices = length(I_concavs_crista_terminalis_sorted);
%This is the approximation of the tot "area under curve." 
I_concavs_crista_terminalis_quotient_within_1SD = (I_concavs_crista_terminalis_dev_avg_pls_1SD_index - I_concavs_crista_terminalisdev_avg_min_1SD_index)/I_concavs_crista_terminalis_tot_number_of_indices;
I_concavs_crista_terminalis_quotient_within_2SD = (I_concavs_crista_terminalis_dev_avg_pls_2SD_index - I_concavs_dev_avg_min_2SD_index)/I_concavs_crista_terminalis_tot_number_of_indices;
I_concavs_crista_terminalis_quotient_within_3SD = (I_concavs_crista_terminalis_dev_avg_pls_3SD_index - I_concavs_crista_terminalis_dev_avg_min_3SD_index)/I_concavs_crista_terminalis_tot_number_of_indices;
if (I_concavs_crista_terminalis_quotient_within_1SD > 0.66 && I_concavs_crista_terminalis_quotient_within_1SD < 0.70) && (I_concavs_crista_terminalis_quotient_within_2SD > 0.93 && I_concavs_crista_terminalis_quotient_within_2SD < 0.97) && (I_concavs_crista_terminalis_quotient_within_3SD > 0.98 && I_concavs_crista_terminalis_quotient_within_3SD < 1)
    isNormal(5,2) = 1;
end
if I_concavs_crista_terminalis_quotient_within_1SD == 0
    isSpike(5,2) = 1;
end
I_concavs_tricuspid_annulus_sorted = sort(I_concavs_tricuspid_annulus);
I_concavs_tricuspid_annulus_dev_avg_min_1SD = concav(length(I_concavs_tricuspid_annulus_sorted),1);
for i = 1:length(I_concavs_tricuspid_annulus_sorted)
    I_concavs_tricuspid_annulus_dev_avg_min_1SD(i) = I_concavs_tricuspid_annulus_sorted(i) - (I_concavs_tricuspid_annulus_average-I_concavs_tricuspid_annulus_SD);
end
I_concavs_tricuspid_annulus_dev_avg_min_1SD = abs(I_concavs_tricuspid_annulus_dev_avg_min_1SD);
for i = 1:length(I_concavs_tricuspid_annulus_sorted)
    if I_concavs_tricuspid_annulus_dev_avg_min_1SD(i) == min(I_concavs_tricuspid_annulus_dev_avg_min_1SD)
        I_concavs_tricuspid_annulus_dev_avg_min_1SD_index = i;
    end
end
I_concavs_tricuspid_annulus_dev_avg_pls_1SD = concav(length(I_concavs_tricuspid_annulus_sorted),1);
for i = 1:length(I_concavs_tricuspid_annulus_sorted)
    I_concavs_tricuspid_annulus_dev_avg_pls_1SD(i) = I_concavs_tricuspid_annulus_sorted(i) - (I_concavs_tricuspid_annulus_average+I_concavs_tricuspid_annulus_SD);
end
I_concavs_tricuspid_annulus_dev_avg_pls_1SD = abs(I_concavs_tricuspid_annulus_dev_avg_pls_1SD);
for i = 1:length(I_concavs_tricuspid_annulus_sorted)
    if I_concavs_tricuspid_annulus_dev_avg_pls_1SD(i) == min(I_concavs_tricuspid_annulus_dev_avg_pls_1SD)
        I_concavs_tricuspid_annulus_dev_avg_pls_1SD_index = i;
    end
end
I_concavs_tricuspid_annulus_dev_avg_min_2SD = concav(length(I_concavs_tricuspid_annulus_sorted),1);
for i = 1:length(I_concavs_tricuspid_annulus_sorted)
    I_concavs_tricuspid_annulus_dev_avg_min_2SD(i) = I_concavs_tricuspid_annulus_sorted(i) - (I_concavs_tricuspid_annulus_average-2*I_concavs_tricuspid_annulus_SD);
end
I_concavs_tricuspid_annulus_dev_avg_min_2SD = abs(I_concavs_tricuspid_annulus_dev_avg_min_2SD);
for i = 1:length(I_concavs_tricuspid_annulus_sorted)
    if I_concavs_tricuspid_annulus_dev_avg_min_2SD(i) == min(I_concavs_tricuspid_annulus_dev_avg_min_2SD)
        I_concavs_dev_avg_min_2SD_index = i;
    end
end
I_concavs_tricuspid_annulus_dev_avg_pls_2SD = concav(length(I_concavs_tricuspid_annulus_sorted),1);
for i = 1:length(I_concavs_tricuspid_annulus_sorted)
    I_concavs_tricuspid_annulus_dev_avg_pls_2SD(i) = I_concavs_tricuspid_annulus_sorted(i) - (I_concavs_tricuspid_annulus_average+2*I_concavs_tricuspid_annulus_SD);
end
I_concavs_tricuspid_annulus_dev_avg_pls_2SD = abs(I_concavs_tricuspid_annulus_dev_avg_pls_2SD);
for i = 1:length(I_concavs_tricuspid_annulus_sorted)
    if I_concavs_tricuspid_annulus_dev_avg_pls_2SD(i) == min(I_concavs_tricuspid_annulus_dev_avg_pls_2SD)
        I_concavs_tricuspid_annulus_dev_avg_pls_2SD_index = i;
    end
end
I_concavs_tricuspid_annulus_dev_avg_min_3SD = concav(length(I_concavs_tricuspid_annulus_sorted),1);
for i = 1:length(I_concavs_tricuspid_annulus_sorted)
    I_concavs_tricuspid_annulus_dev_avg_min_3SD(i) = I_concavs_tricuspid_annulus_sorted(i) - (I_concavs_tricuspid_annulus_average-3*I_concavs_tricuspid_annulus_SD);
end
I_concavs_tricuspid_annulus_dev_avg_min_3SD = abs(I_concavs_tricuspid_annulus_dev_avg_min_3SD);
for i = 1:length(I_concavs_tricuspid_annulus_sorted)
    if I_concavs_tricuspid_annulus_dev_avg_min_3SD(i) == min(I_concavs_tricuspid_annulus_dev_avg_min_3SD)
        I_concavs_tricuspid_annulus_dev_avg_min_3SD_index = i;
    end
end
I_concavs_tricuspid_annulus_dev_avg_pls_3SD = concav(length(I_concavs_tricuspid_annulus_sorted),1);
for i = 1:length(I_concavs_tricuspid_annulus_sorted)
    I_concavs_tricuspid_annulus_dev_avg_pls_3SD(i) = I_concavs_tricuspid_annulus_sorted(i) - (I_concavs_tricuspid_annulus_average+3*I_concavs_tricuspid_annulus_SD);
end
I_concavs_tricuspid_annulus_dev_avg_pls_3SD = abs(I_concavs_tricuspid_annulus_dev_avg_pls_3SD);
for i = 1:length(I_concavs_tricuspid_annulus_sorted)
    if I_concavs_tricuspid_annulus_dev_avg_pls_3SD(i) == min(I_concavs_tricuspid_annulus_dev_avg_pls_3SD)
        I_concavs_tricuspid_annulus_dev_avg_pls_3SD_index = i;
    end
end
I_concavs_tricuspid_annulus_tot_number_of_indices = length(I_concavs_tricuspid_annulus_sorted);
%This is the approximation of the tot "area under curve." 
I_concavs_tricuspid_annulus_quotient_within_1SD = (I_concavs_tricuspid_annulus_dev_avg_pls_1SD_index - I_concavs_tricuspid_annulusdev_avg_min_1SD_index)/I_concavs_tricuspid_annulus_tot_number_of_indices;
I_concavs_tricuspid_annulus_quotient_within_2SD = (I_concavs_tricuspid_annulus_dev_avg_pls_2SD_index - I_concavs_dev_avg_min_2SD_index)/I_concavs_tricuspid_annulus_tot_number_of_indices;
I_concavs_tricuspid_annulus_quotient_within_3SD = (I_concavs_tricuspid_annulus_dev_avg_pls_3SD_index - I_concavs_tricuspid_annulus_dev_avg_min_3SD_index)/I_concavs_tricuspid_annulus_tot_number_of_indices;
if (I_concavs_tricuspid_annulus_quotient_within_1SD > 0.66 && I_concavs_tricuspid_annulus_quotient_within_1SD < 0.70) && (I_concavs_tricuspid_annulus_quotient_within_2SD > 0.93 && I_concavs_tricuspid_annulus_quotient_within_2SD < 0.97) && (I_concavs_tricuspid_annulus_quotient_within_3SD > 0.98 && I_concavs_tricuspid_annulus_quotient_within_3SD < 1)
    isNormal(5,3) = 1;
end
if I_concavs_tricuspid_annulus_quotient_within_1SD == 0
    isSpike(5,3) = 1;
end
I_concavs_coronary_sinus_sorted = sort(I_concavs_coronary_sinus);
I_concavs_coronary_sinus_dev_avg_min_1SD = concav(length(I_concavs_coronary_sinus_sorted),1);
for i = 1:length(I_concavs_coronary_sinus_sorted)
    I_concavs_coronary_sinus_dev_avg_min_1SD(i) = I_concavs_coronary_sinus_sorted(i) - (I_concavs_coronary_sinus_average-I_concavs_coronary_sinus_SD);
end
I_concavs_coronary_sinus_dev_avg_min_1SD = abs(I_concavs_coronary_sinus_dev_avg_min_1SD);
for i = 1:length(I_concavs_coronary_sinus_sorted)
    if I_concavs_coronary_sinus_dev_avg_min_1SD(i) == min(I_concavs_coronary_sinus_dev_avg_min_1SD)
        I_concavs_coronary_sinus_dev_avg_min_1SD_index = i;
    end
end
I_concavs_coronary_sinus_dev_avg_pls_1SD = concav(length(I_concavs_coronary_sinus_sorted),1);
for i = 1:length(I_concavs_coronary_sinus_sorted)
    I_concavs_coronary_sinus_dev_avg_pls_1SD(i) = I_concavs_coronary_sinus_sorted(i) - (I_concavs_coronary_sinus_average+I_concavs_coronary_sinus_SD);
end
I_concavs_coronary_sinus_dev_avg_pls_1SD = abs(I_concavs_coronary_sinus_dev_avg_pls_1SD);
for i = 1:length(I_concavs_coronary_sinus_sorted)
    if I_concavs_coronary_sinus_dev_avg_pls_1SD(i) == min(I_concavs_coronary_sinus_dev_avg_pls_1SD)
        I_concavs_coronary_sinus_dev_avg_pls_1SD_index = i;
    end
end
I_concavs_coronary_sinus_dev_avg_min_2SD = concav(length(I_concavs_coronary_sinus_sorted),1);
for i = 1:length(I_concavs_coronary_sinus_sorted)
    I_concavs_coronary_sinus_dev_avg_min_2SD(i) = I_concavs_coronary_sinus_sorted(i) - (I_concavs_coronary_sinus_average-2*I_concavs_coronary_sinus_SD);
end
I_concavs_coronary_sinus_dev_avg_min_2SD = abs(I_concavs_coronary_sinus_dev_avg_min_2SD);
for i = 1:length(I_concavs_coronary_sinus_sorted)
    if I_concavs_coronary_sinus_dev_avg_min_2SD(i) == min(I_concavs_coronary_sinus_dev_avg_min_2SD)
        I_concavs_dev_avg_min_2SD_index = i;
    end
end
I_concavs_coronary_sinus_dev_avg_pls_2SD = concav(length(I_concavs_coronary_sinus_sorted),1);
for i = 1:length(I_concavs_coronary_sinus_sorted)
    I_concavs_coronary_sinus_dev_avg_pls_2SD(i) = I_concavs_coronary_sinus_sorted(i) - (I_concavs_coronary_sinus_average+2*I_concavs_coronary_sinus_SD);
end
I_concavs_coronary_sinus_dev_avg_pls_2SD = abs(I_concavs_coronary_sinus_dev_avg_pls_2SD);
for i = 1:length(I_concavs_coronary_sinus_sorted)
    if I_concavs_coronary_sinus_dev_avg_pls_2SD(i) == min(I_concavs_coronary_sinus_dev_avg_pls_2SD)
        I_concavs_coronary_sinus_dev_avg_pls_2SD_index = i;
    end
end
I_concavs_coronary_sinus_dev_avg_min_3SD = concav(length(I_concavs_coronary_sinus_sorted),1);
for i = 1:length(I_concavs_coronary_sinus_sorted)
    I_concavs_coronary_sinus_dev_avg_min_3SD(i) = I_concavs_coronary_sinus_sorted(i) - (I_concavs_coronary_sinus_average-3*I_concavs_coronary_sinus_SD);
end
I_concavs_coronary_sinus_dev_avg_min_3SD = abs(I_concavs_coronary_sinus_dev_avg_min_3SD);
for i = 1:length(I_concavs_coronary_sinus_sorted)
    if I_concavs_coronary_sinus_dev_avg_min_3SD(i) == min(I_concavs_coronary_sinus_dev_avg_min_3SD)
        I_concavs_coronary_sinus_dev_avg_min_3SD_index = i;
    end
end
I_concavs_coronary_sinus_dev_avg_pls_3SD = concav(length(I_concavs_coronary_sinus_sorted),1);
for i = 1:length(I_concavs_coronary_sinus_sorted)
    I_concavs_coronary_sinus_dev_avg_pls_3SD(i) = I_concavs_coronary_sinus_sorted(i) - (I_concavs_coronary_sinus_average+3*I_concavs_coronary_sinus_SD);
end
I_concavs_coronary_sinus_dev_avg_pls_3SD = abs(I_concavs_coronary_sinus_dev_avg_pls_3SD);
for i = 1:length(I_concavs_coronary_sinus_sorted)
    if I_concavs_coronary_sinus_dev_avg_pls_3SD(i) == min(I_concavs_coronary_sinus_dev_avg_pls_3SD)
        I_concavs_coronary_sinus_dev_avg_pls_3SD_index = i;
    end
end
I_concavs_coronary_sinus_tot_number_of_indices = length(I_concavs_coronary_sinus_sorted);
%This is the approximation of the tot "area under curve." 
I_concavs_coronary_sinus_quotient_within_1SD = (I_concavs_coronary_sinus_dev_avg_pls_1SD_index - I_concavs_coronary_sinusdev_avg_min_1SD_index)/I_concavs_coronary_sinus_tot_number_of_indices;
I_concavs_coronary_sinus_quotient_within_2SD = (I_concavs_coronary_sinus_dev_avg_pls_2SD_index - I_concavs_dev_avg_min_2SD_index)/I_concavs_coronary_sinus_tot_number_of_indices;
I_concavs_coronary_sinus_quotient_within_3SD = (I_concavs_coronary_sinus_dev_avg_pls_3SD_index - I_concavs_coronary_sinus_dev_avg_min_3SD_index)/I_concavs_coronary_sinus_tot_number_of_indices;
if (I_concavs_coronary_sinus_quotient_within_1SD > 0.66 && I_concavs_coronary_sinus_quotient_within_1SD < 0.70) && (I_concavs_coronary_sinus_quotient_within_2SD > 0.93 && I_concavs_coronary_sinus_quotient_within_2SD < 0.97) && (I_concavs_coronary_sinus_quotient_within_3SD > 0.98 && I_concavs_coronary_sinus_quotient_within_3SD < 1)
    isNormal(5,4) = 1;
end
if I_concavs_coronary_sinus_quotient_within_1SD == 0
    isSpike(5,4) = 1;
end
I_concavs_ostium_sorted = sort(I_concavs_ostium);
I_concavs_ostium_dev_avg_min_1SD = concav(length(I_concavs_ostium_sorted),1);
for i = 1:length(I_concavs_ostium_sorted)
    I_concavs_ostium_dev_avg_min_1SD(i) = I_concavs_ostium_sorted(i) - (I_concavs_ostium_average-I_concavs_ostium_SD);
end
I_concavs_ostium_dev_avg_min_1SD = abs(I_concavs_ostium_dev_avg_min_1SD);
for i = 1:length(I_concavs_ostium_sorted)
    if I_concavs_ostium_dev_avg_min_1SD(i) == min(I_concavs_ostium_dev_avg_min_1SD)
        I_concavs_ostium_dev_avg_min_1SD_index = i;
    end
end
I_concavs_ostium_dev_avg_pls_1SD = concav(length(I_concavs_ostium_sorted),1);
for i = 1:length(I_concavs_ostium_sorted)
    I_concavs_ostium_dev_avg_pls_1SD(i) = I_concavs_ostium_sorted(i) - (I_concavs_ostium_average+I_concavs_ostium_SD);
end
I_concavs_ostium_dev_avg_pls_1SD = abs(I_concavs_ostium_dev_avg_pls_1SD);
for i = 1:length(I_concavs_ostium_sorted)
    if I_concavs_ostium_dev_avg_pls_1SD(i) == min(I_concavs_ostium_dev_avg_pls_1SD)
        I_concavs_ostium_dev_avg_pls_1SD_index = i;
    end
end
I_concavs_ostium_dev_avg_min_2SD = concav(length(I_concavs_ostium_sorted),1);
for i = 1:length(I_concavs_ostium_sorted)
    I_concavs_ostium_dev_avg_min_2SD(i) = I_concavs_ostium_sorted(i) - (I_concavs_ostium_average-2*I_concavs_ostium_SD);
end
I_concavs_ostium_dev_avg_min_2SD = abs(I_concavs_ostium_dev_avg_min_2SD);
for i = 1:length(I_concavs_ostium_sorted)
    if I_concavs_ostium_dev_avg_min_2SD(i) == min(I_concavs_ostium_dev_avg_min_2SD)
        I_concavs_dev_avg_min_2SD_index = i;
    end
end
I_concavs_ostium_dev_avg_pls_2SD = concav(length(I_concavs_ostium_sorted),1);
for i = 1:length(I_concavs_ostium_sorted)
    I_concavs_ostium_dev_avg_pls_2SD(i) = I_concavs_ostium_sorted(i) - (I_concavs_ostium_average+2*I_concavs_ostium_SD);
end
I_concavs_ostium_dev_avg_pls_2SD = abs(I_concavs_ostium_dev_avg_pls_2SD);
for i = 1:length(I_concavs_ostium_sorted)
    if I_concavs_ostium_dev_avg_pls_2SD(i) == min(I_concavs_ostium_dev_avg_pls_2SD)
        I_concavs_ostium_dev_avg_pls_2SD_index = i;
    end
end
I_concavs_ostium_dev_avg_min_3SD = concav(length(I_concavs_ostium_sorted),1);
for i = 1:length(I_concavs_ostium_sorted)
    I_concavs_ostium_dev_avg_min_3SD(i) = I_concavs_ostium_sorted(i) - (I_concavs_ostium_average-3*I_concavs_ostium_SD);
end
I_concavs_ostium_dev_avg_min_3SD = abs(I_concavs_ostium_dev_avg_min_3SD);
for i = 1:length(I_concavs_ostium_sorted)
    if I_concavs_ostium_dev_avg_min_3SD(i) == min(I_concavs_ostium_dev_avg_min_3SD)
        I_concavs_ostium_dev_avg_min_3SD_index = i;
    end
end
I_concavs_ostium_dev_avg_pls_3SD = concav(length(I_concavs_ostium_sorted),1);
for i = 1:length(I_concavs_ostium_sorted)
    I_concavs_ostium_dev_avg_pls_3SD(i) = I_concavs_ostium_sorted(i) - (I_concavs_ostium_average+3*I_concavs_ostium_SD);
end
I_concavs_ostium_dev_avg_pls_3SD = abs(I_concavs_ostium_dev_avg_pls_3SD);
for i = 1:length(I_concavs_ostium_sorted)
    if I_concavs_ostium_dev_avg_pls_3SD(i) == min(I_concavs_ostium_dev_avg_pls_3SD)
        I_concavs_ostium_dev_avg_pls_3SD_index = i;
    end
end
I_concavs_ostium_tot_number_of_indices = length(I_concavs_ostium_sorted);
%This is the approximation of the tot "area under curve." 
I_concavs_ostium_quotient_within_1SD = (I_concavs_ostium_dev_avg_pls_1SD_index - I_concavs_ostiumdev_avg_min_1SD_index)/I_concavs_ostium_tot_number_of_indices;
I_concavs_ostium_quotient_within_2SD = (I_concavs_ostium_dev_avg_pls_2SD_index - I_concavs_dev_avg_min_2SD_index)/I_concavs_ostium_tot_number_of_indices;
I_concavs_ostium_quotient_within_3SD = (I_concavs_ostium_dev_avg_pls_3SD_index - I_concavs_ostium_dev_avg_min_3SD_index)/I_concavs_ostium_tot_number_of_indices;
if (I_concavs_ostium_quotient_within_1SD > 0.66 && I_concavs_ostium_quotient_within_1SD < 0.70) && (I_concavs_ostium_quotient_within_2SD > 0.93 && I_concavs_ostium_quotient_within_2SD < 0.97) && (I_concavs_ostium_quotient_within_3SD > 0.98 && I_concavs_ostium_quotient_within_3SD < 1)
    isNormal(5,5) = 1;
end
if I_concavs_ostium_quotient_within_1SD == 0
    isSpike(5,5) = 1;
end
I_concavs_perinodal_sorted = sort(I_concavs_perinodal);
I_concavs_perinodal_dev_avg_min_1SD = concav(length(I_concavs_perinodal_sorted),1);
for i = 1:length(I_concavs_perinodal_sorted)
    I_concavs_perinodal_dev_avg_min_1SD(i) = I_concavs_perinodal_sorted(i) - (I_concavs_perinodal_average-I_concavs_perinodal_SD);
end
I_concavs_perinodal_dev_avg_min_1SD = abs(I_concavs_perinodal_dev_avg_min_1SD);
for i = 1:length(I_concavs_perinodal_sorted)
    if I_concavs_perinodal_dev_avg_min_1SD(i) == min(I_concavs_perinodal_dev_avg_min_1SD)
        I_concavs_perinodal_dev_avg_min_1SD_index = i;
    end
end
I_concavs_perinodal_dev_avg_pls_1SD = concav(length(I_concavs_perinodal_sorted),1);
for i = 1:length(I_concavs_perinodal_sorted)
    I_concavs_perinodal_dev_avg_pls_1SD(i) = I_concavs_perinodal_sorted(i) - (I_concavs_perinodal_average+I_concavs_perinodal_SD);
end
I_concavs_perinodal_dev_avg_pls_1SD = abs(I_concavs_perinodal_dev_avg_pls_1SD);
for i = 1:length(I_concavs_perinodal_sorted)
    if I_concavs_perinodal_dev_avg_pls_1SD(i) == min(I_concavs_perinodal_dev_avg_pls_1SD)
        I_concavs_perinodal_dev_avg_pls_1SD_index = i;
    end
end
I_concavs_perinodal_dev_avg_min_2SD = concav(length(I_concavs_perinodal_sorted),1);
for i = 1:length(I_concavs_perinodal_sorted)
    I_concavs_perinodal_dev_avg_min_2SD(i) = I_concavs_perinodal_sorted(i) - (I_concavs_perinodal_average-2*I_concavs_perinodal_SD);
end
I_concavs_perinodal_dev_avg_min_2SD = abs(I_concavs_perinodal_dev_avg_min_2SD);
for i = 1:length(I_concavs_perinodal_sorted)
    if I_concavs_perinodal_dev_avg_min_2SD(i) == min(I_concavs_perinodal_dev_avg_min_2SD)
        I_concavs_dev_avg_min_2SD_index = i;
    end
end
I_concavs_perinodal_dev_avg_pls_2SD = concav(length(I_concavs_perinodal_sorted),1);
for i = 1:length(I_concavs_perinodal_sorted)
    I_concavs_perinodal_dev_avg_pls_2SD(i) = I_concavs_perinodal_sorted(i) - (I_concavs_perinodal_average+2*I_concavs_perinodal_SD);
end
I_concavs_perinodal_dev_avg_pls_2SD = abs(I_concavs_perinodal_dev_avg_pls_2SD);
for i = 1:length(I_concavs_perinodal_sorted)
    if I_concavs_perinodal_dev_avg_pls_2SD(i) == min(I_concavs_perinodal_dev_avg_pls_2SD)
        I_concavs_perinodal_dev_avg_pls_2SD_index = i;
    end
end
I_concavs_perinodal_dev_avg_min_3SD = concav(length(I_concavs_perinodal_sorted),1);
for i = 1:length(I_concavs_perinodal_sorted)
    I_concavs_perinodal_dev_avg_min_3SD(i) = I_concavs_perinodal_sorted(i) - (I_concavs_perinodal_average-3*I_concavs_perinodal_SD);
end
I_concavs_perinodal_dev_avg_min_3SD = abs(I_concavs_perinodal_dev_avg_min_3SD);
for i = 1:length(I_concavs_perinodal_sorted)
    if I_concavs_perinodal_dev_avg_min_3SD(i) == min(I_concavs_perinodal_dev_avg_min_3SD)
        I_concavs_perinodal_dev_avg_min_3SD_index = i;
    end
end
I_concavs_perinodal_dev_avg_pls_3SD = concav(length(I_concavs_perinodal_sorted),1);
for i = 1:length(I_concavs_perinodal_sorted)
    I_concavs_perinodal_dev_avg_pls_3SD(i) = I_concavs_perinodal_sorted(i) - (I_concavs_perinodal_average+3*I_concavs_perinodal_SD);
end
I_concavs_perinodal_dev_avg_pls_3SD = abs(I_concavs_perinodal_dev_avg_pls_3SD);
for i = 1:length(I_concavs_perinodal_sorted)
    if I_concavs_perinodal_dev_avg_pls_3SD(i) == min(I_concavs_perinodal_dev_avg_pls_3SD)
        I_concavs_perinodal_dev_avg_pls_3SD_index = i;
    end
end
I_concavs_perinodal_tot_number_of_indices = length(I_concavs_perinodal_sorted);
%This is the approximation of the tot "area under curve." 
I_concavs_perinodal_quotient_within_1SD = (I_concavs_perinodal_dev_avg_pls_1SD_index - I_concavs_perinodaldev_avg_min_1SD_index)/I_concavs_perinodal_tot_number_of_indices;
I_concavs_perinodal_quotient_within_2SD = (I_concavs_perinodal_dev_avg_pls_2SD_index - I_concavs_dev_avg_min_2SD_index)/I_concavs_perinodal_tot_number_of_indices;
I_concavs_perinodal_quotient_within_3SD = (I_concavs_perinodal_dev_avg_pls_3SD_index - I_concavs_perinodal_dev_avg_min_3SD_index)/I_concavs_perinodal_tot_number_of_indices;
if (I_concavs_perinodal_quotient_within_1SD > 0.66 && I_concavs_perinodal_quotient_within_1SD < 0.70) && (I_concavs_perinodal_quotient_within_2SD > 0.93 && I_concavs_perinodal_quotient_within_2SD < 0.97) && (I_concavs_perinodal_quotient_within_3SD > 0.98 && I_concavs_perinodal_quotient_within_3SD < 1)
    isNormal(5,6) = 1;
end
if I_concavs_perinodal_quotient_within_1SD == 0
    isSpike(5,6) = 1;
end
I_concavs_right_septum_sorted = sort(I_concavs_right_septum);
I_concavs_right_septum_dev_avg_min_1SD = concav(length(I_concavs_right_septum_sorted),1);
for i = 1:length(I_concavs_right_septum_sorted)
    I_concavs_right_septum_dev_avg_min_1SD(i) = I_concavs_right_septum_sorted(i) - (I_concavs_right_septum_average-I_concavs_right_septum_SD);
end
I_concavs_right_septum_dev_avg_min_1SD = abs(I_concavs_right_septum_dev_avg_min_1SD);
for i = 1:length(I_concavs_right_septum_sorted)
    if I_concavs_right_septum_dev_avg_min_1SD(i) == min(I_concavs_right_septum_dev_avg_min_1SD)
        I_concavs_right_septum_dev_avg_min_1SD_index = i;
    end
end
I_concavs_right_septum_dev_avg_pls_1SD = concav(length(I_concavs_right_septum_sorted),1);
for i = 1:length(I_concavs_right_septum_sorted)
    I_concavs_right_septum_dev_avg_pls_1SD(i) = I_concavs_right_septum_sorted(i) - (I_concavs_right_septum_average+I_concavs_right_septum_SD);
end
I_concavs_right_septum_dev_avg_pls_1SD = abs(I_concavs_right_septum_dev_avg_pls_1SD);
for i = 1:length(I_concavs_right_septum_sorted)
    if I_concavs_right_septum_dev_avg_pls_1SD(i) == min(I_concavs_right_septum_dev_avg_pls_1SD)
        I_concavs_right_septum_dev_avg_pls_1SD_index = i;
    end
end
I_concavs_right_septum_dev_avg_min_2SD = concav(length(I_concavs_right_septum_sorted),1);
for i = 1:length(I_concavs_right_septum_sorted)
    I_concavs_right_septum_dev_avg_min_2SD(i) = I_concavs_right_septum_sorted(i) - (I_concavs_right_septum_average-2*I_concavs_right_septum_SD);
end
I_concavs_right_septum_dev_avg_min_2SD = abs(I_concavs_right_septum_dev_avg_min_2SD);
for i = 1:length(I_concavs_right_septum_sorted)
    if I_concavs_right_septum_dev_avg_min_2SD(i) == min(I_concavs_right_septum_dev_avg_min_2SD)
        I_concavs_dev_avg_min_2SD_index = i;
    end
end
I_concavs_right_septum_dev_avg_pls_2SD = concav(length(I_concavs_right_septum_sorted),1);
for i = 1:length(I_concavs_right_septum_sorted)
    I_concavs_right_septum_dev_avg_pls_2SD(i) = I_concavs_right_septum_sorted(i) - (I_concavs_right_septum_average+2*I_concavs_right_septum_SD);
end
I_concavs_right_septum_dev_avg_pls_2SD = abs(I_concavs_right_septum_dev_avg_pls_2SD);
for i = 1:length(I_concavs_right_septum_sorted)
    if I_concavs_right_septum_dev_avg_pls_2SD(i) == min(I_concavs_right_septum_dev_avg_pls_2SD)
        I_concavs_right_septum_dev_avg_pls_2SD_index = i;
    end
end
I_concavs_right_septum_dev_avg_min_3SD = concav(length(I_concavs_right_septum_sorted),1);
for i = 1:length(I_concavs_right_septum_sorted)
    I_concavs_right_septum_dev_avg_min_3SD(i) = I_concavs_right_septum_sorted(i) - (I_concavs_right_septum_average-3*I_concavs_right_septum_SD);
end
I_concavs_right_septum_dev_avg_min_3SD = abs(I_concavs_right_septum_dev_avg_min_3SD);
for i = 1:length(I_concavs_right_septum_sorted)
    if I_concavs_right_septum_dev_avg_min_3SD(i) == min(I_concavs_right_septum_dev_avg_min_3SD)
        I_concavs_right_septum_dev_avg_min_3SD_index = i;
    end
end
I_concavs_right_septum_dev_avg_pls_3SD = concav(length(I_concavs_right_septum_sorted),1);
for i = 1:length(I_concavs_right_septum_sorted)
    I_concavs_right_septum_dev_avg_pls_3SD(i) = I_concavs_right_septum_sorted(i) - (I_concavs_right_septum_average+3*I_concavs_right_septum_SD);
end
I_concavs_right_septum_dev_avg_pls_3SD = abs(I_concavs_right_septum_dev_avg_pls_3SD);
for i = 1:length(I_concavs_right_septum_sorted)
    if I_concavs_right_septum_dev_avg_pls_3SD(i) == min(I_concavs_right_septum_dev_avg_pls_3SD)
        I_concavs_right_septum_dev_avg_pls_3SD_index = i;
    end
end
I_concavs_right_septum_tot_number_of_indices = length(I_concavs_right_septum_sorted);
%This is the approximation of the tot "area under curve." 
I_concavs_right_septum_quotient_within_1SD = (I_concavs_right_septum_dev_avg_pls_1SD_index - I_concavs_right_septumdev_avg_min_1SD_index)/I_concavs_right_septum_tot_number_of_indices;
I_concavs_right_septum_quotient_within_2SD = (I_concavs_right_septum_dev_avg_pls_2SD_index - I_concavs_dev_avg_min_2SD_index)/I_concavs_right_septum_tot_number_of_indices;
I_concavs_right_septum_quotient_within_3SD = (I_concavs_right_septum_dev_avg_pls_3SD_index - I_concavs_right_septum_dev_avg_min_3SD_index)/I_concavs_right_septum_tot_number_of_indices;
if (I_concavs_right_septum_quotient_within_1SD > 0.66 && I_concavs_right_septum_quotient_within_1SD < 0.70) && (I_concavs_right_septum_quotient_within_2SD > 0.93 && I_concavs_right_septum_quotient_within_2SD < 0.97) && (I_concavs_right_septum_quotient_within_3SD > 0.98 && I_concavs_right_septum_quotient_within_3SD < 1)
    isNormal(5,7) = 1;
end
if I_concavs_right_septum_quotient_within_1SD == 0
    isSpike(5,7) = 1;
end
I_concavs_right_atrial_appendage_sorted = sort(I_concavs_right_atrial_appendage);
I_concavs_right_atrial_appendage_dev_avg_min_1SD = concav(length(I_concavs_right_atrial_appendage_sorted),1);
for i = 1:length(I_concavs_right_atrial_appendage_sorted)
    I_concavs_right_atrial_appendage_dev_avg_min_1SD(i) = I_concavs_right_atrial_appendage_sorted(i) - (I_concavs_right_atrial_appendage_average-I_concavs_right_atrial_appendage_SD);
end
I_concavs_right_atrial_appendage_dev_avg_min_1SD = abs(I_concavs_right_atrial_appendage_dev_avg_min_1SD);
for i = 1:length(I_concavs_right_atrial_appendage_sorted)
    if I_concavs_right_atrial_appendage_dev_avg_min_1SD(i) == min(I_concavs_right_atrial_appendage_dev_avg_min_1SD)
        I_concavs_right_atrial_appendage_dev_avg_min_1SD_index = i;
    end
end
I_concavs_right_atrial_appendage_dev_avg_pls_1SD = concav(length(I_concavs_right_atrial_appendage_sorted),1);
for i = 1:length(I_concavs_right_atrial_appendage_sorted)
    I_concavs_right_atrial_appendage_dev_avg_pls_1SD(i) = I_concavs_right_atrial_appendage_sorted(i) - (I_concavs_right_atrial_appendage_average+I_concavs_right_atrial_appendage_SD);
end
I_concavs_right_atrial_appendage_dev_avg_pls_1SD = abs(I_concavs_right_atrial_appendage_dev_avg_pls_1SD);
for i = 1:length(I_concavs_right_atrial_appendage_sorted)
    if I_concavs_right_atrial_appendage_dev_avg_pls_1SD(i) == min(I_concavs_right_atrial_appendage_dev_avg_pls_1SD)
        I_concavs_right_atrial_appendage_dev_avg_pls_1SD_index = i;
    end
end
I_concavs_right_atrial_appendage_dev_avg_min_2SD = concav(length(I_concavs_right_atrial_appendage_sorted),1);
for i = 1:length(I_concavs_right_atrial_appendage_sorted)
    I_concavs_right_atrial_appendage_dev_avg_min_2SD(i) = I_concavs_right_atrial_appendage_sorted(i) - (I_concavs_right_atrial_appendage_average-2*I_concavs_right_atrial_appendage_SD);
end
I_concavs_right_atrial_appendage_dev_avg_min_2SD = abs(I_concavs_right_atrial_appendage_dev_avg_min_2SD);
for i = 1:length(I_concavs_right_atrial_appendage_sorted)
    if I_concavs_right_atrial_appendage_dev_avg_min_2SD(i) == min(I_concavs_right_atrial_appendage_dev_avg_min_2SD)
        I_concavs_dev_avg_min_2SD_index = i;
    end
end
I_concavs_right_atrial_appendage_dev_avg_pls_2SD = concav(length(I_concavs_right_atrial_appendage_sorted),1);
for i = 1:length(I_concavs_right_atrial_appendage_sorted)
    I_concavs_right_atrial_appendage_dev_avg_pls_2SD(i) = I_concavs_right_atrial_appendage_sorted(i) - (I_concavs_right_atrial_appendage_average+2*I_concavs_right_atrial_appendage_SD);
end
I_concavs_right_atrial_appendage_dev_avg_pls_2SD = abs(I_concavs_right_atrial_appendage_dev_avg_pls_2SD);
for i = 1:length(I_concavs_right_atrial_appendage_sorted)
    if I_concavs_right_atrial_appendage_dev_avg_pls_2SD(i) == min(I_concavs_right_atrial_appendage_dev_avg_pls_2SD)
        I_concavs_right_atrial_appendage_dev_avg_pls_2SD_index = i;
    end
end
I_concavs_right_atrial_appendage_dev_avg_min_3SD = concav(length(I_concavs_right_atrial_appendage_sorted),1);
for i = 1:length(I_concavs_right_atrial_appendage_sorted)
    I_concavs_right_atrial_appendage_dev_avg_min_3SD(i) = I_concavs_right_atrial_appendage_sorted(i) - (I_concavs_right_atrial_appendage_average-3*I_concavs_right_atrial_appendage_SD);
end
I_concavs_right_atrial_appendage_dev_avg_min_3SD = abs(I_concavs_right_atrial_appendage_dev_avg_min_3SD);
for i = 1:length(I_concavs_right_atrial_appendage_sorted)
    if I_concavs_right_atrial_appendage_dev_avg_min_3SD(i) == min(I_concavs_right_atrial_appendage_dev_avg_min_3SD)
        I_concavs_right_atrial_appendage_dev_avg_min_3SD_index = i;
    end
end
I_concavs_right_atrial_appendage_dev_avg_pls_3SD = concav(length(I_concavs_right_atrial_appendage_sorted),1);
for i = 1:length(I_concavs_right_atrial_appendage_sorted)
    I_concavs_right_atrial_appendage_dev_avg_pls_3SD(i) = I_concavs_right_atrial_appendage_sorted(i) - (I_concavs_right_atrial_appendage_average+3*I_concavs_right_atrial_appendage_SD);
end
I_concavs_right_atrial_appendage_dev_avg_pls_3SD = abs(I_concavs_right_atrial_appendage_dev_avg_pls_3SD);
for i = 1:length(I_concavs_right_atrial_appendage_sorted)
    if I_concavs_right_atrial_appendage_dev_avg_pls_3SD(i) == min(I_concavs_right_atrial_appendage_dev_avg_pls_3SD)
        I_concavs_right_atrial_appendage_dev_avg_pls_3SD_index = i;
    end
end
I_concavs_right_atrial_appendage_tot_number_of_indices = length(I_concavs_right_atrial_appendage_sorted);
%This is the approximation of the tot "area under curve." 
I_concavs_right_atrial_appendage_quotient_within_1SD = (I_concavs_right_atrial_appendage_dev_avg_pls_1SD_index - I_concavs_right_atrial_appendagedev_avg_min_1SD_index)/I_concavs_right_atrial_appendage_tot_number_of_indices;
I_concavs_right_atrial_appendage_quotient_within_2SD = (I_concavs_right_atrial_appendage_dev_avg_pls_2SD_index - I_concavs_dev_avg_min_2SD_index)/I_concavs_right_atrial_appendage_tot_number_of_indices;
I_concavs_right_atrial_appendage_quotient_within_3SD = (I_concavs_right_atrial_appendage_dev_avg_pls_3SD_index - I_concavs_right_atrial_appendage_dev_avg_min_3SD_index)/I_concavs_right_atrial_appendage_tot_number_of_indices;
if (I_concavs_right_atrial_appendage_quotient_within_1SD > 0.66 && I_concavs_right_atrial_appendage_quotient_within_1SD < 0.70) && (I_concavs_right_atrial_appendage_quotient_within_2SD > 0.93 && I_concavs_right_atrial_appendage_quotient_within_2SD < 0.97) && (I_concavs_right_atrial_appendage_quotient_within_3SD > 0.98 && I_concavs_right_atrial_appendage_quotient_within_3SD < 1)
    isNormal(5,8) = 1;
end
if I_concavs_right_atrial_appendage_quotient_within_1SD == 0
    isSpike(5,8) = 1;
end
I_concavs_pulmonary_veins_sorted = sort(I_concavs_pulmonary_veins);
I_concavs_pulmonary_veins_dev_avg_min_1SD = concav(length(I_concavs_pulmonary_veins_sorted),1);
for i = 1:length(I_concavs_pulmonary_veins_sorted)
    I_concavs_pulmonary_veins_dev_avg_min_1SD(i) = I_concavs_pulmonary_veins_sorted(i) - (I_concavs_pulmonary_veins_average-I_concavs_pulmonary_veins_SD);
end
I_concavs_pulmonary_veins_dev_avg_min_1SD = abs(I_concavs_pulmonary_veins_dev_avg_min_1SD);
for i = 1:length(I_concavs_pulmonary_veins_sorted)
    if I_concavs_pulmonary_veins_dev_avg_min_1SD(i) == min(I_concavs_pulmonary_veins_dev_avg_min_1SD)
        I_concavs_pulmonary_veins_dev_avg_min_1SD_index = i;
    end
end
I_concavs_pulmonary_veins_dev_avg_pls_1SD = concav(length(I_concavs_pulmonary_veins_sorted),1);
for i = 1:length(I_concavs_pulmonary_veins_sorted)
    I_concavs_pulmonary_veins_dev_avg_pls_1SD(i) = I_concavs_pulmonary_veins_sorted(i) - (I_concavs_pulmonary_veins_average+I_concavs_pulmonary_veins_SD);
end
I_concavs_pulmonary_veins_dev_avg_pls_1SD = abs(I_concavs_pulmonary_veins_dev_avg_pls_1SD);
for i = 1:length(I_concavs_pulmonary_veins_sorted)
    if I_concavs_pulmonary_veins_dev_avg_pls_1SD(i) == min(I_concavs_pulmonary_veins_dev_avg_pls_1SD)
        I_concavs_pulmonary_veins_dev_avg_pls_1SD_index = i;
    end
end
I_concavs_pulmonary_veins_dev_avg_min_2SD = concav(length(I_concavs_pulmonary_veins_sorted),1);
for i = 1:length(I_concavs_pulmonary_veins_sorted)
    I_concavs_pulmonary_veins_dev_avg_min_2SD(i) = I_concavs_pulmonary_veins_sorted(i) - (I_concavs_pulmonary_veins_average-2*I_concavs_pulmonary_veins_SD);
end
I_concavs_pulmonary_veins_dev_avg_min_2SD = abs(I_concavs_pulmonary_veins_dev_avg_min_2SD);
for i = 1:length(I_concavs_pulmonary_veins_sorted)
    if I_concavs_pulmonary_veins_dev_avg_min_2SD(i) == min(I_concavs_pulmonary_veins_dev_avg_min_2SD)
        I_concavs_dev_avg_min_2SD_index = i;
    end
end
I_concavs_pulmonary_veins_dev_avg_pls_2SD = concav(length(I_concavs_pulmonary_veins_sorted),1);
for i = 1:length(I_concavs_pulmonary_veins_sorted)
    I_concavs_pulmonary_veins_dev_avg_pls_2SD(i) = I_concavs_pulmonary_veins_sorted(i) - (I_concavs_pulmonary_veins_average+2*I_concavs_pulmonary_veins_SD);
end
I_concavs_pulmonary_veins_dev_avg_pls_2SD = abs(I_concavs_pulmonary_veins_dev_avg_pls_2SD);
for i = 1:length(I_concavs_pulmonary_veins_sorted)
    if I_concavs_pulmonary_veins_dev_avg_pls_2SD(i) == min(I_concavs_pulmonary_veins_dev_avg_pls_2SD)
        I_concavs_pulmonary_veins_dev_avg_pls_2SD_index = i;
    end
end
I_concavs_pulmonary_veins_dev_avg_min_3SD = concav(length(I_concavs_pulmonary_veins_sorted),1);
for i = 1:length(I_concavs_pulmonary_veins_sorted)
    I_concavs_pulmonary_veins_dev_avg_min_3SD(i) = I_concavs_pulmonary_veins_sorted(i) - (I_concavs_pulmonary_veins_average-3*I_concavs_pulmonary_veins_SD);
end
I_concavs_pulmonary_veins_dev_avg_min_3SD = abs(I_concavs_pulmonary_veins_dev_avg_min_3SD);
for i = 1:length(I_concavs_pulmonary_veins_sorted)
    if I_concavs_pulmonary_veins_dev_avg_min_3SD(i) == min(I_concavs_pulmonary_veins_dev_avg_min_3SD)
        I_concavs_pulmonary_veins_dev_avg_min_3SD_index = i;
    end
end
I_concavs_pulmonary_veins_dev_avg_pls_3SD = concav(length(I_concavs_pulmonary_veins_sorted),1);
for i = 1:length(I_concavs_pulmonary_veins_sorted)
    I_concavs_pulmonary_veins_dev_avg_pls_3SD(i) = I_concavs_pulmonary_veins_sorted(i) - (I_concavs_pulmonary_veins_average+3*I_concavs_pulmonary_veins_SD);
end
I_concavs_pulmonary_veins_dev_avg_pls_3SD = abs(I_concavs_pulmonary_veins_dev_avg_pls_3SD);
for i = 1:length(I_concavs_pulmonary_veins_sorted)
    if I_concavs_pulmonary_veins_dev_avg_pls_3SD(i) == min(I_concavs_pulmonary_veins_dev_avg_pls_3SD)
        I_concavs_pulmonary_veins_dev_avg_pls_3SD_index = i;
    end
end
I_concavs_pulmonary_veins_tot_number_of_indices = length(I_concavs_pulmonary_veins_sorted);
%This is the approximation of the tot "area under curve." 
I_concavs_pulmonary_veins_quotient_within_1SD = (I_concavs_pulmonary_veins_dev_avg_pls_1SD_index - I_concavs_pulmonary_veinsdev_avg_min_1SD_index)/I_concavs_pulmonary_veins_tot_number_of_indices;
I_concavs_pulmonary_veins_quotient_within_2SD = (I_concavs_pulmonary_veins_dev_avg_pls_2SD_index - I_concavs_dev_avg_min_2SD_index)/I_concavs_pulmonary_veins_tot_number_of_indices;
I_concavs_pulmonary_veins_quotient_within_3SD = (I_concavs_pulmonary_veins_dev_avg_pls_3SD_index - I_concavs_pulmonary_veins_dev_avg_min_3SD_index)/I_concavs_pulmonary_veins_tot_number_of_indices;
if (I_concavs_pulmonary_veins_quotient_within_1SD > 0.66 && I_concavs_pulmonary_veins_quotient_within_1SD < 0.70) && (I_concavs_pulmonary_veins_quotient_within_2SD > 0.93 && I_concavs_pulmonary_veins_quotient_within_2SD < 0.97) && (I_concavs_pulmonary_veins_quotient_within_3SD > 0.98 && I_concavs_pulmonary_veins_quotient_within_3SD < 1)
    isNormal(5,9) = 1;
end
if I_concavs_pulmonary_veins_quotient_within_1SD == 0
    isSpike(5,9) = 1;
end
I_concavs_mitral_annulus_sorted = sort(I_concavs_mitral_annulus);
I_concavs_mitral_annulus_dev_avg_min_1SD = concav(length(I_concavs_mitral_annulus_sorted),1);
for i = 1:length(I_concavs_mitral_annulus_sorted)
    I_concavs_mitral_annulus_dev_avg_min_1SD(i) = I_concavs_mitral_annulus_sorted(i) - (I_concavs_mitral_annulus_average-I_concavs_mitral_annulus_SD);
end
I_concavs_mitral_annulus_dev_avg_min_1SD = abs(I_concavs_mitral_annulus_dev_avg_min_1SD);
for i = 1:length(I_concavs_mitral_annulus_sorted)
    if I_concavs_mitral_annulus_dev_avg_min_1SD(i) == min(I_concavs_mitral_annulus_dev_avg_min_1SD)
        I_concavs_mitral_annulus_dev_avg_min_1SD_index = i;
    end
end
I_concavs_mitral_annulus_dev_avg_pls_1SD = concav(length(I_concavs_mitral_annulus_sorted),1);
for i = 1:length(I_concavs_mitral_annulus_sorted)
    I_concavs_mitral_annulus_dev_avg_pls_1SD(i) = I_concavs_mitral_annulus_sorted(i) - (I_concavs_mitral_annulus_average+I_concavs_mitral_annulus_SD);
end
I_concavs_mitral_annulus_dev_avg_pls_1SD = abs(I_concavs_mitral_annulus_dev_avg_pls_1SD);
for i = 1:length(I_concavs_mitral_annulus_sorted)
    if I_concavs_mitral_annulus_dev_avg_pls_1SD(i) == min(I_concavs_mitral_annulus_dev_avg_pls_1SD)
        I_concavs_mitral_annulus_dev_avg_pls_1SD_index = i;
    end
end
I_concavs_mitral_annulus_dev_avg_min_2SD = concav(length(I_concavs_mitral_annulus_sorted),1);
for i = 1:length(I_concavs_mitral_annulus_sorted)
    I_concavs_mitral_annulus_dev_avg_min_2SD(i) = I_concavs_mitral_annulus_sorted(i) - (I_concavs_mitral_annulus_average-2*I_concavs_mitral_annulus_SD);
end
I_concavs_mitral_annulus_dev_avg_min_2SD = abs(I_concavs_mitral_annulus_dev_avg_min_2SD);
for i = 1:length(I_concavs_mitral_annulus_sorted)
    if I_concavs_mitral_annulus_dev_avg_min_2SD(i) == min(I_concavs_mitral_annulus_dev_avg_min_2SD)
        I_concavs_dev_avg_min_2SD_index = i;
    end
end
I_concavs_mitral_annulus_dev_avg_pls_2SD = concav(length(I_concavs_mitral_annulus_sorted),1);
for i = 1:length(I_concavs_mitral_annulus_sorted)
    I_concavs_mitral_annulus_dev_avg_pls_2SD(i) = I_concavs_mitral_annulus_sorted(i) - (I_concavs_mitral_annulus_average+2*I_concavs_mitral_annulus_SD);
end
I_concavs_mitral_annulus_dev_avg_pls_2SD = abs(I_concavs_mitral_annulus_dev_avg_pls_2SD);
for i = 1:length(I_concavs_mitral_annulus_sorted)
    if I_concavs_mitral_annulus_dev_avg_pls_2SD(i) == min(I_concavs_mitral_annulus_dev_avg_pls_2SD)
        I_concavs_mitral_annulus_dev_avg_pls_2SD_index = i;
    end
end
I_concavs_mitral_annulus_dev_avg_min_3SD = concav(length(I_concavs_mitral_annulus_sorted),1);
for i = 1:length(I_concavs_mitral_annulus_sorted)
    I_concavs_mitral_annulus_dev_avg_min_3SD(i) = I_concavs_mitral_annulus_sorted(i) - (I_concavs_mitral_annulus_average-3*I_concavs_mitral_annulus_SD);
end
I_concavs_mitral_annulus_dev_avg_min_3SD = abs(I_concavs_mitral_annulus_dev_avg_min_3SD);
for i = 1:length(I_concavs_mitral_annulus_sorted)
    if I_concavs_mitral_annulus_dev_avg_min_3SD(i) == min(I_concavs_mitral_annulus_dev_avg_min_3SD)
        I_concavs_mitral_annulus_dev_avg_min_3SD_index = i;
    end
end
I_concavs_mitral_annulus_dev_avg_pls_3SD = concav(length(I_concavs_mitral_annulus_sorted),1);
for i = 1:length(I_concavs_mitral_annulus_sorted)
    I_concavs_mitral_annulus_dev_avg_pls_3SD(i) = I_concavs_mitral_annulus_sorted(i) - (I_concavs_mitral_annulus_average+3*I_concavs_mitral_annulus_SD);
end
I_concavs_mitral_annulus_dev_avg_pls_3SD = abs(I_concavs_mitral_annulus_dev_avg_pls_3SD);
for i = 1:length(I_concavs_mitral_annulus_sorted)
    if I_concavs_mitral_annulus_dev_avg_pls_3SD(i) == min(I_concavs_mitral_annulus_dev_avg_pls_3SD)
        I_concavs_mitral_annulus_dev_avg_pls_3SD_index = i;
    end
end
I_concavs_mitral_annulus_tot_number_of_indices = length(I_concavs_mitral_annulus_sorted);
%This is the approximation of the tot "area under curve." 
I_concavs_mitral_annulus_quotient_within_1SD = (I_concavs_mitral_annulus_dev_avg_pls_1SD_index - I_concavs_mitral_annulusdev_avg_min_1SD_index)/I_concavs_mitral_annulus_tot_number_of_indices;
I_concavs_mitral_annulus_quotient_within_2SD = (I_concavs_mitral_annulus_dev_avg_pls_2SD_index - I_concavs_dev_avg_min_2SD_index)/I_concavs_mitral_annulus_tot_number_of_indices;
I_concavs_mitral_annulus_quotient_within_3SD = (I_concavs_mitral_annulus_dev_avg_pls_3SD_index - I_concavs_mitral_annulus_dev_avg_min_3SD_index)/I_concavs_mitral_annulus_tot_number_of_indices;
if (I_concavs_mitral_annulus_quotient_within_1SD > 0.66 && I_concavs_mitral_annulus_quotient_within_1SD < 0.70) && (I_concavs_mitral_annulus_quotient_within_2SD > 0.93 && I_concavs_mitral_annulus_quotient_within_2SD < 0.97) && (I_concavs_mitral_annulus_quotient_within_3SD > 0.98 && I_concavs_mitral_annulus_quotient_within_3SD < 1)
    isNormal(5,10) = 1;
end
if I_concavs_mitral_annulus_quotient_within_1SD == 0
    isSpike(5,10) = 1;
end
I_concavs_CS_body_sorted = sort(I_concavs_CS_body);
I_concavs_CS_body_dev_avg_min_1SD = concav(length(I_concavs_CS_body_sorted),1);
for i = 1:length(I_concavs_CS_body_sorted)
    I_concavs_CS_body_dev_avg_min_1SD(i) = I_concavs_CS_body_sorted(i) - (I_concavs_CS_body_average-I_concavs_CS_body_SD);
end
I_concavs_CS_body_dev_avg_min_1SD = abs(I_concavs_CS_body_dev_avg_min_1SD);
for i = 1:length(I_concavs_CS_body_sorted)
    if I_concavs_CS_body_dev_avg_min_1SD(i) == min(I_concavs_CS_body_dev_avg_min_1SD)
        I_concavs_CS_body_dev_avg_min_1SD_index = i;
    end
end
I_concavs_CS_body_dev_avg_pls_1SD = concav(length(I_concavs_CS_body_sorted),1);
for i = 1:length(I_concavs_CS_body_sorted)
    I_concavs_CS_body_dev_avg_pls_1SD(i) = I_concavs_CS_body_sorted(i) - (I_concavs_CS_body_average+I_concavs_CS_body_SD);
end
I_concavs_CS_body_dev_avg_pls_1SD = abs(I_concavs_CS_body_dev_avg_pls_1SD);
for i = 1:length(I_concavs_CS_body_sorted)
    if I_concavs_CS_body_dev_avg_pls_1SD(i) == min(I_concavs_CS_body_dev_avg_pls_1SD)
        I_concavs_CS_body_dev_avg_pls_1SD_index = i;
    end
end
I_concavs_CS_body_dev_avg_min_2SD = concav(length(I_concavs_CS_body_sorted),1);
for i = 1:length(I_concavs_CS_body_sorted)
    I_concavs_CS_body_dev_avg_min_2SD(i) = I_concavs_CS_body_sorted(i) - (I_concavs_CS_body_average-2*I_concavs_CS_body_SD);
end
I_concavs_CS_body_dev_avg_min_2SD = abs(I_concavs_CS_body_dev_avg_min_2SD);
for i = 1:length(I_concavs_CS_body_sorted)
    if I_concavs_CS_body_dev_avg_min_2SD(i) == min(I_concavs_CS_body_dev_avg_min_2SD)
        I_concavs_dev_avg_min_2SD_index = i;
    end
end
I_concavs_CS_body_dev_avg_pls_2SD = concav(length(I_concavs_CS_body_sorted),1);
for i = 1:length(I_concavs_CS_body_sorted)
    I_concavs_CS_body_dev_avg_pls_2SD(i) = I_concavs_CS_body_sorted(i) - (I_concavs_CS_body_average+2*I_concavs_CS_body_SD);
end
I_concavs_CS_body_dev_avg_pls_2SD = abs(I_concavs_CS_body_dev_avg_pls_2SD);
for i = 1:length(I_concavs_CS_body_sorted)
    if I_concavs_CS_body_dev_avg_pls_2SD(i) == min(I_concavs_CS_body_dev_avg_pls_2SD)
        I_concavs_CS_body_dev_avg_pls_2SD_index = i;
    end
end
I_concavs_CS_body_dev_avg_min_3SD = concav(length(I_concavs_CS_body_sorted),1);
for i = 1:length(I_concavs_CS_body_sorted)
    I_concavs_CS_body_dev_avg_min_3SD(i) = I_concavs_CS_body_sorted(i) - (I_concavs_CS_body_average-3*I_concavs_CS_body_SD);
end
I_concavs_CS_body_dev_avg_min_3SD = abs(I_concavs_CS_body_dev_avg_min_3SD);
for i = 1:length(I_concavs_CS_body_sorted)
    if I_concavs_CS_body_dev_avg_min_3SD(i) == min(I_concavs_CS_body_dev_avg_min_3SD)
        I_concavs_CS_body_dev_avg_min_3SD_index = i;
    end
end
I_concavs_CS_body_dev_avg_pls_3SD = concav(length(I_concavs_CS_body_sorted),1);
for i = 1:length(I_concavs_CS_body_sorted)
    I_concavs_CS_body_dev_avg_pls_3SD(i) = I_concavs_CS_body_sorted(i) - (I_concavs_CS_body_average+3*I_concavs_CS_body_SD);
end
I_concavs_CS_body_dev_avg_pls_3SD = abs(I_concavs_CS_body_dev_avg_pls_3SD);
for i = 1:length(I_concavs_CS_body_sorted)
    if I_concavs_CS_body_dev_avg_pls_3SD(i) == min(I_concavs_CS_body_dev_avg_pls_3SD)
        I_concavs_CS_body_dev_avg_pls_3SD_index = i;
    end
end
I_concavs_CS_body_tot_number_of_indices = length(I_concavs_CS_body_sorted);
%This is the approximation of the tot "area under curve." 
I_concavs_CS_body_quotient_within_1SD = (I_concavs_CS_body_dev_avg_pls_1SD_index - I_concavs_CS_bodydev_avg_min_1SD_index)/I_concavs_CS_body_tot_number_of_indices;
I_concavs_CS_body_quotient_within_2SD = (I_concavs_CS_body_dev_avg_pls_2SD_index - I_concavs_dev_avg_min_2SD_index)/I_concavs_CS_body_tot_number_of_indices;
I_concavs_CS_body_quotient_within_3SD = (I_concavs_CS_body_dev_avg_pls_3SD_index - I_concavs_CS_body_dev_avg_min_3SD_index)/I_concavs_CS_body_tot_number_of_indices;
if (I_concavs_CS_body_quotient_within_1SD > 0.66 && I_concavs_CS_body_quotient_within_1SD < 0.70) && (I_concavs_CS_body_quotient_within_2SD > 0.93 && I_concavs_CS_body_quotient_within_2SD < 0.97) && (I_concavs_CS_body_quotient_within_3SD > 0.98 && I_concavs_CS_body_quotient_within_3SD < 1)
    isNormal(5,11) = 1;
end
if I_concavs_CS_body_quotient_within_1SD == 0
    isSpike(5,11) = 1;
end
I_concavs_left_septum_sorted = sort(I_concavs_left_septum);
I_concavs_left_septum_dev_avg_min_1SD = concav(length(I_concavs_left_septum_sorted),1);
for i = 1:length(I_concavs_left_septum_sorted)
    I_concavs_left_septum_dev_avg_min_1SD(i) = I_concavs_left_septum_sorted(i) - (I_concavs_left_septum_average-I_concavs_left_septum_SD);
end
I_concavs_left_septum_dev_avg_min_1SD = abs(I_concavs_left_septum_dev_avg_min_1SD);
for i = 1:length(I_concavs_left_septum_sorted)
    if I_concavs_left_septum_dev_avg_min_1SD(i) == min(I_concavs_left_septum_dev_avg_min_1SD)
        I_concavs_left_septum_dev_avg_min_1SD_index = i;
    end
end
I_concavs_left_septum_dev_avg_pls_1SD = concav(length(I_concavs_left_septum_sorted),1);
for i = 1:length(I_concavs_left_septum_sorted)
    I_concavs_left_septum_dev_avg_pls_1SD(i) = I_concavs_left_septum_sorted(i) - (I_concavs_left_septum_average+I_concavs_left_septum_SD);
end
I_concavs_left_septum_dev_avg_pls_1SD = abs(I_concavs_left_septum_dev_avg_pls_1SD);
for i = 1:length(I_concavs_left_septum_sorted)
    if I_concavs_left_septum_dev_avg_pls_1SD(i) == min(I_concavs_left_septum_dev_avg_pls_1SD)
        I_concavs_left_septum_dev_avg_pls_1SD_index = i;
    end
end
I_concavs_left_septum_dev_avg_min_2SD = concav(length(I_concavs_left_septum_sorted),1);
for i = 1:length(I_concavs_left_septum_sorted)
    I_concavs_left_septum_dev_avg_min_2SD(i) = I_concavs_left_septum_sorted(i) - (I_concavs_left_septum_average-2*I_concavs_left_septum_SD);
end
I_concavs_left_septum_dev_avg_min_2SD = abs(I_concavs_left_septum_dev_avg_min_2SD);
for i = 1:length(I_concavs_left_septum_sorted)
    if I_concavs_left_septum_dev_avg_min_2SD(i) == min(I_concavs_left_septum_dev_avg_min_2SD)
        I_concavs_dev_avg_min_2SD_index = i;
    end
end
I_concavs_left_septum_dev_avg_pls_2SD = concav(length(I_concavs_left_septum_sorted),1);
for i = 1:length(I_concavs_left_septum_sorted)
    I_concavs_left_septum_dev_avg_pls_2SD(i) = I_concavs_left_septum_sorted(i) - (I_concavs_left_septum_average+2*I_concavs_left_septum_SD);
end
I_concavs_left_septum_dev_avg_pls_2SD = abs(I_concavs_left_septum_dev_avg_pls_2SD);
for i = 1:length(I_concavs_left_septum_sorted)
    if I_concavs_left_septum_dev_avg_pls_2SD(i) == min(I_concavs_left_septum_dev_avg_pls_2SD)
        I_concavs_left_septum_dev_avg_pls_2SD_index = i;
    end
end
I_concavs_left_septum_dev_avg_min_3SD = concav(length(I_concavs_left_septum_sorted),1);
for i = 1:length(I_concavs_left_septum_sorted)
    I_concavs_left_septum_dev_avg_min_3SD(i) = I_concavs_left_septum_sorted(i) - (I_concavs_left_septum_average-3*I_concavs_left_septum_SD);
end
I_concavs_left_septum_dev_avg_min_3SD = abs(I_concavs_left_septum_dev_avg_min_3SD);
for i = 1:length(I_concavs_left_septum_sorted)
    if I_concavs_left_septum_dev_avg_min_3SD(i) == min(I_concavs_left_septum_dev_avg_min_3SD)
        I_concavs_left_septum_dev_avg_min_3SD_index = i;
    end
end
I_concavs_left_septum_dev_avg_pls_3SD = concav(length(I_concavs_left_septum_sorted),1);
for i = 1:length(I_concavs_left_septum_sorted)
    I_concavs_left_septum_dev_avg_pls_3SD(i) = I_concavs_left_septum_sorted(i) - (I_concavs_left_septum_average+3*I_concavs_left_septum_SD);
end
I_concavs_left_septum_dev_avg_pls_3SD = abs(I_concavs_left_septum_dev_avg_pls_3SD);
for i = 1:length(I_concavs_left_septum_sorted)
    if I_concavs_left_septum_dev_avg_pls_3SD(i) == min(I_concavs_left_septum_dev_avg_pls_3SD)
        I_concavs_left_septum_dev_avg_pls_3SD_index = i;
    end
end
I_concavs_left_septum_tot_number_of_indices = length(I_concavs_left_septum_sorted);
%This is the approximation of the tot "area under curve." 
I_concavs_left_septum_quotient_within_1SD = (I_concavs_left_septum_dev_avg_pls_1SD_index - I_concavs_left_septumdev_avg_min_1SD_index)/I_concavs_left_septum_tot_number_of_indices;
I_concavs_left_septum_quotient_within_2SD = (I_concavs_left_septum_dev_avg_pls_2SD_index - I_concavs_dev_avg_min_2SD_index)/I_concavs_left_septum_tot_number_of_indices;
I_concavs_left_septum_quotient_within_3SD = (I_concavs_left_septum_dev_avg_pls_3SD_index - I_concavs_left_septum_dev_avg_min_3SD_index)/I_concavs_left_septum_tot_number_of_indices;
if (I_concavs_left_septum_quotient_within_1SD > 0.66 && I_concavs_left_septum_quotient_within_1SD < 0.70) && (I_concavs_left_septum_quotient_within_2SD > 0.93 && I_concavs_left_septum_quotient_within_2SD < 0.97) && (I_concavs_left_septum_quotient_within_3SD > 0.98 && I_concavs_left_septum_quotient_within_3SD < 1)
    isNormal(5,12) = 1;
end
if I_concavs_left_septum_quotient_within_1SD == 0
    isSpike(5,12) = 1;
end
I_concavs_left_atrial_appendage_sorted = sort(I_concavs_left_atrial_appendage);
I_concavs_left_atrial_appendage_dev_avg_min_1SD = concav(length(I_concavs_left_atrial_appendage_sorted),1);
for i = 1:length(I_concavs_left_atrial_appendage_sorted)
    I_concavs_left_atrial_appendage_dev_avg_min_1SD(i) = I_concavs_left_atrial_appendage_sorted(i) - (I_concavs_left_atrial_appendage_average-I_concavs_left_atrial_appendage_SD);
end
I_concavs_left_atrial_appendage_dev_avg_min_1SD = abs(I_concavs_left_atrial_appendage_dev_avg_min_1SD);
for i = 1:length(I_concavs_left_atrial_appendage_sorted)
    if I_concavs_left_atrial_appendage_dev_avg_min_1SD(i) == min(I_concavs_left_atrial_appendage_dev_avg_min_1SD)
        I_concavs_left_atrial_appendage_dev_avg_min_1SD_index = i;
    end
end
I_concavs_left_atrial_appendage_dev_avg_pls_1SD = concav(length(I_concavs_left_atrial_appendage_sorted),1);
for i = 1:length(I_concavs_left_atrial_appendage_sorted)
    I_concavs_left_atrial_appendage_dev_avg_pls_1SD(i) = I_concavs_left_atrial_appendage_sorted(i) - (I_concavs_left_atrial_appendage_average+I_concavs_left_atrial_appendage_SD);
end
I_concavs_left_atrial_appendage_dev_avg_pls_1SD = abs(I_concavs_left_atrial_appendage_dev_avg_pls_1SD);
for i = 1:length(I_concavs_left_atrial_appendage_sorted)
    if I_concavs_left_atrial_appendage_dev_avg_pls_1SD(i) == min(I_concavs_left_atrial_appendage_dev_avg_pls_1SD)
        I_concavs_left_atrial_appendage_dev_avg_pls_1SD_index = i;
    end
end
I_concavs_left_atrial_appendage_dev_avg_min_2SD = concav(length(I_concavs_left_atrial_appendage_sorted),1);
for i = 1:length(I_concavs_left_atrial_appendage_sorted)
    I_concavs_left_atrial_appendage_dev_avg_min_2SD(i) = I_concavs_left_atrial_appendage_sorted(i) - (I_concavs_left_atrial_appendage_average-2*I_concavs_left_atrial_appendage_SD);
end
I_concavs_left_atrial_appendage_dev_avg_min_2SD = abs(I_concavs_left_atrial_appendage_dev_avg_min_2SD);
for i = 1:length(I_concavs_left_atrial_appendage_sorted)
    if I_concavs_left_atrial_appendage_dev_avg_min_2SD(i) == min(I_concavs_left_atrial_appendage_dev_avg_min_2SD)
        I_concavs_dev_avg_min_2SD_index = i;
    end
end
I_concavs_left_atrial_appendage_dev_avg_pls_2SD = concav(length(I_concavs_left_atrial_appendage_sorted),1);
for i = 1:length(I_concavs_left_atrial_appendage_sorted)
    I_concavs_left_atrial_appendage_dev_avg_pls_2SD(i) = I_concavs_left_atrial_appendage_sorted(i) - (I_concavs_left_atrial_appendage_average+2*I_concavs_left_atrial_appendage_SD);
end
I_concavs_left_atrial_appendage_dev_avg_pls_2SD = abs(I_concavs_left_atrial_appendage_dev_avg_pls_2SD);
for i = 1:length(I_concavs_left_atrial_appendage_sorted)
    if I_concavs_left_atrial_appendage_dev_avg_pls_2SD(i) == min(I_concavs_left_atrial_appendage_dev_avg_pls_2SD)
        I_concavs_left_atrial_appendage_dev_avg_pls_2SD_index = i;
    end
end
I_concavs_left_atrial_appendage_dev_avg_min_3SD = concav(length(I_concavs_left_atrial_appendage_sorted),1);
for i = 1:length(I_concavs_left_atrial_appendage_sorted)
    I_concavs_left_atrial_appendage_dev_avg_min_3SD(i) = I_concavs_left_atrial_appendage_sorted(i) - (I_concavs_left_atrial_appendage_average-3*I_concavs_left_atrial_appendage_SD);
end
I_concavs_left_atrial_appendage_dev_avg_min_3SD = abs(I_concavs_left_atrial_appendage_dev_avg_min_3SD);
for i = 1:length(I_concavs_left_atrial_appendage_sorted)
    if I_concavs_left_atrial_appendage_dev_avg_min_3SD(i) == min(I_concavs_left_atrial_appendage_dev_avg_min_3SD)
        I_concavs_left_atrial_appendage_dev_avg_min_3SD_index = i;
    end
end
I_concavs_left_atrial_appendage_dev_avg_pls_3SD = concav(length(I_concavs_left_atrial_appendage_sorted),1);
for i = 1:length(I_concavs_left_atrial_appendage_sorted)
    I_concavs_left_atrial_appendage_dev_avg_pls_3SD(i) = I_concavs_left_atrial_appendage_sorted(i) - (I_concavs_left_atrial_appendage_average+3*I_concavs_left_atrial_appendage_SD);
end
I_concavs_left_atrial_appendage_dev_avg_pls_3SD = abs(I_concavs_left_atrial_appendage_dev_avg_pls_3SD);
for i = 1:length(I_concavs_left_atrial_appendage_sorted)
    if I_concavs_left_atrial_appendage_dev_avg_pls_3SD(i) == min(I_concavs_left_atrial_appendage_dev_avg_pls_3SD)
        I_concavs_left_atrial_appendage_dev_avg_pls_3SD_index = i;
    end
end
I_concavs_left_atrial_appendage_tot_number_of_indices = length(I_concavs_left_atrial_appendage_sorted);
%This is the approximation of the tot "area under curve." 
I_concavs_left_atrial_appendage_quotient_within_1SD = (I_concavs_left_atrial_appendage_dev_avg_pls_1SD_index - I_concavs_left_atrial_appendagedev_avg_min_1SD_index)/I_concavs_left_atrial_appendage_tot_number_of_indices;
I_concavs_left_atrial_appendage_quotient_within_2SD = (I_concavs_left_atrial_appendage_dev_avg_pls_2SD_index - I_concavs_dev_avg_min_2SD_index)/I_concavs_left_atrial_appendage_tot_number_of_indices;
I_concavs_left_atrial_appendage_quotient_within_3SD = (I_concavs_left_atrial_appendage_dev_avg_pls_3SD_index - I_concavs_left_atrial_appendage_dev_avg_min_3SD_index)/I_concavs_left_atrial_appendage_tot_number_of_indices;
if (I_concavs_left_atrial_appendage_quotient_within_1SD > 0.66 && I_concavs_left_atrial_appendage_quotient_within_1SD < 0.70) && (I_concavs_left_atrial_appendage_quotient_within_2SD > 0.93 && I_concavs_left_atrial_appendage_quotient_within_2SD < 0.97) && (I_concavs_left_atrial_appendage_quotient_within_3SD > 0.98 && I_concavs_left_atrial_appendage_quotient_within_3SD < 1)
    isNormal(5,13) = 1;
end
if I_concavs_left_atrial_appendage_quotient_within_1SD == 0
    isSpike(5,13) = 1;
end
I_mdslpes_SA_node_sorted = sort(I_mdslpes_SA_node);
I_mdslpes_SA_node_dev_avg_min_1SD = mdslpe(length(I_mdslpes_SA_node_sorted),1);
for i = 1:length(I_mdslpes_SA_node_sorted)
    I_mdslpes_SA_node_dev_avg_min_1SD(i) = I_mdslpes_SA_node_sorted(i) - (I_mdslpes_SA_node_average-I_mdslpes_SA_node_SD);
end
I_mdslpes_SA_node_dev_avg_min_1SD = abs(I_mdslpes_SA_node_dev_avg_min_1SD);
for i = 1:length(I_mdslpes_SA_node_sorted)
    if I_mdslpes_SA_node_dev_avg_min_1SD(i) == min(I_mdslpes_SA_node_dev_avg_min_1SD)
        I_mdslpes_SA_node_dev_avg_min_1SD_index = i;
    end
end
I_mdslpes_SA_node_dev_avg_pls_1SD = mdslpe(length(I_mdslpes_SA_node_sorted),1);
for i = 1:length(I_mdslpes_SA_node_sorted)
    I_mdslpes_SA_node_dev_avg_pls_1SD(i) = I_mdslpes_SA_node_sorted(i) - (I_mdslpes_SA_node_average+I_mdslpes_SA_node_SD);
end
I_mdslpes_SA_node_dev_avg_pls_1SD = abs(I_mdslpes_SA_node_dev_avg_pls_1SD);
for i = 1:length(I_mdslpes_SA_node_sorted)
    if I_mdslpes_SA_node_dev_avg_pls_1SD(i) == min(I_mdslpes_SA_node_dev_avg_pls_1SD)
        I_mdslpes_SA_node_dev_avg_pls_1SD_index = i;
    end
end
I_mdslpes_SA_node_dev_avg_min_2SD = mdslpe(length(I_mdslpes_SA_node_sorted),1);
for i = 1:length(I_mdslpes_SA_node_sorted)
    I_mdslpes_SA_node_dev_avg_min_2SD(i) = I_mdslpes_SA_node_sorted(i) - (I_mdslpes_SA_node_average-2*I_mdslpes_SA_node_SD);
end
I_mdslpes_SA_node_dev_avg_min_2SD = abs(I_mdslpes_SA_node_dev_avg_min_2SD);
for i = 1:length(I_mdslpes_SA_node_sorted)
    if I_mdslpes_SA_node_dev_avg_min_2SD(i) == min(I_mdslpes_SA_node_dev_avg_min_2SD)
        I_mdslpes_dev_avg_min_2SD_index = i;
    end
end
I_mdslpes_SA_node_dev_avg_pls_2SD = mdslpe(length(I_mdslpes_SA_node_sorted),1);
for i = 1:length(I_mdslpes_SA_node_sorted)
    I_mdslpes_SA_node_dev_avg_pls_2SD(i) = I_mdslpes_SA_node_sorted(i) - (I_mdslpes_SA_node_average+2*I_mdslpes_SA_node_SD);
end
I_mdslpes_SA_node_dev_avg_pls_2SD = abs(I_mdslpes_SA_node_dev_avg_pls_2SD);
for i = 1:length(I_mdslpes_SA_node_sorted)
    if I_mdslpes_SA_node_dev_avg_pls_2SD(i) == min(I_mdslpes_SA_node_dev_avg_pls_2SD)
        I_mdslpes_SA_node_dev_avg_pls_2SD_index = i;
    end
end
I_mdslpes_SA_node_dev_avg_min_3SD = mdslpe(length(I_mdslpes_SA_node_sorted),1);
for i = 1:length(I_mdslpes_SA_node_sorted)
    I_mdslpes_SA_node_dev_avg_min_3SD(i) = I_mdslpes_SA_node_sorted(i) - (I_mdslpes_SA_node_average-3*I_mdslpes_SA_node_SD);
end
I_mdslpes_SA_node_dev_avg_min_3SD = abs(I_mdslpes_SA_node_dev_avg_min_3SD);
for i = 1:length(I_mdslpes_SA_node_sorted)
    if I_mdslpes_SA_node_dev_avg_min_3SD(i) == min(I_mdslpes_SA_node_dev_avg_min_3SD)
        I_mdslpes_SA_node_dev_avg_min_3SD_index = i;
    end
end
I_mdslpes_SA_node_dev_avg_pls_3SD = mdslpe(length(I_mdslpes_SA_node_sorted),1);
for i = 1:length(I_mdslpes_SA_node_sorted)
    I_mdslpes_SA_node_dev_avg_pls_3SD(i) = I_mdslpes_SA_node_sorted(i) - (I_mdslpes_SA_node_average+3*I_mdslpes_SA_node_SD);
end
I_mdslpes_SA_node_dev_avg_pls_3SD = abs(I_mdslpes_SA_node_dev_avg_pls_3SD);
for i = 1:length(I_mdslpes_SA_node_sorted)
    if I_mdslpes_SA_node_dev_avg_pls_3SD(i) == min(I_mdslpes_SA_node_dev_avg_pls_3SD)
        I_mdslpes_SA_node_dev_avg_pls_3SD_index = i;
    end
end
I_mdslpes_SA_node_tot_number_of_indices = length(I_mdslpes_SA_node_sorted);
%This is the approximation of the tot "area under curve." 
I_mdslpes_SA_node_quotient_within_1SD = (I_mdslpes_SA_node_dev_avg_pls_1SD_index - I_mdslpes_SA_nodedev_avg_min_1SD_index)/I_mdslpes_SA_node_tot_number_of_indices;
I_mdslpes_SA_node_quotient_within_2SD = (I_mdslpes_SA_node_dev_avg_pls_2SD_index - I_mdslpes_dev_avg_min_2SD_index)/I_mdslpes_SA_node_tot_number_of_indices;
I_mdslpes_SA_node_quotient_within_3SD = (I_mdslpes_SA_node_dev_avg_pls_3SD_index - I_mdslpes_SA_node_dev_avg_min_3SD_index)/I_mdslpes_SA_node_tot_number_of_indices;
if ((I_mdslpes_SA_node_quotient_within_1SD > 0.66 && I_mdslpes_SA_node_quotient_within_1SD < 0.70) && (I_mdslpes_SA_node_quotient_within_2SD > 0.93 && I_mdslpes_SA_node_quotient_within_2SD < 0.97) && (I_mdslpes_SA_node_quotient_within_3SD > 0.98 && I_mdslpes_SA_node_quotient_within_3SD < 1)) 
    isNormal(6,1) = 1;
end
if I_mdslpes_SA_node_quotient_within_1SD == 0
    isSpike(6,1) = 1;
end
I_mdslpes_crista_terminalis_sorted = sort(I_mdslpes_crista_terminalis);
I_mdslpes_crista_terminalis_dev_avg_min_1SD = mdslpe(length(I_mdslpes_crista_terminalis_sorted),1);
for i = 1:length(I_mdslpes_crista_terminalis_sorted)
    I_mdslpes_crista_terminalis_dev_avg_min_1SD(i) = I_mdslpes_crista_terminalis_sorted(i) - (I_mdslpes_crista_terminalis_average-I_mdslpes_crista_terminalis_SD);
end
I_mdslpes_crista_terminalis_dev_avg_min_1SD = abs(I_mdslpes_crista_terminalis_dev_avg_min_1SD);
for i = 1:length(I_mdslpes_crista_terminalis_sorted)
    if I_mdslpes_crista_terminalis_dev_avg_min_1SD(i) == min(I_mdslpes_crista_terminalis_dev_avg_min_1SD)
        I_mdslpes_crista_terminalis_dev_avg_min_1SD_index = i;
    end
end
I_mdslpes_crista_terminalis_dev_avg_pls_1SD = mdslpe(length(I_mdslpes_crista_terminalis_sorted),1);
for i = 1:length(I_mdslpes_crista_terminalis_sorted)
    I_mdslpes_crista_terminalis_dev_avg_pls_1SD(i) = I_mdslpes_crista_terminalis_sorted(i) - (I_mdslpes_crista_terminalis_average+I_mdslpes_crista_terminalis_SD);
end
I_mdslpes_crista_terminalis_dev_avg_pls_1SD = abs(I_mdslpes_crista_terminalis_dev_avg_pls_1SD);
for i = 1:length(I_mdslpes_crista_terminalis_sorted)
    if I_mdslpes_crista_terminalis_dev_avg_pls_1SD(i) == min(I_mdslpes_crista_terminalis_dev_avg_pls_1SD)
        I_mdslpes_crista_terminalis_dev_avg_pls_1SD_index = i;
    end
end
I_mdslpes_crista_terminalis_dev_avg_min_2SD = mdslpe(length(I_mdslpes_crista_terminalis_sorted),1);
for i = 1:length(I_mdslpes_crista_terminalis_sorted)
    I_mdslpes_crista_terminalis_dev_avg_min_2SD(i) = I_mdslpes_crista_terminalis_sorted(i) - (I_mdslpes_crista_terminalis_average-2*I_mdslpes_crista_terminalis_SD);
end
I_mdslpes_crista_terminalis_dev_avg_min_2SD = abs(I_mdslpes_crista_terminalis_dev_avg_min_2SD);
for i = 1:length(I_mdslpes_crista_terminalis_sorted)
    if I_mdslpes_crista_terminalis_dev_avg_min_2SD(i) == min(I_mdslpes_crista_terminalis_dev_avg_min_2SD)
        I_mdslpes_dev_avg_min_2SD_index = i;
    end
end
I_mdslpes_crista_terminalis_dev_avg_pls_2SD = mdslpe(length(I_mdslpes_crista_terminalis_sorted),1);
for i = 1:length(I_mdslpes_crista_terminalis_sorted)
    I_mdslpes_crista_terminalis_dev_avg_pls_2SD(i) = I_mdslpes_crista_terminalis_sorted(i) - (I_mdslpes_crista_terminalis_average+2*I_mdslpes_crista_terminalis_SD);
end
I_mdslpes_crista_terminalis_dev_avg_pls_2SD = abs(I_mdslpes_crista_terminalis_dev_avg_pls_2SD);
for i = 1:length(I_mdslpes_crista_terminalis_sorted)
    if I_mdslpes_crista_terminalis_dev_avg_pls_2SD(i) == min(I_mdslpes_crista_terminalis_dev_avg_pls_2SD)
        I_mdslpes_crista_terminalis_dev_avg_pls_2SD_index = i;
    end
end
I_mdslpes_crista_terminalis_dev_avg_min_3SD = mdslpe(length(I_mdslpes_crista_terminalis_sorted),1);
for i = 1:length(I_mdslpes_crista_terminalis_sorted)
    I_mdslpes_crista_terminalis_dev_avg_min_3SD(i) = I_mdslpes_crista_terminalis_sorted(i) - (I_mdslpes_crista_terminalis_average-3*I_mdslpes_crista_terminalis_SD);
end
I_mdslpes_crista_terminalis_dev_avg_min_3SD = abs(I_mdslpes_crista_terminalis_dev_avg_min_3SD);
for i = 1:length(I_mdslpes_crista_terminalis_sorted)
    if I_mdslpes_crista_terminalis_dev_avg_min_3SD(i) == min(I_mdslpes_crista_terminalis_dev_avg_min_3SD)
        I_mdslpes_crista_terminalis_dev_avg_min_3SD_index = i;
    end
end
I_mdslpes_crista_terminalis_dev_avg_pls_3SD = mdslpe(length(I_mdslpes_crista_terminalis_sorted),1);
for i = 1:length(I_mdslpes_crista_terminalis_sorted)
    I_mdslpes_crista_terminalis_dev_avg_pls_3SD(i) = I_mdslpes_crista_terminalis_sorted(i) - (I_mdslpes_crista_terminalis_average+3*I_mdslpes_crista_terminalis_SD);
end
I_mdslpes_crista_terminalis_dev_avg_pls_3SD = abs(I_mdslpes_crista_terminalis_dev_avg_pls_3SD);
for i = 1:length(I_mdslpes_crista_terminalis_sorted)
    if I_mdslpes_crista_terminalis_dev_avg_pls_3SD(i) == min(I_mdslpes_crista_terminalis_dev_avg_pls_3SD)
        I_mdslpes_crista_terminalis_dev_avg_pls_3SD_index = i;
    end
end
I_mdslpes_crista_terminalis_tot_number_of_indices = length(I_mdslpes_crista_terminalis_sorted);
%This is the approximation of the tot "area under curve." 
I_mdslpes_crista_terminalis_quotient_within_1SD = (I_mdslpes_crista_terminalis_dev_avg_pls_1SD_index - I_mdslpes_crista_terminalisdev_avg_min_1SD_index)/I_mdslpes_crista_terminalis_tot_number_of_indices;
I_mdslpes_crista_terminalis_quotient_within_2SD = (I_mdslpes_crista_terminalis_dev_avg_pls_2SD_index - I_mdslpes_dev_avg_min_2SD_index)/I_mdslpes_crista_terminalis_tot_number_of_indices;
I_mdslpes_crista_terminalis_quotient_within_3SD = (I_mdslpes_crista_terminalis_dev_avg_pls_3SD_index - I_mdslpes_crista_terminalis_dev_avg_min_3SD_index)/I_mdslpes_crista_terminalis_tot_number_of_indices;
if (I_mdslpes_crista_terminalis_quotient_within_1SD > 0.66 && I_mdslpes_crista_terminalis_quotient_within_1SD < 0.70) && (I_mdslpes_crista_terminalis_quotient_within_2SD > 0.93 && I_mdslpes_crista_terminalis_quotient_within_2SD < 0.97) && (I_mdslpes_crista_terminalis_quotient_within_3SD > 0.98 && I_mdslpes_crista_terminalis_quotient_within_3SD < 1)
    isNormal(6,2) = 1;
end
if I_mdslpes_crista_terminalis_quotient_within_1SD == 0
    isSpike(6,2) = 1;
end
I_mdslpes_tricuspid_annulus_sorted = sort(I_mdslpes_tricuspid_annulus);
I_mdslpes_tricuspid_annulus_dev_avg_min_1SD = mdslpe(length(I_mdslpes_tricuspid_annulus_sorted),1);
for i = 1:length(I_mdslpes_tricuspid_annulus_sorted)
    I_mdslpes_tricuspid_annulus_dev_avg_min_1SD(i) = I_mdslpes_tricuspid_annulus_sorted(i) - (I_mdslpes_tricuspid_annulus_average-I_mdslpes_tricuspid_annulus_SD);
end
I_mdslpes_tricuspid_annulus_dev_avg_min_1SD = abs(I_mdslpes_tricuspid_annulus_dev_avg_min_1SD);
for i = 1:length(I_mdslpes_tricuspid_annulus_sorted)
    if I_mdslpes_tricuspid_annulus_dev_avg_min_1SD(i) == min(I_mdslpes_tricuspid_annulus_dev_avg_min_1SD)
        I_mdslpes_tricuspid_annulus_dev_avg_min_1SD_index = i;
    end
end
I_mdslpes_tricuspid_annulus_dev_avg_pls_1SD = mdslpe(length(I_mdslpes_tricuspid_annulus_sorted),1);
for i = 1:length(I_mdslpes_tricuspid_annulus_sorted)
    I_mdslpes_tricuspid_annulus_dev_avg_pls_1SD(i) = I_mdslpes_tricuspid_annulus_sorted(i) - (I_mdslpes_tricuspid_annulus_average+I_mdslpes_tricuspid_annulus_SD);
end
I_mdslpes_tricuspid_annulus_dev_avg_pls_1SD = abs(I_mdslpes_tricuspid_annulus_dev_avg_pls_1SD);
for i = 1:length(I_mdslpes_tricuspid_annulus_sorted)
    if I_mdslpes_tricuspid_annulus_dev_avg_pls_1SD(i) == min(I_mdslpes_tricuspid_annulus_dev_avg_pls_1SD)
        I_mdslpes_tricuspid_annulus_dev_avg_pls_1SD_index = i;
    end
end
I_mdslpes_tricuspid_annulus_dev_avg_min_2SD = mdslpe(length(I_mdslpes_tricuspid_annulus_sorted),1);
for i = 1:length(I_mdslpes_tricuspid_annulus_sorted)
    I_mdslpes_tricuspid_annulus_dev_avg_min_2SD(i) = I_mdslpes_tricuspid_annulus_sorted(i) - (I_mdslpes_tricuspid_annulus_average-2*I_mdslpes_tricuspid_annulus_SD);
end
I_mdslpes_tricuspid_annulus_dev_avg_min_2SD = abs(I_mdslpes_tricuspid_annulus_dev_avg_min_2SD);
for i = 1:length(I_mdslpes_tricuspid_annulus_sorted)
    if I_mdslpes_tricuspid_annulus_dev_avg_min_2SD(i) == min(I_mdslpes_tricuspid_annulus_dev_avg_min_2SD)
        I_mdslpes_dev_avg_min_2SD_index = i;
    end
end
I_mdslpes_tricuspid_annulus_dev_avg_pls_2SD = mdslpe(length(I_mdslpes_tricuspid_annulus_sorted),1);
for i = 1:length(I_mdslpes_tricuspid_annulus_sorted)
    I_mdslpes_tricuspid_annulus_dev_avg_pls_2SD(i) = I_mdslpes_tricuspid_annulus_sorted(i) - (I_mdslpes_tricuspid_annulus_average+2*I_mdslpes_tricuspid_annulus_SD);
end
I_mdslpes_tricuspid_annulus_dev_avg_pls_2SD = abs(I_mdslpes_tricuspid_annulus_dev_avg_pls_2SD);
for i = 1:length(I_mdslpes_tricuspid_annulus_sorted)
    if I_mdslpes_tricuspid_annulus_dev_avg_pls_2SD(i) == min(I_mdslpes_tricuspid_annulus_dev_avg_pls_2SD)
        I_mdslpes_tricuspid_annulus_dev_avg_pls_2SD_index = i;
    end
end
I_mdslpes_tricuspid_annulus_dev_avg_min_3SD = mdslpe(length(I_mdslpes_tricuspid_annulus_sorted),1);
for i = 1:length(I_mdslpes_tricuspid_annulus_sorted)
    I_mdslpes_tricuspid_annulus_dev_avg_min_3SD(i) = I_mdslpes_tricuspid_annulus_sorted(i) - (I_mdslpes_tricuspid_annulus_average-3*I_mdslpes_tricuspid_annulus_SD);
end
I_mdslpes_tricuspid_annulus_dev_avg_min_3SD = abs(I_mdslpes_tricuspid_annulus_dev_avg_min_3SD);
for i = 1:length(I_mdslpes_tricuspid_annulus_sorted)
    if I_mdslpes_tricuspid_annulus_dev_avg_min_3SD(i) == min(I_mdslpes_tricuspid_annulus_dev_avg_min_3SD)
        I_mdslpes_tricuspid_annulus_dev_avg_min_3SD_index = i;
    end
end
I_mdslpes_tricuspid_annulus_dev_avg_pls_3SD = mdslpe(length(I_mdslpes_tricuspid_annulus_sorted),1);
for i = 1:length(I_mdslpes_tricuspid_annulus_sorted)
    I_mdslpes_tricuspid_annulus_dev_avg_pls_3SD(i) = I_mdslpes_tricuspid_annulus_sorted(i) - (I_mdslpes_tricuspid_annulus_average+3*I_mdslpes_tricuspid_annulus_SD);
end
I_mdslpes_tricuspid_annulus_dev_avg_pls_3SD = abs(I_mdslpes_tricuspid_annulus_dev_avg_pls_3SD);
for i = 1:length(I_mdslpes_tricuspid_annulus_sorted)
    if I_mdslpes_tricuspid_annulus_dev_avg_pls_3SD(i) == min(I_mdslpes_tricuspid_annulus_dev_avg_pls_3SD)
        I_mdslpes_tricuspid_annulus_dev_avg_pls_3SD_index = i;
    end
end
I_mdslpes_tricuspid_annulus_tot_number_of_indices = length(I_mdslpes_tricuspid_annulus_sorted);
%This is the approximation of the tot "area under curve." 
I_mdslpes_tricuspid_annulus_quotient_within_1SD = (I_mdslpes_tricuspid_annulus_dev_avg_pls_1SD_index - I_mdslpes_tricuspid_annulusdev_avg_min_1SD_index)/I_mdslpes_tricuspid_annulus_tot_number_of_indices;
I_mdslpes_tricuspid_annulus_quotient_within_2SD = (I_mdslpes_tricuspid_annulus_dev_avg_pls_2SD_index - I_mdslpes_dev_avg_min_2SD_index)/I_mdslpes_tricuspid_annulus_tot_number_of_indices;
I_mdslpes_tricuspid_annulus_quotient_within_3SD = (I_mdslpes_tricuspid_annulus_dev_avg_pls_3SD_index - I_mdslpes_tricuspid_annulus_dev_avg_min_3SD_index)/I_mdslpes_tricuspid_annulus_tot_number_of_indices;
if (I_mdslpes_tricuspid_annulus_quotient_within_1SD > 0.66 && I_mdslpes_tricuspid_annulus_quotient_within_1SD < 0.70) && (I_mdslpes_tricuspid_annulus_quotient_within_2SD > 0.93 && I_mdslpes_tricuspid_annulus_quotient_within_2SD < 0.97) && (I_mdslpes_tricuspid_annulus_quotient_within_3SD > 0.98 && I_mdslpes_tricuspid_annulus_quotient_within_3SD < 1)
    isNormal(6,3) = 1;
end
if I_mdslpes_tricuspid_annulus_quotient_within_1SD == 0
    isSpike(6,3) = 1;
end
I_mdslpes_coronary_sinus_sorted = sort(I_mdslpes_coronary_sinus);
I_mdslpes_coronary_sinus_dev_avg_min_1SD = mdslpe(length(I_mdslpes_coronary_sinus_sorted),1);
for i = 1:length(I_mdslpes_coronary_sinus_sorted)
    I_mdslpes_coronary_sinus_dev_avg_min_1SD(i) = I_mdslpes_coronary_sinus_sorted(i) - (I_mdslpes_coronary_sinus_average-I_mdslpes_coronary_sinus_SD);
end
I_mdslpes_coronary_sinus_dev_avg_min_1SD = abs(I_mdslpes_coronary_sinus_dev_avg_min_1SD);
for i = 1:length(I_mdslpes_coronary_sinus_sorted)
    if I_mdslpes_coronary_sinus_dev_avg_min_1SD(i) == min(I_mdslpes_coronary_sinus_dev_avg_min_1SD)
        I_mdslpes_coronary_sinus_dev_avg_min_1SD_index = i;
    end
end
I_mdslpes_coronary_sinus_dev_avg_pls_1SD = mdslpe(length(I_mdslpes_coronary_sinus_sorted),1);
for i = 1:length(I_mdslpes_coronary_sinus_sorted)
    I_mdslpes_coronary_sinus_dev_avg_pls_1SD(i) = I_mdslpes_coronary_sinus_sorted(i) - (I_mdslpes_coronary_sinus_average+I_mdslpes_coronary_sinus_SD);
end
I_mdslpes_coronary_sinus_dev_avg_pls_1SD = abs(I_mdslpes_coronary_sinus_dev_avg_pls_1SD);
for i = 1:length(I_mdslpes_coronary_sinus_sorted)
    if I_mdslpes_coronary_sinus_dev_avg_pls_1SD(i) == min(I_mdslpes_coronary_sinus_dev_avg_pls_1SD)
        I_mdslpes_coronary_sinus_dev_avg_pls_1SD_index = i;
    end
end
I_mdslpes_coronary_sinus_dev_avg_min_2SD = mdslpe(length(I_mdslpes_coronary_sinus_sorted),1);
for i = 1:length(I_mdslpes_coronary_sinus_sorted)
    I_mdslpes_coronary_sinus_dev_avg_min_2SD(i) = I_mdslpes_coronary_sinus_sorted(i) - (I_mdslpes_coronary_sinus_average-2*I_mdslpes_coronary_sinus_SD);
end
I_mdslpes_coronary_sinus_dev_avg_min_2SD = abs(I_mdslpes_coronary_sinus_dev_avg_min_2SD);
for i = 1:length(I_mdslpes_coronary_sinus_sorted)
    if I_mdslpes_coronary_sinus_dev_avg_min_2SD(i) == min(I_mdslpes_coronary_sinus_dev_avg_min_2SD)
        I_mdslpes_dev_avg_min_2SD_index = i;
    end
end
I_mdslpes_coronary_sinus_dev_avg_pls_2SD = mdslpe(length(I_mdslpes_coronary_sinus_sorted),1);
for i = 1:length(I_mdslpes_coronary_sinus_sorted)
    I_mdslpes_coronary_sinus_dev_avg_pls_2SD(i) = I_mdslpes_coronary_sinus_sorted(i) - (I_mdslpes_coronary_sinus_average+2*I_mdslpes_coronary_sinus_SD);
end
I_mdslpes_coronary_sinus_dev_avg_pls_2SD = abs(I_mdslpes_coronary_sinus_dev_avg_pls_2SD);
for i = 1:length(I_mdslpes_coronary_sinus_sorted)
    if I_mdslpes_coronary_sinus_dev_avg_pls_2SD(i) == min(I_mdslpes_coronary_sinus_dev_avg_pls_2SD)
        I_mdslpes_coronary_sinus_dev_avg_pls_2SD_index = i;
    end
end
I_mdslpes_coronary_sinus_dev_avg_min_3SD = mdslpe(length(I_mdslpes_coronary_sinus_sorted),1);
for i = 1:length(I_mdslpes_coronary_sinus_sorted)
    I_mdslpes_coronary_sinus_dev_avg_min_3SD(i) = I_mdslpes_coronary_sinus_sorted(i) - (I_mdslpes_coronary_sinus_average-3*I_mdslpes_coronary_sinus_SD);
end
I_mdslpes_coronary_sinus_dev_avg_min_3SD = abs(I_mdslpes_coronary_sinus_dev_avg_min_3SD);
for i = 1:length(I_mdslpes_coronary_sinus_sorted)
    if I_mdslpes_coronary_sinus_dev_avg_min_3SD(i) == min(I_mdslpes_coronary_sinus_dev_avg_min_3SD)
        I_mdslpes_coronary_sinus_dev_avg_min_3SD_index = i;
    end
end
I_mdslpes_coronary_sinus_dev_avg_pls_3SD = mdslpe(length(I_mdslpes_coronary_sinus_sorted),1);
for i = 1:length(I_mdslpes_coronary_sinus_sorted)
    I_mdslpes_coronary_sinus_dev_avg_pls_3SD(i) = I_mdslpes_coronary_sinus_sorted(i) - (I_mdslpes_coronary_sinus_average+3*I_mdslpes_coronary_sinus_SD);
end
I_mdslpes_coronary_sinus_dev_avg_pls_3SD = abs(I_mdslpes_coronary_sinus_dev_avg_pls_3SD);
for i = 1:length(I_mdslpes_coronary_sinus_sorted)
    if I_mdslpes_coronary_sinus_dev_avg_pls_3SD(i) == min(I_mdslpes_coronary_sinus_dev_avg_pls_3SD)
        I_mdslpes_coronary_sinus_dev_avg_pls_3SD_index = i;
    end
end
I_mdslpes_coronary_sinus_tot_number_of_indices = length(I_mdslpes_coronary_sinus_sorted);
%This is the approximation of the tot "area under curve." 
I_mdslpes_coronary_sinus_quotient_within_1SD = (I_mdslpes_coronary_sinus_dev_avg_pls_1SD_index - I_mdslpes_coronary_sinusdev_avg_min_1SD_index)/I_mdslpes_coronary_sinus_tot_number_of_indices;
I_mdslpes_coronary_sinus_quotient_within_2SD = (I_mdslpes_coronary_sinus_dev_avg_pls_2SD_index - I_mdslpes_dev_avg_min_2SD_index)/I_mdslpes_coronary_sinus_tot_number_of_indices;
I_mdslpes_coronary_sinus_quotient_within_3SD = (I_mdslpes_coronary_sinus_dev_avg_pls_3SD_index - I_mdslpes_coronary_sinus_dev_avg_min_3SD_index)/I_mdslpes_coronary_sinus_tot_number_of_indices;
if (I_mdslpes_coronary_sinus_quotient_within_1SD > 0.66 && I_mdslpes_coronary_sinus_quotient_within_1SD < 0.70) && (I_mdslpes_coronary_sinus_quotient_within_2SD > 0.93 && I_mdslpes_coronary_sinus_quotient_within_2SD < 0.97) && (I_mdslpes_coronary_sinus_quotient_within_3SD > 0.98 && I_mdslpes_coronary_sinus_quotient_within_3SD < 1)
    isNormal(6,4) = 1;
end
if I_mdslpes_coronary_sinus_quotient_within_1SD == 0
    isSpike(6,4) = 1;
end
I_mdslpes_ostium_sorted = sort(I_mdslpes_ostium);
I_mdslpes_ostium_dev_avg_min_1SD = mdslpe(length(I_mdslpes_ostium_sorted),1);
for i = 1:length(I_mdslpes_ostium_sorted)
    I_mdslpes_ostium_dev_avg_min_1SD(i) = I_mdslpes_ostium_sorted(i) - (I_mdslpes_ostium_average-I_mdslpes_ostium_SD);
end
I_mdslpes_ostium_dev_avg_min_1SD = abs(I_mdslpes_ostium_dev_avg_min_1SD);
for i = 1:length(I_mdslpes_ostium_sorted)
    if I_mdslpes_ostium_dev_avg_min_1SD(i) == min(I_mdslpes_ostium_dev_avg_min_1SD)
        I_mdslpes_ostium_dev_avg_min_1SD_index = i;
    end
end
I_mdslpes_ostium_dev_avg_pls_1SD = mdslpe(length(I_mdslpes_ostium_sorted),1);
for i = 1:length(I_mdslpes_ostium_sorted)
    I_mdslpes_ostium_dev_avg_pls_1SD(i) = I_mdslpes_ostium_sorted(i) - (I_mdslpes_ostium_average+I_mdslpes_ostium_SD);
end
I_mdslpes_ostium_dev_avg_pls_1SD = abs(I_mdslpes_ostium_dev_avg_pls_1SD);
for i = 1:length(I_mdslpes_ostium_sorted)
    if I_mdslpes_ostium_dev_avg_pls_1SD(i) == min(I_mdslpes_ostium_dev_avg_pls_1SD)
        I_mdslpes_ostium_dev_avg_pls_1SD_index = i;
    end
end
I_mdslpes_ostium_dev_avg_min_2SD = mdslpe(length(I_mdslpes_ostium_sorted),1);
for i = 1:length(I_mdslpes_ostium_sorted)
    I_mdslpes_ostium_dev_avg_min_2SD(i) = I_mdslpes_ostium_sorted(i) - (I_mdslpes_ostium_average-2*I_mdslpes_ostium_SD);
end
I_mdslpes_ostium_dev_avg_min_2SD = abs(I_mdslpes_ostium_dev_avg_min_2SD);
for i = 1:length(I_mdslpes_ostium_sorted)
    if I_mdslpes_ostium_dev_avg_min_2SD(i) == min(I_mdslpes_ostium_dev_avg_min_2SD)
        I_mdslpes_dev_avg_min_2SD_index = i;
    end
end
I_mdslpes_ostium_dev_avg_pls_2SD = mdslpe(length(I_mdslpes_ostium_sorted),1);
for i = 1:length(I_mdslpes_ostium_sorted)
    I_mdslpes_ostium_dev_avg_pls_2SD(i) = I_mdslpes_ostium_sorted(i) - (I_mdslpes_ostium_average+2*I_mdslpes_ostium_SD);
end
I_mdslpes_ostium_dev_avg_pls_2SD = abs(I_mdslpes_ostium_dev_avg_pls_2SD);
for i = 1:length(I_mdslpes_ostium_sorted)
    if I_mdslpes_ostium_dev_avg_pls_2SD(i) == min(I_mdslpes_ostium_dev_avg_pls_2SD)
        I_mdslpes_ostium_dev_avg_pls_2SD_index = i;
    end
end
I_mdslpes_ostium_dev_avg_min_3SD = mdslpe(length(I_mdslpes_ostium_sorted),1);
for i = 1:length(I_mdslpes_ostium_sorted)
    I_mdslpes_ostium_dev_avg_min_3SD(i) = I_mdslpes_ostium_sorted(i) - (I_mdslpes_ostium_average-3*I_mdslpes_ostium_SD);
end
I_mdslpes_ostium_dev_avg_min_3SD = abs(I_mdslpes_ostium_dev_avg_min_3SD);
for i = 1:length(I_mdslpes_ostium_sorted)
    if I_mdslpes_ostium_dev_avg_min_3SD(i) == min(I_mdslpes_ostium_dev_avg_min_3SD)
        I_mdslpes_ostium_dev_avg_min_3SD_index = i;
    end
end
I_mdslpes_ostium_dev_avg_pls_3SD = mdslpe(length(I_mdslpes_ostium_sorted),1);
for i = 1:length(I_mdslpes_ostium_sorted)
    I_mdslpes_ostium_dev_avg_pls_3SD(i) = I_mdslpes_ostium_sorted(i) - (I_mdslpes_ostium_average+3*I_mdslpes_ostium_SD);
end
I_mdslpes_ostium_dev_avg_pls_3SD = abs(I_mdslpes_ostium_dev_avg_pls_3SD);
for i = 1:length(I_mdslpes_ostium_sorted)
    if I_mdslpes_ostium_dev_avg_pls_3SD(i) == min(I_mdslpes_ostium_dev_avg_pls_3SD)
        I_mdslpes_ostium_dev_avg_pls_3SD_index = i;
    end
end
I_mdslpes_ostium_tot_number_of_indices = length(I_mdslpes_ostium_sorted);
%This is the approximation of the tot "area under curve." 
I_mdslpes_ostium_quotient_within_1SD = (I_mdslpes_ostium_dev_avg_pls_1SD_index - I_mdslpes_ostiumdev_avg_min_1SD_index)/I_mdslpes_ostium_tot_number_of_indices;
I_mdslpes_ostium_quotient_within_2SD = (I_mdslpes_ostium_dev_avg_pls_2SD_index - I_mdslpes_dev_avg_min_2SD_index)/I_mdslpes_ostium_tot_number_of_indices;
I_mdslpes_ostium_quotient_within_3SD = (I_mdslpes_ostium_dev_avg_pls_3SD_index - I_mdslpes_ostium_dev_avg_min_3SD_index)/I_mdslpes_ostium_tot_number_of_indices;
if (I_mdslpes_ostium_quotient_within_1SD > 0.66 && I_mdslpes_ostium_quotient_within_1SD < 0.70) && (I_mdslpes_ostium_quotient_within_2SD > 0.93 && I_mdslpes_ostium_quotient_within_2SD < 0.97) && (I_mdslpes_ostium_quotient_within_3SD > 0.98 && I_mdslpes_ostium_quotient_within_3SD < 1)
    isNormal(6,5) = 1;
end
if I_mdslpes_ostium_quotient_within_1SD == 0
    isSpike(6,5) = 1;
end
I_mdslpes_perinodal_sorted = sort(I_mdslpes_perinodal);
I_mdslpes_perinodal_dev_avg_min_1SD = mdslpe(length(I_mdslpes_perinodal_sorted),1);
for i = 1:length(I_mdslpes_perinodal_sorted)
    I_mdslpes_perinodal_dev_avg_min_1SD(i) = I_mdslpes_perinodal_sorted(i) - (I_mdslpes_perinodal_average-I_mdslpes_perinodal_SD);
end
I_mdslpes_perinodal_dev_avg_min_1SD = abs(I_mdslpes_perinodal_dev_avg_min_1SD);
for i = 1:length(I_mdslpes_perinodal_sorted)
    if I_mdslpes_perinodal_dev_avg_min_1SD(i) == min(I_mdslpes_perinodal_dev_avg_min_1SD)
        I_mdslpes_perinodal_dev_avg_min_1SD_index = i;
    end
end
I_mdslpes_perinodal_dev_avg_pls_1SD = mdslpe(length(I_mdslpes_perinodal_sorted),1);
for i = 1:length(I_mdslpes_perinodal_sorted)
    I_mdslpes_perinodal_dev_avg_pls_1SD(i) = I_mdslpes_perinodal_sorted(i) - (I_mdslpes_perinodal_average+I_mdslpes_perinodal_SD);
end
I_mdslpes_perinodal_dev_avg_pls_1SD = abs(I_mdslpes_perinodal_dev_avg_pls_1SD);
for i = 1:length(I_mdslpes_perinodal_sorted)
    if I_mdslpes_perinodal_dev_avg_pls_1SD(i) == min(I_mdslpes_perinodal_dev_avg_pls_1SD)
        I_mdslpes_perinodal_dev_avg_pls_1SD_index = i;
    end
end
I_mdslpes_perinodal_dev_avg_min_2SD = mdslpe(length(I_mdslpes_perinodal_sorted),1);
for i = 1:length(I_mdslpes_perinodal_sorted)
    I_mdslpes_perinodal_dev_avg_min_2SD(i) = I_mdslpes_perinodal_sorted(i) - (I_mdslpes_perinodal_average-2*I_mdslpes_perinodal_SD);
end
I_mdslpes_perinodal_dev_avg_min_2SD = abs(I_mdslpes_perinodal_dev_avg_min_2SD);
for i = 1:length(I_mdslpes_perinodal_sorted)
    if I_mdslpes_perinodal_dev_avg_min_2SD(i) == min(I_mdslpes_perinodal_dev_avg_min_2SD)
        I_mdslpes_dev_avg_min_2SD_index = i;
    end
end
I_mdslpes_perinodal_dev_avg_pls_2SD = mdslpe(length(I_mdslpes_perinodal_sorted),1);
for i = 1:length(I_mdslpes_perinodal_sorted)
    I_mdslpes_perinodal_dev_avg_pls_2SD(i) = I_mdslpes_perinodal_sorted(i) - (I_mdslpes_perinodal_average+2*I_mdslpes_perinodal_SD);
end
I_mdslpes_perinodal_dev_avg_pls_2SD = abs(I_mdslpes_perinodal_dev_avg_pls_2SD);
for i = 1:length(I_mdslpes_perinodal_sorted)
    if I_mdslpes_perinodal_dev_avg_pls_2SD(i) == min(I_mdslpes_perinodal_dev_avg_pls_2SD)
        I_mdslpes_perinodal_dev_avg_pls_2SD_index = i;
    end
end
I_mdslpes_perinodal_dev_avg_min_3SD = mdslpe(length(I_mdslpes_perinodal_sorted),1);
for i = 1:length(I_mdslpes_perinodal_sorted)
    I_mdslpes_perinodal_dev_avg_min_3SD(i) = I_mdslpes_perinodal_sorted(i) - (I_mdslpes_perinodal_average-3*I_mdslpes_perinodal_SD);
end
I_mdslpes_perinodal_dev_avg_min_3SD = abs(I_mdslpes_perinodal_dev_avg_min_3SD);
for i = 1:length(I_mdslpes_perinodal_sorted)
    if I_mdslpes_perinodal_dev_avg_min_3SD(i) == min(I_mdslpes_perinodal_dev_avg_min_3SD)
        I_mdslpes_perinodal_dev_avg_min_3SD_index = i;
    end
end
I_mdslpes_perinodal_dev_avg_pls_3SD = mdslpe(length(I_mdslpes_perinodal_sorted),1);
for i = 1:length(I_mdslpes_perinodal_sorted)
    I_mdslpes_perinodal_dev_avg_pls_3SD(i) = I_mdslpes_perinodal_sorted(i) - (I_mdslpes_perinodal_average+3*I_mdslpes_perinodal_SD);
end
I_mdslpes_perinodal_dev_avg_pls_3SD = abs(I_mdslpes_perinodal_dev_avg_pls_3SD);
for i = 1:length(I_mdslpes_perinodal_sorted)
    if I_mdslpes_perinodal_dev_avg_pls_3SD(i) == min(I_mdslpes_perinodal_dev_avg_pls_3SD)
        I_mdslpes_perinodal_dev_avg_pls_3SD_index = i;
    end
end
I_mdslpes_perinodal_tot_number_of_indices = length(I_mdslpes_perinodal_sorted);
%This is the approximation of the tot "area under curve." 
I_mdslpes_perinodal_quotient_within_1SD = (I_mdslpes_perinodal_dev_avg_pls_1SD_index - I_mdslpes_perinodaldev_avg_min_1SD_index)/I_mdslpes_perinodal_tot_number_of_indices;
I_mdslpes_perinodal_quotient_within_2SD = (I_mdslpes_perinodal_dev_avg_pls_2SD_index - I_mdslpes_dev_avg_min_2SD_index)/I_mdslpes_perinodal_tot_number_of_indices;
I_mdslpes_perinodal_quotient_within_3SD = (I_mdslpes_perinodal_dev_avg_pls_3SD_index - I_mdslpes_perinodal_dev_avg_min_3SD_index)/I_mdslpes_perinodal_tot_number_of_indices;
if (I_mdslpes_perinodal_quotient_within_1SD > 0.66 && I_mdslpes_perinodal_quotient_within_1SD < 0.70) && (I_mdslpes_perinodal_quotient_within_2SD > 0.93 && I_mdslpes_perinodal_quotient_within_2SD < 0.97) && (I_mdslpes_perinodal_quotient_within_3SD > 0.98 && I_mdslpes_perinodal_quotient_within_3SD < 1)
    isNormal(6,6) = 1;
end
if I_mdslpes_perinodal_quotient_within_1SD == 0
    isSpike(6,6) = 1;
end
I_mdslpes_right_septum_sorted = sort(I_mdslpes_right_septum);
I_mdslpes_right_septum_dev_avg_min_1SD = mdslpe(length(I_mdslpes_right_septum_sorted),1);
for i = 1:length(I_mdslpes_right_septum_sorted)
    I_mdslpes_right_septum_dev_avg_min_1SD(i) = I_mdslpes_right_septum_sorted(i) - (I_mdslpes_right_septum_average-I_mdslpes_right_septum_SD);
end
I_mdslpes_right_septum_dev_avg_min_1SD = abs(I_mdslpes_right_septum_dev_avg_min_1SD);
for i = 1:length(I_mdslpes_right_septum_sorted)
    if I_mdslpes_right_septum_dev_avg_min_1SD(i) == min(I_mdslpes_right_septum_dev_avg_min_1SD)
        I_mdslpes_right_septum_dev_avg_min_1SD_index = i;
    end
end
I_mdslpes_right_septum_dev_avg_pls_1SD = mdslpe(length(I_mdslpes_right_septum_sorted),1);
for i = 1:length(I_mdslpes_right_septum_sorted)
    I_mdslpes_right_septum_dev_avg_pls_1SD(i) = I_mdslpes_right_septum_sorted(i) - (I_mdslpes_right_septum_average+I_mdslpes_right_septum_SD);
end
I_mdslpes_right_septum_dev_avg_pls_1SD = abs(I_mdslpes_right_septum_dev_avg_pls_1SD);
for i = 1:length(I_mdslpes_right_septum_sorted)
    if I_mdslpes_right_septum_dev_avg_pls_1SD(i) == min(I_mdslpes_right_septum_dev_avg_pls_1SD)
        I_mdslpes_right_septum_dev_avg_pls_1SD_index = i;
    end
end
I_mdslpes_right_septum_dev_avg_min_2SD = mdslpe(length(I_mdslpes_right_septum_sorted),1);
for i = 1:length(I_mdslpes_right_septum_sorted)
    I_mdslpes_right_septum_dev_avg_min_2SD(i) = I_mdslpes_right_septum_sorted(i) - (I_mdslpes_right_septum_average-2*I_mdslpes_right_septum_SD);
end
I_mdslpes_right_septum_dev_avg_min_2SD = abs(I_mdslpes_right_septum_dev_avg_min_2SD);
for i = 1:length(I_mdslpes_right_septum_sorted)
    if I_mdslpes_right_septum_dev_avg_min_2SD(i) == min(I_mdslpes_right_septum_dev_avg_min_2SD)
        I_mdslpes_dev_avg_min_2SD_index = i;
    end
end
I_mdslpes_right_septum_dev_avg_pls_2SD = mdslpe(length(I_mdslpes_right_septum_sorted),1);
for i = 1:length(I_mdslpes_right_septum_sorted)
    I_mdslpes_right_septum_dev_avg_pls_2SD(i) = I_mdslpes_right_septum_sorted(i) - (I_mdslpes_right_septum_average+2*I_mdslpes_right_septum_SD);
end
I_mdslpes_right_septum_dev_avg_pls_2SD = abs(I_mdslpes_right_septum_dev_avg_pls_2SD);
for i = 1:length(I_mdslpes_right_septum_sorted)
    if I_mdslpes_right_septum_dev_avg_pls_2SD(i) == min(I_mdslpes_right_septum_dev_avg_pls_2SD)
        I_mdslpes_right_septum_dev_avg_pls_2SD_index = i;
    end
end
I_mdslpes_right_septum_dev_avg_min_3SD = mdslpe(length(I_mdslpes_right_septum_sorted),1);
for i = 1:length(I_mdslpes_right_septum_sorted)
    I_mdslpes_right_septum_dev_avg_min_3SD(i) = I_mdslpes_right_septum_sorted(i) - (I_mdslpes_right_septum_average-3*I_mdslpes_right_septum_SD);
end
I_mdslpes_right_septum_dev_avg_min_3SD = abs(I_mdslpes_right_septum_dev_avg_min_3SD);
for i = 1:length(I_mdslpes_right_septum_sorted)
    if I_mdslpes_right_septum_dev_avg_min_3SD(i) == min(I_mdslpes_right_septum_dev_avg_min_3SD)
        I_mdslpes_right_septum_dev_avg_min_3SD_index = i;
    end
end
I_mdslpes_right_septum_dev_avg_pls_3SD = mdslpe(length(I_mdslpes_right_septum_sorted),1);
for i = 1:length(I_mdslpes_right_septum_sorted)
    I_mdslpes_right_septum_dev_avg_pls_3SD(i) = I_mdslpes_right_septum_sorted(i) - (I_mdslpes_right_septum_average+3*I_mdslpes_right_septum_SD);
end
I_mdslpes_right_septum_dev_avg_pls_3SD = abs(I_mdslpes_right_septum_dev_avg_pls_3SD);
for i = 1:length(I_mdslpes_right_septum_sorted)
    if I_mdslpes_right_septum_dev_avg_pls_3SD(i) == min(I_mdslpes_right_septum_dev_avg_pls_3SD)
        I_mdslpes_right_septum_dev_avg_pls_3SD_index = i;
    end
end
I_mdslpes_right_septum_tot_number_of_indices = length(I_mdslpes_right_septum_sorted);
%This is the approximation of the tot "area under curve." 
I_mdslpes_right_septum_quotient_within_1SD = (I_mdslpes_right_septum_dev_avg_pls_1SD_index - I_mdslpes_right_septumdev_avg_min_1SD_index)/I_mdslpes_right_septum_tot_number_of_indices;
I_mdslpes_right_septum_quotient_within_2SD = (I_mdslpes_right_septum_dev_avg_pls_2SD_index - I_mdslpes_dev_avg_min_2SD_index)/I_mdslpes_right_septum_tot_number_of_indices;
I_mdslpes_right_septum_quotient_within_3SD = (I_mdslpes_right_septum_dev_avg_pls_3SD_index - I_mdslpes_right_septum_dev_avg_min_3SD_index)/I_mdslpes_right_septum_tot_number_of_indices;
if (I_mdslpes_right_septum_quotient_within_1SD > 0.66 && I_mdslpes_right_septum_quotient_within_1SD < 0.70) && (I_mdslpes_right_septum_quotient_within_2SD > 0.93 && I_mdslpes_right_septum_quotient_within_2SD < 0.97) && (I_mdslpes_right_septum_quotient_within_3SD > 0.98 && I_mdslpes_right_septum_quotient_within_3SD < 1)
    isNormal(6,7) = 1;
end
if I_mdslpes_right_septum_quotient_within_1SD == 0
    isSpike(6,7) = 1;
end
I_mdslpes_right_atrial_appendage_sorted = sort(I_mdslpes_right_atrial_appendage);
I_mdslpes_right_atrial_appendage_dev_avg_min_1SD = mdslpe(length(I_mdslpes_right_atrial_appendage_sorted),1);
for i = 1:length(I_mdslpes_right_atrial_appendage_sorted)
    I_mdslpes_right_atrial_appendage_dev_avg_min_1SD(i) = I_mdslpes_right_atrial_appendage_sorted(i) - (I_mdslpes_right_atrial_appendage_average-I_mdslpes_right_atrial_appendage_SD);
end
I_mdslpes_right_atrial_appendage_dev_avg_min_1SD = abs(I_mdslpes_right_atrial_appendage_dev_avg_min_1SD);
for i = 1:length(I_mdslpes_right_atrial_appendage_sorted)
    if I_mdslpes_right_atrial_appendage_dev_avg_min_1SD(i) == min(I_mdslpes_right_atrial_appendage_dev_avg_min_1SD)
        I_mdslpes_right_atrial_appendage_dev_avg_min_1SD_index = i;
    end
end
I_mdslpes_right_atrial_appendage_dev_avg_pls_1SD = mdslpe(length(I_mdslpes_right_atrial_appendage_sorted),1);
for i = 1:length(I_mdslpes_right_atrial_appendage_sorted)
    I_mdslpes_right_atrial_appendage_dev_avg_pls_1SD(i) = I_mdslpes_right_atrial_appendage_sorted(i) - (I_mdslpes_right_atrial_appendage_average+I_mdslpes_right_atrial_appendage_SD);
end
I_mdslpes_right_atrial_appendage_dev_avg_pls_1SD = abs(I_mdslpes_right_atrial_appendage_dev_avg_pls_1SD);
for i = 1:length(I_mdslpes_right_atrial_appendage_sorted)
    if I_mdslpes_right_atrial_appendage_dev_avg_pls_1SD(i) == min(I_mdslpes_right_atrial_appendage_dev_avg_pls_1SD)
        I_mdslpes_right_atrial_appendage_dev_avg_pls_1SD_index = i;
    end
end
I_mdslpes_right_atrial_appendage_dev_avg_min_2SD = mdslpe(length(I_mdslpes_right_atrial_appendage_sorted),1);
for i = 1:length(I_mdslpes_right_atrial_appendage_sorted)
    I_mdslpes_right_atrial_appendage_dev_avg_min_2SD(i) = I_mdslpes_right_atrial_appendage_sorted(i) - (I_mdslpes_right_atrial_appendage_average-2*I_mdslpes_right_atrial_appendage_SD);
end
I_mdslpes_right_atrial_appendage_dev_avg_min_2SD = abs(I_mdslpes_right_atrial_appendage_dev_avg_min_2SD);
for i = 1:length(I_mdslpes_right_atrial_appendage_sorted)
    if I_mdslpes_right_atrial_appendage_dev_avg_min_2SD(i) == min(I_mdslpes_right_atrial_appendage_dev_avg_min_2SD)
        I_mdslpes_dev_avg_min_2SD_index = i;
    end
end
I_mdslpes_right_atrial_appendage_dev_avg_pls_2SD = mdslpe(length(I_mdslpes_right_atrial_appendage_sorted),1);
for i = 1:length(I_mdslpes_right_atrial_appendage_sorted)
    I_mdslpes_right_atrial_appendage_dev_avg_pls_2SD(i) = I_mdslpes_right_atrial_appendage_sorted(i) - (I_mdslpes_right_atrial_appendage_average+2*I_mdslpes_right_atrial_appendage_SD);
end
I_mdslpes_right_atrial_appendage_dev_avg_pls_2SD = abs(I_mdslpes_right_atrial_appendage_dev_avg_pls_2SD);
for i = 1:length(I_mdslpes_right_atrial_appendage_sorted)
    if I_mdslpes_right_atrial_appendage_dev_avg_pls_2SD(i) == min(I_mdslpes_right_atrial_appendage_dev_avg_pls_2SD)
        I_mdslpes_right_atrial_appendage_dev_avg_pls_2SD_index = i;
    end
end
I_mdslpes_right_atrial_appendage_dev_avg_min_3SD = mdslpe(length(I_mdslpes_right_atrial_appendage_sorted),1);
for i = 1:length(I_mdslpes_right_atrial_appendage_sorted)
    I_mdslpes_right_atrial_appendage_dev_avg_min_3SD(i) = I_mdslpes_right_atrial_appendage_sorted(i) - (I_mdslpes_right_atrial_appendage_average-3*I_mdslpes_right_atrial_appendage_SD);
end
I_mdslpes_right_atrial_appendage_dev_avg_min_3SD = abs(I_mdslpes_right_atrial_appendage_dev_avg_min_3SD);
for i = 1:length(I_mdslpes_right_atrial_appendage_sorted)
    if I_mdslpes_right_atrial_appendage_dev_avg_min_3SD(i) == min(I_mdslpes_right_atrial_appendage_dev_avg_min_3SD)
        I_mdslpes_right_atrial_appendage_dev_avg_min_3SD_index = i;
    end
end
I_mdslpes_right_atrial_appendage_dev_avg_pls_3SD = mdslpe(length(I_mdslpes_right_atrial_appendage_sorted),1);
for i = 1:length(I_mdslpes_right_atrial_appendage_sorted)
    I_mdslpes_right_atrial_appendage_dev_avg_pls_3SD(i) = I_mdslpes_right_atrial_appendage_sorted(i) - (I_mdslpes_right_atrial_appendage_average+3*I_mdslpes_right_atrial_appendage_SD);
end
I_mdslpes_right_atrial_appendage_dev_avg_pls_3SD = abs(I_mdslpes_right_atrial_appendage_dev_avg_pls_3SD);
for i = 1:length(I_mdslpes_right_atrial_appendage_sorted)
    if I_mdslpes_right_atrial_appendage_dev_avg_pls_3SD(i) == min(I_mdslpes_right_atrial_appendage_dev_avg_pls_3SD)
        I_mdslpes_right_atrial_appendage_dev_avg_pls_3SD_index = i;
    end
end
I_mdslpes_right_atrial_appendage_tot_number_of_indices = length(I_mdslpes_right_atrial_appendage_sorted);
%This is the approximation of the tot "area under curve." 
I_mdslpes_right_atrial_appendage_quotient_within_1SD = (I_mdslpes_right_atrial_appendage_dev_avg_pls_1SD_index - I_mdslpes_right_atrial_appendagedev_avg_min_1SD_index)/I_mdslpes_right_atrial_appendage_tot_number_of_indices;
I_mdslpes_right_atrial_appendage_quotient_within_2SD = (I_mdslpes_right_atrial_appendage_dev_avg_pls_2SD_index - I_mdslpes_dev_avg_min_2SD_index)/I_mdslpes_right_atrial_appendage_tot_number_of_indices;
I_mdslpes_right_atrial_appendage_quotient_within_3SD = (I_mdslpes_right_atrial_appendage_dev_avg_pls_3SD_index - I_mdslpes_right_atrial_appendage_dev_avg_min_3SD_index)/I_mdslpes_right_atrial_appendage_tot_number_of_indices;
if (I_mdslpes_right_atrial_appendage_quotient_within_1SD > 0.66 && I_mdslpes_right_atrial_appendage_quotient_within_1SD < 0.70) && (I_mdslpes_right_atrial_appendage_quotient_within_2SD > 0.93 && I_mdslpes_right_atrial_appendage_quotient_within_2SD < 0.97) && (I_mdslpes_right_atrial_appendage_quotient_within_3SD > 0.98 && I_mdslpes_right_atrial_appendage_quotient_within_3SD < 1)
    isNormal(6,8) = 1;
end
if I_mdslpes_right_atrial_appendage_quotient_within_1SD == 0
    isSpike(6,8) = 1;
end
I_mdslpes_pulmonary_veins_sorted = sort(I_mdslpes_pulmonary_veins);
I_mdslpes_pulmonary_veins_dev_avg_min_1SD = mdslpe(length(I_mdslpes_pulmonary_veins_sorted),1);
for i = 1:length(I_mdslpes_pulmonary_veins_sorted)
    I_mdslpes_pulmonary_veins_dev_avg_min_1SD(i) = I_mdslpes_pulmonary_veins_sorted(i) - (I_mdslpes_pulmonary_veins_average-I_mdslpes_pulmonary_veins_SD);
end
I_mdslpes_pulmonary_veins_dev_avg_min_1SD = abs(I_mdslpes_pulmonary_veins_dev_avg_min_1SD);
for i = 1:length(I_mdslpes_pulmonary_veins_sorted)
    if I_mdslpes_pulmonary_veins_dev_avg_min_1SD(i) == min(I_mdslpes_pulmonary_veins_dev_avg_min_1SD)
        I_mdslpes_pulmonary_veins_dev_avg_min_1SD_index = i;
    end
end
I_mdslpes_pulmonary_veins_dev_avg_pls_1SD = mdslpe(length(I_mdslpes_pulmonary_veins_sorted),1);
for i = 1:length(I_mdslpes_pulmonary_veins_sorted)
    I_mdslpes_pulmonary_veins_dev_avg_pls_1SD(i) = I_mdslpes_pulmonary_veins_sorted(i) - (I_mdslpes_pulmonary_veins_average+I_mdslpes_pulmonary_veins_SD);
end
I_mdslpes_pulmonary_veins_dev_avg_pls_1SD = abs(I_mdslpes_pulmonary_veins_dev_avg_pls_1SD);
for i = 1:length(I_mdslpes_pulmonary_veins_sorted)
    if I_mdslpes_pulmonary_veins_dev_avg_pls_1SD(i) == min(I_mdslpes_pulmonary_veins_dev_avg_pls_1SD)
        I_mdslpes_pulmonary_veins_dev_avg_pls_1SD_index = i;
    end
end
I_mdslpes_pulmonary_veins_dev_avg_min_2SD = mdslpe(length(I_mdslpes_pulmonary_veins_sorted),1);
for i = 1:length(I_mdslpes_pulmonary_veins_sorted)
    I_mdslpes_pulmonary_veins_dev_avg_min_2SD(i) = I_mdslpes_pulmonary_veins_sorted(i) - (I_mdslpes_pulmonary_veins_average-2*I_mdslpes_pulmonary_veins_SD);
end
I_mdslpes_pulmonary_veins_dev_avg_min_2SD = abs(I_mdslpes_pulmonary_veins_dev_avg_min_2SD);
for i = 1:length(I_mdslpes_pulmonary_veins_sorted)
    if I_mdslpes_pulmonary_veins_dev_avg_min_2SD(i) == min(I_mdslpes_pulmonary_veins_dev_avg_min_2SD)
        I_mdslpes_dev_avg_min_2SD_index = i;
    end
end
I_mdslpes_pulmonary_veins_dev_avg_pls_2SD = mdslpe(length(I_mdslpes_pulmonary_veins_sorted),1);
for i = 1:length(I_mdslpes_pulmonary_veins_sorted)
    I_mdslpes_pulmonary_veins_dev_avg_pls_2SD(i) = I_mdslpes_pulmonary_veins_sorted(i) - (I_mdslpes_pulmonary_veins_average+2*I_mdslpes_pulmonary_veins_SD);
end
I_mdslpes_pulmonary_veins_dev_avg_pls_2SD = abs(I_mdslpes_pulmonary_veins_dev_avg_pls_2SD);
for i = 1:length(I_mdslpes_pulmonary_veins_sorted)
    if I_mdslpes_pulmonary_veins_dev_avg_pls_2SD(i) == min(I_mdslpes_pulmonary_veins_dev_avg_pls_2SD)
        I_mdslpes_pulmonary_veins_dev_avg_pls_2SD_index = i;
    end
end
I_mdslpes_pulmonary_veins_dev_avg_min_3SD = mdslpe(length(I_mdslpes_pulmonary_veins_sorted),1);
for i = 1:length(I_mdslpes_pulmonary_veins_sorted)
    I_mdslpes_pulmonary_veins_dev_avg_min_3SD(i) = I_mdslpes_pulmonary_veins_sorted(i) - (I_mdslpes_pulmonary_veins_average-3*I_mdslpes_pulmonary_veins_SD);
end
I_mdslpes_pulmonary_veins_dev_avg_min_3SD = abs(I_mdslpes_pulmonary_veins_dev_avg_min_3SD);
for i = 1:length(I_mdslpes_pulmonary_veins_sorted)
    if I_mdslpes_pulmonary_veins_dev_avg_min_3SD(i) == min(I_mdslpes_pulmonary_veins_dev_avg_min_3SD)
        I_mdslpes_pulmonary_veins_dev_avg_min_3SD_index = i;
    end
end
I_mdslpes_pulmonary_veins_dev_avg_pls_3SD = mdslpe(length(I_mdslpes_pulmonary_veins_sorted),1);
for i = 1:length(I_mdslpes_pulmonary_veins_sorted)
    I_mdslpes_pulmonary_veins_dev_avg_pls_3SD(i) = I_mdslpes_pulmonary_veins_sorted(i) - (I_mdslpes_pulmonary_veins_average+3*I_mdslpes_pulmonary_veins_SD);
end
I_mdslpes_pulmonary_veins_dev_avg_pls_3SD = abs(I_mdslpes_pulmonary_veins_dev_avg_pls_3SD);
for i = 1:length(I_mdslpes_pulmonary_veins_sorted)
    if I_mdslpes_pulmonary_veins_dev_avg_pls_3SD(i) == min(I_mdslpes_pulmonary_veins_dev_avg_pls_3SD)
        I_mdslpes_pulmonary_veins_dev_avg_pls_3SD_index = i;
    end
end
I_mdslpes_pulmonary_veins_tot_number_of_indices = length(I_mdslpes_pulmonary_veins_sorted);
%This is the approximation of the tot "area under curve." 
I_mdslpes_pulmonary_veins_quotient_within_1SD = (I_mdslpes_pulmonary_veins_dev_avg_pls_1SD_index - I_mdslpes_pulmonary_veinsdev_avg_min_1SD_index)/I_mdslpes_pulmonary_veins_tot_number_of_indices;
I_mdslpes_pulmonary_veins_quotient_within_2SD = (I_mdslpes_pulmonary_veins_dev_avg_pls_2SD_index - I_mdslpes_dev_avg_min_2SD_index)/I_mdslpes_pulmonary_veins_tot_number_of_indices;
I_mdslpes_pulmonary_veins_quotient_within_3SD = (I_mdslpes_pulmonary_veins_dev_avg_pls_3SD_index - I_mdslpes_pulmonary_veins_dev_avg_min_3SD_index)/I_mdslpes_pulmonary_veins_tot_number_of_indices;
if (I_mdslpes_pulmonary_veins_quotient_within_1SD > 0.66 && I_mdslpes_pulmonary_veins_quotient_within_1SD < 0.70) && (I_mdslpes_pulmonary_veins_quotient_within_2SD > 0.93 && I_mdslpes_pulmonary_veins_quotient_within_2SD < 0.97) && (I_mdslpes_pulmonary_veins_quotient_within_3SD > 0.98 && I_mdslpes_pulmonary_veins_quotient_within_3SD < 1)
    isNormal(6,9) = 1;
end
if I_mdslpes_pulmonary_veins_quotient_within_1SD == 0
    isSpike(6,9) = 1;
end
I_mdslpes_mitral_annulus_sorted = sort(I_mdslpes_mitral_annulus);
I_mdslpes_mitral_annulus_dev_avg_min_1SD = mdslpe(length(I_mdslpes_mitral_annulus_sorted),1);
for i = 1:length(I_mdslpes_mitral_annulus_sorted)
    I_mdslpes_mitral_annulus_dev_avg_min_1SD(i) = I_mdslpes_mitral_annulus_sorted(i) - (I_mdslpes_mitral_annulus_average-I_mdslpes_mitral_annulus_SD);
end
I_mdslpes_mitral_annulus_dev_avg_min_1SD = abs(I_mdslpes_mitral_annulus_dev_avg_min_1SD);
for i = 1:length(I_mdslpes_mitral_annulus_sorted)
    if I_mdslpes_mitral_annulus_dev_avg_min_1SD(i) == min(I_mdslpes_mitral_annulus_dev_avg_min_1SD)
        I_mdslpes_mitral_annulus_dev_avg_min_1SD_index = i;
    end
end
I_mdslpes_mitral_annulus_dev_avg_pls_1SD = mdslpe(length(I_mdslpes_mitral_annulus_sorted),1);
for i = 1:length(I_mdslpes_mitral_annulus_sorted)
    I_mdslpes_mitral_annulus_dev_avg_pls_1SD(i) = I_mdslpes_mitral_annulus_sorted(i) - (I_mdslpes_mitral_annulus_average+I_mdslpes_mitral_annulus_SD);
end
I_mdslpes_mitral_annulus_dev_avg_pls_1SD = abs(I_mdslpes_mitral_annulus_dev_avg_pls_1SD);
for i = 1:length(I_mdslpes_mitral_annulus_sorted)
    if I_mdslpes_mitral_annulus_dev_avg_pls_1SD(i) == min(I_mdslpes_mitral_annulus_dev_avg_pls_1SD)
        I_mdslpes_mitral_annulus_dev_avg_pls_1SD_index = i;
    end
end
I_mdslpes_mitral_annulus_dev_avg_min_2SD = mdslpe(length(I_mdslpes_mitral_annulus_sorted),1);
for i = 1:length(I_mdslpes_mitral_annulus_sorted)
    I_mdslpes_mitral_annulus_dev_avg_min_2SD(i) = I_mdslpes_mitral_annulus_sorted(i) - (I_mdslpes_mitral_annulus_average-2*I_mdslpes_mitral_annulus_SD);
end
I_mdslpes_mitral_annulus_dev_avg_min_2SD = abs(I_mdslpes_mitral_annulus_dev_avg_min_2SD);
for i = 1:length(I_mdslpes_mitral_annulus_sorted)
    if I_mdslpes_mitral_annulus_dev_avg_min_2SD(i) == min(I_mdslpes_mitral_annulus_dev_avg_min_2SD)
        I_mdslpes_dev_avg_min_2SD_index = i;
    end
end
I_mdslpes_mitral_annulus_dev_avg_pls_2SD = mdslpe(length(I_mdslpes_mitral_annulus_sorted),1);
for i = 1:length(I_mdslpes_mitral_annulus_sorted)
    I_mdslpes_mitral_annulus_dev_avg_pls_2SD(i) = I_mdslpes_mitral_annulus_sorted(i) - (I_mdslpes_mitral_annulus_average+2*I_mdslpes_mitral_annulus_SD);
end
I_mdslpes_mitral_annulus_dev_avg_pls_2SD = abs(I_mdslpes_mitral_annulus_dev_avg_pls_2SD);
for i = 1:length(I_mdslpes_mitral_annulus_sorted)
    if I_mdslpes_mitral_annulus_dev_avg_pls_2SD(i) == min(I_mdslpes_mitral_annulus_dev_avg_pls_2SD)
        I_mdslpes_mitral_annulus_dev_avg_pls_2SD_index = i;
    end
end
I_mdslpes_mitral_annulus_dev_avg_min_3SD = mdslpe(length(I_mdslpes_mitral_annulus_sorted),1);
for i = 1:length(I_mdslpes_mitral_annulus_sorted)
    I_mdslpes_mitral_annulus_dev_avg_min_3SD(i) = I_mdslpes_mitral_annulus_sorted(i) - (I_mdslpes_mitral_annulus_average-3*I_mdslpes_mitral_annulus_SD);
end
I_mdslpes_mitral_annulus_dev_avg_min_3SD = abs(I_mdslpes_mitral_annulus_dev_avg_min_3SD);
for i = 1:length(I_mdslpes_mitral_annulus_sorted)
    if I_mdslpes_mitral_annulus_dev_avg_min_3SD(i) == min(I_mdslpes_mitral_annulus_dev_avg_min_3SD)
        I_mdslpes_mitral_annulus_dev_avg_min_3SD_index = i;
    end
end
I_mdslpes_mitral_annulus_dev_avg_pls_3SD = mdslpe(length(I_mdslpes_mitral_annulus_sorted),1);
for i = 1:length(I_mdslpes_mitral_annulus_sorted)
    I_mdslpes_mitral_annulus_dev_avg_pls_3SD(i) = I_mdslpes_mitral_annulus_sorted(i) - (I_mdslpes_mitral_annulus_average+3*I_mdslpes_mitral_annulus_SD);
end
I_mdslpes_mitral_annulus_dev_avg_pls_3SD = abs(I_mdslpes_mitral_annulus_dev_avg_pls_3SD);
for i = 1:length(I_mdslpes_mitral_annulus_sorted)
    if I_mdslpes_mitral_annulus_dev_avg_pls_3SD(i) == min(I_mdslpes_mitral_annulus_dev_avg_pls_3SD)
        I_mdslpes_mitral_annulus_dev_avg_pls_3SD_index = i;
    end
end
I_mdslpes_mitral_annulus_tot_number_of_indices = length(I_mdslpes_mitral_annulus_sorted);
%This is the approximation of the tot "area under curve." 
I_mdslpes_mitral_annulus_quotient_within_1SD = (I_mdslpes_mitral_annulus_dev_avg_pls_1SD_index - I_mdslpes_mitral_annulusdev_avg_min_1SD_index)/I_mdslpes_mitral_annulus_tot_number_of_indices;
I_mdslpes_mitral_annulus_quotient_within_2SD = (I_mdslpes_mitral_annulus_dev_avg_pls_2SD_index - I_mdslpes_dev_avg_min_2SD_index)/I_mdslpes_mitral_annulus_tot_number_of_indices;
I_mdslpes_mitral_annulus_quotient_within_3SD = (I_mdslpes_mitral_annulus_dev_avg_pls_3SD_index - I_mdslpes_mitral_annulus_dev_avg_min_3SD_index)/I_mdslpes_mitral_annulus_tot_number_of_indices;
if (I_mdslpes_mitral_annulus_quotient_within_1SD > 0.66 && I_mdslpes_mitral_annulus_quotient_within_1SD < 0.70) && (I_mdslpes_mitral_annulus_quotient_within_2SD > 0.93 && I_mdslpes_mitral_annulus_quotient_within_2SD < 0.97) && (I_mdslpes_mitral_annulus_quotient_within_3SD > 0.98 && I_mdslpes_mitral_annulus_quotient_within_3SD < 1)
    isNormal(6,10) = 1;
end
if I_mdslpes_mitral_annulus_quotient_within_1SD == 0
    isSpike(6,10) = 1;
end
I_mdslpes_CS_body_sorted = sort(I_mdslpes_CS_body);
I_mdslpes_CS_body_dev_avg_min_1SD = mdslpe(length(I_mdslpes_CS_body_sorted),1);
for i = 1:length(I_mdslpes_CS_body_sorted)
    I_mdslpes_CS_body_dev_avg_min_1SD(i) = I_mdslpes_CS_body_sorted(i) - (I_mdslpes_CS_body_average-I_mdslpes_CS_body_SD);
end
I_mdslpes_CS_body_dev_avg_min_1SD = abs(I_mdslpes_CS_body_dev_avg_min_1SD);
for i = 1:length(I_mdslpes_CS_body_sorted)
    if I_mdslpes_CS_body_dev_avg_min_1SD(i) == min(I_mdslpes_CS_body_dev_avg_min_1SD)
        I_mdslpes_CS_body_dev_avg_min_1SD_index = i;
    end
end
I_mdslpes_CS_body_dev_avg_pls_1SD = mdslpe(length(I_mdslpes_CS_body_sorted),1);
for i = 1:length(I_mdslpes_CS_body_sorted)
    I_mdslpes_CS_body_dev_avg_pls_1SD(i) = I_mdslpes_CS_body_sorted(i) - (I_mdslpes_CS_body_average+I_mdslpes_CS_body_SD);
end
I_mdslpes_CS_body_dev_avg_pls_1SD = abs(I_mdslpes_CS_body_dev_avg_pls_1SD);
for i = 1:length(I_mdslpes_CS_body_sorted)
    if I_mdslpes_CS_body_dev_avg_pls_1SD(i) == min(I_mdslpes_CS_body_dev_avg_pls_1SD)
        I_mdslpes_CS_body_dev_avg_pls_1SD_index = i;
    end
end
I_mdslpes_CS_body_dev_avg_min_2SD = mdslpe(length(I_mdslpes_CS_body_sorted),1);
for i = 1:length(I_mdslpes_CS_body_sorted)
    I_mdslpes_CS_body_dev_avg_min_2SD(i) = I_mdslpes_CS_body_sorted(i) - (I_mdslpes_CS_body_average-2*I_mdslpes_CS_body_SD);
end
I_mdslpes_CS_body_dev_avg_min_2SD = abs(I_mdslpes_CS_body_dev_avg_min_2SD);
for i = 1:length(I_mdslpes_CS_body_sorted)
    if I_mdslpes_CS_body_dev_avg_min_2SD(i) == min(I_mdslpes_CS_body_dev_avg_min_2SD)
        I_mdslpes_dev_avg_min_2SD_index = i;
    end
end
I_mdslpes_CS_body_dev_avg_pls_2SD = mdslpe(length(I_mdslpes_CS_body_sorted),1);
for i = 1:length(I_mdslpes_CS_body_sorted)
    I_mdslpes_CS_body_dev_avg_pls_2SD(i) = I_mdslpes_CS_body_sorted(i) - (I_mdslpes_CS_body_average+2*I_mdslpes_CS_body_SD);
end
I_mdslpes_CS_body_dev_avg_pls_2SD = abs(I_mdslpes_CS_body_dev_avg_pls_2SD);
for i = 1:length(I_mdslpes_CS_body_sorted)
    if I_mdslpes_CS_body_dev_avg_pls_2SD(i) == min(I_mdslpes_CS_body_dev_avg_pls_2SD)
        I_mdslpes_CS_body_dev_avg_pls_2SD_index = i;
    end
end
I_mdslpes_CS_body_dev_avg_min_3SD = mdslpe(length(I_mdslpes_CS_body_sorted),1);
for i = 1:length(I_mdslpes_CS_body_sorted)
    I_mdslpes_CS_body_dev_avg_min_3SD(i) = I_mdslpes_CS_body_sorted(i) - (I_mdslpes_CS_body_average-3*I_mdslpes_CS_body_SD);
end
I_mdslpes_CS_body_dev_avg_min_3SD = abs(I_mdslpes_CS_body_dev_avg_min_3SD);
for i = 1:length(I_mdslpes_CS_body_sorted)
    if I_mdslpes_CS_body_dev_avg_min_3SD(i) == min(I_mdslpes_CS_body_dev_avg_min_3SD)
        I_mdslpes_CS_body_dev_avg_min_3SD_index = i;
    end
end
I_mdslpes_CS_body_dev_avg_pls_3SD = mdslpe(length(I_mdslpes_CS_body_sorted),1);
for i = 1:length(I_mdslpes_CS_body_sorted)
    I_mdslpes_CS_body_dev_avg_pls_3SD(i) = I_mdslpes_CS_body_sorted(i) - (I_mdslpes_CS_body_average+3*I_mdslpes_CS_body_SD);
end
I_mdslpes_CS_body_dev_avg_pls_3SD = abs(I_mdslpes_CS_body_dev_avg_pls_3SD);
for i = 1:length(I_mdslpes_CS_body_sorted)
    if I_mdslpes_CS_body_dev_avg_pls_3SD(i) == min(I_mdslpes_CS_body_dev_avg_pls_3SD)
        I_mdslpes_CS_body_dev_avg_pls_3SD_index = i;
    end
end
I_mdslpes_CS_body_tot_number_of_indices = length(I_mdslpes_CS_body_sorted);
%This is the approximation of the tot "area under curve." 
I_mdslpes_CS_body_quotient_within_1SD = (I_mdslpes_CS_body_dev_avg_pls_1SD_index - I_mdslpes_CS_bodydev_avg_min_1SD_index)/I_mdslpes_CS_body_tot_number_of_indices;
I_mdslpes_CS_body_quotient_within_2SD = (I_mdslpes_CS_body_dev_avg_pls_2SD_index - I_mdslpes_dev_avg_min_2SD_index)/I_mdslpes_CS_body_tot_number_of_indices;
I_mdslpes_CS_body_quotient_within_3SD = (I_mdslpes_CS_body_dev_avg_pls_3SD_index - I_mdslpes_CS_body_dev_avg_min_3SD_index)/I_mdslpes_CS_body_tot_number_of_indices;
if (I_mdslpes_CS_body_quotient_within_1SD > 0.66 && I_mdslpes_CS_body_quotient_within_1SD < 0.70) && (I_mdslpes_CS_body_quotient_within_2SD > 0.93 && I_mdslpes_CS_body_quotient_within_2SD < 0.97) && (I_mdslpes_CS_body_quotient_within_3SD > 0.98 && I_mdslpes_CS_body_quotient_within_3SD < 1)
    isNormal(6,11) = 1;
end
if I_mdslpes_CS_body_quotient_within_1SD == 0
    isSpike(6,11) = 1;
end
I_mdslpes_left_septum_sorted = sort(I_mdslpes_left_septum);
I_mdslpes_left_septum_dev_avg_min_1SD = mdslpe(length(I_mdslpes_left_septum_sorted),1);
for i = 1:length(I_mdslpes_left_septum_sorted)
    I_mdslpes_left_septum_dev_avg_min_1SD(i) = I_mdslpes_left_septum_sorted(i) - (I_mdslpes_left_septum_average-I_mdslpes_left_septum_SD);
end
I_mdslpes_left_septum_dev_avg_min_1SD = abs(I_mdslpes_left_septum_dev_avg_min_1SD);
for i = 1:length(I_mdslpes_left_septum_sorted)
    if I_mdslpes_left_septum_dev_avg_min_1SD(i) == min(I_mdslpes_left_septum_dev_avg_min_1SD)
        I_mdslpes_left_septum_dev_avg_min_1SD_index = i;
    end
end
I_mdslpes_left_septum_dev_avg_pls_1SD = mdslpe(length(I_mdslpes_left_septum_sorted),1);
for i = 1:length(I_mdslpes_left_septum_sorted)
    I_mdslpes_left_septum_dev_avg_pls_1SD(i) = I_mdslpes_left_septum_sorted(i) - (I_mdslpes_left_septum_average+I_mdslpes_left_septum_SD);
end
I_mdslpes_left_septum_dev_avg_pls_1SD = abs(I_mdslpes_left_septum_dev_avg_pls_1SD);
for i = 1:length(I_mdslpes_left_septum_sorted)
    if I_mdslpes_left_septum_dev_avg_pls_1SD(i) == min(I_mdslpes_left_septum_dev_avg_pls_1SD)
        I_mdslpes_left_septum_dev_avg_pls_1SD_index = i;
    end
end
I_mdslpes_left_septum_dev_avg_min_2SD = mdslpe(length(I_mdslpes_left_septum_sorted),1);
for i = 1:length(I_mdslpes_left_septum_sorted)
    I_mdslpes_left_septum_dev_avg_min_2SD(i) = I_mdslpes_left_septum_sorted(i) - (I_mdslpes_left_septum_average-2*I_mdslpes_left_septum_SD);
end
I_mdslpes_left_septum_dev_avg_min_2SD = abs(I_mdslpes_left_septum_dev_avg_min_2SD);
for i = 1:length(I_mdslpes_left_septum_sorted)
    if I_mdslpes_left_septum_dev_avg_min_2SD(i) == min(I_mdslpes_left_septum_dev_avg_min_2SD)
        I_mdslpes_dev_avg_min_2SD_index = i;
    end
end
I_mdslpes_left_septum_dev_avg_pls_2SD = mdslpe(length(I_mdslpes_left_septum_sorted),1);
for i = 1:length(I_mdslpes_left_septum_sorted)
    I_mdslpes_left_septum_dev_avg_pls_2SD(i) = I_mdslpes_left_septum_sorted(i) - (I_mdslpes_left_septum_average+2*I_mdslpes_left_septum_SD);
end
I_mdslpes_left_septum_dev_avg_pls_2SD = abs(I_mdslpes_left_septum_dev_avg_pls_2SD);
for i = 1:length(I_mdslpes_left_septum_sorted)
    if I_mdslpes_left_septum_dev_avg_pls_2SD(i) == min(I_mdslpes_left_septum_dev_avg_pls_2SD)
        I_mdslpes_left_septum_dev_avg_pls_2SD_index = i;
    end
end
I_mdslpes_left_septum_dev_avg_min_3SD = mdslpe(length(I_mdslpes_left_septum_sorted),1);
for i = 1:length(I_mdslpes_left_septum_sorted)
    I_mdslpes_left_septum_dev_avg_min_3SD(i) = I_mdslpes_left_septum_sorted(i) - (I_mdslpes_left_septum_average-3*I_mdslpes_left_septum_SD);
end
I_mdslpes_left_septum_dev_avg_min_3SD = abs(I_mdslpes_left_septum_dev_avg_min_3SD);
for i = 1:length(I_mdslpes_left_septum_sorted)
    if I_mdslpes_left_septum_dev_avg_min_3SD(i) == min(I_mdslpes_left_septum_dev_avg_min_3SD)
        I_mdslpes_left_septum_dev_avg_min_3SD_index = i;
    end
end
I_mdslpes_left_septum_dev_avg_pls_3SD = mdslpe(length(I_mdslpes_left_septum_sorted),1);
for i = 1:length(I_mdslpes_left_septum_sorted)
    I_mdslpes_left_septum_dev_avg_pls_3SD(i) = I_mdslpes_left_septum_sorted(i) - (I_mdslpes_left_septum_average+3*I_mdslpes_left_septum_SD);
end
I_mdslpes_left_septum_dev_avg_pls_3SD = abs(I_mdslpes_left_septum_dev_avg_pls_3SD);
for i = 1:length(I_mdslpes_left_septum_sorted)
    if I_mdslpes_left_septum_dev_avg_pls_3SD(i) == min(I_mdslpes_left_septum_dev_avg_pls_3SD)
        I_mdslpes_left_septum_dev_avg_pls_3SD_index = i;
    end
end
I_mdslpes_left_septum_tot_number_of_indices = length(I_mdslpes_left_septum_sorted);
%This is the approximation of the tot "area under curve." 
I_mdslpes_left_septum_quotient_within_1SD = (I_mdslpes_left_septum_dev_avg_pls_1SD_index - I_mdslpes_left_septumdev_avg_min_1SD_index)/I_mdslpes_left_septum_tot_number_of_indices;
I_mdslpes_left_septum_quotient_within_2SD = (I_mdslpes_left_septum_dev_avg_pls_2SD_index - I_mdslpes_dev_avg_min_2SD_index)/I_mdslpes_left_septum_tot_number_of_indices;
I_mdslpes_left_septum_quotient_within_3SD = (I_mdslpes_left_septum_dev_avg_pls_3SD_index - I_mdslpes_left_septum_dev_avg_min_3SD_index)/I_mdslpes_left_septum_tot_number_of_indices;
if (I_mdslpes_left_septum_quotient_within_1SD > 0.66 && I_mdslpes_left_septum_quotient_within_1SD < 0.70) && (I_mdslpes_left_septum_quotient_within_2SD > 0.93 && I_mdslpes_left_septum_quotient_within_2SD < 0.97) && (I_mdslpes_left_septum_quotient_within_3SD > 0.98 && I_mdslpes_left_septum_quotient_within_3SD < 1)
    isNormal(6,12) = 1;
end
if I_mdslpes_left_septum_quotient_within_1SD == 0
    isSpike(6,12) = 1;
end
I_mdslpes_left_atrial_appendage_sorted = sort(I_mdslpes_left_atrial_appendage);
I_mdslpes_left_atrial_appendage_dev_avg_min_1SD = mdslpe(length(I_mdslpes_left_atrial_appendage_sorted),1);
for i = 1:length(I_mdslpes_left_atrial_appendage_sorted)
    I_mdslpes_left_atrial_appendage_dev_avg_min_1SD(i) = I_mdslpes_left_atrial_appendage_sorted(i) - (I_mdslpes_left_atrial_appendage_average-I_mdslpes_left_atrial_appendage_SD);
end
I_mdslpes_left_atrial_appendage_dev_avg_min_1SD = abs(I_mdslpes_left_atrial_appendage_dev_avg_min_1SD);
for i = 1:length(I_mdslpes_left_atrial_appendage_sorted)
    if I_mdslpes_left_atrial_appendage_dev_avg_min_1SD(i) == min(I_mdslpes_left_atrial_appendage_dev_avg_min_1SD)
        I_mdslpes_left_atrial_appendage_dev_avg_min_1SD_index = i;
    end
end
I_mdslpes_left_atrial_appendage_dev_avg_pls_1SD = mdslpe(length(I_mdslpes_left_atrial_appendage_sorted),1);
for i = 1:length(I_mdslpes_left_atrial_appendage_sorted)
    I_mdslpes_left_atrial_appendage_dev_avg_pls_1SD(i) = I_mdslpes_left_atrial_appendage_sorted(i) - (I_mdslpes_left_atrial_appendage_average+I_mdslpes_left_atrial_appendage_SD);
end
I_mdslpes_left_atrial_appendage_dev_avg_pls_1SD = abs(I_mdslpes_left_atrial_appendage_dev_avg_pls_1SD);
for i = 1:length(I_mdslpes_left_atrial_appendage_sorted)
    if I_mdslpes_left_atrial_appendage_dev_avg_pls_1SD(i) == min(I_mdslpes_left_atrial_appendage_dev_avg_pls_1SD)
        I_mdslpes_left_atrial_appendage_dev_avg_pls_1SD_index = i;
    end
end
I_mdslpes_left_atrial_appendage_dev_avg_min_2SD = mdslpe(length(I_mdslpes_left_atrial_appendage_sorted),1);
for i = 1:length(I_mdslpes_left_atrial_appendage_sorted)
    I_mdslpes_left_atrial_appendage_dev_avg_min_2SD(i) = I_mdslpes_left_atrial_appendage_sorted(i) - (I_mdslpes_left_atrial_appendage_average-2*I_mdslpes_left_atrial_appendage_SD);
end
I_mdslpes_left_atrial_appendage_dev_avg_min_2SD = abs(I_mdslpes_left_atrial_appendage_dev_avg_min_2SD);
for i = 1:length(I_mdslpes_left_atrial_appendage_sorted)
    if I_mdslpes_left_atrial_appendage_dev_avg_min_2SD(i) == min(I_mdslpes_left_atrial_appendage_dev_avg_min_2SD)
        I_mdslpes_dev_avg_min_2SD_index = i;
    end
end
I_mdslpes_left_atrial_appendage_dev_avg_pls_2SD = mdslpe(length(I_mdslpes_left_atrial_appendage_sorted),1);
for i = 1:length(I_mdslpes_left_atrial_appendage_sorted)
    I_mdslpes_left_atrial_appendage_dev_avg_pls_2SD(i) = I_mdslpes_left_atrial_appendage_sorted(i) - (I_mdslpes_left_atrial_appendage_average+2*I_mdslpes_left_atrial_appendage_SD);
end
I_mdslpes_left_atrial_appendage_dev_avg_pls_2SD = abs(I_mdslpes_left_atrial_appendage_dev_avg_pls_2SD);
for i = 1:length(I_mdslpes_left_atrial_appendage_sorted)
    if I_mdslpes_left_atrial_appendage_dev_avg_pls_2SD(i) == min(I_mdslpes_left_atrial_appendage_dev_avg_pls_2SD)
        I_mdslpes_left_atrial_appendage_dev_avg_pls_2SD_index = i;
    end
end
I_mdslpes_left_atrial_appendage_dev_avg_min_3SD = mdslpe(length(I_mdslpes_left_atrial_appendage_sorted),1);
for i = 1:length(I_mdslpes_left_atrial_appendage_sorted)
    I_mdslpes_left_atrial_appendage_dev_avg_min_3SD(i) = I_mdslpes_left_atrial_appendage_sorted(i) - (I_mdslpes_left_atrial_appendage_average-3*I_mdslpes_left_atrial_appendage_SD);
end
I_mdslpes_left_atrial_appendage_dev_avg_min_3SD = abs(I_mdslpes_left_atrial_appendage_dev_avg_min_3SD);
for i = 1:length(I_mdslpes_left_atrial_appendage_sorted)
    if I_mdslpes_left_atrial_appendage_dev_avg_min_3SD(i) == min(I_mdslpes_left_atrial_appendage_dev_avg_min_3SD)
        I_mdslpes_left_atrial_appendage_dev_avg_min_3SD_index = i;
    end
end
I_mdslpes_left_atrial_appendage_dev_avg_pls_3SD = mdslpe(length(I_mdslpes_left_atrial_appendage_sorted),1);
for i = 1:length(I_mdslpes_left_atrial_appendage_sorted)
    I_mdslpes_left_atrial_appendage_dev_avg_pls_3SD(i) = I_mdslpes_left_atrial_appendage_sorted(i) - (I_mdslpes_left_atrial_appendage_average+3*I_mdslpes_left_atrial_appendage_SD);
end
I_mdslpes_left_atrial_appendage_dev_avg_pls_3SD = abs(I_mdslpes_left_atrial_appendage_dev_avg_pls_3SD);
for i = 1:length(I_mdslpes_left_atrial_appendage_sorted)
    if I_mdslpes_left_atrial_appendage_dev_avg_pls_3SD(i) == min(I_mdslpes_left_atrial_appendage_dev_avg_pls_3SD)
        I_mdslpes_left_atrial_appendage_dev_avg_pls_3SD_index = i;
    end
end
I_mdslpes_left_atrial_appendage_tot_number_of_indices = length(I_mdslpes_left_atrial_appendage_sorted);
%This is the approximation of the tot "area under curve." 
I_mdslpes_left_atrial_appendage_quotient_within_1SD = (I_mdslpes_left_atrial_appendage_dev_avg_pls_1SD_index - I_mdslpes_left_atrial_appendagedev_avg_min_1SD_index)/I_mdslpes_left_atrial_appendage_tot_number_of_indices;
I_mdslpes_left_atrial_appendage_quotient_within_2SD = (I_mdslpes_left_atrial_appendage_dev_avg_pls_2SD_index - I_mdslpes_dev_avg_min_2SD_index)/I_mdslpes_left_atrial_appendage_tot_number_of_indices;
I_mdslpes_left_atrial_appendage_quotient_within_3SD = (I_mdslpes_left_atrial_appendage_dev_avg_pls_3SD_index - I_mdslpes_left_atrial_appendage_dev_avg_min_3SD_index)/I_mdslpes_left_atrial_appendage_tot_number_of_indices;
if (I_mdslpes_left_atrial_appendage_quotient_within_1SD > 0.66 && I_mdslpes_left_atrial_appendage_quotient_within_1SD < 0.70) && (I_mdslpes_left_atrial_appendage_quotient_within_2SD > 0.93 && I_mdslpes_left_atrial_appendage_quotient_within_2SD < 0.97) && (I_mdslpes_left_atrial_appendage_quotient_within_3SD > 0.98 && I_mdslpes_left_atrial_appendage_quotient_within_3SD < 1)
    isNormal(6,13) = 1;
end
if I_mdslpes_left_atrial_appendage_quotient_within_1SD == 0
    isSpike(6,13) = 1;
end
I_p_wave_areas_1_SA_node_sorted = sort(I_p_wave_areas_1_SA_node);
I_p_wave_areas_1_SA_node_dev_avg_min_1SD = area_1(length(I_p_wave_areas_1_SA_node_sorted),1);
for i = 1:length(I_p_wave_areas_1_SA_node_sorted)
    I_p_wave_areas_1_SA_node_dev_avg_min_1SD(i) = I_p_wave_areas_1_SA_node_sorted(i) - (I_p_wave_areas_1_SA_node_average-I_p_wave_areas_1_SA_node_SD);
end
I_p_wave_areas_1_SA_node_dev_avg_min_1SD = abs(I_p_wave_areas_1_SA_node_dev_avg_min_1SD);
for i = 1:length(I_p_wave_areas_1_SA_node_sorted)
    if I_p_wave_areas_1_SA_node_dev_avg_min_1SD(i) == min(I_p_wave_areas_1_SA_node_dev_avg_min_1SD)
        I_p_wave_areas_1_SA_node_dev_avg_min_1SD_index = i;
    end
end
I_p_wave_areas_1_SA_node_dev_avg_pls_1SD = area_1(length(I_p_wave_areas_1_SA_node_sorted),1);
for i = 1:length(I_p_wave_areas_1_SA_node_sorted)
    I_p_wave_areas_1_SA_node_dev_avg_pls_1SD(i) = I_p_wave_areas_1_SA_node_sorted(i) - (I_p_wave_areas_1_SA_node_average+I_p_wave_areas_1_SA_node_SD);
end
I_p_wave_areas_1_SA_node_dev_avg_pls_1SD = abs(I_p_wave_areas_1_SA_node_dev_avg_pls_1SD);
for i = 1:length(I_p_wave_areas_1_SA_node_sorted)
    if I_p_wave_areas_1_SA_node_dev_avg_pls_1SD(i) == min(I_p_wave_areas_1_SA_node_dev_avg_pls_1SD)
        I_p_wave_areas_1_SA_node_dev_avg_pls_1SD_index = i;
    end
end
I_p_wave_areas_1_SA_node_dev_avg_min_2SD = area_1(length(I_p_wave_areas_1_SA_node_sorted),1);
for i = 1:length(I_p_wave_areas_1_SA_node_sorted)
    I_p_wave_areas_1_SA_node_dev_avg_min_2SD(i) = I_p_wave_areas_1_SA_node_sorted(i) - (I_p_wave_areas_1_SA_node_average-2*I_p_wave_areas_1_SA_node_SD);
end
I_p_wave_areas_1_SA_node_dev_avg_min_2SD = abs(I_p_wave_areas_1_SA_node_dev_avg_min_2SD);
for i = 1:length(I_p_wave_areas_1_SA_node_sorted)
    if I_p_wave_areas_1_SA_node_dev_avg_min_2SD(i) == min(I_p_wave_areas_1_SA_node_dev_avg_min_2SD)
        I_p_wave_areas_1_dev_avg_min_2SD_index = i;
    end
end
I_p_wave_areas_1_SA_node_dev_avg_pls_2SD = area_1(length(I_p_wave_areas_1_SA_node_sorted),1);
for i = 1:length(I_p_wave_areas_1_SA_node_sorted)
    I_p_wave_areas_1_SA_node_dev_avg_pls_2SD(i) = I_p_wave_areas_1_SA_node_sorted(i) - (I_p_wave_areas_1_SA_node_average+2*I_p_wave_areas_1_SA_node_SD);
end
I_p_wave_areas_1_SA_node_dev_avg_pls_2SD = abs(I_p_wave_areas_1_SA_node_dev_avg_pls_2SD);
for i = 1:length(I_p_wave_areas_1_SA_node_sorted)
    if I_p_wave_areas_1_SA_node_dev_avg_pls_2SD(i) == min(I_p_wave_areas_1_SA_node_dev_avg_pls_2SD)
        I_p_wave_areas_1_SA_node_dev_avg_pls_2SD_index = i;
    end
end
I_p_wave_areas_1_SA_node_dev_avg_min_3SD = area_1(length(I_p_wave_areas_1_SA_node_sorted),1);
for i = 1:length(I_p_wave_areas_1_SA_node_sorted)
    I_p_wave_areas_1_SA_node_dev_avg_min_3SD(i) = I_p_wave_areas_1_SA_node_sorted(i) - (I_p_wave_areas_1_SA_node_average-3*I_p_wave_areas_1_SA_node_SD);
end
I_p_wave_areas_1_SA_node_dev_avg_min_3SD = abs(I_p_wave_areas_1_SA_node_dev_avg_min_3SD);
for i = 1:length(I_p_wave_areas_1_SA_node_sorted)
    if I_p_wave_areas_1_SA_node_dev_avg_min_3SD(i) == min(I_p_wave_areas_1_SA_node_dev_avg_min_3SD)
        I_p_wave_areas_1_SA_node_dev_avg_min_3SD_index = i;
    end
end
I_p_wave_areas_1_SA_node_dev_avg_pls_3SD = area_1(length(I_p_wave_areas_1_SA_node_sorted),1);
for i = 1:length(I_p_wave_areas_1_SA_node_sorted)
    I_p_wave_areas_1_SA_node_dev_avg_pls_3SD(i) = I_p_wave_areas_1_SA_node_sorted(i) - (I_p_wave_areas_1_SA_node_average+3*I_p_wave_areas_1_SA_node_SD);
end
I_p_wave_areas_1_SA_node_dev_avg_pls_3SD = abs(I_p_wave_areas_1_SA_node_dev_avg_pls_3SD);
for i = 1:length(I_p_wave_areas_1_SA_node_sorted)
    if I_p_wave_areas_1_SA_node_dev_avg_pls_3SD(i) == min(I_p_wave_areas_1_SA_node_dev_avg_pls_3SD)
        I_p_wave_areas_1_SA_node_dev_avg_pls_3SD_index = i;
    end
end
I_p_wave_areas_1_SA_node_tot_number_of_indices = length(I_p_wave_areas_1_SA_node_sorted);
%This is the approximation of the tot "area under curve." 
I_p_wave_areas_1_SA_node_quotient_within_1SD = (I_p_wave_areas_1_SA_node_dev_avg_pls_1SD_index - I_p_wave_areas_1_SA_nodedev_avg_min_1SD_index)/I_p_wave_areas_1_SA_node_tot_number_of_indices;
I_p_wave_areas_1_SA_node_quotient_within_2SD = (I_p_wave_areas_1_SA_node_dev_avg_pls_2SD_index - I_p_wave_areas_1_dev_avg_min_2SD_index)/I_p_wave_areas_1_SA_node_tot_number_of_indices;
I_p_wave_areas_1_SA_node_quotient_within_3SD = (I_p_wave_areas_1_SA_node_dev_avg_pls_3SD_index - I_p_wave_areas_1_SA_node_dev_avg_min_3SD_index)/I_p_wave_areas_1_SA_node_tot_number_of_indices;
if ((I_p_wave_areas_1_SA_node_quotient_within_1SD > 0.66 && I_p_wave_areas_1_SA_node_quotient_within_1SD < 0.70) && (I_p_wave_areas_1_SA_node_quotient_within_2SD > 0.93 && I_p_wave_areas_1_SA_node_quotient_within_2SD < 0.97) && (I_p_wave_areas_1_SA_node_quotient_within_3SD > 0.98 && I_p_wave_areas_1_SA_node_quotient_within_3SD < 1)) 
    isNormal(7,1) = 1;
end
if I_p_wave_areas_1_SA_node_quotient_within_1SD == 0
    isSpike(7,1) = 1;
end
I_p_wave_areas_1_crista_terminalis_sorted = sort(I_p_wave_areas_1_crista_terminalis);
I_p_wave_areas_1_crista_terminalis_dev_avg_min_1SD = area_1(length(I_p_wave_areas_1_crista_terminalis_sorted),1);
for i = 1:length(I_p_wave_areas_1_crista_terminalis_sorted)
    I_p_wave_areas_1_crista_terminalis_dev_avg_min_1SD(i) = I_p_wave_areas_1_crista_terminalis_sorted(i) - (I_p_wave_areas_1_crista_terminalis_average-I_p_wave_areas_1_crista_terminalis_SD);
end
I_p_wave_areas_1_crista_terminalis_dev_avg_min_1SD = abs(I_p_wave_areas_1_crista_terminalis_dev_avg_min_1SD);
for i = 1:length(I_p_wave_areas_1_crista_terminalis_sorted)
    if I_p_wave_areas_1_crista_terminalis_dev_avg_min_1SD(i) == min(I_p_wave_areas_1_crista_terminalis_dev_avg_min_1SD)
        I_p_wave_areas_1_crista_terminalis_dev_avg_min_1SD_index = i;
    end
end
I_p_wave_areas_1_crista_terminalis_dev_avg_pls_1SD = area_1(length(I_p_wave_areas_1_crista_terminalis_sorted),1);
for i = 1:length(I_p_wave_areas_1_crista_terminalis_sorted)
    I_p_wave_areas_1_crista_terminalis_dev_avg_pls_1SD(i) = I_p_wave_areas_1_crista_terminalis_sorted(i) - (I_p_wave_areas_1_crista_terminalis_average+I_p_wave_areas_1_crista_terminalis_SD);
end
I_p_wave_areas_1_crista_terminalis_dev_avg_pls_1SD = abs(I_p_wave_areas_1_crista_terminalis_dev_avg_pls_1SD);
for i = 1:length(I_p_wave_areas_1_crista_terminalis_sorted)
    if I_p_wave_areas_1_crista_terminalis_dev_avg_pls_1SD(i) == min(I_p_wave_areas_1_crista_terminalis_dev_avg_pls_1SD)
        I_p_wave_areas_1_crista_terminalis_dev_avg_pls_1SD_index = i;
    end
end
I_p_wave_areas_1_crista_terminalis_dev_avg_min_2SD = area_1(length(I_p_wave_areas_1_crista_terminalis_sorted),1);
for i = 1:length(I_p_wave_areas_1_crista_terminalis_sorted)
    I_p_wave_areas_1_crista_terminalis_dev_avg_min_2SD(i) = I_p_wave_areas_1_crista_terminalis_sorted(i) - (I_p_wave_areas_1_crista_terminalis_average-2*I_p_wave_areas_1_crista_terminalis_SD);
end
I_p_wave_areas_1_crista_terminalis_dev_avg_min_2SD = abs(I_p_wave_areas_1_crista_terminalis_dev_avg_min_2SD);
for i = 1:length(I_p_wave_areas_1_crista_terminalis_sorted)
    if I_p_wave_areas_1_crista_terminalis_dev_avg_min_2SD(i) == min(I_p_wave_areas_1_crista_terminalis_dev_avg_min_2SD)
        I_p_wave_areas_1_dev_avg_min_2SD_index = i;
    end
end
I_p_wave_areas_1_crista_terminalis_dev_avg_pls_2SD = area_1(length(I_p_wave_areas_1_crista_terminalis_sorted),1);
for i = 1:length(I_p_wave_areas_1_crista_terminalis_sorted)
    I_p_wave_areas_1_crista_terminalis_dev_avg_pls_2SD(i) = I_p_wave_areas_1_crista_terminalis_sorted(i) - (I_p_wave_areas_1_crista_terminalis_average+2*I_p_wave_areas_1_crista_terminalis_SD);
end
I_p_wave_areas_1_crista_terminalis_dev_avg_pls_2SD = abs(I_p_wave_areas_1_crista_terminalis_dev_avg_pls_2SD);
for i = 1:length(I_p_wave_areas_1_crista_terminalis_sorted)
    if I_p_wave_areas_1_crista_terminalis_dev_avg_pls_2SD(i) == min(I_p_wave_areas_1_crista_terminalis_dev_avg_pls_2SD)
        I_p_wave_areas_1_crista_terminalis_dev_avg_pls_2SD_index = i;
    end
end
I_p_wave_areas_1_crista_terminalis_dev_avg_min_3SD = area_1(length(I_p_wave_areas_1_crista_terminalis_sorted),1);
for i = 1:length(I_p_wave_areas_1_crista_terminalis_sorted)
    I_p_wave_areas_1_crista_terminalis_dev_avg_min_3SD(i) = I_p_wave_areas_1_crista_terminalis_sorted(i) - (I_p_wave_areas_1_crista_terminalis_average-3*I_p_wave_areas_1_crista_terminalis_SD);
end
I_p_wave_areas_1_crista_terminalis_dev_avg_min_3SD = abs(I_p_wave_areas_1_crista_terminalis_dev_avg_min_3SD);
for i = 1:length(I_p_wave_areas_1_crista_terminalis_sorted)
    if I_p_wave_areas_1_crista_terminalis_dev_avg_min_3SD(i) == min(I_p_wave_areas_1_crista_terminalis_dev_avg_min_3SD)
        I_p_wave_areas_1_crista_terminalis_dev_avg_min_3SD_index = i;
    end
end
I_p_wave_areas_1_crista_terminalis_dev_avg_pls_3SD = area_1(length(I_p_wave_areas_1_crista_terminalis_sorted),1);
for i = 1:length(I_p_wave_areas_1_crista_terminalis_sorted)
    I_p_wave_areas_1_crista_terminalis_dev_avg_pls_3SD(i) = I_p_wave_areas_1_crista_terminalis_sorted(i) - (I_p_wave_areas_1_crista_terminalis_average+3*I_p_wave_areas_1_crista_terminalis_SD);
end
I_p_wave_areas_1_crista_terminalis_dev_avg_pls_3SD = abs(I_p_wave_areas_1_crista_terminalis_dev_avg_pls_3SD);
for i = 1:length(I_p_wave_areas_1_crista_terminalis_sorted)
    if I_p_wave_areas_1_crista_terminalis_dev_avg_pls_3SD(i) == min(I_p_wave_areas_1_crista_terminalis_dev_avg_pls_3SD)
        I_p_wave_areas_1_crista_terminalis_dev_avg_pls_3SD_index = i;
    end
end
I_p_wave_areas_1_crista_terminalis_tot_number_of_indices = length(I_p_wave_areas_1_crista_terminalis_sorted);
%This is the approximation of the tot "area under curve." 
I_p_wave_areas_1_crista_terminalis_quotient_within_1SD = (I_p_wave_areas_1_crista_terminalis_dev_avg_pls_1SD_index - I_p_wave_areas_1_crista_terminalisdev_avg_min_1SD_index)/I_p_wave_areas_1_crista_terminalis_tot_number_of_indices;
I_p_wave_areas_1_crista_terminalis_quotient_within_2SD = (I_p_wave_areas_1_crista_terminalis_dev_avg_pls_2SD_index - I_p_wave_areas_1_dev_avg_min_2SD_index)/I_p_wave_areas_1_crista_terminalis_tot_number_of_indices;
I_p_wave_areas_1_crista_terminalis_quotient_within_3SD = (I_p_wave_areas_1_crista_terminalis_dev_avg_pls_3SD_index - I_p_wave_areas_1_crista_terminalis_dev_avg_min_3SD_index)/I_p_wave_areas_1_crista_terminalis_tot_number_of_indices;
if (I_p_wave_areas_1_crista_terminalis_quotient_within_1SD > 0.66 && I_p_wave_areas_1_crista_terminalis_quotient_within_1SD < 0.70) && (I_p_wave_areas_1_crista_terminalis_quotient_within_2SD > 0.93 && I_p_wave_areas_1_crista_terminalis_quotient_within_2SD < 0.97) && (I_p_wave_areas_1_crista_terminalis_quotient_within_3SD > 0.98 && I_p_wave_areas_1_crista_terminalis_quotient_within_3SD < 1)
    isNormal(7,2) = 1;
end
if I_p_wave_areas_1_crista_terminalis_quotient_within_1SD == 0
    isSpike(7,2) = 1;
end
I_p_wave_areas_1_tricuspid_annulus_sorted = sort(I_p_wave_areas_1_tricuspid_annulus);
I_p_wave_areas_1_tricuspid_annulus_dev_avg_min_1SD = area_1(length(I_p_wave_areas_1_tricuspid_annulus_sorted),1);
for i = 1:length(I_p_wave_areas_1_tricuspid_annulus_sorted)
    I_p_wave_areas_1_tricuspid_annulus_dev_avg_min_1SD(i) = I_p_wave_areas_1_tricuspid_annulus_sorted(i) - (I_p_wave_areas_1_tricuspid_annulus_average-I_p_wave_areas_1_tricuspid_annulus_SD);
end
I_p_wave_areas_1_tricuspid_annulus_dev_avg_min_1SD = abs(I_p_wave_areas_1_tricuspid_annulus_dev_avg_min_1SD);
for i = 1:length(I_p_wave_areas_1_tricuspid_annulus_sorted)
    if I_p_wave_areas_1_tricuspid_annulus_dev_avg_min_1SD(i) == min(I_p_wave_areas_1_tricuspid_annulus_dev_avg_min_1SD)
        I_p_wave_areas_1_tricuspid_annulus_dev_avg_min_1SD_index = i;
    end
end
I_p_wave_areas_1_tricuspid_annulus_dev_avg_pls_1SD = area_1(length(I_p_wave_areas_1_tricuspid_annulus_sorted),1);
for i = 1:length(I_p_wave_areas_1_tricuspid_annulus_sorted)
    I_p_wave_areas_1_tricuspid_annulus_dev_avg_pls_1SD(i) = I_p_wave_areas_1_tricuspid_annulus_sorted(i) - (I_p_wave_areas_1_tricuspid_annulus_average+I_p_wave_areas_1_tricuspid_annulus_SD);
end
I_p_wave_areas_1_tricuspid_annulus_dev_avg_pls_1SD = abs(I_p_wave_areas_1_tricuspid_annulus_dev_avg_pls_1SD);
for i = 1:length(I_p_wave_areas_1_tricuspid_annulus_sorted)
    if I_p_wave_areas_1_tricuspid_annulus_dev_avg_pls_1SD(i) == min(I_p_wave_areas_1_tricuspid_annulus_dev_avg_pls_1SD)
        I_p_wave_areas_1_tricuspid_annulus_dev_avg_pls_1SD_index = i;
    end
end
I_p_wave_areas_1_tricuspid_annulus_dev_avg_min_2SD = area_1(length(I_p_wave_areas_1_tricuspid_annulus_sorted),1);
for i = 1:length(I_p_wave_areas_1_tricuspid_annulus_sorted)
    I_p_wave_areas_1_tricuspid_annulus_dev_avg_min_2SD(i) = I_p_wave_areas_1_tricuspid_annulus_sorted(i) - (I_p_wave_areas_1_tricuspid_annulus_average-2*I_p_wave_areas_1_tricuspid_annulus_SD);
end
I_p_wave_areas_1_tricuspid_annulus_dev_avg_min_2SD = abs(I_p_wave_areas_1_tricuspid_annulus_dev_avg_min_2SD);
for i = 1:length(I_p_wave_areas_1_tricuspid_annulus_sorted)
    if I_p_wave_areas_1_tricuspid_annulus_dev_avg_min_2SD(i) == min(I_p_wave_areas_1_tricuspid_annulus_dev_avg_min_2SD)
        I_p_wave_areas_1_dev_avg_min_2SD_index = i;
    end
end
I_p_wave_areas_1_tricuspid_annulus_dev_avg_pls_2SD = area_1(length(I_p_wave_areas_1_tricuspid_annulus_sorted),1);
for i = 1:length(I_p_wave_areas_1_tricuspid_annulus_sorted)
    I_p_wave_areas_1_tricuspid_annulus_dev_avg_pls_2SD(i) = I_p_wave_areas_1_tricuspid_annulus_sorted(i) - (I_p_wave_areas_1_tricuspid_annulus_average+2*I_p_wave_areas_1_tricuspid_annulus_SD);
end
I_p_wave_areas_1_tricuspid_annulus_dev_avg_pls_2SD = abs(I_p_wave_areas_1_tricuspid_annulus_dev_avg_pls_2SD);
for i = 1:length(I_p_wave_areas_1_tricuspid_annulus_sorted)
    if I_p_wave_areas_1_tricuspid_annulus_dev_avg_pls_2SD(i) == min(I_p_wave_areas_1_tricuspid_annulus_dev_avg_pls_2SD)
        I_p_wave_areas_1_tricuspid_annulus_dev_avg_pls_2SD_index = i;
    end
end
I_p_wave_areas_1_tricuspid_annulus_dev_avg_min_3SD = area_1(length(I_p_wave_areas_1_tricuspid_annulus_sorted),1);
for i = 1:length(I_p_wave_areas_1_tricuspid_annulus_sorted)
    I_p_wave_areas_1_tricuspid_annulus_dev_avg_min_3SD(i) = I_p_wave_areas_1_tricuspid_annulus_sorted(i) - (I_p_wave_areas_1_tricuspid_annulus_average-3*I_p_wave_areas_1_tricuspid_annulus_SD);
end
I_p_wave_areas_1_tricuspid_annulus_dev_avg_min_3SD = abs(I_p_wave_areas_1_tricuspid_annulus_dev_avg_min_3SD);
for i = 1:length(I_p_wave_areas_1_tricuspid_annulus_sorted)
    if I_p_wave_areas_1_tricuspid_annulus_dev_avg_min_3SD(i) == min(I_p_wave_areas_1_tricuspid_annulus_dev_avg_min_3SD)
        I_p_wave_areas_1_tricuspid_annulus_dev_avg_min_3SD_index = i;
    end
end
I_p_wave_areas_1_tricuspid_annulus_dev_avg_pls_3SD = area_1(length(I_p_wave_areas_1_tricuspid_annulus_sorted),1);
for i = 1:length(I_p_wave_areas_1_tricuspid_annulus_sorted)
    I_p_wave_areas_1_tricuspid_annulus_dev_avg_pls_3SD(i) = I_p_wave_areas_1_tricuspid_annulus_sorted(i) - (I_p_wave_areas_1_tricuspid_annulus_average+3*I_p_wave_areas_1_tricuspid_annulus_SD);
end
I_p_wave_areas_1_tricuspid_annulus_dev_avg_pls_3SD = abs(I_p_wave_areas_1_tricuspid_annulus_dev_avg_pls_3SD);
for i = 1:length(I_p_wave_areas_1_tricuspid_annulus_sorted)
    if I_p_wave_areas_1_tricuspid_annulus_dev_avg_pls_3SD(i) == min(I_p_wave_areas_1_tricuspid_annulus_dev_avg_pls_3SD)
        I_p_wave_areas_1_tricuspid_annulus_dev_avg_pls_3SD_index = i;
    end
end
I_p_wave_areas_1_tricuspid_annulus_tot_number_of_indices = length(I_p_wave_areas_1_tricuspid_annulus_sorted);
%This is the approximation of the tot "area under curve." 
I_p_wave_areas_1_tricuspid_annulus_quotient_within_1SD = (I_p_wave_areas_1_tricuspid_annulus_dev_avg_pls_1SD_index - I_p_wave_areas_1_tricuspid_annulusdev_avg_min_1SD_index)/I_p_wave_areas_1_tricuspid_annulus_tot_number_of_indices;
I_p_wave_areas_1_tricuspid_annulus_quotient_within_2SD = (I_p_wave_areas_1_tricuspid_annulus_dev_avg_pls_2SD_index - I_p_wave_areas_1_dev_avg_min_2SD_index)/I_p_wave_areas_1_tricuspid_annulus_tot_number_of_indices;
I_p_wave_areas_1_tricuspid_annulus_quotient_within_3SD = (I_p_wave_areas_1_tricuspid_annulus_dev_avg_pls_3SD_index - I_p_wave_areas_1_tricuspid_annulus_dev_avg_min_3SD_index)/I_p_wave_areas_1_tricuspid_annulus_tot_number_of_indices;
if (I_p_wave_areas_1_tricuspid_annulus_quotient_within_1SD > 0.66 && I_p_wave_areas_1_tricuspid_annulus_quotient_within_1SD < 0.70) && (I_p_wave_areas_1_tricuspid_annulus_quotient_within_2SD > 0.93 && I_p_wave_areas_1_tricuspid_annulus_quotient_within_2SD < 0.97) && (I_p_wave_areas_1_tricuspid_annulus_quotient_within_3SD > 0.98 && I_p_wave_areas_1_tricuspid_annulus_quotient_within_3SD < 1)
    isNormal(7,3) = 1;
end
if I_p_wave_areas_1_tricuspid_annulus_quotient_within_1SD == 0
    isSpike(7,3) = 1;
end
I_p_wave_areas_1_coronary_sinus_sorted = sort(I_p_wave_areas_1_coronary_sinus);
I_p_wave_areas_1_coronary_sinus_dev_avg_min_1SD = area_1(length(I_p_wave_areas_1_coronary_sinus_sorted),1);
for i = 1:length(I_p_wave_areas_1_coronary_sinus_sorted)
    I_p_wave_areas_1_coronary_sinus_dev_avg_min_1SD(i) = I_p_wave_areas_1_coronary_sinus_sorted(i) - (I_p_wave_areas_1_coronary_sinus_average-I_p_wave_areas_1_coronary_sinus_SD);
end
I_p_wave_areas_1_coronary_sinus_dev_avg_min_1SD = abs(I_p_wave_areas_1_coronary_sinus_dev_avg_min_1SD);
for i = 1:length(I_p_wave_areas_1_coronary_sinus_sorted)
    if I_p_wave_areas_1_coronary_sinus_dev_avg_min_1SD(i) == min(I_p_wave_areas_1_coronary_sinus_dev_avg_min_1SD)
        I_p_wave_areas_1_coronary_sinus_dev_avg_min_1SD_index = i;
    end
end
I_p_wave_areas_1_coronary_sinus_dev_avg_pls_1SD = area_1(length(I_p_wave_areas_1_coronary_sinus_sorted),1);
for i = 1:length(I_p_wave_areas_1_coronary_sinus_sorted)
    I_p_wave_areas_1_coronary_sinus_dev_avg_pls_1SD(i) = I_p_wave_areas_1_coronary_sinus_sorted(i) - (I_p_wave_areas_1_coronary_sinus_average+I_p_wave_areas_1_coronary_sinus_SD);
end
I_p_wave_areas_1_coronary_sinus_dev_avg_pls_1SD = abs(I_p_wave_areas_1_coronary_sinus_dev_avg_pls_1SD);
for i = 1:length(I_p_wave_areas_1_coronary_sinus_sorted)
    if I_p_wave_areas_1_coronary_sinus_dev_avg_pls_1SD(i) == min(I_p_wave_areas_1_coronary_sinus_dev_avg_pls_1SD)
        I_p_wave_areas_1_coronary_sinus_dev_avg_pls_1SD_index = i;
    end
end
I_p_wave_areas_1_coronary_sinus_dev_avg_min_2SD = area_1(length(I_p_wave_areas_1_coronary_sinus_sorted),1);
for i = 1:length(I_p_wave_areas_1_coronary_sinus_sorted)
    I_p_wave_areas_1_coronary_sinus_dev_avg_min_2SD(i) = I_p_wave_areas_1_coronary_sinus_sorted(i) - (I_p_wave_areas_1_coronary_sinus_average-2*I_p_wave_areas_1_coronary_sinus_SD);
end
I_p_wave_areas_1_coronary_sinus_dev_avg_min_2SD = abs(I_p_wave_areas_1_coronary_sinus_dev_avg_min_2SD);
for i = 1:length(I_p_wave_areas_1_coronary_sinus_sorted)
    if I_p_wave_areas_1_coronary_sinus_dev_avg_min_2SD(i) == min(I_p_wave_areas_1_coronary_sinus_dev_avg_min_2SD)
        I_p_wave_areas_1_dev_avg_min_2SD_index = i;
    end
end
I_p_wave_areas_1_coronary_sinus_dev_avg_pls_2SD = area_1(length(I_p_wave_areas_1_coronary_sinus_sorted),1);
for i = 1:length(I_p_wave_areas_1_coronary_sinus_sorted)
    I_p_wave_areas_1_coronary_sinus_dev_avg_pls_2SD(i) = I_p_wave_areas_1_coronary_sinus_sorted(i) - (I_p_wave_areas_1_coronary_sinus_average+2*I_p_wave_areas_1_coronary_sinus_SD);
end
I_p_wave_areas_1_coronary_sinus_dev_avg_pls_2SD = abs(I_p_wave_areas_1_coronary_sinus_dev_avg_pls_2SD);
for i = 1:length(I_p_wave_areas_1_coronary_sinus_sorted)
    if I_p_wave_areas_1_coronary_sinus_dev_avg_pls_2SD(i) == min(I_p_wave_areas_1_coronary_sinus_dev_avg_pls_2SD)
        I_p_wave_areas_1_coronary_sinus_dev_avg_pls_2SD_index = i;
    end
end
I_p_wave_areas_1_coronary_sinus_dev_avg_min_3SD = area_1(length(I_p_wave_areas_1_coronary_sinus_sorted),1);
for i = 1:length(I_p_wave_areas_1_coronary_sinus_sorted)
    I_p_wave_areas_1_coronary_sinus_dev_avg_min_3SD(i) = I_p_wave_areas_1_coronary_sinus_sorted(i) - (I_p_wave_areas_1_coronary_sinus_average-3*I_p_wave_areas_1_coronary_sinus_SD);
end
I_p_wave_areas_1_coronary_sinus_dev_avg_min_3SD = abs(I_p_wave_areas_1_coronary_sinus_dev_avg_min_3SD);
for i = 1:length(I_p_wave_areas_1_coronary_sinus_sorted)
    if I_p_wave_areas_1_coronary_sinus_dev_avg_min_3SD(i) == min(I_p_wave_areas_1_coronary_sinus_dev_avg_min_3SD)
        I_p_wave_areas_1_coronary_sinus_dev_avg_min_3SD_index = i;
    end
end
I_p_wave_areas_1_coronary_sinus_dev_avg_pls_3SD = area_1(length(I_p_wave_areas_1_coronary_sinus_sorted),1);
for i = 1:length(I_p_wave_areas_1_coronary_sinus_sorted)
    I_p_wave_areas_1_coronary_sinus_dev_avg_pls_3SD(i) = I_p_wave_areas_1_coronary_sinus_sorted(i) - (I_p_wave_areas_1_coronary_sinus_average+3*I_p_wave_areas_1_coronary_sinus_SD);
end
I_p_wave_areas_1_coronary_sinus_dev_avg_pls_3SD = abs(I_p_wave_areas_1_coronary_sinus_dev_avg_pls_3SD);
for i = 1:length(I_p_wave_areas_1_coronary_sinus_sorted)
    if I_p_wave_areas_1_coronary_sinus_dev_avg_pls_3SD(i) == min(I_p_wave_areas_1_coronary_sinus_dev_avg_pls_3SD)
        I_p_wave_areas_1_coronary_sinus_dev_avg_pls_3SD_index = i;
    end
end
I_p_wave_areas_1_coronary_sinus_tot_number_of_indices = length(I_p_wave_areas_1_coronary_sinus_sorted);
%This is the approximation of the tot "area under curve." 
I_p_wave_areas_1_coronary_sinus_quotient_within_1SD = (I_p_wave_areas_1_coronary_sinus_dev_avg_pls_1SD_index - I_p_wave_areas_1_coronary_sinusdev_avg_min_1SD_index)/I_p_wave_areas_1_coronary_sinus_tot_number_of_indices;
I_p_wave_areas_1_coronary_sinus_quotient_within_2SD = (I_p_wave_areas_1_coronary_sinus_dev_avg_pls_2SD_index - I_p_wave_areas_1_dev_avg_min_2SD_index)/I_p_wave_areas_1_coronary_sinus_tot_number_of_indices;
I_p_wave_areas_1_coronary_sinus_quotient_within_3SD = (I_p_wave_areas_1_coronary_sinus_dev_avg_pls_3SD_index - I_p_wave_areas_1_coronary_sinus_dev_avg_min_3SD_index)/I_p_wave_areas_1_coronary_sinus_tot_number_of_indices;
if (I_p_wave_areas_1_coronary_sinus_quotient_within_1SD > 0.66 && I_p_wave_areas_1_coronary_sinus_quotient_within_1SD < 0.70) && (I_p_wave_areas_1_coronary_sinus_quotient_within_2SD > 0.93 && I_p_wave_areas_1_coronary_sinus_quotient_within_2SD < 0.97) && (I_p_wave_areas_1_coronary_sinus_quotient_within_3SD > 0.98 && I_p_wave_areas_1_coronary_sinus_quotient_within_3SD < 1)
    isNormal(7,4) = 1;
end
if I_p_wave_areas_1_coronary_sinus_quotient_within_1SD == 0
    isSpike(7,4) = 1;
end
I_p_wave_areas_1_ostium_sorted = sort(I_p_wave_areas_1_ostium);
I_p_wave_areas_1_ostium_dev_avg_min_1SD = area_1(length(I_p_wave_areas_1_ostium_sorted),1);
for i = 1:length(I_p_wave_areas_1_ostium_sorted)
    I_p_wave_areas_1_ostium_dev_avg_min_1SD(i) = I_p_wave_areas_1_ostium_sorted(i) - (I_p_wave_areas_1_ostium_average-I_p_wave_areas_1_ostium_SD);
end
I_p_wave_areas_1_ostium_dev_avg_min_1SD = abs(I_p_wave_areas_1_ostium_dev_avg_min_1SD);
for i = 1:length(I_p_wave_areas_1_ostium_sorted)
    if I_p_wave_areas_1_ostium_dev_avg_min_1SD(i) == min(I_p_wave_areas_1_ostium_dev_avg_min_1SD)
        I_p_wave_areas_1_ostium_dev_avg_min_1SD_index = i;
    end
end
I_p_wave_areas_1_ostium_dev_avg_pls_1SD = area_1(length(I_p_wave_areas_1_ostium_sorted),1);
for i = 1:length(I_p_wave_areas_1_ostium_sorted)
    I_p_wave_areas_1_ostium_dev_avg_pls_1SD(i) = I_p_wave_areas_1_ostium_sorted(i) - (I_p_wave_areas_1_ostium_average+I_p_wave_areas_1_ostium_SD);
end
I_p_wave_areas_1_ostium_dev_avg_pls_1SD = abs(I_p_wave_areas_1_ostium_dev_avg_pls_1SD);
for i = 1:length(I_p_wave_areas_1_ostium_sorted)
    if I_p_wave_areas_1_ostium_dev_avg_pls_1SD(i) == min(I_p_wave_areas_1_ostium_dev_avg_pls_1SD)
        I_p_wave_areas_1_ostium_dev_avg_pls_1SD_index = i;
    end
end
I_p_wave_areas_1_ostium_dev_avg_min_2SD = area_1(length(I_p_wave_areas_1_ostium_sorted),1);
for i = 1:length(I_p_wave_areas_1_ostium_sorted)
    I_p_wave_areas_1_ostium_dev_avg_min_2SD(i) = I_p_wave_areas_1_ostium_sorted(i) - (I_p_wave_areas_1_ostium_average-2*I_p_wave_areas_1_ostium_SD);
end
I_p_wave_areas_1_ostium_dev_avg_min_2SD = abs(I_p_wave_areas_1_ostium_dev_avg_min_2SD);
for i = 1:length(I_p_wave_areas_1_ostium_sorted)
    if I_p_wave_areas_1_ostium_dev_avg_min_2SD(i) == min(I_p_wave_areas_1_ostium_dev_avg_min_2SD)
        I_p_wave_areas_1_dev_avg_min_2SD_index = i;
    end
end
I_p_wave_areas_1_ostium_dev_avg_pls_2SD = area_1(length(I_p_wave_areas_1_ostium_sorted),1);
for i = 1:length(I_p_wave_areas_1_ostium_sorted)
    I_p_wave_areas_1_ostium_dev_avg_pls_2SD(i) = I_p_wave_areas_1_ostium_sorted(i) - (I_p_wave_areas_1_ostium_average+2*I_p_wave_areas_1_ostium_SD);
end
I_p_wave_areas_1_ostium_dev_avg_pls_2SD = abs(I_p_wave_areas_1_ostium_dev_avg_pls_2SD);
for i = 1:length(I_p_wave_areas_1_ostium_sorted)
    if I_p_wave_areas_1_ostium_dev_avg_pls_2SD(i) == min(I_p_wave_areas_1_ostium_dev_avg_pls_2SD)
        I_p_wave_areas_1_ostium_dev_avg_pls_2SD_index = i;
    end
end
I_p_wave_areas_1_ostium_dev_avg_min_3SD = area_1(length(I_p_wave_areas_1_ostium_sorted),1);
for i = 1:length(I_p_wave_areas_1_ostium_sorted)
    I_p_wave_areas_1_ostium_dev_avg_min_3SD(i) = I_p_wave_areas_1_ostium_sorted(i) - (I_p_wave_areas_1_ostium_average-3*I_p_wave_areas_1_ostium_SD);
end
I_p_wave_areas_1_ostium_dev_avg_min_3SD = abs(I_p_wave_areas_1_ostium_dev_avg_min_3SD);
for i = 1:length(I_p_wave_areas_1_ostium_sorted)
    if I_p_wave_areas_1_ostium_dev_avg_min_3SD(i) == min(I_p_wave_areas_1_ostium_dev_avg_min_3SD)
        I_p_wave_areas_1_ostium_dev_avg_min_3SD_index = i;
    end
end
I_p_wave_areas_1_ostium_dev_avg_pls_3SD = area_1(length(I_p_wave_areas_1_ostium_sorted),1);
for i = 1:length(I_p_wave_areas_1_ostium_sorted)
    I_p_wave_areas_1_ostium_dev_avg_pls_3SD(i) = I_p_wave_areas_1_ostium_sorted(i) - (I_p_wave_areas_1_ostium_average+3*I_p_wave_areas_1_ostium_SD);
end
I_p_wave_areas_1_ostium_dev_avg_pls_3SD = abs(I_p_wave_areas_1_ostium_dev_avg_pls_3SD);
for i = 1:length(I_p_wave_areas_1_ostium_sorted)
    if I_p_wave_areas_1_ostium_dev_avg_pls_3SD(i) == min(I_p_wave_areas_1_ostium_dev_avg_pls_3SD)
        I_p_wave_areas_1_ostium_dev_avg_pls_3SD_index = i;
    end
end
I_p_wave_areas_1_ostium_tot_number_of_indices = length(I_p_wave_areas_1_ostium_sorted);
%This is the approximation of the tot "area under curve." 
I_p_wave_areas_1_ostium_quotient_within_1SD = (I_p_wave_areas_1_ostium_dev_avg_pls_1SD_index - I_p_wave_areas_1_ostiumdev_avg_min_1SD_index)/I_p_wave_areas_1_ostium_tot_number_of_indices;
I_p_wave_areas_1_ostium_quotient_within_2SD = (I_p_wave_areas_1_ostium_dev_avg_pls_2SD_index - I_p_wave_areas_1_dev_avg_min_2SD_index)/I_p_wave_areas_1_ostium_tot_number_of_indices;
I_p_wave_areas_1_ostium_quotient_within_3SD = (I_p_wave_areas_1_ostium_dev_avg_pls_3SD_index - I_p_wave_areas_1_ostium_dev_avg_min_3SD_index)/I_p_wave_areas_1_ostium_tot_number_of_indices;
if (I_p_wave_areas_1_ostium_quotient_within_1SD > 0.66 && I_p_wave_areas_1_ostium_quotient_within_1SD < 0.70) && (I_p_wave_areas_1_ostium_quotient_within_2SD > 0.93 && I_p_wave_areas_1_ostium_quotient_within_2SD < 0.97) && (I_p_wave_areas_1_ostium_quotient_within_3SD > 0.98 && I_p_wave_areas_1_ostium_quotient_within_3SD < 1)
    isNormal(7,5) = 1;
end
if I_p_wave_areas_1_ostium_quotient_within_1SD == 0
    isSpike(7,5) = 1;
end
I_p_wave_areas_1_perinodal_sorted = sort(I_p_wave_areas_1_perinodal);
I_p_wave_areas_1_perinodal_dev_avg_min_1SD = area_1(length(I_p_wave_areas_1_perinodal_sorted),1);
for i = 1:length(I_p_wave_areas_1_perinodal_sorted)
    I_p_wave_areas_1_perinodal_dev_avg_min_1SD(i) = I_p_wave_areas_1_perinodal_sorted(i) - (I_p_wave_areas_1_perinodal_average-I_p_wave_areas_1_perinodal_SD);
end
I_p_wave_areas_1_perinodal_dev_avg_min_1SD = abs(I_p_wave_areas_1_perinodal_dev_avg_min_1SD);
for i = 1:length(I_p_wave_areas_1_perinodal_sorted)
    if I_p_wave_areas_1_perinodal_dev_avg_min_1SD(i) == min(I_p_wave_areas_1_perinodal_dev_avg_min_1SD)
        I_p_wave_areas_1_perinodal_dev_avg_min_1SD_index = i;
    end
end
I_p_wave_areas_1_perinodal_dev_avg_pls_1SD = area_1(length(I_p_wave_areas_1_perinodal_sorted),1);
for i = 1:length(I_p_wave_areas_1_perinodal_sorted)
    I_p_wave_areas_1_perinodal_dev_avg_pls_1SD(i) = I_p_wave_areas_1_perinodal_sorted(i) - (I_p_wave_areas_1_perinodal_average+I_p_wave_areas_1_perinodal_SD);
end
I_p_wave_areas_1_perinodal_dev_avg_pls_1SD = abs(I_p_wave_areas_1_perinodal_dev_avg_pls_1SD);
for i = 1:length(I_p_wave_areas_1_perinodal_sorted)
    if I_p_wave_areas_1_perinodal_dev_avg_pls_1SD(i) == min(I_p_wave_areas_1_perinodal_dev_avg_pls_1SD)
        I_p_wave_areas_1_perinodal_dev_avg_pls_1SD_index = i;
    end
end
I_p_wave_areas_1_perinodal_dev_avg_min_2SD = area_1(length(I_p_wave_areas_1_perinodal_sorted),1);
for i = 1:length(I_p_wave_areas_1_perinodal_sorted)
    I_p_wave_areas_1_perinodal_dev_avg_min_2SD(i) = I_p_wave_areas_1_perinodal_sorted(i) - (I_p_wave_areas_1_perinodal_average-2*I_p_wave_areas_1_perinodal_SD);
end
I_p_wave_areas_1_perinodal_dev_avg_min_2SD = abs(I_p_wave_areas_1_perinodal_dev_avg_min_2SD);
for i = 1:length(I_p_wave_areas_1_perinodal_sorted)
    if I_p_wave_areas_1_perinodal_dev_avg_min_2SD(i) == min(I_p_wave_areas_1_perinodal_dev_avg_min_2SD)
        I_p_wave_areas_1_dev_avg_min_2SD_index = i;
    end
end
I_p_wave_areas_1_perinodal_dev_avg_pls_2SD = area_1(length(I_p_wave_areas_1_perinodal_sorted),1);
for i = 1:length(I_p_wave_areas_1_perinodal_sorted)
    I_p_wave_areas_1_perinodal_dev_avg_pls_2SD(i) = I_p_wave_areas_1_perinodal_sorted(i) - (I_p_wave_areas_1_perinodal_average+2*I_p_wave_areas_1_perinodal_SD);
end
I_p_wave_areas_1_perinodal_dev_avg_pls_2SD = abs(I_p_wave_areas_1_perinodal_dev_avg_pls_2SD);
for i = 1:length(I_p_wave_areas_1_perinodal_sorted)
    if I_p_wave_areas_1_perinodal_dev_avg_pls_2SD(i) == min(I_p_wave_areas_1_perinodal_dev_avg_pls_2SD)
        I_p_wave_areas_1_perinodal_dev_avg_pls_2SD_index = i;
    end
end
I_p_wave_areas_1_perinodal_dev_avg_min_3SD = area_1(length(I_p_wave_areas_1_perinodal_sorted),1);
for i = 1:length(I_p_wave_areas_1_perinodal_sorted)
    I_p_wave_areas_1_perinodal_dev_avg_min_3SD(i) = I_p_wave_areas_1_perinodal_sorted(i) - (I_p_wave_areas_1_perinodal_average-3*I_p_wave_areas_1_perinodal_SD);
end
I_p_wave_areas_1_perinodal_dev_avg_min_3SD = abs(I_p_wave_areas_1_perinodal_dev_avg_min_3SD);
for i = 1:length(I_p_wave_areas_1_perinodal_sorted)
    if I_p_wave_areas_1_perinodal_dev_avg_min_3SD(i) == min(I_p_wave_areas_1_perinodal_dev_avg_min_3SD)
        I_p_wave_areas_1_perinodal_dev_avg_min_3SD_index = i;
    end
end
I_p_wave_areas_1_perinodal_dev_avg_pls_3SD = area_1(length(I_p_wave_areas_1_perinodal_sorted),1);
for i = 1:length(I_p_wave_areas_1_perinodal_sorted)
    I_p_wave_areas_1_perinodal_dev_avg_pls_3SD(i) = I_p_wave_areas_1_perinodal_sorted(i) - (I_p_wave_areas_1_perinodal_average+3*I_p_wave_areas_1_perinodal_SD);
end
I_p_wave_areas_1_perinodal_dev_avg_pls_3SD = abs(I_p_wave_areas_1_perinodal_dev_avg_pls_3SD);
for i = 1:length(I_p_wave_areas_1_perinodal_sorted)
    if I_p_wave_areas_1_perinodal_dev_avg_pls_3SD(i) == min(I_p_wave_areas_1_perinodal_dev_avg_pls_3SD)
        I_p_wave_areas_1_perinodal_dev_avg_pls_3SD_index = i;
    end
end
I_p_wave_areas_1_perinodal_tot_number_of_indices = length(I_p_wave_areas_1_perinodal_sorted);
%This is the approximation of the tot "area under curve." 
I_p_wave_areas_1_perinodal_quotient_within_1SD = (I_p_wave_areas_1_perinodal_dev_avg_pls_1SD_index - I_p_wave_areas_1_perinodaldev_avg_min_1SD_index)/I_p_wave_areas_1_perinodal_tot_number_of_indices;
I_p_wave_areas_1_perinodal_quotient_within_2SD = (I_p_wave_areas_1_perinodal_dev_avg_pls_2SD_index - I_p_wave_areas_1_dev_avg_min_2SD_index)/I_p_wave_areas_1_perinodal_tot_number_of_indices;
I_p_wave_areas_1_perinodal_quotient_within_3SD = (I_p_wave_areas_1_perinodal_dev_avg_pls_3SD_index - I_p_wave_areas_1_perinodal_dev_avg_min_3SD_index)/I_p_wave_areas_1_perinodal_tot_number_of_indices;
if (I_p_wave_areas_1_perinodal_quotient_within_1SD > 0.66 && I_p_wave_areas_1_perinodal_quotient_within_1SD < 0.70) && (I_p_wave_areas_1_perinodal_quotient_within_2SD > 0.93 && I_p_wave_areas_1_perinodal_quotient_within_2SD < 0.97) && (I_p_wave_areas_1_perinodal_quotient_within_3SD > 0.98 && I_p_wave_areas_1_perinodal_quotient_within_3SD < 1)
    isNormal(7,6) = 1;
end
if I_p_wave_areas_1_perinodal_quotient_within_1SD == 0
    isSpike(7,6) = 1;
end
I_p_wave_areas_1_right_septum_sorted = sort(I_p_wave_areas_1_right_septum);
I_p_wave_areas_1_right_septum_dev_avg_min_1SD = area_1(length(I_p_wave_areas_1_right_septum_sorted),1);
for i = 1:length(I_p_wave_areas_1_right_septum_sorted)
    I_p_wave_areas_1_right_septum_dev_avg_min_1SD(i) = I_p_wave_areas_1_right_septum_sorted(i) - (I_p_wave_areas_1_right_septum_average-I_p_wave_areas_1_right_septum_SD);
end
I_p_wave_areas_1_right_septum_dev_avg_min_1SD = abs(I_p_wave_areas_1_right_septum_dev_avg_min_1SD);
for i = 1:length(I_p_wave_areas_1_right_septum_sorted)
    if I_p_wave_areas_1_right_septum_dev_avg_min_1SD(i) == min(I_p_wave_areas_1_right_septum_dev_avg_min_1SD)
        I_p_wave_areas_1_right_septum_dev_avg_min_1SD_index = i;
    end
end
I_p_wave_areas_1_right_septum_dev_avg_pls_1SD = area_1(length(I_p_wave_areas_1_right_septum_sorted),1);
for i = 1:length(I_p_wave_areas_1_right_septum_sorted)
    I_p_wave_areas_1_right_septum_dev_avg_pls_1SD(i) = I_p_wave_areas_1_right_septum_sorted(i) - (I_p_wave_areas_1_right_septum_average+I_p_wave_areas_1_right_septum_SD);
end
I_p_wave_areas_1_right_septum_dev_avg_pls_1SD = abs(I_p_wave_areas_1_right_septum_dev_avg_pls_1SD);
for i = 1:length(I_p_wave_areas_1_right_septum_sorted)
    if I_p_wave_areas_1_right_septum_dev_avg_pls_1SD(i) == min(I_p_wave_areas_1_right_septum_dev_avg_pls_1SD)
        I_p_wave_areas_1_right_septum_dev_avg_pls_1SD_index = i;
    end
end
I_p_wave_areas_1_right_septum_dev_avg_min_2SD = area_1(length(I_p_wave_areas_1_right_septum_sorted),1);
for i = 1:length(I_p_wave_areas_1_right_septum_sorted)
    I_p_wave_areas_1_right_septum_dev_avg_min_2SD(i) = I_p_wave_areas_1_right_septum_sorted(i) - (I_p_wave_areas_1_right_septum_average-2*I_p_wave_areas_1_right_septum_SD);
end
I_p_wave_areas_1_right_septum_dev_avg_min_2SD = abs(I_p_wave_areas_1_right_septum_dev_avg_min_2SD);
for i = 1:length(I_p_wave_areas_1_right_septum_sorted)
    if I_p_wave_areas_1_right_septum_dev_avg_min_2SD(i) == min(I_p_wave_areas_1_right_septum_dev_avg_min_2SD)
        I_p_wave_areas_1_dev_avg_min_2SD_index = i;
    end
end
I_p_wave_areas_1_right_septum_dev_avg_pls_2SD = area_1(length(I_p_wave_areas_1_right_septum_sorted),1);
for i = 1:length(I_p_wave_areas_1_right_septum_sorted)
    I_p_wave_areas_1_right_septum_dev_avg_pls_2SD(i) = I_p_wave_areas_1_right_septum_sorted(i) - (I_p_wave_areas_1_right_septum_average+2*I_p_wave_areas_1_right_septum_SD);
end
I_p_wave_areas_1_right_septum_dev_avg_pls_2SD = abs(I_p_wave_areas_1_right_septum_dev_avg_pls_2SD);
for i = 1:length(I_p_wave_areas_1_right_septum_sorted)
    if I_p_wave_areas_1_right_septum_dev_avg_pls_2SD(i) == min(I_p_wave_areas_1_right_septum_dev_avg_pls_2SD)
        I_p_wave_areas_1_right_septum_dev_avg_pls_2SD_index = i;
    end
end
I_p_wave_areas_1_right_septum_dev_avg_min_3SD = area_1(length(I_p_wave_areas_1_right_septum_sorted),1);
for i = 1:length(I_p_wave_areas_1_right_septum_sorted)
    I_p_wave_areas_1_right_septum_dev_avg_min_3SD(i) = I_p_wave_areas_1_right_septum_sorted(i) - (I_p_wave_areas_1_right_septum_average-3*I_p_wave_areas_1_right_septum_SD);
end
I_p_wave_areas_1_right_septum_dev_avg_min_3SD = abs(I_p_wave_areas_1_right_septum_dev_avg_min_3SD);
for i = 1:length(I_p_wave_areas_1_right_septum_sorted)
    if I_p_wave_areas_1_right_septum_dev_avg_min_3SD(i) == min(I_p_wave_areas_1_right_septum_dev_avg_min_3SD)
        I_p_wave_areas_1_right_septum_dev_avg_min_3SD_index = i;
    end
end
I_p_wave_areas_1_right_septum_dev_avg_pls_3SD = area_1(length(I_p_wave_areas_1_right_septum_sorted),1);
for i = 1:length(I_p_wave_areas_1_right_septum_sorted)
    I_p_wave_areas_1_right_septum_dev_avg_pls_3SD(i) = I_p_wave_areas_1_right_septum_sorted(i) - (I_p_wave_areas_1_right_septum_average+3*I_p_wave_areas_1_right_septum_SD);
end
I_p_wave_areas_1_right_septum_dev_avg_pls_3SD = abs(I_p_wave_areas_1_right_septum_dev_avg_pls_3SD);
for i = 1:length(I_p_wave_areas_1_right_septum_sorted)
    if I_p_wave_areas_1_right_septum_dev_avg_pls_3SD(i) == min(I_p_wave_areas_1_right_septum_dev_avg_pls_3SD)
        I_p_wave_areas_1_right_septum_dev_avg_pls_3SD_index = i;
    end
end
I_p_wave_areas_1_right_septum_tot_number_of_indices = length(I_p_wave_areas_1_right_septum_sorted);
%This is the approximation of the tot "area under curve." 
I_p_wave_areas_1_right_septum_quotient_within_1SD = (I_p_wave_areas_1_right_septum_dev_avg_pls_1SD_index - I_p_wave_areas_1_right_septumdev_avg_min_1SD_index)/I_p_wave_areas_1_right_septum_tot_number_of_indices;
I_p_wave_areas_1_right_septum_quotient_within_2SD = (I_p_wave_areas_1_right_septum_dev_avg_pls_2SD_index - I_p_wave_areas_1_dev_avg_min_2SD_index)/I_p_wave_areas_1_right_septum_tot_number_of_indices;
I_p_wave_areas_1_right_septum_quotient_within_3SD = (I_p_wave_areas_1_right_septum_dev_avg_pls_3SD_index - I_p_wave_areas_1_right_septum_dev_avg_min_3SD_index)/I_p_wave_areas_1_right_septum_tot_number_of_indices;
if (I_p_wave_areas_1_right_septum_quotient_within_1SD > 0.66 && I_p_wave_areas_1_right_septum_quotient_within_1SD < 0.70) && (I_p_wave_areas_1_right_septum_quotient_within_2SD > 0.93 && I_p_wave_areas_1_right_septum_quotient_within_2SD < 0.97) && (I_p_wave_areas_1_right_septum_quotient_within_3SD > 0.98 && I_p_wave_areas_1_right_septum_quotient_within_3SD < 1)
    isNormal(7,7) = 1;
end
if I_p_wave_areas_1_right_septum_quotient_within_1SD == 0
    isSpike(7,7) = 1;
end
I_p_wave_areas_1_right_atrial_appendage_sorted = sort(I_p_wave_areas_1_right_atrial_appendage);
I_p_wave_areas_1_right_atrial_appendage_dev_avg_min_1SD = area_1(length(I_p_wave_areas_1_right_atrial_appendage_sorted),1);
for i = 1:length(I_p_wave_areas_1_right_atrial_appendage_sorted)
    I_p_wave_areas_1_right_atrial_appendage_dev_avg_min_1SD(i) = I_p_wave_areas_1_right_atrial_appendage_sorted(i) - (I_p_wave_areas_1_right_atrial_appendage_average-I_p_wave_areas_1_right_atrial_appendage_SD);
end
I_p_wave_areas_1_right_atrial_appendage_dev_avg_min_1SD = abs(I_p_wave_areas_1_right_atrial_appendage_dev_avg_min_1SD);
for i = 1:length(I_p_wave_areas_1_right_atrial_appendage_sorted)
    if I_p_wave_areas_1_right_atrial_appendage_dev_avg_min_1SD(i) == min(I_p_wave_areas_1_right_atrial_appendage_dev_avg_min_1SD)
        I_p_wave_areas_1_right_atrial_appendage_dev_avg_min_1SD_index = i;
    end
end
I_p_wave_areas_1_right_atrial_appendage_dev_avg_pls_1SD = area_1(length(I_p_wave_areas_1_right_atrial_appendage_sorted),1);
for i = 1:length(I_p_wave_areas_1_right_atrial_appendage_sorted)
    I_p_wave_areas_1_right_atrial_appendage_dev_avg_pls_1SD(i) = I_p_wave_areas_1_right_atrial_appendage_sorted(i) - (I_p_wave_areas_1_right_atrial_appendage_average+I_p_wave_areas_1_right_atrial_appendage_SD);
end
I_p_wave_areas_1_right_atrial_appendage_dev_avg_pls_1SD = abs(I_p_wave_areas_1_right_atrial_appendage_dev_avg_pls_1SD);
for i = 1:length(I_p_wave_areas_1_right_atrial_appendage_sorted)
    if I_p_wave_areas_1_right_atrial_appendage_dev_avg_pls_1SD(i) == min(I_p_wave_areas_1_right_atrial_appendage_dev_avg_pls_1SD)
        I_p_wave_areas_1_right_atrial_appendage_dev_avg_pls_1SD_index = i;
    end
end
I_p_wave_areas_1_right_atrial_appendage_dev_avg_min_2SD = area_1(length(I_p_wave_areas_1_right_atrial_appendage_sorted),1);
for i = 1:length(I_p_wave_areas_1_right_atrial_appendage_sorted)
    I_p_wave_areas_1_right_atrial_appendage_dev_avg_min_2SD(i) = I_p_wave_areas_1_right_atrial_appendage_sorted(i) - (I_p_wave_areas_1_right_atrial_appendage_average-2*I_p_wave_areas_1_right_atrial_appendage_SD);
end
I_p_wave_areas_1_right_atrial_appendage_dev_avg_min_2SD = abs(I_p_wave_areas_1_right_atrial_appendage_dev_avg_min_2SD);
for i = 1:length(I_p_wave_areas_1_right_atrial_appendage_sorted)
    if I_p_wave_areas_1_right_atrial_appendage_dev_avg_min_2SD(i) == min(I_p_wave_areas_1_right_atrial_appendage_dev_avg_min_2SD)
        I_p_wave_areas_1_dev_avg_min_2SD_index = i;
    end
end
I_p_wave_areas_1_right_atrial_appendage_dev_avg_pls_2SD = area_1(length(I_p_wave_areas_1_right_atrial_appendage_sorted),1);
for i = 1:length(I_p_wave_areas_1_right_atrial_appendage_sorted)
    I_p_wave_areas_1_right_atrial_appendage_dev_avg_pls_2SD(i) = I_p_wave_areas_1_right_atrial_appendage_sorted(i) - (I_p_wave_areas_1_right_atrial_appendage_average+2*I_p_wave_areas_1_right_atrial_appendage_SD);
end
I_p_wave_areas_1_right_atrial_appendage_dev_avg_pls_2SD = abs(I_p_wave_areas_1_right_atrial_appendage_dev_avg_pls_2SD);
for i = 1:length(I_p_wave_areas_1_right_atrial_appendage_sorted)
    if I_p_wave_areas_1_right_atrial_appendage_dev_avg_pls_2SD(i) == min(I_p_wave_areas_1_right_atrial_appendage_dev_avg_pls_2SD)
        I_p_wave_areas_1_right_atrial_appendage_dev_avg_pls_2SD_index = i;
    end
end
I_p_wave_areas_1_right_atrial_appendage_dev_avg_min_3SD = area_1(length(I_p_wave_areas_1_right_atrial_appendage_sorted),1);
for i = 1:length(I_p_wave_areas_1_right_atrial_appendage_sorted)
    I_p_wave_areas_1_right_atrial_appendage_dev_avg_min_3SD(i) = I_p_wave_areas_1_right_atrial_appendage_sorted(i) - (I_p_wave_areas_1_right_atrial_appendage_average-3*I_p_wave_areas_1_right_atrial_appendage_SD);
end
I_p_wave_areas_1_right_atrial_appendage_dev_avg_min_3SD = abs(I_p_wave_areas_1_right_atrial_appendage_dev_avg_min_3SD);
for i = 1:length(I_p_wave_areas_1_right_atrial_appendage_sorted)
    if I_p_wave_areas_1_right_atrial_appendage_dev_avg_min_3SD(i) == min(I_p_wave_areas_1_right_atrial_appendage_dev_avg_min_3SD)
        I_p_wave_areas_1_right_atrial_appendage_dev_avg_min_3SD_index = i;
    end
end
I_p_wave_areas_1_right_atrial_appendage_dev_avg_pls_3SD = area_1(length(I_p_wave_areas_1_right_atrial_appendage_sorted),1);
for i = 1:length(I_p_wave_areas_1_right_atrial_appendage_sorted)
    I_p_wave_areas_1_right_atrial_appendage_dev_avg_pls_3SD(i) = I_p_wave_areas_1_right_atrial_appendage_sorted(i) - (I_p_wave_areas_1_right_atrial_appendage_average+3*I_p_wave_areas_1_right_atrial_appendage_SD);
end
I_p_wave_areas_1_right_atrial_appendage_dev_avg_pls_3SD = abs(I_p_wave_areas_1_right_atrial_appendage_dev_avg_pls_3SD);
for i = 1:length(I_p_wave_areas_1_right_atrial_appendage_sorted)
    if I_p_wave_areas_1_right_atrial_appendage_dev_avg_pls_3SD(i) == min(I_p_wave_areas_1_right_atrial_appendage_dev_avg_pls_3SD)
        I_p_wave_areas_1_right_atrial_appendage_dev_avg_pls_3SD_index = i;
    end
end
I_p_wave_areas_1_right_atrial_appendage_tot_number_of_indices = length(I_p_wave_areas_1_right_atrial_appendage_sorted);
%This is the approximation of the tot "area under curve." 
I_p_wave_areas_1_right_atrial_appendage_quotient_within_1SD = (I_p_wave_areas_1_right_atrial_appendage_dev_avg_pls_1SD_index - I_p_wave_areas_1_right_atrial_appendagedev_avg_min_1SD_index)/I_p_wave_areas_1_right_atrial_appendage_tot_number_of_indices;
I_p_wave_areas_1_right_atrial_appendage_quotient_within_2SD = (I_p_wave_areas_1_right_atrial_appendage_dev_avg_pls_2SD_index - I_p_wave_areas_1_dev_avg_min_2SD_index)/I_p_wave_areas_1_right_atrial_appendage_tot_number_of_indices;
I_p_wave_areas_1_right_atrial_appendage_quotient_within_3SD = (I_p_wave_areas_1_right_atrial_appendage_dev_avg_pls_3SD_index - I_p_wave_areas_1_right_atrial_appendage_dev_avg_min_3SD_index)/I_p_wave_areas_1_right_atrial_appendage_tot_number_of_indices;
if (I_p_wave_areas_1_right_atrial_appendage_quotient_within_1SD > 0.66 && I_p_wave_areas_1_right_atrial_appendage_quotient_within_1SD < 0.70) && (I_p_wave_areas_1_right_atrial_appendage_quotient_within_2SD > 0.93 && I_p_wave_areas_1_right_atrial_appendage_quotient_within_2SD < 0.97) && (I_p_wave_areas_1_right_atrial_appendage_quotient_within_3SD > 0.98 && I_p_wave_areas_1_right_atrial_appendage_quotient_within_3SD < 1)
    isNormal(7,8) = 1;
end
if I_p_wave_areas_1_right_atrial_appendage_quotient_within_1SD == 0
    isSpike(7,8) = 1;
end
I_p_wave_areas_1_pulmonary_veins_sorted = sort(I_p_wave_areas_1_pulmonary_veins);
I_p_wave_areas_1_pulmonary_veins_dev_avg_min_1SD = area_1(length(I_p_wave_areas_1_pulmonary_veins_sorted),1);
for i = 1:length(I_p_wave_areas_1_pulmonary_veins_sorted)
    I_p_wave_areas_1_pulmonary_veins_dev_avg_min_1SD(i) = I_p_wave_areas_1_pulmonary_veins_sorted(i) - (I_p_wave_areas_1_pulmonary_veins_average-I_p_wave_areas_1_pulmonary_veins_SD);
end
I_p_wave_areas_1_pulmonary_veins_dev_avg_min_1SD = abs(I_p_wave_areas_1_pulmonary_veins_dev_avg_min_1SD);
for i = 1:length(I_p_wave_areas_1_pulmonary_veins_sorted)
    if I_p_wave_areas_1_pulmonary_veins_dev_avg_min_1SD(i) == min(I_p_wave_areas_1_pulmonary_veins_dev_avg_min_1SD)
        I_p_wave_areas_1_pulmonary_veins_dev_avg_min_1SD_index = i;
    end
end
I_p_wave_areas_1_pulmonary_veins_dev_avg_pls_1SD = area_1(length(I_p_wave_areas_1_pulmonary_veins_sorted),1);
for i = 1:length(I_p_wave_areas_1_pulmonary_veins_sorted)
    I_p_wave_areas_1_pulmonary_veins_dev_avg_pls_1SD(i) = I_p_wave_areas_1_pulmonary_veins_sorted(i) - (I_p_wave_areas_1_pulmonary_veins_average+I_p_wave_areas_1_pulmonary_veins_SD);
end
I_p_wave_areas_1_pulmonary_veins_dev_avg_pls_1SD = abs(I_p_wave_areas_1_pulmonary_veins_dev_avg_pls_1SD);
for i = 1:length(I_p_wave_areas_1_pulmonary_veins_sorted)
    if I_p_wave_areas_1_pulmonary_veins_dev_avg_pls_1SD(i) == min(I_p_wave_areas_1_pulmonary_veins_dev_avg_pls_1SD)
        I_p_wave_areas_1_pulmonary_veins_dev_avg_pls_1SD_index = i;
    end
end
I_p_wave_areas_1_pulmonary_veins_dev_avg_min_2SD = area_1(length(I_p_wave_areas_1_pulmonary_veins_sorted),1);
for i = 1:length(I_p_wave_areas_1_pulmonary_veins_sorted)
    I_p_wave_areas_1_pulmonary_veins_dev_avg_min_2SD(i) = I_p_wave_areas_1_pulmonary_veins_sorted(i) - (I_p_wave_areas_1_pulmonary_veins_average-2*I_p_wave_areas_1_pulmonary_veins_SD);
end
I_p_wave_areas_1_pulmonary_veins_dev_avg_min_2SD = abs(I_p_wave_areas_1_pulmonary_veins_dev_avg_min_2SD);
for i = 1:length(I_p_wave_areas_1_pulmonary_veins_sorted)
    if I_p_wave_areas_1_pulmonary_veins_dev_avg_min_2SD(i) == min(I_p_wave_areas_1_pulmonary_veins_dev_avg_min_2SD)
        I_p_wave_areas_1_dev_avg_min_2SD_index = i;
    end
end
I_p_wave_areas_1_pulmonary_veins_dev_avg_pls_2SD = area_1(length(I_p_wave_areas_1_pulmonary_veins_sorted),1);
for i = 1:length(I_p_wave_areas_1_pulmonary_veins_sorted)
    I_p_wave_areas_1_pulmonary_veins_dev_avg_pls_2SD(i) = I_p_wave_areas_1_pulmonary_veins_sorted(i) - (I_p_wave_areas_1_pulmonary_veins_average+2*I_p_wave_areas_1_pulmonary_veins_SD);
end
I_p_wave_areas_1_pulmonary_veins_dev_avg_pls_2SD = abs(I_p_wave_areas_1_pulmonary_veins_dev_avg_pls_2SD);
for i = 1:length(I_p_wave_areas_1_pulmonary_veins_sorted)
    if I_p_wave_areas_1_pulmonary_veins_dev_avg_pls_2SD(i) == min(I_p_wave_areas_1_pulmonary_veins_dev_avg_pls_2SD)
        I_p_wave_areas_1_pulmonary_veins_dev_avg_pls_2SD_index = i;
    end
end
I_p_wave_areas_1_pulmonary_veins_dev_avg_min_3SD = area_1(length(I_p_wave_areas_1_pulmonary_veins_sorted),1);
for i = 1:length(I_p_wave_areas_1_pulmonary_veins_sorted)
    I_p_wave_areas_1_pulmonary_veins_dev_avg_min_3SD(i) = I_p_wave_areas_1_pulmonary_veins_sorted(i) - (I_p_wave_areas_1_pulmonary_veins_average-3*I_p_wave_areas_1_pulmonary_veins_SD);
end
I_p_wave_areas_1_pulmonary_veins_dev_avg_min_3SD = abs(I_p_wave_areas_1_pulmonary_veins_dev_avg_min_3SD);
for i = 1:length(I_p_wave_areas_1_pulmonary_veins_sorted)
    if I_p_wave_areas_1_pulmonary_veins_dev_avg_min_3SD(i) == min(I_p_wave_areas_1_pulmonary_veins_dev_avg_min_3SD)
        I_p_wave_areas_1_pulmonary_veins_dev_avg_min_3SD_index = i;
    end
end
I_p_wave_areas_1_pulmonary_veins_dev_avg_pls_3SD = area_1(length(I_p_wave_areas_1_pulmonary_veins_sorted),1);
for i = 1:length(I_p_wave_areas_1_pulmonary_veins_sorted)
    I_p_wave_areas_1_pulmonary_veins_dev_avg_pls_3SD(i) = I_p_wave_areas_1_pulmonary_veins_sorted(i) - (I_p_wave_areas_1_pulmonary_veins_average+3*I_p_wave_areas_1_pulmonary_veins_SD);
end
I_p_wave_areas_1_pulmonary_veins_dev_avg_pls_3SD = abs(I_p_wave_areas_1_pulmonary_veins_dev_avg_pls_3SD);
for i = 1:length(I_p_wave_areas_1_pulmonary_veins_sorted)
    if I_p_wave_areas_1_pulmonary_veins_dev_avg_pls_3SD(i) == min(I_p_wave_areas_1_pulmonary_veins_dev_avg_pls_3SD)
        I_p_wave_areas_1_pulmonary_veins_dev_avg_pls_3SD_index = i;
    end
end
I_p_wave_areas_1_pulmonary_veins_tot_number_of_indices = length(I_p_wave_areas_1_pulmonary_veins_sorted);
%This is the approximation of the tot "area under curve." 
I_p_wave_areas_1_pulmonary_veins_quotient_within_1SD = (I_p_wave_areas_1_pulmonary_veins_dev_avg_pls_1SD_index - I_p_wave_areas_1_pulmonary_veinsdev_avg_min_1SD_index)/I_p_wave_areas_1_pulmonary_veins_tot_number_of_indices;
I_p_wave_areas_1_pulmonary_veins_quotient_within_2SD = (I_p_wave_areas_1_pulmonary_veins_dev_avg_pls_2SD_index - I_p_wave_areas_1_dev_avg_min_2SD_index)/I_p_wave_areas_1_pulmonary_veins_tot_number_of_indices;
I_p_wave_areas_1_pulmonary_veins_quotient_within_3SD = (I_p_wave_areas_1_pulmonary_veins_dev_avg_pls_3SD_index - I_p_wave_areas_1_pulmonary_veins_dev_avg_min_3SD_index)/I_p_wave_areas_1_pulmonary_veins_tot_number_of_indices;
if (I_p_wave_areas_1_pulmonary_veins_quotient_within_1SD > 0.66 && I_p_wave_areas_1_pulmonary_veins_quotient_within_1SD < 0.70) && (I_p_wave_areas_1_pulmonary_veins_quotient_within_2SD > 0.93 && I_p_wave_areas_1_pulmonary_veins_quotient_within_2SD < 0.97) && (I_p_wave_areas_1_pulmonary_veins_quotient_within_3SD > 0.98 && I_p_wave_areas_1_pulmonary_veins_quotient_within_3SD < 1)
    isNormal(7,9) = 1;
end
if I_p_wave_areas_1_pulmonary_veins_quotient_within_1SD == 0
    isSpike(7,9) = 1;
end
I_p_wave_areas_1_mitral_annulus_sorted = sort(I_p_wave_areas_1_mitral_annulus);
I_p_wave_areas_1_mitral_annulus_dev_avg_min_1SD = area_1(length(I_p_wave_areas_1_mitral_annulus_sorted),1);
for i = 1:length(I_p_wave_areas_1_mitral_annulus_sorted)
    I_p_wave_areas_1_mitral_annulus_dev_avg_min_1SD(i) = I_p_wave_areas_1_mitral_annulus_sorted(i) - (I_p_wave_areas_1_mitral_annulus_average-I_p_wave_areas_1_mitral_annulus_SD);
end
I_p_wave_areas_1_mitral_annulus_dev_avg_min_1SD = abs(I_p_wave_areas_1_mitral_annulus_dev_avg_min_1SD);
for i = 1:length(I_p_wave_areas_1_mitral_annulus_sorted)
    if I_p_wave_areas_1_mitral_annulus_dev_avg_min_1SD(i) == min(I_p_wave_areas_1_mitral_annulus_dev_avg_min_1SD)
        I_p_wave_areas_1_mitral_annulus_dev_avg_min_1SD_index = i;
    end
end
I_p_wave_areas_1_mitral_annulus_dev_avg_pls_1SD = area_1(length(I_p_wave_areas_1_mitral_annulus_sorted),1);
for i = 1:length(I_p_wave_areas_1_mitral_annulus_sorted)
    I_p_wave_areas_1_mitral_annulus_dev_avg_pls_1SD(i) = I_p_wave_areas_1_mitral_annulus_sorted(i) - (I_p_wave_areas_1_mitral_annulus_average+I_p_wave_areas_1_mitral_annulus_SD);
end
I_p_wave_areas_1_mitral_annulus_dev_avg_pls_1SD = abs(I_p_wave_areas_1_mitral_annulus_dev_avg_pls_1SD);
for i = 1:length(I_p_wave_areas_1_mitral_annulus_sorted)
    if I_p_wave_areas_1_mitral_annulus_dev_avg_pls_1SD(i) == min(I_p_wave_areas_1_mitral_annulus_dev_avg_pls_1SD)
        I_p_wave_areas_1_mitral_annulus_dev_avg_pls_1SD_index = i;
    end
end
I_p_wave_areas_1_mitral_annulus_dev_avg_min_2SD = area_1(length(I_p_wave_areas_1_mitral_annulus_sorted),1);
for i = 1:length(I_p_wave_areas_1_mitral_annulus_sorted)
    I_p_wave_areas_1_mitral_annulus_dev_avg_min_2SD(i) = I_p_wave_areas_1_mitral_annulus_sorted(i) - (I_p_wave_areas_1_mitral_annulus_average-2*I_p_wave_areas_1_mitral_annulus_SD);
end
I_p_wave_areas_1_mitral_annulus_dev_avg_min_2SD = abs(I_p_wave_areas_1_mitral_annulus_dev_avg_min_2SD);
for i = 1:length(I_p_wave_areas_1_mitral_annulus_sorted)
    if I_p_wave_areas_1_mitral_annulus_dev_avg_min_2SD(i) == min(I_p_wave_areas_1_mitral_annulus_dev_avg_min_2SD)
        I_p_wave_areas_1_dev_avg_min_2SD_index = i;
    end
end
I_p_wave_areas_1_mitral_annulus_dev_avg_pls_2SD = area_1(length(I_p_wave_areas_1_mitral_annulus_sorted),1);
for i = 1:length(I_p_wave_areas_1_mitral_annulus_sorted)
    I_p_wave_areas_1_mitral_annulus_dev_avg_pls_2SD(i) = I_p_wave_areas_1_mitral_annulus_sorted(i) - (I_p_wave_areas_1_mitral_annulus_average+2*I_p_wave_areas_1_mitral_annulus_SD);
end
I_p_wave_areas_1_mitral_annulus_dev_avg_pls_2SD = abs(I_p_wave_areas_1_mitral_annulus_dev_avg_pls_2SD);
for i = 1:length(I_p_wave_areas_1_mitral_annulus_sorted)
    if I_p_wave_areas_1_mitral_annulus_dev_avg_pls_2SD(i) == min(I_p_wave_areas_1_mitral_annulus_dev_avg_pls_2SD)
        I_p_wave_areas_1_mitral_annulus_dev_avg_pls_2SD_index = i;
    end
end
I_p_wave_areas_1_mitral_annulus_dev_avg_min_3SD = area_1(length(I_p_wave_areas_1_mitral_annulus_sorted),1);
for i = 1:length(I_p_wave_areas_1_mitral_annulus_sorted)
    I_p_wave_areas_1_mitral_annulus_dev_avg_min_3SD(i) = I_p_wave_areas_1_mitral_annulus_sorted(i) - (I_p_wave_areas_1_mitral_annulus_average-3*I_p_wave_areas_1_mitral_annulus_SD);
end
I_p_wave_areas_1_mitral_annulus_dev_avg_min_3SD = abs(I_p_wave_areas_1_mitral_annulus_dev_avg_min_3SD);
for i = 1:length(I_p_wave_areas_1_mitral_annulus_sorted)
    if I_p_wave_areas_1_mitral_annulus_dev_avg_min_3SD(i) == min(I_p_wave_areas_1_mitral_annulus_dev_avg_min_3SD)
        I_p_wave_areas_1_mitral_annulus_dev_avg_min_3SD_index = i;
    end
end
I_p_wave_areas_1_mitral_annulus_dev_avg_pls_3SD = area_1(length(I_p_wave_areas_1_mitral_annulus_sorted),1);
for i = 1:length(I_p_wave_areas_1_mitral_annulus_sorted)
    I_p_wave_areas_1_mitral_annulus_dev_avg_pls_3SD(i) = I_p_wave_areas_1_mitral_annulus_sorted(i) - (I_p_wave_areas_1_mitral_annulus_average+3*I_p_wave_areas_1_mitral_annulus_SD);
end
I_p_wave_areas_1_mitral_annulus_dev_avg_pls_3SD = abs(I_p_wave_areas_1_mitral_annulus_dev_avg_pls_3SD);
for i = 1:length(I_p_wave_areas_1_mitral_annulus_sorted)
    if I_p_wave_areas_1_mitral_annulus_dev_avg_pls_3SD(i) == min(I_p_wave_areas_1_mitral_annulus_dev_avg_pls_3SD)
        I_p_wave_areas_1_mitral_annulus_dev_avg_pls_3SD_index = i;
    end
end
I_p_wave_areas_1_mitral_annulus_tot_number_of_indices = length(I_p_wave_areas_1_mitral_annulus_sorted);
%This is the approximation of the tot "area under curve." 
I_p_wave_areas_1_mitral_annulus_quotient_within_1SD = (I_p_wave_areas_1_mitral_annulus_dev_avg_pls_1SD_index - I_p_wave_areas_1_mitral_annulusdev_avg_min_1SD_index)/I_p_wave_areas_1_mitral_annulus_tot_number_of_indices;
I_p_wave_areas_1_mitral_annulus_quotient_within_2SD = (I_p_wave_areas_1_mitral_annulus_dev_avg_pls_2SD_index - I_p_wave_areas_1_dev_avg_min_2SD_index)/I_p_wave_areas_1_mitral_annulus_tot_number_of_indices;
I_p_wave_areas_1_mitral_annulus_quotient_within_3SD = (I_p_wave_areas_1_mitral_annulus_dev_avg_pls_3SD_index - I_p_wave_areas_1_mitral_annulus_dev_avg_min_3SD_index)/I_p_wave_areas_1_mitral_annulus_tot_number_of_indices;
if (I_p_wave_areas_1_mitral_annulus_quotient_within_1SD > 0.66 && I_p_wave_areas_1_mitral_annulus_quotient_within_1SD < 0.70) && (I_p_wave_areas_1_mitral_annulus_quotient_within_2SD > 0.93 && I_p_wave_areas_1_mitral_annulus_quotient_within_2SD < 0.97) && (I_p_wave_areas_1_mitral_annulus_quotient_within_3SD > 0.98 && I_p_wave_areas_1_mitral_annulus_quotient_within_3SD < 1)
    isNormal(7,10) = 1;
end
if I_p_wave_areas_1_mitral_annulus_quotient_within_1SD == 0
    isSpike(7,10) = 1;
end
I_p_wave_areas_1_CS_body_sorted = sort(I_p_wave_areas_1_CS_body);
I_p_wave_areas_1_CS_body_dev_avg_min_1SD = area_1(length(I_p_wave_areas_1_CS_body_sorted),1);
for i = 1:length(I_p_wave_areas_1_CS_body_sorted)
    I_p_wave_areas_1_CS_body_dev_avg_min_1SD(i) = I_p_wave_areas_1_CS_body_sorted(i) - (I_p_wave_areas_1_CS_body_average-I_p_wave_areas_1_CS_body_SD);
end
I_p_wave_areas_1_CS_body_dev_avg_min_1SD = abs(I_p_wave_areas_1_CS_body_dev_avg_min_1SD);
for i = 1:length(I_p_wave_areas_1_CS_body_sorted)
    if I_p_wave_areas_1_CS_body_dev_avg_min_1SD(i) == min(I_p_wave_areas_1_CS_body_dev_avg_min_1SD)
        I_p_wave_areas_1_CS_body_dev_avg_min_1SD_index = i;
    end
end
I_p_wave_areas_1_CS_body_dev_avg_pls_1SD = area_1(length(I_p_wave_areas_1_CS_body_sorted),1);
for i = 1:length(I_p_wave_areas_1_CS_body_sorted)
    I_p_wave_areas_1_CS_body_dev_avg_pls_1SD(i) = I_p_wave_areas_1_CS_body_sorted(i) - (I_p_wave_areas_1_CS_body_average+I_p_wave_areas_1_CS_body_SD);
end
I_p_wave_areas_1_CS_body_dev_avg_pls_1SD = abs(I_p_wave_areas_1_CS_body_dev_avg_pls_1SD);
for i = 1:length(I_p_wave_areas_1_CS_body_sorted)
    if I_p_wave_areas_1_CS_body_dev_avg_pls_1SD(i) == min(I_p_wave_areas_1_CS_body_dev_avg_pls_1SD)
        I_p_wave_areas_1_CS_body_dev_avg_pls_1SD_index = i;
    end
end
I_p_wave_areas_1_CS_body_dev_avg_min_2SD = area_1(length(I_p_wave_areas_1_CS_body_sorted),1);
for i = 1:length(I_p_wave_areas_1_CS_body_sorted)
    I_p_wave_areas_1_CS_body_dev_avg_min_2SD(i) = I_p_wave_areas_1_CS_body_sorted(i) - (I_p_wave_areas_1_CS_body_average-2*I_p_wave_areas_1_CS_body_SD);
end
I_p_wave_areas_1_CS_body_dev_avg_min_2SD = abs(I_p_wave_areas_1_CS_body_dev_avg_min_2SD);
for i = 1:length(I_p_wave_areas_1_CS_body_sorted)
    if I_p_wave_areas_1_CS_body_dev_avg_min_2SD(i) == min(I_p_wave_areas_1_CS_body_dev_avg_min_2SD)
        I_p_wave_areas_1_dev_avg_min_2SD_index = i;
    end
end
I_p_wave_areas_1_CS_body_dev_avg_pls_2SD = area_1(length(I_p_wave_areas_1_CS_body_sorted),1);
for i = 1:length(I_p_wave_areas_1_CS_body_sorted)
    I_p_wave_areas_1_CS_body_dev_avg_pls_2SD(i) = I_p_wave_areas_1_CS_body_sorted(i) - (I_p_wave_areas_1_CS_body_average+2*I_p_wave_areas_1_CS_body_SD);
end
I_p_wave_areas_1_CS_body_dev_avg_pls_2SD = abs(I_p_wave_areas_1_CS_body_dev_avg_pls_2SD);
for i = 1:length(I_p_wave_areas_1_CS_body_sorted)
    if I_p_wave_areas_1_CS_body_dev_avg_pls_2SD(i) == min(I_p_wave_areas_1_CS_body_dev_avg_pls_2SD)
        I_p_wave_areas_1_CS_body_dev_avg_pls_2SD_index = i;
    end
end
I_p_wave_areas_1_CS_body_dev_avg_min_3SD = area_1(length(I_p_wave_areas_1_CS_body_sorted),1);
for i = 1:length(I_p_wave_areas_1_CS_body_sorted)
    I_p_wave_areas_1_CS_body_dev_avg_min_3SD(i) = I_p_wave_areas_1_CS_body_sorted(i) - (I_p_wave_areas_1_CS_body_average-3*I_p_wave_areas_1_CS_body_SD);
end
I_p_wave_areas_1_CS_body_dev_avg_min_3SD = abs(I_p_wave_areas_1_CS_body_dev_avg_min_3SD);
for i = 1:length(I_p_wave_areas_1_CS_body_sorted)
    if I_p_wave_areas_1_CS_body_dev_avg_min_3SD(i) == min(I_p_wave_areas_1_CS_body_dev_avg_min_3SD)
        I_p_wave_areas_1_CS_body_dev_avg_min_3SD_index = i;
    end
end
I_p_wave_areas_1_CS_body_dev_avg_pls_3SD = area_1(length(I_p_wave_areas_1_CS_body_sorted),1);
for i = 1:length(I_p_wave_areas_1_CS_body_sorted)
    I_p_wave_areas_1_CS_body_dev_avg_pls_3SD(i) = I_p_wave_areas_1_CS_body_sorted(i) - (I_p_wave_areas_1_CS_body_average+3*I_p_wave_areas_1_CS_body_SD);
end
I_p_wave_areas_1_CS_body_dev_avg_pls_3SD = abs(I_p_wave_areas_1_CS_body_dev_avg_pls_3SD);
for i = 1:length(I_p_wave_areas_1_CS_body_sorted)
    if I_p_wave_areas_1_CS_body_dev_avg_pls_3SD(i) == min(I_p_wave_areas_1_CS_body_dev_avg_pls_3SD)
        I_p_wave_areas_1_CS_body_dev_avg_pls_3SD_index = i;
    end
end
I_p_wave_areas_1_CS_body_tot_number_of_indices = length(I_p_wave_areas_1_CS_body_sorted);
%This is the approximation of the tot "area under curve." 
I_p_wave_areas_1_CS_body_quotient_within_1SD = (I_p_wave_areas_1_CS_body_dev_avg_pls_1SD_index - I_p_wave_areas_1_CS_bodydev_avg_min_1SD_index)/I_p_wave_areas_1_CS_body_tot_number_of_indices;
I_p_wave_areas_1_CS_body_quotient_within_2SD = (I_p_wave_areas_1_CS_body_dev_avg_pls_2SD_index - I_p_wave_areas_1_dev_avg_min_2SD_index)/I_p_wave_areas_1_CS_body_tot_number_of_indices;
I_p_wave_areas_1_CS_body_quotient_within_3SD = (I_p_wave_areas_1_CS_body_dev_avg_pls_3SD_index - I_p_wave_areas_1_CS_body_dev_avg_min_3SD_index)/I_p_wave_areas_1_CS_body_tot_number_of_indices;
if (I_p_wave_areas_1_CS_body_quotient_within_1SD > 0.66 && I_p_wave_areas_1_CS_body_quotient_within_1SD < 0.70) && (I_p_wave_areas_1_CS_body_quotient_within_2SD > 0.93 && I_p_wave_areas_1_CS_body_quotient_within_2SD < 0.97) && (I_p_wave_areas_1_CS_body_quotient_within_3SD > 0.98 && I_p_wave_areas_1_CS_body_quotient_within_3SD < 1)
    isNormal(7,11) = 1;
end
if I_p_wave_areas_1_CS_body_quotient_within_1SD == 0
    isSpike(7,11) = 1;
end
I_p_wave_areas_1_left_septum_sorted = sort(I_p_wave_areas_1_left_septum);
I_p_wave_areas_1_left_septum_dev_avg_min_1SD = area_1(length(I_p_wave_areas_1_left_septum_sorted),1);
for i = 1:length(I_p_wave_areas_1_left_septum_sorted)
    I_p_wave_areas_1_left_septum_dev_avg_min_1SD(i) = I_p_wave_areas_1_left_septum_sorted(i) - (I_p_wave_areas_1_left_septum_average-I_p_wave_areas_1_left_septum_SD);
end
I_p_wave_areas_1_left_septum_dev_avg_min_1SD = abs(I_p_wave_areas_1_left_septum_dev_avg_min_1SD);
for i = 1:length(I_p_wave_areas_1_left_septum_sorted)
    if I_p_wave_areas_1_left_septum_dev_avg_min_1SD(i) == min(I_p_wave_areas_1_left_septum_dev_avg_min_1SD)
        I_p_wave_areas_1_left_septum_dev_avg_min_1SD_index = i;
    end
end
I_p_wave_areas_1_left_septum_dev_avg_pls_1SD = area_1(length(I_p_wave_areas_1_left_septum_sorted),1);
for i = 1:length(I_p_wave_areas_1_left_septum_sorted)
    I_p_wave_areas_1_left_septum_dev_avg_pls_1SD(i) = I_p_wave_areas_1_left_septum_sorted(i) - (I_p_wave_areas_1_left_septum_average+I_p_wave_areas_1_left_septum_SD);
end
I_p_wave_areas_1_left_septum_dev_avg_pls_1SD = abs(I_p_wave_areas_1_left_septum_dev_avg_pls_1SD);
for i = 1:length(I_p_wave_areas_1_left_septum_sorted)
    if I_p_wave_areas_1_left_septum_dev_avg_pls_1SD(i) == min(I_p_wave_areas_1_left_septum_dev_avg_pls_1SD)
        I_p_wave_areas_1_left_septum_dev_avg_pls_1SD_index = i;
    end
end
I_p_wave_areas_1_left_septum_dev_avg_min_2SD = area_1(length(I_p_wave_areas_1_left_septum_sorted),1);
for i = 1:length(I_p_wave_areas_1_left_septum_sorted)
    I_p_wave_areas_1_left_septum_dev_avg_min_2SD(i) = I_p_wave_areas_1_left_septum_sorted(i) - (I_p_wave_areas_1_left_septum_average-2*I_p_wave_areas_1_left_septum_SD);
end
I_p_wave_areas_1_left_septum_dev_avg_min_2SD = abs(I_p_wave_areas_1_left_septum_dev_avg_min_2SD);
for i = 1:length(I_p_wave_areas_1_left_septum_sorted)
    if I_p_wave_areas_1_left_septum_dev_avg_min_2SD(i) == min(I_p_wave_areas_1_left_septum_dev_avg_min_2SD)
        I_p_wave_areas_1_dev_avg_min_2SD_index = i;
    end
end
I_p_wave_areas_1_left_septum_dev_avg_pls_2SD = area_1(length(I_p_wave_areas_1_left_septum_sorted),1);
for i = 1:length(I_p_wave_areas_1_left_septum_sorted)
    I_p_wave_areas_1_left_septum_dev_avg_pls_2SD(i) = I_p_wave_areas_1_left_septum_sorted(i) - (I_p_wave_areas_1_left_septum_average+2*I_p_wave_areas_1_left_septum_SD);
end
I_p_wave_areas_1_left_septum_dev_avg_pls_2SD = abs(I_p_wave_areas_1_left_septum_dev_avg_pls_2SD);
for i = 1:length(I_p_wave_areas_1_left_septum_sorted)
    if I_p_wave_areas_1_left_septum_dev_avg_pls_2SD(i) == min(I_p_wave_areas_1_left_septum_dev_avg_pls_2SD)
        I_p_wave_areas_1_left_septum_dev_avg_pls_2SD_index = i;
    end
end
I_p_wave_areas_1_left_septum_dev_avg_min_3SD = area_1(length(I_p_wave_areas_1_left_septum_sorted),1);
for i = 1:length(I_p_wave_areas_1_left_septum_sorted)
    I_p_wave_areas_1_left_septum_dev_avg_min_3SD(i) = I_p_wave_areas_1_left_septum_sorted(i) - (I_p_wave_areas_1_left_septum_average-3*I_p_wave_areas_1_left_septum_SD);
end
I_p_wave_areas_1_left_septum_dev_avg_min_3SD = abs(I_p_wave_areas_1_left_septum_dev_avg_min_3SD);
for i = 1:length(I_p_wave_areas_1_left_septum_sorted)
    if I_p_wave_areas_1_left_septum_dev_avg_min_3SD(i) == min(I_p_wave_areas_1_left_septum_dev_avg_min_3SD)
        I_p_wave_areas_1_left_septum_dev_avg_min_3SD_index = i;
    end
end
I_p_wave_areas_1_left_septum_dev_avg_pls_3SD = area_1(length(I_p_wave_areas_1_left_septum_sorted),1);
for i = 1:length(I_p_wave_areas_1_left_septum_sorted)
    I_p_wave_areas_1_left_septum_dev_avg_pls_3SD(i) = I_p_wave_areas_1_left_septum_sorted(i) - (I_p_wave_areas_1_left_septum_average+3*I_p_wave_areas_1_left_septum_SD);
end
I_p_wave_areas_1_left_septum_dev_avg_pls_3SD = abs(I_p_wave_areas_1_left_septum_dev_avg_pls_3SD);
for i = 1:length(I_p_wave_areas_1_left_septum_sorted)
    if I_p_wave_areas_1_left_septum_dev_avg_pls_3SD(i) == min(I_p_wave_areas_1_left_septum_dev_avg_pls_3SD)
        I_p_wave_areas_1_left_septum_dev_avg_pls_3SD_index = i;
    end
end
I_p_wave_areas_1_left_septum_tot_number_of_indices = length(I_p_wave_areas_1_left_septum_sorted);
%This is the approximation of the tot "area under curve." 
I_p_wave_areas_1_left_septum_quotient_within_1SD = (I_p_wave_areas_1_left_septum_dev_avg_pls_1SD_index - I_p_wave_areas_1_left_septumdev_avg_min_1SD_index)/I_p_wave_areas_1_left_septum_tot_number_of_indices;
I_p_wave_areas_1_left_septum_quotient_within_2SD = (I_p_wave_areas_1_left_septum_dev_avg_pls_2SD_index - I_p_wave_areas_1_dev_avg_min_2SD_index)/I_p_wave_areas_1_left_septum_tot_number_of_indices;
I_p_wave_areas_1_left_septum_quotient_within_3SD = (I_p_wave_areas_1_left_septum_dev_avg_pls_3SD_index - I_p_wave_areas_1_left_septum_dev_avg_min_3SD_index)/I_p_wave_areas_1_left_septum_tot_number_of_indices;
if (I_p_wave_areas_1_left_septum_quotient_within_1SD > 0.66 && I_p_wave_areas_1_left_septum_quotient_within_1SD < 0.70) && (I_p_wave_areas_1_left_septum_quotient_within_2SD > 0.93 && I_p_wave_areas_1_left_septum_quotient_within_2SD < 0.97) && (I_p_wave_areas_1_left_septum_quotient_within_3SD > 0.98 && I_p_wave_areas_1_left_septum_quotient_within_3SD < 1)
    isNormal(7,12) = 1;
end
if I_p_wave_areas_1_left_septum_quotient_within_1SD == 0
    isSpike(7,12) = 1;
end
I_p_wave_areas_1_left_atrial_appendage_sorted = sort(I_p_wave_areas_1_left_atrial_appendage);
I_p_wave_areas_1_left_atrial_appendage_dev_avg_min_1SD = area_1(length(I_p_wave_areas_1_left_atrial_appendage_sorted),1);
for i = 1:length(I_p_wave_areas_1_left_atrial_appendage_sorted)
    I_p_wave_areas_1_left_atrial_appendage_dev_avg_min_1SD(i) = I_p_wave_areas_1_left_atrial_appendage_sorted(i) - (I_p_wave_areas_1_left_atrial_appendage_average-I_p_wave_areas_1_left_atrial_appendage_SD);
end
I_p_wave_areas_1_left_atrial_appendage_dev_avg_min_1SD = abs(I_p_wave_areas_1_left_atrial_appendage_dev_avg_min_1SD);
for i = 1:length(I_p_wave_areas_1_left_atrial_appendage_sorted)
    if I_p_wave_areas_1_left_atrial_appendage_dev_avg_min_1SD(i) == min(I_p_wave_areas_1_left_atrial_appendage_dev_avg_min_1SD)
        I_p_wave_areas_1_left_atrial_appendage_dev_avg_min_1SD_index = i;
    end
end
I_p_wave_areas_1_left_atrial_appendage_dev_avg_pls_1SD = area_1(length(I_p_wave_areas_1_left_atrial_appendage_sorted),1);
for i = 1:length(I_p_wave_areas_1_left_atrial_appendage_sorted)
    I_p_wave_areas_1_left_atrial_appendage_dev_avg_pls_1SD(i) = I_p_wave_areas_1_left_atrial_appendage_sorted(i) - (I_p_wave_areas_1_left_atrial_appendage_average+I_p_wave_areas_1_left_atrial_appendage_SD);
end
I_p_wave_areas_1_left_atrial_appendage_dev_avg_pls_1SD = abs(I_p_wave_areas_1_left_atrial_appendage_dev_avg_pls_1SD);
for i = 1:length(I_p_wave_areas_1_left_atrial_appendage_sorted)
    if I_p_wave_areas_1_left_atrial_appendage_dev_avg_pls_1SD(i) == min(I_p_wave_areas_1_left_atrial_appendage_dev_avg_pls_1SD)
        I_p_wave_areas_1_left_atrial_appendage_dev_avg_pls_1SD_index = i;
    end
end
I_p_wave_areas_1_left_atrial_appendage_dev_avg_min_2SD = area_1(length(I_p_wave_areas_1_left_atrial_appendage_sorted),1);
for i = 1:length(I_p_wave_areas_1_left_atrial_appendage_sorted)
    I_p_wave_areas_1_left_atrial_appendage_dev_avg_min_2SD(i) = I_p_wave_areas_1_left_atrial_appendage_sorted(i) - (I_p_wave_areas_1_left_atrial_appendage_average-2*I_p_wave_areas_1_left_atrial_appendage_SD);
end
I_p_wave_areas_1_left_atrial_appendage_dev_avg_min_2SD = abs(I_p_wave_areas_1_left_atrial_appendage_dev_avg_min_2SD);
for i = 1:length(I_p_wave_areas_1_left_atrial_appendage_sorted)
    if I_p_wave_areas_1_left_atrial_appendage_dev_avg_min_2SD(i) == min(I_p_wave_areas_1_left_atrial_appendage_dev_avg_min_2SD)
        I_p_wave_areas_1_dev_avg_min_2SD_index = i;
    end
end
I_p_wave_areas_1_left_atrial_appendage_dev_avg_pls_2SD = area_1(length(I_p_wave_areas_1_left_atrial_appendage_sorted),1);
for i = 1:length(I_p_wave_areas_1_left_atrial_appendage_sorted)
    I_p_wave_areas_1_left_atrial_appendage_dev_avg_pls_2SD(i) = I_p_wave_areas_1_left_atrial_appendage_sorted(i) - (I_p_wave_areas_1_left_atrial_appendage_average+2*I_p_wave_areas_1_left_atrial_appendage_SD);
end
I_p_wave_areas_1_left_atrial_appendage_dev_avg_pls_2SD = abs(I_p_wave_areas_1_left_atrial_appendage_dev_avg_pls_2SD);
for i = 1:length(I_p_wave_areas_1_left_atrial_appendage_sorted)
    if I_p_wave_areas_1_left_atrial_appendage_dev_avg_pls_2SD(i) == min(I_p_wave_areas_1_left_atrial_appendage_dev_avg_pls_2SD)
        I_p_wave_areas_1_left_atrial_appendage_dev_avg_pls_2SD_index = i;
    end
end
I_p_wave_areas_1_left_atrial_appendage_dev_avg_min_3SD = area_1(length(I_p_wave_areas_1_left_atrial_appendage_sorted),1);
for i = 1:length(I_p_wave_areas_1_left_atrial_appendage_sorted)
    I_p_wave_areas_1_left_atrial_appendage_dev_avg_min_3SD(i) = I_p_wave_areas_1_left_atrial_appendage_sorted(i) - (I_p_wave_areas_1_left_atrial_appendage_average-3*I_p_wave_areas_1_left_atrial_appendage_SD);
end
I_p_wave_areas_1_left_atrial_appendage_dev_avg_min_3SD = abs(I_p_wave_areas_1_left_atrial_appendage_dev_avg_min_3SD);
for i = 1:length(I_p_wave_areas_1_left_atrial_appendage_sorted)
    if I_p_wave_areas_1_left_atrial_appendage_dev_avg_min_3SD(i) == min(I_p_wave_areas_1_left_atrial_appendage_dev_avg_min_3SD)
        I_p_wave_areas_1_left_atrial_appendage_dev_avg_min_3SD_index = i;
    end
end
I_p_wave_areas_1_left_atrial_appendage_dev_avg_pls_3SD = area_1(length(I_p_wave_areas_1_left_atrial_appendage_sorted),1);
for i = 1:length(I_p_wave_areas_1_left_atrial_appendage_sorted)
    I_p_wave_areas_1_left_atrial_appendage_dev_avg_pls_3SD(i) = I_p_wave_areas_1_left_atrial_appendage_sorted(i) - (I_p_wave_areas_1_left_atrial_appendage_average+3*I_p_wave_areas_1_left_atrial_appendage_SD);
end
I_p_wave_areas_1_left_atrial_appendage_dev_avg_pls_3SD = abs(I_p_wave_areas_1_left_atrial_appendage_dev_avg_pls_3SD);
for i = 1:length(I_p_wave_areas_1_left_atrial_appendage_sorted)
    if I_p_wave_areas_1_left_atrial_appendage_dev_avg_pls_3SD(i) == min(I_p_wave_areas_1_left_atrial_appendage_dev_avg_pls_3SD)
        I_p_wave_areas_1_left_atrial_appendage_dev_avg_pls_3SD_index = i;
    end
end
I_p_wave_areas_1_left_atrial_appendage_tot_number_of_indices = length(I_p_wave_areas_1_left_atrial_appendage_sorted);
%This is the approximation of the tot "area under curve." 
I_p_wave_areas_1_left_atrial_appendage_quotient_within_1SD = (I_p_wave_areas_1_left_atrial_appendage_dev_avg_pls_1SD_index - I_p_wave_areas_1_left_atrial_appendagedev_avg_min_1SD_index)/I_p_wave_areas_1_left_atrial_appendage_tot_number_of_indices;
I_p_wave_areas_1_left_atrial_appendage_quotient_within_2SD = (I_p_wave_areas_1_left_atrial_appendage_dev_avg_pls_2SD_index - I_p_wave_areas_1_dev_avg_min_2SD_index)/I_p_wave_areas_1_left_atrial_appendage_tot_number_of_indices;
I_p_wave_areas_1_left_atrial_appendage_quotient_within_3SD = (I_p_wave_areas_1_left_atrial_appendage_dev_avg_pls_3SD_index - I_p_wave_areas_1_left_atrial_appendage_dev_avg_min_3SD_index)/I_p_wave_areas_1_left_atrial_appendage_tot_number_of_indices;
if (I_p_wave_areas_1_left_atrial_appendage_quotient_within_1SD > 0.66 && I_p_wave_areas_1_left_atrial_appendage_quotient_within_1SD < 0.70) && (I_p_wave_areas_1_left_atrial_appendage_quotient_within_2SD > 0.93 && I_p_wave_areas_1_left_atrial_appendage_quotient_within_2SD < 0.97) && (I_p_wave_areas_1_left_atrial_appendage_quotient_within_3SD > 0.98 && I_p_wave_areas_1_left_atrial_appendage_quotient_within_3SD < 1)
    isNormal(7,13) = 1;
end
if I_p_wave_areas_1_left_atrial_appendage_quotient_within_1SD == 0
    isSpike(7,13) = 1;
end
I_p_wave_areas_2_SA_node_sorted = sort(I_p_wave_areas_2_SA_node);
I_p_wave_areas_2_SA_node_dev_avg_min_1SD = area_2(length(I_p_wave_areas_2_SA_node_sorted),1);
for i = 1:length(I_p_wave_areas_2_SA_node_sorted)
    I_p_wave_areas_2_SA_node_dev_avg_min_1SD(i) = I_p_wave_areas_2_SA_node_sorted(i) - (I_p_wave_areas_2_SA_node_average-I_p_wave_areas_2_SA_node_SD);
end
I_p_wave_areas_2_SA_node_dev_avg_min_1SD = abs(I_p_wave_areas_2_SA_node_dev_avg_min_1SD);
for i = 1:length(I_p_wave_areas_2_SA_node_sorted)
    if I_p_wave_areas_2_SA_node_dev_avg_min_1SD(i) == min(I_p_wave_areas_2_SA_node_dev_avg_min_1SD)
        I_p_wave_areas_2_SA_node_dev_avg_min_1SD_index = i;
    end
end
I_p_wave_areas_2_SA_node_dev_avg_pls_1SD = area_2(length(I_p_wave_areas_2_SA_node_sorted),1);
for i = 1:length(I_p_wave_areas_2_SA_node_sorted)
    I_p_wave_areas_2_SA_node_dev_avg_pls_1SD(i) = I_p_wave_areas_2_SA_node_sorted(i) - (I_p_wave_areas_2_SA_node_average+I_p_wave_areas_2_SA_node_SD);
end
I_p_wave_areas_2_SA_node_dev_avg_pls_1SD = abs(I_p_wave_areas_2_SA_node_dev_avg_pls_1SD);
for i = 1:length(I_p_wave_areas_2_SA_node_sorted)
    if I_p_wave_areas_2_SA_node_dev_avg_pls_1SD(i) == min(I_p_wave_areas_2_SA_node_dev_avg_pls_1SD)
        I_p_wave_areas_2_SA_node_dev_avg_pls_1SD_index = i;
    end
end
I_p_wave_areas_2_SA_node_dev_avg_min_2SD = area_2(length(I_p_wave_areas_2_SA_node_sorted),1);
for i = 1:length(I_p_wave_areas_2_SA_node_sorted)
    I_p_wave_areas_2_SA_node_dev_avg_min_2SD(i) = I_p_wave_areas_2_SA_node_sorted(i) - (I_p_wave_areas_2_SA_node_average-2*I_p_wave_areas_2_SA_node_SD);
end
I_p_wave_areas_2_SA_node_dev_avg_min_2SD = abs(I_p_wave_areas_2_SA_node_dev_avg_min_2SD);
for i = 1:length(I_p_wave_areas_2_SA_node_sorted)
    if I_p_wave_areas_2_SA_node_dev_avg_min_2SD(i) == min(I_p_wave_areas_2_SA_node_dev_avg_min_2SD)
        I_p_wave_areas_2_dev_avg_min_2SD_index = i;
    end
end
I_p_wave_areas_2_SA_node_dev_avg_pls_2SD = area_2(length(I_p_wave_areas_2_SA_node_sorted),1);
for i = 1:length(I_p_wave_areas_2_SA_node_sorted)
    I_p_wave_areas_2_SA_node_dev_avg_pls_2SD(i) = I_p_wave_areas_2_SA_node_sorted(i) - (I_p_wave_areas_2_SA_node_average+2*I_p_wave_areas_2_SA_node_SD);
end
I_p_wave_areas_2_SA_node_dev_avg_pls_2SD = abs(I_p_wave_areas_2_SA_node_dev_avg_pls_2SD);
for i = 1:length(I_p_wave_areas_2_SA_node_sorted)
    if I_p_wave_areas_2_SA_node_dev_avg_pls_2SD(i) == min(I_p_wave_areas_2_SA_node_dev_avg_pls_2SD)
        I_p_wave_areas_2_SA_node_dev_avg_pls_2SD_index = i;
    end
end
I_p_wave_areas_2_SA_node_dev_avg_min_3SD = area_2(length(I_p_wave_areas_2_SA_node_sorted),1);
for i = 1:length(I_p_wave_areas_2_SA_node_sorted)
    I_p_wave_areas_2_SA_node_dev_avg_min_3SD(i) = I_p_wave_areas_2_SA_node_sorted(i) - (I_p_wave_areas_2_SA_node_average-3*I_p_wave_areas_2_SA_node_SD);
end
I_p_wave_areas_2_SA_node_dev_avg_min_3SD = abs(I_p_wave_areas_2_SA_node_dev_avg_min_3SD);
for i = 1:length(I_p_wave_areas_2_SA_node_sorted)
    if I_p_wave_areas_2_SA_node_dev_avg_min_3SD(i) == min(I_p_wave_areas_2_SA_node_dev_avg_min_3SD)
        I_p_wave_areas_2_SA_node_dev_avg_min_3SD_index = i;
    end
end
I_p_wave_areas_2_SA_node_dev_avg_pls_3SD = area_2(length(I_p_wave_areas_2_SA_node_sorted),1);
for i = 1:length(I_p_wave_areas_2_SA_node_sorted)
    I_p_wave_areas_2_SA_node_dev_avg_pls_3SD(i) = I_p_wave_areas_2_SA_node_sorted(i) - (I_p_wave_areas_2_SA_node_average+3*I_p_wave_areas_2_SA_node_SD);
end
I_p_wave_areas_2_SA_node_dev_avg_pls_3SD = abs(I_p_wave_areas_2_SA_node_dev_avg_pls_3SD);
for i = 1:length(I_p_wave_areas_2_SA_node_sorted)
    if I_p_wave_areas_2_SA_node_dev_avg_pls_3SD(i) == min(I_p_wave_areas_2_SA_node_dev_avg_pls_3SD)
        I_p_wave_areas_2_SA_node_dev_avg_pls_3SD_index = i;
    end
end
I_p_wave_areas_2_SA_node_tot_number_of_indices = length(I_p_wave_areas_2_SA_node_sorted);
%This is the approximation of the tot "area under curve." 
I_p_wave_areas_2_SA_node_quotient_within_1SD = (I_p_wave_areas_2_SA_node_dev_avg_pls_1SD_index - I_p_wave_areas_2_SA_nodedev_avg_min_1SD_index)/I_p_wave_areas_2_SA_node_tot_number_of_indices;
I_p_wave_areas_2_SA_node_quotient_within_2SD = (I_p_wave_areas_2_SA_node_dev_avg_pls_2SD_index - I_p_wave_areas_2_dev_avg_min_2SD_index)/I_p_wave_areas_2_SA_node_tot_number_of_indices;
I_p_wave_areas_2_SA_node_quotient_within_3SD = (I_p_wave_areas_2_SA_node_dev_avg_pls_3SD_index - I_p_wave_areas_2_SA_node_dev_avg_min_3SD_index)/I_p_wave_areas_2_SA_node_tot_number_of_indices;
if ((I_p_wave_areas_2_SA_node_quotient_within_1SD > 0.66 && I_p_wave_areas_2_SA_node_quotient_within_1SD < 0.70) && (I_p_wave_areas_2_SA_node_quotient_within_2SD > 0.93 && I_p_wave_areas_2_SA_node_quotient_within_2SD < 0.97) && (I_p_wave_areas_2_SA_node_quotient_within_3SD > 0.98 && I_p_wave_areas_2_SA_node_quotient_within_3SD < 1)) 
    isNormal(8,1) = 1;
end
if I_p_wave_areas_2_SA_node_quotient_within_1SD == 0
    isSpike(8,1) = 1;
end
I_p_wave_areas_2_crista_terminalis_sorted = sort(I_p_wave_areas_2_crista_terminalis);
I_p_wave_areas_2_crista_terminalis_dev_avg_min_1SD = area_2(length(I_p_wave_areas_2_crista_terminalis_sorted),1);
for i = 1:length(I_p_wave_areas_2_crista_terminalis_sorted)
    I_p_wave_areas_2_crista_terminalis_dev_avg_min_1SD(i) = I_p_wave_areas_2_crista_terminalis_sorted(i) - (I_p_wave_areas_2_crista_terminalis_average-I_p_wave_areas_2_crista_terminalis_SD);
end
I_p_wave_areas_2_crista_terminalis_dev_avg_min_1SD = abs(I_p_wave_areas_2_crista_terminalis_dev_avg_min_1SD);
for i = 1:length(I_p_wave_areas_2_crista_terminalis_sorted)
    if I_p_wave_areas_2_crista_terminalis_dev_avg_min_1SD(i) == min(I_p_wave_areas_2_crista_terminalis_dev_avg_min_1SD)
        I_p_wave_areas_2_crista_terminalis_dev_avg_min_1SD_index = i;
    end
end
I_p_wave_areas_2_crista_terminalis_dev_avg_pls_1SD = area_2(length(I_p_wave_areas_2_crista_terminalis_sorted),1);
for i = 1:length(I_p_wave_areas_2_crista_terminalis_sorted)
    I_p_wave_areas_2_crista_terminalis_dev_avg_pls_1SD(i) = I_p_wave_areas_2_crista_terminalis_sorted(i) - (I_p_wave_areas_2_crista_terminalis_average+I_p_wave_areas_2_crista_terminalis_SD);
end
I_p_wave_areas_2_crista_terminalis_dev_avg_pls_1SD = abs(I_p_wave_areas_2_crista_terminalis_dev_avg_pls_1SD);
for i = 1:length(I_p_wave_areas_2_crista_terminalis_sorted)
    if I_p_wave_areas_2_crista_terminalis_dev_avg_pls_1SD(i) == min(I_p_wave_areas_2_crista_terminalis_dev_avg_pls_1SD)
        I_p_wave_areas_2_crista_terminalis_dev_avg_pls_1SD_index = i;
    end
end
I_p_wave_areas_2_crista_terminalis_dev_avg_min_2SD = area_2(length(I_p_wave_areas_2_crista_terminalis_sorted),1);
for i = 1:length(I_p_wave_areas_2_crista_terminalis_sorted)
    I_p_wave_areas_2_crista_terminalis_dev_avg_min_2SD(i) = I_p_wave_areas_2_crista_terminalis_sorted(i) - (I_p_wave_areas_2_crista_terminalis_average-2*I_p_wave_areas_2_crista_terminalis_SD);
end
I_p_wave_areas_2_crista_terminalis_dev_avg_min_2SD = abs(I_p_wave_areas_2_crista_terminalis_dev_avg_min_2SD);
for i = 1:length(I_p_wave_areas_2_crista_terminalis_sorted)
    if I_p_wave_areas_2_crista_terminalis_dev_avg_min_2SD(i) == min(I_p_wave_areas_2_crista_terminalis_dev_avg_min_2SD)
        I_p_wave_areas_2_dev_avg_min_2SD_index = i;
    end
end
I_p_wave_areas_2_crista_terminalis_dev_avg_pls_2SD = area_2(length(I_p_wave_areas_2_crista_terminalis_sorted),1);
for i = 1:length(I_p_wave_areas_2_crista_terminalis_sorted)
    I_p_wave_areas_2_crista_terminalis_dev_avg_pls_2SD(i) = I_p_wave_areas_2_crista_terminalis_sorted(i) - (I_p_wave_areas_2_crista_terminalis_average+2*I_p_wave_areas_2_crista_terminalis_SD);
end
I_p_wave_areas_2_crista_terminalis_dev_avg_pls_2SD = abs(I_p_wave_areas_2_crista_terminalis_dev_avg_pls_2SD);
for i = 1:length(I_p_wave_areas_2_crista_terminalis_sorted)
    if I_p_wave_areas_2_crista_terminalis_dev_avg_pls_2SD(i) == min(I_p_wave_areas_2_crista_terminalis_dev_avg_pls_2SD)
        I_p_wave_areas_2_crista_terminalis_dev_avg_pls_2SD_index = i;
    end
end
I_p_wave_areas_2_crista_terminalis_dev_avg_min_3SD = area_2(length(I_p_wave_areas_2_crista_terminalis_sorted),1);
for i = 1:length(I_p_wave_areas_2_crista_terminalis_sorted)
    I_p_wave_areas_2_crista_terminalis_dev_avg_min_3SD(i) = I_p_wave_areas_2_crista_terminalis_sorted(i) - (I_p_wave_areas_2_crista_terminalis_average-3*I_p_wave_areas_2_crista_terminalis_SD);
end
I_p_wave_areas_2_crista_terminalis_dev_avg_min_3SD = abs(I_p_wave_areas_2_crista_terminalis_dev_avg_min_3SD);
for i = 1:length(I_p_wave_areas_2_crista_terminalis_sorted)
    if I_p_wave_areas_2_crista_terminalis_dev_avg_min_3SD(i) == min(I_p_wave_areas_2_crista_terminalis_dev_avg_min_3SD)
        I_p_wave_areas_2_crista_terminalis_dev_avg_min_3SD_index = i;
    end
end
I_p_wave_areas_2_crista_terminalis_dev_avg_pls_3SD = area_2(length(I_p_wave_areas_2_crista_terminalis_sorted),1);
for i = 1:length(I_p_wave_areas_2_crista_terminalis_sorted)
    I_p_wave_areas_2_crista_terminalis_dev_avg_pls_3SD(i) = I_p_wave_areas_2_crista_terminalis_sorted(i) - (I_p_wave_areas_2_crista_terminalis_average+3*I_p_wave_areas_2_crista_terminalis_SD);
end
I_p_wave_areas_2_crista_terminalis_dev_avg_pls_3SD = abs(I_p_wave_areas_2_crista_terminalis_dev_avg_pls_3SD);
for i = 1:length(I_p_wave_areas_2_crista_terminalis_sorted)
    if I_p_wave_areas_2_crista_terminalis_dev_avg_pls_3SD(i) == min(I_p_wave_areas_2_crista_terminalis_dev_avg_pls_3SD)
        I_p_wave_areas_2_crista_terminalis_dev_avg_pls_3SD_index = i;
    end
end
I_p_wave_areas_2_crista_terminalis_tot_number_of_indices = length(I_p_wave_areas_2_crista_terminalis_sorted);
%This is the approximation of the tot "area under curve." 
I_p_wave_areas_2_crista_terminalis_quotient_within_1SD = (I_p_wave_areas_2_crista_terminalis_dev_avg_pls_1SD_index - I_p_wave_areas_2_crista_terminalisdev_avg_min_1SD_index)/I_p_wave_areas_2_crista_terminalis_tot_number_of_indices;
I_p_wave_areas_2_crista_terminalis_quotient_within_2SD = (I_p_wave_areas_2_crista_terminalis_dev_avg_pls_2SD_index - I_p_wave_areas_2_dev_avg_min_2SD_index)/I_p_wave_areas_2_crista_terminalis_tot_number_of_indices;
I_p_wave_areas_2_crista_terminalis_quotient_within_3SD = (I_p_wave_areas_2_crista_terminalis_dev_avg_pls_3SD_index - I_p_wave_areas_2_crista_terminalis_dev_avg_min_3SD_index)/I_p_wave_areas_2_crista_terminalis_tot_number_of_indices;
if (I_p_wave_areas_2_crista_terminalis_quotient_within_1SD > 0.66 && I_p_wave_areas_2_crista_terminalis_quotient_within_1SD < 0.70) && (I_p_wave_areas_2_crista_terminalis_quotient_within_2SD > 0.93 && I_p_wave_areas_2_crista_terminalis_quotient_within_2SD < 0.97) && (I_p_wave_areas_2_crista_terminalis_quotient_within_3SD > 0.98 && I_p_wave_areas_2_crista_terminalis_quotient_within_3SD < 1)
    isNormal(8,2) = 1;
end
if I_p_wave_areas_2_crista_terminalis_quotient_within_1SD == 0
    isSpike(8,2) = 1;
end
I_p_wave_areas_2_tricuspid_annulus_sorted = sort(I_p_wave_areas_2_tricuspid_annulus);
I_p_wave_areas_2_tricuspid_annulus_dev_avg_min_1SD = area_2(length(I_p_wave_areas_2_tricuspid_annulus_sorted),1);
for i = 1:length(I_p_wave_areas_2_tricuspid_annulus_sorted)
    I_p_wave_areas_2_tricuspid_annulus_dev_avg_min_1SD(i) = I_p_wave_areas_2_tricuspid_annulus_sorted(i) - (I_p_wave_areas_2_tricuspid_annulus_average-I_p_wave_areas_2_tricuspid_annulus_SD);
end
I_p_wave_areas_2_tricuspid_annulus_dev_avg_min_1SD = abs(I_p_wave_areas_2_tricuspid_annulus_dev_avg_min_1SD);
for i = 1:length(I_p_wave_areas_2_tricuspid_annulus_sorted)
    if I_p_wave_areas_2_tricuspid_annulus_dev_avg_min_1SD(i) == min(I_p_wave_areas_2_tricuspid_annulus_dev_avg_min_1SD)
        I_p_wave_areas_2_tricuspid_annulus_dev_avg_min_1SD_index = i;
    end
end
I_p_wave_areas_2_tricuspid_annulus_dev_avg_pls_1SD = area_2(length(I_p_wave_areas_2_tricuspid_annulus_sorted),1);
for i = 1:length(I_p_wave_areas_2_tricuspid_annulus_sorted)
    I_p_wave_areas_2_tricuspid_annulus_dev_avg_pls_1SD(i) = I_p_wave_areas_2_tricuspid_annulus_sorted(i) - (I_p_wave_areas_2_tricuspid_annulus_average+I_p_wave_areas_2_tricuspid_annulus_SD);
end
I_p_wave_areas_2_tricuspid_annulus_dev_avg_pls_1SD = abs(I_p_wave_areas_2_tricuspid_annulus_dev_avg_pls_1SD);
for i = 1:length(I_p_wave_areas_2_tricuspid_annulus_sorted)
    if I_p_wave_areas_2_tricuspid_annulus_dev_avg_pls_1SD(i) == min(I_p_wave_areas_2_tricuspid_annulus_dev_avg_pls_1SD)
        I_p_wave_areas_2_tricuspid_annulus_dev_avg_pls_1SD_index = i;
    end
end
I_p_wave_areas_2_tricuspid_annulus_dev_avg_min_2SD = area_2(length(I_p_wave_areas_2_tricuspid_annulus_sorted),1);
for i = 1:length(I_p_wave_areas_2_tricuspid_annulus_sorted)
    I_p_wave_areas_2_tricuspid_annulus_dev_avg_min_2SD(i) = I_p_wave_areas_2_tricuspid_annulus_sorted(i) - (I_p_wave_areas_2_tricuspid_annulus_average-2*I_p_wave_areas_2_tricuspid_annulus_SD);
end
I_p_wave_areas_2_tricuspid_annulus_dev_avg_min_2SD = abs(I_p_wave_areas_2_tricuspid_annulus_dev_avg_min_2SD);
for i = 1:length(I_p_wave_areas_2_tricuspid_annulus_sorted)
    if I_p_wave_areas_2_tricuspid_annulus_dev_avg_min_2SD(i) == min(I_p_wave_areas_2_tricuspid_annulus_dev_avg_min_2SD)
        I_p_wave_areas_2_dev_avg_min_2SD_index = i;
    end
end
I_p_wave_areas_2_tricuspid_annulus_dev_avg_pls_2SD = area_2(length(I_p_wave_areas_2_tricuspid_annulus_sorted),1);
for i = 1:length(I_p_wave_areas_2_tricuspid_annulus_sorted)
    I_p_wave_areas_2_tricuspid_annulus_dev_avg_pls_2SD(i) = I_p_wave_areas_2_tricuspid_annulus_sorted(i) - (I_p_wave_areas_2_tricuspid_annulus_average+2*I_p_wave_areas_2_tricuspid_annulus_SD);
end
I_p_wave_areas_2_tricuspid_annulus_dev_avg_pls_2SD = abs(I_p_wave_areas_2_tricuspid_annulus_dev_avg_pls_2SD);
for i = 1:length(I_p_wave_areas_2_tricuspid_annulus_sorted)
    if I_p_wave_areas_2_tricuspid_annulus_dev_avg_pls_2SD(i) == min(I_p_wave_areas_2_tricuspid_annulus_dev_avg_pls_2SD)
        I_p_wave_areas_2_tricuspid_annulus_dev_avg_pls_2SD_index = i;
    end
end
I_p_wave_areas_2_tricuspid_annulus_dev_avg_min_3SD = area_2(length(I_p_wave_areas_2_tricuspid_annulus_sorted),1);
for i = 1:length(I_p_wave_areas_2_tricuspid_annulus_sorted)
    I_p_wave_areas_2_tricuspid_annulus_dev_avg_min_3SD(i) = I_p_wave_areas_2_tricuspid_annulus_sorted(i) - (I_p_wave_areas_2_tricuspid_annulus_average-3*I_p_wave_areas_2_tricuspid_annulus_SD);
end
I_p_wave_areas_2_tricuspid_annulus_dev_avg_min_3SD = abs(I_p_wave_areas_2_tricuspid_annulus_dev_avg_min_3SD);
for i = 1:length(I_p_wave_areas_2_tricuspid_annulus_sorted)
    if I_p_wave_areas_2_tricuspid_annulus_dev_avg_min_3SD(i) == min(I_p_wave_areas_2_tricuspid_annulus_dev_avg_min_3SD)
        I_p_wave_areas_2_tricuspid_annulus_dev_avg_min_3SD_index = i;
    end
end
I_p_wave_areas_2_tricuspid_annulus_dev_avg_pls_3SD = area_2(length(I_p_wave_areas_2_tricuspid_annulus_sorted),1);
for i = 1:length(I_p_wave_areas_2_tricuspid_annulus_sorted)
    I_p_wave_areas_2_tricuspid_annulus_dev_avg_pls_3SD(i) = I_p_wave_areas_2_tricuspid_annulus_sorted(i) - (I_p_wave_areas_2_tricuspid_annulus_average+3*I_p_wave_areas_2_tricuspid_annulus_SD);
end
I_p_wave_areas_2_tricuspid_annulus_dev_avg_pls_3SD = abs(I_p_wave_areas_2_tricuspid_annulus_dev_avg_pls_3SD);
for i = 1:length(I_p_wave_areas_2_tricuspid_annulus_sorted)
    if I_p_wave_areas_2_tricuspid_annulus_dev_avg_pls_3SD(i) == min(I_p_wave_areas_2_tricuspid_annulus_dev_avg_pls_3SD)
        I_p_wave_areas_2_tricuspid_annulus_dev_avg_pls_3SD_index = i;
    end
end
I_p_wave_areas_2_tricuspid_annulus_tot_number_of_indices = length(I_p_wave_areas_2_tricuspid_annulus_sorted);
%This is the approximation of the tot "area under curve." 
I_p_wave_areas_2_tricuspid_annulus_quotient_within_1SD = (I_p_wave_areas_2_tricuspid_annulus_dev_avg_pls_1SD_index - I_p_wave_areas_2_tricuspid_annulusdev_avg_min_1SD_index)/I_p_wave_areas_2_tricuspid_annulus_tot_number_of_indices;
I_p_wave_areas_2_tricuspid_annulus_quotient_within_2SD = (I_p_wave_areas_2_tricuspid_annulus_dev_avg_pls_2SD_index - I_p_wave_areas_2_dev_avg_min_2SD_index)/I_p_wave_areas_2_tricuspid_annulus_tot_number_of_indices;
I_p_wave_areas_2_tricuspid_annulus_quotient_within_3SD = (I_p_wave_areas_2_tricuspid_annulus_dev_avg_pls_3SD_index - I_p_wave_areas_2_tricuspid_annulus_dev_avg_min_3SD_index)/I_p_wave_areas_2_tricuspid_annulus_tot_number_of_indices;
if (I_p_wave_areas_2_tricuspid_annulus_quotient_within_1SD > 0.66 && I_p_wave_areas_2_tricuspid_annulus_quotient_within_1SD < 0.70) && (I_p_wave_areas_2_tricuspid_annulus_quotient_within_2SD > 0.93 && I_p_wave_areas_2_tricuspid_annulus_quotient_within_2SD < 0.97) && (I_p_wave_areas_2_tricuspid_annulus_quotient_within_3SD > 0.98 && I_p_wave_areas_2_tricuspid_annulus_quotient_within_3SD < 1)
    isNormal(8,3) = 1;
end
if I_p_wave_areas_2_tricuspid_annulus_quotient_within_1SD == 0
    isSpike(8,3) = 1;
end
I_p_wave_areas_2_coronary_sinus_sorted = sort(I_p_wave_areas_2_coronary_sinus);
I_p_wave_areas_2_coronary_sinus_dev_avg_min_1SD = area_2(length(I_p_wave_areas_2_coronary_sinus_sorted),1);
for i = 1:length(I_p_wave_areas_2_coronary_sinus_sorted)
    I_p_wave_areas_2_coronary_sinus_dev_avg_min_1SD(i) = I_p_wave_areas_2_coronary_sinus_sorted(i) - (I_p_wave_areas_2_coronary_sinus_average-I_p_wave_areas_2_coronary_sinus_SD);
end
I_p_wave_areas_2_coronary_sinus_dev_avg_min_1SD = abs(I_p_wave_areas_2_coronary_sinus_dev_avg_min_1SD);
for i = 1:length(I_p_wave_areas_2_coronary_sinus_sorted)
    if I_p_wave_areas_2_coronary_sinus_dev_avg_min_1SD(i) == min(I_p_wave_areas_2_coronary_sinus_dev_avg_min_1SD)
        I_p_wave_areas_2_coronary_sinus_dev_avg_min_1SD_index = i;
    end
end
I_p_wave_areas_2_coronary_sinus_dev_avg_pls_1SD = area_2(length(I_p_wave_areas_2_coronary_sinus_sorted),1);
for i = 1:length(I_p_wave_areas_2_coronary_sinus_sorted)
    I_p_wave_areas_2_coronary_sinus_dev_avg_pls_1SD(i) = I_p_wave_areas_2_coronary_sinus_sorted(i) - (I_p_wave_areas_2_coronary_sinus_average+I_p_wave_areas_2_coronary_sinus_SD);
end
I_p_wave_areas_2_coronary_sinus_dev_avg_pls_1SD = abs(I_p_wave_areas_2_coronary_sinus_dev_avg_pls_1SD);
for i = 1:length(I_p_wave_areas_2_coronary_sinus_sorted)
    if I_p_wave_areas_2_coronary_sinus_dev_avg_pls_1SD(i) == min(I_p_wave_areas_2_coronary_sinus_dev_avg_pls_1SD)
        I_p_wave_areas_2_coronary_sinus_dev_avg_pls_1SD_index = i;
    end
end
I_p_wave_areas_2_coronary_sinus_dev_avg_min_2SD = area_2(length(I_p_wave_areas_2_coronary_sinus_sorted),1);
for i = 1:length(I_p_wave_areas_2_coronary_sinus_sorted)
    I_p_wave_areas_2_coronary_sinus_dev_avg_min_2SD(i) = I_p_wave_areas_2_coronary_sinus_sorted(i) - (I_p_wave_areas_2_coronary_sinus_average-2*I_p_wave_areas_2_coronary_sinus_SD);
end
I_p_wave_areas_2_coronary_sinus_dev_avg_min_2SD = abs(I_p_wave_areas_2_coronary_sinus_dev_avg_min_2SD);
for i = 1:length(I_p_wave_areas_2_coronary_sinus_sorted)
    if I_p_wave_areas_2_coronary_sinus_dev_avg_min_2SD(i) == min(I_p_wave_areas_2_coronary_sinus_dev_avg_min_2SD)
        I_p_wave_areas_2_dev_avg_min_2SD_index = i;
    end
end
I_p_wave_areas_2_coronary_sinus_dev_avg_pls_2SD = area_2(length(I_p_wave_areas_2_coronary_sinus_sorted),1);
for i = 1:length(I_p_wave_areas_2_coronary_sinus_sorted)
    I_p_wave_areas_2_coronary_sinus_dev_avg_pls_2SD(i) = I_p_wave_areas_2_coronary_sinus_sorted(i) - (I_p_wave_areas_2_coronary_sinus_average+2*I_p_wave_areas_2_coronary_sinus_SD);
end
I_p_wave_areas_2_coronary_sinus_dev_avg_pls_2SD = abs(I_p_wave_areas_2_coronary_sinus_dev_avg_pls_2SD);
for i = 1:length(I_p_wave_areas_2_coronary_sinus_sorted)
    if I_p_wave_areas_2_coronary_sinus_dev_avg_pls_2SD(i) == min(I_p_wave_areas_2_coronary_sinus_dev_avg_pls_2SD)
        I_p_wave_areas_2_coronary_sinus_dev_avg_pls_2SD_index = i;
    end
end
I_p_wave_areas_2_coronary_sinus_dev_avg_min_3SD = area_2(length(I_p_wave_areas_2_coronary_sinus_sorted),1);
for i = 1:length(I_p_wave_areas_2_coronary_sinus_sorted)
    I_p_wave_areas_2_coronary_sinus_dev_avg_min_3SD(i) = I_p_wave_areas_2_coronary_sinus_sorted(i) - (I_p_wave_areas_2_coronary_sinus_average-3*I_p_wave_areas_2_coronary_sinus_SD);
end
I_p_wave_areas_2_coronary_sinus_dev_avg_min_3SD = abs(I_p_wave_areas_2_coronary_sinus_dev_avg_min_3SD);
for i = 1:length(I_p_wave_areas_2_coronary_sinus_sorted)
    if I_p_wave_areas_2_coronary_sinus_dev_avg_min_3SD(i) == min(I_p_wave_areas_2_coronary_sinus_dev_avg_min_3SD)
        I_p_wave_areas_2_coronary_sinus_dev_avg_min_3SD_index = i;
    end
end
I_p_wave_areas_2_coronary_sinus_dev_avg_pls_3SD = area_2(length(I_p_wave_areas_2_coronary_sinus_sorted),1);
for i = 1:length(I_p_wave_areas_2_coronary_sinus_sorted)
    I_p_wave_areas_2_coronary_sinus_dev_avg_pls_3SD(i) = I_p_wave_areas_2_coronary_sinus_sorted(i) - (I_p_wave_areas_2_coronary_sinus_average+3*I_p_wave_areas_2_coronary_sinus_SD);
end
I_p_wave_areas_2_coronary_sinus_dev_avg_pls_3SD = abs(I_p_wave_areas_2_coronary_sinus_dev_avg_pls_3SD);
for i = 1:length(I_p_wave_areas_2_coronary_sinus_sorted)
    if I_p_wave_areas_2_coronary_sinus_dev_avg_pls_3SD(i) == min(I_p_wave_areas_2_coronary_sinus_dev_avg_pls_3SD)
        I_p_wave_areas_2_coronary_sinus_dev_avg_pls_3SD_index = i;
    end
end
I_p_wave_areas_2_coronary_sinus_tot_number_of_indices = length(I_p_wave_areas_2_coronary_sinus_sorted);
%This is the approximation of the tot "area under curve." 
I_p_wave_areas_2_coronary_sinus_quotient_within_1SD = (I_p_wave_areas_2_coronary_sinus_dev_avg_pls_1SD_index - I_p_wave_areas_2_coronary_sinusdev_avg_min_1SD_index)/I_p_wave_areas_2_coronary_sinus_tot_number_of_indices;
I_p_wave_areas_2_coronary_sinus_quotient_within_2SD = (I_p_wave_areas_2_coronary_sinus_dev_avg_pls_2SD_index - I_p_wave_areas_2_dev_avg_min_2SD_index)/I_p_wave_areas_2_coronary_sinus_tot_number_of_indices;
I_p_wave_areas_2_coronary_sinus_quotient_within_3SD = (I_p_wave_areas_2_coronary_sinus_dev_avg_pls_3SD_index - I_p_wave_areas_2_coronary_sinus_dev_avg_min_3SD_index)/I_p_wave_areas_2_coronary_sinus_tot_number_of_indices;
if (I_p_wave_areas_2_coronary_sinus_quotient_within_1SD > 0.66 && I_p_wave_areas_2_coronary_sinus_quotient_within_1SD < 0.70) && (I_p_wave_areas_2_coronary_sinus_quotient_within_2SD > 0.93 && I_p_wave_areas_2_coronary_sinus_quotient_within_2SD < 0.97) && (I_p_wave_areas_2_coronary_sinus_quotient_within_3SD > 0.98 && I_p_wave_areas_2_coronary_sinus_quotient_within_3SD < 1)
    isNormal(8,4) = 1;
end
if I_p_wave_areas_2_coronary_sinus_quotient_within_1SD == 0
    isSpike(8,4) = 1;
end
I_p_wave_areas_2_ostium_sorted = sort(I_p_wave_areas_2_ostium);
I_p_wave_areas_2_ostium_dev_avg_min_1SD = area_2(length(I_p_wave_areas_2_ostium_sorted),1);
for i = 1:length(I_p_wave_areas_2_ostium_sorted)
    I_p_wave_areas_2_ostium_dev_avg_min_1SD(i) = I_p_wave_areas_2_ostium_sorted(i) - (I_p_wave_areas_2_ostium_average-I_p_wave_areas_2_ostium_SD);
end
I_p_wave_areas_2_ostium_dev_avg_min_1SD = abs(I_p_wave_areas_2_ostium_dev_avg_min_1SD);
for i = 1:length(I_p_wave_areas_2_ostium_sorted)
    if I_p_wave_areas_2_ostium_dev_avg_min_1SD(i) == min(I_p_wave_areas_2_ostium_dev_avg_min_1SD)
        I_p_wave_areas_2_ostium_dev_avg_min_1SD_index = i;
    end
end
I_p_wave_areas_2_ostium_dev_avg_pls_1SD = area_2(length(I_p_wave_areas_2_ostium_sorted),1);
for i = 1:length(I_p_wave_areas_2_ostium_sorted)
    I_p_wave_areas_2_ostium_dev_avg_pls_1SD(i) = I_p_wave_areas_2_ostium_sorted(i) - (I_p_wave_areas_2_ostium_average+I_p_wave_areas_2_ostium_SD);
end
I_p_wave_areas_2_ostium_dev_avg_pls_1SD = abs(I_p_wave_areas_2_ostium_dev_avg_pls_1SD);
for i = 1:length(I_p_wave_areas_2_ostium_sorted)
    if I_p_wave_areas_2_ostium_dev_avg_pls_1SD(i) == min(I_p_wave_areas_2_ostium_dev_avg_pls_1SD)
        I_p_wave_areas_2_ostium_dev_avg_pls_1SD_index = i;
    end
end
I_p_wave_areas_2_ostium_dev_avg_min_2SD = area_2(length(I_p_wave_areas_2_ostium_sorted),1);
for i = 1:length(I_p_wave_areas_2_ostium_sorted)
    I_p_wave_areas_2_ostium_dev_avg_min_2SD(i) = I_p_wave_areas_2_ostium_sorted(i) - (I_p_wave_areas_2_ostium_average-2*I_p_wave_areas_2_ostium_SD);
end
I_p_wave_areas_2_ostium_dev_avg_min_2SD = abs(I_p_wave_areas_2_ostium_dev_avg_min_2SD);
for i = 1:length(I_p_wave_areas_2_ostium_sorted)
    if I_p_wave_areas_2_ostium_dev_avg_min_2SD(i) == min(I_p_wave_areas_2_ostium_dev_avg_min_2SD)
        I_p_wave_areas_2_dev_avg_min_2SD_index = i;
    end
end
I_p_wave_areas_2_ostium_dev_avg_pls_2SD = area_2(length(I_p_wave_areas_2_ostium_sorted),1);
for i = 1:length(I_p_wave_areas_2_ostium_sorted)
    I_p_wave_areas_2_ostium_dev_avg_pls_2SD(i) = I_p_wave_areas_2_ostium_sorted(i) - (I_p_wave_areas_2_ostium_average+2*I_p_wave_areas_2_ostium_SD);
end
I_p_wave_areas_2_ostium_dev_avg_pls_2SD = abs(I_p_wave_areas_2_ostium_dev_avg_pls_2SD);
for i = 1:length(I_p_wave_areas_2_ostium_sorted)
    if I_p_wave_areas_2_ostium_dev_avg_pls_2SD(i) == min(I_p_wave_areas_2_ostium_dev_avg_pls_2SD)
        I_p_wave_areas_2_ostium_dev_avg_pls_2SD_index = i;
    end
end
I_p_wave_areas_2_ostium_dev_avg_min_3SD = area_2(length(I_p_wave_areas_2_ostium_sorted),1);
for i = 1:length(I_p_wave_areas_2_ostium_sorted)
    I_p_wave_areas_2_ostium_dev_avg_min_3SD(i) = I_p_wave_areas_2_ostium_sorted(i) - (I_p_wave_areas_2_ostium_average-3*I_p_wave_areas_2_ostium_SD);
end
I_p_wave_areas_2_ostium_dev_avg_min_3SD = abs(I_p_wave_areas_2_ostium_dev_avg_min_3SD);
for i = 1:length(I_p_wave_areas_2_ostium_sorted)
    if I_p_wave_areas_2_ostium_dev_avg_min_3SD(i) == min(I_p_wave_areas_2_ostium_dev_avg_min_3SD)
        I_p_wave_areas_2_ostium_dev_avg_min_3SD_index = i;
    end
end
I_p_wave_areas_2_ostium_dev_avg_pls_3SD = area_2(length(I_p_wave_areas_2_ostium_sorted),1);
for i = 1:length(I_p_wave_areas_2_ostium_sorted)
    I_p_wave_areas_2_ostium_dev_avg_pls_3SD(i) = I_p_wave_areas_2_ostium_sorted(i) - (I_p_wave_areas_2_ostium_average+3*I_p_wave_areas_2_ostium_SD);
end
I_p_wave_areas_2_ostium_dev_avg_pls_3SD = abs(I_p_wave_areas_2_ostium_dev_avg_pls_3SD);
for i = 1:length(I_p_wave_areas_2_ostium_sorted)
    if I_p_wave_areas_2_ostium_dev_avg_pls_3SD(i) == min(I_p_wave_areas_2_ostium_dev_avg_pls_3SD)
        I_p_wave_areas_2_ostium_dev_avg_pls_3SD_index = i;
    end
end
I_p_wave_areas_2_ostium_tot_number_of_indices = length(I_p_wave_areas_2_ostium_sorted);
%This is the approximation of the tot "area under curve." 
I_p_wave_areas_2_ostium_quotient_within_1SD = (I_p_wave_areas_2_ostium_dev_avg_pls_1SD_index - I_p_wave_areas_2_ostiumdev_avg_min_1SD_index)/I_p_wave_areas_2_ostium_tot_number_of_indices;
I_p_wave_areas_2_ostium_quotient_within_2SD = (I_p_wave_areas_2_ostium_dev_avg_pls_2SD_index - I_p_wave_areas_2_dev_avg_min_2SD_index)/I_p_wave_areas_2_ostium_tot_number_of_indices;
I_p_wave_areas_2_ostium_quotient_within_3SD = (I_p_wave_areas_2_ostium_dev_avg_pls_3SD_index - I_p_wave_areas_2_ostium_dev_avg_min_3SD_index)/I_p_wave_areas_2_ostium_tot_number_of_indices;
if (I_p_wave_areas_2_ostium_quotient_within_1SD > 0.66 && I_p_wave_areas_2_ostium_quotient_within_1SD < 0.70) && (I_p_wave_areas_2_ostium_quotient_within_2SD > 0.93 && I_p_wave_areas_2_ostium_quotient_within_2SD < 0.97) && (I_p_wave_areas_2_ostium_quotient_within_3SD > 0.98 && I_p_wave_areas_2_ostium_quotient_within_3SD < 1)
    isNormal(8,5) = 1;
end
if I_p_wave_areas_2_ostium_quotient_within_1SD == 0
    isSpike(8,5) = 1;
end
I_p_wave_areas_2_perinodal_sorted = sort(I_p_wave_areas_2_perinodal);
I_p_wave_areas_2_perinodal_dev_avg_min_1SD = area_2(length(I_p_wave_areas_2_perinodal_sorted),1);
for i = 1:length(I_p_wave_areas_2_perinodal_sorted)
    I_p_wave_areas_2_perinodal_dev_avg_min_1SD(i) = I_p_wave_areas_2_perinodal_sorted(i) - (I_p_wave_areas_2_perinodal_average-I_p_wave_areas_2_perinodal_SD);
end
I_p_wave_areas_2_perinodal_dev_avg_min_1SD = abs(I_p_wave_areas_2_perinodal_dev_avg_min_1SD);
for i = 1:length(I_p_wave_areas_2_perinodal_sorted)
    if I_p_wave_areas_2_perinodal_dev_avg_min_1SD(i) == min(I_p_wave_areas_2_perinodal_dev_avg_min_1SD)
        I_p_wave_areas_2_perinodal_dev_avg_min_1SD_index = i;
    end
end
I_p_wave_areas_2_perinodal_dev_avg_pls_1SD = area_2(length(I_p_wave_areas_2_perinodal_sorted),1);
for i = 1:length(I_p_wave_areas_2_perinodal_sorted)
    I_p_wave_areas_2_perinodal_dev_avg_pls_1SD(i) = I_p_wave_areas_2_perinodal_sorted(i) - (I_p_wave_areas_2_perinodal_average+I_p_wave_areas_2_perinodal_SD);
end
I_p_wave_areas_2_perinodal_dev_avg_pls_1SD = abs(I_p_wave_areas_2_perinodal_dev_avg_pls_1SD);
for i = 1:length(I_p_wave_areas_2_perinodal_sorted)
    if I_p_wave_areas_2_perinodal_dev_avg_pls_1SD(i) == min(I_p_wave_areas_2_perinodal_dev_avg_pls_1SD)
        I_p_wave_areas_2_perinodal_dev_avg_pls_1SD_index = i;
    end
end
I_p_wave_areas_2_perinodal_dev_avg_min_2SD = area_2(length(I_p_wave_areas_2_perinodal_sorted),1);
for i = 1:length(I_p_wave_areas_2_perinodal_sorted)
    I_p_wave_areas_2_perinodal_dev_avg_min_2SD(i) = I_p_wave_areas_2_perinodal_sorted(i) - (I_p_wave_areas_2_perinodal_average-2*I_p_wave_areas_2_perinodal_SD);
end
I_p_wave_areas_2_perinodal_dev_avg_min_2SD = abs(I_p_wave_areas_2_perinodal_dev_avg_min_2SD);
for i = 1:length(I_p_wave_areas_2_perinodal_sorted)
    if I_p_wave_areas_2_perinodal_dev_avg_min_2SD(i) == min(I_p_wave_areas_2_perinodal_dev_avg_min_2SD)
        I_p_wave_areas_2_dev_avg_min_2SD_index = i;
    end
end
I_p_wave_areas_2_perinodal_dev_avg_pls_2SD = area_2(length(I_p_wave_areas_2_perinodal_sorted),1);
for i = 1:length(I_p_wave_areas_2_perinodal_sorted)
    I_p_wave_areas_2_perinodal_dev_avg_pls_2SD(i) = I_p_wave_areas_2_perinodal_sorted(i) - (I_p_wave_areas_2_perinodal_average+2*I_p_wave_areas_2_perinodal_SD);
end
I_p_wave_areas_2_perinodal_dev_avg_pls_2SD = abs(I_p_wave_areas_2_perinodal_dev_avg_pls_2SD);
for i = 1:length(I_p_wave_areas_2_perinodal_sorted)
    if I_p_wave_areas_2_perinodal_dev_avg_pls_2SD(i) == min(I_p_wave_areas_2_perinodal_dev_avg_pls_2SD)
        I_p_wave_areas_2_perinodal_dev_avg_pls_2SD_index = i;
    end
end
I_p_wave_areas_2_perinodal_dev_avg_min_3SD = area_2(length(I_p_wave_areas_2_perinodal_sorted),1);
for i = 1:length(I_p_wave_areas_2_perinodal_sorted)
    I_p_wave_areas_2_perinodal_dev_avg_min_3SD(i) = I_p_wave_areas_2_perinodal_sorted(i) - (I_p_wave_areas_2_perinodal_average-3*I_p_wave_areas_2_perinodal_SD);
end
I_p_wave_areas_2_perinodal_dev_avg_min_3SD = abs(I_p_wave_areas_2_perinodal_dev_avg_min_3SD);
for i = 1:length(I_p_wave_areas_2_perinodal_sorted)
    if I_p_wave_areas_2_perinodal_dev_avg_min_3SD(i) == min(I_p_wave_areas_2_perinodal_dev_avg_min_3SD)
        I_p_wave_areas_2_perinodal_dev_avg_min_3SD_index = i;
    end
end
I_p_wave_areas_2_perinodal_dev_avg_pls_3SD = area_2(length(I_p_wave_areas_2_perinodal_sorted),1);
for i = 1:length(I_p_wave_areas_2_perinodal_sorted)
    I_p_wave_areas_2_perinodal_dev_avg_pls_3SD(i) = I_p_wave_areas_2_perinodal_sorted(i) - (I_p_wave_areas_2_perinodal_average+3*I_p_wave_areas_2_perinodal_SD);
end
I_p_wave_areas_2_perinodal_dev_avg_pls_3SD = abs(I_p_wave_areas_2_perinodal_dev_avg_pls_3SD);
for i = 1:length(I_p_wave_areas_2_perinodal_sorted)
    if I_p_wave_areas_2_perinodal_dev_avg_pls_3SD(i) == min(I_p_wave_areas_2_perinodal_dev_avg_pls_3SD)
        I_p_wave_areas_2_perinodal_dev_avg_pls_3SD_index = i;
    end
end
I_p_wave_areas_2_perinodal_tot_number_of_indices = length(I_p_wave_areas_2_perinodal_sorted);
%This is the approximation of the tot "area under curve." 
I_p_wave_areas_2_perinodal_quotient_within_1SD = (I_p_wave_areas_2_perinodal_dev_avg_pls_1SD_index - I_p_wave_areas_2_perinodaldev_avg_min_1SD_index)/I_p_wave_areas_2_perinodal_tot_number_of_indices;
I_p_wave_areas_2_perinodal_quotient_within_2SD = (I_p_wave_areas_2_perinodal_dev_avg_pls_2SD_index - I_p_wave_areas_2_dev_avg_min_2SD_index)/I_p_wave_areas_2_perinodal_tot_number_of_indices;
I_p_wave_areas_2_perinodal_quotient_within_3SD = (I_p_wave_areas_2_perinodal_dev_avg_pls_3SD_index - I_p_wave_areas_2_perinodal_dev_avg_min_3SD_index)/I_p_wave_areas_2_perinodal_tot_number_of_indices;
if (I_p_wave_areas_2_perinodal_quotient_within_1SD > 0.66 && I_p_wave_areas_2_perinodal_quotient_within_1SD < 0.70) && (I_p_wave_areas_2_perinodal_quotient_within_2SD > 0.93 && I_p_wave_areas_2_perinodal_quotient_within_2SD < 0.97) && (I_p_wave_areas_2_perinodal_quotient_within_3SD > 0.98 && I_p_wave_areas_2_perinodal_quotient_within_3SD < 1)
    isNormal(8,6) = 1;
end
if I_p_wave_areas_2_perinodal_quotient_within_1SD == 0
    isSpike(8,6) = 1;
end
I_p_wave_areas_2_right_septum_sorted = sort(I_p_wave_areas_2_right_septum);
I_p_wave_areas_2_right_septum_dev_avg_min_1SD = area_2(length(I_p_wave_areas_2_right_septum_sorted),1);
for i = 1:length(I_p_wave_areas_2_right_septum_sorted)
    I_p_wave_areas_2_right_septum_dev_avg_min_1SD(i) = I_p_wave_areas_2_right_septum_sorted(i) - (I_p_wave_areas_2_right_septum_average-I_p_wave_areas_2_right_septum_SD);
end
I_p_wave_areas_2_right_septum_dev_avg_min_1SD = abs(I_p_wave_areas_2_right_septum_dev_avg_min_1SD);
for i = 1:length(I_p_wave_areas_2_right_septum_sorted)
    if I_p_wave_areas_2_right_septum_dev_avg_min_1SD(i) == min(I_p_wave_areas_2_right_septum_dev_avg_min_1SD)
        I_p_wave_areas_2_right_septum_dev_avg_min_1SD_index = i;
    end
end
I_p_wave_areas_2_right_septum_dev_avg_pls_1SD = area_2(length(I_p_wave_areas_2_right_septum_sorted),1);
for i = 1:length(I_p_wave_areas_2_right_septum_sorted)
    I_p_wave_areas_2_right_septum_dev_avg_pls_1SD(i) = I_p_wave_areas_2_right_septum_sorted(i) - (I_p_wave_areas_2_right_septum_average+I_p_wave_areas_2_right_septum_SD);
end
I_p_wave_areas_2_right_septum_dev_avg_pls_1SD = abs(I_p_wave_areas_2_right_septum_dev_avg_pls_1SD);
for i = 1:length(I_p_wave_areas_2_right_septum_sorted)
    if I_p_wave_areas_2_right_septum_dev_avg_pls_1SD(i) == min(I_p_wave_areas_2_right_septum_dev_avg_pls_1SD)
        I_p_wave_areas_2_right_septum_dev_avg_pls_1SD_index = i;
    end
end
I_p_wave_areas_2_right_septum_dev_avg_min_2SD = area_2(length(I_p_wave_areas_2_right_septum_sorted),1);
for i = 1:length(I_p_wave_areas_2_right_septum_sorted)
    I_p_wave_areas_2_right_septum_dev_avg_min_2SD(i) = I_p_wave_areas_2_right_septum_sorted(i) - (I_p_wave_areas_2_right_septum_average-2*I_p_wave_areas_2_right_septum_SD);
end
I_p_wave_areas_2_right_septum_dev_avg_min_2SD = abs(I_p_wave_areas_2_right_septum_dev_avg_min_2SD);
for i = 1:length(I_p_wave_areas_2_right_septum_sorted)
    if I_p_wave_areas_2_right_septum_dev_avg_min_2SD(i) == min(I_p_wave_areas_2_right_septum_dev_avg_min_2SD)
        I_p_wave_areas_2_dev_avg_min_2SD_index = i;
    end
end
I_p_wave_areas_2_right_septum_dev_avg_pls_2SD = area_2(length(I_p_wave_areas_2_right_septum_sorted),1);
for i = 1:length(I_p_wave_areas_2_right_septum_sorted)
    I_p_wave_areas_2_right_septum_dev_avg_pls_2SD(i) = I_p_wave_areas_2_right_septum_sorted(i) - (I_p_wave_areas_2_right_septum_average+2*I_p_wave_areas_2_right_septum_SD);
end
I_p_wave_areas_2_right_septum_dev_avg_pls_2SD = abs(I_p_wave_areas_2_right_septum_dev_avg_pls_2SD);
for i = 1:length(I_p_wave_areas_2_right_septum_sorted)
    if I_p_wave_areas_2_right_septum_dev_avg_pls_2SD(i) == min(I_p_wave_areas_2_right_septum_dev_avg_pls_2SD)
        I_p_wave_areas_2_right_septum_dev_avg_pls_2SD_index = i;
    end
end
I_p_wave_areas_2_right_septum_dev_avg_min_3SD = area_2(length(I_p_wave_areas_2_right_septum_sorted),1);
for i = 1:length(I_p_wave_areas_2_right_septum_sorted)
    I_p_wave_areas_2_right_septum_dev_avg_min_3SD(i) = I_p_wave_areas_2_right_septum_sorted(i) - (I_p_wave_areas_2_right_septum_average-3*I_p_wave_areas_2_right_septum_SD);
end
I_p_wave_areas_2_right_septum_dev_avg_min_3SD = abs(I_p_wave_areas_2_right_septum_dev_avg_min_3SD);
for i = 1:length(I_p_wave_areas_2_right_septum_sorted)
    if I_p_wave_areas_2_right_septum_dev_avg_min_3SD(i) == min(I_p_wave_areas_2_right_septum_dev_avg_min_3SD)
        I_p_wave_areas_2_right_septum_dev_avg_min_3SD_index = i;
    end
end
I_p_wave_areas_2_right_septum_dev_avg_pls_3SD = area_2(length(I_p_wave_areas_2_right_septum_sorted),1);
for i = 1:length(I_p_wave_areas_2_right_septum_sorted)
    I_p_wave_areas_2_right_septum_dev_avg_pls_3SD(i) = I_p_wave_areas_2_right_septum_sorted(i) - (I_p_wave_areas_2_right_septum_average+3*I_p_wave_areas_2_right_septum_SD);
end
I_p_wave_areas_2_right_septum_dev_avg_pls_3SD = abs(I_p_wave_areas_2_right_septum_dev_avg_pls_3SD);
for i = 1:length(I_p_wave_areas_2_right_septum_sorted)
    if I_p_wave_areas_2_right_septum_dev_avg_pls_3SD(i) == min(I_p_wave_areas_2_right_septum_dev_avg_pls_3SD)
        I_p_wave_areas_2_right_septum_dev_avg_pls_3SD_index = i;
    end
end
I_p_wave_areas_2_right_septum_tot_number_of_indices = length(I_p_wave_areas_2_right_septum_sorted);
%This is the approximation of the tot "area under curve." 
I_p_wave_areas_2_right_septum_quotient_within_1SD = (I_p_wave_areas_2_right_septum_dev_avg_pls_1SD_index - I_p_wave_areas_2_right_septumdev_avg_min_1SD_index)/I_p_wave_areas_2_right_septum_tot_number_of_indices;
I_p_wave_areas_2_right_septum_quotient_within_2SD = (I_p_wave_areas_2_right_septum_dev_avg_pls_2SD_index - I_p_wave_areas_2_dev_avg_min_2SD_index)/I_p_wave_areas_2_right_septum_tot_number_of_indices;
I_p_wave_areas_2_right_septum_quotient_within_3SD = (I_p_wave_areas_2_right_septum_dev_avg_pls_3SD_index - I_p_wave_areas_2_right_septum_dev_avg_min_3SD_index)/I_p_wave_areas_2_right_septum_tot_number_of_indices;
if (I_p_wave_areas_2_right_septum_quotient_within_1SD > 0.66 && I_p_wave_areas_2_right_septum_quotient_within_1SD < 0.70) && (I_p_wave_areas_2_right_septum_quotient_within_2SD > 0.93 && I_p_wave_areas_2_right_septum_quotient_within_2SD < 0.97) && (I_p_wave_areas_2_right_septum_quotient_within_3SD > 0.98 && I_p_wave_areas_2_right_septum_quotient_within_3SD < 1)
    isNormal(8,7) = 1;
end
if I_p_wave_areas_2_right_septum_quotient_within_1SD == 0
    isSpike(8,7) = 1;
end
I_p_wave_areas_2_right_atrial_appendage_sorted = sort(I_p_wave_areas_2_right_atrial_appendage);
I_p_wave_areas_2_right_atrial_appendage_dev_avg_min_1SD = area_2(length(I_p_wave_areas_2_right_atrial_appendage_sorted),1);
for i = 1:length(I_p_wave_areas_2_right_atrial_appendage_sorted)
    I_p_wave_areas_2_right_atrial_appendage_dev_avg_min_1SD(i) = I_p_wave_areas_2_right_atrial_appendage_sorted(i) - (I_p_wave_areas_2_right_atrial_appendage_average-I_p_wave_areas_2_right_atrial_appendage_SD);
end
I_p_wave_areas_2_right_atrial_appendage_dev_avg_min_1SD = abs(I_p_wave_areas_2_right_atrial_appendage_dev_avg_min_1SD);
for i = 1:length(I_p_wave_areas_2_right_atrial_appendage_sorted)
    if I_p_wave_areas_2_right_atrial_appendage_dev_avg_min_1SD(i) == min(I_p_wave_areas_2_right_atrial_appendage_dev_avg_min_1SD)
        I_p_wave_areas_2_right_atrial_appendage_dev_avg_min_1SD_index = i;
    end
end
I_p_wave_areas_2_right_atrial_appendage_dev_avg_pls_1SD = area_2(length(I_p_wave_areas_2_right_atrial_appendage_sorted),1);
for i = 1:length(I_p_wave_areas_2_right_atrial_appendage_sorted)
    I_p_wave_areas_2_right_atrial_appendage_dev_avg_pls_1SD(i) = I_p_wave_areas_2_right_atrial_appendage_sorted(i) - (I_p_wave_areas_2_right_atrial_appendage_average+I_p_wave_areas_2_right_atrial_appendage_SD);
end
I_p_wave_areas_2_right_atrial_appendage_dev_avg_pls_1SD = abs(I_p_wave_areas_2_right_atrial_appendage_dev_avg_pls_1SD);
for i = 1:length(I_p_wave_areas_2_right_atrial_appendage_sorted)
    if I_p_wave_areas_2_right_atrial_appendage_dev_avg_pls_1SD(i) == min(I_p_wave_areas_2_right_atrial_appendage_dev_avg_pls_1SD)
        I_p_wave_areas_2_right_atrial_appendage_dev_avg_pls_1SD_index = i;
    end
end
I_p_wave_areas_2_right_atrial_appendage_dev_avg_min_2SD = area_2(length(I_p_wave_areas_2_right_atrial_appendage_sorted),1);
for i = 1:length(I_p_wave_areas_2_right_atrial_appendage_sorted)
    I_p_wave_areas_2_right_atrial_appendage_dev_avg_min_2SD(i) = I_p_wave_areas_2_right_atrial_appendage_sorted(i) - (I_p_wave_areas_2_right_atrial_appendage_average-2*I_p_wave_areas_2_right_atrial_appendage_SD);
end
I_p_wave_areas_2_right_atrial_appendage_dev_avg_min_2SD = abs(I_p_wave_areas_2_right_atrial_appendage_dev_avg_min_2SD);
for i = 1:length(I_p_wave_areas_2_right_atrial_appendage_sorted)
    if I_p_wave_areas_2_right_atrial_appendage_dev_avg_min_2SD(i) == min(I_p_wave_areas_2_right_atrial_appendage_dev_avg_min_2SD)
        I_p_wave_areas_2_dev_avg_min_2SD_index = i;
    end
end
I_p_wave_areas_2_right_atrial_appendage_dev_avg_pls_2SD = area_2(length(I_p_wave_areas_2_right_atrial_appendage_sorted),1);
for i = 1:length(I_p_wave_areas_2_right_atrial_appendage_sorted)
    I_p_wave_areas_2_right_atrial_appendage_dev_avg_pls_2SD(i) = I_p_wave_areas_2_right_atrial_appendage_sorted(i) - (I_p_wave_areas_2_right_atrial_appendage_average+2*I_p_wave_areas_2_right_atrial_appendage_SD);
end
I_p_wave_areas_2_right_atrial_appendage_dev_avg_pls_2SD = abs(I_p_wave_areas_2_right_atrial_appendage_dev_avg_pls_2SD);
for i = 1:length(I_p_wave_areas_2_right_atrial_appendage_sorted)
    if I_p_wave_areas_2_right_atrial_appendage_dev_avg_pls_2SD(i) == min(I_p_wave_areas_2_right_atrial_appendage_dev_avg_pls_2SD)
        I_p_wave_areas_2_right_atrial_appendage_dev_avg_pls_2SD_index = i;
    end
end
I_p_wave_areas_2_right_atrial_appendage_dev_avg_min_3SD = area_2(length(I_p_wave_areas_2_right_atrial_appendage_sorted),1);
for i = 1:length(I_p_wave_areas_2_right_atrial_appendage_sorted)
    I_p_wave_areas_2_right_atrial_appendage_dev_avg_min_3SD(i) = I_p_wave_areas_2_right_atrial_appendage_sorted(i) - (I_p_wave_areas_2_right_atrial_appendage_average-3*I_p_wave_areas_2_right_atrial_appendage_SD);
end
I_p_wave_areas_2_right_atrial_appendage_dev_avg_min_3SD = abs(I_p_wave_areas_2_right_atrial_appendage_dev_avg_min_3SD);
for i = 1:length(I_p_wave_areas_2_right_atrial_appendage_sorted)
    if I_p_wave_areas_2_right_atrial_appendage_dev_avg_min_3SD(i) == min(I_p_wave_areas_2_right_atrial_appendage_dev_avg_min_3SD)
        I_p_wave_areas_2_right_atrial_appendage_dev_avg_min_3SD_index = i;
    end
end
I_p_wave_areas_2_right_atrial_appendage_dev_avg_pls_3SD = area_2(length(I_p_wave_areas_2_right_atrial_appendage_sorted),1);
for i = 1:length(I_p_wave_areas_2_right_atrial_appendage_sorted)
    I_p_wave_areas_2_right_atrial_appendage_dev_avg_pls_3SD(i) = I_p_wave_areas_2_right_atrial_appendage_sorted(i) - (I_p_wave_areas_2_right_atrial_appendage_average+3*I_p_wave_areas_2_right_atrial_appendage_SD);
end
I_p_wave_areas_2_right_atrial_appendage_dev_avg_pls_3SD = abs(I_p_wave_areas_2_right_atrial_appendage_dev_avg_pls_3SD);
for i = 1:length(I_p_wave_areas_2_right_atrial_appendage_sorted)
    if I_p_wave_areas_2_right_atrial_appendage_dev_avg_pls_3SD(i) == min(I_p_wave_areas_2_right_atrial_appendage_dev_avg_pls_3SD)
        I_p_wave_areas_2_right_atrial_appendage_dev_avg_pls_3SD_index = i;
    end
end
I_p_wave_areas_2_right_atrial_appendage_tot_number_of_indices = length(I_p_wave_areas_2_right_atrial_appendage_sorted);
%This is the approximation of the tot "area under curve." 
I_p_wave_areas_2_right_atrial_appendage_quotient_within_1SD = (I_p_wave_areas_2_right_atrial_appendage_dev_avg_pls_1SD_index - I_p_wave_areas_2_right_atrial_appendagedev_avg_min_1SD_index)/I_p_wave_areas_2_right_atrial_appendage_tot_number_of_indices;
I_p_wave_areas_2_right_atrial_appendage_quotient_within_2SD = (I_p_wave_areas_2_right_atrial_appendage_dev_avg_pls_2SD_index - I_p_wave_areas_2_dev_avg_min_2SD_index)/I_p_wave_areas_2_right_atrial_appendage_tot_number_of_indices;
I_p_wave_areas_2_right_atrial_appendage_quotient_within_3SD = (I_p_wave_areas_2_right_atrial_appendage_dev_avg_pls_3SD_index - I_p_wave_areas_2_right_atrial_appendage_dev_avg_min_3SD_index)/I_p_wave_areas_2_right_atrial_appendage_tot_number_of_indices;
if (I_p_wave_areas_2_right_atrial_appendage_quotient_within_1SD > 0.66 && I_p_wave_areas_2_right_atrial_appendage_quotient_within_1SD < 0.70) && (I_p_wave_areas_2_right_atrial_appendage_quotient_within_2SD > 0.93 && I_p_wave_areas_2_right_atrial_appendage_quotient_within_2SD < 0.97) && (I_p_wave_areas_2_right_atrial_appendage_quotient_within_3SD > 0.98 && I_p_wave_areas_2_right_atrial_appendage_quotient_within_3SD < 1)
    isNormal(8,8) = 1;
end
if I_p_wave_areas_2_right_atrial_appendage_quotient_within_1SD == 0
    isSpike(8,8) = 1;
end
I_p_wave_areas_2_pulmonary_veins_sorted = sort(I_p_wave_areas_2_pulmonary_veins);
I_p_wave_areas_2_pulmonary_veins_dev_avg_min_1SD = area_2(length(I_p_wave_areas_2_pulmonary_veins_sorted),1);
for i = 1:length(I_p_wave_areas_2_pulmonary_veins_sorted)
    I_p_wave_areas_2_pulmonary_veins_dev_avg_min_1SD(i) = I_p_wave_areas_2_pulmonary_veins_sorted(i) - (I_p_wave_areas_2_pulmonary_veins_average-I_p_wave_areas_2_pulmonary_veins_SD);
end
I_p_wave_areas_2_pulmonary_veins_dev_avg_min_1SD = abs(I_p_wave_areas_2_pulmonary_veins_dev_avg_min_1SD);
for i = 1:length(I_p_wave_areas_2_pulmonary_veins_sorted)
    if I_p_wave_areas_2_pulmonary_veins_dev_avg_min_1SD(i) == min(I_p_wave_areas_2_pulmonary_veins_dev_avg_min_1SD)
        I_p_wave_areas_2_pulmonary_veins_dev_avg_min_1SD_index = i;
    end
end
I_p_wave_areas_2_pulmonary_veins_dev_avg_pls_1SD = area_2(length(I_p_wave_areas_2_pulmonary_veins_sorted),1);
for i = 1:length(I_p_wave_areas_2_pulmonary_veins_sorted)
    I_p_wave_areas_2_pulmonary_veins_dev_avg_pls_1SD(i) = I_p_wave_areas_2_pulmonary_veins_sorted(i) - (I_p_wave_areas_2_pulmonary_veins_average+I_p_wave_areas_2_pulmonary_veins_SD);
end
I_p_wave_areas_2_pulmonary_veins_dev_avg_pls_1SD = abs(I_p_wave_areas_2_pulmonary_veins_dev_avg_pls_1SD);
for i = 1:length(I_p_wave_areas_2_pulmonary_veins_sorted)
    if I_p_wave_areas_2_pulmonary_veins_dev_avg_pls_1SD(i) == min(I_p_wave_areas_2_pulmonary_veins_dev_avg_pls_1SD)
        I_p_wave_areas_2_pulmonary_veins_dev_avg_pls_1SD_index = i;
    end
end
I_p_wave_areas_2_pulmonary_veins_dev_avg_min_2SD = area_2(length(I_p_wave_areas_2_pulmonary_veins_sorted),1);
for i = 1:length(I_p_wave_areas_2_pulmonary_veins_sorted)
    I_p_wave_areas_2_pulmonary_veins_dev_avg_min_2SD(i) = I_p_wave_areas_2_pulmonary_veins_sorted(i) - (I_p_wave_areas_2_pulmonary_veins_average-2*I_p_wave_areas_2_pulmonary_veins_SD);
end
I_p_wave_areas_2_pulmonary_veins_dev_avg_min_2SD = abs(I_p_wave_areas_2_pulmonary_veins_dev_avg_min_2SD);
for i = 1:length(I_p_wave_areas_2_pulmonary_veins_sorted)
    if I_p_wave_areas_2_pulmonary_veins_dev_avg_min_2SD(i) == min(I_p_wave_areas_2_pulmonary_veins_dev_avg_min_2SD)
        I_p_wave_areas_2_dev_avg_min_2SD_index = i;
    end
end
I_p_wave_areas_2_pulmonary_veins_dev_avg_pls_2SD = area_2(length(I_p_wave_areas_2_pulmonary_veins_sorted),1);
for i = 1:length(I_p_wave_areas_2_pulmonary_veins_sorted)
    I_p_wave_areas_2_pulmonary_veins_dev_avg_pls_2SD(i) = I_p_wave_areas_2_pulmonary_veins_sorted(i) - (I_p_wave_areas_2_pulmonary_veins_average+2*I_p_wave_areas_2_pulmonary_veins_SD);
end
I_p_wave_areas_2_pulmonary_veins_dev_avg_pls_2SD = abs(I_p_wave_areas_2_pulmonary_veins_dev_avg_pls_2SD);
for i = 1:length(I_p_wave_areas_2_pulmonary_veins_sorted)
    if I_p_wave_areas_2_pulmonary_veins_dev_avg_pls_2SD(i) == min(I_p_wave_areas_2_pulmonary_veins_dev_avg_pls_2SD)
        I_p_wave_areas_2_pulmonary_veins_dev_avg_pls_2SD_index = i;
    end
end
I_p_wave_areas_2_pulmonary_veins_dev_avg_min_3SD = area_2(length(I_p_wave_areas_2_pulmonary_veins_sorted),1);
for i = 1:length(I_p_wave_areas_2_pulmonary_veins_sorted)
    I_p_wave_areas_2_pulmonary_veins_dev_avg_min_3SD(i) = I_p_wave_areas_2_pulmonary_veins_sorted(i) - (I_p_wave_areas_2_pulmonary_veins_average-3*I_p_wave_areas_2_pulmonary_veins_SD);
end
I_p_wave_areas_2_pulmonary_veins_dev_avg_min_3SD = abs(I_p_wave_areas_2_pulmonary_veins_dev_avg_min_3SD);
for i = 1:length(I_p_wave_areas_2_pulmonary_veins_sorted)
    if I_p_wave_areas_2_pulmonary_veins_dev_avg_min_3SD(i) == min(I_p_wave_areas_2_pulmonary_veins_dev_avg_min_3SD)
        I_p_wave_areas_2_pulmonary_veins_dev_avg_min_3SD_index = i;
    end
end
I_p_wave_areas_2_pulmonary_veins_dev_avg_pls_3SD = area_2(length(I_p_wave_areas_2_pulmonary_veins_sorted),1);
for i = 1:length(I_p_wave_areas_2_pulmonary_veins_sorted)
    I_p_wave_areas_2_pulmonary_veins_dev_avg_pls_3SD(i) = I_p_wave_areas_2_pulmonary_veins_sorted(i) - (I_p_wave_areas_2_pulmonary_veins_average+3*I_p_wave_areas_2_pulmonary_veins_SD);
end
I_p_wave_areas_2_pulmonary_veins_dev_avg_pls_3SD = abs(I_p_wave_areas_2_pulmonary_veins_dev_avg_pls_3SD);
for i = 1:length(I_p_wave_areas_2_pulmonary_veins_sorted)
    if I_p_wave_areas_2_pulmonary_veins_dev_avg_pls_3SD(i) == min(I_p_wave_areas_2_pulmonary_veins_dev_avg_pls_3SD)
        I_p_wave_areas_2_pulmonary_veins_dev_avg_pls_3SD_index = i;
    end
end
I_p_wave_areas_2_pulmonary_veins_tot_number_of_indices = length(I_p_wave_areas_2_pulmonary_veins_sorted);
%This is the approximation of the tot "area under curve." 
I_p_wave_areas_2_pulmonary_veins_quotient_within_1SD = (I_p_wave_areas_2_pulmonary_veins_dev_avg_pls_1SD_index - I_p_wave_areas_2_pulmonary_veinsdev_avg_min_1SD_index)/I_p_wave_areas_2_pulmonary_veins_tot_number_of_indices;
I_p_wave_areas_2_pulmonary_veins_quotient_within_2SD = (I_p_wave_areas_2_pulmonary_veins_dev_avg_pls_2SD_index - I_p_wave_areas_2_dev_avg_min_2SD_index)/I_p_wave_areas_2_pulmonary_veins_tot_number_of_indices;
I_p_wave_areas_2_pulmonary_veins_quotient_within_3SD = (I_p_wave_areas_2_pulmonary_veins_dev_avg_pls_3SD_index - I_p_wave_areas_2_pulmonary_veins_dev_avg_min_3SD_index)/I_p_wave_areas_2_pulmonary_veins_tot_number_of_indices;
if (I_p_wave_areas_2_pulmonary_veins_quotient_within_1SD > 0.66 && I_p_wave_areas_2_pulmonary_veins_quotient_within_1SD < 0.70) && (I_p_wave_areas_2_pulmonary_veins_quotient_within_2SD > 0.93 && I_p_wave_areas_2_pulmonary_veins_quotient_within_2SD < 0.97) && (I_p_wave_areas_2_pulmonary_veins_quotient_within_3SD > 0.98 && I_p_wave_areas_2_pulmonary_veins_quotient_within_3SD < 1)
    isNormal(8,9) = 1;
end
if I_p_wave_areas_2_pulmonary_veins_quotient_within_1SD == 0
    isSpike(8,9) = 1;
end
I_p_wave_areas_2_mitral_annulus_sorted = sort(I_p_wave_areas_2_mitral_annulus);
I_p_wave_areas_2_mitral_annulus_dev_avg_min_1SD = area_2(length(I_p_wave_areas_2_mitral_annulus_sorted),1);
for i = 1:length(I_p_wave_areas_2_mitral_annulus_sorted)
    I_p_wave_areas_2_mitral_annulus_dev_avg_min_1SD(i) = I_p_wave_areas_2_mitral_annulus_sorted(i) - (I_p_wave_areas_2_mitral_annulus_average-I_p_wave_areas_2_mitral_annulus_SD);
end
I_p_wave_areas_2_mitral_annulus_dev_avg_min_1SD = abs(I_p_wave_areas_2_mitral_annulus_dev_avg_min_1SD);
for i = 1:length(I_p_wave_areas_2_mitral_annulus_sorted)
    if I_p_wave_areas_2_mitral_annulus_dev_avg_min_1SD(i) == min(I_p_wave_areas_2_mitral_annulus_dev_avg_min_1SD)
        I_p_wave_areas_2_mitral_annulus_dev_avg_min_1SD_index = i;
    end
end
I_p_wave_areas_2_mitral_annulus_dev_avg_pls_1SD = area_2(length(I_p_wave_areas_2_mitral_annulus_sorted),1);
for i = 1:length(I_p_wave_areas_2_mitral_annulus_sorted)
    I_p_wave_areas_2_mitral_annulus_dev_avg_pls_1SD(i) = I_p_wave_areas_2_mitral_annulus_sorted(i) - (I_p_wave_areas_2_mitral_annulus_average+I_p_wave_areas_2_mitral_annulus_SD);
end
I_p_wave_areas_2_mitral_annulus_dev_avg_pls_1SD = abs(I_p_wave_areas_2_mitral_annulus_dev_avg_pls_1SD);
for i = 1:length(I_p_wave_areas_2_mitral_annulus_sorted)
    if I_p_wave_areas_2_mitral_annulus_dev_avg_pls_1SD(i) == min(I_p_wave_areas_2_mitral_annulus_dev_avg_pls_1SD)
        I_p_wave_areas_2_mitral_annulus_dev_avg_pls_1SD_index = i;
    end
end
I_p_wave_areas_2_mitral_annulus_dev_avg_min_2SD = area_2(length(I_p_wave_areas_2_mitral_annulus_sorted),1);
for i = 1:length(I_p_wave_areas_2_mitral_annulus_sorted)
    I_p_wave_areas_2_mitral_annulus_dev_avg_min_2SD(i) = I_p_wave_areas_2_mitral_annulus_sorted(i) - (I_p_wave_areas_2_mitral_annulus_average-2*I_p_wave_areas_2_mitral_annulus_SD);
end
I_p_wave_areas_2_mitral_annulus_dev_avg_min_2SD = abs(I_p_wave_areas_2_mitral_annulus_dev_avg_min_2SD);
for i = 1:length(I_p_wave_areas_2_mitral_annulus_sorted)
    if I_p_wave_areas_2_mitral_annulus_dev_avg_min_2SD(i) == min(I_p_wave_areas_2_mitral_annulus_dev_avg_min_2SD)
        I_p_wave_areas_2_dev_avg_min_2SD_index = i;
    end
end
I_p_wave_areas_2_mitral_annulus_dev_avg_pls_2SD = area_2(length(I_p_wave_areas_2_mitral_annulus_sorted),1);
for i = 1:length(I_p_wave_areas_2_mitral_annulus_sorted)
    I_p_wave_areas_2_mitral_annulus_dev_avg_pls_2SD(i) = I_p_wave_areas_2_mitral_annulus_sorted(i) - (I_p_wave_areas_2_mitral_annulus_average+2*I_p_wave_areas_2_mitral_annulus_SD);
end
I_p_wave_areas_2_mitral_annulus_dev_avg_pls_2SD = abs(I_p_wave_areas_2_mitral_annulus_dev_avg_pls_2SD);
for i = 1:length(I_p_wave_areas_2_mitral_annulus_sorted)
    if I_p_wave_areas_2_mitral_annulus_dev_avg_pls_2SD(i) == min(I_p_wave_areas_2_mitral_annulus_dev_avg_pls_2SD)
        I_p_wave_areas_2_mitral_annulus_dev_avg_pls_2SD_index = i;
    end
end
I_p_wave_areas_2_mitral_annulus_dev_avg_min_3SD = area_2(length(I_p_wave_areas_2_mitral_annulus_sorted),1);
for i = 1:length(I_p_wave_areas_2_mitral_annulus_sorted)
    I_p_wave_areas_2_mitral_annulus_dev_avg_min_3SD(i) = I_p_wave_areas_2_mitral_annulus_sorted(i) - (I_p_wave_areas_2_mitral_annulus_average-3*I_p_wave_areas_2_mitral_annulus_SD);
end
I_p_wave_areas_2_mitral_annulus_dev_avg_min_3SD = abs(I_p_wave_areas_2_mitral_annulus_dev_avg_min_3SD);
for i = 1:length(I_p_wave_areas_2_mitral_annulus_sorted)
    if I_p_wave_areas_2_mitral_annulus_dev_avg_min_3SD(i) == min(I_p_wave_areas_2_mitral_annulus_dev_avg_min_3SD)
        I_p_wave_areas_2_mitral_annulus_dev_avg_min_3SD_index = i;
    end
end
I_p_wave_areas_2_mitral_annulus_dev_avg_pls_3SD = area_2(length(I_p_wave_areas_2_mitral_annulus_sorted),1);
for i = 1:length(I_p_wave_areas_2_mitral_annulus_sorted)
    I_p_wave_areas_2_mitral_annulus_dev_avg_pls_3SD(i) = I_p_wave_areas_2_mitral_annulus_sorted(i) - (I_p_wave_areas_2_mitral_annulus_average+3*I_p_wave_areas_2_mitral_annulus_SD);
end
I_p_wave_areas_2_mitral_annulus_dev_avg_pls_3SD = abs(I_p_wave_areas_2_mitral_annulus_dev_avg_pls_3SD);
for i = 1:length(I_p_wave_areas_2_mitral_annulus_sorted)
    if I_p_wave_areas_2_mitral_annulus_dev_avg_pls_3SD(i) == min(I_p_wave_areas_2_mitral_annulus_dev_avg_pls_3SD)
        I_p_wave_areas_2_mitral_annulus_dev_avg_pls_3SD_index = i;
    end
end
I_p_wave_areas_2_mitral_annulus_tot_number_of_indices = length(I_p_wave_areas_2_mitral_annulus_sorted);
%This is the approximation of the tot "area under curve." 
I_p_wave_areas_2_mitral_annulus_quotient_within_1SD = (I_p_wave_areas_2_mitral_annulus_dev_avg_pls_1SD_index - I_p_wave_areas_2_mitral_annulusdev_avg_min_1SD_index)/I_p_wave_areas_2_mitral_annulus_tot_number_of_indices;
I_p_wave_areas_2_mitral_annulus_quotient_within_2SD = (I_p_wave_areas_2_mitral_annulus_dev_avg_pls_2SD_index - I_p_wave_areas_2_dev_avg_min_2SD_index)/I_p_wave_areas_2_mitral_annulus_tot_number_of_indices;
I_p_wave_areas_2_mitral_annulus_quotient_within_3SD = (I_p_wave_areas_2_mitral_annulus_dev_avg_pls_3SD_index - I_p_wave_areas_2_mitral_annulus_dev_avg_min_3SD_index)/I_p_wave_areas_2_mitral_annulus_tot_number_of_indices;
if (I_p_wave_areas_2_mitral_annulus_quotient_within_1SD > 0.66 && I_p_wave_areas_2_mitral_annulus_quotient_within_1SD < 0.70) && (I_p_wave_areas_2_mitral_annulus_quotient_within_2SD > 0.93 && I_p_wave_areas_2_mitral_annulus_quotient_within_2SD < 0.97) && (I_p_wave_areas_2_mitral_annulus_quotient_within_3SD > 0.98 && I_p_wave_areas_2_mitral_annulus_quotient_within_3SD < 1)
    isNormal(8,10) = 1;
end
if I_p_wave_areas_2_mitral_annulus_quotient_within_1SD == 0
    isSpike(8,10) = 1;
end
I_p_wave_areas_2_CS_body_sorted = sort(I_p_wave_areas_2_CS_body);
I_p_wave_areas_2_CS_body_dev_avg_min_1SD = area_2(length(I_p_wave_areas_2_CS_body_sorted),1);
for i = 1:length(I_p_wave_areas_2_CS_body_sorted)
    I_p_wave_areas_2_CS_body_dev_avg_min_1SD(i) = I_p_wave_areas_2_CS_body_sorted(i) - (I_p_wave_areas_2_CS_body_average-I_p_wave_areas_2_CS_body_SD);
end
I_p_wave_areas_2_CS_body_dev_avg_min_1SD = abs(I_p_wave_areas_2_CS_body_dev_avg_min_1SD);
for i = 1:length(I_p_wave_areas_2_CS_body_sorted)
    if I_p_wave_areas_2_CS_body_dev_avg_min_1SD(i) == min(I_p_wave_areas_2_CS_body_dev_avg_min_1SD)
        I_p_wave_areas_2_CS_body_dev_avg_min_1SD_index = i;
    end
end
I_p_wave_areas_2_CS_body_dev_avg_pls_1SD = area_2(length(I_p_wave_areas_2_CS_body_sorted),1);
for i = 1:length(I_p_wave_areas_2_CS_body_sorted)
    I_p_wave_areas_2_CS_body_dev_avg_pls_1SD(i) = I_p_wave_areas_2_CS_body_sorted(i) - (I_p_wave_areas_2_CS_body_average+I_p_wave_areas_2_CS_body_SD);
end
I_p_wave_areas_2_CS_body_dev_avg_pls_1SD = abs(I_p_wave_areas_2_CS_body_dev_avg_pls_1SD);
for i = 1:length(I_p_wave_areas_2_CS_body_sorted)
    if I_p_wave_areas_2_CS_body_dev_avg_pls_1SD(i) == min(I_p_wave_areas_2_CS_body_dev_avg_pls_1SD)
        I_p_wave_areas_2_CS_body_dev_avg_pls_1SD_index = i;
    end
end
I_p_wave_areas_2_CS_body_dev_avg_min_2SD = area_2(length(I_p_wave_areas_2_CS_body_sorted),1);
for i = 1:length(I_p_wave_areas_2_CS_body_sorted)
    I_p_wave_areas_2_CS_body_dev_avg_min_2SD(i) = I_p_wave_areas_2_CS_body_sorted(i) - (I_p_wave_areas_2_CS_body_average-2*I_p_wave_areas_2_CS_body_SD);
end
I_p_wave_areas_2_CS_body_dev_avg_min_2SD = abs(I_p_wave_areas_2_CS_body_dev_avg_min_2SD);
for i = 1:length(I_p_wave_areas_2_CS_body_sorted)
    if I_p_wave_areas_2_CS_body_dev_avg_min_2SD(i) == min(I_p_wave_areas_2_CS_body_dev_avg_min_2SD)
        I_p_wave_areas_2_dev_avg_min_2SD_index = i;
    end
end
I_p_wave_areas_2_CS_body_dev_avg_pls_2SD = area_2(length(I_p_wave_areas_2_CS_body_sorted),1);
for i = 1:length(I_p_wave_areas_2_CS_body_sorted)
    I_p_wave_areas_2_CS_body_dev_avg_pls_2SD(i) = I_p_wave_areas_2_CS_body_sorted(i) - (I_p_wave_areas_2_CS_body_average+2*I_p_wave_areas_2_CS_body_SD);
end
I_p_wave_areas_2_CS_body_dev_avg_pls_2SD = abs(I_p_wave_areas_2_CS_body_dev_avg_pls_2SD);
for i = 1:length(I_p_wave_areas_2_CS_body_sorted)
    if I_p_wave_areas_2_CS_body_dev_avg_pls_2SD(i) == min(I_p_wave_areas_2_CS_body_dev_avg_pls_2SD)
        I_p_wave_areas_2_CS_body_dev_avg_pls_2SD_index = i;
    end
end
I_p_wave_areas_2_CS_body_dev_avg_min_3SD = area_2(length(I_p_wave_areas_2_CS_body_sorted),1);
for i = 1:length(I_p_wave_areas_2_CS_body_sorted)
    I_p_wave_areas_2_CS_body_dev_avg_min_3SD(i) = I_p_wave_areas_2_CS_body_sorted(i) - (I_p_wave_areas_2_CS_body_average-3*I_p_wave_areas_2_CS_body_SD);
end
I_p_wave_areas_2_CS_body_dev_avg_min_3SD = abs(I_p_wave_areas_2_CS_body_dev_avg_min_3SD);
for i = 1:length(I_p_wave_areas_2_CS_body_sorted)
    if I_p_wave_areas_2_CS_body_dev_avg_min_3SD(i) == min(I_p_wave_areas_2_CS_body_dev_avg_min_3SD)
        I_p_wave_areas_2_CS_body_dev_avg_min_3SD_index = i;
    end
end
I_p_wave_areas_2_CS_body_dev_avg_pls_3SD = area_2(length(I_p_wave_areas_2_CS_body_sorted),1);
for i = 1:length(I_p_wave_areas_2_CS_body_sorted)
    I_p_wave_areas_2_CS_body_dev_avg_pls_3SD(i) = I_p_wave_areas_2_CS_body_sorted(i) - (I_p_wave_areas_2_CS_body_average+3*I_p_wave_areas_2_CS_body_SD);
end
I_p_wave_areas_2_CS_body_dev_avg_pls_3SD = abs(I_p_wave_areas_2_CS_body_dev_avg_pls_3SD);
for i = 1:length(I_p_wave_areas_2_CS_body_sorted)
    if I_p_wave_areas_2_CS_body_dev_avg_pls_3SD(i) == min(I_p_wave_areas_2_CS_body_dev_avg_pls_3SD)
        I_p_wave_areas_2_CS_body_dev_avg_pls_3SD_index = i;
    end
end
I_p_wave_areas_2_CS_body_tot_number_of_indices = length(I_p_wave_areas_2_CS_body_sorted);
%This is the approximation of the tot "area under curve." 
I_p_wave_areas_2_CS_body_quotient_within_1SD = (I_p_wave_areas_2_CS_body_dev_avg_pls_1SD_index - I_p_wave_areas_2_CS_bodydev_avg_min_1SD_index)/I_p_wave_areas_2_CS_body_tot_number_of_indices;
I_p_wave_areas_2_CS_body_quotient_within_2SD = (I_p_wave_areas_2_CS_body_dev_avg_pls_2SD_index - I_p_wave_areas_2_dev_avg_min_2SD_index)/I_p_wave_areas_2_CS_body_tot_number_of_indices;
I_p_wave_areas_2_CS_body_quotient_within_3SD = (I_p_wave_areas_2_CS_body_dev_avg_pls_3SD_index - I_p_wave_areas_2_CS_body_dev_avg_min_3SD_index)/I_p_wave_areas_2_CS_body_tot_number_of_indices;
if (I_p_wave_areas_2_CS_body_quotient_within_1SD > 0.66 && I_p_wave_areas_2_CS_body_quotient_within_1SD < 0.70) && (I_p_wave_areas_2_CS_body_quotient_within_2SD > 0.93 && I_p_wave_areas_2_CS_body_quotient_within_2SD < 0.97) && (I_p_wave_areas_2_CS_body_quotient_within_3SD > 0.98 && I_p_wave_areas_2_CS_body_quotient_within_3SD < 1)
    isNormal(8,11) = 1;
end
if I_p_wave_areas_2_CS_body_quotient_within_1SD == 0
    isSpike(8,11) = 1;
end
I_p_wave_areas_2_left_septum_sorted = sort(I_p_wave_areas_2_left_septum);
I_p_wave_areas_2_left_septum_dev_avg_min_1SD = area_2(length(I_p_wave_areas_2_left_septum_sorted),1);
for i = 1:length(I_p_wave_areas_2_left_septum_sorted)
    I_p_wave_areas_2_left_septum_dev_avg_min_1SD(i) = I_p_wave_areas_2_left_septum_sorted(i) - (I_p_wave_areas_2_left_septum_average-I_p_wave_areas_2_left_septum_SD);
end
I_p_wave_areas_2_left_septum_dev_avg_min_1SD = abs(I_p_wave_areas_2_left_septum_dev_avg_min_1SD);
for i = 1:length(I_p_wave_areas_2_left_septum_sorted)
    if I_p_wave_areas_2_left_septum_dev_avg_min_1SD(i) == min(I_p_wave_areas_2_left_septum_dev_avg_min_1SD)
        I_p_wave_areas_2_left_septum_dev_avg_min_1SD_index = i;
    end
end
I_p_wave_areas_2_left_septum_dev_avg_pls_1SD = area_2(length(I_p_wave_areas_2_left_septum_sorted),1);
for i = 1:length(I_p_wave_areas_2_left_septum_sorted)
    I_p_wave_areas_2_left_septum_dev_avg_pls_1SD(i) = I_p_wave_areas_2_left_septum_sorted(i) - (I_p_wave_areas_2_left_septum_average+I_p_wave_areas_2_left_septum_SD);
end
I_p_wave_areas_2_left_septum_dev_avg_pls_1SD = abs(I_p_wave_areas_2_left_septum_dev_avg_pls_1SD);
for i = 1:length(I_p_wave_areas_2_left_septum_sorted)
    if I_p_wave_areas_2_left_septum_dev_avg_pls_1SD(i) == min(I_p_wave_areas_2_left_septum_dev_avg_pls_1SD)
        I_p_wave_areas_2_left_septum_dev_avg_pls_1SD_index = i;
    end
end
I_p_wave_areas_2_left_septum_dev_avg_min_2SD = area_2(length(I_p_wave_areas_2_left_septum_sorted),1);
for i = 1:length(I_p_wave_areas_2_left_septum_sorted)
    I_p_wave_areas_2_left_septum_dev_avg_min_2SD(i) = I_p_wave_areas_2_left_septum_sorted(i) - (I_p_wave_areas_2_left_septum_average-2*I_p_wave_areas_2_left_septum_SD);
end
I_p_wave_areas_2_left_septum_dev_avg_min_2SD = abs(I_p_wave_areas_2_left_septum_dev_avg_min_2SD);
for i = 1:length(I_p_wave_areas_2_left_septum_sorted)
    if I_p_wave_areas_2_left_septum_dev_avg_min_2SD(i) == min(I_p_wave_areas_2_left_septum_dev_avg_min_2SD)
        I_p_wave_areas_2_dev_avg_min_2SD_index = i;
    end
end
I_p_wave_areas_2_left_septum_dev_avg_pls_2SD = area_2(length(I_p_wave_areas_2_left_septum_sorted),1);
for i = 1:length(I_p_wave_areas_2_left_septum_sorted)
    I_p_wave_areas_2_left_septum_dev_avg_pls_2SD(i) = I_p_wave_areas_2_left_septum_sorted(i) - (I_p_wave_areas_2_left_septum_average+2*I_p_wave_areas_2_left_septum_SD);
end
I_p_wave_areas_2_left_septum_dev_avg_pls_2SD = abs(I_p_wave_areas_2_left_septum_dev_avg_pls_2SD);
for i = 1:length(I_p_wave_areas_2_left_septum_sorted)
    if I_p_wave_areas_2_left_septum_dev_avg_pls_2SD(i) == min(I_p_wave_areas_2_left_septum_dev_avg_pls_2SD)
        I_p_wave_areas_2_left_septum_dev_avg_pls_2SD_index = i;
    end
end
I_p_wave_areas_2_left_septum_dev_avg_min_3SD = area_2(length(I_p_wave_areas_2_left_septum_sorted),1);
for i = 1:length(I_p_wave_areas_2_left_septum_sorted)
    I_p_wave_areas_2_left_septum_dev_avg_min_3SD(i) = I_p_wave_areas_2_left_septum_sorted(i) - (I_p_wave_areas_2_left_septum_average-3*I_p_wave_areas_2_left_septum_SD);
end
I_p_wave_areas_2_left_septum_dev_avg_min_3SD = abs(I_p_wave_areas_2_left_septum_dev_avg_min_3SD);
for i = 1:length(I_p_wave_areas_2_left_septum_sorted)
    if I_p_wave_areas_2_left_septum_dev_avg_min_3SD(i) == min(I_p_wave_areas_2_left_septum_dev_avg_min_3SD)
        I_p_wave_areas_2_left_septum_dev_avg_min_3SD_index = i;
    end
end
I_p_wave_areas_2_left_septum_dev_avg_pls_3SD = area_2(length(I_p_wave_areas_2_left_septum_sorted),1);
for i = 1:length(I_p_wave_areas_2_left_septum_sorted)
    I_p_wave_areas_2_left_septum_dev_avg_pls_3SD(i) = I_p_wave_areas_2_left_septum_sorted(i) - (I_p_wave_areas_2_left_septum_average+3*I_p_wave_areas_2_left_septum_SD);
end
I_p_wave_areas_2_left_septum_dev_avg_pls_3SD = abs(I_p_wave_areas_2_left_septum_dev_avg_pls_3SD);
for i = 1:length(I_p_wave_areas_2_left_septum_sorted)
    if I_p_wave_areas_2_left_septum_dev_avg_pls_3SD(i) == min(I_p_wave_areas_2_left_septum_dev_avg_pls_3SD)
        I_p_wave_areas_2_left_septum_dev_avg_pls_3SD_index = i;
    end
end
I_p_wave_areas_2_left_septum_tot_number_of_indices = length(I_p_wave_areas_2_left_septum_sorted);
%This is the approximation of the tot "area under curve." 
I_p_wave_areas_2_left_septum_quotient_within_1SD = (I_p_wave_areas_2_left_septum_dev_avg_pls_1SD_index - I_p_wave_areas_2_left_septumdev_avg_min_1SD_index)/I_p_wave_areas_2_left_septum_tot_number_of_indices;
I_p_wave_areas_2_left_septum_quotient_within_2SD = (I_p_wave_areas_2_left_septum_dev_avg_pls_2SD_index - I_p_wave_areas_2_dev_avg_min_2SD_index)/I_p_wave_areas_2_left_septum_tot_number_of_indices;
I_p_wave_areas_2_left_septum_quotient_within_3SD = (I_p_wave_areas_2_left_septum_dev_avg_pls_3SD_index - I_p_wave_areas_2_left_septum_dev_avg_min_3SD_index)/I_p_wave_areas_2_left_septum_tot_number_of_indices;
if (I_p_wave_areas_2_left_septum_quotient_within_1SD > 0.66 && I_p_wave_areas_2_left_septum_quotient_within_1SD < 0.70) && (I_p_wave_areas_2_left_septum_quotient_within_2SD > 0.93 && I_p_wave_areas_2_left_septum_quotient_within_2SD < 0.97) && (I_p_wave_areas_2_left_septum_quotient_within_3SD > 0.98 && I_p_wave_areas_2_left_septum_quotient_within_3SD < 1)
    isNormal(8,12) = 1;
end
if I_p_wave_areas_2_left_septum_quotient_within_1SD == 0
    isSpike(8,12) = 1;
end
I_p_wave_areas_2_left_atrial_appendage_sorted = sort(I_p_wave_areas_2_left_atrial_appendage);
I_p_wave_areas_2_left_atrial_appendage_dev_avg_min_1SD = area_2(length(I_p_wave_areas_2_left_atrial_appendage_sorted),1);
for i = 1:length(I_p_wave_areas_2_left_atrial_appendage_sorted)
    I_p_wave_areas_2_left_atrial_appendage_dev_avg_min_1SD(i) = I_p_wave_areas_2_left_atrial_appendage_sorted(i) - (I_p_wave_areas_2_left_atrial_appendage_average-I_p_wave_areas_2_left_atrial_appendage_SD);
end
I_p_wave_areas_2_left_atrial_appendage_dev_avg_min_1SD = abs(I_p_wave_areas_2_left_atrial_appendage_dev_avg_min_1SD);
for i = 1:length(I_p_wave_areas_2_left_atrial_appendage_sorted)
    if I_p_wave_areas_2_left_atrial_appendage_dev_avg_min_1SD(i) == min(I_p_wave_areas_2_left_atrial_appendage_dev_avg_min_1SD)
        I_p_wave_areas_2_left_atrial_appendage_dev_avg_min_1SD_index = i;
    end
end
I_p_wave_areas_2_left_atrial_appendage_dev_avg_pls_1SD = area_2(length(I_p_wave_areas_2_left_atrial_appendage_sorted),1);
for i = 1:length(I_p_wave_areas_2_left_atrial_appendage_sorted)
    I_p_wave_areas_2_left_atrial_appendage_dev_avg_pls_1SD(i) = I_p_wave_areas_2_left_atrial_appendage_sorted(i) - (I_p_wave_areas_2_left_atrial_appendage_average+I_p_wave_areas_2_left_atrial_appendage_SD);
end
I_p_wave_areas_2_left_atrial_appendage_dev_avg_pls_1SD = abs(I_p_wave_areas_2_left_atrial_appendage_dev_avg_pls_1SD);
for i = 1:length(I_p_wave_areas_2_left_atrial_appendage_sorted)
    if I_p_wave_areas_2_left_atrial_appendage_dev_avg_pls_1SD(i) == min(I_p_wave_areas_2_left_atrial_appendage_dev_avg_pls_1SD)
        I_p_wave_areas_2_left_atrial_appendage_dev_avg_pls_1SD_index = i;
    end
end
I_p_wave_areas_2_left_atrial_appendage_dev_avg_min_2SD = area_2(length(I_p_wave_areas_2_left_atrial_appendage_sorted),1);
for i = 1:length(I_p_wave_areas_2_left_atrial_appendage_sorted)
    I_p_wave_areas_2_left_atrial_appendage_dev_avg_min_2SD(i) = I_p_wave_areas_2_left_atrial_appendage_sorted(i) - (I_p_wave_areas_2_left_atrial_appendage_average-2*I_p_wave_areas_2_left_atrial_appendage_SD);
end
I_p_wave_areas_2_left_atrial_appendage_dev_avg_min_2SD = abs(I_p_wave_areas_2_left_atrial_appendage_dev_avg_min_2SD);
for i = 1:length(I_p_wave_areas_2_left_atrial_appendage_sorted)
    if I_p_wave_areas_2_left_atrial_appendage_dev_avg_min_2SD(i) == min(I_p_wave_areas_2_left_atrial_appendage_dev_avg_min_2SD)
        I_p_wave_areas_2_dev_avg_min_2SD_index = i;
    end
end
I_p_wave_areas_2_left_atrial_appendage_dev_avg_pls_2SD = area_2(length(I_p_wave_areas_2_left_atrial_appendage_sorted),1);
for i = 1:length(I_p_wave_areas_2_left_atrial_appendage_sorted)
    I_p_wave_areas_2_left_atrial_appendage_dev_avg_pls_2SD(i) = I_p_wave_areas_2_left_atrial_appendage_sorted(i) - (I_p_wave_areas_2_left_atrial_appendage_average+2*I_p_wave_areas_2_left_atrial_appendage_SD);
end
I_p_wave_areas_2_left_atrial_appendage_dev_avg_pls_2SD = abs(I_p_wave_areas_2_left_atrial_appendage_dev_avg_pls_2SD);
for i = 1:length(I_p_wave_areas_2_left_atrial_appendage_sorted)
    if I_p_wave_areas_2_left_atrial_appendage_dev_avg_pls_2SD(i) == min(I_p_wave_areas_2_left_atrial_appendage_dev_avg_pls_2SD)
        I_p_wave_areas_2_left_atrial_appendage_dev_avg_pls_2SD_index = i;
    end
end
I_p_wave_areas_2_left_atrial_appendage_dev_avg_min_3SD = area_2(length(I_p_wave_areas_2_left_atrial_appendage_sorted),1);
for i = 1:length(I_p_wave_areas_2_left_atrial_appendage_sorted)
    I_p_wave_areas_2_left_atrial_appendage_dev_avg_min_3SD(i) = I_p_wave_areas_2_left_atrial_appendage_sorted(i) - (I_p_wave_areas_2_left_atrial_appendage_average-3*I_p_wave_areas_2_left_atrial_appendage_SD);
end
I_p_wave_areas_2_left_atrial_appendage_dev_avg_min_3SD = abs(I_p_wave_areas_2_left_atrial_appendage_dev_avg_min_3SD);
for i = 1:length(I_p_wave_areas_2_left_atrial_appendage_sorted)
    if I_p_wave_areas_2_left_atrial_appendage_dev_avg_min_3SD(i) == min(I_p_wave_areas_2_left_atrial_appendage_dev_avg_min_3SD)
        I_p_wave_areas_2_left_atrial_appendage_dev_avg_min_3SD_index = i;
    end
end
I_p_wave_areas_2_left_atrial_appendage_dev_avg_pls_3SD = area_2(length(I_p_wave_areas_2_left_atrial_appendage_sorted),1);
for i = 1:length(I_p_wave_areas_2_left_atrial_appendage_sorted)
    I_p_wave_areas_2_left_atrial_appendage_dev_avg_pls_3SD(i) = I_p_wave_areas_2_left_atrial_appendage_sorted(i) - (I_p_wave_areas_2_left_atrial_appendage_average+3*I_p_wave_areas_2_left_atrial_appendage_SD);
end
I_p_wave_areas_2_left_atrial_appendage_dev_avg_pls_3SD = abs(I_p_wave_areas_2_left_atrial_appendage_dev_avg_pls_3SD);
for i = 1:length(I_p_wave_areas_2_left_atrial_appendage_sorted)
    if I_p_wave_areas_2_left_atrial_appendage_dev_avg_pls_3SD(i) == min(I_p_wave_areas_2_left_atrial_appendage_dev_avg_pls_3SD)
        I_p_wave_areas_2_left_atrial_appendage_dev_avg_pls_3SD_index = i;
    end
end
I_p_wave_areas_2_left_atrial_appendage_tot_number_of_indices = length(I_p_wave_areas_2_left_atrial_appendage_sorted);
%This is the approximation of the tot "area under curve." 
I_p_wave_areas_2_left_atrial_appendage_quotient_within_1SD = (I_p_wave_areas_2_left_atrial_appendage_dev_avg_pls_1SD_index - I_p_wave_areas_2_left_atrial_appendagedev_avg_min_1SD_index)/I_p_wave_areas_2_left_atrial_appendage_tot_number_of_indices;
I_p_wave_areas_2_left_atrial_appendage_quotient_within_2SD = (I_p_wave_areas_2_left_atrial_appendage_dev_avg_pls_2SD_index - I_p_wave_areas_2_dev_avg_min_2SD_index)/I_p_wave_areas_2_left_atrial_appendage_tot_number_of_indices;
I_p_wave_areas_2_left_atrial_appendage_quotient_within_3SD = (I_p_wave_areas_2_left_atrial_appendage_dev_avg_pls_3SD_index - I_p_wave_areas_2_left_atrial_appendage_dev_avg_min_3SD_index)/I_p_wave_areas_2_left_atrial_appendage_tot_number_of_indices;
if (I_p_wave_areas_2_left_atrial_appendage_quotient_within_1SD > 0.66 && I_p_wave_areas_2_left_atrial_appendage_quotient_within_1SD < 0.70) && (I_p_wave_areas_2_left_atrial_appendage_quotient_within_2SD > 0.93 && I_p_wave_areas_2_left_atrial_appendage_quotient_within_2SD < 0.97) && (I_p_wave_areas_2_left_atrial_appendage_quotient_within_3SD > 0.98 && I_p_wave_areas_2_left_atrial_appendage_quotient_within_3SD < 1)
    isNormal(8,13) = 1;
end
if I_p_wave_areas_2_left_atrial_appendage_quotient_within_1SD == 0
    isSpike(8,13) = 1;
end
I_p_wave_areas_3_SA_node_sorted = sort(I_p_wave_areas_3_SA_node);
I_p_wave_areas_3_SA_node_dev_avg_min_1SD = area_3(length(I_p_wave_areas_3_SA_node_sorted),1);
for i = 1:length(I_p_wave_areas_3_SA_node_sorted)
    I_p_wave_areas_3_SA_node_dev_avg_min_1SD(i) = I_p_wave_areas_3_SA_node_sorted(i) - (I_p_wave_areas_3_SA_node_average-I_p_wave_areas_3_SA_node_SD);
end
I_p_wave_areas_3_SA_node_dev_avg_min_1SD = abs(I_p_wave_areas_3_SA_node_dev_avg_min_1SD);
for i = 1:length(I_p_wave_areas_3_SA_node_sorted)
    if I_p_wave_areas_3_SA_node_dev_avg_min_1SD(i) == min(I_p_wave_areas_3_SA_node_dev_avg_min_1SD)
        I_p_wave_areas_3_SA_node_dev_avg_min_1SD_index = i;
    end
end
I_p_wave_areas_3_SA_node_dev_avg_pls_1SD = area_3(length(I_p_wave_areas_3_SA_node_sorted),1);
for i = 1:length(I_p_wave_areas_3_SA_node_sorted)
    I_p_wave_areas_3_SA_node_dev_avg_pls_1SD(i) = I_p_wave_areas_3_SA_node_sorted(i) - (I_p_wave_areas_3_SA_node_average+I_p_wave_areas_3_SA_node_SD);
end
I_p_wave_areas_3_SA_node_dev_avg_pls_1SD = abs(I_p_wave_areas_3_SA_node_dev_avg_pls_1SD);
for i = 1:length(I_p_wave_areas_3_SA_node_sorted)
    if I_p_wave_areas_3_SA_node_dev_avg_pls_1SD(i) == min(I_p_wave_areas_3_SA_node_dev_avg_pls_1SD)
        I_p_wave_areas_3_SA_node_dev_avg_pls_1SD_index = i;
    end
end
I_p_wave_areas_3_SA_node_dev_avg_min_2SD = area_3(length(I_p_wave_areas_3_SA_node_sorted),1);
for i = 1:length(I_p_wave_areas_3_SA_node_sorted)
    I_p_wave_areas_3_SA_node_dev_avg_min_2SD(i) = I_p_wave_areas_3_SA_node_sorted(i) - (I_p_wave_areas_3_SA_node_average-2*I_p_wave_areas_3_SA_node_SD);
end
I_p_wave_areas_3_SA_node_dev_avg_min_2SD = abs(I_p_wave_areas_3_SA_node_dev_avg_min_2SD);
for i = 1:length(I_p_wave_areas_3_SA_node_sorted)
    if I_p_wave_areas_3_SA_node_dev_avg_min_2SD(i) == min(I_p_wave_areas_3_SA_node_dev_avg_min_2SD)
        I_p_wave_areas_3_dev_avg_min_2SD_index = i;
    end
end
I_p_wave_areas_3_SA_node_dev_avg_pls_2SD = area_3(length(I_p_wave_areas_3_SA_node_sorted),1);
for i = 1:length(I_p_wave_areas_3_SA_node_sorted)
    I_p_wave_areas_3_SA_node_dev_avg_pls_2SD(i) = I_p_wave_areas_3_SA_node_sorted(i) - (I_p_wave_areas_3_SA_node_average+2*I_p_wave_areas_3_SA_node_SD);
end
I_p_wave_areas_3_SA_node_dev_avg_pls_2SD = abs(I_p_wave_areas_3_SA_node_dev_avg_pls_2SD);
for i = 1:length(I_p_wave_areas_3_SA_node_sorted)
    if I_p_wave_areas_3_SA_node_dev_avg_pls_2SD(i) == min(I_p_wave_areas_3_SA_node_dev_avg_pls_2SD)
        I_p_wave_areas_3_SA_node_dev_avg_pls_2SD_index = i;
    end
end
I_p_wave_areas_3_SA_node_dev_avg_min_3SD = area_3(length(I_p_wave_areas_3_SA_node_sorted),1);
for i = 1:length(I_p_wave_areas_3_SA_node_sorted)
    I_p_wave_areas_3_SA_node_dev_avg_min_3SD(i) = I_p_wave_areas_3_SA_node_sorted(i) - (I_p_wave_areas_3_SA_node_average-3*I_p_wave_areas_3_SA_node_SD);
end
I_p_wave_areas_3_SA_node_dev_avg_min_3SD = abs(I_p_wave_areas_3_SA_node_dev_avg_min_3SD);
for i = 1:length(I_p_wave_areas_3_SA_node_sorted)
    if I_p_wave_areas_3_SA_node_dev_avg_min_3SD(i) == min(I_p_wave_areas_3_SA_node_dev_avg_min_3SD)
        I_p_wave_areas_3_SA_node_dev_avg_min_3SD_index = i;
    end
end
I_p_wave_areas_3_SA_node_dev_avg_pls_3SD = area_3(length(I_p_wave_areas_3_SA_node_sorted),1);
for i = 1:length(I_p_wave_areas_3_SA_node_sorted)
    I_p_wave_areas_3_SA_node_dev_avg_pls_3SD(i) = I_p_wave_areas_3_SA_node_sorted(i) - (I_p_wave_areas_3_SA_node_average+3*I_p_wave_areas_3_SA_node_SD);
end
I_p_wave_areas_3_SA_node_dev_avg_pls_3SD = abs(I_p_wave_areas_3_SA_node_dev_avg_pls_3SD);
for i = 1:length(I_p_wave_areas_3_SA_node_sorted)
    if I_p_wave_areas_3_SA_node_dev_avg_pls_3SD(i) == min(I_p_wave_areas_3_SA_node_dev_avg_pls_3SD)
        I_p_wave_areas_3_SA_node_dev_avg_pls_3SD_index = i;
    end
end
I_p_wave_areas_3_SA_node_tot_number_of_indices = length(I_p_wave_areas_3_SA_node_sorted);
%This is the approximation of the tot "area under curve." 
I_p_wave_areas_3_SA_node_quotient_within_1SD = (I_p_wave_areas_3_SA_node_dev_avg_pls_1SD_index - I_p_wave_areas_3_SA_nodedev_avg_min_1SD_index)/I_p_wave_areas_3_SA_node_tot_number_of_indices;
I_p_wave_areas_3_SA_node_quotient_within_2SD = (I_p_wave_areas_3_SA_node_dev_avg_pls_2SD_index - I_p_wave_areas_3_dev_avg_min_2SD_index)/I_p_wave_areas_3_SA_node_tot_number_of_indices;
I_p_wave_areas_3_SA_node_quotient_within_3SD = (I_p_wave_areas_3_SA_node_dev_avg_pls_3SD_index - I_p_wave_areas_3_SA_node_dev_avg_min_3SD_index)/I_p_wave_areas_3_SA_node_tot_number_of_indices;
if ((I_p_wave_areas_3_SA_node_quotient_within_1SD > 0.66 && I_p_wave_areas_3_SA_node_quotient_within_1SD < 0.70) && (I_p_wave_areas_3_SA_node_quotient_within_2SD > 0.93 && I_p_wave_areas_3_SA_node_quotient_within_2SD < 0.97) && (I_p_wave_areas_3_SA_node_quotient_within_3SD > 0.98 && I_p_wave_areas_3_SA_node_quotient_within_3SD < 1)) 
    isNormal(9,1) = 1;
end
if I_p_wave_areas_3_SA_node_quotient_within_1SD == 0
    isSpike(9,1) = 1;
end
I_p_wave_areas_3_crista_terminalis_sorted = sort(I_p_wave_areas_3_crista_terminalis);
I_p_wave_areas_3_crista_terminalis_dev_avg_min_1SD = area_3(length(I_p_wave_areas_3_crista_terminalis_sorted),1);
for i = 1:length(I_p_wave_areas_3_crista_terminalis_sorted)
    I_p_wave_areas_3_crista_terminalis_dev_avg_min_1SD(i) = I_p_wave_areas_3_crista_terminalis_sorted(i) - (I_p_wave_areas_3_crista_terminalis_average-I_p_wave_areas_3_crista_terminalis_SD);
end
I_p_wave_areas_3_crista_terminalis_dev_avg_min_1SD = abs(I_p_wave_areas_3_crista_terminalis_dev_avg_min_1SD);
for i = 1:length(I_p_wave_areas_3_crista_terminalis_sorted)
    if I_p_wave_areas_3_crista_terminalis_dev_avg_min_1SD(i) == min(I_p_wave_areas_3_crista_terminalis_dev_avg_min_1SD)
        I_p_wave_areas_3_crista_terminalis_dev_avg_min_1SD_index = i;
    end
end
I_p_wave_areas_3_crista_terminalis_dev_avg_pls_1SD = area_3(length(I_p_wave_areas_3_crista_terminalis_sorted),1);
for i = 1:length(I_p_wave_areas_3_crista_terminalis_sorted)
    I_p_wave_areas_3_crista_terminalis_dev_avg_pls_1SD(i) = I_p_wave_areas_3_crista_terminalis_sorted(i) - (I_p_wave_areas_3_crista_terminalis_average+I_p_wave_areas_3_crista_terminalis_SD);
end
I_p_wave_areas_3_crista_terminalis_dev_avg_pls_1SD = abs(I_p_wave_areas_3_crista_terminalis_dev_avg_pls_1SD);
for i = 1:length(I_p_wave_areas_3_crista_terminalis_sorted)
    if I_p_wave_areas_3_crista_terminalis_dev_avg_pls_1SD(i) == min(I_p_wave_areas_3_crista_terminalis_dev_avg_pls_1SD)
        I_p_wave_areas_3_crista_terminalis_dev_avg_pls_1SD_index = i;
    end
end
I_p_wave_areas_3_crista_terminalis_dev_avg_min_2SD = area_3(length(I_p_wave_areas_3_crista_terminalis_sorted),1);
for i = 1:length(I_p_wave_areas_3_crista_terminalis_sorted)
    I_p_wave_areas_3_crista_terminalis_dev_avg_min_2SD(i) = I_p_wave_areas_3_crista_terminalis_sorted(i) - (I_p_wave_areas_3_crista_terminalis_average-2*I_p_wave_areas_3_crista_terminalis_SD);
end
I_p_wave_areas_3_crista_terminalis_dev_avg_min_2SD = abs(I_p_wave_areas_3_crista_terminalis_dev_avg_min_2SD);
for i = 1:length(I_p_wave_areas_3_crista_terminalis_sorted)
    if I_p_wave_areas_3_crista_terminalis_dev_avg_min_2SD(i) == min(I_p_wave_areas_3_crista_terminalis_dev_avg_min_2SD)
        I_p_wave_areas_3_dev_avg_min_2SD_index = i;
    end
end
I_p_wave_areas_3_crista_terminalis_dev_avg_pls_2SD = area_3(length(I_p_wave_areas_3_crista_terminalis_sorted),1);
for i = 1:length(I_p_wave_areas_3_crista_terminalis_sorted)
    I_p_wave_areas_3_crista_terminalis_dev_avg_pls_2SD(i) = I_p_wave_areas_3_crista_terminalis_sorted(i) - (I_p_wave_areas_3_crista_terminalis_average+2*I_p_wave_areas_3_crista_terminalis_SD);
end
I_p_wave_areas_3_crista_terminalis_dev_avg_pls_2SD = abs(I_p_wave_areas_3_crista_terminalis_dev_avg_pls_2SD);
for i = 1:length(I_p_wave_areas_3_crista_terminalis_sorted)
    if I_p_wave_areas_3_crista_terminalis_dev_avg_pls_2SD(i) == min(I_p_wave_areas_3_crista_terminalis_dev_avg_pls_2SD)
        I_p_wave_areas_3_crista_terminalis_dev_avg_pls_2SD_index = i;
    end
end
I_p_wave_areas_3_crista_terminalis_dev_avg_min_3SD = area_3(length(I_p_wave_areas_3_crista_terminalis_sorted),1);
for i = 1:length(I_p_wave_areas_3_crista_terminalis_sorted)
    I_p_wave_areas_3_crista_terminalis_dev_avg_min_3SD(i) = I_p_wave_areas_3_crista_terminalis_sorted(i) - (I_p_wave_areas_3_crista_terminalis_average-3*I_p_wave_areas_3_crista_terminalis_SD);
end
I_p_wave_areas_3_crista_terminalis_dev_avg_min_3SD = abs(I_p_wave_areas_3_crista_terminalis_dev_avg_min_3SD);
for i = 1:length(I_p_wave_areas_3_crista_terminalis_sorted)
    if I_p_wave_areas_3_crista_terminalis_dev_avg_min_3SD(i) == min(I_p_wave_areas_3_crista_terminalis_dev_avg_min_3SD)
        I_p_wave_areas_3_crista_terminalis_dev_avg_min_3SD_index = i;
    end
end
I_p_wave_areas_3_crista_terminalis_dev_avg_pls_3SD = area_3(length(I_p_wave_areas_3_crista_terminalis_sorted),1);
for i = 1:length(I_p_wave_areas_3_crista_terminalis_sorted)
    I_p_wave_areas_3_crista_terminalis_dev_avg_pls_3SD(i) = I_p_wave_areas_3_crista_terminalis_sorted(i) - (I_p_wave_areas_3_crista_terminalis_average+3*I_p_wave_areas_3_crista_terminalis_SD);
end
I_p_wave_areas_3_crista_terminalis_dev_avg_pls_3SD = abs(I_p_wave_areas_3_crista_terminalis_dev_avg_pls_3SD);
for i = 1:length(I_p_wave_areas_3_crista_terminalis_sorted)
    if I_p_wave_areas_3_crista_terminalis_dev_avg_pls_3SD(i) == min(I_p_wave_areas_3_crista_terminalis_dev_avg_pls_3SD)
        I_p_wave_areas_3_crista_terminalis_dev_avg_pls_3SD_index = i;
    end
end
I_p_wave_areas_3_crista_terminalis_tot_number_of_indices = length(I_p_wave_areas_3_crista_terminalis_sorted);
%This is the approximation of the tot "area under curve." 
I_p_wave_areas_3_crista_terminalis_quotient_within_1SD = (I_p_wave_areas_3_crista_terminalis_dev_avg_pls_1SD_index - I_p_wave_areas_3_crista_terminalisdev_avg_min_1SD_index)/I_p_wave_areas_3_crista_terminalis_tot_number_of_indices;
I_p_wave_areas_3_crista_terminalis_quotient_within_2SD = (I_p_wave_areas_3_crista_terminalis_dev_avg_pls_2SD_index - I_p_wave_areas_3_dev_avg_min_2SD_index)/I_p_wave_areas_3_crista_terminalis_tot_number_of_indices;
I_p_wave_areas_3_crista_terminalis_quotient_within_3SD = (I_p_wave_areas_3_crista_terminalis_dev_avg_pls_3SD_index - I_p_wave_areas_3_crista_terminalis_dev_avg_min_3SD_index)/I_p_wave_areas_3_crista_terminalis_tot_number_of_indices;
if (I_p_wave_areas_3_crista_terminalis_quotient_within_1SD > 0.66 && I_p_wave_areas_3_crista_terminalis_quotient_within_1SD < 0.70) && (I_p_wave_areas_3_crista_terminalis_quotient_within_2SD > 0.93 && I_p_wave_areas_3_crista_terminalis_quotient_within_2SD < 0.97) && (I_p_wave_areas_3_crista_terminalis_quotient_within_3SD > 0.98 && I_p_wave_areas_3_crista_terminalis_quotient_within_3SD < 1)
    isNormal(9,2) = 1;
end
if I_p_wave_areas_3_crista_terminalis_quotient_within_1SD == 0
    isSpike(9,2) = 1;
end
I_p_wave_areas_3_tricuspid_annulus_sorted = sort(I_p_wave_areas_3_tricuspid_annulus);
I_p_wave_areas_3_tricuspid_annulus_dev_avg_min_1SD = area_3(length(I_p_wave_areas_3_tricuspid_annulus_sorted),1);
for i = 1:length(I_p_wave_areas_3_tricuspid_annulus_sorted)
    I_p_wave_areas_3_tricuspid_annulus_dev_avg_min_1SD(i) = I_p_wave_areas_3_tricuspid_annulus_sorted(i) - (I_p_wave_areas_3_tricuspid_annulus_average-I_p_wave_areas_3_tricuspid_annulus_SD);
end
I_p_wave_areas_3_tricuspid_annulus_dev_avg_min_1SD = abs(I_p_wave_areas_3_tricuspid_annulus_dev_avg_min_1SD);
for i = 1:length(I_p_wave_areas_3_tricuspid_annulus_sorted)
    if I_p_wave_areas_3_tricuspid_annulus_dev_avg_min_1SD(i) == min(I_p_wave_areas_3_tricuspid_annulus_dev_avg_min_1SD)
        I_p_wave_areas_3_tricuspid_annulus_dev_avg_min_1SD_index = i;
    end
end
I_p_wave_areas_3_tricuspid_annulus_dev_avg_pls_1SD = area_3(length(I_p_wave_areas_3_tricuspid_annulus_sorted),1);
for i = 1:length(I_p_wave_areas_3_tricuspid_annulus_sorted)
    I_p_wave_areas_3_tricuspid_annulus_dev_avg_pls_1SD(i) = I_p_wave_areas_3_tricuspid_annulus_sorted(i) - (I_p_wave_areas_3_tricuspid_annulus_average+I_p_wave_areas_3_tricuspid_annulus_SD);
end
I_p_wave_areas_3_tricuspid_annulus_dev_avg_pls_1SD = abs(I_p_wave_areas_3_tricuspid_annulus_dev_avg_pls_1SD);
for i = 1:length(I_p_wave_areas_3_tricuspid_annulus_sorted)
    if I_p_wave_areas_3_tricuspid_annulus_dev_avg_pls_1SD(i) == min(I_p_wave_areas_3_tricuspid_annulus_dev_avg_pls_1SD)
        I_p_wave_areas_3_tricuspid_annulus_dev_avg_pls_1SD_index = i;
    end
end
I_p_wave_areas_3_tricuspid_annulus_dev_avg_min_2SD = area_3(length(I_p_wave_areas_3_tricuspid_annulus_sorted),1);
for i = 1:length(I_p_wave_areas_3_tricuspid_annulus_sorted)
    I_p_wave_areas_3_tricuspid_annulus_dev_avg_min_2SD(i) = I_p_wave_areas_3_tricuspid_annulus_sorted(i) - (I_p_wave_areas_3_tricuspid_annulus_average-2*I_p_wave_areas_3_tricuspid_annulus_SD);
end
I_p_wave_areas_3_tricuspid_annulus_dev_avg_min_2SD = abs(I_p_wave_areas_3_tricuspid_annulus_dev_avg_min_2SD);
for i = 1:length(I_p_wave_areas_3_tricuspid_annulus_sorted)
    if I_p_wave_areas_3_tricuspid_annulus_dev_avg_min_2SD(i) == min(I_p_wave_areas_3_tricuspid_annulus_dev_avg_min_2SD)
        I_p_wave_areas_3_dev_avg_min_2SD_index = i;
    end
end
I_p_wave_areas_3_tricuspid_annulus_dev_avg_pls_2SD = area_3(length(I_p_wave_areas_3_tricuspid_annulus_sorted),1);
for i = 1:length(I_p_wave_areas_3_tricuspid_annulus_sorted)
    I_p_wave_areas_3_tricuspid_annulus_dev_avg_pls_2SD(i) = I_p_wave_areas_3_tricuspid_annulus_sorted(i) - (I_p_wave_areas_3_tricuspid_annulus_average+2*I_p_wave_areas_3_tricuspid_annulus_SD);
end
I_p_wave_areas_3_tricuspid_annulus_dev_avg_pls_2SD = abs(I_p_wave_areas_3_tricuspid_annulus_dev_avg_pls_2SD);
for i = 1:length(I_p_wave_areas_3_tricuspid_annulus_sorted)
    if I_p_wave_areas_3_tricuspid_annulus_dev_avg_pls_2SD(i) == min(I_p_wave_areas_3_tricuspid_annulus_dev_avg_pls_2SD)
        I_p_wave_areas_3_tricuspid_annulus_dev_avg_pls_2SD_index = i;
    end
end
I_p_wave_areas_3_tricuspid_annulus_dev_avg_min_3SD = area_3(length(I_p_wave_areas_3_tricuspid_annulus_sorted),1);
for i = 1:length(I_p_wave_areas_3_tricuspid_annulus_sorted)
    I_p_wave_areas_3_tricuspid_annulus_dev_avg_min_3SD(i) = I_p_wave_areas_3_tricuspid_annulus_sorted(i) - (I_p_wave_areas_3_tricuspid_annulus_average-3*I_p_wave_areas_3_tricuspid_annulus_SD);
end
I_p_wave_areas_3_tricuspid_annulus_dev_avg_min_3SD = abs(I_p_wave_areas_3_tricuspid_annulus_dev_avg_min_3SD);
for i = 1:length(I_p_wave_areas_3_tricuspid_annulus_sorted)
    if I_p_wave_areas_3_tricuspid_annulus_dev_avg_min_3SD(i) == min(I_p_wave_areas_3_tricuspid_annulus_dev_avg_min_3SD)
        I_p_wave_areas_3_tricuspid_annulus_dev_avg_min_3SD_index = i;
    end
end
I_p_wave_areas_3_tricuspid_annulus_dev_avg_pls_3SD = area_3(length(I_p_wave_areas_3_tricuspid_annulus_sorted),1);
for i = 1:length(I_p_wave_areas_3_tricuspid_annulus_sorted)
    I_p_wave_areas_3_tricuspid_annulus_dev_avg_pls_3SD(i) = I_p_wave_areas_3_tricuspid_annulus_sorted(i) - (I_p_wave_areas_3_tricuspid_annulus_average+3*I_p_wave_areas_3_tricuspid_annulus_SD);
end
I_p_wave_areas_3_tricuspid_annulus_dev_avg_pls_3SD = abs(I_p_wave_areas_3_tricuspid_annulus_dev_avg_pls_3SD);
for i = 1:length(I_p_wave_areas_3_tricuspid_annulus_sorted)
    if I_p_wave_areas_3_tricuspid_annulus_dev_avg_pls_3SD(i) == min(I_p_wave_areas_3_tricuspid_annulus_dev_avg_pls_3SD)
        I_p_wave_areas_3_tricuspid_annulus_dev_avg_pls_3SD_index = i;
    end
end
I_p_wave_areas_3_tricuspid_annulus_tot_number_of_indices = length(I_p_wave_areas_3_tricuspid_annulus_sorted);
%This is the approximation of the tot "area under curve." 
I_p_wave_areas_3_tricuspid_annulus_quotient_within_1SD = (I_p_wave_areas_3_tricuspid_annulus_dev_avg_pls_1SD_index - I_p_wave_areas_3_tricuspid_annulusdev_avg_min_1SD_index)/I_p_wave_areas_3_tricuspid_annulus_tot_number_of_indices;
I_p_wave_areas_3_tricuspid_annulus_quotient_within_2SD = (I_p_wave_areas_3_tricuspid_annulus_dev_avg_pls_2SD_index - I_p_wave_areas_3_dev_avg_min_2SD_index)/I_p_wave_areas_3_tricuspid_annulus_tot_number_of_indices;
I_p_wave_areas_3_tricuspid_annulus_quotient_within_3SD = (I_p_wave_areas_3_tricuspid_annulus_dev_avg_pls_3SD_index - I_p_wave_areas_3_tricuspid_annulus_dev_avg_min_3SD_index)/I_p_wave_areas_3_tricuspid_annulus_tot_number_of_indices;
if (I_p_wave_areas_3_tricuspid_annulus_quotient_within_1SD > 0.66 && I_p_wave_areas_3_tricuspid_annulus_quotient_within_1SD < 0.70) && (I_p_wave_areas_3_tricuspid_annulus_quotient_within_2SD > 0.93 && I_p_wave_areas_3_tricuspid_annulus_quotient_within_2SD < 0.97) && (I_p_wave_areas_3_tricuspid_annulus_quotient_within_3SD > 0.98 && I_p_wave_areas_3_tricuspid_annulus_quotient_within_3SD < 1)
    isNormal(9,3) = 1;
end
if I_p_wave_areas_3_tricuspid_annulus_quotient_within_1SD == 0
    isSpike(9,3) = 1;
end
I_p_wave_areas_3_coronary_sinus_sorted = sort(I_p_wave_areas_3_coronary_sinus);
I_p_wave_areas_3_coronary_sinus_dev_avg_min_1SD = area_3(length(I_p_wave_areas_3_coronary_sinus_sorted),1);
for i = 1:length(I_p_wave_areas_3_coronary_sinus_sorted)
    I_p_wave_areas_3_coronary_sinus_dev_avg_min_1SD(i) = I_p_wave_areas_3_coronary_sinus_sorted(i) - (I_p_wave_areas_3_coronary_sinus_average-I_p_wave_areas_3_coronary_sinus_SD);
end
I_p_wave_areas_3_coronary_sinus_dev_avg_min_1SD = abs(I_p_wave_areas_3_coronary_sinus_dev_avg_min_1SD);
for i = 1:length(I_p_wave_areas_3_coronary_sinus_sorted)
    if I_p_wave_areas_3_coronary_sinus_dev_avg_min_1SD(i) == min(I_p_wave_areas_3_coronary_sinus_dev_avg_min_1SD)
        I_p_wave_areas_3_coronary_sinus_dev_avg_min_1SD_index = i;
    end
end
I_p_wave_areas_3_coronary_sinus_dev_avg_pls_1SD = area_3(length(I_p_wave_areas_3_coronary_sinus_sorted),1);
for i = 1:length(I_p_wave_areas_3_coronary_sinus_sorted)
    I_p_wave_areas_3_coronary_sinus_dev_avg_pls_1SD(i) = I_p_wave_areas_3_coronary_sinus_sorted(i) - (I_p_wave_areas_3_coronary_sinus_average+I_p_wave_areas_3_coronary_sinus_SD);
end
I_p_wave_areas_3_coronary_sinus_dev_avg_pls_1SD = abs(I_p_wave_areas_3_coronary_sinus_dev_avg_pls_1SD);
for i = 1:length(I_p_wave_areas_3_coronary_sinus_sorted)
    if I_p_wave_areas_3_coronary_sinus_dev_avg_pls_1SD(i) == min(I_p_wave_areas_3_coronary_sinus_dev_avg_pls_1SD)
        I_p_wave_areas_3_coronary_sinus_dev_avg_pls_1SD_index = i;
    end
end
I_p_wave_areas_3_coronary_sinus_dev_avg_min_2SD = area_3(length(I_p_wave_areas_3_coronary_sinus_sorted),1);
for i = 1:length(I_p_wave_areas_3_coronary_sinus_sorted)
    I_p_wave_areas_3_coronary_sinus_dev_avg_min_2SD(i) = I_p_wave_areas_3_coronary_sinus_sorted(i) - (I_p_wave_areas_3_coronary_sinus_average-2*I_p_wave_areas_3_coronary_sinus_SD);
end
I_p_wave_areas_3_coronary_sinus_dev_avg_min_2SD = abs(I_p_wave_areas_3_coronary_sinus_dev_avg_min_2SD);
for i = 1:length(I_p_wave_areas_3_coronary_sinus_sorted)
    if I_p_wave_areas_3_coronary_sinus_dev_avg_min_2SD(i) == min(I_p_wave_areas_3_coronary_sinus_dev_avg_min_2SD)
        I_p_wave_areas_3_dev_avg_min_2SD_index = i;
    end
end
I_p_wave_areas_3_coronary_sinus_dev_avg_pls_2SD = area_3(length(I_p_wave_areas_3_coronary_sinus_sorted),1);
for i = 1:length(I_p_wave_areas_3_coronary_sinus_sorted)
    I_p_wave_areas_3_coronary_sinus_dev_avg_pls_2SD(i) = I_p_wave_areas_3_coronary_sinus_sorted(i) - (I_p_wave_areas_3_coronary_sinus_average+2*I_p_wave_areas_3_coronary_sinus_SD);
end
I_p_wave_areas_3_coronary_sinus_dev_avg_pls_2SD = abs(I_p_wave_areas_3_coronary_sinus_dev_avg_pls_2SD);
for i = 1:length(I_p_wave_areas_3_coronary_sinus_sorted)
    if I_p_wave_areas_3_coronary_sinus_dev_avg_pls_2SD(i) == min(I_p_wave_areas_3_coronary_sinus_dev_avg_pls_2SD)
        I_p_wave_areas_3_coronary_sinus_dev_avg_pls_2SD_index = i;
    end
end
I_p_wave_areas_3_coronary_sinus_dev_avg_min_3SD = area_3(length(I_p_wave_areas_3_coronary_sinus_sorted),1);
for i = 1:length(I_p_wave_areas_3_coronary_sinus_sorted)
    I_p_wave_areas_3_coronary_sinus_dev_avg_min_3SD(i) = I_p_wave_areas_3_coronary_sinus_sorted(i) - (I_p_wave_areas_3_coronary_sinus_average-3*I_p_wave_areas_3_coronary_sinus_SD);
end
I_p_wave_areas_3_coronary_sinus_dev_avg_min_3SD = abs(I_p_wave_areas_3_coronary_sinus_dev_avg_min_3SD);
for i = 1:length(I_p_wave_areas_3_coronary_sinus_sorted)
    if I_p_wave_areas_3_coronary_sinus_dev_avg_min_3SD(i) == min(I_p_wave_areas_3_coronary_sinus_dev_avg_min_3SD)
        I_p_wave_areas_3_coronary_sinus_dev_avg_min_3SD_index = i;
    end
end
I_p_wave_areas_3_coronary_sinus_dev_avg_pls_3SD = area_3(length(I_p_wave_areas_3_coronary_sinus_sorted),1);
for i = 1:length(I_p_wave_areas_3_coronary_sinus_sorted)
    I_p_wave_areas_3_coronary_sinus_dev_avg_pls_3SD(i) = I_p_wave_areas_3_coronary_sinus_sorted(i) - (I_p_wave_areas_3_coronary_sinus_average+3*I_p_wave_areas_3_coronary_sinus_SD);
end
I_p_wave_areas_3_coronary_sinus_dev_avg_pls_3SD = abs(I_p_wave_areas_3_coronary_sinus_dev_avg_pls_3SD);
for i = 1:length(I_p_wave_areas_3_coronary_sinus_sorted)
    if I_p_wave_areas_3_coronary_sinus_dev_avg_pls_3SD(i) == min(I_p_wave_areas_3_coronary_sinus_dev_avg_pls_3SD)
        I_p_wave_areas_3_coronary_sinus_dev_avg_pls_3SD_index = i;
    end
end
I_p_wave_areas_3_coronary_sinus_tot_number_of_indices = length(I_p_wave_areas_3_coronary_sinus_sorted);
%This is the approximation of the tot "area under curve." 
I_p_wave_areas_3_coronary_sinus_quotient_within_1SD = (I_p_wave_areas_3_coronary_sinus_dev_avg_pls_1SD_index - I_p_wave_areas_3_coronary_sinusdev_avg_min_1SD_index)/I_p_wave_areas_3_coronary_sinus_tot_number_of_indices;
I_p_wave_areas_3_coronary_sinus_quotient_within_2SD = (I_p_wave_areas_3_coronary_sinus_dev_avg_pls_2SD_index - I_p_wave_areas_3_dev_avg_min_2SD_index)/I_p_wave_areas_3_coronary_sinus_tot_number_of_indices;
I_p_wave_areas_3_coronary_sinus_quotient_within_3SD = (I_p_wave_areas_3_coronary_sinus_dev_avg_pls_3SD_index - I_p_wave_areas_3_coronary_sinus_dev_avg_min_3SD_index)/I_p_wave_areas_3_coronary_sinus_tot_number_of_indices;
if (I_p_wave_areas_3_coronary_sinus_quotient_within_1SD > 0.66 && I_p_wave_areas_3_coronary_sinus_quotient_within_1SD < 0.70) && (I_p_wave_areas_3_coronary_sinus_quotient_within_2SD > 0.93 && I_p_wave_areas_3_coronary_sinus_quotient_within_2SD < 0.97) && (I_p_wave_areas_3_coronary_sinus_quotient_within_3SD > 0.98 && I_p_wave_areas_3_coronary_sinus_quotient_within_3SD < 1)
    isNormal(9,4) = 1;
end
if I_p_wave_areas_3_coronary_sinus_quotient_within_1SD == 0
    isSpike(9,4) = 1;
end
I_p_wave_areas_3_ostium_sorted = sort(I_p_wave_areas_3_ostium);
I_p_wave_areas_3_ostium_dev_avg_min_1SD = area_3(length(I_p_wave_areas_3_ostium_sorted),1);
for i = 1:length(I_p_wave_areas_3_ostium_sorted)
    I_p_wave_areas_3_ostium_dev_avg_min_1SD(i) = I_p_wave_areas_3_ostium_sorted(i) - (I_p_wave_areas_3_ostium_average-I_p_wave_areas_3_ostium_SD);
end
I_p_wave_areas_3_ostium_dev_avg_min_1SD = abs(I_p_wave_areas_3_ostium_dev_avg_min_1SD);
for i = 1:length(I_p_wave_areas_3_ostium_sorted)
    if I_p_wave_areas_3_ostium_dev_avg_min_1SD(i) == min(I_p_wave_areas_3_ostium_dev_avg_min_1SD)
        I_p_wave_areas_3_ostium_dev_avg_min_1SD_index = i;
    end
end
I_p_wave_areas_3_ostium_dev_avg_pls_1SD = area_3(length(I_p_wave_areas_3_ostium_sorted),1);
for i = 1:length(I_p_wave_areas_3_ostium_sorted)
    I_p_wave_areas_3_ostium_dev_avg_pls_1SD(i) = I_p_wave_areas_3_ostium_sorted(i) - (I_p_wave_areas_3_ostium_average+I_p_wave_areas_3_ostium_SD);
end
I_p_wave_areas_3_ostium_dev_avg_pls_1SD = abs(I_p_wave_areas_3_ostium_dev_avg_pls_1SD);
for i = 1:length(I_p_wave_areas_3_ostium_sorted)
    if I_p_wave_areas_3_ostium_dev_avg_pls_1SD(i) == min(I_p_wave_areas_3_ostium_dev_avg_pls_1SD)
        I_p_wave_areas_3_ostium_dev_avg_pls_1SD_index = i;
    end
end
I_p_wave_areas_3_ostium_dev_avg_min_2SD = area_3(length(I_p_wave_areas_3_ostium_sorted),1);
for i = 1:length(I_p_wave_areas_3_ostium_sorted)
    I_p_wave_areas_3_ostium_dev_avg_min_2SD(i) = I_p_wave_areas_3_ostium_sorted(i) - (I_p_wave_areas_3_ostium_average-2*I_p_wave_areas_3_ostium_SD);
end
I_p_wave_areas_3_ostium_dev_avg_min_2SD = abs(I_p_wave_areas_3_ostium_dev_avg_min_2SD);
for i = 1:length(I_p_wave_areas_3_ostium_sorted)
    if I_p_wave_areas_3_ostium_dev_avg_min_2SD(i) == min(I_p_wave_areas_3_ostium_dev_avg_min_2SD)
        I_p_wave_areas_3_dev_avg_min_2SD_index = i;
    end
end
I_p_wave_areas_3_ostium_dev_avg_pls_2SD = area_3(length(I_p_wave_areas_3_ostium_sorted),1);
for i = 1:length(I_p_wave_areas_3_ostium_sorted)
    I_p_wave_areas_3_ostium_dev_avg_pls_2SD(i) = I_p_wave_areas_3_ostium_sorted(i) - (I_p_wave_areas_3_ostium_average+2*I_p_wave_areas_3_ostium_SD);
end
I_p_wave_areas_3_ostium_dev_avg_pls_2SD = abs(I_p_wave_areas_3_ostium_dev_avg_pls_2SD);
for i = 1:length(I_p_wave_areas_3_ostium_sorted)
    if I_p_wave_areas_3_ostium_dev_avg_pls_2SD(i) == min(I_p_wave_areas_3_ostium_dev_avg_pls_2SD)
        I_p_wave_areas_3_ostium_dev_avg_pls_2SD_index = i;
    end
end
I_p_wave_areas_3_ostium_dev_avg_min_3SD = area_3(length(I_p_wave_areas_3_ostium_sorted),1);
for i = 1:length(I_p_wave_areas_3_ostium_sorted)
    I_p_wave_areas_3_ostium_dev_avg_min_3SD(i) = I_p_wave_areas_3_ostium_sorted(i) - (I_p_wave_areas_3_ostium_average-3*I_p_wave_areas_3_ostium_SD);
end
I_p_wave_areas_3_ostium_dev_avg_min_3SD = abs(I_p_wave_areas_3_ostium_dev_avg_min_3SD);
for i = 1:length(I_p_wave_areas_3_ostium_sorted)
    if I_p_wave_areas_3_ostium_dev_avg_min_3SD(i) == min(I_p_wave_areas_3_ostium_dev_avg_min_3SD)
        I_p_wave_areas_3_ostium_dev_avg_min_3SD_index = i;
    end
end
I_p_wave_areas_3_ostium_dev_avg_pls_3SD = area_3(length(I_p_wave_areas_3_ostium_sorted),1);
for i = 1:length(I_p_wave_areas_3_ostium_sorted)
    I_p_wave_areas_3_ostium_dev_avg_pls_3SD(i) = I_p_wave_areas_3_ostium_sorted(i) - (I_p_wave_areas_3_ostium_average+3*I_p_wave_areas_3_ostium_SD);
end
I_p_wave_areas_3_ostium_dev_avg_pls_3SD = abs(I_p_wave_areas_3_ostium_dev_avg_pls_3SD);
for i = 1:length(I_p_wave_areas_3_ostium_sorted)
    if I_p_wave_areas_3_ostium_dev_avg_pls_3SD(i) == min(I_p_wave_areas_3_ostium_dev_avg_pls_3SD)
        I_p_wave_areas_3_ostium_dev_avg_pls_3SD_index = i;
    end
end
I_p_wave_areas_3_ostium_tot_number_of_indices = length(I_p_wave_areas_3_ostium_sorted);
%This is the approximation of the tot "area under curve." 
I_p_wave_areas_3_ostium_quotient_within_1SD = (I_p_wave_areas_3_ostium_dev_avg_pls_1SD_index - I_p_wave_areas_3_ostiumdev_avg_min_1SD_index)/I_p_wave_areas_3_ostium_tot_number_of_indices;
I_p_wave_areas_3_ostium_quotient_within_2SD = (I_p_wave_areas_3_ostium_dev_avg_pls_2SD_index - I_p_wave_areas_3_dev_avg_min_2SD_index)/I_p_wave_areas_3_ostium_tot_number_of_indices;
I_p_wave_areas_3_ostium_quotient_within_3SD = (I_p_wave_areas_3_ostium_dev_avg_pls_3SD_index - I_p_wave_areas_3_ostium_dev_avg_min_3SD_index)/I_p_wave_areas_3_ostium_tot_number_of_indices;
if (I_p_wave_areas_3_ostium_quotient_within_1SD > 0.66 && I_p_wave_areas_3_ostium_quotient_within_1SD < 0.70) && (I_p_wave_areas_3_ostium_quotient_within_2SD > 0.93 && I_p_wave_areas_3_ostium_quotient_within_2SD < 0.97) && (I_p_wave_areas_3_ostium_quotient_within_3SD > 0.98 && I_p_wave_areas_3_ostium_quotient_within_3SD < 1)
    isNormal(9,5) = 1;
end
if I_p_wave_areas_3_ostium_quotient_within_1SD == 0
    isSpike(9,5) = 1;
end
I_p_wave_areas_3_perinodal_sorted = sort(I_p_wave_areas_3_perinodal);
I_p_wave_areas_3_perinodal_dev_avg_min_1SD = area_3(length(I_p_wave_areas_3_perinodal_sorted),1);
for i = 1:length(I_p_wave_areas_3_perinodal_sorted)
    I_p_wave_areas_3_perinodal_dev_avg_min_1SD(i) = I_p_wave_areas_3_perinodal_sorted(i) - (I_p_wave_areas_3_perinodal_average-I_p_wave_areas_3_perinodal_SD);
end
I_p_wave_areas_3_perinodal_dev_avg_min_1SD = abs(I_p_wave_areas_3_perinodal_dev_avg_min_1SD);
for i = 1:length(I_p_wave_areas_3_perinodal_sorted)
    if I_p_wave_areas_3_perinodal_dev_avg_min_1SD(i) == min(I_p_wave_areas_3_perinodal_dev_avg_min_1SD)
        I_p_wave_areas_3_perinodal_dev_avg_min_1SD_index = i;
    end
end
I_p_wave_areas_3_perinodal_dev_avg_pls_1SD = area_3(length(I_p_wave_areas_3_perinodal_sorted),1);
for i = 1:length(I_p_wave_areas_3_perinodal_sorted)
    I_p_wave_areas_3_perinodal_dev_avg_pls_1SD(i) = I_p_wave_areas_3_perinodal_sorted(i) - (I_p_wave_areas_3_perinodal_average+I_p_wave_areas_3_perinodal_SD);
end
I_p_wave_areas_3_perinodal_dev_avg_pls_1SD = abs(I_p_wave_areas_3_perinodal_dev_avg_pls_1SD);
for i = 1:length(I_p_wave_areas_3_perinodal_sorted)
    if I_p_wave_areas_3_perinodal_dev_avg_pls_1SD(i) == min(I_p_wave_areas_3_perinodal_dev_avg_pls_1SD)
        I_p_wave_areas_3_perinodal_dev_avg_pls_1SD_index = i;
    end
end
I_p_wave_areas_3_perinodal_dev_avg_min_2SD = area_3(length(I_p_wave_areas_3_perinodal_sorted),1);
for i = 1:length(I_p_wave_areas_3_perinodal_sorted)
    I_p_wave_areas_3_perinodal_dev_avg_min_2SD(i) = I_p_wave_areas_3_perinodal_sorted(i) - (I_p_wave_areas_3_perinodal_average-2*I_p_wave_areas_3_perinodal_SD);
end
I_p_wave_areas_3_perinodal_dev_avg_min_2SD = abs(I_p_wave_areas_3_perinodal_dev_avg_min_2SD);
for i = 1:length(I_p_wave_areas_3_perinodal_sorted)
    if I_p_wave_areas_3_perinodal_dev_avg_min_2SD(i) == min(I_p_wave_areas_3_perinodal_dev_avg_min_2SD)
        I_p_wave_areas_3_dev_avg_min_2SD_index = i;
    end
end
I_p_wave_areas_3_perinodal_dev_avg_pls_2SD = area_3(length(I_p_wave_areas_3_perinodal_sorted),1);
for i = 1:length(I_p_wave_areas_3_perinodal_sorted)
    I_p_wave_areas_3_perinodal_dev_avg_pls_2SD(i) = I_p_wave_areas_3_perinodal_sorted(i) - (I_p_wave_areas_3_perinodal_average+2*I_p_wave_areas_3_perinodal_SD);
end
I_p_wave_areas_3_perinodal_dev_avg_pls_2SD = abs(I_p_wave_areas_3_perinodal_dev_avg_pls_2SD);
for i = 1:length(I_p_wave_areas_3_perinodal_sorted)
    if I_p_wave_areas_3_perinodal_dev_avg_pls_2SD(i) == min(I_p_wave_areas_3_perinodal_dev_avg_pls_2SD)
        I_p_wave_areas_3_perinodal_dev_avg_pls_2SD_index = i;
    end
end
I_p_wave_areas_3_perinodal_dev_avg_min_3SD = area_3(length(I_p_wave_areas_3_perinodal_sorted),1);
for i = 1:length(I_p_wave_areas_3_perinodal_sorted)
    I_p_wave_areas_3_perinodal_dev_avg_min_3SD(i) = I_p_wave_areas_3_perinodal_sorted(i) - (I_p_wave_areas_3_perinodal_average-3*I_p_wave_areas_3_perinodal_SD);
end
I_p_wave_areas_3_perinodal_dev_avg_min_3SD = abs(I_p_wave_areas_3_perinodal_dev_avg_min_3SD);
for i = 1:length(I_p_wave_areas_3_perinodal_sorted)
    if I_p_wave_areas_3_perinodal_dev_avg_min_3SD(i) == min(I_p_wave_areas_3_perinodal_dev_avg_min_3SD)
        I_p_wave_areas_3_perinodal_dev_avg_min_3SD_index = i;
    end
end
I_p_wave_areas_3_perinodal_dev_avg_pls_3SD = area_3(length(I_p_wave_areas_3_perinodal_sorted),1);
for i = 1:length(I_p_wave_areas_3_perinodal_sorted)
    I_p_wave_areas_3_perinodal_dev_avg_pls_3SD(i) = I_p_wave_areas_3_perinodal_sorted(i) - (I_p_wave_areas_3_perinodal_average+3*I_p_wave_areas_3_perinodal_SD);
end
I_p_wave_areas_3_perinodal_dev_avg_pls_3SD = abs(I_p_wave_areas_3_perinodal_dev_avg_pls_3SD);
for i = 1:length(I_p_wave_areas_3_perinodal_sorted)
    if I_p_wave_areas_3_perinodal_dev_avg_pls_3SD(i) == min(I_p_wave_areas_3_perinodal_dev_avg_pls_3SD)
        I_p_wave_areas_3_perinodal_dev_avg_pls_3SD_index = i;
    end
end
I_p_wave_areas_3_perinodal_tot_number_of_indices = length(I_p_wave_areas_3_perinodal_sorted);
%This is the approximation of the tot "area under curve." 
I_p_wave_areas_3_perinodal_quotient_within_1SD = (I_p_wave_areas_3_perinodal_dev_avg_pls_1SD_index - I_p_wave_areas_3_perinodaldev_avg_min_1SD_index)/I_p_wave_areas_3_perinodal_tot_number_of_indices;
I_p_wave_areas_3_perinodal_quotient_within_2SD = (I_p_wave_areas_3_perinodal_dev_avg_pls_2SD_index - I_p_wave_areas_3_dev_avg_min_2SD_index)/I_p_wave_areas_3_perinodal_tot_number_of_indices;
I_p_wave_areas_3_perinodal_quotient_within_3SD = (I_p_wave_areas_3_perinodal_dev_avg_pls_3SD_index - I_p_wave_areas_3_perinodal_dev_avg_min_3SD_index)/I_p_wave_areas_3_perinodal_tot_number_of_indices;
if (I_p_wave_areas_3_perinodal_quotient_within_1SD > 0.66 && I_p_wave_areas_3_perinodal_quotient_within_1SD < 0.70) && (I_p_wave_areas_3_perinodal_quotient_within_2SD > 0.93 && I_p_wave_areas_3_perinodal_quotient_within_2SD < 0.97) && (I_p_wave_areas_3_perinodal_quotient_within_3SD > 0.98 && I_p_wave_areas_3_perinodal_quotient_within_3SD < 1)
    isNormal(9,6) = 1;
end
if I_p_wave_areas_3_perinodal_quotient_within_1SD == 0
    isSpike(9,6) = 1;
end
I_p_wave_areas_3_right_septum_sorted = sort(I_p_wave_areas_3_right_septum);
I_p_wave_areas_3_right_septum_dev_avg_min_1SD = area_3(length(I_p_wave_areas_3_right_septum_sorted),1);
for i = 1:length(I_p_wave_areas_3_right_septum_sorted)
    I_p_wave_areas_3_right_septum_dev_avg_min_1SD(i) = I_p_wave_areas_3_right_septum_sorted(i) - (I_p_wave_areas_3_right_septum_average-I_p_wave_areas_3_right_septum_SD);
end
I_p_wave_areas_3_right_septum_dev_avg_min_1SD = abs(I_p_wave_areas_3_right_septum_dev_avg_min_1SD);
for i = 1:length(I_p_wave_areas_3_right_septum_sorted)
    if I_p_wave_areas_3_right_septum_dev_avg_min_1SD(i) == min(I_p_wave_areas_3_right_septum_dev_avg_min_1SD)
        I_p_wave_areas_3_right_septum_dev_avg_min_1SD_index = i;
    end
end
I_p_wave_areas_3_right_septum_dev_avg_pls_1SD = area_3(length(I_p_wave_areas_3_right_septum_sorted),1);
for i = 1:length(I_p_wave_areas_3_right_septum_sorted)
    I_p_wave_areas_3_right_septum_dev_avg_pls_1SD(i) = I_p_wave_areas_3_right_septum_sorted(i) - (I_p_wave_areas_3_right_septum_average+I_p_wave_areas_3_right_septum_SD);
end
I_p_wave_areas_3_right_septum_dev_avg_pls_1SD = abs(I_p_wave_areas_3_right_septum_dev_avg_pls_1SD);
for i = 1:length(I_p_wave_areas_3_right_septum_sorted)
    if I_p_wave_areas_3_right_septum_dev_avg_pls_1SD(i) == min(I_p_wave_areas_3_right_septum_dev_avg_pls_1SD)
        I_p_wave_areas_3_right_septum_dev_avg_pls_1SD_index = i;
    end
end
I_p_wave_areas_3_right_septum_dev_avg_min_2SD = area_3(length(I_p_wave_areas_3_right_septum_sorted),1);
for i = 1:length(I_p_wave_areas_3_right_septum_sorted)
    I_p_wave_areas_3_right_septum_dev_avg_min_2SD(i) = I_p_wave_areas_3_right_septum_sorted(i) - (I_p_wave_areas_3_right_septum_average-2*I_p_wave_areas_3_right_septum_SD);
end
I_p_wave_areas_3_right_septum_dev_avg_min_2SD = abs(I_p_wave_areas_3_right_septum_dev_avg_min_2SD);
for i = 1:length(I_p_wave_areas_3_right_septum_sorted)
    if I_p_wave_areas_3_right_septum_dev_avg_min_2SD(i) == min(I_p_wave_areas_3_right_septum_dev_avg_min_2SD)
        I_p_wave_areas_3_dev_avg_min_2SD_index = i;
    end
end
I_p_wave_areas_3_right_septum_dev_avg_pls_2SD = area_3(length(I_p_wave_areas_3_right_septum_sorted),1);
for i = 1:length(I_p_wave_areas_3_right_septum_sorted)
    I_p_wave_areas_3_right_septum_dev_avg_pls_2SD(i) = I_p_wave_areas_3_right_septum_sorted(i) - (I_p_wave_areas_3_right_septum_average+2*I_p_wave_areas_3_right_septum_SD);
end
I_p_wave_areas_3_right_septum_dev_avg_pls_2SD = abs(I_p_wave_areas_3_right_septum_dev_avg_pls_2SD);
for i = 1:length(I_p_wave_areas_3_right_septum_sorted)
    if I_p_wave_areas_3_right_septum_dev_avg_pls_2SD(i) == min(I_p_wave_areas_3_right_septum_dev_avg_pls_2SD)
        I_p_wave_areas_3_right_septum_dev_avg_pls_2SD_index = i;
    end
end
I_p_wave_areas_3_right_septum_dev_avg_min_3SD = area_3(length(I_p_wave_areas_3_right_septum_sorted),1);
for i = 1:length(I_p_wave_areas_3_right_septum_sorted)
    I_p_wave_areas_3_right_septum_dev_avg_min_3SD(i) = I_p_wave_areas_3_right_septum_sorted(i) - (I_p_wave_areas_3_right_septum_average-3*I_p_wave_areas_3_right_septum_SD);
end
I_p_wave_areas_3_right_septum_dev_avg_min_3SD = abs(I_p_wave_areas_3_right_septum_dev_avg_min_3SD);
for i = 1:length(I_p_wave_areas_3_right_septum_sorted)
    if I_p_wave_areas_3_right_septum_dev_avg_min_3SD(i) == min(I_p_wave_areas_3_right_septum_dev_avg_min_3SD)
        I_p_wave_areas_3_right_septum_dev_avg_min_3SD_index = i;
    end
end
I_p_wave_areas_3_right_septum_dev_avg_pls_3SD = area_3(length(I_p_wave_areas_3_right_septum_sorted),1);
for i = 1:length(I_p_wave_areas_3_right_septum_sorted)
    I_p_wave_areas_3_right_septum_dev_avg_pls_3SD(i) = I_p_wave_areas_3_right_septum_sorted(i) - (I_p_wave_areas_3_right_septum_average+3*I_p_wave_areas_3_right_septum_SD);
end
I_p_wave_areas_3_right_septum_dev_avg_pls_3SD = abs(I_p_wave_areas_3_right_septum_dev_avg_pls_3SD);
for i = 1:length(I_p_wave_areas_3_right_septum_sorted)
    if I_p_wave_areas_3_right_septum_dev_avg_pls_3SD(i) == min(I_p_wave_areas_3_right_septum_dev_avg_pls_3SD)
        I_p_wave_areas_3_right_septum_dev_avg_pls_3SD_index = i;
    end
end
I_p_wave_areas_3_right_septum_tot_number_of_indices = length(I_p_wave_areas_3_right_septum_sorted);
%This is the approximation of the tot "area under curve." 
I_p_wave_areas_3_right_septum_quotient_within_1SD = (I_p_wave_areas_3_right_septum_dev_avg_pls_1SD_index - I_p_wave_areas_3_right_septumdev_avg_min_1SD_index)/I_p_wave_areas_3_right_septum_tot_number_of_indices;
I_p_wave_areas_3_right_septum_quotient_within_2SD = (I_p_wave_areas_3_right_septum_dev_avg_pls_2SD_index - I_p_wave_areas_3_dev_avg_min_2SD_index)/I_p_wave_areas_3_right_septum_tot_number_of_indices;
I_p_wave_areas_3_right_septum_quotient_within_3SD = (I_p_wave_areas_3_right_septum_dev_avg_pls_3SD_index - I_p_wave_areas_3_right_septum_dev_avg_min_3SD_index)/I_p_wave_areas_3_right_septum_tot_number_of_indices;
if (I_p_wave_areas_3_right_septum_quotient_within_1SD > 0.66 && I_p_wave_areas_3_right_septum_quotient_within_1SD < 0.70) && (I_p_wave_areas_3_right_septum_quotient_within_2SD > 0.93 && I_p_wave_areas_3_right_septum_quotient_within_2SD < 0.97) && (I_p_wave_areas_3_right_septum_quotient_within_3SD > 0.98 && I_p_wave_areas_3_right_septum_quotient_within_3SD < 1)
    isNormal(9,7) = 1;
end
if I_p_wave_areas_3_right_septum_quotient_within_1SD == 0
    isSpike(9,7) = 1;
end
I_p_wave_areas_3_right_atrial_appendage_sorted = sort(I_p_wave_areas_3_right_atrial_appendage);
I_p_wave_areas_3_right_atrial_appendage_dev_avg_min_1SD = area_3(length(I_p_wave_areas_3_right_atrial_appendage_sorted),1);
for i = 1:length(I_p_wave_areas_3_right_atrial_appendage_sorted)
    I_p_wave_areas_3_right_atrial_appendage_dev_avg_min_1SD(i) = I_p_wave_areas_3_right_atrial_appendage_sorted(i) - (I_p_wave_areas_3_right_atrial_appendage_average-I_p_wave_areas_3_right_atrial_appendage_SD);
end
I_p_wave_areas_3_right_atrial_appendage_dev_avg_min_1SD = abs(I_p_wave_areas_3_right_atrial_appendage_dev_avg_min_1SD);
for i = 1:length(I_p_wave_areas_3_right_atrial_appendage_sorted)
    if I_p_wave_areas_3_right_atrial_appendage_dev_avg_min_1SD(i) == min(I_p_wave_areas_3_right_atrial_appendage_dev_avg_min_1SD)
        I_p_wave_areas_3_right_atrial_appendage_dev_avg_min_1SD_index = i;
    end
end
I_p_wave_areas_3_right_atrial_appendage_dev_avg_pls_1SD = area_3(length(I_p_wave_areas_3_right_atrial_appendage_sorted),1);
for i = 1:length(I_p_wave_areas_3_right_atrial_appendage_sorted)
    I_p_wave_areas_3_right_atrial_appendage_dev_avg_pls_1SD(i) = I_p_wave_areas_3_right_atrial_appendage_sorted(i) - (I_p_wave_areas_3_right_atrial_appendage_average+I_p_wave_areas_3_right_atrial_appendage_SD);
end
I_p_wave_areas_3_right_atrial_appendage_dev_avg_pls_1SD = abs(I_p_wave_areas_3_right_atrial_appendage_dev_avg_pls_1SD);
for i = 1:length(I_p_wave_areas_3_right_atrial_appendage_sorted)
    if I_p_wave_areas_3_right_atrial_appendage_dev_avg_pls_1SD(i) == min(I_p_wave_areas_3_right_atrial_appendage_dev_avg_pls_1SD)
        I_p_wave_areas_3_right_atrial_appendage_dev_avg_pls_1SD_index = i;
    end
end
I_p_wave_areas_3_right_atrial_appendage_dev_avg_min_2SD = area_3(length(I_p_wave_areas_3_right_atrial_appendage_sorted),1);
for i = 1:length(I_p_wave_areas_3_right_atrial_appendage_sorted)
    I_p_wave_areas_3_right_atrial_appendage_dev_avg_min_2SD(i) = I_p_wave_areas_3_right_atrial_appendage_sorted(i) - (I_p_wave_areas_3_right_atrial_appendage_average-2*I_p_wave_areas_3_right_atrial_appendage_SD);
end
I_p_wave_areas_3_right_atrial_appendage_dev_avg_min_2SD = abs(I_p_wave_areas_3_right_atrial_appendage_dev_avg_min_2SD);
for i = 1:length(I_p_wave_areas_3_right_atrial_appendage_sorted)
    if I_p_wave_areas_3_right_atrial_appendage_dev_avg_min_2SD(i) == min(I_p_wave_areas_3_right_atrial_appendage_dev_avg_min_2SD)
        I_p_wave_areas_3_dev_avg_min_2SD_index = i;
    end
end
I_p_wave_areas_3_right_atrial_appendage_dev_avg_pls_2SD = area_3(length(I_p_wave_areas_3_right_atrial_appendage_sorted),1);
for i = 1:length(I_p_wave_areas_3_right_atrial_appendage_sorted)
    I_p_wave_areas_3_right_atrial_appendage_dev_avg_pls_2SD(i) = I_p_wave_areas_3_right_atrial_appendage_sorted(i) - (I_p_wave_areas_3_right_atrial_appendage_average+2*I_p_wave_areas_3_right_atrial_appendage_SD);
end
I_p_wave_areas_3_right_atrial_appendage_dev_avg_pls_2SD = abs(I_p_wave_areas_3_right_atrial_appendage_dev_avg_pls_2SD);
for i = 1:length(I_p_wave_areas_3_right_atrial_appendage_sorted)
    if I_p_wave_areas_3_right_atrial_appendage_dev_avg_pls_2SD(i) == min(I_p_wave_areas_3_right_atrial_appendage_dev_avg_pls_2SD)
        I_p_wave_areas_3_right_atrial_appendage_dev_avg_pls_2SD_index = i;
    end
end
I_p_wave_areas_3_right_atrial_appendage_dev_avg_min_3SD = area_3(length(I_p_wave_areas_3_right_atrial_appendage_sorted),1);
for i = 1:length(I_p_wave_areas_3_right_atrial_appendage_sorted)
    I_p_wave_areas_3_right_atrial_appendage_dev_avg_min_3SD(i) = I_p_wave_areas_3_right_atrial_appendage_sorted(i) - (I_p_wave_areas_3_right_atrial_appendage_average-3*I_p_wave_areas_3_right_atrial_appendage_SD);
end
I_p_wave_areas_3_right_atrial_appendage_dev_avg_min_3SD = abs(I_p_wave_areas_3_right_atrial_appendage_dev_avg_min_3SD);
for i = 1:length(I_p_wave_areas_3_right_atrial_appendage_sorted)
    if I_p_wave_areas_3_right_atrial_appendage_dev_avg_min_3SD(i) == min(I_p_wave_areas_3_right_atrial_appendage_dev_avg_min_3SD)
        I_p_wave_areas_3_right_atrial_appendage_dev_avg_min_3SD_index = i;
    end
end
I_p_wave_areas_3_right_atrial_appendage_dev_avg_pls_3SD = area_3(length(I_p_wave_areas_3_right_atrial_appendage_sorted),1);
for i = 1:length(I_p_wave_areas_3_right_atrial_appendage_sorted)
    I_p_wave_areas_3_right_atrial_appendage_dev_avg_pls_3SD(i) = I_p_wave_areas_3_right_atrial_appendage_sorted(i) - (I_p_wave_areas_3_right_atrial_appendage_average+3*I_p_wave_areas_3_right_atrial_appendage_SD);
end
I_p_wave_areas_3_right_atrial_appendage_dev_avg_pls_3SD = abs(I_p_wave_areas_3_right_atrial_appendage_dev_avg_pls_3SD);
for i = 1:length(I_p_wave_areas_3_right_atrial_appendage_sorted)
    if I_p_wave_areas_3_right_atrial_appendage_dev_avg_pls_3SD(i) == min(I_p_wave_areas_3_right_atrial_appendage_dev_avg_pls_3SD)
        I_p_wave_areas_3_right_atrial_appendage_dev_avg_pls_3SD_index = i;
    end
end
I_p_wave_areas_3_right_atrial_appendage_tot_number_of_indices = length(I_p_wave_areas_3_right_atrial_appendage_sorted);
%This is the approximation of the tot "area under curve." 
I_p_wave_areas_3_right_atrial_appendage_quotient_within_1SD = (I_p_wave_areas_3_right_atrial_appendage_dev_avg_pls_1SD_index - I_p_wave_areas_3_right_atrial_appendagedev_avg_min_1SD_index)/I_p_wave_areas_3_right_atrial_appendage_tot_number_of_indices;
I_p_wave_areas_3_right_atrial_appendage_quotient_within_2SD = (I_p_wave_areas_3_right_atrial_appendage_dev_avg_pls_2SD_index - I_p_wave_areas_3_dev_avg_min_2SD_index)/I_p_wave_areas_3_right_atrial_appendage_tot_number_of_indices;
I_p_wave_areas_3_right_atrial_appendage_quotient_within_3SD = (I_p_wave_areas_3_right_atrial_appendage_dev_avg_pls_3SD_index - I_p_wave_areas_3_right_atrial_appendage_dev_avg_min_3SD_index)/I_p_wave_areas_3_right_atrial_appendage_tot_number_of_indices;
if (I_p_wave_areas_3_right_atrial_appendage_quotient_within_1SD > 0.66 && I_p_wave_areas_3_right_atrial_appendage_quotient_within_1SD < 0.70) && (I_p_wave_areas_3_right_atrial_appendage_quotient_within_2SD > 0.93 && I_p_wave_areas_3_right_atrial_appendage_quotient_within_2SD < 0.97) && (I_p_wave_areas_3_right_atrial_appendage_quotient_within_3SD > 0.98 && I_p_wave_areas_3_right_atrial_appendage_quotient_within_3SD < 1)
    isNormal(9,8) = 1;
end
if I_p_wave_areas_3_right_atrial_appendage_quotient_within_1SD == 0
    isSpike(9,8) = 1;
end
I_p_wave_areas_3_pulmonary_veins_sorted = sort(I_p_wave_areas_3_pulmonary_veins);
I_p_wave_areas_3_pulmonary_veins_dev_avg_min_1SD = area_3(length(I_p_wave_areas_3_pulmonary_veins_sorted),1);
for i = 1:length(I_p_wave_areas_3_pulmonary_veins_sorted)
    I_p_wave_areas_3_pulmonary_veins_dev_avg_min_1SD(i) = I_p_wave_areas_3_pulmonary_veins_sorted(i) - (I_p_wave_areas_3_pulmonary_veins_average-I_p_wave_areas_3_pulmonary_veins_SD);
end
I_p_wave_areas_3_pulmonary_veins_dev_avg_min_1SD = abs(I_p_wave_areas_3_pulmonary_veins_dev_avg_min_1SD);
for i = 1:length(I_p_wave_areas_3_pulmonary_veins_sorted)
    if I_p_wave_areas_3_pulmonary_veins_dev_avg_min_1SD(i) == min(I_p_wave_areas_3_pulmonary_veins_dev_avg_min_1SD)
        I_p_wave_areas_3_pulmonary_veins_dev_avg_min_1SD_index = i;
    end
end
I_p_wave_areas_3_pulmonary_veins_dev_avg_pls_1SD = area_3(length(I_p_wave_areas_3_pulmonary_veins_sorted),1);
for i = 1:length(I_p_wave_areas_3_pulmonary_veins_sorted)
    I_p_wave_areas_3_pulmonary_veins_dev_avg_pls_1SD(i) = I_p_wave_areas_3_pulmonary_veins_sorted(i) - (I_p_wave_areas_3_pulmonary_veins_average+I_p_wave_areas_3_pulmonary_veins_SD);
end
I_p_wave_areas_3_pulmonary_veins_dev_avg_pls_1SD = abs(I_p_wave_areas_3_pulmonary_veins_dev_avg_pls_1SD);
for i = 1:length(I_p_wave_areas_3_pulmonary_veins_sorted)
    if I_p_wave_areas_3_pulmonary_veins_dev_avg_pls_1SD(i) == min(I_p_wave_areas_3_pulmonary_veins_dev_avg_pls_1SD)
        I_p_wave_areas_3_pulmonary_veins_dev_avg_pls_1SD_index = i;
    end
end
I_p_wave_areas_3_pulmonary_veins_dev_avg_min_2SD = area_3(length(I_p_wave_areas_3_pulmonary_veins_sorted),1);
for i = 1:length(I_p_wave_areas_3_pulmonary_veins_sorted)
    I_p_wave_areas_3_pulmonary_veins_dev_avg_min_2SD(i) = I_p_wave_areas_3_pulmonary_veins_sorted(i) - (I_p_wave_areas_3_pulmonary_veins_average-2*I_p_wave_areas_3_pulmonary_veins_SD);
end
I_p_wave_areas_3_pulmonary_veins_dev_avg_min_2SD = abs(I_p_wave_areas_3_pulmonary_veins_dev_avg_min_2SD);
for i = 1:length(I_p_wave_areas_3_pulmonary_veins_sorted)
    if I_p_wave_areas_3_pulmonary_veins_dev_avg_min_2SD(i) == min(I_p_wave_areas_3_pulmonary_veins_dev_avg_min_2SD)
        I_p_wave_areas_3_dev_avg_min_2SD_index = i;
    end
end
I_p_wave_areas_3_pulmonary_veins_dev_avg_pls_2SD = area_3(length(I_p_wave_areas_3_pulmonary_veins_sorted),1);
for i = 1:length(I_p_wave_areas_3_pulmonary_veins_sorted)
    I_p_wave_areas_3_pulmonary_veins_dev_avg_pls_2SD(i) = I_p_wave_areas_3_pulmonary_veins_sorted(i) - (I_p_wave_areas_3_pulmonary_veins_average+2*I_p_wave_areas_3_pulmonary_veins_SD);
end
I_p_wave_areas_3_pulmonary_veins_dev_avg_pls_2SD = abs(I_p_wave_areas_3_pulmonary_veins_dev_avg_pls_2SD);
for i = 1:length(I_p_wave_areas_3_pulmonary_veins_sorted)
    if I_p_wave_areas_3_pulmonary_veins_dev_avg_pls_2SD(i) == min(I_p_wave_areas_3_pulmonary_veins_dev_avg_pls_2SD)
        I_p_wave_areas_3_pulmonary_veins_dev_avg_pls_2SD_index = i;
    end
end
I_p_wave_areas_3_pulmonary_veins_dev_avg_min_3SD = area_3(length(I_p_wave_areas_3_pulmonary_veins_sorted),1);
for i = 1:length(I_p_wave_areas_3_pulmonary_veins_sorted)
    I_p_wave_areas_3_pulmonary_veins_dev_avg_min_3SD(i) = I_p_wave_areas_3_pulmonary_veins_sorted(i) - (I_p_wave_areas_3_pulmonary_veins_average-3*I_p_wave_areas_3_pulmonary_veins_SD);
end
I_p_wave_areas_3_pulmonary_veins_dev_avg_min_3SD = abs(I_p_wave_areas_3_pulmonary_veins_dev_avg_min_3SD);
for i = 1:length(I_p_wave_areas_3_pulmonary_veins_sorted)
    if I_p_wave_areas_3_pulmonary_veins_dev_avg_min_3SD(i) == min(I_p_wave_areas_3_pulmonary_veins_dev_avg_min_3SD)
        I_p_wave_areas_3_pulmonary_veins_dev_avg_min_3SD_index = i;
    end
end
I_p_wave_areas_3_pulmonary_veins_dev_avg_pls_3SD = area_3(length(I_p_wave_areas_3_pulmonary_veins_sorted),1);
for i = 1:length(I_p_wave_areas_3_pulmonary_veins_sorted)
    I_p_wave_areas_3_pulmonary_veins_dev_avg_pls_3SD(i) = I_p_wave_areas_3_pulmonary_veins_sorted(i) - (I_p_wave_areas_3_pulmonary_veins_average+3*I_p_wave_areas_3_pulmonary_veins_SD);
end
I_p_wave_areas_3_pulmonary_veins_dev_avg_pls_3SD = abs(I_p_wave_areas_3_pulmonary_veins_dev_avg_pls_3SD);
for i = 1:length(I_p_wave_areas_3_pulmonary_veins_sorted)
    if I_p_wave_areas_3_pulmonary_veins_dev_avg_pls_3SD(i) == min(I_p_wave_areas_3_pulmonary_veins_dev_avg_pls_3SD)
        I_p_wave_areas_3_pulmonary_veins_dev_avg_pls_3SD_index = i;
    end
end
I_p_wave_areas_3_pulmonary_veins_tot_number_of_indices = length(I_p_wave_areas_3_pulmonary_veins_sorted);
%This is the approximation of the tot "area under curve." 
I_p_wave_areas_3_pulmonary_veins_quotient_within_1SD = (I_p_wave_areas_3_pulmonary_veins_dev_avg_pls_1SD_index - I_p_wave_areas_3_pulmonary_veinsdev_avg_min_1SD_index)/I_p_wave_areas_3_pulmonary_veins_tot_number_of_indices;
I_p_wave_areas_3_pulmonary_veins_quotient_within_2SD = (I_p_wave_areas_3_pulmonary_veins_dev_avg_pls_2SD_index - I_p_wave_areas_3_dev_avg_min_2SD_index)/I_p_wave_areas_3_pulmonary_veins_tot_number_of_indices;
I_p_wave_areas_3_pulmonary_veins_quotient_within_3SD = (I_p_wave_areas_3_pulmonary_veins_dev_avg_pls_3SD_index - I_p_wave_areas_3_pulmonary_veins_dev_avg_min_3SD_index)/I_p_wave_areas_3_pulmonary_veins_tot_number_of_indices;
if (I_p_wave_areas_3_pulmonary_veins_quotient_within_1SD > 0.66 && I_p_wave_areas_3_pulmonary_veins_quotient_within_1SD < 0.70) && (I_p_wave_areas_3_pulmonary_veins_quotient_within_2SD > 0.93 && I_p_wave_areas_3_pulmonary_veins_quotient_within_2SD < 0.97) && (I_p_wave_areas_3_pulmonary_veins_quotient_within_3SD > 0.98 && I_p_wave_areas_3_pulmonary_veins_quotient_within_3SD < 1)
    isNormal(9,9) = 1;
end
if I_p_wave_areas_3_pulmonary_veins_quotient_within_1SD == 0
    isSpike(9,9) = 1;
end
I_p_wave_areas_3_mitral_annulus_sorted = sort(I_p_wave_areas_3_mitral_annulus);
I_p_wave_areas_3_mitral_annulus_dev_avg_min_1SD = area_3(length(I_p_wave_areas_3_mitral_annulus_sorted),1);
for i = 1:length(I_p_wave_areas_3_mitral_annulus_sorted)
    I_p_wave_areas_3_mitral_annulus_dev_avg_min_1SD(i) = I_p_wave_areas_3_mitral_annulus_sorted(i) - (I_p_wave_areas_3_mitral_annulus_average-I_p_wave_areas_3_mitral_annulus_SD);
end
I_p_wave_areas_3_mitral_annulus_dev_avg_min_1SD = abs(I_p_wave_areas_3_mitral_annulus_dev_avg_min_1SD);
for i = 1:length(I_p_wave_areas_3_mitral_annulus_sorted)
    if I_p_wave_areas_3_mitral_annulus_dev_avg_min_1SD(i) == min(I_p_wave_areas_3_mitral_annulus_dev_avg_min_1SD)
        I_p_wave_areas_3_mitral_annulus_dev_avg_min_1SD_index = i;
    end
end
I_p_wave_areas_3_mitral_annulus_dev_avg_pls_1SD = area_3(length(I_p_wave_areas_3_mitral_annulus_sorted),1);
for i = 1:length(I_p_wave_areas_3_mitral_annulus_sorted)
    I_p_wave_areas_3_mitral_annulus_dev_avg_pls_1SD(i) = I_p_wave_areas_3_mitral_annulus_sorted(i) - (I_p_wave_areas_3_mitral_annulus_average+I_p_wave_areas_3_mitral_annulus_SD);
end
I_p_wave_areas_3_mitral_annulus_dev_avg_pls_1SD = abs(I_p_wave_areas_3_mitral_annulus_dev_avg_pls_1SD);
for i = 1:length(I_p_wave_areas_3_mitral_annulus_sorted)
    if I_p_wave_areas_3_mitral_annulus_dev_avg_pls_1SD(i) == min(I_p_wave_areas_3_mitral_annulus_dev_avg_pls_1SD)
        I_p_wave_areas_3_mitral_annulus_dev_avg_pls_1SD_index = i;
    end
end
I_p_wave_areas_3_mitral_annulus_dev_avg_min_2SD = area_3(length(I_p_wave_areas_3_mitral_annulus_sorted),1);
for i = 1:length(I_p_wave_areas_3_mitral_annulus_sorted)
    I_p_wave_areas_3_mitral_annulus_dev_avg_min_2SD(i) = I_p_wave_areas_3_mitral_annulus_sorted(i) - (I_p_wave_areas_3_mitral_annulus_average-2*I_p_wave_areas_3_mitral_annulus_SD);
end
I_p_wave_areas_3_mitral_annulus_dev_avg_min_2SD = abs(I_p_wave_areas_3_mitral_annulus_dev_avg_min_2SD);
for i = 1:length(I_p_wave_areas_3_mitral_annulus_sorted)
    if I_p_wave_areas_3_mitral_annulus_dev_avg_min_2SD(i) == min(I_p_wave_areas_3_mitral_annulus_dev_avg_min_2SD)
        I_p_wave_areas_3_dev_avg_min_2SD_index = i;
    end
end
I_p_wave_areas_3_mitral_annulus_dev_avg_pls_2SD = area_3(length(I_p_wave_areas_3_mitral_annulus_sorted),1);
for i = 1:length(I_p_wave_areas_3_mitral_annulus_sorted)
    I_p_wave_areas_3_mitral_annulus_dev_avg_pls_2SD(i) = I_p_wave_areas_3_mitral_annulus_sorted(i) - (I_p_wave_areas_3_mitral_annulus_average+2*I_p_wave_areas_3_mitral_annulus_SD);
end
I_p_wave_areas_3_mitral_annulus_dev_avg_pls_2SD = abs(I_p_wave_areas_3_mitral_annulus_dev_avg_pls_2SD);
for i = 1:length(I_p_wave_areas_3_mitral_annulus_sorted)
    if I_p_wave_areas_3_mitral_annulus_dev_avg_pls_2SD(i) == min(I_p_wave_areas_3_mitral_annulus_dev_avg_pls_2SD)
        I_p_wave_areas_3_mitral_annulus_dev_avg_pls_2SD_index = i;
    end
end
I_p_wave_areas_3_mitral_annulus_dev_avg_min_3SD = area_3(length(I_p_wave_areas_3_mitral_annulus_sorted),1);
for i = 1:length(I_p_wave_areas_3_mitral_annulus_sorted)
    I_p_wave_areas_3_mitral_annulus_dev_avg_min_3SD(i) = I_p_wave_areas_3_mitral_annulus_sorted(i) - (I_p_wave_areas_3_mitral_annulus_average-3*I_p_wave_areas_3_mitral_annulus_SD);
end
I_p_wave_areas_3_mitral_annulus_dev_avg_min_3SD = abs(I_p_wave_areas_3_mitral_annulus_dev_avg_min_3SD);
for i = 1:length(I_p_wave_areas_3_mitral_annulus_sorted)
    if I_p_wave_areas_3_mitral_annulus_dev_avg_min_3SD(i) == min(I_p_wave_areas_3_mitral_annulus_dev_avg_min_3SD)
        I_p_wave_areas_3_mitral_annulus_dev_avg_min_3SD_index = i;
    end
end
I_p_wave_areas_3_mitral_annulus_dev_avg_pls_3SD = area_3(length(I_p_wave_areas_3_mitral_annulus_sorted),1);
for i = 1:length(I_p_wave_areas_3_mitral_annulus_sorted)
    I_p_wave_areas_3_mitral_annulus_dev_avg_pls_3SD(i) = I_p_wave_areas_3_mitral_annulus_sorted(i) - (I_p_wave_areas_3_mitral_annulus_average+3*I_p_wave_areas_3_mitral_annulus_SD);
end
I_p_wave_areas_3_mitral_annulus_dev_avg_pls_3SD = abs(I_p_wave_areas_3_mitral_annulus_dev_avg_pls_3SD);
for i = 1:length(I_p_wave_areas_3_mitral_annulus_sorted)
    if I_p_wave_areas_3_mitral_annulus_dev_avg_pls_3SD(i) == min(I_p_wave_areas_3_mitral_annulus_dev_avg_pls_3SD)
        I_p_wave_areas_3_mitral_annulus_dev_avg_pls_3SD_index = i;
    end
end
I_p_wave_areas_3_mitral_annulus_tot_number_of_indices = length(I_p_wave_areas_3_mitral_annulus_sorted);
%This is the approximation of the tot "area under curve." 
I_p_wave_areas_3_mitral_annulus_quotient_within_1SD = (I_p_wave_areas_3_mitral_annulus_dev_avg_pls_1SD_index - I_p_wave_areas_3_mitral_annulusdev_avg_min_1SD_index)/I_p_wave_areas_3_mitral_annulus_tot_number_of_indices;
I_p_wave_areas_3_mitral_annulus_quotient_within_2SD = (I_p_wave_areas_3_mitral_annulus_dev_avg_pls_2SD_index - I_p_wave_areas_3_dev_avg_min_2SD_index)/I_p_wave_areas_3_mitral_annulus_tot_number_of_indices;
I_p_wave_areas_3_mitral_annulus_quotient_within_3SD = (I_p_wave_areas_3_mitral_annulus_dev_avg_pls_3SD_index - I_p_wave_areas_3_mitral_annulus_dev_avg_min_3SD_index)/I_p_wave_areas_3_mitral_annulus_tot_number_of_indices;
if (I_p_wave_areas_3_mitral_annulus_quotient_within_1SD > 0.66 && I_p_wave_areas_3_mitral_annulus_quotient_within_1SD < 0.70) && (I_p_wave_areas_3_mitral_annulus_quotient_within_2SD > 0.93 && I_p_wave_areas_3_mitral_annulus_quotient_within_2SD < 0.97) && (I_p_wave_areas_3_mitral_annulus_quotient_within_3SD > 0.98 && I_p_wave_areas_3_mitral_annulus_quotient_within_3SD < 1)
    isNormal(9,10) = 1;
end
if I_p_wave_areas_3_mitral_annulus_quotient_within_1SD == 0
    isSpike(9,10) = 1;
end
I_p_wave_areas_3_CS_body_sorted = sort(I_p_wave_areas_3_CS_body);
I_p_wave_areas_3_CS_body_dev_avg_min_1SD = area_3(length(I_p_wave_areas_3_CS_body_sorted),1);
for i = 1:length(I_p_wave_areas_3_CS_body_sorted)
    I_p_wave_areas_3_CS_body_dev_avg_min_1SD(i) = I_p_wave_areas_3_CS_body_sorted(i) - (I_p_wave_areas_3_CS_body_average-I_p_wave_areas_3_CS_body_SD);
end
I_p_wave_areas_3_CS_body_dev_avg_min_1SD = abs(I_p_wave_areas_3_CS_body_dev_avg_min_1SD);
for i = 1:length(I_p_wave_areas_3_CS_body_sorted)
    if I_p_wave_areas_3_CS_body_dev_avg_min_1SD(i) == min(I_p_wave_areas_3_CS_body_dev_avg_min_1SD)
        I_p_wave_areas_3_CS_body_dev_avg_min_1SD_index = i;
    end
end
I_p_wave_areas_3_CS_body_dev_avg_pls_1SD = area_3(length(I_p_wave_areas_3_CS_body_sorted),1);
for i = 1:length(I_p_wave_areas_3_CS_body_sorted)
    I_p_wave_areas_3_CS_body_dev_avg_pls_1SD(i) = I_p_wave_areas_3_CS_body_sorted(i) - (I_p_wave_areas_3_CS_body_average+I_p_wave_areas_3_CS_body_SD);
end
I_p_wave_areas_3_CS_body_dev_avg_pls_1SD = abs(I_p_wave_areas_3_CS_body_dev_avg_pls_1SD);
for i = 1:length(I_p_wave_areas_3_CS_body_sorted)
    if I_p_wave_areas_3_CS_body_dev_avg_pls_1SD(i) == min(I_p_wave_areas_3_CS_body_dev_avg_pls_1SD)
        I_p_wave_areas_3_CS_body_dev_avg_pls_1SD_index = i;
    end
end
I_p_wave_areas_3_CS_body_dev_avg_min_2SD = area_3(length(I_p_wave_areas_3_CS_body_sorted),1);
for i = 1:length(I_p_wave_areas_3_CS_body_sorted)
    I_p_wave_areas_3_CS_body_dev_avg_min_2SD(i) = I_p_wave_areas_3_CS_body_sorted(i) - (I_p_wave_areas_3_CS_body_average-2*I_p_wave_areas_3_CS_body_SD);
end
I_p_wave_areas_3_CS_body_dev_avg_min_2SD = abs(I_p_wave_areas_3_CS_body_dev_avg_min_2SD);
for i = 1:length(I_p_wave_areas_3_CS_body_sorted)
    if I_p_wave_areas_3_CS_body_dev_avg_min_2SD(i) == min(I_p_wave_areas_3_CS_body_dev_avg_min_2SD)
        I_p_wave_areas_3_dev_avg_min_2SD_index = i;
    end
end
I_p_wave_areas_3_CS_body_dev_avg_pls_2SD = area_3(length(I_p_wave_areas_3_CS_body_sorted),1);
for i = 1:length(I_p_wave_areas_3_CS_body_sorted)
    I_p_wave_areas_3_CS_body_dev_avg_pls_2SD(i) = I_p_wave_areas_3_CS_body_sorted(i) - (I_p_wave_areas_3_CS_body_average+2*I_p_wave_areas_3_CS_body_SD);
end
I_p_wave_areas_3_CS_body_dev_avg_pls_2SD = abs(I_p_wave_areas_3_CS_body_dev_avg_pls_2SD);
for i = 1:length(I_p_wave_areas_3_CS_body_sorted)
    if I_p_wave_areas_3_CS_body_dev_avg_pls_2SD(i) == min(I_p_wave_areas_3_CS_body_dev_avg_pls_2SD)
        I_p_wave_areas_3_CS_body_dev_avg_pls_2SD_index = i;
    end
end
I_p_wave_areas_3_CS_body_dev_avg_min_3SD = area_3(length(I_p_wave_areas_3_CS_body_sorted),1);
for i = 1:length(I_p_wave_areas_3_CS_body_sorted)
    I_p_wave_areas_3_CS_body_dev_avg_min_3SD(i) = I_p_wave_areas_3_CS_body_sorted(i) - (I_p_wave_areas_3_CS_body_average-3*I_p_wave_areas_3_CS_body_SD);
end
I_p_wave_areas_3_CS_body_dev_avg_min_3SD = abs(I_p_wave_areas_3_CS_body_dev_avg_min_3SD);
for i = 1:length(I_p_wave_areas_3_CS_body_sorted)
    if I_p_wave_areas_3_CS_body_dev_avg_min_3SD(i) == min(I_p_wave_areas_3_CS_body_dev_avg_min_3SD)
        I_p_wave_areas_3_CS_body_dev_avg_min_3SD_index = i;
    end
end
I_p_wave_areas_3_CS_body_dev_avg_pls_3SD = area_3(length(I_p_wave_areas_3_CS_body_sorted),1);
for i = 1:length(I_p_wave_areas_3_CS_body_sorted)
    I_p_wave_areas_3_CS_body_dev_avg_pls_3SD(i) = I_p_wave_areas_3_CS_body_sorted(i) - (I_p_wave_areas_3_CS_body_average+3*I_p_wave_areas_3_CS_body_SD);
end
I_p_wave_areas_3_CS_body_dev_avg_pls_3SD = abs(I_p_wave_areas_3_CS_body_dev_avg_pls_3SD);
for i = 1:length(I_p_wave_areas_3_CS_body_sorted)
    if I_p_wave_areas_3_CS_body_dev_avg_pls_3SD(i) == min(I_p_wave_areas_3_CS_body_dev_avg_pls_3SD)
        I_p_wave_areas_3_CS_body_dev_avg_pls_3SD_index = i;
    end
end
I_p_wave_areas_3_CS_body_tot_number_of_indices = length(I_p_wave_areas_3_CS_body_sorted);
%This is the approximation of the tot "area under curve." 
I_p_wave_areas_3_CS_body_quotient_within_1SD = (I_p_wave_areas_3_CS_body_dev_avg_pls_1SD_index - I_p_wave_areas_3_CS_bodydev_avg_min_1SD_index)/I_p_wave_areas_3_CS_body_tot_number_of_indices;
I_p_wave_areas_3_CS_body_quotient_within_2SD = (I_p_wave_areas_3_CS_body_dev_avg_pls_2SD_index - I_p_wave_areas_3_dev_avg_min_2SD_index)/I_p_wave_areas_3_CS_body_tot_number_of_indices;
I_p_wave_areas_3_CS_body_quotient_within_3SD = (I_p_wave_areas_3_CS_body_dev_avg_pls_3SD_index - I_p_wave_areas_3_CS_body_dev_avg_min_3SD_index)/I_p_wave_areas_3_CS_body_tot_number_of_indices;
if (I_p_wave_areas_3_CS_body_quotient_within_1SD > 0.66 && I_p_wave_areas_3_CS_body_quotient_within_1SD < 0.70) && (I_p_wave_areas_3_CS_body_quotient_within_2SD > 0.93 && I_p_wave_areas_3_CS_body_quotient_within_2SD < 0.97) && (I_p_wave_areas_3_CS_body_quotient_within_3SD > 0.98 && I_p_wave_areas_3_CS_body_quotient_within_3SD < 1)
    isNormal(9,11) = 1;
end
if I_p_wave_areas_3_CS_body_quotient_within_1SD == 0
    isSpike(9,11) = 1;
end
I_p_wave_areas_3_left_septum_sorted = sort(I_p_wave_areas_3_left_septum);
I_p_wave_areas_3_left_septum_dev_avg_min_1SD = area_3(length(I_p_wave_areas_3_left_septum_sorted),1);
for i = 1:length(I_p_wave_areas_3_left_septum_sorted)
    I_p_wave_areas_3_left_septum_dev_avg_min_1SD(i) = I_p_wave_areas_3_left_septum_sorted(i) - (I_p_wave_areas_3_left_septum_average-I_p_wave_areas_3_left_septum_SD);
end
I_p_wave_areas_3_left_septum_dev_avg_min_1SD = abs(I_p_wave_areas_3_left_septum_dev_avg_min_1SD);
for i = 1:length(I_p_wave_areas_3_left_septum_sorted)
    if I_p_wave_areas_3_left_septum_dev_avg_min_1SD(i) == min(I_p_wave_areas_3_left_septum_dev_avg_min_1SD)
        I_p_wave_areas_3_left_septum_dev_avg_min_1SD_index = i;
    end
end
I_p_wave_areas_3_left_septum_dev_avg_pls_1SD = area_3(length(I_p_wave_areas_3_left_septum_sorted),1);
for i = 1:length(I_p_wave_areas_3_left_septum_sorted)
    I_p_wave_areas_3_left_septum_dev_avg_pls_1SD(i) = I_p_wave_areas_3_left_septum_sorted(i) - (I_p_wave_areas_3_left_septum_average+I_p_wave_areas_3_left_septum_SD);
end
I_p_wave_areas_3_left_septum_dev_avg_pls_1SD = abs(I_p_wave_areas_3_left_septum_dev_avg_pls_1SD);
for i = 1:length(I_p_wave_areas_3_left_septum_sorted)
    if I_p_wave_areas_3_left_septum_dev_avg_pls_1SD(i) == min(I_p_wave_areas_3_left_septum_dev_avg_pls_1SD)
        I_p_wave_areas_3_left_septum_dev_avg_pls_1SD_index = i;
    end
end
I_p_wave_areas_3_left_septum_dev_avg_min_2SD = area_3(length(I_p_wave_areas_3_left_septum_sorted),1);
for i = 1:length(I_p_wave_areas_3_left_septum_sorted)
    I_p_wave_areas_3_left_septum_dev_avg_min_2SD(i) = I_p_wave_areas_3_left_septum_sorted(i) - (I_p_wave_areas_3_left_septum_average-2*I_p_wave_areas_3_left_septum_SD);
end
I_p_wave_areas_3_left_septum_dev_avg_min_2SD = abs(I_p_wave_areas_3_left_septum_dev_avg_min_2SD);
for i = 1:length(I_p_wave_areas_3_left_septum_sorted)
    if I_p_wave_areas_3_left_septum_dev_avg_min_2SD(i) == min(I_p_wave_areas_3_left_septum_dev_avg_min_2SD)
        I_p_wave_areas_3_dev_avg_min_2SD_index = i;
    end
end
I_p_wave_areas_3_left_septum_dev_avg_pls_2SD = area_3(length(I_p_wave_areas_3_left_septum_sorted),1);
for i = 1:length(I_p_wave_areas_3_left_septum_sorted)
    I_p_wave_areas_3_left_septum_dev_avg_pls_2SD(i) = I_p_wave_areas_3_left_septum_sorted(i) - (I_p_wave_areas_3_left_septum_average+2*I_p_wave_areas_3_left_septum_SD);
end
I_p_wave_areas_3_left_septum_dev_avg_pls_2SD = abs(I_p_wave_areas_3_left_septum_dev_avg_pls_2SD);
for i = 1:length(I_p_wave_areas_3_left_septum_sorted)
    if I_p_wave_areas_3_left_septum_dev_avg_pls_2SD(i) == min(I_p_wave_areas_3_left_septum_dev_avg_pls_2SD)
        I_p_wave_areas_3_left_septum_dev_avg_pls_2SD_index = i;
    end
end
I_p_wave_areas_3_left_septum_dev_avg_min_3SD = area_3(length(I_p_wave_areas_3_left_septum_sorted),1);
for i = 1:length(I_p_wave_areas_3_left_septum_sorted)
    I_p_wave_areas_3_left_septum_dev_avg_min_3SD(i) = I_p_wave_areas_3_left_septum_sorted(i) - (I_p_wave_areas_3_left_septum_average-3*I_p_wave_areas_3_left_septum_SD);
end
I_p_wave_areas_3_left_septum_dev_avg_min_3SD = abs(I_p_wave_areas_3_left_septum_dev_avg_min_3SD);
for i = 1:length(I_p_wave_areas_3_left_septum_sorted)
    if I_p_wave_areas_3_left_septum_dev_avg_min_3SD(i) == min(I_p_wave_areas_3_left_septum_dev_avg_min_3SD)
        I_p_wave_areas_3_left_septum_dev_avg_min_3SD_index = i;
    end
end
I_p_wave_areas_3_left_septum_dev_avg_pls_3SD = area_3(length(I_p_wave_areas_3_left_septum_sorted),1);
for i = 1:length(I_p_wave_areas_3_left_septum_sorted)
    I_p_wave_areas_3_left_septum_dev_avg_pls_3SD(i) = I_p_wave_areas_3_left_septum_sorted(i) - (I_p_wave_areas_3_left_septum_average+3*I_p_wave_areas_3_left_septum_SD);
end
I_p_wave_areas_3_left_septum_dev_avg_pls_3SD = abs(I_p_wave_areas_3_left_septum_dev_avg_pls_3SD);
for i = 1:length(I_p_wave_areas_3_left_septum_sorted)
    if I_p_wave_areas_3_left_septum_dev_avg_pls_3SD(i) == min(I_p_wave_areas_3_left_septum_dev_avg_pls_3SD)
        I_p_wave_areas_3_left_septum_dev_avg_pls_3SD_index = i;
    end
end
I_p_wave_areas_3_left_septum_tot_number_of_indices = length(I_p_wave_areas_3_left_septum_sorted);
%This is the approximation of the tot "area under curve." 
I_p_wave_areas_3_left_septum_quotient_within_1SD = (I_p_wave_areas_3_left_septum_dev_avg_pls_1SD_index - I_p_wave_areas_3_left_septumdev_avg_min_1SD_index)/I_p_wave_areas_3_left_septum_tot_number_of_indices;
I_p_wave_areas_3_left_septum_quotient_within_2SD = (I_p_wave_areas_3_left_septum_dev_avg_pls_2SD_index - I_p_wave_areas_3_dev_avg_min_2SD_index)/I_p_wave_areas_3_left_septum_tot_number_of_indices;
I_p_wave_areas_3_left_septum_quotient_within_3SD = (I_p_wave_areas_3_left_septum_dev_avg_pls_3SD_index - I_p_wave_areas_3_left_septum_dev_avg_min_3SD_index)/I_p_wave_areas_3_left_septum_tot_number_of_indices;
if (I_p_wave_areas_3_left_septum_quotient_within_1SD > 0.66 && I_p_wave_areas_3_left_septum_quotient_within_1SD < 0.70) && (I_p_wave_areas_3_left_septum_quotient_within_2SD > 0.93 && I_p_wave_areas_3_left_septum_quotient_within_2SD < 0.97) && (I_p_wave_areas_3_left_septum_quotient_within_3SD > 0.98 && I_p_wave_areas_3_left_septum_quotient_within_3SD < 1)
    isNormal(9,12) = 1;
end
if I_p_wave_areas_3_left_septum_quotient_within_1SD == 0
    isSpike(9,12) = 1;
end
I_p_wave_areas_3_left_atrial_appendage_sorted = sort(I_p_wave_areas_3_left_atrial_appendage);
I_p_wave_areas_3_left_atrial_appendage_dev_avg_min_1SD = area_3(length(I_p_wave_areas_3_left_atrial_appendage_sorted),1);
for i = 1:length(I_p_wave_areas_3_left_atrial_appendage_sorted)
    I_p_wave_areas_3_left_atrial_appendage_dev_avg_min_1SD(i) = I_p_wave_areas_3_left_atrial_appendage_sorted(i) - (I_p_wave_areas_3_left_atrial_appendage_average-I_p_wave_areas_3_left_atrial_appendage_SD);
end
I_p_wave_areas_3_left_atrial_appendage_dev_avg_min_1SD = abs(I_p_wave_areas_3_left_atrial_appendage_dev_avg_min_1SD);
for i = 1:length(I_p_wave_areas_3_left_atrial_appendage_sorted)
    if I_p_wave_areas_3_left_atrial_appendage_dev_avg_min_1SD(i) == min(I_p_wave_areas_3_left_atrial_appendage_dev_avg_min_1SD)
        I_p_wave_areas_3_left_atrial_appendage_dev_avg_min_1SD_index = i;
    end
end
I_p_wave_areas_3_left_atrial_appendage_dev_avg_pls_1SD = area_3(length(I_p_wave_areas_3_left_atrial_appendage_sorted),1);
for i = 1:length(I_p_wave_areas_3_left_atrial_appendage_sorted)
    I_p_wave_areas_3_left_atrial_appendage_dev_avg_pls_1SD(i) = I_p_wave_areas_3_left_atrial_appendage_sorted(i) - (I_p_wave_areas_3_left_atrial_appendage_average+I_p_wave_areas_3_left_atrial_appendage_SD);
end
I_p_wave_areas_3_left_atrial_appendage_dev_avg_pls_1SD = abs(I_p_wave_areas_3_left_atrial_appendage_dev_avg_pls_1SD);
for i = 1:length(I_p_wave_areas_3_left_atrial_appendage_sorted)
    if I_p_wave_areas_3_left_atrial_appendage_dev_avg_pls_1SD(i) == min(I_p_wave_areas_3_left_atrial_appendage_dev_avg_pls_1SD)
        I_p_wave_areas_3_left_atrial_appendage_dev_avg_pls_1SD_index = i;
    end
end
I_p_wave_areas_3_left_atrial_appendage_dev_avg_min_2SD = area_3(length(I_p_wave_areas_3_left_atrial_appendage_sorted),1);
for i = 1:length(I_p_wave_areas_3_left_atrial_appendage_sorted)
    I_p_wave_areas_3_left_atrial_appendage_dev_avg_min_2SD(i) = I_p_wave_areas_3_left_atrial_appendage_sorted(i) - (I_p_wave_areas_3_left_atrial_appendage_average-2*I_p_wave_areas_3_left_atrial_appendage_SD);
end
I_p_wave_areas_3_left_atrial_appendage_dev_avg_min_2SD = abs(I_p_wave_areas_3_left_atrial_appendage_dev_avg_min_2SD);
for i = 1:length(I_p_wave_areas_3_left_atrial_appendage_sorted)
    if I_p_wave_areas_3_left_atrial_appendage_dev_avg_min_2SD(i) == min(I_p_wave_areas_3_left_atrial_appendage_dev_avg_min_2SD)
        I_p_wave_areas_3_dev_avg_min_2SD_index = i;
    end
end
I_p_wave_areas_3_left_atrial_appendage_dev_avg_pls_2SD = area_3(length(I_p_wave_areas_3_left_atrial_appendage_sorted),1);
for i = 1:length(I_p_wave_areas_3_left_atrial_appendage_sorted)
    I_p_wave_areas_3_left_atrial_appendage_dev_avg_pls_2SD(i) = I_p_wave_areas_3_left_atrial_appendage_sorted(i) - (I_p_wave_areas_3_left_atrial_appendage_average+2*I_p_wave_areas_3_left_atrial_appendage_SD);
end
I_p_wave_areas_3_left_atrial_appendage_dev_avg_pls_2SD = abs(I_p_wave_areas_3_left_atrial_appendage_dev_avg_pls_2SD);
for i = 1:length(I_p_wave_areas_3_left_atrial_appendage_sorted)
    if I_p_wave_areas_3_left_atrial_appendage_dev_avg_pls_2SD(i) == min(I_p_wave_areas_3_left_atrial_appendage_dev_avg_pls_2SD)
        I_p_wave_areas_3_left_atrial_appendage_dev_avg_pls_2SD_index = i;
    end
end
I_p_wave_areas_3_left_atrial_appendage_dev_avg_min_3SD = area_3(length(I_p_wave_areas_3_left_atrial_appendage_sorted),1);
for i = 1:length(I_p_wave_areas_3_left_atrial_appendage_sorted)
    I_p_wave_areas_3_left_atrial_appendage_dev_avg_min_3SD(i) = I_p_wave_areas_3_left_atrial_appendage_sorted(i) - (I_p_wave_areas_3_left_atrial_appendage_average-3*I_p_wave_areas_3_left_atrial_appendage_SD);
end
I_p_wave_areas_3_left_atrial_appendage_dev_avg_min_3SD = abs(I_p_wave_areas_3_left_atrial_appendage_dev_avg_min_3SD);
for i = 1:length(I_p_wave_areas_3_left_atrial_appendage_sorted)
    if I_p_wave_areas_3_left_atrial_appendage_dev_avg_min_3SD(i) == min(I_p_wave_areas_3_left_atrial_appendage_dev_avg_min_3SD)
        I_p_wave_areas_3_left_atrial_appendage_dev_avg_min_3SD_index = i;
    end
end
I_p_wave_areas_3_left_atrial_appendage_dev_avg_pls_3SD = area_3(length(I_p_wave_areas_3_left_atrial_appendage_sorted),1);
for i = 1:length(I_p_wave_areas_3_left_atrial_appendage_sorted)
    I_p_wave_areas_3_left_atrial_appendage_dev_avg_pls_3SD(i) = I_p_wave_areas_3_left_atrial_appendage_sorted(i) - (I_p_wave_areas_3_left_atrial_appendage_average+3*I_p_wave_areas_3_left_atrial_appendage_SD);
end
I_p_wave_areas_3_left_atrial_appendage_dev_avg_pls_3SD = abs(I_p_wave_areas_3_left_atrial_appendage_dev_avg_pls_3SD);
for i = 1:length(I_p_wave_areas_3_left_atrial_appendage_sorted)
    if I_p_wave_areas_3_left_atrial_appendage_dev_avg_pls_3SD(i) == min(I_p_wave_areas_3_left_atrial_appendage_dev_avg_pls_3SD)
        I_p_wave_areas_3_left_atrial_appendage_dev_avg_pls_3SD_index = i;
    end
end
I_p_wave_areas_3_left_atrial_appendage_tot_number_of_indices = length(I_p_wave_areas_3_left_atrial_appendage_sorted);
%This is the approximation of the tot "area under curve." 
I_p_wave_areas_3_left_atrial_appendage_quotient_within_1SD = (I_p_wave_areas_3_left_atrial_appendage_dev_avg_pls_1SD_index - I_p_wave_areas_3_left_atrial_appendagedev_avg_min_1SD_index)/I_p_wave_areas_3_left_atrial_appendage_tot_number_of_indices;
I_p_wave_areas_3_left_atrial_appendage_quotient_within_2SD = (I_p_wave_areas_3_left_atrial_appendage_dev_avg_pls_2SD_index - I_p_wave_areas_3_dev_avg_min_2SD_index)/I_p_wave_areas_3_left_atrial_appendage_tot_number_of_indices;
I_p_wave_areas_3_left_atrial_appendage_quotient_within_3SD = (I_p_wave_areas_3_left_atrial_appendage_dev_avg_pls_3SD_index - I_p_wave_areas_3_left_atrial_appendage_dev_avg_min_3SD_index)/I_p_wave_areas_3_left_atrial_appendage_tot_number_of_indices;
if (I_p_wave_areas_3_left_atrial_appendage_quotient_within_1SD > 0.66 && I_p_wave_areas_3_left_atrial_appendage_quotient_within_1SD < 0.70) && (I_p_wave_areas_3_left_atrial_appendage_quotient_within_2SD > 0.93 && I_p_wave_areas_3_left_atrial_appendage_quotient_within_2SD < 0.97) && (I_p_wave_areas_3_left_atrial_appendage_quotient_within_3SD > 0.98 && I_p_wave_areas_3_left_atrial_appendage_quotient_within_3SD < 1)
    isNormal(9,13) = 1;
end
if I_p_wave_areas_3_left_atrial_appendage_quotient_within_1SD == 0
    isSpike(9,13) = 1;
end
I_p_wave_areas_4_SA_node_sorted = sort(I_p_wave_areas_4_SA_node);
I_p_wave_areas_4_SA_node_dev_avg_min_1SD = area_4(length(I_p_wave_areas_4_SA_node_sorted),1);
for i = 1:length(I_p_wave_areas_4_SA_node_sorted)
    I_p_wave_areas_4_SA_node_dev_avg_min_1SD(i) = I_p_wave_areas_4_SA_node_sorted(i) - (I_p_wave_areas_4_SA_node_average-I_p_wave_areas_4_SA_node_SD);
end
I_p_wave_areas_4_SA_node_dev_avg_min_1SD = abs(I_p_wave_areas_4_SA_node_dev_avg_min_1SD);
for i = 1:length(I_p_wave_areas_4_SA_node_sorted)
    if I_p_wave_areas_4_SA_node_dev_avg_min_1SD(i) == min(I_p_wave_areas_4_SA_node_dev_avg_min_1SD)
        I_p_wave_areas_4_SA_node_dev_avg_min_1SD_index = i;
    end
end
I_p_wave_areas_4_SA_node_dev_avg_pls_1SD = area_4(length(I_p_wave_areas_4_SA_node_sorted),1);
for i = 1:length(I_p_wave_areas_4_SA_node_sorted)
    I_p_wave_areas_4_SA_node_dev_avg_pls_1SD(i) = I_p_wave_areas_4_SA_node_sorted(i) - (I_p_wave_areas_4_SA_node_average+I_p_wave_areas_4_SA_node_SD);
end
I_p_wave_areas_4_SA_node_dev_avg_pls_1SD = abs(I_p_wave_areas_4_SA_node_dev_avg_pls_1SD);
for i = 1:length(I_p_wave_areas_4_SA_node_sorted)
    if I_p_wave_areas_4_SA_node_dev_avg_pls_1SD(i) == min(I_p_wave_areas_4_SA_node_dev_avg_pls_1SD)
        I_p_wave_areas_4_SA_node_dev_avg_pls_1SD_index = i;
    end
end
I_p_wave_areas_4_SA_node_dev_avg_min_2SD = area_4(length(I_p_wave_areas_4_SA_node_sorted),1);
for i = 1:length(I_p_wave_areas_4_SA_node_sorted)
    I_p_wave_areas_4_SA_node_dev_avg_min_2SD(i) = I_p_wave_areas_4_SA_node_sorted(i) - (I_p_wave_areas_4_SA_node_average-2*I_p_wave_areas_4_SA_node_SD);
end
I_p_wave_areas_4_SA_node_dev_avg_min_2SD = abs(I_p_wave_areas_4_SA_node_dev_avg_min_2SD);
for i = 1:length(I_p_wave_areas_4_SA_node_sorted)
    if I_p_wave_areas_4_SA_node_dev_avg_min_2SD(i) == min(I_p_wave_areas_4_SA_node_dev_avg_min_2SD)
        I_p_wave_areas_4_dev_avg_min_2SD_index = i;
    end
end
I_p_wave_areas_4_SA_node_dev_avg_pls_2SD = area_4(length(I_p_wave_areas_4_SA_node_sorted),1);
for i = 1:length(I_p_wave_areas_4_SA_node_sorted)
    I_p_wave_areas_4_SA_node_dev_avg_pls_2SD(i) = I_p_wave_areas_4_SA_node_sorted(i) - (I_p_wave_areas_4_SA_node_average+2*I_p_wave_areas_4_SA_node_SD);
end
I_p_wave_areas_4_SA_node_dev_avg_pls_2SD = abs(I_p_wave_areas_4_SA_node_dev_avg_pls_2SD);
for i = 1:length(I_p_wave_areas_4_SA_node_sorted)
    if I_p_wave_areas_4_SA_node_dev_avg_pls_2SD(i) == min(I_p_wave_areas_4_SA_node_dev_avg_pls_2SD)
        I_p_wave_areas_4_SA_node_dev_avg_pls_2SD_index = i;
    end
end
I_p_wave_areas_4_SA_node_dev_avg_min_3SD = area_4(length(I_p_wave_areas_4_SA_node_sorted),1);
for i = 1:length(I_p_wave_areas_4_SA_node_sorted)
    I_p_wave_areas_4_SA_node_dev_avg_min_3SD(i) = I_p_wave_areas_4_SA_node_sorted(i) - (I_p_wave_areas_4_SA_node_average-3*I_p_wave_areas_4_SA_node_SD);
end
I_p_wave_areas_4_SA_node_dev_avg_min_3SD = abs(I_p_wave_areas_4_SA_node_dev_avg_min_3SD);
for i = 1:length(I_p_wave_areas_4_SA_node_sorted)
    if I_p_wave_areas_4_SA_node_dev_avg_min_3SD(i) == min(I_p_wave_areas_4_SA_node_dev_avg_min_3SD)
        I_p_wave_areas_4_SA_node_dev_avg_min_3SD_index = i;
    end
end
I_p_wave_areas_4_SA_node_dev_avg_pls_3SD = area_4(length(I_p_wave_areas_4_SA_node_sorted),1);
for i = 1:length(I_p_wave_areas_4_SA_node_sorted)
    I_p_wave_areas_4_SA_node_dev_avg_pls_3SD(i) = I_p_wave_areas_4_SA_node_sorted(i) - (I_p_wave_areas_4_SA_node_average+3*I_p_wave_areas_4_SA_node_SD);
end
I_p_wave_areas_4_SA_node_dev_avg_pls_3SD = abs(I_p_wave_areas_4_SA_node_dev_avg_pls_3SD);
for i = 1:length(I_p_wave_areas_4_SA_node_sorted)
    if I_p_wave_areas_4_SA_node_dev_avg_pls_3SD(i) == min(I_p_wave_areas_4_SA_node_dev_avg_pls_3SD)
        I_p_wave_areas_4_SA_node_dev_avg_pls_3SD_index = i;
    end
end
I_p_wave_areas_4_SA_node_tot_number_of_indices = length(I_p_wave_areas_4_SA_node_sorted);
%This is the approximation of the tot "area under curve." 
I_p_wave_areas_4_SA_node_quotient_within_1SD = (I_p_wave_areas_4_SA_node_dev_avg_pls_1SD_index - I_p_wave_areas_4_SA_nodedev_avg_min_1SD_index)/I_p_wave_areas_4_SA_node_tot_number_of_indices;
I_p_wave_areas_4_SA_node_quotient_within_2SD = (I_p_wave_areas_4_SA_node_dev_avg_pls_2SD_index - I_p_wave_areas_4_dev_avg_min_2SD_index)/I_p_wave_areas_4_SA_node_tot_number_of_indices;
I_p_wave_areas_4_SA_node_quotient_within_3SD = (I_p_wave_areas_4_SA_node_dev_avg_pls_3SD_index - I_p_wave_areas_4_SA_node_dev_avg_min_3SD_index)/I_p_wave_areas_4_SA_node_tot_number_of_indices;
if ((I_p_wave_areas_4_SA_node_quotient_within_1SD > 0.66 && I_p_wave_areas_4_SA_node_quotient_within_1SD < 0.70) && (I_p_wave_areas_4_SA_node_quotient_within_2SD > 0.93 && I_p_wave_areas_4_SA_node_quotient_within_2SD < 0.97) && (I_p_wave_areas_4_SA_node_quotient_within_3SD > 0.98 && I_p_wave_areas_4_SA_node_quotient_within_3SD < 1)) 
    isNormal(10,1) = 1;
end
if I_p_wave_areas_4_SA_node_quotient_within_1SD == 0
    isSpike(10,1) = 1;
end
I_p_wave_areas_4_crista_terminalis_sorted = sort(I_p_wave_areas_4_crista_terminalis);
I_p_wave_areas_4_crista_terminalis_dev_avg_min_1SD = area_4(length(I_p_wave_areas_4_crista_terminalis_sorted),1);
for i = 1:length(I_p_wave_areas_4_crista_terminalis_sorted)
    I_p_wave_areas_4_crista_terminalis_dev_avg_min_1SD(i) = I_p_wave_areas_4_crista_terminalis_sorted(i) - (I_p_wave_areas_4_crista_terminalis_average-I_p_wave_areas_4_crista_terminalis_SD);
end
I_p_wave_areas_4_crista_terminalis_dev_avg_min_1SD = abs(I_p_wave_areas_4_crista_terminalis_dev_avg_min_1SD);
for i = 1:length(I_p_wave_areas_4_crista_terminalis_sorted)
    if I_p_wave_areas_4_crista_terminalis_dev_avg_min_1SD(i) == min(I_p_wave_areas_4_crista_terminalis_dev_avg_min_1SD)
        I_p_wave_areas_4_crista_terminalis_dev_avg_min_1SD_index = i;
    end
end
I_p_wave_areas_4_crista_terminalis_dev_avg_pls_1SD = area_4(length(I_p_wave_areas_4_crista_terminalis_sorted),1);
for i = 1:length(I_p_wave_areas_4_crista_terminalis_sorted)
    I_p_wave_areas_4_crista_terminalis_dev_avg_pls_1SD(i) = I_p_wave_areas_4_crista_terminalis_sorted(i) - (I_p_wave_areas_4_crista_terminalis_average+I_p_wave_areas_4_crista_terminalis_SD);
end
I_p_wave_areas_4_crista_terminalis_dev_avg_pls_1SD = abs(I_p_wave_areas_4_crista_terminalis_dev_avg_pls_1SD);
for i = 1:length(I_p_wave_areas_4_crista_terminalis_sorted)
    if I_p_wave_areas_4_crista_terminalis_dev_avg_pls_1SD(i) == min(I_p_wave_areas_4_crista_terminalis_dev_avg_pls_1SD)
        I_p_wave_areas_4_crista_terminalis_dev_avg_pls_1SD_index = i;
    end
end
I_p_wave_areas_4_crista_terminalis_dev_avg_min_2SD = area_4(length(I_p_wave_areas_4_crista_terminalis_sorted),1);
for i = 1:length(I_p_wave_areas_4_crista_terminalis_sorted)
    I_p_wave_areas_4_crista_terminalis_dev_avg_min_2SD(i) = I_p_wave_areas_4_crista_terminalis_sorted(i) - (I_p_wave_areas_4_crista_terminalis_average-2*I_p_wave_areas_4_crista_terminalis_SD);
end
I_p_wave_areas_4_crista_terminalis_dev_avg_min_2SD = abs(I_p_wave_areas_4_crista_terminalis_dev_avg_min_2SD);
for i = 1:length(I_p_wave_areas_4_crista_terminalis_sorted)
    if I_p_wave_areas_4_crista_terminalis_dev_avg_min_2SD(i) == min(I_p_wave_areas_4_crista_terminalis_dev_avg_min_2SD)
        I_p_wave_areas_4_dev_avg_min_2SD_index = i;
    end
end
I_p_wave_areas_4_crista_terminalis_dev_avg_pls_2SD = area_4(length(I_p_wave_areas_4_crista_terminalis_sorted),1);
for i = 1:length(I_p_wave_areas_4_crista_terminalis_sorted)
    I_p_wave_areas_4_crista_terminalis_dev_avg_pls_2SD(i) = I_p_wave_areas_4_crista_terminalis_sorted(i) - (I_p_wave_areas_4_crista_terminalis_average+2*I_p_wave_areas_4_crista_terminalis_SD);
end
I_p_wave_areas_4_crista_terminalis_dev_avg_pls_2SD = abs(I_p_wave_areas_4_crista_terminalis_dev_avg_pls_2SD);
for i = 1:length(I_p_wave_areas_4_crista_terminalis_sorted)
    if I_p_wave_areas_4_crista_terminalis_dev_avg_pls_2SD(i) == min(I_p_wave_areas_4_crista_terminalis_dev_avg_pls_2SD)
        I_p_wave_areas_4_crista_terminalis_dev_avg_pls_2SD_index = i;
    end
end
I_p_wave_areas_4_crista_terminalis_dev_avg_min_3SD = area_4(length(I_p_wave_areas_4_crista_terminalis_sorted),1);
for i = 1:length(I_p_wave_areas_4_crista_terminalis_sorted)
    I_p_wave_areas_4_crista_terminalis_dev_avg_min_3SD(i) = I_p_wave_areas_4_crista_terminalis_sorted(i) - (I_p_wave_areas_4_crista_terminalis_average-3*I_p_wave_areas_4_crista_terminalis_SD);
end
I_p_wave_areas_4_crista_terminalis_dev_avg_min_3SD = abs(I_p_wave_areas_4_crista_terminalis_dev_avg_min_3SD);
for i = 1:length(I_p_wave_areas_4_crista_terminalis_sorted)
    if I_p_wave_areas_4_crista_terminalis_dev_avg_min_3SD(i) == min(I_p_wave_areas_4_crista_terminalis_dev_avg_min_3SD)
        I_p_wave_areas_4_crista_terminalis_dev_avg_min_3SD_index = i;
    end
end
I_p_wave_areas_4_crista_terminalis_dev_avg_pls_3SD = area_4(length(I_p_wave_areas_4_crista_terminalis_sorted),1);
for i = 1:length(I_p_wave_areas_4_crista_terminalis_sorted)
    I_p_wave_areas_4_crista_terminalis_dev_avg_pls_3SD(i) = I_p_wave_areas_4_crista_terminalis_sorted(i) - (I_p_wave_areas_4_crista_terminalis_average+3*I_p_wave_areas_4_crista_terminalis_SD);
end
I_p_wave_areas_4_crista_terminalis_dev_avg_pls_3SD = abs(I_p_wave_areas_4_crista_terminalis_dev_avg_pls_3SD);
for i = 1:length(I_p_wave_areas_4_crista_terminalis_sorted)
    if I_p_wave_areas_4_crista_terminalis_dev_avg_pls_3SD(i) == min(I_p_wave_areas_4_crista_terminalis_dev_avg_pls_3SD)
        I_p_wave_areas_4_crista_terminalis_dev_avg_pls_3SD_index = i;
    end
end
I_p_wave_areas_4_crista_terminalis_tot_number_of_indices = length(I_p_wave_areas_4_crista_terminalis_sorted);
%This is the approximation of the tot "area under curve." 
I_p_wave_areas_4_crista_terminalis_quotient_within_1SD = (I_p_wave_areas_4_crista_terminalis_dev_avg_pls_1SD_index - I_p_wave_areas_4_crista_terminalisdev_avg_min_1SD_index)/I_p_wave_areas_4_crista_terminalis_tot_number_of_indices;
I_p_wave_areas_4_crista_terminalis_quotient_within_2SD = (I_p_wave_areas_4_crista_terminalis_dev_avg_pls_2SD_index - I_p_wave_areas_4_dev_avg_min_2SD_index)/I_p_wave_areas_4_crista_terminalis_tot_number_of_indices;
I_p_wave_areas_4_crista_terminalis_quotient_within_3SD = (I_p_wave_areas_4_crista_terminalis_dev_avg_pls_3SD_index - I_p_wave_areas_4_crista_terminalis_dev_avg_min_3SD_index)/I_p_wave_areas_4_crista_terminalis_tot_number_of_indices;
if (I_p_wave_areas_4_crista_terminalis_quotient_within_1SD > 0.66 && I_p_wave_areas_4_crista_terminalis_quotient_within_1SD < 0.70) && (I_p_wave_areas_4_crista_terminalis_quotient_within_2SD > 0.93 && I_p_wave_areas_4_crista_terminalis_quotient_within_2SD < 0.97) && (I_p_wave_areas_4_crista_terminalis_quotient_within_3SD > 0.98 && I_p_wave_areas_4_crista_terminalis_quotient_within_3SD < 1)
    isNormal(10,2) = 1;
end
if I_p_wave_areas_4_crista_terminalis_quotient_within_1SD == 0
    isSpike(10,2) = 1;
end
I_p_wave_areas_4_tricuspid_annulus_sorted = sort(I_p_wave_areas_4_tricuspid_annulus);
I_p_wave_areas_4_tricuspid_annulus_dev_avg_min_1SD = area_4(length(I_p_wave_areas_4_tricuspid_annulus_sorted),1);
for i = 1:length(I_p_wave_areas_4_tricuspid_annulus_sorted)
    I_p_wave_areas_4_tricuspid_annulus_dev_avg_min_1SD(i) = I_p_wave_areas_4_tricuspid_annulus_sorted(i) - (I_p_wave_areas_4_tricuspid_annulus_average-I_p_wave_areas_4_tricuspid_annulus_SD);
end
I_p_wave_areas_4_tricuspid_annulus_dev_avg_min_1SD = abs(I_p_wave_areas_4_tricuspid_annulus_dev_avg_min_1SD);
for i = 1:length(I_p_wave_areas_4_tricuspid_annulus_sorted)
    if I_p_wave_areas_4_tricuspid_annulus_dev_avg_min_1SD(i) == min(I_p_wave_areas_4_tricuspid_annulus_dev_avg_min_1SD)
        I_p_wave_areas_4_tricuspid_annulus_dev_avg_min_1SD_index = i;
    end
end
I_p_wave_areas_4_tricuspid_annulus_dev_avg_pls_1SD = area_4(length(I_p_wave_areas_4_tricuspid_annulus_sorted),1);
for i = 1:length(I_p_wave_areas_4_tricuspid_annulus_sorted)
    I_p_wave_areas_4_tricuspid_annulus_dev_avg_pls_1SD(i) = I_p_wave_areas_4_tricuspid_annulus_sorted(i) - (I_p_wave_areas_4_tricuspid_annulus_average+I_p_wave_areas_4_tricuspid_annulus_SD);
end
I_p_wave_areas_4_tricuspid_annulus_dev_avg_pls_1SD = abs(I_p_wave_areas_4_tricuspid_annulus_dev_avg_pls_1SD);
for i = 1:length(I_p_wave_areas_4_tricuspid_annulus_sorted)
    if I_p_wave_areas_4_tricuspid_annulus_dev_avg_pls_1SD(i) == min(I_p_wave_areas_4_tricuspid_annulus_dev_avg_pls_1SD)
        I_p_wave_areas_4_tricuspid_annulus_dev_avg_pls_1SD_index = i;
    end
end
I_p_wave_areas_4_tricuspid_annulus_dev_avg_min_2SD = area_4(length(I_p_wave_areas_4_tricuspid_annulus_sorted),1);
for i = 1:length(I_p_wave_areas_4_tricuspid_annulus_sorted)
    I_p_wave_areas_4_tricuspid_annulus_dev_avg_min_2SD(i) = I_p_wave_areas_4_tricuspid_annulus_sorted(i) - (I_p_wave_areas_4_tricuspid_annulus_average-2*I_p_wave_areas_4_tricuspid_annulus_SD);
end
I_p_wave_areas_4_tricuspid_annulus_dev_avg_min_2SD = abs(I_p_wave_areas_4_tricuspid_annulus_dev_avg_min_2SD);
for i = 1:length(I_p_wave_areas_4_tricuspid_annulus_sorted)
    if I_p_wave_areas_4_tricuspid_annulus_dev_avg_min_2SD(i) == min(I_p_wave_areas_4_tricuspid_annulus_dev_avg_min_2SD)
        I_p_wave_areas_4_dev_avg_min_2SD_index = i;
    end
end
I_p_wave_areas_4_tricuspid_annulus_dev_avg_pls_2SD = area_4(length(I_p_wave_areas_4_tricuspid_annulus_sorted),1);
for i = 1:length(I_p_wave_areas_4_tricuspid_annulus_sorted)
    I_p_wave_areas_4_tricuspid_annulus_dev_avg_pls_2SD(i) = I_p_wave_areas_4_tricuspid_annulus_sorted(i) - (I_p_wave_areas_4_tricuspid_annulus_average+2*I_p_wave_areas_4_tricuspid_annulus_SD);
end
I_p_wave_areas_4_tricuspid_annulus_dev_avg_pls_2SD = abs(I_p_wave_areas_4_tricuspid_annulus_dev_avg_pls_2SD);
for i = 1:length(I_p_wave_areas_4_tricuspid_annulus_sorted)
    if I_p_wave_areas_4_tricuspid_annulus_dev_avg_pls_2SD(i) == min(I_p_wave_areas_4_tricuspid_annulus_dev_avg_pls_2SD)
        I_p_wave_areas_4_tricuspid_annulus_dev_avg_pls_2SD_index = i;
    end
end
I_p_wave_areas_4_tricuspid_annulus_dev_avg_min_3SD = area_4(length(I_p_wave_areas_4_tricuspid_annulus_sorted),1);
for i = 1:length(I_p_wave_areas_4_tricuspid_annulus_sorted)
    I_p_wave_areas_4_tricuspid_annulus_dev_avg_min_3SD(i) = I_p_wave_areas_4_tricuspid_annulus_sorted(i) - (I_p_wave_areas_4_tricuspid_annulus_average-3*I_p_wave_areas_4_tricuspid_annulus_SD);
end
I_p_wave_areas_4_tricuspid_annulus_dev_avg_min_3SD = abs(I_p_wave_areas_4_tricuspid_annulus_dev_avg_min_3SD);
for i = 1:length(I_p_wave_areas_4_tricuspid_annulus_sorted)
    if I_p_wave_areas_4_tricuspid_annulus_dev_avg_min_3SD(i) == min(I_p_wave_areas_4_tricuspid_annulus_dev_avg_min_3SD)
        I_p_wave_areas_4_tricuspid_annulus_dev_avg_min_3SD_index = i;
    end
end
I_p_wave_areas_4_tricuspid_annulus_dev_avg_pls_3SD = area_4(length(I_p_wave_areas_4_tricuspid_annulus_sorted),1);
for i = 1:length(I_p_wave_areas_4_tricuspid_annulus_sorted)
    I_p_wave_areas_4_tricuspid_annulus_dev_avg_pls_3SD(i) = I_p_wave_areas_4_tricuspid_annulus_sorted(i) - (I_p_wave_areas_4_tricuspid_annulus_average+3*I_p_wave_areas_4_tricuspid_annulus_SD);
end
I_p_wave_areas_4_tricuspid_annulus_dev_avg_pls_3SD = abs(I_p_wave_areas_4_tricuspid_annulus_dev_avg_pls_3SD);
for i = 1:length(I_p_wave_areas_4_tricuspid_annulus_sorted)
    if I_p_wave_areas_4_tricuspid_annulus_dev_avg_pls_3SD(i) == min(I_p_wave_areas_4_tricuspid_annulus_dev_avg_pls_3SD)
        I_p_wave_areas_4_tricuspid_annulus_dev_avg_pls_3SD_index = i;
    end
end
I_p_wave_areas_4_tricuspid_annulus_tot_number_of_indices = length(I_p_wave_areas_4_tricuspid_annulus_sorted);
%This is the approximation of the tot "area under curve." 
I_p_wave_areas_4_tricuspid_annulus_quotient_within_1SD = (I_p_wave_areas_4_tricuspid_annulus_dev_avg_pls_1SD_index - I_p_wave_areas_4_tricuspid_annulusdev_avg_min_1SD_index)/I_p_wave_areas_4_tricuspid_annulus_tot_number_of_indices;
I_p_wave_areas_4_tricuspid_annulus_quotient_within_2SD = (I_p_wave_areas_4_tricuspid_annulus_dev_avg_pls_2SD_index - I_p_wave_areas_4_dev_avg_min_2SD_index)/I_p_wave_areas_4_tricuspid_annulus_tot_number_of_indices;
I_p_wave_areas_4_tricuspid_annulus_quotient_within_3SD = (I_p_wave_areas_4_tricuspid_annulus_dev_avg_pls_3SD_index - I_p_wave_areas_4_tricuspid_annulus_dev_avg_min_3SD_index)/I_p_wave_areas_4_tricuspid_annulus_tot_number_of_indices;
if (I_p_wave_areas_4_tricuspid_annulus_quotient_within_1SD > 0.66 && I_p_wave_areas_4_tricuspid_annulus_quotient_within_1SD < 0.70) && (I_p_wave_areas_4_tricuspid_annulus_quotient_within_2SD > 0.93 && I_p_wave_areas_4_tricuspid_annulus_quotient_within_2SD < 0.97) && (I_p_wave_areas_4_tricuspid_annulus_quotient_within_3SD > 0.98 && I_p_wave_areas_4_tricuspid_annulus_quotient_within_3SD < 1)
    isNormal(10,3) = 1;
end
if I_p_wave_areas_4_tricuspid_annulus_quotient_within_1SD == 0
    isSpike(10,3) = 1;
end
I_p_wave_areas_4_coronary_sinus_sorted = sort(I_p_wave_areas_4_coronary_sinus);
I_p_wave_areas_4_coronary_sinus_dev_avg_min_1SD = area_4(length(I_p_wave_areas_4_coronary_sinus_sorted),1);
for i = 1:length(I_p_wave_areas_4_coronary_sinus_sorted)
    I_p_wave_areas_4_coronary_sinus_dev_avg_min_1SD(i) = I_p_wave_areas_4_coronary_sinus_sorted(i) - (I_p_wave_areas_4_coronary_sinus_average-I_p_wave_areas_4_coronary_sinus_SD);
end
I_p_wave_areas_4_coronary_sinus_dev_avg_min_1SD = abs(I_p_wave_areas_4_coronary_sinus_dev_avg_min_1SD);
for i = 1:length(I_p_wave_areas_4_coronary_sinus_sorted)
    if I_p_wave_areas_4_coronary_sinus_dev_avg_min_1SD(i) == min(I_p_wave_areas_4_coronary_sinus_dev_avg_min_1SD)
        I_p_wave_areas_4_coronary_sinus_dev_avg_min_1SD_index = i;
    end
end
I_p_wave_areas_4_coronary_sinus_dev_avg_pls_1SD = area_4(length(I_p_wave_areas_4_coronary_sinus_sorted),1);
for i = 1:length(I_p_wave_areas_4_coronary_sinus_sorted)
    I_p_wave_areas_4_coronary_sinus_dev_avg_pls_1SD(i) = I_p_wave_areas_4_coronary_sinus_sorted(i) - (I_p_wave_areas_4_coronary_sinus_average+I_p_wave_areas_4_coronary_sinus_SD);
end
I_p_wave_areas_4_coronary_sinus_dev_avg_pls_1SD = abs(I_p_wave_areas_4_coronary_sinus_dev_avg_pls_1SD);
for i = 1:length(I_p_wave_areas_4_coronary_sinus_sorted)
    if I_p_wave_areas_4_coronary_sinus_dev_avg_pls_1SD(i) == min(I_p_wave_areas_4_coronary_sinus_dev_avg_pls_1SD)
        I_p_wave_areas_4_coronary_sinus_dev_avg_pls_1SD_index = i;
    end
end
I_p_wave_areas_4_coronary_sinus_dev_avg_min_2SD = area_4(length(I_p_wave_areas_4_coronary_sinus_sorted),1);
for i = 1:length(I_p_wave_areas_4_coronary_sinus_sorted)
    I_p_wave_areas_4_coronary_sinus_dev_avg_min_2SD(i) = I_p_wave_areas_4_coronary_sinus_sorted(i) - (I_p_wave_areas_4_coronary_sinus_average-2*I_p_wave_areas_4_coronary_sinus_SD);
end
I_p_wave_areas_4_coronary_sinus_dev_avg_min_2SD = abs(I_p_wave_areas_4_coronary_sinus_dev_avg_min_2SD);
for i = 1:length(I_p_wave_areas_4_coronary_sinus_sorted)
    if I_p_wave_areas_4_coronary_sinus_dev_avg_min_2SD(i) == min(I_p_wave_areas_4_coronary_sinus_dev_avg_min_2SD)
        I_p_wave_areas_4_dev_avg_min_2SD_index = i;
    end
end
I_p_wave_areas_4_coronary_sinus_dev_avg_pls_2SD = area_4(length(I_p_wave_areas_4_coronary_sinus_sorted),1);
for i = 1:length(I_p_wave_areas_4_coronary_sinus_sorted)
    I_p_wave_areas_4_coronary_sinus_dev_avg_pls_2SD(i) = I_p_wave_areas_4_coronary_sinus_sorted(i) - (I_p_wave_areas_4_coronary_sinus_average+2*I_p_wave_areas_4_coronary_sinus_SD);
end
I_p_wave_areas_4_coronary_sinus_dev_avg_pls_2SD = abs(I_p_wave_areas_4_coronary_sinus_dev_avg_pls_2SD);
for i = 1:length(I_p_wave_areas_4_coronary_sinus_sorted)
    if I_p_wave_areas_4_coronary_sinus_dev_avg_pls_2SD(i) == min(I_p_wave_areas_4_coronary_sinus_dev_avg_pls_2SD)
        I_p_wave_areas_4_coronary_sinus_dev_avg_pls_2SD_index = i;
    end
end
I_p_wave_areas_4_coronary_sinus_dev_avg_min_3SD = area_4(length(I_p_wave_areas_4_coronary_sinus_sorted),1);
for i = 1:length(I_p_wave_areas_4_coronary_sinus_sorted)
    I_p_wave_areas_4_coronary_sinus_dev_avg_min_3SD(i) = I_p_wave_areas_4_coronary_sinus_sorted(i) - (I_p_wave_areas_4_coronary_sinus_average-3*I_p_wave_areas_4_coronary_sinus_SD);
end
I_p_wave_areas_4_coronary_sinus_dev_avg_min_3SD = abs(I_p_wave_areas_4_coronary_sinus_dev_avg_min_3SD);
for i = 1:length(I_p_wave_areas_4_coronary_sinus_sorted)
    if I_p_wave_areas_4_coronary_sinus_dev_avg_min_3SD(i) == min(I_p_wave_areas_4_coronary_sinus_dev_avg_min_3SD)
        I_p_wave_areas_4_coronary_sinus_dev_avg_min_3SD_index = i;
    end
end
I_p_wave_areas_4_coronary_sinus_dev_avg_pls_3SD = area_4(length(I_p_wave_areas_4_coronary_sinus_sorted),1);
for i = 1:length(I_p_wave_areas_4_coronary_sinus_sorted)
    I_p_wave_areas_4_coronary_sinus_dev_avg_pls_3SD(i) = I_p_wave_areas_4_coronary_sinus_sorted(i) - (I_p_wave_areas_4_coronary_sinus_average+3*I_p_wave_areas_4_coronary_sinus_SD);
end
I_p_wave_areas_4_coronary_sinus_dev_avg_pls_3SD = abs(I_p_wave_areas_4_coronary_sinus_dev_avg_pls_3SD);
for i = 1:length(I_p_wave_areas_4_coronary_sinus_sorted)
    if I_p_wave_areas_4_coronary_sinus_dev_avg_pls_3SD(i) == min(I_p_wave_areas_4_coronary_sinus_dev_avg_pls_3SD)
        I_p_wave_areas_4_coronary_sinus_dev_avg_pls_3SD_index = i;
    end
end
I_p_wave_areas_4_coronary_sinus_tot_number_of_indices = length(I_p_wave_areas_4_coronary_sinus_sorted);
%This is the approximation of the tot "area under curve." 
I_p_wave_areas_4_coronary_sinus_quotient_within_1SD = (I_p_wave_areas_4_coronary_sinus_dev_avg_pls_1SD_index - I_p_wave_areas_4_coronary_sinusdev_avg_min_1SD_index)/I_p_wave_areas_4_coronary_sinus_tot_number_of_indices;
I_p_wave_areas_4_coronary_sinus_quotient_within_2SD = (I_p_wave_areas_4_coronary_sinus_dev_avg_pls_2SD_index - I_p_wave_areas_4_dev_avg_min_2SD_index)/I_p_wave_areas_4_coronary_sinus_tot_number_of_indices;
I_p_wave_areas_4_coronary_sinus_quotient_within_3SD = (I_p_wave_areas_4_coronary_sinus_dev_avg_pls_3SD_index - I_p_wave_areas_4_coronary_sinus_dev_avg_min_3SD_index)/I_p_wave_areas_4_coronary_sinus_tot_number_of_indices;
if (I_p_wave_areas_4_coronary_sinus_quotient_within_1SD > 0.66 && I_p_wave_areas_4_coronary_sinus_quotient_within_1SD < 0.70) && (I_p_wave_areas_4_coronary_sinus_quotient_within_2SD > 0.93 && I_p_wave_areas_4_coronary_sinus_quotient_within_2SD < 0.97) && (I_p_wave_areas_4_coronary_sinus_quotient_within_3SD > 0.98 && I_p_wave_areas_4_coronary_sinus_quotient_within_3SD < 1)
    isNormal(10,4) = 1;
end
if I_p_wave_areas_4_coronary_sinus_quotient_within_1SD == 0
    isSpike(10,4) = 1;
end
I_p_wave_areas_4_ostium_sorted = sort(I_p_wave_areas_4_ostium);
I_p_wave_areas_4_ostium_dev_avg_min_1SD = area_4(length(I_p_wave_areas_4_ostium_sorted),1);
for i = 1:length(I_p_wave_areas_4_ostium_sorted)
    I_p_wave_areas_4_ostium_dev_avg_min_1SD(i) = I_p_wave_areas_4_ostium_sorted(i) - (I_p_wave_areas_4_ostium_average-I_p_wave_areas_4_ostium_SD);
end
I_p_wave_areas_4_ostium_dev_avg_min_1SD = abs(I_p_wave_areas_4_ostium_dev_avg_min_1SD);
for i = 1:length(I_p_wave_areas_4_ostium_sorted)
    if I_p_wave_areas_4_ostium_dev_avg_min_1SD(i) == min(I_p_wave_areas_4_ostium_dev_avg_min_1SD)
        I_p_wave_areas_4_ostium_dev_avg_min_1SD_index = i;
    end
end
I_p_wave_areas_4_ostium_dev_avg_pls_1SD = area_4(length(I_p_wave_areas_4_ostium_sorted),1);
for i = 1:length(I_p_wave_areas_4_ostium_sorted)
    I_p_wave_areas_4_ostium_dev_avg_pls_1SD(i) = I_p_wave_areas_4_ostium_sorted(i) - (I_p_wave_areas_4_ostium_average+I_p_wave_areas_4_ostium_SD);
end
I_p_wave_areas_4_ostium_dev_avg_pls_1SD = abs(I_p_wave_areas_4_ostium_dev_avg_pls_1SD);
for i = 1:length(I_p_wave_areas_4_ostium_sorted)
    if I_p_wave_areas_4_ostium_dev_avg_pls_1SD(i) == min(I_p_wave_areas_4_ostium_dev_avg_pls_1SD)
        I_p_wave_areas_4_ostium_dev_avg_pls_1SD_index = i;
    end
end
I_p_wave_areas_4_ostium_dev_avg_min_2SD = area_4(length(I_p_wave_areas_4_ostium_sorted),1);
for i = 1:length(I_p_wave_areas_4_ostium_sorted)
    I_p_wave_areas_4_ostium_dev_avg_min_2SD(i) = I_p_wave_areas_4_ostium_sorted(i) - (I_p_wave_areas_4_ostium_average-2*I_p_wave_areas_4_ostium_SD);
end
I_p_wave_areas_4_ostium_dev_avg_min_2SD = abs(I_p_wave_areas_4_ostium_dev_avg_min_2SD);
for i = 1:length(I_p_wave_areas_4_ostium_sorted)
    if I_p_wave_areas_4_ostium_dev_avg_min_2SD(i) == min(I_p_wave_areas_4_ostium_dev_avg_min_2SD)
        I_p_wave_areas_4_dev_avg_min_2SD_index = i;
    end
end
I_p_wave_areas_4_ostium_dev_avg_pls_2SD = area_4(length(I_p_wave_areas_4_ostium_sorted),1);
for i = 1:length(I_p_wave_areas_4_ostium_sorted)
    I_p_wave_areas_4_ostium_dev_avg_pls_2SD(i) = I_p_wave_areas_4_ostium_sorted(i) - (I_p_wave_areas_4_ostium_average+2*I_p_wave_areas_4_ostium_SD);
end
I_p_wave_areas_4_ostium_dev_avg_pls_2SD = abs(I_p_wave_areas_4_ostium_dev_avg_pls_2SD);
for i = 1:length(I_p_wave_areas_4_ostium_sorted)
    if I_p_wave_areas_4_ostium_dev_avg_pls_2SD(i) == min(I_p_wave_areas_4_ostium_dev_avg_pls_2SD)
        I_p_wave_areas_4_ostium_dev_avg_pls_2SD_index = i;
    end
end
I_p_wave_areas_4_ostium_dev_avg_min_3SD = area_4(length(I_p_wave_areas_4_ostium_sorted),1);
for i = 1:length(I_p_wave_areas_4_ostium_sorted)
    I_p_wave_areas_4_ostium_dev_avg_min_3SD(i) = I_p_wave_areas_4_ostium_sorted(i) - (I_p_wave_areas_4_ostium_average-3*I_p_wave_areas_4_ostium_SD);
end
I_p_wave_areas_4_ostium_dev_avg_min_3SD = abs(I_p_wave_areas_4_ostium_dev_avg_min_3SD);
for i = 1:length(I_p_wave_areas_4_ostium_sorted)
    if I_p_wave_areas_4_ostium_dev_avg_min_3SD(i) == min(I_p_wave_areas_4_ostium_dev_avg_min_3SD)
        I_p_wave_areas_4_ostium_dev_avg_min_3SD_index = i;
    end
end
I_p_wave_areas_4_ostium_dev_avg_pls_3SD = area_4(length(I_p_wave_areas_4_ostium_sorted),1);
for i = 1:length(I_p_wave_areas_4_ostium_sorted)
    I_p_wave_areas_4_ostium_dev_avg_pls_3SD(i) = I_p_wave_areas_4_ostium_sorted(i) - (I_p_wave_areas_4_ostium_average+3*I_p_wave_areas_4_ostium_SD);
end
I_p_wave_areas_4_ostium_dev_avg_pls_3SD = abs(I_p_wave_areas_4_ostium_dev_avg_pls_3SD);
for i = 1:length(I_p_wave_areas_4_ostium_sorted)
    if I_p_wave_areas_4_ostium_dev_avg_pls_3SD(i) == min(I_p_wave_areas_4_ostium_dev_avg_pls_3SD)
        I_p_wave_areas_4_ostium_dev_avg_pls_3SD_index = i;
    end
end
I_p_wave_areas_4_ostium_tot_number_of_indices = length(I_p_wave_areas_4_ostium_sorted);
%This is the approximation of the tot "area under curve." 
I_p_wave_areas_4_ostium_quotient_within_1SD = (I_p_wave_areas_4_ostium_dev_avg_pls_1SD_index - I_p_wave_areas_4_ostiumdev_avg_min_1SD_index)/I_p_wave_areas_4_ostium_tot_number_of_indices;
I_p_wave_areas_4_ostium_quotient_within_2SD = (I_p_wave_areas_4_ostium_dev_avg_pls_2SD_index - I_p_wave_areas_4_dev_avg_min_2SD_index)/I_p_wave_areas_4_ostium_tot_number_of_indices;
I_p_wave_areas_4_ostium_quotient_within_3SD = (I_p_wave_areas_4_ostium_dev_avg_pls_3SD_index - I_p_wave_areas_4_ostium_dev_avg_min_3SD_index)/I_p_wave_areas_4_ostium_tot_number_of_indices;
if (I_p_wave_areas_4_ostium_quotient_within_1SD > 0.66 && I_p_wave_areas_4_ostium_quotient_within_1SD < 0.70) && (I_p_wave_areas_4_ostium_quotient_within_2SD > 0.93 && I_p_wave_areas_4_ostium_quotient_within_2SD < 0.97) && (I_p_wave_areas_4_ostium_quotient_within_3SD > 0.98 && I_p_wave_areas_4_ostium_quotient_within_3SD < 1)
    isNormal(10,5) = 1;
end
if I_p_wave_areas_4_ostium_quotient_within_1SD == 0
    isSpike(10,5) = 1;
end
I_p_wave_areas_4_perinodal_sorted = sort(I_p_wave_areas_4_perinodal);
I_p_wave_areas_4_perinodal_dev_avg_min_1SD = area_4(length(I_p_wave_areas_4_perinodal_sorted),1);
for i = 1:length(I_p_wave_areas_4_perinodal_sorted)
    I_p_wave_areas_4_perinodal_dev_avg_min_1SD(i) = I_p_wave_areas_4_perinodal_sorted(i) - (I_p_wave_areas_4_perinodal_average-I_p_wave_areas_4_perinodal_SD);
end
I_p_wave_areas_4_perinodal_dev_avg_min_1SD = abs(I_p_wave_areas_4_perinodal_dev_avg_min_1SD);
for i = 1:length(I_p_wave_areas_4_perinodal_sorted)
    if I_p_wave_areas_4_perinodal_dev_avg_min_1SD(i) == min(I_p_wave_areas_4_perinodal_dev_avg_min_1SD)
        I_p_wave_areas_4_perinodal_dev_avg_min_1SD_index = i;
    end
end
I_p_wave_areas_4_perinodal_dev_avg_pls_1SD = area_4(length(I_p_wave_areas_4_perinodal_sorted),1);
for i = 1:length(I_p_wave_areas_4_perinodal_sorted)
    I_p_wave_areas_4_perinodal_dev_avg_pls_1SD(i) = I_p_wave_areas_4_perinodal_sorted(i) - (I_p_wave_areas_4_perinodal_average+I_p_wave_areas_4_perinodal_SD);
end
I_p_wave_areas_4_perinodal_dev_avg_pls_1SD = abs(I_p_wave_areas_4_perinodal_dev_avg_pls_1SD);
for i = 1:length(I_p_wave_areas_4_perinodal_sorted)
    if I_p_wave_areas_4_perinodal_dev_avg_pls_1SD(i) == min(I_p_wave_areas_4_perinodal_dev_avg_pls_1SD)
        I_p_wave_areas_4_perinodal_dev_avg_pls_1SD_index = i;
    end
end
I_p_wave_areas_4_perinodal_dev_avg_min_2SD = area_4(length(I_p_wave_areas_4_perinodal_sorted),1);
for i = 1:length(I_p_wave_areas_4_perinodal_sorted)
    I_p_wave_areas_4_perinodal_dev_avg_min_2SD(i) = I_p_wave_areas_4_perinodal_sorted(i) - (I_p_wave_areas_4_perinodal_average-2*I_p_wave_areas_4_perinodal_SD);
end
I_p_wave_areas_4_perinodal_dev_avg_min_2SD = abs(I_p_wave_areas_4_perinodal_dev_avg_min_2SD);
for i = 1:length(I_p_wave_areas_4_perinodal_sorted)
    if I_p_wave_areas_4_perinodal_dev_avg_min_2SD(i) == min(I_p_wave_areas_4_perinodal_dev_avg_min_2SD)
        I_p_wave_areas_4_dev_avg_min_2SD_index = i;
    end
end
I_p_wave_areas_4_perinodal_dev_avg_pls_2SD = area_4(length(I_p_wave_areas_4_perinodal_sorted),1);
for i = 1:length(I_p_wave_areas_4_perinodal_sorted)
    I_p_wave_areas_4_perinodal_dev_avg_pls_2SD(i) = I_p_wave_areas_4_perinodal_sorted(i) - (I_p_wave_areas_4_perinodal_average+2*I_p_wave_areas_4_perinodal_SD);
end
I_p_wave_areas_4_perinodal_dev_avg_pls_2SD = abs(I_p_wave_areas_4_perinodal_dev_avg_pls_2SD);
for i = 1:length(I_p_wave_areas_4_perinodal_sorted)
    if I_p_wave_areas_4_perinodal_dev_avg_pls_2SD(i) == min(I_p_wave_areas_4_perinodal_dev_avg_pls_2SD)
        I_p_wave_areas_4_perinodal_dev_avg_pls_2SD_index = i;
    end
end
I_p_wave_areas_4_perinodal_dev_avg_min_3SD = area_4(length(I_p_wave_areas_4_perinodal_sorted),1);
for i = 1:length(I_p_wave_areas_4_perinodal_sorted)
    I_p_wave_areas_4_perinodal_dev_avg_min_3SD(i) = I_p_wave_areas_4_perinodal_sorted(i) - (I_p_wave_areas_4_perinodal_average-3*I_p_wave_areas_4_perinodal_SD);
end
I_p_wave_areas_4_perinodal_dev_avg_min_3SD = abs(I_p_wave_areas_4_perinodal_dev_avg_min_3SD);
for i = 1:length(I_p_wave_areas_4_perinodal_sorted)
    if I_p_wave_areas_4_perinodal_dev_avg_min_3SD(i) == min(I_p_wave_areas_4_perinodal_dev_avg_min_3SD)
        I_p_wave_areas_4_perinodal_dev_avg_min_3SD_index = i;
    end
end
I_p_wave_areas_4_perinodal_dev_avg_pls_3SD = area_4(length(I_p_wave_areas_4_perinodal_sorted),1);
for i = 1:length(I_p_wave_areas_4_perinodal_sorted)
    I_p_wave_areas_4_perinodal_dev_avg_pls_3SD(i) = I_p_wave_areas_4_perinodal_sorted(i) - (I_p_wave_areas_4_perinodal_average+3*I_p_wave_areas_4_perinodal_SD);
end
I_p_wave_areas_4_perinodal_dev_avg_pls_3SD = abs(I_p_wave_areas_4_perinodal_dev_avg_pls_3SD);
for i = 1:length(I_p_wave_areas_4_perinodal_sorted)
    if I_p_wave_areas_4_perinodal_dev_avg_pls_3SD(i) == min(I_p_wave_areas_4_perinodal_dev_avg_pls_3SD)
        I_p_wave_areas_4_perinodal_dev_avg_pls_3SD_index = i;
    end
end
I_p_wave_areas_4_perinodal_tot_number_of_indices = length(I_p_wave_areas_4_perinodal_sorted);
%This is the approximation of the tot "area under curve." 
I_p_wave_areas_4_perinodal_quotient_within_1SD = (I_p_wave_areas_4_perinodal_dev_avg_pls_1SD_index - I_p_wave_areas_4_perinodaldev_avg_min_1SD_index)/I_p_wave_areas_4_perinodal_tot_number_of_indices;
I_p_wave_areas_4_perinodal_quotient_within_2SD = (I_p_wave_areas_4_perinodal_dev_avg_pls_2SD_index - I_p_wave_areas_4_dev_avg_min_2SD_index)/I_p_wave_areas_4_perinodal_tot_number_of_indices;
I_p_wave_areas_4_perinodal_quotient_within_3SD = (I_p_wave_areas_4_perinodal_dev_avg_pls_3SD_index - I_p_wave_areas_4_perinodal_dev_avg_min_3SD_index)/I_p_wave_areas_4_perinodal_tot_number_of_indices;
if (I_p_wave_areas_4_perinodal_quotient_within_1SD > 0.66 && I_p_wave_areas_4_perinodal_quotient_within_1SD < 0.70) && (I_p_wave_areas_4_perinodal_quotient_within_2SD > 0.93 && I_p_wave_areas_4_perinodal_quotient_within_2SD < 0.97) && (I_p_wave_areas_4_perinodal_quotient_within_3SD > 0.98 && I_p_wave_areas_4_perinodal_quotient_within_3SD < 1)
    isNormal(10,6) = 1;
end
if I_p_wave_areas_4_perinodal_quotient_within_1SD == 0
    isSpike(10,6) = 1;
end
I_p_wave_areas_4_right_septum_sorted = sort(I_p_wave_areas_4_right_septum);
I_p_wave_areas_4_right_septum_dev_avg_min_1SD = area_4(length(I_p_wave_areas_4_right_septum_sorted),1);
for i = 1:length(I_p_wave_areas_4_right_septum_sorted)
    I_p_wave_areas_4_right_septum_dev_avg_min_1SD(i) = I_p_wave_areas_4_right_septum_sorted(i) - (I_p_wave_areas_4_right_septum_average-I_p_wave_areas_4_right_septum_SD);
end
I_p_wave_areas_4_right_septum_dev_avg_min_1SD = abs(I_p_wave_areas_4_right_septum_dev_avg_min_1SD);
for i = 1:length(I_p_wave_areas_4_right_septum_sorted)
    if I_p_wave_areas_4_right_septum_dev_avg_min_1SD(i) == min(I_p_wave_areas_4_right_septum_dev_avg_min_1SD)
        I_p_wave_areas_4_right_septum_dev_avg_min_1SD_index = i;
    end
end
I_p_wave_areas_4_right_septum_dev_avg_pls_1SD = area_4(length(I_p_wave_areas_4_right_septum_sorted),1);
for i = 1:length(I_p_wave_areas_4_right_septum_sorted)
    I_p_wave_areas_4_right_septum_dev_avg_pls_1SD(i) = I_p_wave_areas_4_right_septum_sorted(i) - (I_p_wave_areas_4_right_septum_average+I_p_wave_areas_4_right_septum_SD);
end
I_p_wave_areas_4_right_septum_dev_avg_pls_1SD = abs(I_p_wave_areas_4_right_septum_dev_avg_pls_1SD);
for i = 1:length(I_p_wave_areas_4_right_septum_sorted)
    if I_p_wave_areas_4_right_septum_dev_avg_pls_1SD(i) == min(I_p_wave_areas_4_right_septum_dev_avg_pls_1SD)
        I_p_wave_areas_4_right_septum_dev_avg_pls_1SD_index = i;
    end
end
I_p_wave_areas_4_right_septum_dev_avg_min_2SD = area_4(length(I_p_wave_areas_4_right_septum_sorted),1);
for i = 1:length(I_p_wave_areas_4_right_septum_sorted)
    I_p_wave_areas_4_right_septum_dev_avg_min_2SD(i) = I_p_wave_areas_4_right_septum_sorted(i) - (I_p_wave_areas_4_right_septum_average-2*I_p_wave_areas_4_right_septum_SD);
end
I_p_wave_areas_4_right_septum_dev_avg_min_2SD = abs(I_p_wave_areas_4_right_septum_dev_avg_min_2SD);
for i = 1:length(I_p_wave_areas_4_right_septum_sorted)
    if I_p_wave_areas_4_right_septum_dev_avg_min_2SD(i) == min(I_p_wave_areas_4_right_septum_dev_avg_min_2SD)
        I_p_wave_areas_4_dev_avg_min_2SD_index = i;
    end
end
I_p_wave_areas_4_right_septum_dev_avg_pls_2SD = area_4(length(I_p_wave_areas_4_right_septum_sorted),1);
for i = 1:length(I_p_wave_areas_4_right_septum_sorted)
    I_p_wave_areas_4_right_septum_dev_avg_pls_2SD(i) = I_p_wave_areas_4_right_septum_sorted(i) - (I_p_wave_areas_4_right_septum_average+2*I_p_wave_areas_4_right_septum_SD);
end
I_p_wave_areas_4_right_septum_dev_avg_pls_2SD = abs(I_p_wave_areas_4_right_septum_dev_avg_pls_2SD);
for i = 1:length(I_p_wave_areas_4_right_septum_sorted)
    if I_p_wave_areas_4_right_septum_dev_avg_pls_2SD(i) == min(I_p_wave_areas_4_right_septum_dev_avg_pls_2SD)
        I_p_wave_areas_4_right_septum_dev_avg_pls_2SD_index = i;
    end
end
I_p_wave_areas_4_right_septum_dev_avg_min_3SD = area_4(length(I_p_wave_areas_4_right_septum_sorted),1);
for i = 1:length(I_p_wave_areas_4_right_septum_sorted)
    I_p_wave_areas_4_right_septum_dev_avg_min_3SD(i) = I_p_wave_areas_4_right_septum_sorted(i) - (I_p_wave_areas_4_right_septum_average-3*I_p_wave_areas_4_right_septum_SD);
end
I_p_wave_areas_4_right_septum_dev_avg_min_3SD = abs(I_p_wave_areas_4_right_septum_dev_avg_min_3SD);
for i = 1:length(I_p_wave_areas_4_right_septum_sorted)
    if I_p_wave_areas_4_right_septum_dev_avg_min_3SD(i) == min(I_p_wave_areas_4_right_septum_dev_avg_min_3SD)
        I_p_wave_areas_4_right_septum_dev_avg_min_3SD_index = i;
    end
end
I_p_wave_areas_4_right_septum_dev_avg_pls_3SD = area_4(length(I_p_wave_areas_4_right_septum_sorted),1);
for i = 1:length(I_p_wave_areas_4_right_septum_sorted)
    I_p_wave_areas_4_right_septum_dev_avg_pls_3SD(i) = I_p_wave_areas_4_right_septum_sorted(i) - (I_p_wave_areas_4_right_septum_average+3*I_p_wave_areas_4_right_septum_SD);
end
I_p_wave_areas_4_right_septum_dev_avg_pls_3SD = abs(I_p_wave_areas_4_right_septum_dev_avg_pls_3SD);
for i = 1:length(I_p_wave_areas_4_right_septum_sorted)
    if I_p_wave_areas_4_right_septum_dev_avg_pls_3SD(i) == min(I_p_wave_areas_4_right_septum_dev_avg_pls_3SD)
        I_p_wave_areas_4_right_septum_dev_avg_pls_3SD_index = i;
    end
end
I_p_wave_areas_4_right_septum_tot_number_of_indices = length(I_p_wave_areas_4_right_septum_sorted);
%This is the approximation of the tot "area under curve." 
I_p_wave_areas_4_right_septum_quotient_within_1SD = (I_p_wave_areas_4_right_septum_dev_avg_pls_1SD_index - I_p_wave_areas_4_right_septumdev_avg_min_1SD_index)/I_p_wave_areas_4_right_septum_tot_number_of_indices;
I_p_wave_areas_4_right_septum_quotient_within_2SD = (I_p_wave_areas_4_right_septum_dev_avg_pls_2SD_index - I_p_wave_areas_4_dev_avg_min_2SD_index)/I_p_wave_areas_4_right_septum_tot_number_of_indices;
I_p_wave_areas_4_right_septum_quotient_within_3SD = (I_p_wave_areas_4_right_septum_dev_avg_pls_3SD_index - I_p_wave_areas_4_right_septum_dev_avg_min_3SD_index)/I_p_wave_areas_4_right_septum_tot_number_of_indices;
if (I_p_wave_areas_4_right_septum_quotient_within_1SD > 0.66 && I_p_wave_areas_4_right_septum_quotient_within_1SD < 0.70) && (I_p_wave_areas_4_right_septum_quotient_within_2SD > 0.93 && I_p_wave_areas_4_right_septum_quotient_within_2SD < 0.97) && (I_p_wave_areas_4_right_septum_quotient_within_3SD > 0.98 && I_p_wave_areas_4_right_septum_quotient_within_3SD < 1)
    isNormal(10,7) = 1;
end
if I_p_wave_areas_4_right_septum_quotient_within_1SD == 0
    isSpike(10,7) = 1;
end
I_p_wave_areas_4_right_atrial_appendage_sorted = sort(I_p_wave_areas_4_right_atrial_appendage);
I_p_wave_areas_4_right_atrial_appendage_dev_avg_min_1SD = area_4(length(I_p_wave_areas_4_right_atrial_appendage_sorted),1);
for i = 1:length(I_p_wave_areas_4_right_atrial_appendage_sorted)
    I_p_wave_areas_4_right_atrial_appendage_dev_avg_min_1SD(i) = I_p_wave_areas_4_right_atrial_appendage_sorted(i) - (I_p_wave_areas_4_right_atrial_appendage_average-I_p_wave_areas_4_right_atrial_appendage_SD);
end
I_p_wave_areas_4_right_atrial_appendage_dev_avg_min_1SD = abs(I_p_wave_areas_4_right_atrial_appendage_dev_avg_min_1SD);
for i = 1:length(I_p_wave_areas_4_right_atrial_appendage_sorted)
    if I_p_wave_areas_4_right_atrial_appendage_dev_avg_min_1SD(i) == min(I_p_wave_areas_4_right_atrial_appendage_dev_avg_min_1SD)
        I_p_wave_areas_4_right_atrial_appendage_dev_avg_min_1SD_index = i;
    end
end
I_p_wave_areas_4_right_atrial_appendage_dev_avg_pls_1SD = area_4(length(I_p_wave_areas_4_right_atrial_appendage_sorted),1);
for i = 1:length(I_p_wave_areas_4_right_atrial_appendage_sorted)
    I_p_wave_areas_4_right_atrial_appendage_dev_avg_pls_1SD(i) = I_p_wave_areas_4_right_atrial_appendage_sorted(i) - (I_p_wave_areas_4_right_atrial_appendage_average+I_p_wave_areas_4_right_atrial_appendage_SD);
end
I_p_wave_areas_4_right_atrial_appendage_dev_avg_pls_1SD = abs(I_p_wave_areas_4_right_atrial_appendage_dev_avg_pls_1SD);
for i = 1:length(I_p_wave_areas_4_right_atrial_appendage_sorted)
    if I_p_wave_areas_4_right_atrial_appendage_dev_avg_pls_1SD(i) == min(I_p_wave_areas_4_right_atrial_appendage_dev_avg_pls_1SD)
        I_p_wave_areas_4_right_atrial_appendage_dev_avg_pls_1SD_index = i;
    end
end
I_p_wave_areas_4_right_atrial_appendage_dev_avg_min_2SD = area_4(length(I_p_wave_areas_4_right_atrial_appendage_sorted),1);
for i = 1:length(I_p_wave_areas_4_right_atrial_appendage_sorted)
    I_p_wave_areas_4_right_atrial_appendage_dev_avg_min_2SD(i) = I_p_wave_areas_4_right_atrial_appendage_sorted(i) - (I_p_wave_areas_4_right_atrial_appendage_average-2*I_p_wave_areas_4_right_atrial_appendage_SD);
end
I_p_wave_areas_4_right_atrial_appendage_dev_avg_min_2SD = abs(I_p_wave_areas_4_right_atrial_appendage_dev_avg_min_2SD);
for i = 1:length(I_p_wave_areas_4_right_atrial_appendage_sorted)
    if I_p_wave_areas_4_right_atrial_appendage_dev_avg_min_2SD(i) == min(I_p_wave_areas_4_right_atrial_appendage_dev_avg_min_2SD)
        I_p_wave_areas_4_dev_avg_min_2SD_index = i;
    end
end
I_p_wave_areas_4_right_atrial_appendage_dev_avg_pls_2SD = area_4(length(I_p_wave_areas_4_right_atrial_appendage_sorted),1);
for i = 1:length(I_p_wave_areas_4_right_atrial_appendage_sorted)
    I_p_wave_areas_4_right_atrial_appendage_dev_avg_pls_2SD(i) = I_p_wave_areas_4_right_atrial_appendage_sorted(i) - (I_p_wave_areas_4_right_atrial_appendage_average+2*I_p_wave_areas_4_right_atrial_appendage_SD);
end
I_p_wave_areas_4_right_atrial_appendage_dev_avg_pls_2SD = abs(I_p_wave_areas_4_right_atrial_appendage_dev_avg_pls_2SD);
for i = 1:length(I_p_wave_areas_4_right_atrial_appendage_sorted)
    if I_p_wave_areas_4_right_atrial_appendage_dev_avg_pls_2SD(i) == min(I_p_wave_areas_4_right_atrial_appendage_dev_avg_pls_2SD)
        I_p_wave_areas_4_right_atrial_appendage_dev_avg_pls_2SD_index = i;
    end
end
I_p_wave_areas_4_right_atrial_appendage_dev_avg_min_3SD = area_4(length(I_p_wave_areas_4_right_atrial_appendage_sorted),1);
for i = 1:length(I_p_wave_areas_4_right_atrial_appendage_sorted)
    I_p_wave_areas_4_right_atrial_appendage_dev_avg_min_3SD(i) = I_p_wave_areas_4_right_atrial_appendage_sorted(i) - (I_p_wave_areas_4_right_atrial_appendage_average-3*I_p_wave_areas_4_right_atrial_appendage_SD);
end
I_p_wave_areas_4_right_atrial_appendage_dev_avg_min_3SD = abs(I_p_wave_areas_4_right_atrial_appendage_dev_avg_min_3SD);
for i = 1:length(I_p_wave_areas_4_right_atrial_appendage_sorted)
    if I_p_wave_areas_4_right_atrial_appendage_dev_avg_min_3SD(i) == min(I_p_wave_areas_4_right_atrial_appendage_dev_avg_min_3SD)
        I_p_wave_areas_4_right_atrial_appendage_dev_avg_min_3SD_index = i;
    end
end
I_p_wave_areas_4_right_atrial_appendage_dev_avg_pls_3SD = area_4(length(I_p_wave_areas_4_right_atrial_appendage_sorted),1);
for i = 1:length(I_p_wave_areas_4_right_atrial_appendage_sorted)
    I_p_wave_areas_4_right_atrial_appendage_dev_avg_pls_3SD(i) = I_p_wave_areas_4_right_atrial_appendage_sorted(i) - (I_p_wave_areas_4_right_atrial_appendage_average+3*I_p_wave_areas_4_right_atrial_appendage_SD);
end
I_p_wave_areas_4_right_atrial_appendage_dev_avg_pls_3SD = abs(I_p_wave_areas_4_right_atrial_appendage_dev_avg_pls_3SD);
for i = 1:length(I_p_wave_areas_4_right_atrial_appendage_sorted)
    if I_p_wave_areas_4_right_atrial_appendage_dev_avg_pls_3SD(i) == min(I_p_wave_areas_4_right_atrial_appendage_dev_avg_pls_3SD)
        I_p_wave_areas_4_right_atrial_appendage_dev_avg_pls_3SD_index = i;
    end
end
I_p_wave_areas_4_right_atrial_appendage_tot_number_of_indices = length(I_p_wave_areas_4_right_atrial_appendage_sorted);
%This is the approximation of the tot "area under curve." 
I_p_wave_areas_4_right_atrial_appendage_quotient_within_1SD = (I_p_wave_areas_4_right_atrial_appendage_dev_avg_pls_1SD_index - I_p_wave_areas_4_right_atrial_appendagedev_avg_min_1SD_index)/I_p_wave_areas_4_right_atrial_appendage_tot_number_of_indices;
I_p_wave_areas_4_right_atrial_appendage_quotient_within_2SD = (I_p_wave_areas_4_right_atrial_appendage_dev_avg_pls_2SD_index - I_p_wave_areas_4_dev_avg_min_2SD_index)/I_p_wave_areas_4_right_atrial_appendage_tot_number_of_indices;
I_p_wave_areas_4_right_atrial_appendage_quotient_within_3SD = (I_p_wave_areas_4_right_atrial_appendage_dev_avg_pls_3SD_index - I_p_wave_areas_4_right_atrial_appendage_dev_avg_min_3SD_index)/I_p_wave_areas_4_right_atrial_appendage_tot_number_of_indices;
if (I_p_wave_areas_4_right_atrial_appendage_quotient_within_1SD > 0.66 && I_p_wave_areas_4_right_atrial_appendage_quotient_within_1SD < 0.70) && (I_p_wave_areas_4_right_atrial_appendage_quotient_within_2SD > 0.93 && I_p_wave_areas_4_right_atrial_appendage_quotient_within_2SD < 0.97) && (I_p_wave_areas_4_right_atrial_appendage_quotient_within_3SD > 0.98 && I_p_wave_areas_4_right_atrial_appendage_quotient_within_3SD < 1)
    isNormal(10,8) = 1;
end
if I_p_wave_areas_4_right_atrial_appendage_quotient_within_1SD == 0
    isSpike(10,8) = 1;
end
I_p_wave_areas_4_pulmonary_veins_sorted = sort(I_p_wave_areas_4_pulmonary_veins);
I_p_wave_areas_4_pulmonary_veins_dev_avg_min_1SD = area_4(length(I_p_wave_areas_4_pulmonary_veins_sorted),1);
for i = 1:length(I_p_wave_areas_4_pulmonary_veins_sorted)
    I_p_wave_areas_4_pulmonary_veins_dev_avg_min_1SD(i) = I_p_wave_areas_4_pulmonary_veins_sorted(i) - (I_p_wave_areas_4_pulmonary_veins_average-I_p_wave_areas_4_pulmonary_veins_SD);
end
I_p_wave_areas_4_pulmonary_veins_dev_avg_min_1SD = abs(I_p_wave_areas_4_pulmonary_veins_dev_avg_min_1SD);
for i = 1:length(I_p_wave_areas_4_pulmonary_veins_sorted)
    if I_p_wave_areas_4_pulmonary_veins_dev_avg_min_1SD(i) == min(I_p_wave_areas_4_pulmonary_veins_dev_avg_min_1SD)
        I_p_wave_areas_4_pulmonary_veins_dev_avg_min_1SD_index = i;
    end
end
I_p_wave_areas_4_pulmonary_veins_dev_avg_pls_1SD = area_4(length(I_p_wave_areas_4_pulmonary_veins_sorted),1);
for i = 1:length(I_p_wave_areas_4_pulmonary_veins_sorted)
    I_p_wave_areas_4_pulmonary_veins_dev_avg_pls_1SD(i) = I_p_wave_areas_4_pulmonary_veins_sorted(i) - (I_p_wave_areas_4_pulmonary_veins_average+I_p_wave_areas_4_pulmonary_veins_SD);
end
I_p_wave_areas_4_pulmonary_veins_dev_avg_pls_1SD = abs(I_p_wave_areas_4_pulmonary_veins_dev_avg_pls_1SD);
for i = 1:length(I_p_wave_areas_4_pulmonary_veins_sorted)
    if I_p_wave_areas_4_pulmonary_veins_dev_avg_pls_1SD(i) == min(I_p_wave_areas_4_pulmonary_veins_dev_avg_pls_1SD)
        I_p_wave_areas_4_pulmonary_veins_dev_avg_pls_1SD_index = i;
    end
end
I_p_wave_areas_4_pulmonary_veins_dev_avg_min_2SD = area_4(length(I_p_wave_areas_4_pulmonary_veins_sorted),1);
for i = 1:length(I_p_wave_areas_4_pulmonary_veins_sorted)
    I_p_wave_areas_4_pulmonary_veins_dev_avg_min_2SD(i) = I_p_wave_areas_4_pulmonary_veins_sorted(i) - (I_p_wave_areas_4_pulmonary_veins_average-2*I_p_wave_areas_4_pulmonary_veins_SD);
end
I_p_wave_areas_4_pulmonary_veins_dev_avg_min_2SD = abs(I_p_wave_areas_4_pulmonary_veins_dev_avg_min_2SD);
for i = 1:length(I_p_wave_areas_4_pulmonary_veins_sorted)
    if I_p_wave_areas_4_pulmonary_veins_dev_avg_min_2SD(i) == min(I_p_wave_areas_4_pulmonary_veins_dev_avg_min_2SD)
        I_p_wave_areas_4_dev_avg_min_2SD_index = i;
    end
end
I_p_wave_areas_4_pulmonary_veins_dev_avg_pls_2SD = area_4(length(I_p_wave_areas_4_pulmonary_veins_sorted),1);
for i = 1:length(I_p_wave_areas_4_pulmonary_veins_sorted)
    I_p_wave_areas_4_pulmonary_veins_dev_avg_pls_2SD(i) = I_p_wave_areas_4_pulmonary_veins_sorted(i) - (I_p_wave_areas_4_pulmonary_veins_average+2*I_p_wave_areas_4_pulmonary_veins_SD);
end
I_p_wave_areas_4_pulmonary_veins_dev_avg_pls_2SD = abs(I_p_wave_areas_4_pulmonary_veins_dev_avg_pls_2SD);
for i = 1:length(I_p_wave_areas_4_pulmonary_veins_sorted)
    if I_p_wave_areas_4_pulmonary_veins_dev_avg_pls_2SD(i) == min(I_p_wave_areas_4_pulmonary_veins_dev_avg_pls_2SD)
        I_p_wave_areas_4_pulmonary_veins_dev_avg_pls_2SD_index = i;
    end
end
I_p_wave_areas_4_pulmonary_veins_dev_avg_min_3SD = area_4(length(I_p_wave_areas_4_pulmonary_veins_sorted),1);
for i = 1:length(I_p_wave_areas_4_pulmonary_veins_sorted)
    I_p_wave_areas_4_pulmonary_veins_dev_avg_min_3SD(i) = I_p_wave_areas_4_pulmonary_veins_sorted(i) - (I_p_wave_areas_4_pulmonary_veins_average-3*I_p_wave_areas_4_pulmonary_veins_SD);
end
I_p_wave_areas_4_pulmonary_veins_dev_avg_min_3SD = abs(I_p_wave_areas_4_pulmonary_veins_dev_avg_min_3SD);
for i = 1:length(I_p_wave_areas_4_pulmonary_veins_sorted)
    if I_p_wave_areas_4_pulmonary_veins_dev_avg_min_3SD(i) == min(I_p_wave_areas_4_pulmonary_veins_dev_avg_min_3SD)
        I_p_wave_areas_4_pulmonary_veins_dev_avg_min_3SD_index = i;
    end
end
I_p_wave_areas_4_pulmonary_veins_dev_avg_pls_3SD = area_4(length(I_p_wave_areas_4_pulmonary_veins_sorted),1);
for i = 1:length(I_p_wave_areas_4_pulmonary_veins_sorted)
    I_p_wave_areas_4_pulmonary_veins_dev_avg_pls_3SD(i) = I_p_wave_areas_4_pulmonary_veins_sorted(i) - (I_p_wave_areas_4_pulmonary_veins_average+3*I_p_wave_areas_4_pulmonary_veins_SD);
end
I_p_wave_areas_4_pulmonary_veins_dev_avg_pls_3SD = abs(I_p_wave_areas_4_pulmonary_veins_dev_avg_pls_3SD);
for i = 1:length(I_p_wave_areas_4_pulmonary_veins_sorted)
    if I_p_wave_areas_4_pulmonary_veins_dev_avg_pls_3SD(i) == min(I_p_wave_areas_4_pulmonary_veins_dev_avg_pls_3SD)
        I_p_wave_areas_4_pulmonary_veins_dev_avg_pls_3SD_index = i;
    end
end
I_p_wave_areas_4_pulmonary_veins_tot_number_of_indices = length(I_p_wave_areas_4_pulmonary_veins_sorted);
%This is the approximation of the tot "area under curve." 
I_p_wave_areas_4_pulmonary_veins_quotient_within_1SD = (I_p_wave_areas_4_pulmonary_veins_dev_avg_pls_1SD_index - I_p_wave_areas_4_pulmonary_veinsdev_avg_min_1SD_index)/I_p_wave_areas_4_pulmonary_veins_tot_number_of_indices;
I_p_wave_areas_4_pulmonary_veins_quotient_within_2SD = (I_p_wave_areas_4_pulmonary_veins_dev_avg_pls_2SD_index - I_p_wave_areas_4_dev_avg_min_2SD_index)/I_p_wave_areas_4_pulmonary_veins_tot_number_of_indices;
I_p_wave_areas_4_pulmonary_veins_quotient_within_3SD = (I_p_wave_areas_4_pulmonary_veins_dev_avg_pls_3SD_index - I_p_wave_areas_4_pulmonary_veins_dev_avg_min_3SD_index)/I_p_wave_areas_4_pulmonary_veins_tot_number_of_indices;
if (I_p_wave_areas_4_pulmonary_veins_quotient_within_1SD > 0.66 && I_p_wave_areas_4_pulmonary_veins_quotient_within_1SD < 0.70) && (I_p_wave_areas_4_pulmonary_veins_quotient_within_2SD > 0.93 && I_p_wave_areas_4_pulmonary_veins_quotient_within_2SD < 0.97) && (I_p_wave_areas_4_pulmonary_veins_quotient_within_3SD > 0.98 && I_p_wave_areas_4_pulmonary_veins_quotient_within_3SD < 1)
    isNormal(10,9) = 1;
end
if I_p_wave_areas_4_pulmonary_veins_quotient_within_1SD == 0
    isSpike(10,9) = 1;
end
I_p_wave_areas_4_mitral_annulus_sorted = sort(I_p_wave_areas_4_mitral_annulus);
I_p_wave_areas_4_mitral_annulus_dev_avg_min_1SD = area_4(length(I_p_wave_areas_4_mitral_annulus_sorted),1);
for i = 1:length(I_p_wave_areas_4_mitral_annulus_sorted)
    I_p_wave_areas_4_mitral_annulus_dev_avg_min_1SD(i) = I_p_wave_areas_4_mitral_annulus_sorted(i) - (I_p_wave_areas_4_mitral_annulus_average-I_p_wave_areas_4_mitral_annulus_SD);
end
I_p_wave_areas_4_mitral_annulus_dev_avg_min_1SD = abs(I_p_wave_areas_4_mitral_annulus_dev_avg_min_1SD);
for i = 1:length(I_p_wave_areas_4_mitral_annulus_sorted)
    if I_p_wave_areas_4_mitral_annulus_dev_avg_min_1SD(i) == min(I_p_wave_areas_4_mitral_annulus_dev_avg_min_1SD)
        I_p_wave_areas_4_mitral_annulus_dev_avg_min_1SD_index = i;
    end
end
I_p_wave_areas_4_mitral_annulus_dev_avg_pls_1SD = area_4(length(I_p_wave_areas_4_mitral_annulus_sorted),1);
for i = 1:length(I_p_wave_areas_4_mitral_annulus_sorted)
    I_p_wave_areas_4_mitral_annulus_dev_avg_pls_1SD(i) = I_p_wave_areas_4_mitral_annulus_sorted(i) - (I_p_wave_areas_4_mitral_annulus_average+I_p_wave_areas_4_mitral_annulus_SD);
end
I_p_wave_areas_4_mitral_annulus_dev_avg_pls_1SD = abs(I_p_wave_areas_4_mitral_annulus_dev_avg_pls_1SD);
for i = 1:length(I_p_wave_areas_4_mitral_annulus_sorted)
    if I_p_wave_areas_4_mitral_annulus_dev_avg_pls_1SD(i) == min(I_p_wave_areas_4_mitral_annulus_dev_avg_pls_1SD)
        I_p_wave_areas_4_mitral_annulus_dev_avg_pls_1SD_index = i;
    end
end
I_p_wave_areas_4_mitral_annulus_dev_avg_min_2SD = area_4(length(I_p_wave_areas_4_mitral_annulus_sorted),1);
for i = 1:length(I_p_wave_areas_4_mitral_annulus_sorted)
    I_p_wave_areas_4_mitral_annulus_dev_avg_min_2SD(i) = I_p_wave_areas_4_mitral_annulus_sorted(i) - (I_p_wave_areas_4_mitral_annulus_average-2*I_p_wave_areas_4_mitral_annulus_SD);
end
I_p_wave_areas_4_mitral_annulus_dev_avg_min_2SD = abs(I_p_wave_areas_4_mitral_annulus_dev_avg_min_2SD);
for i = 1:length(I_p_wave_areas_4_mitral_annulus_sorted)
    if I_p_wave_areas_4_mitral_annulus_dev_avg_min_2SD(i) == min(I_p_wave_areas_4_mitral_annulus_dev_avg_min_2SD)
        I_p_wave_areas_4_dev_avg_min_2SD_index = i;
    end
end
I_p_wave_areas_4_mitral_annulus_dev_avg_pls_2SD = area_4(length(I_p_wave_areas_4_mitral_annulus_sorted),1);
for i = 1:length(I_p_wave_areas_4_mitral_annulus_sorted)
    I_p_wave_areas_4_mitral_annulus_dev_avg_pls_2SD(i) = I_p_wave_areas_4_mitral_annulus_sorted(i) - (I_p_wave_areas_4_mitral_annulus_average+2*I_p_wave_areas_4_mitral_annulus_SD);
end
I_p_wave_areas_4_mitral_annulus_dev_avg_pls_2SD = abs(I_p_wave_areas_4_mitral_annulus_dev_avg_pls_2SD);
for i = 1:length(I_p_wave_areas_4_mitral_annulus_sorted)
    if I_p_wave_areas_4_mitral_annulus_dev_avg_pls_2SD(i) == min(I_p_wave_areas_4_mitral_annulus_dev_avg_pls_2SD)
        I_p_wave_areas_4_mitral_annulus_dev_avg_pls_2SD_index = i;
    end
end
I_p_wave_areas_4_mitral_annulus_dev_avg_min_3SD = area_4(length(I_p_wave_areas_4_mitral_annulus_sorted),1);
for i = 1:length(I_p_wave_areas_4_mitral_annulus_sorted)
    I_p_wave_areas_4_mitral_annulus_dev_avg_min_3SD(i) = I_p_wave_areas_4_mitral_annulus_sorted(i) - (I_p_wave_areas_4_mitral_annulus_average-3*I_p_wave_areas_4_mitral_annulus_SD);
end
I_p_wave_areas_4_mitral_annulus_dev_avg_min_3SD = abs(I_p_wave_areas_4_mitral_annulus_dev_avg_min_3SD);
for i = 1:length(I_p_wave_areas_4_mitral_annulus_sorted)
    if I_p_wave_areas_4_mitral_annulus_dev_avg_min_3SD(i) == min(I_p_wave_areas_4_mitral_annulus_dev_avg_min_3SD)
        I_p_wave_areas_4_mitral_annulus_dev_avg_min_3SD_index = i;
    end
end
I_p_wave_areas_4_mitral_annulus_dev_avg_pls_3SD = area_4(length(I_p_wave_areas_4_mitral_annulus_sorted),1);
for i = 1:length(I_p_wave_areas_4_mitral_annulus_sorted)
    I_p_wave_areas_4_mitral_annulus_dev_avg_pls_3SD(i) = I_p_wave_areas_4_mitral_annulus_sorted(i) - (I_p_wave_areas_4_mitral_annulus_average+3*I_p_wave_areas_4_mitral_annulus_SD);
end
I_p_wave_areas_4_mitral_annulus_dev_avg_pls_3SD = abs(I_p_wave_areas_4_mitral_annulus_dev_avg_pls_3SD);
for i = 1:length(I_p_wave_areas_4_mitral_annulus_sorted)
    if I_p_wave_areas_4_mitral_annulus_dev_avg_pls_3SD(i) == min(I_p_wave_areas_4_mitral_annulus_dev_avg_pls_3SD)
        I_p_wave_areas_4_mitral_annulus_dev_avg_pls_3SD_index = i;
    end
end
I_p_wave_areas_4_mitral_annulus_tot_number_of_indices = length(I_p_wave_areas_4_mitral_annulus_sorted);
%This is the approximation of the tot "area under curve." 
I_p_wave_areas_4_mitral_annulus_quotient_within_1SD = (I_p_wave_areas_4_mitral_annulus_dev_avg_pls_1SD_index - I_p_wave_areas_4_mitral_annulusdev_avg_min_1SD_index)/I_p_wave_areas_4_mitral_annulus_tot_number_of_indices;
I_p_wave_areas_4_mitral_annulus_quotient_within_2SD = (I_p_wave_areas_4_mitral_annulus_dev_avg_pls_2SD_index - I_p_wave_areas_4_dev_avg_min_2SD_index)/I_p_wave_areas_4_mitral_annulus_tot_number_of_indices;
I_p_wave_areas_4_mitral_annulus_quotient_within_3SD = (I_p_wave_areas_4_mitral_annulus_dev_avg_pls_3SD_index - I_p_wave_areas_4_mitral_annulus_dev_avg_min_3SD_index)/I_p_wave_areas_4_mitral_annulus_tot_number_of_indices;
if (I_p_wave_areas_4_mitral_annulus_quotient_within_1SD > 0.66 && I_p_wave_areas_4_mitral_annulus_quotient_within_1SD < 0.70) && (I_p_wave_areas_4_mitral_annulus_quotient_within_2SD > 0.93 && I_p_wave_areas_4_mitral_annulus_quotient_within_2SD < 0.97) && (I_p_wave_areas_4_mitral_annulus_quotient_within_3SD > 0.98 && I_p_wave_areas_4_mitral_annulus_quotient_within_3SD < 1)
    isNormal(10,10) = 1;
end
if I_p_wave_areas_4_mitral_annulus_quotient_within_1SD == 0
    isSpike(10,10) = 1;
end
I_p_wave_areas_4_CS_body_sorted = sort(I_p_wave_areas_4_CS_body);
I_p_wave_areas_4_CS_body_dev_avg_min_1SD = area_4(length(I_p_wave_areas_4_CS_body_sorted),1);
for i = 1:length(I_p_wave_areas_4_CS_body_sorted)
    I_p_wave_areas_4_CS_body_dev_avg_min_1SD(i) = I_p_wave_areas_4_CS_body_sorted(i) - (I_p_wave_areas_4_CS_body_average-I_p_wave_areas_4_CS_body_SD);
end
I_p_wave_areas_4_CS_body_dev_avg_min_1SD = abs(I_p_wave_areas_4_CS_body_dev_avg_min_1SD);
for i = 1:length(I_p_wave_areas_4_CS_body_sorted)
    if I_p_wave_areas_4_CS_body_dev_avg_min_1SD(i) == min(I_p_wave_areas_4_CS_body_dev_avg_min_1SD)
        I_p_wave_areas_4_CS_body_dev_avg_min_1SD_index = i;
    end
end
I_p_wave_areas_4_CS_body_dev_avg_pls_1SD = area_4(length(I_p_wave_areas_4_CS_body_sorted),1);
for i = 1:length(I_p_wave_areas_4_CS_body_sorted)
    I_p_wave_areas_4_CS_body_dev_avg_pls_1SD(i) = I_p_wave_areas_4_CS_body_sorted(i) - (I_p_wave_areas_4_CS_body_average+I_p_wave_areas_4_CS_body_SD);
end
I_p_wave_areas_4_CS_body_dev_avg_pls_1SD = abs(I_p_wave_areas_4_CS_body_dev_avg_pls_1SD);
for i = 1:length(I_p_wave_areas_4_CS_body_sorted)
    if I_p_wave_areas_4_CS_body_dev_avg_pls_1SD(i) == min(I_p_wave_areas_4_CS_body_dev_avg_pls_1SD)
        I_p_wave_areas_4_CS_body_dev_avg_pls_1SD_index = i;
    end
end
I_p_wave_areas_4_CS_body_dev_avg_min_2SD = area_4(length(I_p_wave_areas_4_CS_body_sorted),1);
for i = 1:length(I_p_wave_areas_4_CS_body_sorted)
    I_p_wave_areas_4_CS_body_dev_avg_min_2SD(i) = I_p_wave_areas_4_CS_body_sorted(i) - (I_p_wave_areas_4_CS_body_average-2*I_p_wave_areas_4_CS_body_SD);
end
I_p_wave_areas_4_CS_body_dev_avg_min_2SD = abs(I_p_wave_areas_4_CS_body_dev_avg_min_2SD);
for i = 1:length(I_p_wave_areas_4_CS_body_sorted)
    if I_p_wave_areas_4_CS_body_dev_avg_min_2SD(i) == min(I_p_wave_areas_4_CS_body_dev_avg_min_2SD)
        I_p_wave_areas_4_dev_avg_min_2SD_index = i;
    end
end
I_p_wave_areas_4_CS_body_dev_avg_pls_2SD = area_4(length(I_p_wave_areas_4_CS_body_sorted),1);
for i = 1:length(I_p_wave_areas_4_CS_body_sorted)
    I_p_wave_areas_4_CS_body_dev_avg_pls_2SD(i) = I_p_wave_areas_4_CS_body_sorted(i) - (I_p_wave_areas_4_CS_body_average+2*I_p_wave_areas_4_CS_body_SD);
end
I_p_wave_areas_4_CS_body_dev_avg_pls_2SD = abs(I_p_wave_areas_4_CS_body_dev_avg_pls_2SD);
for i = 1:length(I_p_wave_areas_4_CS_body_sorted)
    if I_p_wave_areas_4_CS_body_dev_avg_pls_2SD(i) == min(I_p_wave_areas_4_CS_body_dev_avg_pls_2SD)
        I_p_wave_areas_4_CS_body_dev_avg_pls_2SD_index = i;
    end
end
I_p_wave_areas_4_CS_body_dev_avg_min_3SD = area_4(length(I_p_wave_areas_4_CS_body_sorted),1);
for i = 1:length(I_p_wave_areas_4_CS_body_sorted)
    I_p_wave_areas_4_CS_body_dev_avg_min_3SD(i) = I_p_wave_areas_4_CS_body_sorted(i) - (I_p_wave_areas_4_CS_body_average-3*I_p_wave_areas_4_CS_body_SD);
end
I_p_wave_areas_4_CS_body_dev_avg_min_3SD = abs(I_p_wave_areas_4_CS_body_dev_avg_min_3SD);
for i = 1:length(I_p_wave_areas_4_CS_body_sorted)
    if I_p_wave_areas_4_CS_body_dev_avg_min_3SD(i) == min(I_p_wave_areas_4_CS_body_dev_avg_min_3SD)
        I_p_wave_areas_4_CS_body_dev_avg_min_3SD_index = i;
    end
end
I_p_wave_areas_4_CS_body_dev_avg_pls_3SD = area_4(length(I_p_wave_areas_4_CS_body_sorted),1);
for i = 1:length(I_p_wave_areas_4_CS_body_sorted)
    I_p_wave_areas_4_CS_body_dev_avg_pls_3SD(i) = I_p_wave_areas_4_CS_body_sorted(i) - (I_p_wave_areas_4_CS_body_average+3*I_p_wave_areas_4_CS_body_SD);
end
I_p_wave_areas_4_CS_body_dev_avg_pls_3SD = abs(I_p_wave_areas_4_CS_body_dev_avg_pls_3SD);
for i = 1:length(I_p_wave_areas_4_CS_body_sorted)
    if I_p_wave_areas_4_CS_body_dev_avg_pls_3SD(i) == min(I_p_wave_areas_4_CS_body_dev_avg_pls_3SD)
        I_p_wave_areas_4_CS_body_dev_avg_pls_3SD_index = i;
    end
end
I_p_wave_areas_4_CS_body_tot_number_of_indices = length(I_p_wave_areas_4_CS_body_sorted);
%This is the approximation of the tot "area under curve." 
I_p_wave_areas_4_CS_body_quotient_within_1SD = (I_p_wave_areas_4_CS_body_dev_avg_pls_1SD_index - I_p_wave_areas_4_CS_bodydev_avg_min_1SD_index)/I_p_wave_areas_4_CS_body_tot_number_of_indices;
I_p_wave_areas_4_CS_body_quotient_within_2SD = (I_p_wave_areas_4_CS_body_dev_avg_pls_2SD_index - I_p_wave_areas_4_dev_avg_min_2SD_index)/I_p_wave_areas_4_CS_body_tot_number_of_indices;
I_p_wave_areas_4_CS_body_quotient_within_3SD = (I_p_wave_areas_4_CS_body_dev_avg_pls_3SD_index - I_p_wave_areas_4_CS_body_dev_avg_min_3SD_index)/I_p_wave_areas_4_CS_body_tot_number_of_indices;
if (I_p_wave_areas_4_CS_body_quotient_within_1SD > 0.66 && I_p_wave_areas_4_CS_body_quotient_within_1SD < 0.70) && (I_p_wave_areas_4_CS_body_quotient_within_2SD > 0.93 && I_p_wave_areas_4_CS_body_quotient_within_2SD < 0.97) && (I_p_wave_areas_4_CS_body_quotient_within_3SD > 0.98 && I_p_wave_areas_4_CS_body_quotient_within_3SD < 1)
    isNormal(10,11) = 1;
end
if I_p_wave_areas_4_CS_body_quotient_within_1SD == 0
    isSpike(10,11) = 1;
end
I_p_wave_areas_4_left_septum_sorted = sort(I_p_wave_areas_4_left_septum);
I_p_wave_areas_4_left_septum_dev_avg_min_1SD = area_4(length(I_p_wave_areas_4_left_septum_sorted),1);
for i = 1:length(I_p_wave_areas_4_left_septum_sorted)
    I_p_wave_areas_4_left_septum_dev_avg_min_1SD(i) = I_p_wave_areas_4_left_septum_sorted(i) - (I_p_wave_areas_4_left_septum_average-I_p_wave_areas_4_left_septum_SD);
end
I_p_wave_areas_4_left_septum_dev_avg_min_1SD = abs(I_p_wave_areas_4_left_septum_dev_avg_min_1SD);
for i = 1:length(I_p_wave_areas_4_left_septum_sorted)
    if I_p_wave_areas_4_left_septum_dev_avg_min_1SD(i) == min(I_p_wave_areas_4_left_septum_dev_avg_min_1SD)
        I_p_wave_areas_4_left_septum_dev_avg_min_1SD_index = i;
    end
end
I_p_wave_areas_4_left_septum_dev_avg_pls_1SD = area_4(length(I_p_wave_areas_4_left_septum_sorted),1);
for i = 1:length(I_p_wave_areas_4_left_septum_sorted)
    I_p_wave_areas_4_left_septum_dev_avg_pls_1SD(i) = I_p_wave_areas_4_left_septum_sorted(i) - (I_p_wave_areas_4_left_septum_average+I_p_wave_areas_4_left_septum_SD);
end
I_p_wave_areas_4_left_septum_dev_avg_pls_1SD = abs(I_p_wave_areas_4_left_septum_dev_avg_pls_1SD);
for i = 1:length(I_p_wave_areas_4_left_septum_sorted)
    if I_p_wave_areas_4_left_septum_dev_avg_pls_1SD(i) == min(I_p_wave_areas_4_left_septum_dev_avg_pls_1SD)
        I_p_wave_areas_4_left_septum_dev_avg_pls_1SD_index = i;
    end
end
I_p_wave_areas_4_left_septum_dev_avg_min_2SD = area_4(length(I_p_wave_areas_4_left_septum_sorted),1);
for i = 1:length(I_p_wave_areas_4_left_septum_sorted)
    I_p_wave_areas_4_left_septum_dev_avg_min_2SD(i) = I_p_wave_areas_4_left_septum_sorted(i) - (I_p_wave_areas_4_left_septum_average-2*I_p_wave_areas_4_left_septum_SD);
end
I_p_wave_areas_4_left_septum_dev_avg_min_2SD = abs(I_p_wave_areas_4_left_septum_dev_avg_min_2SD);
for i = 1:length(I_p_wave_areas_4_left_septum_sorted)
    if I_p_wave_areas_4_left_septum_dev_avg_min_2SD(i) == min(I_p_wave_areas_4_left_septum_dev_avg_min_2SD)
        I_p_wave_areas_4_dev_avg_min_2SD_index = i;
    end
end
I_p_wave_areas_4_left_septum_dev_avg_pls_2SD = area_4(length(I_p_wave_areas_4_left_septum_sorted),1);
for i = 1:length(I_p_wave_areas_4_left_septum_sorted)
    I_p_wave_areas_4_left_septum_dev_avg_pls_2SD(i) = I_p_wave_areas_4_left_septum_sorted(i) - (I_p_wave_areas_4_left_septum_average+2*I_p_wave_areas_4_left_septum_SD);
end
I_p_wave_areas_4_left_septum_dev_avg_pls_2SD = abs(I_p_wave_areas_4_left_septum_dev_avg_pls_2SD);
for i = 1:length(I_p_wave_areas_4_left_septum_sorted)
    if I_p_wave_areas_4_left_septum_dev_avg_pls_2SD(i) == min(I_p_wave_areas_4_left_septum_dev_avg_pls_2SD)
        I_p_wave_areas_4_left_septum_dev_avg_pls_2SD_index = i;
    end
end
I_p_wave_areas_4_left_septum_dev_avg_min_3SD = area_4(length(I_p_wave_areas_4_left_septum_sorted),1);
for i = 1:length(I_p_wave_areas_4_left_septum_sorted)
    I_p_wave_areas_4_left_septum_dev_avg_min_3SD(i) = I_p_wave_areas_4_left_septum_sorted(i) - (I_p_wave_areas_4_left_septum_average-3*I_p_wave_areas_4_left_septum_SD);
end
I_p_wave_areas_4_left_septum_dev_avg_min_3SD = abs(I_p_wave_areas_4_left_septum_dev_avg_min_3SD);
for i = 1:length(I_p_wave_areas_4_left_septum_sorted)
    if I_p_wave_areas_4_left_septum_dev_avg_min_3SD(i) == min(I_p_wave_areas_4_left_septum_dev_avg_min_3SD)
        I_p_wave_areas_4_left_septum_dev_avg_min_3SD_index = i;
    end
end
I_p_wave_areas_4_left_septum_dev_avg_pls_3SD = area_4(length(I_p_wave_areas_4_left_septum_sorted),1);
for i = 1:length(I_p_wave_areas_4_left_septum_sorted)
    I_p_wave_areas_4_left_septum_dev_avg_pls_3SD(i) = I_p_wave_areas_4_left_septum_sorted(i) - (I_p_wave_areas_4_left_septum_average+3*I_p_wave_areas_4_left_septum_SD);
end
I_p_wave_areas_4_left_septum_dev_avg_pls_3SD = abs(I_p_wave_areas_4_left_septum_dev_avg_pls_3SD);
for i = 1:length(I_p_wave_areas_4_left_septum_sorted)
    if I_p_wave_areas_4_left_septum_dev_avg_pls_3SD(i) == min(I_p_wave_areas_4_left_septum_dev_avg_pls_3SD)
        I_p_wave_areas_4_left_septum_dev_avg_pls_3SD_index = i;
    end
end
I_p_wave_areas_4_left_septum_tot_number_of_indices = length(I_p_wave_areas_4_left_septum_sorted);
%This is the approximation of the tot "area under curve." 
I_p_wave_areas_4_left_septum_quotient_within_1SD = (I_p_wave_areas_4_left_septum_dev_avg_pls_1SD_index - I_p_wave_areas_4_left_septumdev_avg_min_1SD_index)/I_p_wave_areas_4_left_septum_tot_number_of_indices;
I_p_wave_areas_4_left_septum_quotient_within_2SD = (I_p_wave_areas_4_left_septum_dev_avg_pls_2SD_index - I_p_wave_areas_4_dev_avg_min_2SD_index)/I_p_wave_areas_4_left_septum_tot_number_of_indices;
I_p_wave_areas_4_left_septum_quotient_within_3SD = (I_p_wave_areas_4_left_septum_dev_avg_pls_3SD_index - I_p_wave_areas_4_left_septum_dev_avg_min_3SD_index)/I_p_wave_areas_4_left_septum_tot_number_of_indices;
if (I_p_wave_areas_4_left_septum_quotient_within_1SD > 0.66 && I_p_wave_areas_4_left_septum_quotient_within_1SD < 0.70) && (I_p_wave_areas_4_left_septum_quotient_within_2SD > 0.93 && I_p_wave_areas_4_left_septum_quotient_within_2SD < 0.97) && (I_p_wave_areas_4_left_septum_quotient_within_3SD > 0.98 && I_p_wave_areas_4_left_septum_quotient_within_3SD < 1)
    isNormal(10,12) = 1;
end
if I_p_wave_areas_4_left_septum_quotient_within_1SD == 0
    isSpike(10,12) = 1;
end
I_p_wave_areas_4_left_atrial_appendage_sorted = sort(I_p_wave_areas_4_left_atrial_appendage);
I_p_wave_areas_4_left_atrial_appendage_dev_avg_min_1SD = area_4(length(I_p_wave_areas_4_left_atrial_appendage_sorted),1);
for i = 1:length(I_p_wave_areas_4_left_atrial_appendage_sorted)
    I_p_wave_areas_4_left_atrial_appendage_dev_avg_min_1SD(i) = I_p_wave_areas_4_left_atrial_appendage_sorted(i) - (I_p_wave_areas_4_left_atrial_appendage_average-I_p_wave_areas_4_left_atrial_appendage_SD);
end
I_p_wave_areas_4_left_atrial_appendage_dev_avg_min_1SD = abs(I_p_wave_areas_4_left_atrial_appendage_dev_avg_min_1SD);
for i = 1:length(I_p_wave_areas_4_left_atrial_appendage_sorted)
    if I_p_wave_areas_4_left_atrial_appendage_dev_avg_min_1SD(i) == min(I_p_wave_areas_4_left_atrial_appendage_dev_avg_min_1SD)
        I_p_wave_areas_4_left_atrial_appendage_dev_avg_min_1SD_index = i;
    end
end
I_p_wave_areas_4_left_atrial_appendage_dev_avg_pls_1SD = area_4(length(I_p_wave_areas_4_left_atrial_appendage_sorted),1);
for i = 1:length(I_p_wave_areas_4_left_atrial_appendage_sorted)
    I_p_wave_areas_4_left_atrial_appendage_dev_avg_pls_1SD(i) = I_p_wave_areas_4_left_atrial_appendage_sorted(i) - (I_p_wave_areas_4_left_atrial_appendage_average+I_p_wave_areas_4_left_atrial_appendage_SD);
end
I_p_wave_areas_4_left_atrial_appendage_dev_avg_pls_1SD = abs(I_p_wave_areas_4_left_atrial_appendage_dev_avg_pls_1SD);
for i = 1:length(I_p_wave_areas_4_left_atrial_appendage_sorted)
    if I_p_wave_areas_4_left_atrial_appendage_dev_avg_pls_1SD(i) == min(I_p_wave_areas_4_left_atrial_appendage_dev_avg_pls_1SD)
        I_p_wave_areas_4_left_atrial_appendage_dev_avg_pls_1SD_index = i;
    end
end
I_p_wave_areas_4_left_atrial_appendage_dev_avg_min_2SD = area_4(length(I_p_wave_areas_4_left_atrial_appendage_sorted),1);
for i = 1:length(I_p_wave_areas_4_left_atrial_appendage_sorted)
    I_p_wave_areas_4_left_atrial_appendage_dev_avg_min_2SD(i) = I_p_wave_areas_4_left_atrial_appendage_sorted(i) - (I_p_wave_areas_4_left_atrial_appendage_average-2*I_p_wave_areas_4_left_atrial_appendage_SD);
end
I_p_wave_areas_4_left_atrial_appendage_dev_avg_min_2SD = abs(I_p_wave_areas_4_left_atrial_appendage_dev_avg_min_2SD);
for i = 1:length(I_p_wave_areas_4_left_atrial_appendage_sorted)
    if I_p_wave_areas_4_left_atrial_appendage_dev_avg_min_2SD(i) == min(I_p_wave_areas_4_left_atrial_appendage_dev_avg_min_2SD)
        I_p_wave_areas_4_dev_avg_min_2SD_index = i;
    end
end
I_p_wave_areas_4_left_atrial_appendage_dev_avg_pls_2SD = area_4(length(I_p_wave_areas_4_left_atrial_appendage_sorted),1);
for i = 1:length(I_p_wave_areas_4_left_atrial_appendage_sorted)
    I_p_wave_areas_4_left_atrial_appendage_dev_avg_pls_2SD(i) = I_p_wave_areas_4_left_atrial_appendage_sorted(i) - (I_p_wave_areas_4_left_atrial_appendage_average+2*I_p_wave_areas_4_left_atrial_appendage_SD);
end
I_p_wave_areas_4_left_atrial_appendage_dev_avg_pls_2SD = abs(I_p_wave_areas_4_left_atrial_appendage_dev_avg_pls_2SD);
for i = 1:length(I_p_wave_areas_4_left_atrial_appendage_sorted)
    if I_p_wave_areas_4_left_atrial_appendage_dev_avg_pls_2SD(i) == min(I_p_wave_areas_4_left_atrial_appendage_dev_avg_pls_2SD)
        I_p_wave_areas_4_left_atrial_appendage_dev_avg_pls_2SD_index = i;
    end
end
I_p_wave_areas_4_left_atrial_appendage_dev_avg_min_3SD = area_4(length(I_p_wave_areas_4_left_atrial_appendage_sorted),1);
for i = 1:length(I_p_wave_areas_4_left_atrial_appendage_sorted)
    I_p_wave_areas_4_left_atrial_appendage_dev_avg_min_3SD(i) = I_p_wave_areas_4_left_atrial_appendage_sorted(i) - (I_p_wave_areas_4_left_atrial_appendage_average-3*I_p_wave_areas_4_left_atrial_appendage_SD);
end
I_p_wave_areas_4_left_atrial_appendage_dev_avg_min_3SD = abs(I_p_wave_areas_4_left_atrial_appendage_dev_avg_min_3SD);
for i = 1:length(I_p_wave_areas_4_left_atrial_appendage_sorted)
    if I_p_wave_areas_4_left_atrial_appendage_dev_avg_min_3SD(i) == min(I_p_wave_areas_4_left_atrial_appendage_dev_avg_min_3SD)
        I_p_wave_areas_4_left_atrial_appendage_dev_avg_min_3SD_index = i;
    end
end
I_p_wave_areas_4_left_atrial_appendage_dev_avg_pls_3SD = area_4(length(I_p_wave_areas_4_left_atrial_appendage_sorted),1);
for i = 1:length(I_p_wave_areas_4_left_atrial_appendage_sorted)
    I_p_wave_areas_4_left_atrial_appendage_dev_avg_pls_3SD(i) = I_p_wave_areas_4_left_atrial_appendage_sorted(i) - (I_p_wave_areas_4_left_atrial_appendage_average+3*I_p_wave_areas_4_left_atrial_appendage_SD);
end
I_p_wave_areas_4_left_atrial_appendage_dev_avg_pls_3SD = abs(I_p_wave_areas_4_left_atrial_appendage_dev_avg_pls_3SD);
for i = 1:length(I_p_wave_areas_4_left_atrial_appendage_sorted)
    if I_p_wave_areas_4_left_atrial_appendage_dev_avg_pls_3SD(i) == min(I_p_wave_areas_4_left_atrial_appendage_dev_avg_pls_3SD)
        I_p_wave_areas_4_left_atrial_appendage_dev_avg_pls_3SD_index = i;
    end
end
I_p_wave_areas_4_left_atrial_appendage_tot_number_of_indices = length(I_p_wave_areas_4_left_atrial_appendage_sorted);
%This is the approximation of the tot "area under curve." 
I_p_wave_areas_4_left_atrial_appendage_quotient_within_1SD = (I_p_wave_areas_4_left_atrial_appendage_dev_avg_pls_1SD_index - I_p_wave_areas_4_left_atrial_appendagedev_avg_min_1SD_index)/I_p_wave_areas_4_left_atrial_appendage_tot_number_of_indices;
I_p_wave_areas_4_left_atrial_appendage_quotient_within_2SD = (I_p_wave_areas_4_left_atrial_appendage_dev_avg_pls_2SD_index - I_p_wave_areas_4_dev_avg_min_2SD_index)/I_p_wave_areas_4_left_atrial_appendage_tot_number_of_indices;
I_p_wave_areas_4_left_atrial_appendage_quotient_within_3SD = (I_p_wave_areas_4_left_atrial_appendage_dev_avg_pls_3SD_index - I_p_wave_areas_4_left_atrial_appendage_dev_avg_min_3SD_index)/I_p_wave_areas_4_left_atrial_appendage_tot_number_of_indices;
if (I_p_wave_areas_4_left_atrial_appendage_quotient_within_1SD > 0.66 && I_p_wave_areas_4_left_atrial_appendage_quotient_within_1SD < 0.70) && (I_p_wave_areas_4_left_atrial_appendage_quotient_within_2SD > 0.93 && I_p_wave_areas_4_left_atrial_appendage_quotient_within_2SD < 0.97) && (I_p_wave_areas_4_left_atrial_appendage_quotient_within_3SD > 0.98 && I_p_wave_areas_4_left_atrial_appendage_quotient_within_3SD < 1)
    isNormal(10,13) = 1;
end
if I_p_wave_areas_4_left_atrial_appendage_quotient_within_1SD == 0
    isSpike(10,13) = 1;
end
I_p_wave_areas_5_SA_node_sorted = sort(I_p_wave_areas_5_SA_node);
I_p_wave_areas_5_SA_node_dev_avg_min_1SD = area_5(length(I_p_wave_areas_5_SA_node_sorted),1);
for i = 1:length(I_p_wave_areas_5_SA_node_sorted)
    I_p_wave_areas_5_SA_node_dev_avg_min_1SD(i) = I_p_wave_areas_5_SA_node_sorted(i) - (I_p_wave_areas_5_SA_node_average-I_p_wave_areas_5_SA_node_SD);
end
I_p_wave_areas_5_SA_node_dev_avg_min_1SD = abs(I_p_wave_areas_5_SA_node_dev_avg_min_1SD);
for i = 1:length(I_p_wave_areas_5_SA_node_sorted)
    if I_p_wave_areas_5_SA_node_dev_avg_min_1SD(i) == min(I_p_wave_areas_5_SA_node_dev_avg_min_1SD)
        I_p_wave_areas_5_SA_node_dev_avg_min_1SD_index = i;
    end
end
I_p_wave_areas_5_SA_node_dev_avg_pls_1SD = area_5(length(I_p_wave_areas_5_SA_node_sorted),1);
for i = 1:length(I_p_wave_areas_5_SA_node_sorted)
    I_p_wave_areas_5_SA_node_dev_avg_pls_1SD(i) = I_p_wave_areas_5_SA_node_sorted(i) - (I_p_wave_areas_5_SA_node_average+I_p_wave_areas_5_SA_node_SD);
end
I_p_wave_areas_5_SA_node_dev_avg_pls_1SD = abs(I_p_wave_areas_5_SA_node_dev_avg_pls_1SD);
for i = 1:length(I_p_wave_areas_5_SA_node_sorted)
    if I_p_wave_areas_5_SA_node_dev_avg_pls_1SD(i) == min(I_p_wave_areas_5_SA_node_dev_avg_pls_1SD)
        I_p_wave_areas_5_SA_node_dev_avg_pls_1SD_index = i;
    end
end
I_p_wave_areas_5_SA_node_dev_avg_min_2SD = area_5(length(I_p_wave_areas_5_SA_node_sorted),1);
for i = 1:length(I_p_wave_areas_5_SA_node_sorted)
    I_p_wave_areas_5_SA_node_dev_avg_min_2SD(i) = I_p_wave_areas_5_SA_node_sorted(i) - (I_p_wave_areas_5_SA_node_average-2*I_p_wave_areas_5_SA_node_SD);
end
I_p_wave_areas_5_SA_node_dev_avg_min_2SD = abs(I_p_wave_areas_5_SA_node_dev_avg_min_2SD);
for i = 1:length(I_p_wave_areas_5_SA_node_sorted)
    if I_p_wave_areas_5_SA_node_dev_avg_min_2SD(i) == min(I_p_wave_areas_5_SA_node_dev_avg_min_2SD)
        I_p_wave_areas_5_dev_avg_min_2SD_index = i;
    end
end
I_p_wave_areas_5_SA_node_dev_avg_pls_2SD = area_5(length(I_p_wave_areas_5_SA_node_sorted),1);
for i = 1:length(I_p_wave_areas_5_SA_node_sorted)
    I_p_wave_areas_5_SA_node_dev_avg_pls_2SD(i) = I_p_wave_areas_5_SA_node_sorted(i) - (I_p_wave_areas_5_SA_node_average+2*I_p_wave_areas_5_SA_node_SD);
end
I_p_wave_areas_5_SA_node_dev_avg_pls_2SD = abs(I_p_wave_areas_5_SA_node_dev_avg_pls_2SD);
for i = 1:length(I_p_wave_areas_5_SA_node_sorted)
    if I_p_wave_areas_5_SA_node_dev_avg_pls_2SD(i) == min(I_p_wave_areas_5_SA_node_dev_avg_pls_2SD)
        I_p_wave_areas_5_SA_node_dev_avg_pls_2SD_index = i;
    end
end
I_p_wave_areas_5_SA_node_dev_avg_min_3SD = area_5(length(I_p_wave_areas_5_SA_node_sorted),1);
for i = 1:length(I_p_wave_areas_5_SA_node_sorted)
    I_p_wave_areas_5_SA_node_dev_avg_min_3SD(i) = I_p_wave_areas_5_SA_node_sorted(i) - (I_p_wave_areas_5_SA_node_average-3*I_p_wave_areas_5_SA_node_SD);
end
I_p_wave_areas_5_SA_node_dev_avg_min_3SD = abs(I_p_wave_areas_5_SA_node_dev_avg_min_3SD);
for i = 1:length(I_p_wave_areas_5_SA_node_sorted)
    if I_p_wave_areas_5_SA_node_dev_avg_min_3SD(i) == min(I_p_wave_areas_5_SA_node_dev_avg_min_3SD)
        I_p_wave_areas_5_SA_node_dev_avg_min_3SD_index = i;
    end
end
I_p_wave_areas_5_SA_node_dev_avg_pls_3SD = area_5(length(I_p_wave_areas_5_SA_node_sorted),1);
for i = 1:length(I_p_wave_areas_5_SA_node_sorted)
    I_p_wave_areas_5_SA_node_dev_avg_pls_3SD(i) = I_p_wave_areas_5_SA_node_sorted(i) - (I_p_wave_areas_5_SA_node_average+3*I_p_wave_areas_5_SA_node_SD);
end
I_p_wave_areas_5_SA_node_dev_avg_pls_3SD = abs(I_p_wave_areas_5_SA_node_dev_avg_pls_3SD);
for i = 1:length(I_p_wave_areas_5_SA_node_sorted)
    if I_p_wave_areas_5_SA_node_dev_avg_pls_3SD(i) == min(I_p_wave_areas_5_SA_node_dev_avg_pls_3SD)
        I_p_wave_areas_5_SA_node_dev_avg_pls_3SD_index = i;
    end
end
I_p_wave_areas_5_SA_node_tot_number_of_indices = length(I_p_wave_areas_5_SA_node_sorted);
%This is the approximation of the tot "area under curve." 
I_p_wave_areas_5_SA_node_quotient_within_1SD = (I_p_wave_areas_5_SA_node_dev_avg_pls_1SD_index - I_p_wave_areas_5_SA_nodedev_avg_min_1SD_index)/I_p_wave_areas_5_SA_node_tot_number_of_indices;
I_p_wave_areas_5_SA_node_quotient_within_2SD = (I_p_wave_areas_5_SA_node_dev_avg_pls_2SD_index - I_p_wave_areas_5_dev_avg_min_2SD_index)/I_p_wave_areas_5_SA_node_tot_number_of_indices;
I_p_wave_areas_5_SA_node_quotient_within_3SD = (I_p_wave_areas_5_SA_node_dev_avg_pls_3SD_index - I_p_wave_areas_5_SA_node_dev_avg_min_3SD_index)/I_p_wave_areas_5_SA_node_tot_number_of_indices;
if ((I_p_wave_areas_5_SA_node_quotient_within_1SD > 0.66 && I_p_wave_areas_5_SA_node_quotient_within_1SD < 0.70) && (I_p_wave_areas_5_SA_node_quotient_within_2SD > 0.93 && I_p_wave_areas_5_SA_node_quotient_within_2SD < 0.97) && (I_p_wave_areas_5_SA_node_quotient_within_3SD > 0.98 && I_p_wave_areas_5_SA_node_quotient_within_3SD < 1)) 
    isNormal(11,1) = 1;
end
if I_p_wave_areas_5_SA_node_quotient_within_1SD == 0
    isSpike(11,1) = 1;
end
I_p_wave_areas_5_crista_terminalis_sorted = sort(I_p_wave_areas_5_crista_terminalis);
I_p_wave_areas_5_crista_terminalis_dev_avg_min_1SD = area_5(length(I_p_wave_areas_5_crista_terminalis_sorted),1);
for i = 1:length(I_p_wave_areas_5_crista_terminalis_sorted)
    I_p_wave_areas_5_crista_terminalis_dev_avg_min_1SD(i) = I_p_wave_areas_5_crista_terminalis_sorted(i) - (I_p_wave_areas_5_crista_terminalis_average-I_p_wave_areas_5_crista_terminalis_SD);
end
I_p_wave_areas_5_crista_terminalis_dev_avg_min_1SD = abs(I_p_wave_areas_5_crista_terminalis_dev_avg_min_1SD);
for i = 1:length(I_p_wave_areas_5_crista_terminalis_sorted)
    if I_p_wave_areas_5_crista_terminalis_dev_avg_min_1SD(i) == min(I_p_wave_areas_5_crista_terminalis_dev_avg_min_1SD)
        I_p_wave_areas_5_crista_terminalis_dev_avg_min_1SD_index = i;
    end
end
I_p_wave_areas_5_crista_terminalis_dev_avg_pls_1SD = area_5(length(I_p_wave_areas_5_crista_terminalis_sorted),1);
for i = 1:length(I_p_wave_areas_5_crista_terminalis_sorted)
    I_p_wave_areas_5_crista_terminalis_dev_avg_pls_1SD(i) = I_p_wave_areas_5_crista_terminalis_sorted(i) - (I_p_wave_areas_5_crista_terminalis_average+I_p_wave_areas_5_crista_terminalis_SD);
end
I_p_wave_areas_5_crista_terminalis_dev_avg_pls_1SD = abs(I_p_wave_areas_5_crista_terminalis_dev_avg_pls_1SD);
for i = 1:length(I_p_wave_areas_5_crista_terminalis_sorted)
    if I_p_wave_areas_5_crista_terminalis_dev_avg_pls_1SD(i) == min(I_p_wave_areas_5_crista_terminalis_dev_avg_pls_1SD)
        I_p_wave_areas_5_crista_terminalis_dev_avg_pls_1SD_index = i;
    end
end
I_p_wave_areas_5_crista_terminalis_dev_avg_min_2SD = area_5(length(I_p_wave_areas_5_crista_terminalis_sorted),1);
for i = 1:length(I_p_wave_areas_5_crista_terminalis_sorted)
    I_p_wave_areas_5_crista_terminalis_dev_avg_min_2SD(i) = I_p_wave_areas_5_crista_terminalis_sorted(i) - (I_p_wave_areas_5_crista_terminalis_average-2*I_p_wave_areas_5_crista_terminalis_SD);
end
I_p_wave_areas_5_crista_terminalis_dev_avg_min_2SD = abs(I_p_wave_areas_5_crista_terminalis_dev_avg_min_2SD);
for i = 1:length(I_p_wave_areas_5_crista_terminalis_sorted)
    if I_p_wave_areas_5_crista_terminalis_dev_avg_min_2SD(i) == min(I_p_wave_areas_5_crista_terminalis_dev_avg_min_2SD)
        I_p_wave_areas_5_dev_avg_min_2SD_index = i;
    end
end
I_p_wave_areas_5_crista_terminalis_dev_avg_pls_2SD = area_5(length(I_p_wave_areas_5_crista_terminalis_sorted),1);
for i = 1:length(I_p_wave_areas_5_crista_terminalis_sorted)
    I_p_wave_areas_5_crista_terminalis_dev_avg_pls_2SD(i) = I_p_wave_areas_5_crista_terminalis_sorted(i) - (I_p_wave_areas_5_crista_terminalis_average+2*I_p_wave_areas_5_crista_terminalis_SD);
end
I_p_wave_areas_5_crista_terminalis_dev_avg_pls_2SD = abs(I_p_wave_areas_5_crista_terminalis_dev_avg_pls_2SD);
for i = 1:length(I_p_wave_areas_5_crista_terminalis_sorted)
    if I_p_wave_areas_5_crista_terminalis_dev_avg_pls_2SD(i) == min(I_p_wave_areas_5_crista_terminalis_dev_avg_pls_2SD)
        I_p_wave_areas_5_crista_terminalis_dev_avg_pls_2SD_index = i;
    end
end
I_p_wave_areas_5_crista_terminalis_dev_avg_min_3SD = area_5(length(I_p_wave_areas_5_crista_terminalis_sorted),1);
for i = 1:length(I_p_wave_areas_5_crista_terminalis_sorted)
    I_p_wave_areas_5_crista_terminalis_dev_avg_min_3SD(i) = I_p_wave_areas_5_crista_terminalis_sorted(i) - (I_p_wave_areas_5_crista_terminalis_average-3*I_p_wave_areas_5_crista_terminalis_SD);
end
I_p_wave_areas_5_crista_terminalis_dev_avg_min_3SD = abs(I_p_wave_areas_5_crista_terminalis_dev_avg_min_3SD);
for i = 1:length(I_p_wave_areas_5_crista_terminalis_sorted)
    if I_p_wave_areas_5_crista_terminalis_dev_avg_min_3SD(i) == min(I_p_wave_areas_5_crista_terminalis_dev_avg_min_3SD)
        I_p_wave_areas_5_crista_terminalis_dev_avg_min_3SD_index = i;
    end
end
I_p_wave_areas_5_crista_terminalis_dev_avg_pls_3SD = area_5(length(I_p_wave_areas_5_crista_terminalis_sorted),1);
for i = 1:length(I_p_wave_areas_5_crista_terminalis_sorted)
    I_p_wave_areas_5_crista_terminalis_dev_avg_pls_3SD(i) = I_p_wave_areas_5_crista_terminalis_sorted(i) - (I_p_wave_areas_5_crista_terminalis_average+3*I_p_wave_areas_5_crista_terminalis_SD);
end
I_p_wave_areas_5_crista_terminalis_dev_avg_pls_3SD = abs(I_p_wave_areas_5_crista_terminalis_dev_avg_pls_3SD);
for i = 1:length(I_p_wave_areas_5_crista_terminalis_sorted)
    if I_p_wave_areas_5_crista_terminalis_dev_avg_pls_3SD(i) == min(I_p_wave_areas_5_crista_terminalis_dev_avg_pls_3SD)
        I_p_wave_areas_5_crista_terminalis_dev_avg_pls_3SD_index = i;
    end
end
I_p_wave_areas_5_crista_terminalis_tot_number_of_indices = length(I_p_wave_areas_5_crista_terminalis_sorted);
%This is the approximation of the tot "area under curve." 
I_p_wave_areas_5_crista_terminalis_quotient_within_1SD = (I_p_wave_areas_5_crista_terminalis_dev_avg_pls_1SD_index - I_p_wave_areas_5_crista_terminalisdev_avg_min_1SD_index)/I_p_wave_areas_5_crista_terminalis_tot_number_of_indices;
I_p_wave_areas_5_crista_terminalis_quotient_within_2SD = (I_p_wave_areas_5_crista_terminalis_dev_avg_pls_2SD_index - I_p_wave_areas_5_dev_avg_min_2SD_index)/I_p_wave_areas_5_crista_terminalis_tot_number_of_indices;
I_p_wave_areas_5_crista_terminalis_quotient_within_3SD = (I_p_wave_areas_5_crista_terminalis_dev_avg_pls_3SD_index - I_p_wave_areas_5_crista_terminalis_dev_avg_min_3SD_index)/I_p_wave_areas_5_crista_terminalis_tot_number_of_indices;
if (I_p_wave_areas_5_crista_terminalis_quotient_within_1SD > 0.66 && I_p_wave_areas_5_crista_terminalis_quotient_within_1SD < 0.70) && (I_p_wave_areas_5_crista_terminalis_quotient_within_2SD > 0.93 && I_p_wave_areas_5_crista_terminalis_quotient_within_2SD < 0.97) && (I_p_wave_areas_5_crista_terminalis_quotient_within_3SD > 0.98 && I_p_wave_areas_5_crista_terminalis_quotient_within_3SD < 1)
    isNormal(11,2) = 1;
end
if I_p_wave_areas_5_crista_terminalis_quotient_within_1SD == 0
    isSpike(11,2) = 1;
end
I_p_wave_areas_5_tricuspid_annulus_sorted = sort(I_p_wave_areas_5_tricuspid_annulus);
I_p_wave_areas_5_tricuspid_annulus_dev_avg_min_1SD = area_5(length(I_p_wave_areas_5_tricuspid_annulus_sorted),1);
for i = 1:length(I_p_wave_areas_5_tricuspid_annulus_sorted)
    I_p_wave_areas_5_tricuspid_annulus_dev_avg_min_1SD(i) = I_p_wave_areas_5_tricuspid_annulus_sorted(i) - (I_p_wave_areas_5_tricuspid_annulus_average-I_p_wave_areas_5_tricuspid_annulus_SD);
end
I_p_wave_areas_5_tricuspid_annulus_dev_avg_min_1SD = abs(I_p_wave_areas_5_tricuspid_annulus_dev_avg_min_1SD);
for i = 1:length(I_p_wave_areas_5_tricuspid_annulus_sorted)
    if I_p_wave_areas_5_tricuspid_annulus_dev_avg_min_1SD(i) == min(I_p_wave_areas_5_tricuspid_annulus_dev_avg_min_1SD)
        I_p_wave_areas_5_tricuspid_annulus_dev_avg_min_1SD_index = i;
    end
end
I_p_wave_areas_5_tricuspid_annulus_dev_avg_pls_1SD = area_5(length(I_p_wave_areas_5_tricuspid_annulus_sorted),1);
for i = 1:length(I_p_wave_areas_5_tricuspid_annulus_sorted)
    I_p_wave_areas_5_tricuspid_annulus_dev_avg_pls_1SD(i) = I_p_wave_areas_5_tricuspid_annulus_sorted(i) - (I_p_wave_areas_5_tricuspid_annulus_average+I_p_wave_areas_5_tricuspid_annulus_SD);
end
I_p_wave_areas_5_tricuspid_annulus_dev_avg_pls_1SD = abs(I_p_wave_areas_5_tricuspid_annulus_dev_avg_pls_1SD);
for i = 1:length(I_p_wave_areas_5_tricuspid_annulus_sorted)
    if I_p_wave_areas_5_tricuspid_annulus_dev_avg_pls_1SD(i) == min(I_p_wave_areas_5_tricuspid_annulus_dev_avg_pls_1SD)
        I_p_wave_areas_5_tricuspid_annulus_dev_avg_pls_1SD_index = i;
    end
end
I_p_wave_areas_5_tricuspid_annulus_dev_avg_min_2SD = area_5(length(I_p_wave_areas_5_tricuspid_annulus_sorted),1);
for i = 1:length(I_p_wave_areas_5_tricuspid_annulus_sorted)
    I_p_wave_areas_5_tricuspid_annulus_dev_avg_min_2SD(i) = I_p_wave_areas_5_tricuspid_annulus_sorted(i) - (I_p_wave_areas_5_tricuspid_annulus_average-2*I_p_wave_areas_5_tricuspid_annulus_SD);
end
I_p_wave_areas_5_tricuspid_annulus_dev_avg_min_2SD = abs(I_p_wave_areas_5_tricuspid_annulus_dev_avg_min_2SD);
for i = 1:length(I_p_wave_areas_5_tricuspid_annulus_sorted)
    if I_p_wave_areas_5_tricuspid_annulus_dev_avg_min_2SD(i) == min(I_p_wave_areas_5_tricuspid_annulus_dev_avg_min_2SD)
        I_p_wave_areas_5_dev_avg_min_2SD_index = i;
    end
end
I_p_wave_areas_5_tricuspid_annulus_dev_avg_pls_2SD = area_5(length(I_p_wave_areas_5_tricuspid_annulus_sorted),1);
for i = 1:length(I_p_wave_areas_5_tricuspid_annulus_sorted)
    I_p_wave_areas_5_tricuspid_annulus_dev_avg_pls_2SD(i) = I_p_wave_areas_5_tricuspid_annulus_sorted(i) - (I_p_wave_areas_5_tricuspid_annulus_average+2*I_p_wave_areas_5_tricuspid_annulus_SD);
end
I_p_wave_areas_5_tricuspid_annulus_dev_avg_pls_2SD = abs(I_p_wave_areas_5_tricuspid_annulus_dev_avg_pls_2SD);
for i = 1:length(I_p_wave_areas_5_tricuspid_annulus_sorted)
    if I_p_wave_areas_5_tricuspid_annulus_dev_avg_pls_2SD(i) == min(I_p_wave_areas_5_tricuspid_annulus_dev_avg_pls_2SD)
        I_p_wave_areas_5_tricuspid_annulus_dev_avg_pls_2SD_index = i;
    end
end
I_p_wave_areas_5_tricuspid_annulus_dev_avg_min_3SD = area_5(length(I_p_wave_areas_5_tricuspid_annulus_sorted),1);
for i = 1:length(I_p_wave_areas_5_tricuspid_annulus_sorted)
    I_p_wave_areas_5_tricuspid_annulus_dev_avg_min_3SD(i) = I_p_wave_areas_5_tricuspid_annulus_sorted(i) - (I_p_wave_areas_5_tricuspid_annulus_average-3*I_p_wave_areas_5_tricuspid_annulus_SD);
end
I_p_wave_areas_5_tricuspid_annulus_dev_avg_min_3SD = abs(I_p_wave_areas_5_tricuspid_annulus_dev_avg_min_3SD);
for i = 1:length(I_p_wave_areas_5_tricuspid_annulus_sorted)
    if I_p_wave_areas_5_tricuspid_annulus_dev_avg_min_3SD(i) == min(I_p_wave_areas_5_tricuspid_annulus_dev_avg_min_3SD)
        I_p_wave_areas_5_tricuspid_annulus_dev_avg_min_3SD_index = i;
    end
end
I_p_wave_areas_5_tricuspid_annulus_dev_avg_pls_3SD = area_5(length(I_p_wave_areas_5_tricuspid_annulus_sorted),1);
for i = 1:length(I_p_wave_areas_5_tricuspid_annulus_sorted)
    I_p_wave_areas_5_tricuspid_annulus_dev_avg_pls_3SD(i) = I_p_wave_areas_5_tricuspid_annulus_sorted(i) - (I_p_wave_areas_5_tricuspid_annulus_average+3*I_p_wave_areas_5_tricuspid_annulus_SD);
end
I_p_wave_areas_5_tricuspid_annulus_dev_avg_pls_3SD = abs(I_p_wave_areas_5_tricuspid_annulus_dev_avg_pls_3SD);
for i = 1:length(I_p_wave_areas_5_tricuspid_annulus_sorted)
    if I_p_wave_areas_5_tricuspid_annulus_dev_avg_pls_3SD(i) == min(I_p_wave_areas_5_tricuspid_annulus_dev_avg_pls_3SD)
        I_p_wave_areas_5_tricuspid_annulus_dev_avg_pls_3SD_index = i;
    end
end
I_p_wave_areas_5_tricuspid_annulus_tot_number_of_indices = length(I_p_wave_areas_5_tricuspid_annulus_sorted);
%This is the approximation of the tot "area under curve." 
I_p_wave_areas_5_tricuspid_annulus_quotient_within_1SD = (I_p_wave_areas_5_tricuspid_annulus_dev_avg_pls_1SD_index - I_p_wave_areas_5_tricuspid_annulusdev_avg_min_1SD_index)/I_p_wave_areas_5_tricuspid_annulus_tot_number_of_indices;
I_p_wave_areas_5_tricuspid_annulus_quotient_within_2SD = (I_p_wave_areas_5_tricuspid_annulus_dev_avg_pls_2SD_index - I_p_wave_areas_5_dev_avg_min_2SD_index)/I_p_wave_areas_5_tricuspid_annulus_tot_number_of_indices;
I_p_wave_areas_5_tricuspid_annulus_quotient_within_3SD = (I_p_wave_areas_5_tricuspid_annulus_dev_avg_pls_3SD_index - I_p_wave_areas_5_tricuspid_annulus_dev_avg_min_3SD_index)/I_p_wave_areas_5_tricuspid_annulus_tot_number_of_indices;
if (I_p_wave_areas_5_tricuspid_annulus_quotient_within_1SD > 0.66 && I_p_wave_areas_5_tricuspid_annulus_quotient_within_1SD < 0.70) && (I_p_wave_areas_5_tricuspid_annulus_quotient_within_2SD > 0.93 && I_p_wave_areas_5_tricuspid_annulus_quotient_within_2SD < 0.97) && (I_p_wave_areas_5_tricuspid_annulus_quotient_within_3SD > 0.98 && I_p_wave_areas_5_tricuspid_annulus_quotient_within_3SD < 1)
    isNormal(11,3) = 1;
end
if I_p_wave_areas_5_tricuspid_annulus_quotient_within_1SD == 0
    isSpike(11,3) = 1;
end
I_p_wave_areas_5_coronary_sinus_sorted = sort(I_p_wave_areas_5_coronary_sinus);
I_p_wave_areas_5_coronary_sinus_dev_avg_min_1SD = area_5(length(I_p_wave_areas_5_coronary_sinus_sorted),1);
for i = 1:length(I_p_wave_areas_5_coronary_sinus_sorted)
    I_p_wave_areas_5_coronary_sinus_dev_avg_min_1SD(i) = I_p_wave_areas_5_coronary_sinus_sorted(i) - (I_p_wave_areas_5_coronary_sinus_average-I_p_wave_areas_5_coronary_sinus_SD);
end
I_p_wave_areas_5_coronary_sinus_dev_avg_min_1SD = abs(I_p_wave_areas_5_coronary_sinus_dev_avg_min_1SD);
for i = 1:length(I_p_wave_areas_5_coronary_sinus_sorted)
    if I_p_wave_areas_5_coronary_sinus_dev_avg_min_1SD(i) == min(I_p_wave_areas_5_coronary_sinus_dev_avg_min_1SD)
        I_p_wave_areas_5_coronary_sinus_dev_avg_min_1SD_index = i;
    end
end
I_p_wave_areas_5_coronary_sinus_dev_avg_pls_1SD = area_5(length(I_p_wave_areas_5_coronary_sinus_sorted),1);
for i = 1:length(I_p_wave_areas_5_coronary_sinus_sorted)
    I_p_wave_areas_5_coronary_sinus_dev_avg_pls_1SD(i) = I_p_wave_areas_5_coronary_sinus_sorted(i) - (I_p_wave_areas_5_coronary_sinus_average+I_p_wave_areas_5_coronary_sinus_SD);
end
I_p_wave_areas_5_coronary_sinus_dev_avg_pls_1SD = abs(I_p_wave_areas_5_coronary_sinus_dev_avg_pls_1SD);
for i = 1:length(I_p_wave_areas_5_coronary_sinus_sorted)
    if I_p_wave_areas_5_coronary_sinus_dev_avg_pls_1SD(i) == min(I_p_wave_areas_5_coronary_sinus_dev_avg_pls_1SD)
        I_p_wave_areas_5_coronary_sinus_dev_avg_pls_1SD_index = i;
    end
end
I_p_wave_areas_5_coronary_sinus_dev_avg_min_2SD = area_5(length(I_p_wave_areas_5_coronary_sinus_sorted),1);
for i = 1:length(I_p_wave_areas_5_coronary_sinus_sorted)
    I_p_wave_areas_5_coronary_sinus_dev_avg_min_2SD(i) = I_p_wave_areas_5_coronary_sinus_sorted(i) - (I_p_wave_areas_5_coronary_sinus_average-2*I_p_wave_areas_5_coronary_sinus_SD);
end
I_p_wave_areas_5_coronary_sinus_dev_avg_min_2SD = abs(I_p_wave_areas_5_coronary_sinus_dev_avg_min_2SD);
for i = 1:length(I_p_wave_areas_5_coronary_sinus_sorted)
    if I_p_wave_areas_5_coronary_sinus_dev_avg_min_2SD(i) == min(I_p_wave_areas_5_coronary_sinus_dev_avg_min_2SD)
        I_p_wave_areas_5_dev_avg_min_2SD_index = i;
    end
end
I_p_wave_areas_5_coronary_sinus_dev_avg_pls_2SD = area_5(length(I_p_wave_areas_5_coronary_sinus_sorted),1);
for i = 1:length(I_p_wave_areas_5_coronary_sinus_sorted)
    I_p_wave_areas_5_coronary_sinus_dev_avg_pls_2SD(i) = I_p_wave_areas_5_coronary_sinus_sorted(i) - (I_p_wave_areas_5_coronary_sinus_average+2*I_p_wave_areas_5_coronary_sinus_SD);
end
I_p_wave_areas_5_coronary_sinus_dev_avg_pls_2SD = abs(I_p_wave_areas_5_coronary_sinus_dev_avg_pls_2SD);
for i = 1:length(I_p_wave_areas_5_coronary_sinus_sorted)
    if I_p_wave_areas_5_coronary_sinus_dev_avg_pls_2SD(i) == min(I_p_wave_areas_5_coronary_sinus_dev_avg_pls_2SD)
        I_p_wave_areas_5_coronary_sinus_dev_avg_pls_2SD_index = i;
    end
end
I_p_wave_areas_5_coronary_sinus_dev_avg_min_3SD = area_5(length(I_p_wave_areas_5_coronary_sinus_sorted),1);
for i = 1:length(I_p_wave_areas_5_coronary_sinus_sorted)
    I_p_wave_areas_5_coronary_sinus_dev_avg_min_3SD(i) = I_p_wave_areas_5_coronary_sinus_sorted(i) - (I_p_wave_areas_5_coronary_sinus_average-3*I_p_wave_areas_5_coronary_sinus_SD);
end
I_p_wave_areas_5_coronary_sinus_dev_avg_min_3SD = abs(I_p_wave_areas_5_coronary_sinus_dev_avg_min_3SD);
for i = 1:length(I_p_wave_areas_5_coronary_sinus_sorted)
    if I_p_wave_areas_5_coronary_sinus_dev_avg_min_3SD(i) == min(I_p_wave_areas_5_coronary_sinus_dev_avg_min_3SD)
        I_p_wave_areas_5_coronary_sinus_dev_avg_min_3SD_index = i;
    end
end
I_p_wave_areas_5_coronary_sinus_dev_avg_pls_3SD = area_5(length(I_p_wave_areas_5_coronary_sinus_sorted),1);
for i = 1:length(I_p_wave_areas_5_coronary_sinus_sorted)
    I_p_wave_areas_5_coronary_sinus_dev_avg_pls_3SD(i) = I_p_wave_areas_5_coronary_sinus_sorted(i) - (I_p_wave_areas_5_coronary_sinus_average+3*I_p_wave_areas_5_coronary_sinus_SD);
end
I_p_wave_areas_5_coronary_sinus_dev_avg_pls_3SD = abs(I_p_wave_areas_5_coronary_sinus_dev_avg_pls_3SD);
for i = 1:length(I_p_wave_areas_5_coronary_sinus_sorted)
    if I_p_wave_areas_5_coronary_sinus_dev_avg_pls_3SD(i) == min(I_p_wave_areas_5_coronary_sinus_dev_avg_pls_3SD)
        I_p_wave_areas_5_coronary_sinus_dev_avg_pls_3SD_index = i;
    end
end
I_p_wave_areas_5_coronary_sinus_tot_number_of_indices = length(I_p_wave_areas_5_coronary_sinus_sorted);
%This is the approximation of the tot "area under curve." 
I_p_wave_areas_5_coronary_sinus_quotient_within_1SD = (I_p_wave_areas_5_coronary_sinus_dev_avg_pls_1SD_index - I_p_wave_areas_5_coronary_sinusdev_avg_min_1SD_index)/I_p_wave_areas_5_coronary_sinus_tot_number_of_indices;
I_p_wave_areas_5_coronary_sinus_quotient_within_2SD = (I_p_wave_areas_5_coronary_sinus_dev_avg_pls_2SD_index - I_p_wave_areas_5_dev_avg_min_2SD_index)/I_p_wave_areas_5_coronary_sinus_tot_number_of_indices;
I_p_wave_areas_5_coronary_sinus_quotient_within_3SD = (I_p_wave_areas_5_coronary_sinus_dev_avg_pls_3SD_index - I_p_wave_areas_5_coronary_sinus_dev_avg_min_3SD_index)/I_p_wave_areas_5_coronary_sinus_tot_number_of_indices;
if (I_p_wave_areas_5_coronary_sinus_quotient_within_1SD > 0.66 && I_p_wave_areas_5_coronary_sinus_quotient_within_1SD < 0.70) && (I_p_wave_areas_5_coronary_sinus_quotient_within_2SD > 0.93 && I_p_wave_areas_5_coronary_sinus_quotient_within_2SD < 0.97) && (I_p_wave_areas_5_coronary_sinus_quotient_within_3SD > 0.98 && I_p_wave_areas_5_coronary_sinus_quotient_within_3SD < 1)
    isNormal(11,4) = 1;
end
if I_p_wave_areas_5_coronary_sinus_quotient_within_1SD == 0
    isSpike(11,4) = 1;
end
I_p_wave_areas_5_ostium_sorted = sort(I_p_wave_areas_5_ostium);
I_p_wave_areas_5_ostium_dev_avg_min_1SD = area_5(length(I_p_wave_areas_5_ostium_sorted),1);
for i = 1:length(I_p_wave_areas_5_ostium_sorted)
    I_p_wave_areas_5_ostium_dev_avg_min_1SD(i) = I_p_wave_areas_5_ostium_sorted(i) - (I_p_wave_areas_5_ostium_average-I_p_wave_areas_5_ostium_SD);
end
I_p_wave_areas_5_ostium_dev_avg_min_1SD = abs(I_p_wave_areas_5_ostium_dev_avg_min_1SD);
for i = 1:length(I_p_wave_areas_5_ostium_sorted)
    if I_p_wave_areas_5_ostium_dev_avg_min_1SD(i) == min(I_p_wave_areas_5_ostium_dev_avg_min_1SD)
        I_p_wave_areas_5_ostium_dev_avg_min_1SD_index = i;
    end
end
I_p_wave_areas_5_ostium_dev_avg_pls_1SD = area_5(length(I_p_wave_areas_5_ostium_sorted),1);
for i = 1:length(I_p_wave_areas_5_ostium_sorted)
    I_p_wave_areas_5_ostium_dev_avg_pls_1SD(i) = I_p_wave_areas_5_ostium_sorted(i) - (I_p_wave_areas_5_ostium_average+I_p_wave_areas_5_ostium_SD);
end
I_p_wave_areas_5_ostium_dev_avg_pls_1SD = abs(I_p_wave_areas_5_ostium_dev_avg_pls_1SD);
for i = 1:length(I_p_wave_areas_5_ostium_sorted)
    if I_p_wave_areas_5_ostium_dev_avg_pls_1SD(i) == min(I_p_wave_areas_5_ostium_dev_avg_pls_1SD)
        I_p_wave_areas_5_ostium_dev_avg_pls_1SD_index = i;
    end
end
I_p_wave_areas_5_ostium_dev_avg_min_2SD = area_5(length(I_p_wave_areas_5_ostium_sorted),1);
for i = 1:length(I_p_wave_areas_5_ostium_sorted)
    I_p_wave_areas_5_ostium_dev_avg_min_2SD(i) = I_p_wave_areas_5_ostium_sorted(i) - (I_p_wave_areas_5_ostium_average-2*I_p_wave_areas_5_ostium_SD);
end
I_p_wave_areas_5_ostium_dev_avg_min_2SD = abs(I_p_wave_areas_5_ostium_dev_avg_min_2SD);
for i = 1:length(I_p_wave_areas_5_ostium_sorted)
    if I_p_wave_areas_5_ostium_dev_avg_min_2SD(i) == min(I_p_wave_areas_5_ostium_dev_avg_min_2SD)
        I_p_wave_areas_5_dev_avg_min_2SD_index = i;
    end
end
I_p_wave_areas_5_ostium_dev_avg_pls_2SD = area_5(length(I_p_wave_areas_5_ostium_sorted),1);
for i = 1:length(I_p_wave_areas_5_ostium_sorted)
    I_p_wave_areas_5_ostium_dev_avg_pls_2SD(i) = I_p_wave_areas_5_ostium_sorted(i) - (I_p_wave_areas_5_ostium_average+2*I_p_wave_areas_5_ostium_SD);
end
I_p_wave_areas_5_ostium_dev_avg_pls_2SD = abs(I_p_wave_areas_5_ostium_dev_avg_pls_2SD);
for i = 1:length(I_p_wave_areas_5_ostium_sorted)
    if I_p_wave_areas_5_ostium_dev_avg_pls_2SD(i) == min(I_p_wave_areas_5_ostium_dev_avg_pls_2SD)
        I_p_wave_areas_5_ostium_dev_avg_pls_2SD_index = i;
    end
end
I_p_wave_areas_5_ostium_dev_avg_min_3SD = area_5(length(I_p_wave_areas_5_ostium_sorted),1);
for i = 1:length(I_p_wave_areas_5_ostium_sorted)
    I_p_wave_areas_5_ostium_dev_avg_min_3SD(i) = I_p_wave_areas_5_ostium_sorted(i) - (I_p_wave_areas_5_ostium_average-3*I_p_wave_areas_5_ostium_SD);
end
I_p_wave_areas_5_ostium_dev_avg_min_3SD = abs(I_p_wave_areas_5_ostium_dev_avg_min_3SD);
for i = 1:length(I_p_wave_areas_5_ostium_sorted)
    if I_p_wave_areas_5_ostium_dev_avg_min_3SD(i) == min(I_p_wave_areas_5_ostium_dev_avg_min_3SD)
        I_p_wave_areas_5_ostium_dev_avg_min_3SD_index = i;
    end
end
I_p_wave_areas_5_ostium_dev_avg_pls_3SD = area_5(length(I_p_wave_areas_5_ostium_sorted),1);
for i = 1:length(I_p_wave_areas_5_ostium_sorted)
    I_p_wave_areas_5_ostium_dev_avg_pls_3SD(i) = I_p_wave_areas_5_ostium_sorted(i) - (I_p_wave_areas_5_ostium_average+3*I_p_wave_areas_5_ostium_SD);
end
I_p_wave_areas_5_ostium_dev_avg_pls_3SD = abs(I_p_wave_areas_5_ostium_dev_avg_pls_3SD);
for i = 1:length(I_p_wave_areas_5_ostium_sorted)
    if I_p_wave_areas_5_ostium_dev_avg_pls_3SD(i) == min(I_p_wave_areas_5_ostium_dev_avg_pls_3SD)
        I_p_wave_areas_5_ostium_dev_avg_pls_3SD_index = i;
    end
end
I_p_wave_areas_5_ostium_tot_number_of_indices = length(I_p_wave_areas_5_ostium_sorted);
%This is the approximation of the tot "area under curve." 
I_p_wave_areas_5_ostium_quotient_within_1SD = (I_p_wave_areas_5_ostium_dev_avg_pls_1SD_index - I_p_wave_areas_5_ostiumdev_avg_min_1SD_index)/I_p_wave_areas_5_ostium_tot_number_of_indices;
I_p_wave_areas_5_ostium_quotient_within_2SD = (I_p_wave_areas_5_ostium_dev_avg_pls_2SD_index - I_p_wave_areas_5_dev_avg_min_2SD_index)/I_p_wave_areas_5_ostium_tot_number_of_indices;
I_p_wave_areas_5_ostium_quotient_within_3SD = (I_p_wave_areas_5_ostium_dev_avg_pls_3SD_index - I_p_wave_areas_5_ostium_dev_avg_min_3SD_index)/I_p_wave_areas_5_ostium_tot_number_of_indices;
if (I_p_wave_areas_5_ostium_quotient_within_1SD > 0.66 && I_p_wave_areas_5_ostium_quotient_within_1SD < 0.70) && (I_p_wave_areas_5_ostium_quotient_within_2SD > 0.93 && I_p_wave_areas_5_ostium_quotient_within_2SD < 0.97) && (I_p_wave_areas_5_ostium_quotient_within_3SD > 0.98 && I_p_wave_areas_5_ostium_quotient_within_3SD < 1)
    isNormal(11,5) = 1;
end
if I_p_wave_areas_5_ostium_quotient_within_1SD == 0
    isSpike(11,5) = 1;
end
I_p_wave_areas_5_perinodal_sorted = sort(I_p_wave_areas_5_perinodal);
I_p_wave_areas_5_perinodal_dev_avg_min_1SD = area_5(length(I_p_wave_areas_5_perinodal_sorted),1);
for i = 1:length(I_p_wave_areas_5_perinodal_sorted)
    I_p_wave_areas_5_perinodal_dev_avg_min_1SD(i) = I_p_wave_areas_5_perinodal_sorted(i) - (I_p_wave_areas_5_perinodal_average-I_p_wave_areas_5_perinodal_SD);
end
I_p_wave_areas_5_perinodal_dev_avg_min_1SD = abs(I_p_wave_areas_5_perinodal_dev_avg_min_1SD);
for i = 1:length(I_p_wave_areas_5_perinodal_sorted)
    if I_p_wave_areas_5_perinodal_dev_avg_min_1SD(i) == min(I_p_wave_areas_5_perinodal_dev_avg_min_1SD)
        I_p_wave_areas_5_perinodal_dev_avg_min_1SD_index = i;
    end
end
I_p_wave_areas_5_perinodal_dev_avg_pls_1SD = area_5(length(I_p_wave_areas_5_perinodal_sorted),1);
for i = 1:length(I_p_wave_areas_5_perinodal_sorted)
    I_p_wave_areas_5_perinodal_dev_avg_pls_1SD(i) = I_p_wave_areas_5_perinodal_sorted(i) - (I_p_wave_areas_5_perinodal_average+I_p_wave_areas_5_perinodal_SD);
end
I_p_wave_areas_5_perinodal_dev_avg_pls_1SD = abs(I_p_wave_areas_5_perinodal_dev_avg_pls_1SD);
for i = 1:length(I_p_wave_areas_5_perinodal_sorted)
    if I_p_wave_areas_5_perinodal_dev_avg_pls_1SD(i) == min(I_p_wave_areas_5_perinodal_dev_avg_pls_1SD)
        I_p_wave_areas_5_perinodal_dev_avg_pls_1SD_index = i;
    end
end
I_p_wave_areas_5_perinodal_dev_avg_min_2SD = area_5(length(I_p_wave_areas_5_perinodal_sorted),1);
for i = 1:length(I_p_wave_areas_5_perinodal_sorted)
    I_p_wave_areas_5_perinodal_dev_avg_min_2SD(i) = I_p_wave_areas_5_perinodal_sorted(i) - (I_p_wave_areas_5_perinodal_average-2*I_p_wave_areas_5_perinodal_SD);
end
I_p_wave_areas_5_perinodal_dev_avg_min_2SD = abs(I_p_wave_areas_5_perinodal_dev_avg_min_2SD);
for i = 1:length(I_p_wave_areas_5_perinodal_sorted)
    if I_p_wave_areas_5_perinodal_dev_avg_min_2SD(i) == min(I_p_wave_areas_5_perinodal_dev_avg_min_2SD)
        I_p_wave_areas_5_dev_avg_min_2SD_index = i;
    end
end
I_p_wave_areas_5_perinodal_dev_avg_pls_2SD = area_5(length(I_p_wave_areas_5_perinodal_sorted),1);
for i = 1:length(I_p_wave_areas_5_perinodal_sorted)
    I_p_wave_areas_5_perinodal_dev_avg_pls_2SD(i) = I_p_wave_areas_5_perinodal_sorted(i) - (I_p_wave_areas_5_perinodal_average+2*I_p_wave_areas_5_perinodal_SD);
end
I_p_wave_areas_5_perinodal_dev_avg_pls_2SD = abs(I_p_wave_areas_5_perinodal_dev_avg_pls_2SD);
for i = 1:length(I_p_wave_areas_5_perinodal_sorted)
    if I_p_wave_areas_5_perinodal_dev_avg_pls_2SD(i) == min(I_p_wave_areas_5_perinodal_dev_avg_pls_2SD)
        I_p_wave_areas_5_perinodal_dev_avg_pls_2SD_index = i;
    end
end
I_p_wave_areas_5_perinodal_dev_avg_min_3SD = area_5(length(I_p_wave_areas_5_perinodal_sorted),1);
for i = 1:length(I_p_wave_areas_5_perinodal_sorted)
    I_p_wave_areas_5_perinodal_dev_avg_min_3SD(i) = I_p_wave_areas_5_perinodal_sorted(i) - (I_p_wave_areas_5_perinodal_average-3*I_p_wave_areas_5_perinodal_SD);
end
I_p_wave_areas_5_perinodal_dev_avg_min_3SD = abs(I_p_wave_areas_5_perinodal_dev_avg_min_3SD);
for i = 1:length(I_p_wave_areas_5_perinodal_sorted)
    if I_p_wave_areas_5_perinodal_dev_avg_min_3SD(i) == min(I_p_wave_areas_5_perinodal_dev_avg_min_3SD)
        I_p_wave_areas_5_perinodal_dev_avg_min_3SD_index = i;
    end
end
I_p_wave_areas_5_perinodal_dev_avg_pls_3SD = area_5(length(I_p_wave_areas_5_perinodal_sorted),1);
for i = 1:length(I_p_wave_areas_5_perinodal_sorted)
    I_p_wave_areas_5_perinodal_dev_avg_pls_3SD(i) = I_p_wave_areas_5_perinodal_sorted(i) - (I_p_wave_areas_5_perinodal_average+3*I_p_wave_areas_5_perinodal_SD);
end
I_p_wave_areas_5_perinodal_dev_avg_pls_3SD = abs(I_p_wave_areas_5_perinodal_dev_avg_pls_3SD);
for i = 1:length(I_p_wave_areas_5_perinodal_sorted)
    if I_p_wave_areas_5_perinodal_dev_avg_pls_3SD(i) == min(I_p_wave_areas_5_perinodal_dev_avg_pls_3SD)
        I_p_wave_areas_5_perinodal_dev_avg_pls_3SD_index = i;
    end
end
I_p_wave_areas_5_perinodal_tot_number_of_indices = length(I_p_wave_areas_5_perinodal_sorted);
%This is the approximation of the tot "area under curve." 
I_p_wave_areas_5_perinodal_quotient_within_1SD = (I_p_wave_areas_5_perinodal_dev_avg_pls_1SD_index - I_p_wave_areas_5_perinodaldev_avg_min_1SD_index)/I_p_wave_areas_5_perinodal_tot_number_of_indices;
I_p_wave_areas_5_perinodal_quotient_within_2SD = (I_p_wave_areas_5_perinodal_dev_avg_pls_2SD_index - I_p_wave_areas_5_dev_avg_min_2SD_index)/I_p_wave_areas_5_perinodal_tot_number_of_indices;
I_p_wave_areas_5_perinodal_quotient_within_3SD = (I_p_wave_areas_5_perinodal_dev_avg_pls_3SD_index - I_p_wave_areas_5_perinodal_dev_avg_min_3SD_index)/I_p_wave_areas_5_perinodal_tot_number_of_indices;
if (I_p_wave_areas_5_perinodal_quotient_within_1SD > 0.66 && I_p_wave_areas_5_perinodal_quotient_within_1SD < 0.70) && (I_p_wave_areas_5_perinodal_quotient_within_2SD > 0.93 && I_p_wave_areas_5_perinodal_quotient_within_2SD < 0.97) && (I_p_wave_areas_5_perinodal_quotient_within_3SD > 0.98 && I_p_wave_areas_5_perinodal_quotient_within_3SD < 1)
    isNormal(11,6) = 1;
end
if I_p_wave_areas_5_perinodal_quotient_within_1SD == 0
    isSpike(11,6) = 1;
end
I_p_wave_areas_5_right_septum_sorted = sort(I_p_wave_areas_5_right_septum);
I_p_wave_areas_5_right_septum_dev_avg_min_1SD = area_5(length(I_p_wave_areas_5_right_septum_sorted),1);
for i = 1:length(I_p_wave_areas_5_right_septum_sorted)
    I_p_wave_areas_5_right_septum_dev_avg_min_1SD(i) = I_p_wave_areas_5_right_septum_sorted(i) - (I_p_wave_areas_5_right_septum_average-I_p_wave_areas_5_right_septum_SD);
end
I_p_wave_areas_5_right_septum_dev_avg_min_1SD = abs(I_p_wave_areas_5_right_septum_dev_avg_min_1SD);
for i = 1:length(I_p_wave_areas_5_right_septum_sorted)
    if I_p_wave_areas_5_right_septum_dev_avg_min_1SD(i) == min(I_p_wave_areas_5_right_septum_dev_avg_min_1SD)
        I_p_wave_areas_5_right_septum_dev_avg_min_1SD_index = i;
    end
end
I_p_wave_areas_5_right_septum_dev_avg_pls_1SD = area_5(length(I_p_wave_areas_5_right_septum_sorted),1);
for i = 1:length(I_p_wave_areas_5_right_septum_sorted)
    I_p_wave_areas_5_right_septum_dev_avg_pls_1SD(i) = I_p_wave_areas_5_right_septum_sorted(i) - (I_p_wave_areas_5_right_septum_average+I_p_wave_areas_5_right_septum_SD);
end
I_p_wave_areas_5_right_septum_dev_avg_pls_1SD = abs(I_p_wave_areas_5_right_septum_dev_avg_pls_1SD);
for i = 1:length(I_p_wave_areas_5_right_septum_sorted)
    if I_p_wave_areas_5_right_septum_dev_avg_pls_1SD(i) == min(I_p_wave_areas_5_right_septum_dev_avg_pls_1SD)
        I_p_wave_areas_5_right_septum_dev_avg_pls_1SD_index = i;
    end
end
I_p_wave_areas_5_right_septum_dev_avg_min_2SD = area_5(length(I_p_wave_areas_5_right_septum_sorted),1);
for i = 1:length(I_p_wave_areas_5_right_septum_sorted)
    I_p_wave_areas_5_right_septum_dev_avg_min_2SD(i) = I_p_wave_areas_5_right_septum_sorted(i) - (I_p_wave_areas_5_right_septum_average-2*I_p_wave_areas_5_right_septum_SD);
end
I_p_wave_areas_5_right_septum_dev_avg_min_2SD = abs(I_p_wave_areas_5_right_septum_dev_avg_min_2SD);
for i = 1:length(I_p_wave_areas_5_right_septum_sorted)
    if I_p_wave_areas_5_right_septum_dev_avg_min_2SD(i) == min(I_p_wave_areas_5_right_septum_dev_avg_min_2SD)
        I_p_wave_areas_5_dev_avg_min_2SD_index = i;
    end
end
I_p_wave_areas_5_right_septum_dev_avg_pls_2SD = area_5(length(I_p_wave_areas_5_right_septum_sorted),1);
for i = 1:length(I_p_wave_areas_5_right_septum_sorted)
    I_p_wave_areas_5_right_septum_dev_avg_pls_2SD(i) = I_p_wave_areas_5_right_septum_sorted(i) - (I_p_wave_areas_5_right_septum_average+2*I_p_wave_areas_5_right_septum_SD);
end
I_p_wave_areas_5_right_septum_dev_avg_pls_2SD = abs(I_p_wave_areas_5_right_septum_dev_avg_pls_2SD);
for i = 1:length(I_p_wave_areas_5_right_septum_sorted)
    if I_p_wave_areas_5_right_septum_dev_avg_pls_2SD(i) == min(I_p_wave_areas_5_right_septum_dev_avg_pls_2SD)
        I_p_wave_areas_5_right_septum_dev_avg_pls_2SD_index = i;
    end
end
I_p_wave_areas_5_right_septum_dev_avg_min_3SD = area_5(length(I_p_wave_areas_5_right_septum_sorted),1);
for i = 1:length(I_p_wave_areas_5_right_septum_sorted)
    I_p_wave_areas_5_right_septum_dev_avg_min_3SD(i) = I_p_wave_areas_5_right_septum_sorted(i) - (I_p_wave_areas_5_right_septum_average-3*I_p_wave_areas_5_right_septum_SD);
end
I_p_wave_areas_5_right_septum_dev_avg_min_3SD = abs(I_p_wave_areas_5_right_septum_dev_avg_min_3SD);
for i = 1:length(I_p_wave_areas_5_right_septum_sorted)
    if I_p_wave_areas_5_right_septum_dev_avg_min_3SD(i) == min(I_p_wave_areas_5_right_septum_dev_avg_min_3SD)
        I_p_wave_areas_5_right_septum_dev_avg_min_3SD_index = i;
    end
end
I_p_wave_areas_5_right_septum_dev_avg_pls_3SD = area_5(length(I_p_wave_areas_5_right_septum_sorted),1);
for i = 1:length(I_p_wave_areas_5_right_septum_sorted)
    I_p_wave_areas_5_right_septum_dev_avg_pls_3SD(i) = I_p_wave_areas_5_right_septum_sorted(i) - (I_p_wave_areas_5_right_septum_average+3*I_p_wave_areas_5_right_septum_SD);
end
I_p_wave_areas_5_right_septum_dev_avg_pls_3SD = abs(I_p_wave_areas_5_right_septum_dev_avg_pls_3SD);
for i = 1:length(I_p_wave_areas_5_right_septum_sorted)
    if I_p_wave_areas_5_right_septum_dev_avg_pls_3SD(i) == min(I_p_wave_areas_5_right_septum_dev_avg_pls_3SD)
        I_p_wave_areas_5_right_septum_dev_avg_pls_3SD_index = i;
    end
end
I_p_wave_areas_5_right_septum_tot_number_of_indices = length(I_p_wave_areas_5_right_septum_sorted);
%This is the approximation of the tot "area under curve." 
I_p_wave_areas_5_right_septum_quotient_within_1SD = (I_p_wave_areas_5_right_septum_dev_avg_pls_1SD_index - I_p_wave_areas_5_right_septumdev_avg_min_1SD_index)/I_p_wave_areas_5_right_septum_tot_number_of_indices;
I_p_wave_areas_5_right_septum_quotient_within_2SD = (I_p_wave_areas_5_right_septum_dev_avg_pls_2SD_index - I_p_wave_areas_5_dev_avg_min_2SD_index)/I_p_wave_areas_5_right_septum_tot_number_of_indices;
I_p_wave_areas_5_right_septum_quotient_within_3SD = (I_p_wave_areas_5_right_septum_dev_avg_pls_3SD_index - I_p_wave_areas_5_right_septum_dev_avg_min_3SD_index)/I_p_wave_areas_5_right_septum_tot_number_of_indices;
if (I_p_wave_areas_5_right_septum_quotient_within_1SD > 0.66 && I_p_wave_areas_5_right_septum_quotient_within_1SD < 0.70) && (I_p_wave_areas_5_right_septum_quotient_within_2SD > 0.93 && I_p_wave_areas_5_right_septum_quotient_within_2SD < 0.97) && (I_p_wave_areas_5_right_septum_quotient_within_3SD > 0.98 && I_p_wave_areas_5_right_septum_quotient_within_3SD < 1)
    isNormal(11,7) = 1;
end
if I_p_wave_areas_5_right_septum_quotient_within_1SD == 0
    isSpike(11,7) = 1;
end
I_p_wave_areas_5_right_atrial_appendage_sorted = sort(I_p_wave_areas_5_right_atrial_appendage);
I_p_wave_areas_5_right_atrial_appendage_dev_avg_min_1SD = area_5(length(I_p_wave_areas_5_right_atrial_appendage_sorted),1);
for i = 1:length(I_p_wave_areas_5_right_atrial_appendage_sorted)
    I_p_wave_areas_5_right_atrial_appendage_dev_avg_min_1SD(i) = I_p_wave_areas_5_right_atrial_appendage_sorted(i) - (I_p_wave_areas_5_right_atrial_appendage_average-I_p_wave_areas_5_right_atrial_appendage_SD);
end
I_p_wave_areas_5_right_atrial_appendage_dev_avg_min_1SD = abs(I_p_wave_areas_5_right_atrial_appendage_dev_avg_min_1SD);
for i = 1:length(I_p_wave_areas_5_right_atrial_appendage_sorted)
    if I_p_wave_areas_5_right_atrial_appendage_dev_avg_min_1SD(i) == min(I_p_wave_areas_5_right_atrial_appendage_dev_avg_min_1SD)
        I_p_wave_areas_5_right_atrial_appendage_dev_avg_min_1SD_index = i;
    end
end
I_p_wave_areas_5_right_atrial_appendage_dev_avg_pls_1SD = area_5(length(I_p_wave_areas_5_right_atrial_appendage_sorted),1);
for i = 1:length(I_p_wave_areas_5_right_atrial_appendage_sorted)
    I_p_wave_areas_5_right_atrial_appendage_dev_avg_pls_1SD(i) = I_p_wave_areas_5_right_atrial_appendage_sorted(i) - (I_p_wave_areas_5_right_atrial_appendage_average+I_p_wave_areas_5_right_atrial_appendage_SD);
end
I_p_wave_areas_5_right_atrial_appendage_dev_avg_pls_1SD = abs(I_p_wave_areas_5_right_atrial_appendage_dev_avg_pls_1SD);
for i = 1:length(I_p_wave_areas_5_right_atrial_appendage_sorted)
    if I_p_wave_areas_5_right_atrial_appendage_dev_avg_pls_1SD(i) == min(I_p_wave_areas_5_right_atrial_appendage_dev_avg_pls_1SD)
        I_p_wave_areas_5_right_atrial_appendage_dev_avg_pls_1SD_index = i;
    end
end
I_p_wave_areas_5_right_atrial_appendage_dev_avg_min_2SD = area_5(length(I_p_wave_areas_5_right_atrial_appendage_sorted),1);
for i = 1:length(I_p_wave_areas_5_right_atrial_appendage_sorted)
    I_p_wave_areas_5_right_atrial_appendage_dev_avg_min_2SD(i) = I_p_wave_areas_5_right_atrial_appendage_sorted(i) - (I_p_wave_areas_5_right_atrial_appendage_average-2*I_p_wave_areas_5_right_atrial_appendage_SD);
end
I_p_wave_areas_5_right_atrial_appendage_dev_avg_min_2SD = abs(I_p_wave_areas_5_right_atrial_appendage_dev_avg_min_2SD);
for i = 1:length(I_p_wave_areas_5_right_atrial_appendage_sorted)
    if I_p_wave_areas_5_right_atrial_appendage_dev_avg_min_2SD(i) == min(I_p_wave_areas_5_right_atrial_appendage_dev_avg_min_2SD)
        I_p_wave_areas_5_dev_avg_min_2SD_index = i;
    end
end
I_p_wave_areas_5_right_atrial_appendage_dev_avg_pls_2SD = area_5(length(I_p_wave_areas_5_right_atrial_appendage_sorted),1);
for i = 1:length(I_p_wave_areas_5_right_atrial_appendage_sorted)
    I_p_wave_areas_5_right_atrial_appendage_dev_avg_pls_2SD(i) = I_p_wave_areas_5_right_atrial_appendage_sorted(i) - (I_p_wave_areas_5_right_atrial_appendage_average+2*I_p_wave_areas_5_right_atrial_appendage_SD);
end
I_p_wave_areas_5_right_atrial_appendage_dev_avg_pls_2SD = abs(I_p_wave_areas_5_right_atrial_appendage_dev_avg_pls_2SD);
for i = 1:length(I_p_wave_areas_5_right_atrial_appendage_sorted)
    if I_p_wave_areas_5_right_atrial_appendage_dev_avg_pls_2SD(i) == min(I_p_wave_areas_5_right_atrial_appendage_dev_avg_pls_2SD)
        I_p_wave_areas_5_right_atrial_appendage_dev_avg_pls_2SD_index = i;
    end
end
I_p_wave_areas_5_right_atrial_appendage_dev_avg_min_3SD = area_5(length(I_p_wave_areas_5_right_atrial_appendage_sorted),1);
for i = 1:length(I_p_wave_areas_5_right_atrial_appendage_sorted)
    I_p_wave_areas_5_right_atrial_appendage_dev_avg_min_3SD(i) = I_p_wave_areas_5_right_atrial_appendage_sorted(i) - (I_p_wave_areas_5_right_atrial_appendage_average-3*I_p_wave_areas_5_right_atrial_appendage_SD);
end
I_p_wave_areas_5_right_atrial_appendage_dev_avg_min_3SD = abs(I_p_wave_areas_5_right_atrial_appendage_dev_avg_min_3SD);
for i = 1:length(I_p_wave_areas_5_right_atrial_appendage_sorted)
    if I_p_wave_areas_5_right_atrial_appendage_dev_avg_min_3SD(i) == min(I_p_wave_areas_5_right_atrial_appendage_dev_avg_min_3SD)
        I_p_wave_areas_5_right_atrial_appendage_dev_avg_min_3SD_index = i;
    end
end
I_p_wave_areas_5_right_atrial_appendage_dev_avg_pls_3SD = area_5(length(I_p_wave_areas_5_right_atrial_appendage_sorted),1);
for i = 1:length(I_p_wave_areas_5_right_atrial_appendage_sorted)
    I_p_wave_areas_5_right_atrial_appendage_dev_avg_pls_3SD(i) = I_p_wave_areas_5_right_atrial_appendage_sorted(i) - (I_p_wave_areas_5_right_atrial_appendage_average+3*I_p_wave_areas_5_right_atrial_appendage_SD);
end
I_p_wave_areas_5_right_atrial_appendage_dev_avg_pls_3SD = abs(I_p_wave_areas_5_right_atrial_appendage_dev_avg_pls_3SD);
for i = 1:length(I_p_wave_areas_5_right_atrial_appendage_sorted)
    if I_p_wave_areas_5_right_atrial_appendage_dev_avg_pls_3SD(i) == min(I_p_wave_areas_5_right_atrial_appendage_dev_avg_pls_3SD)
        I_p_wave_areas_5_right_atrial_appendage_dev_avg_pls_3SD_index = i;
    end
end
I_p_wave_areas_5_right_atrial_appendage_tot_number_of_indices = length(I_p_wave_areas_5_right_atrial_appendage_sorted);
%This is the approximation of the tot "area under curve." 
I_p_wave_areas_5_right_atrial_appendage_quotient_within_1SD = (I_p_wave_areas_5_right_atrial_appendage_dev_avg_pls_1SD_index - I_p_wave_areas_5_right_atrial_appendagedev_avg_min_1SD_index)/I_p_wave_areas_5_right_atrial_appendage_tot_number_of_indices;
I_p_wave_areas_5_right_atrial_appendage_quotient_within_2SD = (I_p_wave_areas_5_right_atrial_appendage_dev_avg_pls_2SD_index - I_p_wave_areas_5_dev_avg_min_2SD_index)/I_p_wave_areas_5_right_atrial_appendage_tot_number_of_indices;
I_p_wave_areas_5_right_atrial_appendage_quotient_within_3SD = (I_p_wave_areas_5_right_atrial_appendage_dev_avg_pls_3SD_index - I_p_wave_areas_5_right_atrial_appendage_dev_avg_min_3SD_index)/I_p_wave_areas_5_right_atrial_appendage_tot_number_of_indices;
if (I_p_wave_areas_5_right_atrial_appendage_quotient_within_1SD > 0.66 && I_p_wave_areas_5_right_atrial_appendage_quotient_within_1SD < 0.70) && (I_p_wave_areas_5_right_atrial_appendage_quotient_within_2SD > 0.93 && I_p_wave_areas_5_right_atrial_appendage_quotient_within_2SD < 0.97) && (I_p_wave_areas_5_right_atrial_appendage_quotient_within_3SD > 0.98 && I_p_wave_areas_5_right_atrial_appendage_quotient_within_3SD < 1)
    isNormal(11,8) = 1;
end
if I_p_wave_areas_5_right_atrial_appendage_quotient_within_1SD == 0
    isSpike(11,8) = 1;
end
I_p_wave_areas_5_pulmonary_veins_sorted = sort(I_p_wave_areas_5_pulmonary_veins);
I_p_wave_areas_5_pulmonary_veins_dev_avg_min_1SD = area_5(length(I_p_wave_areas_5_pulmonary_veins_sorted),1);
for i = 1:length(I_p_wave_areas_5_pulmonary_veins_sorted)
    I_p_wave_areas_5_pulmonary_veins_dev_avg_min_1SD(i) = I_p_wave_areas_5_pulmonary_veins_sorted(i) - (I_p_wave_areas_5_pulmonary_veins_average-I_p_wave_areas_5_pulmonary_veins_SD);
end
I_p_wave_areas_5_pulmonary_veins_dev_avg_min_1SD = abs(I_p_wave_areas_5_pulmonary_veins_dev_avg_min_1SD);
for i = 1:length(I_p_wave_areas_5_pulmonary_veins_sorted)
    if I_p_wave_areas_5_pulmonary_veins_dev_avg_min_1SD(i) == min(I_p_wave_areas_5_pulmonary_veins_dev_avg_min_1SD)
        I_p_wave_areas_5_pulmonary_veins_dev_avg_min_1SD_index = i;
    end
end
I_p_wave_areas_5_pulmonary_veins_dev_avg_pls_1SD = area_5(length(I_p_wave_areas_5_pulmonary_veins_sorted),1);
for i = 1:length(I_p_wave_areas_5_pulmonary_veins_sorted)
    I_p_wave_areas_5_pulmonary_veins_dev_avg_pls_1SD(i) = I_p_wave_areas_5_pulmonary_veins_sorted(i) - (I_p_wave_areas_5_pulmonary_veins_average+I_p_wave_areas_5_pulmonary_veins_SD);
end
I_p_wave_areas_5_pulmonary_veins_dev_avg_pls_1SD = abs(I_p_wave_areas_5_pulmonary_veins_dev_avg_pls_1SD);
for i = 1:length(I_p_wave_areas_5_pulmonary_veins_sorted)
    if I_p_wave_areas_5_pulmonary_veins_dev_avg_pls_1SD(i) == min(I_p_wave_areas_5_pulmonary_veins_dev_avg_pls_1SD)
        I_p_wave_areas_5_pulmonary_veins_dev_avg_pls_1SD_index = i;
    end
end
I_p_wave_areas_5_pulmonary_veins_dev_avg_min_2SD = area_5(length(I_p_wave_areas_5_pulmonary_veins_sorted),1);
for i = 1:length(I_p_wave_areas_5_pulmonary_veins_sorted)
    I_p_wave_areas_5_pulmonary_veins_dev_avg_min_2SD(i) = I_p_wave_areas_5_pulmonary_veins_sorted(i) - (I_p_wave_areas_5_pulmonary_veins_average-2*I_p_wave_areas_5_pulmonary_veins_SD);
end
I_p_wave_areas_5_pulmonary_veins_dev_avg_min_2SD = abs(I_p_wave_areas_5_pulmonary_veins_dev_avg_min_2SD);
for i = 1:length(I_p_wave_areas_5_pulmonary_veins_sorted)
    if I_p_wave_areas_5_pulmonary_veins_dev_avg_min_2SD(i) == min(I_p_wave_areas_5_pulmonary_veins_dev_avg_min_2SD)
        I_p_wave_areas_5_dev_avg_min_2SD_index = i;
    end
end
I_p_wave_areas_5_pulmonary_veins_dev_avg_pls_2SD = area_5(length(I_p_wave_areas_5_pulmonary_veins_sorted),1);
for i = 1:length(I_p_wave_areas_5_pulmonary_veins_sorted)
    I_p_wave_areas_5_pulmonary_veins_dev_avg_pls_2SD(i) = I_p_wave_areas_5_pulmonary_veins_sorted(i) - (I_p_wave_areas_5_pulmonary_veins_average+2*I_p_wave_areas_5_pulmonary_veins_SD);
end
I_p_wave_areas_5_pulmonary_veins_dev_avg_pls_2SD = abs(I_p_wave_areas_5_pulmonary_veins_dev_avg_pls_2SD);
for i = 1:length(I_p_wave_areas_5_pulmonary_veins_sorted)
    if I_p_wave_areas_5_pulmonary_veins_dev_avg_pls_2SD(i) == min(I_p_wave_areas_5_pulmonary_veins_dev_avg_pls_2SD)
        I_p_wave_areas_5_pulmonary_veins_dev_avg_pls_2SD_index = i;
    end
end
I_p_wave_areas_5_pulmonary_veins_dev_avg_min_3SD = area_5(length(I_p_wave_areas_5_pulmonary_veins_sorted),1);
for i = 1:length(I_p_wave_areas_5_pulmonary_veins_sorted)
    I_p_wave_areas_5_pulmonary_veins_dev_avg_min_3SD(i) = I_p_wave_areas_5_pulmonary_veins_sorted(i) - (I_p_wave_areas_5_pulmonary_veins_average-3*I_p_wave_areas_5_pulmonary_veins_SD);
end
I_p_wave_areas_5_pulmonary_veins_dev_avg_min_3SD = abs(I_p_wave_areas_5_pulmonary_veins_dev_avg_min_3SD);
for i = 1:length(I_p_wave_areas_5_pulmonary_veins_sorted)
    if I_p_wave_areas_5_pulmonary_veins_dev_avg_min_3SD(i) == min(I_p_wave_areas_5_pulmonary_veins_dev_avg_min_3SD)
        I_p_wave_areas_5_pulmonary_veins_dev_avg_min_3SD_index = i;
    end
end
I_p_wave_areas_5_pulmonary_veins_dev_avg_pls_3SD = area_5(length(I_p_wave_areas_5_pulmonary_veins_sorted),1);
for i = 1:length(I_p_wave_areas_5_pulmonary_veins_sorted)
    I_p_wave_areas_5_pulmonary_veins_dev_avg_pls_3SD(i) = I_p_wave_areas_5_pulmonary_veins_sorted(i) - (I_p_wave_areas_5_pulmonary_veins_average+3*I_p_wave_areas_5_pulmonary_veins_SD);
end
I_p_wave_areas_5_pulmonary_veins_dev_avg_pls_3SD = abs(I_p_wave_areas_5_pulmonary_veins_dev_avg_pls_3SD);
for i = 1:length(I_p_wave_areas_5_pulmonary_veins_sorted)
    if I_p_wave_areas_5_pulmonary_veins_dev_avg_pls_3SD(i) == min(I_p_wave_areas_5_pulmonary_veins_dev_avg_pls_3SD)
        I_p_wave_areas_5_pulmonary_veins_dev_avg_pls_3SD_index = i;
    end
end
I_p_wave_areas_5_pulmonary_veins_tot_number_of_indices = length(I_p_wave_areas_5_pulmonary_veins_sorted);
%This is the approximation of the tot "area under curve." 
I_p_wave_areas_5_pulmonary_veins_quotient_within_1SD = (I_p_wave_areas_5_pulmonary_veins_dev_avg_pls_1SD_index - I_p_wave_areas_5_pulmonary_veinsdev_avg_min_1SD_index)/I_p_wave_areas_5_pulmonary_veins_tot_number_of_indices;
I_p_wave_areas_5_pulmonary_veins_quotient_within_2SD = (I_p_wave_areas_5_pulmonary_veins_dev_avg_pls_2SD_index - I_p_wave_areas_5_dev_avg_min_2SD_index)/I_p_wave_areas_5_pulmonary_veins_tot_number_of_indices;
I_p_wave_areas_5_pulmonary_veins_quotient_within_3SD = (I_p_wave_areas_5_pulmonary_veins_dev_avg_pls_3SD_index - I_p_wave_areas_5_pulmonary_veins_dev_avg_min_3SD_index)/I_p_wave_areas_5_pulmonary_veins_tot_number_of_indices;
if (I_p_wave_areas_5_pulmonary_veins_quotient_within_1SD > 0.66 && I_p_wave_areas_5_pulmonary_veins_quotient_within_1SD < 0.70) && (I_p_wave_areas_5_pulmonary_veins_quotient_within_2SD > 0.93 && I_p_wave_areas_5_pulmonary_veins_quotient_within_2SD < 0.97) && (I_p_wave_areas_5_pulmonary_veins_quotient_within_3SD > 0.98 && I_p_wave_areas_5_pulmonary_veins_quotient_within_3SD < 1)
    isNormal(11,9) = 1;
end
if I_p_wave_areas_5_pulmonary_veins_quotient_within_1SD == 0
    isSpike(11,9) = 1;
end
I_p_wave_areas_5_mitral_annulus_sorted = sort(I_p_wave_areas_5_mitral_annulus);
I_p_wave_areas_5_mitral_annulus_dev_avg_min_1SD = area_5(length(I_p_wave_areas_5_mitral_annulus_sorted),1);
for i = 1:length(I_p_wave_areas_5_mitral_annulus_sorted)
    I_p_wave_areas_5_mitral_annulus_dev_avg_min_1SD(i) = I_p_wave_areas_5_mitral_annulus_sorted(i) - (I_p_wave_areas_5_mitral_annulus_average-I_p_wave_areas_5_mitral_annulus_SD);
end
I_p_wave_areas_5_mitral_annulus_dev_avg_min_1SD = abs(I_p_wave_areas_5_mitral_annulus_dev_avg_min_1SD);
for i = 1:length(I_p_wave_areas_5_mitral_annulus_sorted)
    if I_p_wave_areas_5_mitral_annulus_dev_avg_min_1SD(i) == min(I_p_wave_areas_5_mitral_annulus_dev_avg_min_1SD)
        I_p_wave_areas_5_mitral_annulus_dev_avg_min_1SD_index = i;
    end
end
I_p_wave_areas_5_mitral_annulus_dev_avg_pls_1SD = area_5(length(I_p_wave_areas_5_mitral_annulus_sorted),1);
for i = 1:length(I_p_wave_areas_5_mitral_annulus_sorted)
    I_p_wave_areas_5_mitral_annulus_dev_avg_pls_1SD(i) = I_p_wave_areas_5_mitral_annulus_sorted(i) - (I_p_wave_areas_5_mitral_annulus_average+I_p_wave_areas_5_mitral_annulus_SD);
end
I_p_wave_areas_5_mitral_annulus_dev_avg_pls_1SD = abs(I_p_wave_areas_5_mitral_annulus_dev_avg_pls_1SD);
for i = 1:length(I_p_wave_areas_5_mitral_annulus_sorted)
    if I_p_wave_areas_5_mitral_annulus_dev_avg_pls_1SD(i) == min(I_p_wave_areas_5_mitral_annulus_dev_avg_pls_1SD)
        I_p_wave_areas_5_mitral_annulus_dev_avg_pls_1SD_index = i;
    end
end
I_p_wave_areas_5_mitral_annulus_dev_avg_min_2SD = area_5(length(I_p_wave_areas_5_mitral_annulus_sorted),1);
for i = 1:length(I_p_wave_areas_5_mitral_annulus_sorted)
    I_p_wave_areas_5_mitral_annulus_dev_avg_min_2SD(i) = I_p_wave_areas_5_mitral_annulus_sorted(i) - (I_p_wave_areas_5_mitral_annulus_average-2*I_p_wave_areas_5_mitral_annulus_SD);
end
I_p_wave_areas_5_mitral_annulus_dev_avg_min_2SD = abs(I_p_wave_areas_5_mitral_annulus_dev_avg_min_2SD);
for i = 1:length(I_p_wave_areas_5_mitral_annulus_sorted)
    if I_p_wave_areas_5_mitral_annulus_dev_avg_min_2SD(i) == min(I_p_wave_areas_5_mitral_annulus_dev_avg_min_2SD)
        I_p_wave_areas_5_dev_avg_min_2SD_index = i;
    end
end
I_p_wave_areas_5_mitral_annulus_dev_avg_pls_2SD = area_5(length(I_p_wave_areas_5_mitral_annulus_sorted),1);
for i = 1:length(I_p_wave_areas_5_mitral_annulus_sorted)
    I_p_wave_areas_5_mitral_annulus_dev_avg_pls_2SD(i) = I_p_wave_areas_5_mitral_annulus_sorted(i) - (I_p_wave_areas_5_mitral_annulus_average+2*I_p_wave_areas_5_mitral_annulus_SD);
end
I_p_wave_areas_5_mitral_annulus_dev_avg_pls_2SD = abs(I_p_wave_areas_5_mitral_annulus_dev_avg_pls_2SD);
for i = 1:length(I_p_wave_areas_5_mitral_annulus_sorted)
    if I_p_wave_areas_5_mitral_annulus_dev_avg_pls_2SD(i) == min(I_p_wave_areas_5_mitral_annulus_dev_avg_pls_2SD)
        I_p_wave_areas_5_mitral_annulus_dev_avg_pls_2SD_index = i;
    end
end
I_p_wave_areas_5_mitral_annulus_dev_avg_min_3SD = area_5(length(I_p_wave_areas_5_mitral_annulus_sorted),1);
for i = 1:length(I_p_wave_areas_5_mitral_annulus_sorted)
    I_p_wave_areas_5_mitral_annulus_dev_avg_min_3SD(i) = I_p_wave_areas_5_mitral_annulus_sorted(i) - (I_p_wave_areas_5_mitral_annulus_average-3*I_p_wave_areas_5_mitral_annulus_SD);
end
I_p_wave_areas_5_mitral_annulus_dev_avg_min_3SD = abs(I_p_wave_areas_5_mitral_annulus_dev_avg_min_3SD);
for i = 1:length(I_p_wave_areas_5_mitral_annulus_sorted)
    if I_p_wave_areas_5_mitral_annulus_dev_avg_min_3SD(i) == min(I_p_wave_areas_5_mitral_annulus_dev_avg_min_3SD)
        I_p_wave_areas_5_mitral_annulus_dev_avg_min_3SD_index = i;
    end
end
I_p_wave_areas_5_mitral_annulus_dev_avg_pls_3SD = area_5(length(I_p_wave_areas_5_mitral_annulus_sorted),1);
for i = 1:length(I_p_wave_areas_5_mitral_annulus_sorted)
    I_p_wave_areas_5_mitral_annulus_dev_avg_pls_3SD(i) = I_p_wave_areas_5_mitral_annulus_sorted(i) - (I_p_wave_areas_5_mitral_annulus_average+3*I_p_wave_areas_5_mitral_annulus_SD);
end
I_p_wave_areas_5_mitral_annulus_dev_avg_pls_3SD = abs(I_p_wave_areas_5_mitral_annulus_dev_avg_pls_3SD);
for i = 1:length(I_p_wave_areas_5_mitral_annulus_sorted)
    if I_p_wave_areas_5_mitral_annulus_dev_avg_pls_3SD(i) == min(I_p_wave_areas_5_mitral_annulus_dev_avg_pls_3SD)
        I_p_wave_areas_5_mitral_annulus_dev_avg_pls_3SD_index = i;
    end
end
I_p_wave_areas_5_mitral_annulus_tot_number_of_indices = length(I_p_wave_areas_5_mitral_annulus_sorted);
%This is the approximation of the tot "area under curve." 
I_p_wave_areas_5_mitral_annulus_quotient_within_1SD = (I_p_wave_areas_5_mitral_annulus_dev_avg_pls_1SD_index - I_p_wave_areas_5_mitral_annulusdev_avg_min_1SD_index)/I_p_wave_areas_5_mitral_annulus_tot_number_of_indices;
I_p_wave_areas_5_mitral_annulus_quotient_within_2SD = (I_p_wave_areas_5_mitral_annulus_dev_avg_pls_2SD_index - I_p_wave_areas_5_dev_avg_min_2SD_index)/I_p_wave_areas_5_mitral_annulus_tot_number_of_indices;
I_p_wave_areas_5_mitral_annulus_quotient_within_3SD = (I_p_wave_areas_5_mitral_annulus_dev_avg_pls_3SD_index - I_p_wave_areas_5_mitral_annulus_dev_avg_min_3SD_index)/I_p_wave_areas_5_mitral_annulus_tot_number_of_indices;
if (I_p_wave_areas_5_mitral_annulus_quotient_within_1SD > 0.66 && I_p_wave_areas_5_mitral_annulus_quotient_within_1SD < 0.70) && (I_p_wave_areas_5_mitral_annulus_quotient_within_2SD > 0.93 && I_p_wave_areas_5_mitral_annulus_quotient_within_2SD < 0.97) && (I_p_wave_areas_5_mitral_annulus_quotient_within_3SD > 0.98 && I_p_wave_areas_5_mitral_annulus_quotient_within_3SD < 1)
    isNormal(11,10) = 1;
end
if I_p_wave_areas_5_mitral_annulus_quotient_within_1SD == 0
    isSpike(11,10) = 1;
end
I_p_wave_areas_5_CS_body_sorted = sort(I_p_wave_areas_5_CS_body);
I_p_wave_areas_5_CS_body_dev_avg_min_1SD = area_5(length(I_p_wave_areas_5_CS_body_sorted),1);
for i = 1:length(I_p_wave_areas_5_CS_body_sorted)
    I_p_wave_areas_5_CS_body_dev_avg_min_1SD(i) = I_p_wave_areas_5_CS_body_sorted(i) - (I_p_wave_areas_5_CS_body_average-I_p_wave_areas_5_CS_body_SD);
end
I_p_wave_areas_5_CS_body_dev_avg_min_1SD = abs(I_p_wave_areas_5_CS_body_dev_avg_min_1SD);
for i = 1:length(I_p_wave_areas_5_CS_body_sorted)
    if I_p_wave_areas_5_CS_body_dev_avg_min_1SD(i) == min(I_p_wave_areas_5_CS_body_dev_avg_min_1SD)
        I_p_wave_areas_5_CS_body_dev_avg_min_1SD_index = i;
    end
end
I_p_wave_areas_5_CS_body_dev_avg_pls_1SD = area_5(length(I_p_wave_areas_5_CS_body_sorted),1);
for i = 1:length(I_p_wave_areas_5_CS_body_sorted)
    I_p_wave_areas_5_CS_body_dev_avg_pls_1SD(i) = I_p_wave_areas_5_CS_body_sorted(i) - (I_p_wave_areas_5_CS_body_average+I_p_wave_areas_5_CS_body_SD);
end
I_p_wave_areas_5_CS_body_dev_avg_pls_1SD = abs(I_p_wave_areas_5_CS_body_dev_avg_pls_1SD);
for i = 1:length(I_p_wave_areas_5_CS_body_sorted)
    if I_p_wave_areas_5_CS_body_dev_avg_pls_1SD(i) == min(I_p_wave_areas_5_CS_body_dev_avg_pls_1SD)
        I_p_wave_areas_5_CS_body_dev_avg_pls_1SD_index = i;
    end
end
I_p_wave_areas_5_CS_body_dev_avg_min_2SD = area_5(length(I_p_wave_areas_5_CS_body_sorted),1);
for i = 1:length(I_p_wave_areas_5_CS_body_sorted)
    I_p_wave_areas_5_CS_body_dev_avg_min_2SD(i) = I_p_wave_areas_5_CS_body_sorted(i) - (I_p_wave_areas_5_CS_body_average-2*I_p_wave_areas_5_CS_body_SD);
end
I_p_wave_areas_5_CS_body_dev_avg_min_2SD = abs(I_p_wave_areas_5_CS_body_dev_avg_min_2SD);
for i = 1:length(I_p_wave_areas_5_CS_body_sorted)
    if I_p_wave_areas_5_CS_body_dev_avg_min_2SD(i) == min(I_p_wave_areas_5_CS_body_dev_avg_min_2SD)
        I_p_wave_areas_5_dev_avg_min_2SD_index = i;
    end
end
I_p_wave_areas_5_CS_body_dev_avg_pls_2SD = area_5(length(I_p_wave_areas_5_CS_body_sorted),1);
for i = 1:length(I_p_wave_areas_5_CS_body_sorted)
    I_p_wave_areas_5_CS_body_dev_avg_pls_2SD(i) = I_p_wave_areas_5_CS_body_sorted(i) - (I_p_wave_areas_5_CS_body_average+2*I_p_wave_areas_5_CS_body_SD);
end
I_p_wave_areas_5_CS_body_dev_avg_pls_2SD = abs(I_p_wave_areas_5_CS_body_dev_avg_pls_2SD);
for i = 1:length(I_p_wave_areas_5_CS_body_sorted)
    if I_p_wave_areas_5_CS_body_dev_avg_pls_2SD(i) == min(I_p_wave_areas_5_CS_body_dev_avg_pls_2SD)
        I_p_wave_areas_5_CS_body_dev_avg_pls_2SD_index = i;
    end
end
I_p_wave_areas_5_CS_body_dev_avg_min_3SD = area_5(length(I_p_wave_areas_5_CS_body_sorted),1);
for i = 1:length(I_p_wave_areas_5_CS_body_sorted)
    I_p_wave_areas_5_CS_body_dev_avg_min_3SD(i) = I_p_wave_areas_5_CS_body_sorted(i) - (I_p_wave_areas_5_CS_body_average-3*I_p_wave_areas_5_CS_body_SD);
end
I_p_wave_areas_5_CS_body_dev_avg_min_3SD = abs(I_p_wave_areas_5_CS_body_dev_avg_min_3SD);
for i = 1:length(I_p_wave_areas_5_CS_body_sorted)
    if I_p_wave_areas_5_CS_body_dev_avg_min_3SD(i) == min(I_p_wave_areas_5_CS_body_dev_avg_min_3SD)
        I_p_wave_areas_5_CS_body_dev_avg_min_3SD_index = i;
    end
end
I_p_wave_areas_5_CS_body_dev_avg_pls_3SD = area_5(length(I_p_wave_areas_5_CS_body_sorted),1);
for i = 1:length(I_p_wave_areas_5_CS_body_sorted)
    I_p_wave_areas_5_CS_body_dev_avg_pls_3SD(i) = I_p_wave_areas_5_CS_body_sorted(i) - (I_p_wave_areas_5_CS_body_average+3*I_p_wave_areas_5_CS_body_SD);
end
I_p_wave_areas_5_CS_body_dev_avg_pls_3SD = abs(I_p_wave_areas_5_CS_body_dev_avg_pls_3SD);
for i = 1:length(I_p_wave_areas_5_CS_body_sorted)
    if I_p_wave_areas_5_CS_body_dev_avg_pls_3SD(i) == min(I_p_wave_areas_5_CS_body_dev_avg_pls_3SD)
        I_p_wave_areas_5_CS_body_dev_avg_pls_3SD_index = i;
    end
end
I_p_wave_areas_5_CS_body_tot_number_of_indices = length(I_p_wave_areas_5_CS_body_sorted);
%This is the approximation of the tot "area under curve." 
I_p_wave_areas_5_CS_body_quotient_within_1SD = (I_p_wave_areas_5_CS_body_dev_avg_pls_1SD_index - I_p_wave_areas_5_CS_bodydev_avg_min_1SD_index)/I_p_wave_areas_5_CS_body_tot_number_of_indices;
I_p_wave_areas_5_CS_body_quotient_within_2SD = (I_p_wave_areas_5_CS_body_dev_avg_pls_2SD_index - I_p_wave_areas_5_dev_avg_min_2SD_index)/I_p_wave_areas_5_CS_body_tot_number_of_indices;
I_p_wave_areas_5_CS_body_quotient_within_3SD = (I_p_wave_areas_5_CS_body_dev_avg_pls_3SD_index - I_p_wave_areas_5_CS_body_dev_avg_min_3SD_index)/I_p_wave_areas_5_CS_body_tot_number_of_indices;
if (I_p_wave_areas_5_CS_body_quotient_within_1SD > 0.66 && I_p_wave_areas_5_CS_body_quotient_within_1SD < 0.70) && (I_p_wave_areas_5_CS_body_quotient_within_2SD > 0.93 && I_p_wave_areas_5_CS_body_quotient_within_2SD < 0.97) && (I_p_wave_areas_5_CS_body_quotient_within_3SD > 0.98 && I_p_wave_areas_5_CS_body_quotient_within_3SD < 1)
    isNormal(11,11) = 1;
end
if I_p_wave_areas_5_CS_body_quotient_within_1SD == 0
    isSpike(11,11) = 1;
end
I_p_wave_areas_5_left_septum_sorted = sort(I_p_wave_areas_5_left_septum);
I_p_wave_areas_5_left_septum_dev_avg_min_1SD = area_5(length(I_p_wave_areas_5_left_septum_sorted),1);
for i = 1:length(I_p_wave_areas_5_left_septum_sorted)
    I_p_wave_areas_5_left_septum_dev_avg_min_1SD(i) = I_p_wave_areas_5_left_septum_sorted(i) - (I_p_wave_areas_5_left_septum_average-I_p_wave_areas_5_left_septum_SD);
end
I_p_wave_areas_5_left_septum_dev_avg_min_1SD = abs(I_p_wave_areas_5_left_septum_dev_avg_min_1SD);
for i = 1:length(I_p_wave_areas_5_left_septum_sorted)
    if I_p_wave_areas_5_left_septum_dev_avg_min_1SD(i) == min(I_p_wave_areas_5_left_septum_dev_avg_min_1SD)
        I_p_wave_areas_5_left_septum_dev_avg_min_1SD_index = i;
    end
end
I_p_wave_areas_5_left_septum_dev_avg_pls_1SD = area_5(length(I_p_wave_areas_5_left_septum_sorted),1);
for i = 1:length(I_p_wave_areas_5_left_septum_sorted)
    I_p_wave_areas_5_left_septum_dev_avg_pls_1SD(i) = I_p_wave_areas_5_left_septum_sorted(i) - (I_p_wave_areas_5_left_septum_average+I_p_wave_areas_5_left_septum_SD);
end
I_p_wave_areas_5_left_septum_dev_avg_pls_1SD = abs(I_p_wave_areas_5_left_septum_dev_avg_pls_1SD);
for i = 1:length(I_p_wave_areas_5_left_septum_sorted)
    if I_p_wave_areas_5_left_septum_dev_avg_pls_1SD(i) == min(I_p_wave_areas_5_left_septum_dev_avg_pls_1SD)
        I_p_wave_areas_5_left_septum_dev_avg_pls_1SD_index = i;
    end
end
I_p_wave_areas_5_left_septum_dev_avg_min_2SD = area_5(length(I_p_wave_areas_5_left_septum_sorted),1);
for i = 1:length(I_p_wave_areas_5_left_septum_sorted)
    I_p_wave_areas_5_left_septum_dev_avg_min_2SD(i) = I_p_wave_areas_5_left_septum_sorted(i) - (I_p_wave_areas_5_left_septum_average-2*I_p_wave_areas_5_left_septum_SD);
end
I_p_wave_areas_5_left_septum_dev_avg_min_2SD = abs(I_p_wave_areas_5_left_septum_dev_avg_min_2SD);
for i = 1:length(I_p_wave_areas_5_left_septum_sorted)
    if I_p_wave_areas_5_left_septum_dev_avg_min_2SD(i) == min(I_p_wave_areas_5_left_septum_dev_avg_min_2SD)
        I_p_wave_areas_5_dev_avg_min_2SD_index = i;
    end
end
I_p_wave_areas_5_left_septum_dev_avg_pls_2SD = area_5(length(I_p_wave_areas_5_left_septum_sorted),1);
for i = 1:length(I_p_wave_areas_5_left_septum_sorted)
    I_p_wave_areas_5_left_septum_dev_avg_pls_2SD(i) = I_p_wave_areas_5_left_septum_sorted(i) - (I_p_wave_areas_5_left_septum_average+2*I_p_wave_areas_5_left_septum_SD);
end
I_p_wave_areas_5_left_septum_dev_avg_pls_2SD = abs(I_p_wave_areas_5_left_septum_dev_avg_pls_2SD);
for i = 1:length(I_p_wave_areas_5_left_septum_sorted)
    if I_p_wave_areas_5_left_septum_dev_avg_pls_2SD(i) == min(I_p_wave_areas_5_left_septum_dev_avg_pls_2SD)
        I_p_wave_areas_5_left_septum_dev_avg_pls_2SD_index = i;
    end
end
I_p_wave_areas_5_left_septum_dev_avg_min_3SD = area_5(length(I_p_wave_areas_5_left_septum_sorted),1);
for i = 1:length(I_p_wave_areas_5_left_septum_sorted)
    I_p_wave_areas_5_left_septum_dev_avg_min_3SD(i) = I_p_wave_areas_5_left_septum_sorted(i) - (I_p_wave_areas_5_left_septum_average-3*I_p_wave_areas_5_left_septum_SD);
end
I_p_wave_areas_5_left_septum_dev_avg_min_3SD = abs(I_p_wave_areas_5_left_septum_dev_avg_min_3SD);
for i = 1:length(I_p_wave_areas_5_left_septum_sorted)
    if I_p_wave_areas_5_left_septum_dev_avg_min_3SD(i) == min(I_p_wave_areas_5_left_septum_dev_avg_min_3SD)
        I_p_wave_areas_5_left_septum_dev_avg_min_3SD_index = i;
    end
end
I_p_wave_areas_5_left_septum_dev_avg_pls_3SD = area_5(length(I_p_wave_areas_5_left_septum_sorted),1);
for i = 1:length(I_p_wave_areas_5_left_septum_sorted)
    I_p_wave_areas_5_left_septum_dev_avg_pls_3SD(i) = I_p_wave_areas_5_left_septum_sorted(i) - (I_p_wave_areas_5_left_septum_average+3*I_p_wave_areas_5_left_septum_SD);
end
I_p_wave_areas_5_left_septum_dev_avg_pls_3SD = abs(I_p_wave_areas_5_left_septum_dev_avg_pls_3SD);
for i = 1:length(I_p_wave_areas_5_left_septum_sorted)
    if I_p_wave_areas_5_left_septum_dev_avg_pls_3SD(i) == min(I_p_wave_areas_5_left_septum_dev_avg_pls_3SD)
        I_p_wave_areas_5_left_septum_dev_avg_pls_3SD_index = i;
    end
end
I_p_wave_areas_5_left_septum_tot_number_of_indices = length(I_p_wave_areas_5_left_septum_sorted);
%This is the approximation of the tot "area under curve." 
I_p_wave_areas_5_left_septum_quotient_within_1SD = (I_p_wave_areas_5_left_septum_dev_avg_pls_1SD_index - I_p_wave_areas_5_left_septumdev_avg_min_1SD_index)/I_p_wave_areas_5_left_septum_tot_number_of_indices;
I_p_wave_areas_5_left_septum_quotient_within_2SD = (I_p_wave_areas_5_left_septum_dev_avg_pls_2SD_index - I_p_wave_areas_5_dev_avg_min_2SD_index)/I_p_wave_areas_5_left_septum_tot_number_of_indices;
I_p_wave_areas_5_left_septum_quotient_within_3SD = (I_p_wave_areas_5_left_septum_dev_avg_pls_3SD_index - I_p_wave_areas_5_left_septum_dev_avg_min_3SD_index)/I_p_wave_areas_5_left_septum_tot_number_of_indices;
if (I_p_wave_areas_5_left_septum_quotient_within_1SD > 0.66 && I_p_wave_areas_5_left_septum_quotient_within_1SD < 0.70) && (I_p_wave_areas_5_left_septum_quotient_within_2SD > 0.93 && I_p_wave_areas_5_left_septum_quotient_within_2SD < 0.97) && (I_p_wave_areas_5_left_septum_quotient_within_3SD > 0.98 && I_p_wave_areas_5_left_septum_quotient_within_3SD < 1)
    isNormal(11,12) = 1;
end
if I_p_wave_areas_5_left_septum_quotient_within_1SD == 0
    isSpike(11,12) = 1;
end
I_p_wave_areas_5_left_atrial_appendage_sorted = sort(I_p_wave_areas_5_left_atrial_appendage);
I_p_wave_areas_5_left_atrial_appendage_dev_avg_min_1SD = area_5(length(I_p_wave_areas_5_left_atrial_appendage_sorted),1);
for i = 1:length(I_p_wave_areas_5_left_atrial_appendage_sorted)
    I_p_wave_areas_5_left_atrial_appendage_dev_avg_min_1SD(i) = I_p_wave_areas_5_left_atrial_appendage_sorted(i) - (I_p_wave_areas_5_left_atrial_appendage_average-I_p_wave_areas_5_left_atrial_appendage_SD);
end
I_p_wave_areas_5_left_atrial_appendage_dev_avg_min_1SD = abs(I_p_wave_areas_5_left_atrial_appendage_dev_avg_min_1SD);
for i = 1:length(I_p_wave_areas_5_left_atrial_appendage_sorted)
    if I_p_wave_areas_5_left_atrial_appendage_dev_avg_min_1SD(i) == min(I_p_wave_areas_5_left_atrial_appendage_dev_avg_min_1SD)
        I_p_wave_areas_5_left_atrial_appendage_dev_avg_min_1SD_index = i;
    end
end
I_p_wave_areas_5_left_atrial_appendage_dev_avg_pls_1SD = area_5(length(I_p_wave_areas_5_left_atrial_appendage_sorted),1);
for i = 1:length(I_p_wave_areas_5_left_atrial_appendage_sorted)
    I_p_wave_areas_5_left_atrial_appendage_dev_avg_pls_1SD(i) = I_p_wave_areas_5_left_atrial_appendage_sorted(i) - (I_p_wave_areas_5_left_atrial_appendage_average+I_p_wave_areas_5_left_atrial_appendage_SD);
end
I_p_wave_areas_5_left_atrial_appendage_dev_avg_pls_1SD = abs(I_p_wave_areas_5_left_atrial_appendage_dev_avg_pls_1SD);
for i = 1:length(I_p_wave_areas_5_left_atrial_appendage_sorted)
    if I_p_wave_areas_5_left_atrial_appendage_dev_avg_pls_1SD(i) == min(I_p_wave_areas_5_left_atrial_appendage_dev_avg_pls_1SD)
        I_p_wave_areas_5_left_atrial_appendage_dev_avg_pls_1SD_index = i;
    end
end
I_p_wave_areas_5_left_atrial_appendage_dev_avg_min_2SD = area_5(length(I_p_wave_areas_5_left_atrial_appendage_sorted),1);
for i = 1:length(I_p_wave_areas_5_left_atrial_appendage_sorted)
    I_p_wave_areas_5_left_atrial_appendage_dev_avg_min_2SD(i) = I_p_wave_areas_5_left_atrial_appendage_sorted(i) - (I_p_wave_areas_5_left_atrial_appendage_average-2*I_p_wave_areas_5_left_atrial_appendage_SD);
end
I_p_wave_areas_5_left_atrial_appendage_dev_avg_min_2SD = abs(I_p_wave_areas_5_left_atrial_appendage_dev_avg_min_2SD);
for i = 1:length(I_p_wave_areas_5_left_atrial_appendage_sorted)
    if I_p_wave_areas_5_left_atrial_appendage_dev_avg_min_2SD(i) == min(I_p_wave_areas_5_left_atrial_appendage_dev_avg_min_2SD)
        I_p_wave_areas_5_dev_avg_min_2SD_index = i;
    end
end
I_p_wave_areas_5_left_atrial_appendage_dev_avg_pls_2SD = area_5(length(I_p_wave_areas_5_left_atrial_appendage_sorted),1);
for i = 1:length(I_p_wave_areas_5_left_atrial_appendage_sorted)
    I_p_wave_areas_5_left_atrial_appendage_dev_avg_pls_2SD(i) = I_p_wave_areas_5_left_atrial_appendage_sorted(i) - (I_p_wave_areas_5_left_atrial_appendage_average+2*I_p_wave_areas_5_left_atrial_appendage_SD);
end
I_p_wave_areas_5_left_atrial_appendage_dev_avg_pls_2SD = abs(I_p_wave_areas_5_left_atrial_appendage_dev_avg_pls_2SD);
for i = 1:length(I_p_wave_areas_5_left_atrial_appendage_sorted)
    if I_p_wave_areas_5_left_atrial_appendage_dev_avg_pls_2SD(i) == min(I_p_wave_areas_5_left_atrial_appendage_dev_avg_pls_2SD)
        I_p_wave_areas_5_left_atrial_appendage_dev_avg_pls_2SD_index = i;
    end
end
I_p_wave_areas_5_left_atrial_appendage_dev_avg_min_3SD = area_5(length(I_p_wave_areas_5_left_atrial_appendage_sorted),1);
for i = 1:length(I_p_wave_areas_5_left_atrial_appendage_sorted)
    I_p_wave_areas_5_left_atrial_appendage_dev_avg_min_3SD(i) = I_p_wave_areas_5_left_atrial_appendage_sorted(i) - (I_p_wave_areas_5_left_atrial_appendage_average-3*I_p_wave_areas_5_left_atrial_appendage_SD);
end
I_p_wave_areas_5_left_atrial_appendage_dev_avg_min_3SD = abs(I_p_wave_areas_5_left_atrial_appendage_dev_avg_min_3SD);
for i = 1:length(I_p_wave_areas_5_left_atrial_appendage_sorted)
    if I_p_wave_areas_5_left_atrial_appendage_dev_avg_min_3SD(i) == min(I_p_wave_areas_5_left_atrial_appendage_dev_avg_min_3SD)
        I_p_wave_areas_5_left_atrial_appendage_dev_avg_min_3SD_index = i;
    end
end
I_p_wave_areas_5_left_atrial_appendage_dev_avg_pls_3SD = area_5(length(I_p_wave_areas_5_left_atrial_appendage_sorted),1);
for i = 1:length(I_p_wave_areas_5_left_atrial_appendage_sorted)
    I_p_wave_areas_5_left_atrial_appendage_dev_avg_pls_3SD(i) = I_p_wave_areas_5_left_atrial_appendage_sorted(i) - (I_p_wave_areas_5_left_atrial_appendage_average+3*I_p_wave_areas_5_left_atrial_appendage_SD);
end
I_p_wave_areas_5_left_atrial_appendage_dev_avg_pls_3SD = abs(I_p_wave_areas_5_left_atrial_appendage_dev_avg_pls_3SD);
for i = 1:length(I_p_wave_areas_5_left_atrial_appendage_sorted)
    if I_p_wave_areas_5_left_atrial_appendage_dev_avg_pls_3SD(i) == min(I_p_wave_areas_5_left_atrial_appendage_dev_avg_pls_3SD)
        I_p_wave_areas_5_left_atrial_appendage_dev_avg_pls_3SD_index = i;
    end
end
I_p_wave_areas_5_left_atrial_appendage_tot_number_of_indices = length(I_p_wave_areas_5_left_atrial_appendage_sorted);
%This is the approximation of the tot "area under curve." 
I_p_wave_areas_5_left_atrial_appendage_quotient_within_1SD = (I_p_wave_areas_5_left_atrial_appendage_dev_avg_pls_1SD_index - I_p_wave_areas_5_left_atrial_appendagedev_avg_min_1SD_index)/I_p_wave_areas_5_left_atrial_appendage_tot_number_of_indices;
I_p_wave_areas_5_left_atrial_appendage_quotient_within_2SD = (I_p_wave_areas_5_left_atrial_appendage_dev_avg_pls_2SD_index - I_p_wave_areas_5_dev_avg_min_2SD_index)/I_p_wave_areas_5_left_atrial_appendage_tot_number_of_indices;
I_p_wave_areas_5_left_atrial_appendage_quotient_within_3SD = (I_p_wave_areas_5_left_atrial_appendage_dev_avg_pls_3SD_index - I_p_wave_areas_5_left_atrial_appendage_dev_avg_min_3SD_index)/I_p_wave_areas_5_left_atrial_appendage_tot_number_of_indices;
if (I_p_wave_areas_5_left_atrial_appendage_quotient_within_1SD > 0.66 && I_p_wave_areas_5_left_atrial_appendage_quotient_within_1SD < 0.70) && (I_p_wave_areas_5_left_atrial_appendage_quotient_within_2SD > 0.93 && I_p_wave_areas_5_left_atrial_appendage_quotient_within_2SD < 0.97) && (I_p_wave_areas_5_left_atrial_appendage_quotient_within_3SD > 0.98 && I_p_wave_areas_5_left_atrial_appendage_quotient_within_3SD < 1)
    isNormal(11,13) = 1;
end
if I_p_wave_areas_5_left_atrial_appendage_quotient_within_1SD == 0
    isSpike(11,13) = 1;
end
number_of_p_waves = length(p_wave_start_indices);
for i = 1:number_of_p_waves
    if I_num_zeros_SA_node_SD == 0 && I_number_of_zeros_array(i) == I_num_zeros_SA_node_average 
        I_num_zeros_z_scores_SA_node(i) = 0;
    else
        if I_num_zeros_SA_node_SD == 0 && ~(I_number_of_zeros_array(i) == I_num_zeros_SA_node_average)
            I_num_zeros_z_scores_SA_node(i) = 10*(abs(I_num_zeros_z_scores_SA_node(i)-I_num_zeros_SA_node_average));
        end
    end
end
for i = 1:number_of_p_waves
    if I_maxima_SA_node_SD == 0 && I_maxima(i) == I_maxima_SA_node_average 
        I_maxima_z_scores_SA_node(i) = 0;
    else
        if I_maxima_SA_node_SD == 0 && ~(I_maxima(i) == I_maxima_SA_node_average)
            I_maxima_z_scores_SA_node(i) = 10*(abs(I_maxima(i)-I_maxima_SA_node_average));
        end
    end
end
for i = 1:number_of_p_waves
    if I_minima_SA_node_SD == 0 && I_minima(i) == I_minima_SA_node_average 
        I_minima_z_scores_SA_node(i) = 0;
    else
        if I_minima_SA_node_SD == 0 && ~(I_minima(i) == I_minima_SA_node_average)
            I_minima_z_scores_SA_node(i) = 10*(abs(I_minima(i)-I_minima_SA_node_average));
        end
    end
end
for i = 1:number_of_p_waves
    if I_sym_inds_SA_node_SD == 0 && I_sym_inds(i) == I_sym_inds_SA_node_average 
        I_sym_inds_z_scores_SA_node(i) = 0;
    else
        if I_sym_inds_SA_node_SD == 0 && ~(I_sym_inds(i) == I_sym_inds_SA_node_average)
            I_sym_inds_z_scores_SA_node(i) = 10*(abs(I_sym_inds(i)-I_sym_inds_SA_node_average));
        end
    end
end
for i = 1:number_of_p_waves
    if I_concavs_SA_node_SD == 0 && I_concavs(i) == I_concavs_SA_node_average 
        I_concavs_z_scores_SA_node(i) = 0;
    else
        if I_concavs_SA_node_SD == 0 && ~(I_concavs(i) == I_concavs_SA_node_average)
            I_concavs_z_scores_SA_node(i) = 10*(abs(I_concavs(i)-I_concavs_SA_node_average));
        end
    end
end
for i = 1:number_of_p_waves
    if I_mdslpes_SA_node_SD == 0 && I_mdslpes(i) == I_mdslpes_SA_node_average 
        I_mdslpes_z_scores_SA_node(i) = 0;
    else
        if I_mdslpes_SA_node_SD == 0 && ~(I_mdslpes(i) == I_mdslpes_SA_node_average)
            I_mdslpes_z_scores_SA_node(i) = 10*(abs(I_mdslpes(i)-I_mdslpes_SA_node_average));
        end
    end
end
for i = 1:number_of_p_waves
    if I_p_wave_areas_1_SA_node_SD == 0 && I_p_wave_areas_1(i) == I_p_wave_areas_1_SA_node_average 
        I_p_wave_areas_1_z_scores_SA_node(i) = 0;
    else
        if I_p_wave_areas_1_SA_node_SD == 0 && ~(I_p_wave_areas_1(i) == I_p_wave_areas_1_SA_node_average)
            I_p_wave_areas_1_z_scores_SA_node(i) = 10*(abs(I_p_wave_areas_1(i)-I_p_wave_areas_1_SA_node_average));
        end
    end
end
for i = 1:number_of_p_waves
    if I_p_wave_areas_2_SA_node_SD == 0 && I_p_wave_areas_2(i) == I_p_wave_areas_2_SA_node_average 
        I_p_wave_areas_2_z_scores_SA_node(i) = 0;
    else
        if I_p_wave_areas_2_SA_node_SD == 0 && ~(I_p_wave_areas_2(i) == I_p_wave_areas_2_SA_node_average)
            I_p_wave_areas_2_z_scores_SA_node(i) = 10*(abs(I_p_wave_areas_2(i)-I_p_wave_areas_2_SA_node_average));
        end
    end
end
for i = 1:number_of_p_waves
    if I_p_wave_areas_3_SA_node_SD == 0 && I_p_wave_areas_3(i) == I_p_wave_areas_3_SA_node_average 
        I_p_wave_areas_3_z_scores_SA_node(i) = 0;
    else
        if I_p_wave_areas_3_SA_node_SD == 0 && ~(I_p_wave_areas_3(i) == I_p_wave_areas_3_SA_node_average)
            I_p_wave_areas_3_z_scores_SA_node(i) = 10*(abs(I_p_wave_areas_3(i)-I_p_wave_areas_3_SA_node_average));
        end
    end
end
for i = 1:number_of_p_waves
    if I_p_wave_areas_4_SA_node_SD == 0 && I_p_wave_areas_4(i) == I_p_wave_areas_4_SA_node_average 
        I_p_wave_areas_4_z_scores_SA_node(i) = 0;
    else
        if I_p_wave_areas_4_SA_node_SD == 0 && ~(I_p_wave_areas_4(i) == I_p_wave_areas_4_SA_node_average)
            I_p_wave_areas_4_z_scores_SA_node(i) = 10*(abs(I_p_wave_areas_4(i)-I_p_wave_areas_4_SA_node_average));
        end
    end
end
for i = 1:number_of_p_waves
    if I_p_wave_areas_5_SA_node_SD == 0 && I_p_wave_areas_5(i) == I_p_wave_areas_5_SA_node_average 
        I_p_wave_areas_5_z_scores_SA_node(i) = 0;
    else
        if I_p_wave_areas_5_SA_node_SD == 0 && ~(I_p_wave_areas_5(i) == I_p_wave_areas_5_SA_node_average)
            I_p_wave_areas_5_z_scores_SA_node(i) = 10*(abs(I_p_wave_areas_5(i)-I_p_wave_areas_5_SA_node_average));
        end
    end
end
for i = 1:number_of_p_waves
    if I_num_zeros_crista_terminalis_SD == 0 && I_number_of_zeros_array(i) == I_num_zeros_crista_terminalis_average 
        I_num_zeros_z_scores_crista_terminalis(i) = 0;
    else
        if I_num_zeros_crista_terminalis_SD == 0 && ~(I_number_of_zeros_array(i) == I_num_zeros_crista_terminalis_average)
            I_num_zeros_z_scores_crista_terminalis(i) = 10*(abs(I_num_zeros_z_scores_crista_terminalis(i)-I_num_zeros_crista_terminalis_average));
        end
    end
end
for i = 1:number_of_p_waves
    if I_maxima_crista_terminalis_SD == 0 && I_maxima(i) == I_maxima_crista_terminalis_average 
        I_maxima_z_scores_crista_terminalis(i) = 0;
    else
        if I_maxima_crista_terminalis_SD == 0 && ~(I_maxima(i) == I_maxima_crista_terminalis_average)
            I_maxima_z_scores_crista_terminalis(i) = 10*(abs(I_maxima(i)-I_maxima_crista_terminalis_average));
        end
    end
end
for i = 1:number_of_p_waves
    if I_minima_crista_terminalis_SD == 0 && I_minima(i) == I_minima_crista_terminalis_average 
        I_minima_z_scores_crista_terminalis(i) = 0;
    else
        if I_minima_crista_terminalis_SD == 0 && ~(I_minima(i) == I_minima_crista_terminalis_average)
            I_minima_z_scores_crista_terminalis(i) = 10*(abs(I_minima(i)-I_minima_crista_terminalis_average));
        end
    end
end
for i = 1:number_of_p_waves
    if I_sym_inds_crista_terminalis_SD == 0 && I_sym_inds(i) == I_sym_inds_crista_terminalis_average 
        I_sym_inds_z_scores_crista_terminalis(i) = 0;
    else
        if I_sym_inds_crista_terminalis_SD == 0 && ~(I_sym_inds(i) == I_sym_inds_crista_terminalis_average)
            I_sym_inds_z_scores_crista_terminalis(i) = 10*(abs(I_sym_inds(i)-I_sym_inds_crista_terminalis_average));
        end
    end
end
for i = 1:number_of_p_waves
    if I_concavs_crista_terminalis_SD == 0 && I_concavs(i) == I_concavs_crista_terminalis_average 
        I_concavs_z_scores_crista_terminalis(i) = 0;
    else
        if I_concavs_crista_terminalis_SD == 0 && ~(I_concavs(i) == I_concavs_crista_terminalis_average)
            I_concavs_z_scores_crista_terminalis(i) = 10*(abs(I_concavs(i)-I_concavs_crista_terminalis_average));
        end
    end
end
for i = 1:number_of_p_waves
    if I_mdslpes_crista_terminalis_SD == 0 && I_mdslpes(i) == I_mdslpes_crista_terminalis_average 
        I_mdslpes_z_scores_crista_terminalis(i) = 0;
    else
        if I_mdslpes_crista_terminalis_SD == 0 && ~(I_mdslpes(i) == I_mdslpes_crista_terminalis_average)
            I_mdslpes_z_scores_crista_terminalis(i) = 10*(abs(I_mdslpes(i)-I_mdslpes_crista_terminalis_average));
        end
    end
end
for i = 1:number_of_p_waves
    if I_p_wave_areas_1_crista_terminalis_SD == 0 && I_p_wave_areas_1(i) == I_p_wave_areas_1_crista_terminalis_average 
        I_p_wave_areas_1_z_scores_crista_terminalis(i) = 0;
    else
        if I_p_wave_areas_1_crista_terminalis_SD == 0 && ~(I_p_wave_areas_1(i) == I_p_wave_areas_1_crista_terminalis_average)
            I_p_wave_areas_1_z_scores_crista_terminalis(i) = 10*(abs(I_p_wave_areas_1(i)-I_p_wave_areas_1_crista_terminalis_average));
        end
    end
end
for i = 1:number_of_p_waves
    if I_p_wave_areas_2_crista_terminalis_SD == 0 && I_p_wave_areas_2(i) == I_p_wave_areas_2_crista_terminalis_average 
        I_p_wave_areas_2_z_scores_crista_terminalis(i) = 0;
    else
        if I_p_wave_areas_2_crista_terminalis_SD == 0 && ~(I_p_wave_areas_2(i) == I_p_wave_areas_2_crista_terminalis_average)
            I_p_wave_areas_2_z_scores_crista_terminalis(i) = 10*(abs(I_p_wave_areas_2(i)-I_p_wave_areas_2_crista_terminalis_average));
        end
    end
end
for i = 1:number_of_p_waves
    if I_p_wave_areas_3_crista_terminalis_SD == 0 && I_p_wave_areas_3(i) == I_p_wave_areas_3_crista_terminalis_average 
        I_p_wave_areas_3_z_scores_crista_terminalis(i) = 0;
    else
        if I_p_wave_areas_3_crista_terminalis_SD == 0 && ~(I_p_wave_areas_3(i) == I_p_wave_areas_3_crista_terminalis_average)
            I_p_wave_areas_3_z_scores_crista_terminalis(i) = 10*(abs(I_p_wave_areas_3(i)-I_p_wave_areas_3_crista_terminalis_average));
        end
    end
end
for i = 1:number_of_p_waves
    if I_p_wave_areas_4_crista_terminalis_SD == 0 && I_p_wave_areas_4(i) == I_p_wave_areas_4_crista_terminalis_average 
        I_p_wave_areas_4_z_scores_crista_terminalis(i) = 0;
    else
        if I_p_wave_areas_4_crista_terminalis_SD == 0 && ~(I_p_wave_areas_4(i) == I_p_wave_areas_4_crista_terminalis_average)
            I_p_wave_areas_4_z_scores_crista_terminalis(i) = 10*(abs(I_p_wave_areas_4(i)-I_p_wave_areas_4_crista_terminalis_average));
        end
    end
end
for i = 1:number_of_p_waves
    if I_p_wave_areas_5_crista_terminalis_SD == 0 && I_p_wave_areas_5(i) == I_p_wave_areas_5_crista_terminalis_average 
        I_p_wave_areas_5_z_scores_crista_terminalis(i) = 0;
    else
        if I_p_wave_areas_5_crista_terminalis_SD == 0 && ~(I_p_wave_areas_5(i) == I_p_wave_areas_5_crista_terminalis_average)
            I_p_wave_areas_5_z_scores_crista_terminalis(i) = 10*(abs(I_p_wave_areas_5(i)-I_p_wave_areas_5_crista_terminalis_average));
        end
    end
end
for i = 1:number_of_p_waves
    if I_num_zeros_tricuspid_annulus_SD == 0 && I_number_of_zeros_array(i) == I_num_zeros_tricuspid_annulus_average 
        I_num_zeros_z_scores_tricuspid_annulus(i) = 0;
    else
        if I_num_zeros_tricuspid_annulus_SD == 0 && ~(I_number_of_zeros_array(i) == I_num_zeros_tricuspid_annulus_average)
            I_num_zeros_z_scores_tricuspid_annulus(i) = 10*(abs(I_num_zeros_z_scores_tricuspid_annulus(i)-I_num_zeros_tricuspid_annulus_average));
        end
    end
end
for i = 1:number_of_p_waves
    if I_maxima_tricuspid_annulus_SD == 0 && I_maxima(i) == I_maxima_tricuspid_annulus_average 
        I_maxima_z_scores_tricuspid_annulus(i) = 0;
    else
        if I_maxima_tricuspid_annulus_SD == 0 && ~(I_maxima(i) == I_maxima_tricuspid_annulus_average)
            I_maxima_z_scores_tricuspid_annulus(i) = 10*(abs(I_maxima(i)-I_maxima_tricuspid_annulus_average));
        end
    end
end
for i = 1:number_of_p_waves
    if I_minima_tricuspid_annulus_SD == 0 && I_minima(i) == I_minima_tricuspid_annulus_average 
        I_minima_z_scores_tricuspid_annulus(i) = 0;
    else
        if I_minima_tricuspid_annulus_SD == 0 && ~(I_minima(i) == I_minima_tricuspid_annulus_average)
            I_minima_z_scores_tricuspid_annulus(i) = 10*(abs(I_minima(i)-I_minima_tricuspid_annulus_average));
        end
    end
end
for i = 1:number_of_p_waves
    if I_sym_inds_tricuspid_annulus_SD == 0 && I_sym_inds(i) == I_sym_inds_tricuspid_annulus_average 
        I_sym_inds_z_scores_tricuspid_annulus(i) = 0;
    else
        if I_sym_inds_tricuspid_annulus_SD == 0 && ~(I_sym_inds(i) == I_sym_inds_tricuspid_annulus_average)
            I_sym_inds_z_scores_tricuspid_annulus(i) = 10*(abs(I_sym_inds(i)-I_sym_inds_tricuspid_annulus_average));
        end
    end
end
for i = 1:number_of_p_waves
    if I_concavs_tricuspid_annulus_SD == 0 && I_concavs(i) == I_concavs_tricuspid_annulus_average 
        I_concavs_z_scores_tricuspid_annulus(i) = 0;
    else
        if I_concavs_tricuspid_annulus_SD == 0 && ~(I_concavs(i) == I_concavs_tricuspid_annulus_average)
            I_concavs_z_scores_tricuspid_annulus(i) = 10*(abs(I_concavs(i)-I_concavs_tricuspid_annulus_average));
        end
    end
end
for i = 1:number_of_p_waves
    if I_mdslpes_tricuspid_annulus_SD == 0 && I_mdslpes(i) == I_mdslpes_tricuspid_annulus_average 
        I_mdslpes_z_scores_tricuspid_annulus(i) = 0;
    else
        if I_mdslpes_tricuspid_annulus_SD == 0 && ~(I_mdslpes(i) == I_mdslpes_tricuspid_annulus_average)
            I_mdslpes_z_scores_tricuspid_annulus(i) = 10*(abs(I_mdslpes(i)-I_mdslpes_tricuspid_annulus_average));
        end
    end
end
for i = 1:number_of_p_waves
    if I_p_wave_areas_1_tricuspid_annulus_SD == 0 && I_p_wave_areas_1(i) == I_p_wave_areas_1_tricuspid_annulus_average 
        I_p_wave_areas_1_z_scores_tricuspid_annulus(i) = 0;
    else
        if I_p_wave_areas_1_tricuspid_annulus_SD == 0 && ~(I_p_wave_areas_1(i) == I_p_wave_areas_1_tricuspid_annulus_average)
            I_p_wave_areas_1_z_scores_tricuspid_annulus(i) = 10*(abs(I_p_wave_areas_1(i)-I_p_wave_areas_1_tricuspid_annulus_average));
        end
    end
end
for i = 1:number_of_p_waves
    if I_p_wave_areas_2_tricuspid_annulus_SD == 0 && I_p_wave_areas_2(i) == I_p_wave_areas_2_tricuspid_annulus_average 
        I_p_wave_areas_2_z_scores_tricuspid_annulus(i) = 0;
    else
        if I_p_wave_areas_2_tricuspid_annulus_SD == 0 && ~(I_p_wave_areas_2(i) == I_p_wave_areas_2_tricuspid_annulus_average)
            I_p_wave_areas_2_z_scores_tricuspid_annulus(i) = 10*(abs(I_p_wave_areas_2(i)-I_p_wave_areas_2_tricuspid_annulus_average));
        end
    end
end
for i = 1:number_of_p_waves
    if I_p_wave_areas_3_tricuspid_annulus_SD == 0 && I_p_wave_areas_3(i) == I_p_wave_areas_3_tricuspid_annulus_average 
        I_p_wave_areas_3_z_scores_tricuspid_annulus(i) = 0;
    else
        if I_p_wave_areas_3_tricuspid_annulus_SD == 0 && ~(I_p_wave_areas_3(i) == I_p_wave_areas_3_tricuspid_annulus_average)
            I_p_wave_areas_3_z_scores_tricuspid_annulus(i) = 10*(abs(I_p_wave_areas_3(i)-I_p_wave_areas_3_tricuspid_annulus_average));
        end
    end
end
for i = 1:number_of_p_waves
    if I_p_wave_areas_4_tricuspid_annulus_SD == 0 && I_p_wave_areas_4(i) == I_p_wave_areas_4_tricuspid_annulus_average 
        I_p_wave_areas_4_z_scores_tricuspid_annulus(i) = 0;
    else
        if I_p_wave_areas_4_tricuspid_annulus_SD == 0 && ~(I_p_wave_areas_4(i) == I_p_wave_areas_4_tricuspid_annulus_average)
            I_p_wave_areas_4_z_scores_tricuspid_annulus(i) = 10*(abs(I_p_wave_areas_4(i)-I_p_wave_areas_4_tricuspid_annulus_average));
        end
    end
end
for i = 1:number_of_p_waves
    if I_p_wave_areas_5_tricuspid_annulus_SD == 0 && I_p_wave_areas_5(i) == I_p_wave_areas_5_tricuspid_annulus_average 
        I_p_wave_areas_5_z_scores_tricuspid_annulus(i) = 0;
    else
        if I_p_wave_areas_5_tricuspid_annulus_SD == 0 && ~(I_p_wave_areas_5(i) == I_p_wave_areas_5_tricuspid_annulus_average)
            I_p_wave_areas_5_z_scores_tricuspid_annulus(i) = 10*(abs(I_p_wave_areas_5(i)-I_p_wave_areas_5_tricuspid_annulus_average));
        end
    end
end
for i = 1:number_of_p_waves
    if I_num_zeros_coronary_sinus_SD == 0 && I_number_of_zeros_array(i) == I_num_zeros_coronary_sinus_average 
        I_num_zeros_z_scores_coronary_sinus(i) = 0;
    else
        if I_num_zeros_coronary_sinus_SD == 0 && ~(I_number_of_zeros_array(i) == I_num_zeros_coronary_sinus_average)
            I_num_zeros_z_scores_coronary_sinus(i) = 10*(abs(I_num_zeros_z_scores_coronary_sinus(i)-I_num_zeros_coronary_sinus_average));
        end
    end
end
for i = 1:number_of_p_waves
    if I_maxima_coronary_sinus_SD == 0 && I_maxima(i) == I_maxima_coronary_sinus_average 
        I_maxima_z_scores_coronary_sinus(i) = 0;
    else
        if I_maxima_coronary_sinus_SD == 0 && ~(I_maxima(i) == I_maxima_coronary_sinus_average)
            I_maxima_z_scores_coronary_sinus(i) = 10*(abs(I_maxima(i)-I_maxima_coronary_sinus_average));
        end
    end
end
for i = 1:number_of_p_waves
    if I_minima_coronary_sinus_SD == 0 && I_minima(i) == I_minima_coronary_sinus_average 
        I_minima_z_scores_coronary_sinus(i) = 0;
    else
        if I_minima_coronary_sinus_SD == 0 && ~(I_minima(i) == I_minima_coronary_sinus_average)
            I_minima_z_scores_coronary_sinus(i) = 10*(abs(I_minima(i)-I_minima_coronary_sinus_average));
        end
    end
end
for i = 1:number_of_p_waves
    if I_sym_inds_coronary_sinus_SD == 0 && I_sym_inds(i) == I_sym_inds_coronary_sinus_average 
        I_sym_inds_z_scores_coronary_sinus(i) = 0;
    else
        if I_sym_inds_coronary_sinus_SD == 0 && ~(I_sym_inds(i) == I_sym_inds_coronary_sinus_average)
            I_sym_inds_z_scores_coronary_sinus(i) = 10*(abs(I_sym_inds(i)-I_sym_inds_coronary_sinus_average));
        end
    end
end
for i = 1:number_of_p_waves
    if I_concavs_coronary_sinus_SD == 0 && I_concavs(i) == I_concavs_coronary_sinus_average 
        I_concavs_z_scores_coronary_sinus(i) = 0;
    else
        if I_concavs_coronary_sinus_SD == 0 && ~(I_concavs(i) == I_concavs_coronary_sinus_average)
            I_concavs_z_scores_coronary_sinus(i) = 10*(abs(I_concavs(i)-I_concavs_coronary_sinus_average));
        end
    end
end
for i = 1:number_of_p_waves
    if I_mdslpes_coronary_sinus_SD == 0 && I_mdslpes(i) == I_mdslpes_coronary_sinus_average 
        I_mdslpes_z_scores_coronary_sinus(i) = 0;
    else
        if I_mdslpes_coronary_sinus_SD == 0 && ~(I_mdslpes(i) == I_mdslpes_coronary_sinus_average)
            I_mdslpes_z_scores_coronary_sinus(i) = 10*(abs(I_mdslpes(i)-I_mdslpes_coronary_sinus_average));
        end
    end
end
for i = 1:number_of_p_waves
    if I_p_wave_areas_1_coronary_sinus_SD == 0 && I_p_wave_areas_1(i) == I_p_wave_areas_1_coronary_sinus_average 
        I_p_wave_areas_1_z_scores_coronary_sinus(i) = 0;
    else
        if I_p_wave_areas_1_coronary_sinus_SD == 0 && ~(I_p_wave_areas_1(i) == I_p_wave_areas_1_coronary_sinus_average)
            I_p_wave_areas_1_z_scores_coronary_sinus(i) = 10*(abs(I_p_wave_areas_1(i)-I_p_wave_areas_1_coronary_sinus_average));
        end
    end
end
for i = 1:number_of_p_waves
    if I_p_wave_areas_2_coronary_sinus_SD == 0 && I_p_wave_areas_2(i) == I_p_wave_areas_2_coronary_sinus_average 
        I_p_wave_areas_2_z_scores_coronary_sinus(i) = 0;
    else
        if I_p_wave_areas_2_coronary_sinus_SD == 0 && ~(I_p_wave_areas_2(i) == I_p_wave_areas_2_coronary_sinus_average)
            I_p_wave_areas_2_z_scores_coronary_sinus(i) = 10*(abs(I_p_wave_areas_2(i)-I_p_wave_areas_2_coronary_sinus_average));
        end
    end
end
for i = 1:number_of_p_waves
    if I_p_wave_areas_3_coronary_sinus_SD == 0 && I_p_wave_areas_3(i) == I_p_wave_areas_3_coronary_sinus_average 
        I_p_wave_areas_3_z_scores_coronary_sinus(i) = 0;
    else
        if I_p_wave_areas_3_coronary_sinus_SD == 0 && ~(I_p_wave_areas_3(i) == I_p_wave_areas_3_coronary_sinus_average)
            I_p_wave_areas_3_z_scores_coronary_sinus(i) = 10*(abs(I_p_wave_areas_3(i)-I_p_wave_areas_3_coronary_sinus_average));
        end
    end
end
for i = 1:number_of_p_waves
    if I_p_wave_areas_4_coronary_sinus_SD == 0 && I_p_wave_areas_4(i) == I_p_wave_areas_4_coronary_sinus_average 
        I_p_wave_areas_4_z_scores_coronary_sinus(i) = 0;
    else
        if I_p_wave_areas_4_coronary_sinus_SD == 0 && ~(I_p_wave_areas_4(i) == I_p_wave_areas_4_coronary_sinus_average)
            I_p_wave_areas_4_z_scores_coronary_sinus(i) = 10*(abs(I_p_wave_areas_4(i)-I_p_wave_areas_4_coronary_sinus_average));
        end
    end
end
for i = 1:number_of_p_waves
    if I_p_wave_areas_5_coronary_sinus_SD == 0 && I_p_wave_areas_5(i) == I_p_wave_areas_5_coronary_sinus_average 
        I_p_wave_areas_5_z_scores_coronary_sinus(i) = 0;
    else
        if I_p_wave_areas_5_coronary_sinus_SD == 0 && ~(I_p_wave_areas_5(i) == I_p_wave_areas_5_coronary_sinus_average)
            I_p_wave_areas_5_z_scores_coronary_sinus(i) = 10*(abs(I_p_wave_areas_5(i)-I_p_wave_areas_5_coronary_sinus_average));
        end
    end
end
for i = 1:number_of_p_waves
    if I_num_zeros_ostium_SD == 0 && I_number_of_zeros_array(i) == I_num_zeros_ostium_average 
        I_num_zeros_z_scores_ostium(i) = 0;
    else
        if I_num_zeros_ostium_SD == 0 && ~(I_number_of_zeros_array(i) == I_num_zeros_ostium_average)
            I_num_zeros_z_scores_ostium(i) = 10*(abs(I_num_zeros_z_scores_ostium(i)-I_num_zeros_ostium_average));
        end
    end
end
for i = 1:number_of_p_waves
    if I_maxima_ostium_SD == 0 && I_maxima(i) == I_maxima_ostium_average 
        I_maxima_z_scores_ostium(i) = 0;
    else
        if I_maxima_ostium_SD == 0 && ~(I_maxima(i) == I_maxima_ostium_average)
            I_maxima_z_scores_ostium(i) = 10*(abs(I_maxima(i)-I_maxima_ostium_average));
        end
    end
end
for i = 1:number_of_p_waves
    if I_minima_ostium_SD == 0 && I_minima(i) == I_minima_ostium_average 
        I_minima_z_scores_ostium(i) = 0;
    else
        if I_minima_ostium_SD == 0 && ~(I_minima(i) == I_minima_ostium_average)
            I_minima_z_scores_ostium(i) = 10*(abs(I_minima(i)-I_minima_ostium_average));
        end
    end
end
for i = 1:number_of_p_waves
    if I_sym_inds_ostium_SD == 0 && I_sym_inds(i) == I_sym_inds_ostium_average 
        I_sym_inds_z_scores_ostium(i) = 0;
    else
        if I_sym_inds_ostium_SD == 0 && ~(I_sym_inds(i) == I_sym_inds_ostium_average)
            I_sym_inds_z_scores_ostium(i) = 10*(abs(I_sym_inds(i)-I_sym_inds_ostium_average));
        end
    end
end
for i = 1:number_of_p_waves
    if I_concavs_ostium_SD == 0 && I_concavs(i) == I_concavs_ostium_average 
        I_concavs_z_scores_ostium(i) = 0;
    else
        if I_concavs_ostium_SD == 0 && ~(I_concavs(i) == I_concavs_ostium_average)
            I_concavs_z_scores_ostium(i) = 10*(abs(I_concavs(i)-I_concavs_ostium_average));
        end
    end
end
for i = 1:number_of_p_waves
    if I_mdslpes_ostium_SD == 0 && I_mdslpes(i) == I_mdslpes_ostium_average 
        I_mdslpes_z_scores_ostium(i) = 0;
    else
        if I_mdslpes_ostium_SD == 0 && ~(I_mdslpes(i) == I_mdslpes_ostium_average)
            I_mdslpes_z_scores_ostium(i) = 10*(abs(I_mdslpes(i)-I_mdslpes_ostium_average));
        end
    end
end
for i = 1:number_of_p_waves
    if I_p_wave_areas_1_ostium_SD == 0 && I_p_wave_areas_1(i) == I_p_wave_areas_1_ostium_average 
        I_p_wave_areas_1_z_scores_ostium(i) = 0;
    else
        if I_p_wave_areas_1_ostium_SD == 0 && ~(I_p_wave_areas_1(i) == I_p_wave_areas_1_ostium_average)
            I_p_wave_areas_1_z_scores_ostium(i) = 10*(abs(I_p_wave_areas_1(i)-I_p_wave_areas_1_ostium_average));
        end
    end
end
for i = 1:number_of_p_waves
    if I_p_wave_areas_2_ostium_SD == 0 && I_p_wave_areas_2(i) == I_p_wave_areas_2_ostium_average 
        I_p_wave_areas_2_z_scores_ostium(i) = 0;
    else
        if I_p_wave_areas_2_ostium_SD == 0 && ~(I_p_wave_areas_2(i) == I_p_wave_areas_2_ostium_average)
            I_p_wave_areas_2_z_scores_ostium(i) = 10*(abs(I_p_wave_areas_2(i)-I_p_wave_areas_2_ostium_average));
        end
    end
end
for i = 1:number_of_p_waves
    if I_p_wave_areas_3_ostium_SD == 0 && I_p_wave_areas_3(i) == I_p_wave_areas_3_ostium_average 
        I_p_wave_areas_3_z_scores_ostium(i) = 0;
    else
        if I_p_wave_areas_3_ostium_SD == 0 && ~(I_p_wave_areas_3(i) == I_p_wave_areas_3_ostium_average)
            I_p_wave_areas_3_z_scores_ostium(i) = 10*(abs(I_p_wave_areas_3(i)-I_p_wave_areas_3_ostium_average));
        end
    end
end
for i = 1:number_of_p_waves
    if I_p_wave_areas_4_ostium_SD == 0 && I_p_wave_areas_4(i) == I_p_wave_areas_4_ostium_average 
        I_p_wave_areas_4_z_scores_ostium(i) = 0;
    else
        if I_p_wave_areas_4_ostium_SD == 0 && ~(I_p_wave_areas_4(i) == I_p_wave_areas_4_ostium_average)
            I_p_wave_areas_4_z_scores_ostium(i) = 10*(abs(I_p_wave_areas_4(i)-I_p_wave_areas_4_ostium_average));
        end
    end
end
for i = 1:number_of_p_waves
    if I_p_wave_areas_5_ostium_SD == 0 && I_p_wave_areas_5(i) == I_p_wave_areas_5_ostium_average 
        I_p_wave_areas_5_z_scores_ostium(i) = 0;
    else
        if I_p_wave_areas_5_ostium_SD == 0 && ~(I_p_wave_areas_5(i) == I_p_wave_areas_5_ostium_average)
            I_p_wave_areas_5_z_scores_ostium(i) = 10*(abs(I_p_wave_areas_5(i)-I_p_wave_areas_5_ostium_average));
        end
    end
end
for i = 1:number_of_p_waves
    if I_num_zeros_perinodal_SD == 0 && I_number_of_zeros_array(i) == I_num_zeros_perinodal_average 
        I_num_zeros_z_scores_perinodal(i) = 0;
    else
        if I_num_zeros_perinodal_SD == 0 && ~(I_number_of_zeros_array(i) == I_num_zeros_perinodal_average)
            I_num_zeros_z_scores_perinodal(i) = 10*(abs(I_num_zeros_z_scores_perinodal(i)-I_num_zeros_perinodal_average));
        end
    end
end
for i = 1:number_of_p_waves
    if I_maxima_perinodal_SD == 0 && I_maxima(i) == I_maxima_perinodal_average 
        I_maxima_z_scores_perinodal(i) = 0;
    else
        if I_maxima_perinodal_SD == 0 && ~(I_maxima(i) == I_maxima_perinodal_average)
            I_maxima_z_scores_perinodal(i) = 10*(abs(I_maxima(i)-I_maxima_perinodal_average));
        end
    end
end
for i = 1:number_of_p_waves
    if I_minima_perinodal_SD == 0 && I_minima(i) == I_minima_perinodal_average 
        I_minima_z_scores_perinodal(i) = 0;
    else
        if I_minima_perinodal_SD == 0 && ~(I_minima(i) == I_minima_perinodal_average)
            I_minima_z_scores_perinodal(i) = 10*(abs(I_minima(i)-I_minima_perinodal_average));
        end
    end
end
for i = 1:number_of_p_waves
    if I_sym_inds_perinodal_SD == 0 && I_sym_inds(i) == I_sym_inds_perinodal_average 
        I_sym_inds_z_scores_perinodal(i) = 0;
    else
        if I_sym_inds_perinodal_SD == 0 && ~(I_sym_inds(i) == I_sym_inds_perinodal_average)
            I_sym_inds_z_scores_perinodal(i) = 10*(abs(I_sym_inds(i)-I_sym_inds_perinodal_average));
        end
    end
end
for i = 1:number_of_p_waves
    if I_concavs_perinodal_SD == 0 && I_concavs(i) == I_concavs_perinodal_average 
        I_concavs_z_scores_perinodal(i) = 0;
    else
        if I_concavs_perinodal_SD == 0 && ~(I_concavs(i) == I_concavs_perinodal_average)
            I_concavs_z_scores_perinodal(i) = 10*(abs(I_concavs(i)-I_concavs_perinodal_average));
        end
    end
end
for i = 1:number_of_p_waves
    if I_mdslpes_perinodal_SD == 0 && I_mdslpes(i) == I_mdslpes_perinodal_average 
        I_mdslpes_z_scores_perinodal(i) = 0;
    else
        if I_mdslpes_perinodal_SD == 0 && ~(I_mdslpes(i) == I_mdslpes_perinodal_average)
            I_mdslpes_z_scores_perinodal(i) = 10*(abs(I_mdslpes(i)-I_mdslpes_perinodal_average));
        end
    end
end
for i = 1:number_of_p_waves
    if I_p_wave_areas_1_perinodal_SD == 0 && I_p_wave_areas_1(i) == I_p_wave_areas_1_perinodal_average 
        I_p_wave_areas_1_z_scores_perinodal(i) = 0;
    else
        if I_p_wave_areas_1_perinodal_SD == 0 && ~(I_p_wave_areas_1(i) == I_p_wave_areas_1_perinodal_average)
            I_p_wave_areas_1_z_scores_perinodal(i) = 10*(abs(I_p_wave_areas_1(i)-I_p_wave_areas_1_perinodal_average));
        end
    end
end
for i = 1:number_of_p_waves
    if I_p_wave_areas_2_perinodal_SD == 0 && I_p_wave_areas_2(i) == I_p_wave_areas_2_perinodal_average 
        I_p_wave_areas_2_z_scores_perinodal(i) = 0;
    else
        if I_p_wave_areas_2_perinodal_SD == 0 && ~(I_p_wave_areas_2(i) == I_p_wave_areas_2_perinodal_average)
            I_p_wave_areas_2_z_scores_perinodal(i) = 10*(abs(I_p_wave_areas_2(i)-I_p_wave_areas_2_perinodal_average));
        end
    end
end
for i = 1:number_of_p_waves
    if I_p_wave_areas_3_perinodal_SD == 0 && I_p_wave_areas_3(i) == I_p_wave_areas_3_perinodal_average 
        I_p_wave_areas_3_z_scores_perinodal(i) = 0;
    else
        if I_p_wave_areas_3_perinodal_SD == 0 && ~(I_p_wave_areas_3(i) == I_p_wave_areas_3_perinodal_average)
            I_p_wave_areas_3_z_scores_perinodal(i) = 10*(abs(I_p_wave_areas_3(i)-I_p_wave_areas_3_perinodal_average));
        end
    end
end
for i = 1:number_of_p_waves
    if I_p_wave_areas_4_perinodal_SD == 0 && I_p_wave_areas_4(i) == I_p_wave_areas_4_perinodal_average 
        I_p_wave_areas_4_z_scores_perinodal(i) = 0;
    else
        if I_p_wave_areas_4_perinodal_SD == 0 && ~(I_p_wave_areas_4(i) == I_p_wave_areas_4_perinodal_average)
            I_p_wave_areas_4_z_scores_perinodal(i) = 10*(abs(I_p_wave_areas_4(i)-I_p_wave_areas_4_perinodal_average));
        end
    end
end
for i = 1:number_of_p_waves
    if I_p_wave_areas_5_perinodal_SD == 0 && I_p_wave_areas_5(i) == I_p_wave_areas_5_perinodal_average 
        I_p_wave_areas_5_z_scores_perinodal(i) = 0;
    else
        if I_p_wave_areas_5_perinodal_SD == 0 && ~(I_p_wave_areas_5(i) == I_p_wave_areas_5_perinodal_average)
            I_p_wave_areas_5_z_scores_perinodal(i) = 10*(abs(I_p_wave_areas_5(i)-I_p_wave_areas_5_perinodal_average));
        end
    end
end
for i = 1:number_of_p_waves
    if I_num_zeros_right_septum_SD == 0 && I_number_of_zeros_array(i) == I_num_zeros_right_septum_average 
        I_num_zeros_z_scores_right_septum(i) = 0;
    else
        if I_num_zeros_right_septum_SD == 0 && ~(I_number_of_zeros_array(i) == I_num_zeros_right_septum_average)
            I_num_zeros_z_scores_right_septum(i) = 10*(abs(I_num_zeros_z_scores_right_septum(i)-I_num_zeros_right_septum_average));
        end
    end
end
for i = 1:number_of_p_waves
    if I_maxima_right_septum_SD == 0 && I_maxima(i) == I_maxima_right_septum_average 
        I_maxima_z_scores_right_septum(i) = 0;
    else
        if I_maxima_right_septum_SD == 0 && ~(I_maxima(i) == I_maxima_right_septum_average)
            I_maxima_z_scores_right_septum(i) = 10*(abs(I_maxima(i)-I_maxima_right_septum_average));
        end
    end
end
for i = 1:number_of_p_waves
    if I_minima_right_septum_SD == 0 && I_minima(i) == I_minima_right_septum_average 
        I_minima_z_scores_right_septum(i) = 0;
    else
        if I_minima_right_septum_SD == 0 && ~(I_minima(i) == I_minima_right_septum_average)
            I_minima_z_scores_right_septum(i) = 10*(abs(I_minima(i)-I_minima_right_septum_average));
        end
    end
end
for i = 1:number_of_p_waves
    if I_sym_inds_right_septum_SD == 0 && I_sym_inds(i) == I_sym_inds_right_septum_average 
        I_sym_inds_z_scores_right_septum(i) = 0;
    else
        if I_sym_inds_right_septum_SD == 0 && ~(I_sym_inds(i) == I_sym_inds_right_septum_average)
            I_sym_inds_z_scores_right_septum(i) = 10*(abs(I_sym_inds(i)-I_sym_inds_right_septum_average));
        end
    end
end
for i = 1:number_of_p_waves
    if I_concavs_right_septum_SD == 0 && I_concavs(i) == I_concavs_right_septum_average 
        I_concavs_z_scores_right_septum(i) = 0;
    else
        if I_concavs_right_septum_SD == 0 && ~(I_concavs(i) == I_concavs_right_septum_average)
            I_concavs_z_scores_right_septum(i) = 10*(abs(I_concavs(i)-I_concavs_right_septum_average));
        end
    end
end
for i = 1:number_of_p_waves
    if I_mdslpes_right_septum_SD == 0 && I_mdslpes(i) == I_mdslpes_right_septum_average 
        I_mdslpes_z_scores_right_septum(i) = 0;
    else
        if I_mdslpes_right_septum_SD == 0 && ~(I_mdslpes(i) == I_mdslpes_right_septum_average)
            I_mdslpes_z_scores_right_septum(i) = 10*(abs(I_mdslpes(i)-I_mdslpes_right_septum_average));
        end
    end
end
for i = 1:number_of_p_waves
    if I_p_wave_areas_1_right_septum_SD == 0 && I_p_wave_areas_1(i) == I_p_wave_areas_1_right_septum_average 
        I_p_wave_areas_1_z_scores_right_septum(i) = 0;
    else
        if I_p_wave_areas_1_right_septum_SD == 0 && ~(I_p_wave_areas_1(i) == I_p_wave_areas_1_right_septum_average)
            I_p_wave_areas_1_z_scores_right_septum(i) = 10*(abs(I_p_wave_areas_1(i)-I_p_wave_areas_1_right_septum_average));
        end
    end
end
for i = 1:number_of_p_waves
    if I_p_wave_areas_2_right_septum_SD == 0 && I_p_wave_areas_2(i) == I_p_wave_areas_2_right_septum_average 
        I_p_wave_areas_2_z_scores_right_septum(i) = 0;
    else
        if I_p_wave_areas_2_right_septum_SD == 0 && ~(I_p_wave_areas_2(i) == I_p_wave_areas_2_right_septum_average)
            I_p_wave_areas_2_z_scores_right_septum(i) = 10*(abs(I_p_wave_areas_2(i)-I_p_wave_areas_2_right_septum_average));
        end
    end
end
for i = 1:number_of_p_waves
    if I_p_wave_areas_3_right_septum_SD == 0 && I_p_wave_areas_3(i) == I_p_wave_areas_3_right_septum_average 
        I_p_wave_areas_3_z_scores_right_septum(i) = 0;
    else
        if I_p_wave_areas_3_right_septum_SD == 0 && ~(I_p_wave_areas_3(i) == I_p_wave_areas_3_right_septum_average)
            I_p_wave_areas_3_z_scores_right_septum(i) = 10*(abs(I_p_wave_areas_3(i)-I_p_wave_areas_3_right_septum_average));
        end
    end
end
for i = 1:number_of_p_waves
    if I_p_wave_areas_4_right_septum_SD == 0 && I_p_wave_areas_4(i) == I_p_wave_areas_4_right_septum_average 
        I_p_wave_areas_4_z_scores_right_septum(i) = 0;
    else
        if I_p_wave_areas_4_right_septum_SD == 0 && ~(I_p_wave_areas_4(i) == I_p_wave_areas_4_right_septum_average)
            I_p_wave_areas_4_z_scores_right_septum(i) = 10*(abs(I_p_wave_areas_4(i)-I_p_wave_areas_4_right_septum_average));
        end
    end
end
for i = 1:number_of_p_waves
    if I_p_wave_areas_5_right_septum_SD == 0 && I_p_wave_areas_5(i) == I_p_wave_areas_5_right_septum_average 
        I_p_wave_areas_5_z_scores_right_septum(i) = 0;
    else
        if I_p_wave_areas_5_right_septum_SD == 0 && ~(I_p_wave_areas_5(i) == I_p_wave_areas_5_right_septum_average)
            I_p_wave_areas_5_z_scores_right_septum(i) = 10*(abs(I_p_wave_areas_5(i)-I_p_wave_areas_5_right_septum_average));
        end
    end
end
for i = 1:number_of_p_waves
    if I_num_zeros_right_atrial_appendage_SD == 0 && I_number_of_zeros_array(i) == I_num_zeros_right_atrial_appendage_average 
        I_num_zeros_z_scores_right_atrial_appendage(i) = 0;
    else
        if I_num_zeros_right_atrial_appendage_SD == 0 && ~(I_number_of_zeros_array(i) == I_num_zeros_right_atrial_appendage_average)
            I_num_zeros_z_scores_right_atrial_appendage(i) = 10*(abs(I_num_zeros_z_scores_right_atrial_appendage(i)-I_num_zeros_right_atrial_appendage_average));
        end
    end
end
for i = 1:number_of_p_waves
    if I_maxima_right_atrial_appendage_SD == 0 && I_maxima(i) == I_maxima_right_atrial_appendage_average 
        I_maxima_z_scores_right_atrial_appendage(i) = 0;
    else
        if I_maxima_right_atrial_appendage_SD == 0 && ~(I_maxima(i) == I_maxima_right_atrial_appendage_average)
            I_maxima_z_scores_right_atrial_appendage(i) = 10*(abs(I_maxima(i)-I_maxima_right_atrial_appendage_average));
        end
    end
end
for i = 1:number_of_p_waves
    if I_minima_right_atrial_appendage_SD == 0 && I_minima(i) == I_minima_right_atrial_appendage_average 
        I_minima_z_scores_right_atrial_appendage(i) = 0;
    else
        if I_minima_right_atrial_appendage_SD == 0 && ~(I_minima(i) == I_minima_right_atrial_appendage_average)
            I_minima_z_scores_right_atrial_appendage(i) = 10*(abs(I_minima(i)-I_minima_right_atrial_appendage_average));
        end
    end
end
for i = 1:number_of_p_waves
    if I_sym_inds_right_atrial_appendage_SD == 0 && I_sym_inds(i) == I_sym_inds_right_atrial_appendage_average 
        I_sym_inds_z_scores_right_atrial_appendage(i) = 0;
    else
        if I_sym_inds_right_atrial_appendage_SD == 0 && ~(I_sym_inds(i) == I_sym_inds_right_atrial_appendage_average)
            I_sym_inds_z_scores_right_atrial_appendage(i) = 10*(abs(I_sym_inds(i)-I_sym_inds_right_atrial_appendage_average));
        end
    end
end
for i = 1:number_of_p_waves
    if I_concavs_right_atrial_appendage_SD == 0 && I_concavs(i) == I_concavs_right_atrial_appendage_average 
        I_concavs_z_scores_right_atrial_appendage(i) = 0;
    else
        if I_concavs_right_atrial_appendage_SD == 0 && ~(I_concavs(i) == I_concavs_right_atrial_appendage_average)
            I_concavs_z_scores_right_atrial_appendage(i) = 10*(abs(I_concavs(i)-I_concavs_right_atrial_appendage_average));
        end
    end
end
for i = 1:number_of_p_waves
    if I_mdslpes_right_atrial_appendage_SD == 0 && I_mdslpes(i) == I_mdslpes_right_atrial_appendage_average 
        I_mdslpes_z_scores_right_atrial_appendage(i) = 0;
    else
        if I_mdslpes_right_atrial_appendage_SD == 0 && ~(I_mdslpes(i) == I_mdslpes_right_atrial_appendage_average)
            I_mdslpes_z_scores_right_atrial_appendage(i) = 10*(abs(I_mdslpes(i)-I_mdslpes_right_atrial_appendage_average));
        end
    end
end
for i = 1:number_of_p_waves
    if I_p_wave_areas_1_right_atrial_appendage_SD == 0 && I_p_wave_areas_1(i) == I_p_wave_areas_1_right_atrial_appendage_average 
        I_p_wave_areas_1_z_scores_right_atrial_appendage(i) = 0;
    else
        if I_p_wave_areas_1_right_atrial_appendage_SD == 0 && ~(I_p_wave_areas_1(i) == I_p_wave_areas_1_right_atrial_appendage_average)
            I_p_wave_areas_1_z_scores_right_atrial_appendage(i) = 10*(abs(I_p_wave_areas_1(i)-I_p_wave_areas_1_right_atrial_appendage_average));
        end
    end
end
for i = 1:number_of_p_waves
    if I_p_wave_areas_2_right_atrial_appendage_SD == 0 && I_p_wave_areas_2(i) == I_p_wave_areas_2_right_atrial_appendage_average 
        I_p_wave_areas_2_z_scores_right_atrial_appendage(i) = 0;
    else
        if I_p_wave_areas_2_right_atrial_appendage_SD == 0 && ~(I_p_wave_areas_2(i) == I_p_wave_areas_2_right_atrial_appendage_average)
            I_p_wave_areas_2_z_scores_right_atrial_appendage(i) = 10*(abs(I_p_wave_areas_2(i)-I_p_wave_areas_2_right_atrial_appendage_average));
        end
    end
end
for i = 1:number_of_p_waves
    if I_p_wave_areas_3_right_atrial_appendage_SD == 0 && I_p_wave_areas_3(i) == I_p_wave_areas_3_right_atrial_appendage_average 
        I_p_wave_areas_3_z_scores_right_atrial_appendage(i) = 0;
    else
        if I_p_wave_areas_3_right_atrial_appendage_SD == 0 && ~(I_p_wave_areas_3(i) == I_p_wave_areas_3_right_atrial_appendage_average)
            I_p_wave_areas_3_z_scores_right_atrial_appendage(i) = 10*(abs(I_p_wave_areas_3(i)-I_p_wave_areas_3_right_atrial_appendage_average));
        end
    end
end
for i = 1:number_of_p_waves
    if I_p_wave_areas_4_right_atrial_appendage_SD == 0 && I_p_wave_areas_4(i) == I_p_wave_areas_4_right_atrial_appendage_average 
        I_p_wave_areas_4_z_scores_right_atrial_appendage(i) = 0;
    else
        if I_p_wave_areas_4_right_atrial_appendage_SD == 0 && ~(I_p_wave_areas_4(i) == I_p_wave_areas_4_right_atrial_appendage_average)
            I_p_wave_areas_4_z_scores_right_atrial_appendage(i) = 10*(abs(I_p_wave_areas_4(i)-I_p_wave_areas_4_right_atrial_appendage_average));
        end
    end
end
for i = 1:number_of_p_waves
    if I_p_wave_areas_5_right_atrial_appendage_SD == 0 && I_p_wave_areas_5(i) == I_p_wave_areas_5_right_atrial_appendage_average 
        I_p_wave_areas_5_z_scores_right_atrial_appendage(i) = 0;
    else
        if I_p_wave_areas_5_right_atrial_appendage_SD == 0 && ~(I_p_wave_areas_5(i) == I_p_wave_areas_5_right_atrial_appendage_average)
            I_p_wave_areas_5_z_scores_right_atrial_appendage(i) = 10*(abs(I_p_wave_areas_5(i)-I_p_wave_areas_5_right_atrial_appendage_average));
        end
    end
end
for i = 1:number_of_p_waves
    if I_num_zeros_pulmonary_veins_SD == 0 && I_number_of_zeros_array(i) == I_num_zeros_pulmonary_veins_average 
        I_num_zeros_z_scores_pulmonary_veins(i) = 0;
    else
        if I_num_zeros_pulmonary_veins_SD == 0 && ~(I_number_of_zeros_array(i) == I_num_zeros_pulmonary_veins_average)
            I_num_zeros_z_scores_pulmonary_veins(i) = 10*(abs(I_num_zeros_z_scores_pulmonary_veins(i)-I_num_zeros_pulmonary_veins_average));
        end
    end
end
for i = 1:number_of_p_waves
    if I_maxima_pulmonary_veins_SD == 0 && I_maxima(i) == I_maxima_pulmonary_veins_average 
        I_maxima_z_scores_pulmonary_veins(i) = 0;
    else
        if I_maxima_pulmonary_veins_SD == 0 && ~(I_maxima(i) == I_maxima_pulmonary_veins_average)
            I_maxima_z_scores_pulmonary_veins(i) = 10*(abs(I_maxima(i)-I_maxima_pulmonary_veins_average));
        end
    end
end
for i = 1:number_of_p_waves
    if I_minima_pulmonary_veins_SD == 0 && I_minima(i) == I_minima_pulmonary_veins_average 
        I_minima_z_scores_pulmonary_veins(i) = 0;
    else
        if I_minima_pulmonary_veins_SD == 0 && ~(I_minima(i) == I_minima_pulmonary_veins_average)
            I_minima_z_scores_pulmonary_veins(i) = 10*(abs(I_minima(i)-I_minima_pulmonary_veins_average));
        end
    end
end
for i = 1:number_of_p_waves
    if I_sym_inds_pulmonary_veins_SD == 0 && I_sym_inds(i) == I_sym_inds_pulmonary_veins_average 
        I_sym_inds_z_scores_pulmonary_veins(i) = 0;
    else
        if I_sym_inds_pulmonary_veins_SD == 0 && ~(I_sym_inds(i) == I_sym_inds_pulmonary_veins_average)
            I_sym_inds_z_scores_pulmonary_veins(i) = 10*(abs(I_sym_inds(i)-I_sym_inds_pulmonary_veins_average));
        end
    end
end
for i = 1:number_of_p_waves
    if I_concavs_pulmonary_veins_SD == 0 && I_concavs(i) == I_concavs_pulmonary_veins_average 
        I_concavs_z_scores_pulmonary_veins(i) = 0;
    else
        if I_concavs_pulmonary_veins_SD == 0 && ~(I_concavs(i) == I_concavs_pulmonary_veins_average)
            I_concavs_z_scores_pulmonary_veins(i) = 10*(abs(I_concavs(i)-I_concavs_pulmonary_veins_average));
        end
    end
end
for i = 1:number_of_p_waves
    if I_mdslpes_pulmonary_veins_SD == 0 && I_mdslpes(i) == I_mdslpes_pulmonary_veins_average 
        I_mdslpes_z_scores_pulmonary_veins(i) = 0;
    else
        if I_mdslpes_pulmonary_veins_SD == 0 && ~(I_mdslpes(i) == I_mdslpes_pulmonary_veins_average)
            I_mdslpes_z_scores_pulmonary_veins(i) = 10*(abs(I_mdslpes(i)-I_mdslpes_pulmonary_veins_average));
        end
    end
end
for i = 1:number_of_p_waves
    if I_p_wave_areas_1_pulmonary_veins_SD == 0 && I_p_wave_areas_1(i) == I_p_wave_areas_1_pulmonary_veins_average 
        I_p_wave_areas_1_z_scores_pulmonary_veins(i) = 0;
    else
        if I_p_wave_areas_1_pulmonary_veins_SD == 0 && ~(I_p_wave_areas_1(i) == I_p_wave_areas_1_pulmonary_veins_average)
            I_p_wave_areas_1_z_scores_pulmonary_veins(i) = 10*(abs(I_p_wave_areas_1(i)-I_p_wave_areas_1_pulmonary_veins_average));
        end
    end
end
for i = 1:number_of_p_waves
    if I_p_wave_areas_2_pulmonary_veins_SD == 0 && I_p_wave_areas_2(i) == I_p_wave_areas_2_pulmonary_veins_average 
        I_p_wave_areas_2_z_scores_pulmonary_veins(i) = 0;
    else
        if I_p_wave_areas_2_pulmonary_veins_SD == 0 && ~(I_p_wave_areas_2(i) == I_p_wave_areas_2_pulmonary_veins_average)
            I_p_wave_areas_2_z_scores_pulmonary_veins(i) = 10*(abs(I_p_wave_areas_2(i)-I_p_wave_areas_2_pulmonary_veins_average));
        end
    end
end
for i = 1:number_of_p_waves
    if I_p_wave_areas_3_pulmonary_veins_SD == 0 && I_p_wave_areas_3(i) == I_p_wave_areas_3_pulmonary_veins_average 
        I_p_wave_areas_3_z_scores_pulmonary_veins(i) = 0;
    else
        if I_p_wave_areas_3_pulmonary_veins_SD == 0 && ~(I_p_wave_areas_3(i) == I_p_wave_areas_3_pulmonary_veins_average)
            I_p_wave_areas_3_z_scores_pulmonary_veins(i) = 10*(abs(I_p_wave_areas_3(i)-I_p_wave_areas_3_pulmonary_veins_average));
        end
    end
end
for i = 1:number_of_p_waves
    if I_p_wave_areas_4_pulmonary_veins_SD == 0 && I_p_wave_areas_4(i) == I_p_wave_areas_4_pulmonary_veins_average 
        I_p_wave_areas_4_z_scores_pulmonary_veins(i) = 0;
    else
        if I_p_wave_areas_4_pulmonary_veins_SD == 0 && ~(I_p_wave_areas_4(i) == I_p_wave_areas_4_pulmonary_veins_average)
            I_p_wave_areas_4_z_scores_pulmonary_veins(i) = 10*(abs(I_p_wave_areas_4(i)-I_p_wave_areas_4_pulmonary_veins_average));
        end
    end
end
for i = 1:number_of_p_waves
    if I_p_wave_areas_5_pulmonary_veins_SD == 0 && I_p_wave_areas_5(i) == I_p_wave_areas_5_pulmonary_veins_average 
        I_p_wave_areas_5_z_scores_pulmonary_veins(i) = 0;
    else
        if I_p_wave_areas_5_pulmonary_veins_SD == 0 && ~(I_p_wave_areas_5(i) == I_p_wave_areas_5_pulmonary_veins_average)
            I_p_wave_areas_5_z_scores_pulmonary_veins(i) = 10*(abs(I_p_wave_areas_5(i)-I_p_wave_areas_5_pulmonary_veins_average));
        end
    end
end
for i = 1:number_of_p_waves
    if I_num_zeros_mitral_annulus_SD == 0 && I_number_of_zeros_array(i) == I_num_zeros_mitral_annulus_average 
        I_num_zeros_z_scores_mitral_annulus(i) = 0;
    else
        if I_num_zeros_mitral_annulus_SD == 0 && ~(I_number_of_zeros_array(i) == I_num_zeros_mitral_annulus_average)
            I_num_zeros_z_scores_mitral_annulus(i) = 10*(abs(I_num_zeros_z_scores_mitral_annulus(i)-I_num_zeros_mitral_annulus_average));
        end
    end
end
for i = 1:number_of_p_waves
    if I_maxima_mitral_annulus_SD == 0 && I_maxima(i) == I_maxima_mitral_annulus_average 
        I_maxima_z_scores_mitral_annulus(i) = 0;
    else
        if I_maxima_mitral_annulus_SD == 0 && ~(I_maxima(i) == I_maxima_mitral_annulus_average)
            I_maxima_z_scores_mitral_annulus(i) = 10*(abs(I_maxima(i)-I_maxima_mitral_annulus_average));
        end
    end
end
for i = 1:number_of_p_waves
    if I_minima_mitral_annulus_SD == 0 && I_minima(i) == I_minima_mitral_annulus_average 
        I_minima_z_scores_mitral_annulus(i) = 0;
    else
        if I_minima_mitral_annulus_SD == 0 && ~(I_minima(i) == I_minima_mitral_annulus_average)
            I_minima_z_scores_mitral_annulus(i) = 10*(abs(I_minima(i)-I_minima_mitral_annulus_average));
        end
    end
end
for i = 1:number_of_p_waves
    if I_sym_inds_mitral_annulus_SD == 0 && I_sym_inds(i) == I_sym_inds_mitral_annulus_average 
        I_sym_inds_z_scores_mitral_annulus(i) = 0;
    else
        if I_sym_inds_mitral_annulus_SD == 0 && ~(I_sym_inds(i) == I_sym_inds_mitral_annulus_average)
            I_sym_inds_z_scores_mitral_annulus(i) = 10*(abs(I_sym_inds(i)-I_sym_inds_mitral_annulus_average));
        end
    end
end
for i = 1:number_of_p_waves
    if I_concavs_mitral_annulus_SD == 0 && I_concavs(i) == I_concavs_mitral_annulus_average 
        I_concavs_z_scores_mitral_annulus(i) = 0;
    else
        if I_concavs_mitral_annulus_SD == 0 && ~(I_concavs(i) == I_concavs_mitral_annulus_average)
            I_concavs_z_scores_mitral_annulus(i) = 10*(abs(I_concavs(i)-I_concavs_mitral_annulus_average));
        end
    end
end
for i = 1:number_of_p_waves
    if I_mdslpes_mitral_annulus_SD == 0 && I_mdslpes(i) == I_mdslpes_mitral_annulus_average 
        I_mdslpes_z_scores_mitral_annulus(i) = 0;
    else
        if I_mdslpes_mitral_annulus_SD == 0 && ~(I_mdslpes(i) == I_mdslpes_mitral_annulus_average)
            I_mdslpes_z_scores_mitral_annulus(i) = 10*(abs(I_mdslpes(i)-I_mdslpes_mitral_annulus_average));
        end
    end
end
for i = 1:number_of_p_waves
    if I_p_wave_areas_1_mitral_annulus_SD == 0 && I_p_wave_areas_1(i) == I_p_wave_areas_1_mitral_annulus_average 
        I_p_wave_areas_1_z_scores_mitral_annulus(i) = 0;
    else
        if I_p_wave_areas_1_mitral_annulus_SD == 0 && ~(I_p_wave_areas_1(i) == I_p_wave_areas_1_mitral_annulus_average)
            I_p_wave_areas_1_z_scores_mitral_annulus(i) = 10*(abs(I_p_wave_areas_1(i)-I_p_wave_areas_1_mitral_annulus_average));
        end
    end
end
for i = 1:number_of_p_waves
    if I_p_wave_areas_2_mitral_annulus_SD == 0 && I_p_wave_areas_2(i) == I_p_wave_areas_2_mitral_annulus_average 
        I_p_wave_areas_2_z_scores_mitral_annulus(i) = 0;
    else
        if I_p_wave_areas_2_mitral_annulus_SD == 0 && ~(I_p_wave_areas_2(i) == I_p_wave_areas_2_mitral_annulus_average)
            I_p_wave_areas_2_z_scores_mitral_annulus(i) = 10*(abs(I_p_wave_areas_2(i)-I_p_wave_areas_2_mitral_annulus_average));
        end
    end
end
for i = 1:number_of_p_waves
    if I_p_wave_areas_3_mitral_annulus_SD == 0 && I_p_wave_areas_3(i) == I_p_wave_areas_3_mitral_annulus_average 
        I_p_wave_areas_3_z_scores_mitral_annulus(i) = 0;
    else
        if I_p_wave_areas_3_mitral_annulus_SD == 0 && ~(I_p_wave_areas_3(i) == I_p_wave_areas_3_mitral_annulus_average)
            I_p_wave_areas_3_z_scores_mitral_annulus(i) = 10*(abs(I_p_wave_areas_3(i)-I_p_wave_areas_3_mitral_annulus_average));
        end
    end
end
for i = 1:number_of_p_waves
    if I_p_wave_areas_4_mitral_annulus_SD == 0 && I_p_wave_areas_4(i) == I_p_wave_areas_4_mitral_annulus_average 
        I_p_wave_areas_4_z_scores_mitral_annulus(i) = 0;
    else
        if I_p_wave_areas_4_mitral_annulus_SD == 0 && ~(I_p_wave_areas_4(i) == I_p_wave_areas_4_mitral_annulus_average)
            I_p_wave_areas_4_z_scores_mitral_annulus(i) = 10*(abs(I_p_wave_areas_4(i)-I_p_wave_areas_4_mitral_annulus_average));
        end
    end
end
for i = 1:number_of_p_waves
    if I_p_wave_areas_5_mitral_annulus_SD == 0 && I_p_wave_areas_5(i) == I_p_wave_areas_5_mitral_annulus_average 
        I_p_wave_areas_5_z_scores_mitral_annulus(i) = 0;
    else
        if I_p_wave_areas_5_mitral_annulus_SD == 0 && ~(I_p_wave_areas_5(i) == I_p_wave_areas_5_mitral_annulus_average)
            I_p_wave_areas_5_z_scores_mitral_annulus(i) = 10*(abs(I_p_wave_areas_5(i)-I_p_wave_areas_5_mitral_annulus_average));
        end
    end
end
for i = 1:number_of_p_waves
    if I_num_zeros_CS_body_SD == 0 && I_number_of_zeros_array(i) == I_num_zeros_CS_body_average 
        I_num_zeros_z_scores_CS_body(i) = 0;
    else
        if I_num_zeros_CS_body_SD == 0 && ~(I_number_of_zeros_array(i) == I_num_zeros_CS_body_average)
            I_num_zeros_z_scores_CS_body(i) = 10*(abs(I_num_zeros_z_scores_CS_body(i)-I_num_zeros_CS_body_average));
        end
    end
end
for i = 1:number_of_p_waves
    if I_maxima_CS_body_SD == 0 && I_maxima(i) == I_maxima_CS_body_average 
        I_maxima_z_scores_CS_body(i) = 0;
    else
        if I_maxima_CS_body_SD == 0 && ~(I_maxima(i) == I_maxima_CS_body_average)
            I_maxima_z_scores_CS_body(i) = 10*(abs(I_maxima(i)-I_maxima_CS_body_average));
        end
    end
end
for i = 1:number_of_p_waves
    if I_minima_CS_body_SD == 0 && I_minima(i) == I_minima_CS_body_average 
        I_minima_z_scores_CS_body(i) = 0;
    else
        if I_minima_CS_body_SD == 0 && ~(I_minima(i) == I_minima_CS_body_average)
            I_minima_z_scores_CS_body(i) = 10*(abs(I_minima(i)-I_minima_CS_body_average));
        end
    end
end
for i = 1:number_of_p_waves
    if I_sym_inds_CS_body_SD == 0 && I_sym_inds(i) == I_sym_inds_CS_body_average 
        I_sym_inds_z_scores_CS_body(i) = 0;
    else
        if I_sym_inds_CS_body_SD == 0 && ~(I_sym_inds(i) == I_sym_inds_CS_body_average)
            I_sym_inds_z_scores_CS_body(i) = 10*(abs(I_sym_inds(i)-I_sym_inds_CS_body_average));
        end
    end
end
for i = 1:number_of_p_waves
    if I_concavs_CS_body_SD == 0 && I_concavs(i) == I_concavs_CS_body_average 
        I_concavs_z_scores_CS_body(i) = 0;
    else
        if I_concavs_CS_body_SD == 0 && ~(I_concavs(i) == I_concavs_CS_body_average)
            I_concavs_z_scores_CS_body(i) = 10*(abs(I_concavs(i)-I_concavs_CS_body_average));
        end
    end
end
for i = 1:number_of_p_waves
    if I_mdslpes_CS_body_SD == 0 && I_mdslpes(i) == I_mdslpes_CS_body_average 
        I_mdslpes_z_scores_CS_body(i) = 0;
    else
        if I_mdslpes_CS_body_SD == 0 && ~(I_mdslpes(i) == I_mdslpes_CS_body_average)
            I_mdslpes_z_scores_CS_body(i) = 10*(abs(I_mdslpes(i)-I_mdslpes_CS_body_average));
        end
    end
end
for i = 1:number_of_p_waves
    if I_p_wave_areas_1_CS_body_SD == 0 && I_p_wave_areas_1(i) == I_p_wave_areas_1_CS_body_average 
        I_p_wave_areas_1_z_scores_CS_body(i) = 0;
    else
        if I_p_wave_areas_1_CS_body_SD == 0 && ~(I_p_wave_areas_1(i) == I_p_wave_areas_1_CS_body_average)
            I_p_wave_areas_1_z_scores_CS_body(i) = 10*(abs(I_p_wave_areas_1(i)-I_p_wave_areas_1_CS_body_average));
        end
    end
end
for i = 1:number_of_p_waves
    if I_p_wave_areas_2_CS_body_SD == 0 && I_p_wave_areas_2(i) == I_p_wave_areas_2_CS_body_average 
        I_p_wave_areas_2_z_scores_CS_body(i) = 0;
    else
        if I_p_wave_areas_2_CS_body_SD == 0 && ~(I_p_wave_areas_2(i) == I_p_wave_areas_2_CS_body_average)
            I_p_wave_areas_2_z_scores_CS_body(i) = 10*(abs(I_p_wave_areas_2(i)-I_p_wave_areas_2_CS_body_average));
        end
    end
end
for i = 1:number_of_p_waves
    if I_p_wave_areas_3_CS_body_SD == 0 && I_p_wave_areas_3(i) == I_p_wave_areas_3_CS_body_average 
        I_p_wave_areas_3_z_scores_CS_body(i) = 0;
    else
        if I_p_wave_areas_3_CS_body_SD == 0 && ~(I_p_wave_areas_3(i) == I_p_wave_areas_3_CS_body_average)
            I_p_wave_areas_3_z_scores_CS_body(i) = 10*(abs(I_p_wave_areas_3(i)-I_p_wave_areas_3_CS_body_average));
        end
    end
end
for i = 1:number_of_p_waves
    if I_p_wave_areas_4_CS_body_SD == 0 && I_p_wave_areas_4(i) == I_p_wave_areas_4_CS_body_average 
        I_p_wave_areas_4_z_scores_CS_body(i) = 0;
    else
        if I_p_wave_areas_4_CS_body_SD == 0 && ~(I_p_wave_areas_4(i) == I_p_wave_areas_4_CS_body_average)
            I_p_wave_areas_4_z_scores_CS_body(i) = 10*(abs(I_p_wave_areas_4(i)-I_p_wave_areas_4_CS_body_average));
        end
    end
end
for i = 1:number_of_p_waves
    if I_p_wave_areas_5_CS_body_SD == 0 && I_p_wave_areas_5(i) == I_p_wave_areas_5_CS_body_average 
        I_p_wave_areas_5_z_scores_CS_body(i) = 0;
    else
        if I_p_wave_areas_5_CS_body_SD == 0 && ~(I_p_wave_areas_5(i) == I_p_wave_areas_5_CS_body_average)
            I_p_wave_areas_5_z_scores_CS_body(i) = 10*(abs(I_p_wave_areas_5(i)-I_p_wave_areas_5_CS_body_average));
        end
    end
end
for i = 1:number_of_p_waves
    if I_num_zeros_left_septum_SD == 0 && I_number_of_zeros_array(i) == I_num_zeros_left_septum_average 
        I_num_zeros_z_scores_left_septum(i) = 0;
    else
        if I_num_zeros_left_septum_SD == 0 && ~(I_number_of_zeros_array(i) == I_num_zeros_left_septum_average)
            I_num_zeros_z_scores_left_septum(i) = 10*(abs(I_num_zeros_z_scores_left_septum(i)-I_num_zeros_left_septum_average));
        end
    end
end
for i = 1:number_of_p_waves
    if I_maxima_left_septum_SD == 0 && I_maxima(i) == I_maxima_left_septum_average 
        I_maxima_z_scores_left_septum(i) = 0;
    else
        if I_maxima_left_septum_SD == 0 && ~(I_maxima(i) == I_maxima_left_septum_average)
            I_maxima_z_scores_left_septum(i) = 10*(abs(I_maxima(i)-I_maxima_left_septum_average));
        end
    end
end
for i = 1:number_of_p_waves
    if I_minima_left_septum_SD == 0 && I_minima(i) == I_minima_left_septum_average 
        I_minima_z_scores_left_septum(i) = 0;
    else
        if I_minima_left_septum_SD == 0 && ~(I_minima(i) == I_minima_left_septum_average)
            I_minima_z_scores_left_septum(i) = 10*(abs(I_minima(i)-I_minima_left_septum_average));
        end
    end
end
for i = 1:number_of_p_waves
    if I_sym_inds_left_septum_SD == 0 && I_sym_inds(i) == I_sym_inds_left_septum_average 
        I_sym_inds_z_scores_left_septum(i) = 0;
    else
        if I_sym_inds_left_septum_SD == 0 && ~(I_sym_inds(i) == I_sym_inds_left_septum_average)
            I_sym_inds_z_scores_left_septum(i) = 10*(abs(I_sym_inds(i)-I_sym_inds_left_septum_average));
        end
    end
end
for i = 1:number_of_p_waves
    if I_concavs_left_septum_SD == 0 && I_concavs(i) == I_concavs_left_septum_average 
        I_concavs_z_scores_left_septum(i) = 0;
    else
        if I_concavs_left_septum_SD == 0 && ~(I_concavs(i) == I_concavs_left_septum_average)
            I_concavs_z_scores_left_septum(i) = 10*(abs(I_concavs(i)-I_concavs_left_septum_average));
        end
    end
end
for i = 1:number_of_p_waves
    if I_mdslpes_left_septum_SD == 0 && I_mdslpes(i) == I_mdslpes_left_septum_average 
        I_mdslpes_z_scores_left_septum(i) = 0;
    else
        if I_mdslpes_left_septum_SD == 0 && ~(I_mdslpes(i) == I_mdslpes_left_septum_average)
            I_mdslpes_z_scores_left_septum(i) = 10*(abs(I_mdslpes(i)-I_mdslpes_left_septum_average));
        end
    end
end
for i = 1:number_of_p_waves
    if I_p_wave_areas_1_left_septum_SD == 0 && I_p_wave_areas_1(i) == I_p_wave_areas_1_left_septum_average 
        I_p_wave_areas_1_z_scores_left_septum(i) = 0;
    else
        if I_p_wave_areas_1_left_septum_SD == 0 && ~(I_p_wave_areas_1(i) == I_p_wave_areas_1_left_septum_average)
            I_p_wave_areas_1_z_scores_left_septum(i) = 10*(abs(I_p_wave_areas_1(i)-I_p_wave_areas_1_left_septum_average));
        end
    end
end
for i = 1:number_of_p_waves
    if I_p_wave_areas_2_left_septum_SD == 0 && I_p_wave_areas_2(i) == I_p_wave_areas_2_left_septum_average 
        I_p_wave_areas_2_z_scores_left_septum(i) = 0;
    else
        if I_p_wave_areas_2_left_septum_SD == 0 && ~(I_p_wave_areas_2(i) == I_p_wave_areas_2_left_septum_average)
            I_p_wave_areas_2_z_scores_left_septum(i) = 10*(abs(I_p_wave_areas_2(i)-I_p_wave_areas_2_left_septum_average));
        end
    end
end
for i = 1:number_of_p_waves
    if I_p_wave_areas_3_left_septum_SD == 0 && I_p_wave_areas_3(i) == I_p_wave_areas_3_left_septum_average 
        I_p_wave_areas_3_z_scores_left_septum(i) = 0;
    else
        if I_p_wave_areas_3_left_septum_SD == 0 && ~(I_p_wave_areas_3(i) == I_p_wave_areas_3_left_septum_average)
            I_p_wave_areas_3_z_scores_left_septum(i) = 10*(abs(I_p_wave_areas_3(i)-I_p_wave_areas_3_left_septum_average));
        end
    end
end
for i = 1:number_of_p_waves
    if I_p_wave_areas_4_left_septum_SD == 0 && I_p_wave_areas_4(i) == I_p_wave_areas_4_left_septum_average 
        I_p_wave_areas_4_z_scores_left_septum(i) = 0;
    else
        if I_p_wave_areas_4_left_septum_SD == 0 && ~(I_p_wave_areas_4(i) == I_p_wave_areas_4_left_septum_average)
            I_p_wave_areas_4_z_scores_left_septum(i) = 10*(abs(I_p_wave_areas_4(i)-I_p_wave_areas_4_left_septum_average));
        end
    end
end
for i = 1:number_of_p_waves
    if I_p_wave_areas_5_left_septum_SD == 0 && I_p_wave_areas_5(i) == I_p_wave_areas_5_left_septum_average 
        I_p_wave_areas_5_z_scores_left_septum(i) = 0;
    else
        if I_p_wave_areas_5_left_septum_SD == 0 && ~(I_p_wave_areas_5(i) == I_p_wave_areas_5_left_septum_average)
            I_p_wave_areas_5_z_scores_left_septum(i) = 10*(abs(I_p_wave_areas_5(i)-I_p_wave_areas_5_left_septum_average));
        end
    end
end
for i = 1:number_of_p_waves
    if I_num_zeros_left_atrial_appendage_SD == 0 && I_number_of_zeros_array(i) == I_num_zeros_left_atrial_appendage_average 
        I_num_zeros_z_scores_left_atrial_appendage(i) = 0;
    else
        if I_num_zeros_left_atrial_appendage_SD == 0 && ~(I_number_of_zeros_array(i) == I_num_zeros_left_atrial_appendage_average)
            I_num_zeros_z_scores_left_atrial_appendage(i) = 10*(abs(I_num_zeros_z_scores_left_atrial_appendage(i)-I_num_zeros_left_atrial_appendage_average));
        end
    end
end
for i = 1:number_of_p_waves
    if I_maxima_left_atrial_appendage_SD == 0 && I_maxima(i) == I_maxima_left_atrial_appendage_average 
        I_maxima_z_scores_left_atrial_appendage(i) = 0;
    else
        if I_maxima_left_atrial_appendage_SD == 0 && ~(I_maxima(i) == I_maxima_left_atrial_appendage_average)
            I_maxima_z_scores_left_atrial_appendage(i) = 10*(abs(I_maxima(i)-I_maxima_left_atrial_appendage_average));
        end
    end
end
for i = 1:number_of_p_waves
    if I_minima_left_atrial_appendage_SD == 0 && I_minima(i) == I_minima_left_atrial_appendage_average 
        I_minima_z_scores_left_atrial_appendage(i) = 0;
    else
        if I_minima_left_atrial_appendage_SD == 0 && ~(I_minima(i) == I_minima_left_atrial_appendage_average)
            I_minima_z_scores_left_atrial_appendage(i) = 10*(abs(I_minima(i)-I_minima_left_atrial_appendage_average));
        end
    end
end
for i = 1:number_of_p_waves
    if I_sym_inds_left_atrial_appendage_SD == 0 && I_sym_inds(i) == I_sym_inds_left_atrial_appendage_average 
        I_sym_inds_z_scores_left_atrial_appendage(i) = 0;
    else
        if I_sym_inds_left_atrial_appendage_SD == 0 && ~(I_sym_inds(i) == I_sym_inds_left_atrial_appendage_average)
            I_sym_inds_z_scores_left_atrial_appendage(i) = 10*(abs(I_sym_inds(i)-I_sym_inds_left_atrial_appendage_average));
        end
    end
end
for i = 1:number_of_p_waves
    if I_concavs_left_atrial_appendage_SD == 0 && I_concavs(i) == I_concavs_left_atrial_appendage_average 
        I_concavs_z_scores_left_atrial_appendage(i) = 0;
    else
        if I_concavs_left_atrial_appendage_SD == 0 && ~(I_concavs(i) == I_concavs_left_atrial_appendage_average)
            I_concavs_z_scores_left_atrial_appendage(i) = 10*(abs(I_concavs(i)-I_concavs_left_atrial_appendage_average));
        end
    end
end
for i = 1:number_of_p_waves
    if I_mdslpes_left_atrial_appendage_SD == 0 && I_mdslpes(i) == I_mdslpes_left_atrial_appendage_average 
        I_mdslpes_z_scores_left_atrial_appendage(i) = 0;
    else
        if I_mdslpes_left_atrial_appendage_SD == 0 && ~(I_mdslpes(i) == I_mdslpes_left_atrial_appendage_average)
            I_mdslpes_z_scores_left_atrial_appendage(i) = 10*(abs(I_mdslpes(i)-I_mdslpes_left_atrial_appendage_average));
        end
    end
end
for i = 1:number_of_p_waves
    if I_p_wave_areas_1_left_atrial_appendage_SD == 0 && I_p_wave_areas_1(i) == I_p_wave_areas_1_left_atrial_appendage_average 
        I_p_wave_areas_1_z_scores_left_atrial_appendage(i) = 0;
    else
        if I_p_wave_areas_1_left_atrial_appendage_SD == 0 && ~(I_p_wave_areas_1(i) == I_p_wave_areas_1_left_atrial_appendage_average)
            I_p_wave_areas_1_z_scores_left_atrial_appendage(i) = 10*(abs(I_p_wave_areas_1(i)-I_p_wave_areas_1_left_atrial_appendage_average));
        end
    end
end
for i = 1:number_of_p_waves
    if I_p_wave_areas_2_left_atrial_appendage_SD == 0 && I_p_wave_areas_2(i) == I_p_wave_areas_2_left_atrial_appendage_average 
        I_p_wave_areas_2_z_scores_left_atrial_appendage(i) = 0;
    else
        if I_p_wave_areas_2_left_atrial_appendage_SD == 0 && ~(I_p_wave_areas_2(i) == I_p_wave_areas_2_left_atrial_appendage_average)
            I_p_wave_areas_2_z_scores_left_atrial_appendage(i) = 10*(abs(I_p_wave_areas_2(i)-I_p_wave_areas_2_left_atrial_appendage_average));
        end
    end
end
for i = 1:number_of_p_waves
    if I_p_wave_areas_3_left_atrial_appendage_SD == 0 && I_p_wave_areas_3(i) == I_p_wave_areas_3_left_atrial_appendage_average 
        I_p_wave_areas_3_z_scores_left_atrial_appendage(i) = 0;
    else
        if I_p_wave_areas_3_left_atrial_appendage_SD == 0 && ~(I_p_wave_areas_3(i) == I_p_wave_areas_3_left_atrial_appendage_average)
            I_p_wave_areas_3_z_scores_left_atrial_appendage(i) = 10*(abs(I_p_wave_areas_3(i)-I_p_wave_areas_3_left_atrial_appendage_average));
        end
    end
end
for i = 1:number_of_p_waves
    if I_p_wave_areas_4_left_atrial_appendage_SD == 0 && I_p_wave_areas_4(i) == I_p_wave_areas_4_left_atrial_appendage_average 
        I_p_wave_areas_4_z_scores_left_atrial_appendage(i) = 0;
    else
        if I_p_wave_areas_4_left_atrial_appendage_SD == 0 && ~(I_p_wave_areas_4(i) == I_p_wave_areas_4_left_atrial_appendage_average)
            I_p_wave_areas_4_z_scores_left_atrial_appendage(i) = 10*(abs(I_p_wave_areas_4(i)-I_p_wave_areas_4_left_atrial_appendage_average));
        end
    end
end
for i = 1:number_of_p_waves
    if I_p_wave_areas_5_left_atrial_appendage_SD == 0 && I_p_wave_areas_5(i) == I_p_wave_areas_5_left_atrial_appendage_average 
        I_p_wave_areas_5_z_scores_left_atrial_appendage(i) = 0;
    else
        if I_p_wave_areas_5_left_atrial_appendage_SD == 0 && ~(I_p_wave_areas_5(i) == I_p_wave_areas_5_left_atrial_appendage_average)
            I_p_wave_areas_5_z_scores_left_atrial_appendage(i) = 10*(abs(I_p_wave_areas_5(i)-I_p_wave_areas_5_left_atrial_appendage_average));
        end
    end
end
row_vec_sum_isNormal = sum(isNormal);
row_vec_sum_isSpike = sum(isSpike);
sum_isNormal = sum(row_vec_sum_isNormal);
sum_isSpike = sum(row_vec_sum_isSpike);
fraction_of_dist_normal = sum_isNormal/numel(isNormal);
fraction_of_dist_spike = sum_isSpike/numel(isSpike);
%Below are arrays that determine the location of the ectopic pacemaker
%according to each parameter.
ectopic_pacemaker_location_indices_num_zeros = zeros(length(p_wave_start_indices),1);
ectopic_pacemaker_location_indices_maxima = zeros(length(p_wave_start_indices),1);
ectopic_pacemaker_location_indices_minima = zeros(length(p_wave_start_indices),1);
ectopic_pacemaker_location_indices_sym_inds = zeros(length(p_wave_start_indices),1);
ectopic_pacemaker_location_indices_concavs = zeros(length(p_wave_start_indices),1);
ectopic_pacemaker_location_indices_mdslpes = zeros(length(p_wave_start_indices),1);
ectopic_pacemaker_location_indices_p_wave_areas_1 = zeros(length(p_wave_start_indices),1);
ectopic_pacemaker_location_indices_p_wave_areas_2 = zeros(length(p_wave_start_indices),1);
ectopic_pacemaker_location_indices_p_wave_areas_3 = zeros(length(p_wave_start_indices),1);
ectopic_pacemaker_location_indices_p_wave_areas_4 = zeros(length(p_wave_start_indices),1);
ectopic_pacemaker_location_indices_p_wave_areas_5 = zeros(length(p_wave_start_indices),1);
%Figure out what to name these. Figure out a plan before coding
location_counter_region1 = zeros(length(p_wave_start_indices),1);
location_counter_region2 = zeros(length(p_wave_start_indices),1);
location_counter_region3 = zeros(length(p_wave_start_indices),1);
location_counter_region4 = zeros(length(p_wave_start_indices),1);
location_counter_region5 = zeros(length(p_wave_start_indices),1);
location_counter_region6 = zeros(length(p_wave_start_indices),1);
location_counter_region7 = zeros(length(p_wave_start_indices),1);
location_counter_region8 = zeros(length(p_wave_start_indices),1);
location_counter_region9 = zeros(length(p_wave_start_indices),1);
location_counter_region10 = zeros(length(p_wave_start_indices),1);
location_counter_region11 = zeros(length(p_wave_start_indices),1);
location_different_regions_counter = zeros(length(p_wave_start_indices),1);
%See above.
%Answer: Create arrays with the lengths of the number of p-waves and put a
%counter for each region.
%The latter of these two counters determines how many regions the algorithm
%simultaneously determines to be the region of origin for a given p-wave.
%The former determines how many measurements indicate the likely region of
%origin.
number_of_p_waves = length(p_wave_start_indices);
number_of_parameters = 11;
number_of_regions = 13;
isNormal = zeros(number_of_parameters, number_of_regions);
z_score_threshold_difference = 1;
significance_level = 2;
number_of_possible_regions = 13;
for i = 1:number_of_p_waves
    if ~(isNormal(1,1)) && ~(isSpike(1,1))
        I_num_zeros_z_scores_SA_node(i) = 1000;
    end
    if ~(isNormal(1,2)) && ~(isSpike(1,2))
        I_num_zeros_z_scores_crista_terminalis(i) = 1000;
    end
    if ~(isNormal(1,3)) && ~(isSpike(1,3))
        I_num_zeros_z_scores_tricuspid_annulus(i) = 1000;
    end
    if ~(isNormal(1,4)) && ~(isSpike(1,4))
        I_num_zeros_z_scores_coronary_sinus(i) = 1000;
    end
    if ~(isNormal(1,5)) && ~(isSpike(1,5))
        I_num_zeros_z_scores_ostium(i) = 1000;
    end
    if ~(isNormal(1,6)) && ~(isSpike(1,6))
        I_num_zeros_z_scores_perinodal(i) = 1000;
    end
    if ~(isNormal(1,7)) && ~(isSpike(1,7))
        I_num_zeros_z_scores_right_septum(i) = 1000;
    end
    if ~(isNormal(1,8)) && ~(isSpike(1,8))
        I_num_zeros_z_scores_right_atrial_appendage(i) = 1000;
    end
    if ~(isNormal(1,9)) && ~(isSpike(1,9))
        I_num_zeros_z_scores_pulmonary_veins(i) = 1000;
    end
    if ~(isNormal(1,10)) && ~(isSpike(1,10))
        I_num_zeros_z_scores_mitral_annulus(i) = 1000;
    end
    if ~(isNormal(1,11)) && ~(isSpike(1,11))
        I_num_zeros_z_scores_CS_body(i) = 1000;
    end
    if ~(isNormal(1,12)) && ~(isSpike(1,12))
        I_num_zeros_z_scores_left_septum(i) = 1000;
    end
    if ~(isNormal(1,13)) && ~(isSpike(1,13))
        I_num_zeros_z_scores_left_atrial_appendage(i) = 1000;
    end
    if ~(isNormal(2,1)) 
        I_maxima_z_scores_SA_node(i) = 1000;
    end
    if ~(isNormal(2,2))
        I_maxima_z_scores_crista_terminalis(i) = 1000;
    end
    if ~(isNormal(2,3))
        I_maxima_z_scores_tricuspid_annulus(i) = 1000;
    end
    if ~(isNormal(2,4))
        I_maxima_z_scores_coronary_sinus(i) = 1000;
    end
    if ~(isNormal(2,5))
        I_maxima_z_scores_ostium(i) = 1000;
    end
    if ~(isNormal(2,6))
        I_maxima_z_scores_perinodal(i) = 1000;
    end
    if ~(isNormal(2,7))
        I_maxima_z_scores_right_septum(i) = 1000;
    end
    if ~(isNormal(2,8))
        I_maxima_z_scores_right_atrial_appendage(i) = 1000;
    end
    if ~(isNormal(2,9))
        I_maxima_z_scores_pulmonary_veins(i) = 1000;
    end
    if ~(isNormal(2,10))
        I_maxima_z_scores_mitral_annulus(i) = 1000;
    end
    if ~(isNormal(2,11))
        I_maxima_z_scores_CS_body(i) = 1000;
    end
    if ~(isNormal(2,12))
        I_maxima_z_scores_left_septum(i) = 1000;
    end
    if ~(isNormal(2,13))
        I_maxima_z_scores_left_atrial_appendage(i) = 1000;
    end
    if ~(isNormal(3,1))
        I_minima_z_scores_SA_node(i) = 1000;
    end
    if ~(isNormal(3,2))
        I_minima_z_scores_crista_terminalis(i) = 1000;
    end
    if ~(isNormal(3,3))
        I_minima_z_scores_tricuspid_annulus(i) = 1000;
    end
    if ~(isNormal(3,4))
        I_minima_z_scores_coronary_sinus(i) = 1000;
    end
    if ~(isNormal(3,5))
        I_minima_z_scores_ostium(i) = 1000;
    end
    if ~(isNormal(3,6))
        I_num_zeros_z_scores_perinodal(i) = 1000;
    end
    if ~(isNormal(3,7))
        I_minima_z_scores_right_septum(i) = 1000;
    end
    if ~(isNormal(3,8))
        I_minima_z_scores_right_atrial_appendage(i) = 1000;
    end
    if ~(isNormal(3,9))
        I_minima_z_scores_pulmonary_veins(i) = 1000;
    end
    if ~(isNormal(3,10))
        I_minima_z_scores_mitral_annulus(i) = 1000;
    end
    if ~(isNormal(3,11))
        I_minima_z_scores_CS_body(i) = 1000;
    end
    if ~(isNormal(3,12))
        I_minima_z_scores_left_septum(i) = 1000;
    end
    if ~(isNormal(3,13))
        I_minima_z_scores_left_atrial_appendage(i) = 1000;
    end
    if ~(isNormal(4,1))
        I_sym_inds_z_scores_SA_node(i) = 1000;
    end
    if ~(isNormal(4,2))
        I_sym_inds_z_scores_crista_terminalis(i) = 1000;
    end
    if ~(isNormal(4,3))
        I_sym_inds_z_scores_tricuspid_annulus(i) = 1000;
    end
    if ~(isNormal(4,4))
        I_sym_inds_z_scores_coronary_sinus(i) = 1000;
    end
    if ~(isNormal(4,5))
        I_sym_inds_z_scores_ostium(i) = 1000;
    end
    if ~(isNormal(4,6))
        I_sym_inds_z_scores_perinodal(i) = 1000;
    end
    if ~(isNormal(4,7))
        I_sym_inds_z_scores_right_septum(i) = 1000;
    end
    if ~(isNormal(4,8))
        I_sym_inds_z_scores_right_atrial_appendage(i) = 1000;
    end
    if ~(isNormal(4,9))
        I_sym_inds_z_scores_pulmonary_veins(i) = 1000;
    end
    if ~(isNormal(4,10))
        I_sym_inds_z_scores_mitral_annulus(i) = 1000;
    end
    if ~(isNormal(4,11))
        I_sym_inds_z_scores_CS_body(i) = 1000;
    end
    if ~(isNormal(4,12))
        I_sym_inds_z_scores_left_septum(i) = 1000;
    end
    if ~(isNormal(4,13))
        I_sym_inds_z_scores_left_atrial_appendage(i) = 1000;
    end
    if ~(isNormal(5,1))
        I_concavs_z_scores_SA_node(i) = 1000;
    end
    if ~(isNormal(5,2))
        I_concavs_z_scores_crista_terminalis(i) = 1000;
    end
    if ~(isNormal(5,3))
        I_concavs_z_scores_tricuspid_annulus(i) = 1000;
    end
    if ~(isNormal(5,4))
        I_concavs_z_scores_coronary_sinus(i) = 1000;
    end
    if ~(isNormal(5,5))
        I_concavs_z_scores_ostium(i) = 1000;
    end
    if ~(isNormal(5,6))
        I_concavs_z_scores_perinodal(i) = 1000;
    end
    if ~(isNormal(5,7))
        I_concavs_z_scores_right_septum(i) = 1000;
    end
    if ~(isNormal(5,8))
        I_concavs_z_scores_right_atrial_appendage(i) = 1000;
    end
    if ~(isNormal(5,9))
        I_concavs_z_scores_pulmonary_veins(i) = 1000;
    end
    if ~(isNormal(5,10))
        I_concavs_z_scores_mitral_annulus(i) = 1000;
    end
    if ~(isNormal(5,11))
        I_concavs_z_scores_CS_body(i) = 1000;
    end
    if ~(isNormal(5,12))
        I_concavs_z_scores_left_septum(i) = 1000;
    end
    if ~(isNormal(5,13))
        I_concavs_z_scores_left_atrial_appendage(i) = 1000;
    end
    if ~(isNormal(6,1))
        I_mdslpes_z_scores_SA_node(i) = 1000;
    end
    if ~(isNormal(6,2))
        I_mdslpes_z_scores_crista_terminalis(i) = 1000;
    end
    if ~(isNormal(6,3))
        I_mdslpes_z_scores_tricuspid_annulus(i) = 1000;
    end
    if ~(isNormal(6,4))
        I_mdslpes_z_scores_coronary_sinus(i) = 1000;
    end
    if ~(isNormal(6,5))
        I_mdslpes_z_scores_ostium(i) = 1000;
    end
    if ~(isNormal(6,6))
        I_mdslpes_z_scores_perinodal(i) = 1000;
    end
    if ~(isNormal(6,7))
        I_mdslpes_z_scores_right_septum(i) = 1000;
    end
    if ~(isNormal(6,8))
        I_mdslpes_z_scores_right_atrial_appendage(i) = 1000;
    end
    if ~(isNormal(6,9))
        I_mdslpes_z_scores_pulmonary_veins(i) = 1000;
    end
    if ~(isNormal(6,10))
        I_mdslpes_z_scores_mitral_annulus(i) = 1000;
    end
    if ~(isNormal(6,11))
        I_mdslpes_z_scores_CS_body(i) = 1000;
    end
    if ~(isNormal(6,12))
        I_mdslpes_z_scores_left_septum(i) = 1000;
    end
    if ~(isNormal(6,13))
        I_mdslpes_z_scores_left_atrial_appendage(i) = 1000;
    end
    if ~(isNormal(7,1))
        I_p_wave_areas_1_z_scores_SA_node(i) = 1000;
    end
    if ~(isNormal(7,2))
        I_p_wave_areas_1_z_scores_crista_terminalis(i) = 1000;
    end
    if ~(isNormal(7,3))
        I_p_wave_areas_1_z_scores_tricuspid_annulus(i) = 1000;
    end
    if ~(isNormal(7,4))
        I_p_wave_areas_1_z_scores_coronary_sinus(i) = 1000;
    end
    if ~(isNormal(7,5))
        I_p_wave_areas_1_z_scores_ostium(i) = 1000;
    end
    if ~(isNormal(7,6))
        I_p_wave_areas_1_z_scores_perinodal(i) = 1000;
    end
    if ~(isNormal(7,7))
        I_p_wave_areas_1_z_scores_right_septum(i) = 1000;
    end
    if ~(isNormal(7,8))
        I_p_wave_areas_1_z_scores_right_atrial_appendage(i) = 1000;
    end
    if ~(isNormal(7,9))
        I_p_wave_areas_1_z_scores_pulmonary_veins(i) = 1000;
    end
    if ~(isNormal(7,10))
        I_p_wave_areas_1_z_scores_mitral_annulus(i) = 1000;
    end
    if ~(isNormal(7,11))
        I_p_wave_areas_1_z_scores_CS_body(i) = 1000;
    end
    if ~(isNormal(7,12))
        I_p_wave_areas_1_z_scores_left_septum(i) = 1000;
    end
    if ~(isNormal(7,13))
        I_p_wave_areas_1_z_scores_left_atrial_appendage(i) = 1000;
    end
    if ~(isNormal(8,1))
        I_p_wave_areas_2_z_scores_SA_node(i) = 1000;
    end
    if ~(isNormal(8,2))
        I_p_wave_areas_2_z_scores_crista_terminalis(i) = 1000;
    end
    if ~(isNormal(8,3))
        I_p_wave_areas_2_z_scores_tricuspid_annulus(i) = 1000;
    end
    if ~(isNormal(8,4))
        I_p_wave_areas_2_z_scores_coronary_sinus(i) = 1000;
    end
    if ~(isNormal(8,5))
        I_p_wave_areas_2_z_scores_ostium(i) = 1000;
    end
    if ~(isNormal(8,6))
        I_p_wave_areas_2_z_scores_perinodal(i) = 1000;
    end
    if ~(isNormal(8,7))
        I_p_wave_areas_2_z_scores_right_septum(i) = 1000;
    end
    if ~(isNormal(8,8))
        I_p_wave_areas_2_z_scores_right_atrial_appendage(i) = 1000;
    end
    if ~(isNormal(8,9))
        I_p_wave_areas_2_z_scores_pulmonary_veins(i) = 1000;
    end
    if ~(isNormal(8,10))
        I_p_wave_areas_2_z_scores_mitral_annulus(i) = 1000;
    end
    if ~(isNormal(8,11))
        I_p_wave_areas_2_z_scores_CS_body(i) = 1000;
    end
    if ~(isNormal(8,12))
        I_p_wave_areas_2_z_scores_left_septum(i) = 1000;
    end
    if ~(isNormal(8,13))
        I_p_wave_areas_2_z_scores_left_atrial_appendage(i) = 1000;
    end
    if ~(isNormal(9,1))
        I_p_wave_areas_3_z_scores_SA_node(i) = 1000;
    end
    if ~(isNormal(9,2))
        I_p_wave_areas_3_z_scores_crista_terminalis(i) = 1000;
    end
    if ~(isNormal(9,3))
        I_p_wave_areas_3_z_scores_tricuspid_annulus(i) = 1000;
    end
    if ~(isNormal(9,4))
        I_p_wave_areas_3_z_scores_coronary_sinus(i) = 1000;
    end
    if ~(isNormal(9,5))
        I_p_wave_areas_3_z_scores_ostium(i) = 1000;
    end
    if ~(isNormal(9,6))
        I_p_wave_areas_3_z_scores_perinodal(i) = 1000;
    end
    if ~(isNormal(9,7))
        I_p_wave_areas_3_z_scores_right_septum(i) = 1000;
    end
    if ~(isNormal(9,8))
        I_p_wave_areas_3_z_scores_right_atrial_appendage(i) = 1000;
    end
    if ~(isNormal(9,9))
        I_p_wave_areas_3_z_scores_pulmonary_veins(i) = 1000;
    end
    if ~(isNormal(9,10))
        I_p_wave_areas_3_z_scores_mitral_annulus(i) = 1000;
    end
    if ~(isNormal(9,11))
        I_p_wave_areas_3_z_scores_CS_body(i) = 1000;
    end
    if ~(isNormal(9,12))
        I_p_wave_areas_3_z_scores_left_septum(i) = 1000;
    end
    if ~(isNormal(9,13))
        I_p_wave_areas_3_z_scores_left_atrial_appendage(i) = 1000;
    end
    if ~(isNormal(10,1))
        I_p_wave_areas_4_z_scores_SA_node(i) = 1000;
    end
    if ~(isNormal(10,2))
        I_p_wave_areas_4_z_scores_crista_terminalis(i) = 1000;
    end
    if ~(isNormal(10,3))
        I_p_wave_areas_4_z_scores_tricuspid_annulus(i) = 1000;
    end
    if ~(isNormal(10,4))
        I_p_wave_areas_4_z_scores_coronary_sinus(i) = 1000;
    end
    if ~(isNormal(10,5))
        I_p_wave_areas_4_z_scores_ostium(i) = 1000;
    end
    if ~(isNormal(10,6))
        I_p_wave_areas_4_z_scores_perinodal(i) = 1000;
    end
    if ~(isNormal(10,7))
        I_p_wave_areas_4_z_scores_right_septum(i) = 1000;
    end
    if ~(isNormal(10,8))
        I_p_wave_areas_4_z_scores_right_atrial_appendage(i) = 1000;
    end
    if ~(isNormal(10,9))
        I_p_wave_areas_4_z_scores_pulmonary_veins(i) = 1000;
    end
    if ~(isNormal(10,10))
        I_p_wave_areas_4_z_scores_mitral_annulus(i) = 1000;
    end
    if ~(isNormal(10,11))
        I_p_wave_areas_4_z_scores_CS_body(i) = 1000;
    end
    if ~(isNormal(10,12))
        I_p_wave_areas_4_z_scores_left_septum(i) = 1000;
    end
    if ~(isNormal(10,13))
        I_p_wave_areas_4_z_scores_left_atrial_appendage(i) = 1000;
    end
    if ~(isNormal(11,1))
        I_p_wave_areas_5_z_scores_SA_node(i) = 1000;
    end
    if ~(isNormal(11,2))
        I_p_wave_areas_5_z_scores_crista_terminalis(i) = 1000;
    end
    if ~(isNormal(11,3))
        I_p_wave_areas_5_z_scores_tricuspid_annulus(i) = 1000;
    end
    if ~(isNormal(11,4))
        I_p_wave_areas_5_z_scores_coronary_sinus(i) = 1000;
    end
    if ~(isNormal(11,5))
        I_p_wave_areas_5_z_scores_ostium(i) = 1000;
    end
    if ~(isNormal(11,6))
        I_p_wave_areas_5_z_scores_perinodal(i) = 1000;
    end
    if ~(isNormal(11,7))
        I_p_wave_areas_5_z_scores_right_septum(i) = 1000;
    end
    if ~(isNormal(11,8))
        I_p_wave_areas_5_z_scores_right_atrial_appendage(i) = 1000;
    end
    if ~(isNormal(11,9))
        I_p_wave_areas_5_z_scores_pulmonary_veins(i) = 1000;
    end
    if ~(isNormal(11,10))
        I_p_wave_areas_5_z_scores_mitral_annulus(i) = 1000;
    end
    if ~(isNormal(11,11))
        I_p_wave_areas_5_z_scores_CS_body(i) = 1000;
    end
    if ~(isNormal(11,12))
        I_p_wave_areas_5_z_scores_left_septum(i) = 1000;
    end
    if ~(isNormal(11,13))
        I_p_wave_areas_5_z_scores_left_atrial_appendage(i) = 1000;
    end
end
for i = 1:number_of_p_waves
    I_num_zeros_z_scores_regions_array = [I_num_zeros_z_scores_SA_node(i),I_num_zeros_z_scores_crista_terminalis(i),I_num_zeros_z_scores_tricuspid_annulus(i),I_num_zeros_z_scores_coronary_sinus(i),I_num_zeros_z_scores_ostium(i),I_num_zeros_z_scores_perinodal(i),I_num_zeros_z_scores_right_septum(i),I_num_zeros_z_scores_right_atrial_appendage(i),I_num_zeros_z_scores_pulmonary_veins(i),I_num_zeros_z_scores_mitral_annulus(i),I_num_zeros_z_scores_CS_body(i),I_num_zeros_z_scores_left_septum(i),I_num_zeros_z_scores_left_atrial_appendage(i)];
    I_maxima_z_scores_regions_array = [I_maxima_z_scores_SA_node(i),I_maxima_z_scores_crista_terminalis(i),I_maxima_z_scores_tricuspid_annulus(i),I_maxima_z_scores_coronary_sinus(i),I_maxima_z_scores_ostium(i),I_maxima_z_scores_perinodal(i),I_maxima_z_scores_right_septum(i),I_maxima_z_scores_right_atrial_appendage(i),I_maxima_z_scores_pulmonary_veins(i),I_maxima_z_scores_mitral_annulus(i),I_maxima_z_scores_CS_body(i),I_maxima_z_scores_left_septum(i),I_maxima_z_scores_left_atrial_appendage(i)];
    I_minima_z_scores_regions_array = [I_minima_z_scores_SA_node(i),I_minima_z_scores_crista_terminalis(i),I_minima_z_scores_tricuspid_annulus(i),I_minima_z_scores_coronary_sinus(i),I_minima_z_scores_ostium(i),I_minima_z_scores_perinodal(i),I_minima_z_scores_right_septum(i),I_minima_z_scores_right_atrial_appendage(i),I_minima_z_scores_pulmonary_veins(i),I_minima_z_scores_mitral_annulus(i),I_minima_z_scores_CS_body(i),I_minima_z_scores_left_septum(i),I_minima_z_scores_left_atrial_appendage(i)];
    I_sym_inds_z_scores_regions_array = [I_sym_inds_z_scores_SA_node(i),I_sym_inds_z_scores_crista_terminalis(i),I_sym_inds_z_scores_tricuspid_annulus(i),I_sym_inds_z_scores_coronary_sinus(i),I_sym_inds_z_scores_ostium(i),I_sym_inds_z_scores_perinodal(i),I_sym_inds_z_scores_right_septum(i),I_sym_inds_z_scores_right_atrial_appendage(i),I_sym_inds_z_scores_pulmonary_veins(i),I_sym_inds_z_scores_mitral_annulus(i),I_sym_inds_z_scores_CS_body(i),I_sym_inds_z_scores_left_septum(i),I_sym_inds_z_scores_left_atrial_appendage(i)];
    I_concavs_z_scores_regions_array = [I_concavs_z_scores_SA_node(i),I_concavs_z_scores_crista_terminalis(i),I_concavs_z_scores_tricuspid_annulus(i),I_concavs_z_scores_coronary_sinus(i),I_concavs_z_scores_ostium(i),I_concavs_z_scores_perinodal(i),I_concavs_z_scores_right_septum(i),I_concavs_z_scores_right_atrial_appendage(i),I_concavs_z_scores_pulmonary_veins(i),I_concavs_z_scores_mitral_annulus(i),I_concavs_z_scores_CS_body(i),I_concavs_z_scores_left_septum(i),I_concavs_z_scores_left_atrial_appendage(i)];
    I_mdslpes_z_scores_regions_array = [I_mdslpes_z_scores_SA_node(i),I_mdslpes_z_scores_crista_terminalis(i),I_mdslpes_z_scores_tricuspid_annulus(i),I_mdslpes_z_scores_coronary_sinus(i),I_mdslpes_z_scores_ostium(i),I_mdslpes_z_scores_perinodal(i),I_mdslpes_z_scores_right_septum(i),I_mdslpes_z_scores_right_atrial_appendage(i),I_mdslpes_z_scores_pulmonary_veins(i),I_mdslpes_z_scores_mitral_annulus(i),I_mdslpes_z_scores_CS_body(i),I_mdslpes_z_scores_left_septum(i),I_mdslpes_z_scores_left_atrial_appendage(i)];
    I_p_wave_areas_1_z_scores_regions_array = [I_p_wave_areas_1_z_scores_SA_node(i),I_p_wave_areas_1_z_scores_crista_terminalis(i),I_p_wave_areas_1_z_scores_tricuspid_annulus(i),I_p_wave_areas_1_z_scores_coronary_sinus(i),I_p_wave_areas_1_z_scores_ostium(i),I_p_wave_areas_1_z_scores_perinodal(i),I_p_wave_areas_1_z_scores_right_septum(i),I_p_wave_areas_1_z_scores_right_atrial_appendage(i),I_p_wave_areas_1_z_scores_pulmonary_veins(i),I_p_wave_areas_1_z_scores_mitral_annulus(i),I_p_wave_areas_1_z_scores_CS_body(i),I_p_wave_areas_1_z_scores_left_septum(i),I_p_wave_areas_1_z_scores_left_atrial_appendage(i)];
    I_p_wave_areas_2_z_scores_regions_array = [I_p_wave_areas_2_z_scores_SA_node(i),I_p_wave_areas_2_z_scores_crista_terminalis(i),I_p_wave_areas_2_z_scores_tricuspid_annulus(i),I_p_wave_areas_2_z_scores_coronary_sinus(i),I_p_wave_areas_2_z_scores_ostium(i),I_p_wave_areas_2_z_scores_perinodal(i),I_p_wave_areas_2_z_scores_right_septum(i),I_p_wave_areas_2_z_scores_right_atrial_appendage(i),I_p_wave_areas_2_z_scores_pulmonary_veins(i),I_p_wave_areas_2_z_scores_mitral_annulus(i),I_p_wave_areas_2_z_scores_CS_body(i),I_p_wave_areas_2_z_scores_left_septum(i),I_p_wave_areas_2_z_scores_left_atrial_appendage(i)];
    I_p_wave_areas_3_z_scores_regions_array = [I_p_wave_areas_3_z_scores_SA_node(i),I_p_wave_areas_3_z_scores_crista_terminalis(i),I_p_wave_areas_3_z_scores_tricuspid_annulus(i),I_p_wave_areas_3_z_scores_coronary_sinus(i),I_p_wave_areas_3_z_scores_ostium(i),I_p_wave_areas_3_z_scores_perinodal(i),I_p_wave_areas_3_z_scores_right_septum(i),I_p_wave_areas_3_z_scores_right_atrial_appendage(i),I_p_wave_areas_3_z_scores_pulmonary_veins(i),I_p_wave_areas_3_z_scores_mitral_annulus(i),I_p_wave_areas_3_z_scores_CS_body(i),I_p_wave_areas_3_z_scores_left_septum(i),I_p_wave_areas_3_z_scores_left_atrial_appendage(i)];
    I_p_wave_areas_4_z_scores_regions_array = [I_p_wave_areas_4_z_scores_SA_node(i),I_p_wave_areas_4_z_scores_crista_terminalis(i),I_p_wave_areas_4_z_scores_tricuspid_annulus(i),I_p_wave_areas_4_z_scores_coronary_sinus(i),I_p_wave_areas_4_z_scores_ostium(i),I_p_wave_areas_4_z_scores_perinodal(i),I_p_wave_areas_4_z_scores_right_septum(i),I_p_wave_areas_4_z_scores_right_atrial_appendage(i),I_p_wave_areas_4_z_scores_pulmonary_veins(i),I_p_wave_areas_4_z_scores_mitral_annulus(i),I_p_wave_areas_4_z_scores_CS_body(i),I_p_wave_areas_4_z_scores_left_septum(i),I_p_wave_areas_4_z_scores_left_atrial_appendage(i)];
    I_p_wave_areas_5_z_scores_regions_array = [I_p_wave_areas_5_z_scores_SA_node(i),I_p_wave_areas_5_z_scores_crista_terminalis(i),I_p_wave_areas_5_z_scores_tricuspid_annulus(i),I_p_wave_areas_5_z_scores_coronary_sinus(i),I_p_wave_areas_5_z_scores_ostium(i),I_p_wave_areas_5_z_scores_perinodal(i),I_p_wave_areas_5_z_scores_right_septum(i),I_p_wave_areas_5_z_scores_right_atrial_appendage(i),I_p_wave_areas_5_z_scores_pulmonary_veins(i),I_p_wave_areas_5_z_scores_mitral_annulus(i),I_p_wave_areas_5_z_scores_CS_body(i),I_p_wave_areas_5_z_scores_left_septum(i),I_p_wave_areas_5_z_scores_left_atrial_appendage(i)];
    I_num_zeros_z_scores_regions_array_abs = abs(I_num_zeros_z_scores_regions_array);
    I_maxima_z_scores_regions_array_abs = abs(I_maxima_z_scores_regions_array);
    I_minima_z_scores_regions_array_abs = abs(I_minima_z_scores_regions_array);
    I_sym_inds_z_scores_regions_array_abs = abs(I_sym_inds_z_scores_regions_array);
    I_concavs_z_scores_regions_array_abs = abs(I_concavs_z_scores_regions_array);
    I_mdslpes_z_scores_regions_array_abs = abs(I_mdslpes_z_scores_regions_array);
    I_p_wave_areas_1_z_scores_regions_array_abs = abs(I_p_wave_areas_1_z_scores_regions_array);
    I_p_wave_areas_2_z_scores_regions_array_abs = abs(I_p_wave_areas_2_z_scores_regions_array);
    I_p_wave_areas_3_z_scores_regions_array_abs = abs(I_p_wave_areas_3_z_scores_regions_array);
    I_p_wave_areas_4_z_scores_regions_array_abs = abs(I_p_wave_areas_4_z_scores_regions_array);
    I_p_wave_areas_5_z_scores_regions_array_abs = abs(I_p_wave_areas_5_z_scores_regions_array);
    I_num_zeros_z_scores_regions_array_sorted = sort(I_num_zeros_z_scores_regions_array_abs);
    I_maxima_z_scores_regions_array_sorted = sort(I_maxima_z_scores_regions_array_abs);
    I_minima_z_scores_regions_array_sorted = sort(I_minima_z_scores_regions_array_abs);
    I_sym_inds_z_scores_regions_array_sorted = sort(I_sym_inds_z_scores_regions_array_abs);
    I_concavs_z_scores_regions_array_sorted = sort(I_concavs_z_scores_regions_array_abs);
    I_mdslpes_z_scores_regions_array_sorted = sort(I_mdslpes_z_scores_regions_array_abs);
    I_p_wave_areas_1_z_scores_regions_array_sorted = sort(I_p_wave_areas_1_z_scores_regions_array_abs);
    I_p_wave_areas_2_z_scores_regions_array_sorted = sort(I_p_wave_areas_2_z_scores_regions_array_abs);
    I_p_wave_areas_3_z_scores_regions_array_sorted = sort(I_p_wave_areas_3_z_scores_regions_array_abs);
    I_p_wave_areas_4_z_scores_regions_array_sorted = sort(I_p_wave_areas_4_z_scores_regions_array_abs);
    I_p_wave_areas_5_z_scores_regions_array_sorted = sort(I_p_wave_areas_5_z_scores_regions_array_abs);
    if (I_num_zeros_z_scores_regions_array_sorted(1) < significance_level) && (I_num_zeros_z_scores_regions_array_sorted(2) > significance_level) && ((I_num_zeros_z_scores_regions_array_sorted(2) - I_num_zeros_z_scores_regions_array_sorted(2)) > z_score_threshold_difference)
        for j = 1:number_of_possible_regions
            if I_num_zeros_z_scores_regions_array_sorted(1) == I_num_zeros_z_scores_regions_array_abs(j)
                ectopic_pacemaker_location_indices_num_zeros(i) = j;
            end
            %This for loop determines the region number for this
            %parameter/lead
        end
        location_counter_region1(i) = location_counter_region1(i) + 1;
    else
        ectopic_pacemaker_location_indices_num_zeros(i) = 0;
        %this means that the first metric could not be used
    end
    if (I_maxima_z_scores_regions_array_sorted(1) < significance_level) && (I_maxima_z_scores_regions_array_sorted(2) > significance_level) && ((I_maxima_z_scores_regions_array_sorted(2) - I_maxima_z_scores_regions_array_sorted(2)) > z_score_threshold_difference)
        for j = 1:number_of_possible_regions
            if I_maxima_z_scores_regions_array_sorted(1) == I_maxima_z_scores_regions_array_abs(j)
                ectopic_pacemaker_location_indices_maxima(i) = j;
            end
            %This for loop determines the region number
        end
        if ectopic_pacemaker_location_indices_num_zeros(i) == ectopic_pacemaker_location_indices_maxima(i)
            location_counter_region1(i) = location_counter_region1(i) + 1;
        else
            location_counter_region2(i) = location_counter_region2(i) + 1;
            location_different_regions_counter(i) = location_different_regions_counter(i) + 1;
        end
    else
        ectopic_pacemaker_location_indices_maxima(i) = 0;
        %this means that the first metric could not be used
    end
    if (I_minima_z_scores_regions_array_sorted(1) < significance_level) && (I_minima_z_scores_regions_array_sorted(2) > significance_level) && ((I_minima_z_scores_regions_array_sorted(2) - I_minima_z_scores_regions_array_sorted(2)) > z_score_threshold_difference)
        for j = 1:number_of_possible_regions
            if I_minima_z_scores_regions_array_sorted(1) == I_minima_z_scores_regions_array_abs(j)
                ectopic_pacemaker_location_indices_minima(i) = j;
            end
            %This for loop determines the region number
        end
        if ectopic_pacemaker_location_indices_num_zeros(i) == ectopic_pacemaker_location_indices_minima(i)
            location_counter_region1(i) = location_counter_region1(i) + 1;
        else
            if ectopic_pacemaker_location_indices_maxima(i) == ectopic_pacemaker_location_indices_minima(i)
                location_counter_region2(i) = location_counter_region2(i) + 1;
            else
                
                location_counter_region3(i) = location_counter_region3(i) + 1;
                location_different_regions_counter(i) = location_different_regions_counter(i) + 1;
            end
        end
    else
        ectopic_pacemaker_location_indices_minima(i) = 0;
        %this means that the first metric could not be used
    end
    if (I_sym_inds_z_scores_regions_array_sorted(1) < significance_level) && (I_sym_inds_z_scores_regions_array_sorted(2) > significance_level) && ((I_sym_inds_z_scores_regions_array_sorted(2) - I_sym_inds_z_scores_regions_array_sorted(2)) > z_score_threshold_difference)
        for j = 1:number_of_possible_regions
            if I_sym_inds_z_scores_regions_array_sorted(1) == I_sym_inds_z_scores_regions_array_abs(j)
                ectopic_pacemaker_location_indices_sym_inds(i) = j;
            end
            %This for loop determines the region number
        end
        if ectopic_pacemaker_location_indices_num_zeros(i) == ectopic_pacemaker_location_indices_sym_inds(i)
            location_counter_region1(i) = location_counter_region1(i) + 1;
        else
            if ectopic_pacemaker_location_indices_maxima(i) == ectopic_pacemaker_location_indices_sym_inds(i)
                location_counter_region2(i) = location_counter_region2(i) + 1;
            else
                if ectopic_pacemaker_location_indices_minima(i) == ectopic_pacemaker_location_indices_sym_inds(i)
                    location_counter_region3(i) = location_counter_region3(i) + 1;
                else
                    location_counter_region4(i) = location_counter_region4(i) + 1;
                    location_different_regions_counter(i) = location_different_regions_counter(i) + 1;
                end
            end
        end
    else
        ectopic_pacemaker_location_indices_sym_inds(i) = 0;
        %this means that the first metric could not be used
    end
    if (I_concavs_z_scores_regions_array_sorted(1) < significance_level) && (I_concavs_z_scores_regions_array_sorted(2) > significance_level) && ((I_concavs_z_scores_regions_array_sorted(2) - I_concavs_z_scores_regions_array_sorted(2)) > z_score_threshold_difference)
        for j = 1:number_of_possible_regions
            if I_concavs_z_scores_regions_array_sorted(1) == I_concavs_z_scores_regions_array_abs(j)
                ectopic_pacemaker_location_indices_concavs(i) = j;
            end
            %This for loop determines the region number
        end
        if ectopic_pacemaker_location_indices_num_zeros(i) == ectopic_pacemaker_location_indices_concavs(i)
            location_counter_region1(i) = location_counter_region1(i) + 1;
        else
            if ectopic_pacemaker_location_indices_maxima(i) == ectopic_pacemaker_location_indices_concavs(i)
                location_counter_region2(i) = location_counter_region2(i) + 1;
            else
                if ectopic_pacemaker_location_indices_minima(i) == ectopic_pacemaker_location_indices_concavs(i)
                    location_counter_region3(i) = location_counter_region3(i) + 1;
                else
                    if ectopic_pacemaker_location_indices_sym_inds(i) == ectopic_pacemaker_location_indices_concavs(i)
                        location_counter_region4(i) = location_counter_region4(i) + 1;
                    else
                        location_different_regions_counter(i) = location_different_regions_counter(i) + 1;
                        location_counter_region5(i) = location_counter_region5(i) + 1;
                    end
                end
            end
        end
    else
        ectopic_pacemaker_location_indices_concavs(i) = 0;
        %this means that the first metric could not be used
    end
    if (I_mdslpes_z_scores_regions_array_sorted(1) < significance_level) && (I_mdslpes_z_scores_regions_array_sorted(2) > significance_level) && ((I_mdslpes_z_scores_regions_array_sorted(2) - I_mdslpes_z_scores_regions_array_sorted(2)) > z_score_threshold_difference)
        for j = 1:number_of_possible_regions
            if I_mdslpes_z_scores_regions_array_sorted(1) == I_mdslpes_z_scores_regions_array_abs(j)
                ectopic_pacemaker_location_indices_mdslpes(i) = j;
            end
            %This for loop determines the region number
        end
        if ectopic_pacemaker_location_indices_num_zeros(i) == ectopic_pacemaker_location_indices_mdslpes(i)
            location_counter_region1(i) = location_counter_region1(i) + 1;
        else
            if ectopic_pacemaker_location_indices_maxima(i) == ectopic_pacemaker_location_indices_mdslpes(i)
                location_counter_region2(i) = location_counter_region2(i) + 1;
            else
                if ectopic_pacemaker_location_indices_minima(i) == ectopic_pacemaker_location_indices_mdslpes(i)
                    location_counter_region3(i) = location_counter_region3(i) + 1;
                else
                    if ectopic_pacemaker_location_indices_sym_inds(i) == ectopic_pacemaker_location_indices_mdslpes(i)
                        location_counter_region4(i) = location_counter_region4(i) + 1;
                    else
                        if ectopic_pacemaker_location_indices_concavs(i) == ectopic_pacemaker_location_indices_mdslpes(i)
                            location_counter_region5(i) = location_counter_region5(i) + 1;
                        else
                            location_counter_region6(i) = location_counter_region6(i) + 1;
                            location_different_regions_counter(i) = location_different_regions_counter(i) + 1;
                        end
                    end
                end
            end
        end
    else
        ectopic_pacemaker_location_indices_mdslpes(i) = 0;
        %this means that the first metric could not be used
    end
    if (I_p_wave_areas_1_z_scores_regions_array_sorted(1) < significance_level) && (I_p_wave_areas_1_z_scores_regions_array_sorted(2) > significance_level) && ((I_p_wave_areas_1_z_scores_regions_array_sorted(2) - I_p_wave_areas_1_z_scores_regions_array_sorted(2)) > z_score_threshold_difference)
        for j = 1:number_of_possible_regions
            if I_p_wave_areas_1_z_scores_regions_array_sorted(1) == I_p_wave_areas_1_z_scores_regions_array_abs(j)
                ectopic_pacemaker_location_indices_p_wave_areas_1(i) = j;
            end
            %This for loop determines the region number
        end
        if ectopic_pacemaker_location_indices_num_zeros(i) == ectopic_pacemaker_location_indices_p_wave_areas_1(i)
            location_counter_region1(i) = location_counter_region1(i) + 1;
        else
            if ectopic_pacemaker_location_indices_maxima(i) == ectopic_pacemaker_location_indices_p_wave_areas_1(i)
                location_counter_region2(i) = location_counter_region2(i) + 1;
            else
                if ectopic_pacemaker_location_indices_minima(i) == ectopic_pacemaker_location_indices_p_wave_areas_1(i)
                    location_counter_region3(i) = location_counter_region3(i) + 1;
                else
                    if ectopic_pacemaker_location_indices_sym_inds(i) == ectopic_pacemaker_location_indices_p_wave_areas_1(i)
                        location_counter_region4(i) = location_counter_region4(i) + 1;
                    else
                        if ectopic_pacemaker_location_indices_concavs(i) == ectopic_pacemaker_location_indices_p_wave_areas_1(i)
                            location_counter_region5(i) = location_counter_region5(i) + 1;
                        else
                            if ectopic_pacemaker_location_indices_mdslpes(i) == ectopic_pacemaker_location_indices_p_wave_areas_1(i)
                                location_counter_region6(i) = location_counter_region6(i) + 1;
                            else
                                location_counter_region7(i) = location_counter_region7(i) + 1;
                                location_different_regions_counter(i) = location_different_regions_counter(i) + 1;
                            end
                        end
                    end
                end
            end
        end
    else
        ectopic_pacemaker_location_indices_p_wave_areas_1(i) = 0;
        %this means that the first metric could not be used
    end
    if (I_p_wave_areas_2_z_scores_regions_array_sorted(1) < significance_level) && (I_p_wave_areas_2_z_scores_regions_array_sorted(2) > significance_level) && ((I_p_wave_areas_2_z_scores_regions_array_sorted(2) - I_p_wave_areas_2_z_scores_regions_array_sorted(2)) > z_score_threshold_difference)
        for j = 1:number_of_possible_regions
            if I_p_wave_areas_2_z_scores_regions_array_sorted(1) == I_p_wave_areas_2_z_scores_regions_array_abs(j)
                ectopic_pacemaker_location_indices_p_wave_areas_2(i) = j;
            end
            %This for loop determines the region number
        end
        if ectopic_pacemaker_location_indices_num_zeros(i) == ectopic_pacemaker_location_indices_p_wave_areas_2(i)
            location_counter_region1(i) = location_counter_region1(i) + 1;
        else
            if ectopic_pacemaker_location_indices_maxima(i) == ectopic_pacemaker_location_indices_p_wave_areas_2(i)
                location_counter_region2(i) = location_counter_region2(i) + 1;
            else
                if ectopic_pacemaker_location_indices_minima(i) == ectopic_pacemaker_location_indices_p_wave_areas_2(i)
                    location_counter_region3(i) = location_counter_region3(i) + 1;
                else
                    if ectopic_pacemaker_location_indices_sym_inds(i) == ectopic_pacemaker_location_indices_p_wave_areas_2(i)
                        location_counter_region4(i) = location_counter_region4(i) + 1;
                    else
                        if ectopic_pacemaker_location_indices_concavs(i) == ectopic_pacemaker_location_indices_p_wave_areas_2(i)
                            location_counter_region5(i) = location_counter_region5(i) + 1;
                        else
                            if ectopic_pacemaker_location_indices_mdslpes(i) == ectopic_pacemaker_location_indices_p_wave_areas_2(i)
                                location_counter_region6(i) = location_counter_region6(i) + 1;
                            else
                                if ectopic_pacemaker_location_indices_p_wave_areas_1(i) == ectopic_pacemaker_location_indices_p_wave_areas_2(i)
                                    location_counter_region7(i) = location_counter_region7(i) + 1;
                                else
                                    location_counter_region8(i) = location_counter_region8(i) + 1;
                                    location_different_regions_counter(i) = location_different_regions_counter(i) + 1;
                                end
                            end
                        end
                    end
                end
            end
        end
    else
        ectopic_pacemaker_location_indices_p_wave_areas_2(i) = 0;
        %this means that the first metric could not be used
    end
    if (I_p_wave_areas_4_z_scores_regions_array_sorted(1) < significance_level) && (I_p_wave_areas_4_z_scores_regions_array_sorted(2) > significance_level) && ((I_p_wave_areas_4_z_scores_regions_array_sorted(2) - I_p_wave_areas_4_z_scores_regions_array_sorted(2)) > z_score_threshold_difference)
        for j = 1:number_of_possible_regions
            if I_p_wave_areas_4_z_scores_regions_array_sorted(1) == I_p_wave_areas_4_z_scores_regions_array_abs(j)
                ectopic_pacemaker_location_indices_p_wave_areas_4(i) = j;
            end
            %This for loop determines the region number
        end
        if ectopic_pacemaker_location_indices_num_zeros(i) == ectopic_pacemaker_location_indices_p_wave_areas_3(i)
            location_counter_region1(i) = location_counter_region1(i) + 1;
        else
            if ectopic_pacemaker_location_indices_maxima(i) == ectopic_pacemaker_location_indices_p_wave_areas_3(i)
                location_counter_region2(i) = location_counter_region2(i) + 1;
            else
                if ectopic_pacemaker_location_indices_minima(i) == ectopic_pacemaker_location_indices_p_wave_areas_3(i)
                    location_counter_region3(i) = location_counter_region3(i) + 1;
                else
                    if ectopic_pacemaker_location_indices_sym_inds(i) == ectopic_pacemaker_location_indices_p_wave_areas_3(i)
                        location_counter_region4(i) = location_counter_region4(i) + 1;
                    else
                        if ectopic_pacemaker_location_indices_concavs(i) == ectopic_pacemaker_location_indices_p_wave_areas_3(i)
                            location_counter_region5(i) = location_counter_region5(i) + 1;
                        else
                            if ectopic_pacemaker_location_indices_mdslpes(i) == ectopic_pacemaker_location_indices_p_wave_areas_3(i)
                                location_counter_region6(i) = location_counter_region6(i) + 1;
                            else
                                if ectopic_pacemaker_location_indices_p_wave_areas_1(i) == ectopic_pacemaker_location_indices_p_wave_areas_3(i)
                                    location_counter_region7(i) = location_counter_region7(i) + 1;
                                else
                                    if ectopic_pacemaker_location_indices_p_wave_areas_2(i) == ectopic_pacemaker_location_indices_p_wave_areas_3(i)
                                        location_counter_region8(i) = location_counter_region8(i) + 1;
                                    else
                                        location_counter_region9(i) = location_counter_region9(i) + 1;
                                        location_different_regions_counter(i) = location_different_regions_counter(i) + 1;
                                    end
                                end
                            end
                        end
                    end
                end
            end
        end
    else
        ectopic_pacemaker_location_indices_p_wave_areas_3(i) = 0;
        %this means that the first metric could not be used
    end
    if (I_p_wave_areas_4_z_scores_regions_array_sorted(1) < significance_level) && (I_p_wave_areas_4_z_scores_regions_array_sorted(2) > significance_level) && ((I_p_wave_areas_4_z_scores_regions_array_sorted(2) - I_p_wave_areas_4_z_scores_regions_array_sorted(2)) > z_score_threshold_difference)
        for j = 1:number_of_possible_regions
            if I_p_wave_areas_4_z_scores_regions_array_sorted(1) == I_p_wave_areas_4_z_scores_regions_array_abs(j)
                ectopic_pacemaker_location_indices_p_wave_areas_3(i) = j;
            end
            %This for loop determines the region number
        end
        if ectopic_pacemaker_location_indices_num_zeros(i) == ectopic_pacemaker_location_indices_p_wave_areas_4(i)
            location_counter_region1(i) = location_counter_region1(i) + 1;
        else
            if ectopic_pacemaker_location_indices_maxima(i) == ectopic_pacemaker_location_indices_p_wave_areas_4(i)
                location_counter_region2(i) = location_counter_region2(i) + 1;
            else
                if ectopic_pacemaker_location_indices_minima(i) == ectopic_pacemaker_location_indices_p_wave_areas_4(i)
                    location_counter_region3(i) = location_counter_region3(i) + 1;
                else
                    if ectopic_pacemaker_location_indices_sym_inds(i) == ectopic_pacemaker_location_indices_p_wave_areas_4(i)
                        location_counter_region4(i) = location_counter_region4(i) + 1;
                    else
                        if ectopic_pacemaker_location_indices_concavs(i) == ectopic_pacemaker_location_indices_p_wave_areas_4(i)
                            location_counter_region5(i) = location_counter_region5(i) + 1;
                        else
                            if ectopic_pacemaker_location_indices_mdslpes(i) == ectopic_pacemaker_location_indices_p_wave_areas_4(i)
                                location_counter_region6(i) = location_counter_region6(i) + 1;
                            else
                                if ectopic_pacemaker_location_indices_p_wave_areas_1(i) == ectopic_pacemaker_location_indices_p_wave_areas_4(i)
                                    location_counter_region7(i) = location_counter_region7(i) + 1;
                                else
                                    if ectopic_pacemaker_location_indices_p_wave_areas_2(i) == ectopic_pacemaker_location_indices_p_wave_areas_4(i)
                                        location_counter_region8(i) = location_counter_region8(i) + 1;
                                    else
                                        if ectopic_pacemaker_location_indices_p_wave_areas_3(i) == ectopic_pacemaker_location_indices_p_wave_areas_4(i)
                                            location_counter_region9(i) = location_counter_region9(i) + 1;
                                        else
                                            location_counter_region10(i) = location_counter_region10(i) + 1;
                                            location_different_regions_counter(i) = location_different_regions_counter(i) + 1;
                                        end
                                    end
                                end
                            end
                        end
                    end
                end
            end
        end
    else
        ectopic_pacemaker_location_indices_p_wave_areas_4(i) = 0;
        %this means that the first metric could not be used
    end
    if (I_p_wave_areas_5_z_scores_regions_array_sorted(1) < significance_level) && (I_p_wave_areas_5_z_scores_regions_array_sorted(2) > significance_level) && ((I_p_wave_areas_5_z_scores_regions_array_sorted(2) - I_p_wave_areas_5_z_scores_regions_array_sorted(2)) > z_score_threshold_difference)
        for j = 1:number_of_possible_regions
            if I_p_wave_areas_5_z_scores_regions_array_sorted(1) == I_p_wave_areas_5_z_scores_regions_array_abs(j)
                ectopic_pacemaker_location_indices_p_wave_areas_3(i) = j;
            end
            %This for loop determines the region number
        end
        if ectopic_pacemaker_location_indices_num_zeros(i) == ectopic_pacemaker_location_indices_p_wave_areas_5(i)
            location_counter_region1(i) = location_counter_region1(i) + 1;
        else
            if ectopic_pacemaker_location_indices_maxima(i) == ectopic_pacemaker_location_indices_p_wave_areas_5(i)
                location_counter_region2(i) = location_counter_region2(i) + 1;
            else
                if ectopic_pacemaker_location_indices_minima(i) == ectopic_pacemaker_location_indices_p_wave_areas_5(i)
                    location_counter_region3(i) = location_counter_region3(i) + 1;
                else
                    if ectopic_pacemaker_location_indices_sym_inds(i) == ectopic_pacemaker_location_indices_p_wave_areas_5(i)
                        location_counter_region4(i) = location_counter_region4(i) + 1;
                    else
                        if ectopic_pacemaker_location_indices_concavs(i) == ectopic_pacemaker_location_indices_p_wave_areas_5(i)
                            location_counter_region5(i) = location_counter_region5(i) + 1;
                        else
                            if ectopic_pacemaker_location_indices_mdslpes(i) == ectopic_pacemaker_location_indices_p_wave_areas_5(i)
                                location_counter_region6(i) = location_counter_region6(i) + 1;
                            else
                                if ectopic_pacemaker_location_indices_p_wave_areas_1(i) == ectopic_pacemaker_location_indices_p_wave_areas_5(i)
                                    location_counter_region7(i) = location_counter_region7(i) + 1;
                                else
                                    if ectopic_pacemaker_location_indices_p_wave_areas_2(i) == ectopic_pacemaker_location_indices_p_wave_areas_5(i)
                                        location_counter_region8(i) = location_counter_region8(i) + 1;
                                    else
                                        if ectopic_pacemaker_location_indices_p_wave_areas_3(i) == ectopic_pacemaker_location_indices_p_wave_areas_5(i)
                                            location_counter_region9(i) = location_counter_region9(i) + 1;
                                        else
                                            if ectopic_pacemaker_location_indices_p_wave_areas_4(i) == ectopic_pacemaker_location_indices_p_wave_areas_5(i)
                                                location_counter_region10(i) = location_counter_region10(i) + 1;
                                            else
                                                location_counter_region11(i) = location_counter_region11(i) + 1;
                                                location_different_regions_counter(i) = location_different_regions_counter(i) + 1;
                                            end
                                        end
                                    end
                                end
                            end
                        end
                    end
                end
            end
        end
    else
        ectopic_pacemaker_location_indices_p_wave_areas_5(i) = 0;
        %this means that the first metric could not be used
    end
end
ectopic_pacemaker_locations_indices = zeros(number_of_p_waves, number_of_parameters);
ectopic_pacemaker_locations_counters = zeros(number_of_p_waves, number_of_regions);
for i = 1:number_of_p_waves
    ectopic_pacemaker_locations_indices(i,1) = ectopic_pacemaker_location_indices_num_zeros(i);
    ectopic_pacemaker_locations_indices(i,2) = ectopic_pacemaker_location_indices_maxima(i);
    ectopic_pacemaker_locations_indices(i,3) = ectopic_pacemaker_location_indices_minima(i);
    ectopic_pacemaker_locations_indices(i,4) = ectopic_pacemaker_location_indices_sym_inds(i);
    ectopic_pacemaker_locations_indices(i,5) = ectopic_pacemaker_location_indices_concavs(i);
    ectopic_pacemaker_locations_indices(i,6) = ectopic_pacemaker_location_indices_mdlpes(i);
    ectopic_pacemaker_locations_indices(i,7) = ectopic_pacemaker_location_indices_p_wave_areas_1(i);
    ectopic_pacemaker_locations_indices(i,8) = ectopic_pacemaker_location_indices_p_wave_areas_2(i);
    ectopic_pacemaker_locations_indices(i,9) = ectopic_pacemaker_location_indices_p_wave_areas_3(i);
    ectopic_pacemaker_locations_indices(i,10) = ectopic_pacemaker_location_indices_p_wave_areas_4(i);
    ectopic_pacemaker_locations_indices(i,11) = ectopic_pacemaker_location_indices_p_wave_areas_5(i);
    ectopic_pacemaker_locations_counters(i,1) = location_counter_region1(i);
    ectopic_pacemaker_locations_counters(i,2) = location_counter_region2(i);
    ectopic_pacemaker_locations_counters(i,3) = location_counter_region3(i);
    ectopic_pacemaker_locations_counters(i,4) = location_counter_region4(i);
    ectopic_pacemaker_locations_counters(i,5) = location_counter_region5(i);
    ectopic_pacemaker_locations_counters(i,6) = location_counter_region6(i);
    ectopic_pacemaker_locations_counters(i,7) = location_counter_region7(i);
    ectopic_pacemaker_locations_counters(i,8) = location_counter_region8(i);
    ectopic_pacemaker_locations_counters(i,9) = location_counter_region9(i);
    ectopic_pacemaker_locations_counters(i,10) = location_counter_region10(i);
    ectopic_pacemaker_locations_counters(i,11) = location_counter_region11(i);
end
ectopic_pacemaker_location_indices_overall_determination = zeros(number_of_p_waves,1);
regions_counter_quotient_threshold = 0.7;
different_regions_quotient_threshold = 0.3;
for i = 1:number_of_p_waves
        %things you'll need
        ectopic_pacemaker_locations_counters_row = ectopic_pacemaker_locations_counters(i,:);
        ectopic_pacemaker_locations_indices_row = ectopic_pacemaker_locations_indices(i,:);
        ectopic_pacemaker_locations_counters_row_sorted = sort(ectopic_pacemaker_locations_counters_row);
        ectopic_pacemaker_locations_indices_row_sorted = sort(ectopic_pacemaker_locations_indices_row);
        most_common_next_most_common_region_diff = ectopic_pacemaker_locations_counters_row_sorted(end) - ectopic_pacemaker_locations_counters_row_sorted(end-1);
        most_common_next_most_common_region_tot_regions_quotient = most_common_next_most_common_region_diff/sum(ectopic_pacemaker_locations_counters_row);
        different_regions_quotient = location_different_regions_counter(i)/number_of_regions;
        condition1_4_region_determination = ectopic_pacemaker_locations_counters(i,j) == max(ectopic_pacemaker_locations_counters(i,:));
        condition2_4_region_determination = most_common_next_most_common_region_tot_regions_quotient > regions_counter_quotient_threshold;
        condition3_4_region_determination = different_regions_quotient < different_regions_quotient;
    for j = 1:number_of_parameters
        if condition1_4_region_determination && condition2_4_region_determination && condition3_4_region_determination
            ectopic_pacemaker_location_indices_overall_determination(i) = j;
        end
    end
end
for i = 1:length(ectopic_pacemaker_location_indices_overall_determination)
    ectopic_pacemaker_location_indices_overall_determination(i) = ectopic_pacemaker_location_indices_overall_determination(i) + 1;
end
regions = ['Cannot locate region  ','SA Node               '; 'Crista Terminalis     '; 'Tricuspid Annulus     '; 'Coronary Sinus        '; 'Ostium                '; 'Perinodal             '; 'Right Septum          '; 'Right Atrial Appendage'; 'Pulmonary Veins       '; 'Mitral Annulus        '; 'CS Body               '; 'Left Septum           '; 'Left Atrial Appendage '];
number_of_characters_req = 22;
%Above array has spaces in order to ensure that all region strings have the
%same number of elements
ectopic_pacemaker_region_locations = zeros(number_of_p_waves, number_of_characters_req);
for i = 1:number_of_p_waves
    for j =1:number_of_characters_req
        ectopic_pacemaker_region_locations(i,j) = 32;
    end
end
ectopic_pacemaker_region_locations = char(ectopic_pacemaker_region_locations);
for i = 1:number_of_p_waves
    ectopic_pacemaker_region_locations(i,:) = regions(ectopic_pacemaker_location_indices_overall_determination(i),:);
end

else
    
    disp('Signal is too noisy to make valid inferences about cardiac activity.')
    
end

