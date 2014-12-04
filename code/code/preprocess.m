% Preprocess the data by down-sampling

num_preictal_segments = 42;
num_interictal_segments = 500;

sampling_rate = 60;

base_file_pre = './Dog_2/Dog_2_preictal_segment_00';
base_file_inter = './Dog_2/Dog_2_interictal_segment_0';


for i = 1:num_preictal_segments
    i
    if i < 10
        data = load(strcat(base_file_pre,'0',int2str(i)));
    else
        data = load(strcat(base_file_pre, int2str(i)));
    end
    data = getfield(data, strcat('preictal_segment_', int2str(i)));
    data = data.data(:, 1:sampling_rate:end);
    
    save(strcat('./processed_dog_2/pre_', int2str(i)), 'data');    
     
end

for i = 1:num_interictal_segments
    i
    if i < 10
        data = load(strcat(base_file_inter,'00',int2str(i)));
    elseif i < 100
        data = load(strcat(base_file_inter,'0',int2str(i)));
    else
        data = load(strcat(base_file_inter, int2str(i)));
    end
    data = getfield(data, strcat('interictal_segment_', int2str(i)));
    data = data.data(:, 1:sampling_rate:end);
    
    save(strcat('./processed_dog_2/inter_', int2str(i)), 'data');    
     
end