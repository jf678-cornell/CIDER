function [newdata,stdev] = averageEvery(sampling_length,start_index,data)
data_size = size(data);
timetotal = data_size(3);
newset = start_index:sampling_length:timetotal;
newdata_size = data_size;
if newset(end)+sampling_length-1>timetotal
    newset = newset(1:end-1);
end


newdata_size(3) = length(newset);



newdata = zeros(newdata_size);
stdev = zeros(newdata_size);

for i = 1:length(newset)
    newdata(:,:,i,:) = mean(data(:,:,newset(i):(newset(i)-1+sampling_length),:),3);
    stdev(:,:,i,:) = std(data(:,:,newset(i):(newset(i)-1+sampling_length),:),0,3);
end



end