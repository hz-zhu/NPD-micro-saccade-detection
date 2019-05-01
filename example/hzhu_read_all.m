function data = hzhu_read_all(path,name)
if isempty(path)
    path = pwd;
end

if isempty(name)
    name = '.csv';
end

listing = dir([path,'\']);
counter = 0;
data = [];
for i = 1:length(listing)
    if contains(listing(i).name,name)
        local_name = listing(i).name;
        counter = counter+1;
        eval(['data.',local_name(1:end-4),'=csvread(''',path,'\',local_name,''');']);
        disp(['File ',listing(i).name,' loaded'])
    end
    data.n = counter;
end
end

