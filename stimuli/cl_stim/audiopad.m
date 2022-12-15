function MATpad

file = dir('*.mat');
Fs = 44100;
for i = 1:length(file)
    file(i).name
    load(file(i).name); 
    pad = zeros(1, round((5*Fs-length(y))/2));     
    y = [pad';y;pad'];
    save(file(i).name,'y');
end