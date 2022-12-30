
function SR=get_sample_rate(filename)
if exist(filename)~=2
    filename=strcat(['raw/' filename]);
end
if exist(filename)==2
fid = fopen(filename);
bytes_per_header = 2^14; %default for all Neuralynx files
status = fseek(fid,bytes_per_header+8+4,'bof');
SR = fread(fid,1,'uint32');
fclose(fid);
else
   SR=[]; 
end
end

