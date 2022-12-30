function [data,timestamps,info]=read_neuralynx_file(filename, varargin)

if exist('raw')==7
    filename=strcat(['raw/' filename]);
end
data=[]; timestamps=[]; info=[];

if exist(filename)==2
    [~,~,filetype] = fileparts(filename);
    if ~any(strcmp(filetype,{'.ncs','.nvt','.ntt','.nev'}))
        error('File extension not recognized');
    end
   
      if ~isempty(varargin)
         parameters=varargin;
      else
         parameters={[1 1 1 1 1],1,1, []};
      end
     
    if strcmp(computer,'MACI64')  %if using a mac
        switch filetype
            case '.ncs'   %CSC file
                    [timestamps, ~, ~, data.NumberOfValidSamples,...
                    data.Samples, info.Header] = Nlx2MatCSC_v3(filename, parameters{:});
            case '.nvt'  %video file
                [timestamps, data.ExtractedX, data.ExtractedY, ~, data.Targets,  ~, info.Header] = Nlx2MatVT_v3(filename, [1 1 1 1 1 1],1, 1, []);

            case '.ntt'  %spike file
                [timestamps, ~, ~, data.Params, ~,info.header] = Nlx2MatSpike_v3(filename, [1 1 1 1 1],1,1, []);

            case '.nev'  %event string
                [timestamps, data.Event_ID, data.ttl, data.Extras, data.EventStrings, info.Header] = Nlx2MatEV_v3(filename,  [1 1 1 1 1],1,1, []);
        end
    else  %if using a PC
        switch filetype
            case '.ncs'   %CSC file
                [timestamps, ~, ~, data.NumberOfValidSamples,...
                    data.Samples, info.Header] = Nlx2MatCSC(filename,  parameters{:});
               
            case '.nvt'  %video file
                [timestamps, data.ExtractedX, data.ExtractedY, ~, data.Targets, ...
                    ~, info.Header] = Nlx2MatVT(filename,  [1 1 1 1 1 1],1, 1, []);

            case '.ntt'  %spike file

                [timestamps, ~, ~, data.Params, ~,info.header] = Nlx2MatSpike(filename,  [1 1 1 1 1],1,1, []);

            case '.nev'  %event string
                 [timestamps, data.Event_ID, data.ttl, data.Extras, data.EventStrings, info.Header]  = Nlx2MatEV(filename, [1 1 1 1 1],1,1, []);
               
        end
        
        
    end
    
else
    disp('file does not exist')    
end
