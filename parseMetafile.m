function [phd, ic1, eng, expt]=parseMetafile(fileName, modeN)
% [phd, ic1, eng, expt]=parseMetafile(fileName, modeN)
phd=-1;
ic1=-1;
eng=-1;
expt=-1;
fid=fopen(fileName);
expt = -1.;
exptset = -1.;
if (nargin == 1)
    modeN = -1; % photodiode number is correctly handled
end

if(fid ~= -1)
try
    while feof(fid)==0
        line=fgetl(fid);
        if (length(line) > 8 & line(1:7) =='% IC1 :')
                ic1= str2num(strtrim(line(9:end)));
                if ic1 == 0
                    ic1 =1.0;
                end    
                continue;
        end            
%         if (length(line) > 8 & line(1:7) =='% IC1 :')
%                 ic1= str2num(strtrim(line(9:end)));
%                 if ic1 == 0
%                     ic1 =1.0;
%                 end    
%         end    
        
        switch modeN
            case -1
                if length(line) > 14 & lower(line(1:14)) =='% photodiode :'
                        phd= str2num(strtrim(line(15:end)));
                        if phd==0.
                            phd = 1.0;
                        end    
                        continue;
                end                    
%                 if length(line) > 14 & line(1:14) =='% photodiode :'
%                         phd= str2num(strtrim(line(15:end)));
%                         if phd==0.
%                             phd = 1.0;
%                         end    
%                 end    
            case 0   % SAXS
                if length(line) > 10 & line(1:10) =='% SAXSBS :'
                        phd= str2num(strtrim(line(11:end)));
                        if phd==0.
                            phd = 1.0;
                        end  
                end
                
            case 1  % GISAXS
                if length(line) > 12 & line(1:12) =='% GISAXSBS :'
                    phd= str2num(strtrim(line(13:end)));
                    if phd==0.
                        phd = 1.0;
                    end
                end
            case 2  % GISAXS
                if length(line) > 12 & line(1:12) =='% GIWAXSBS :'
                    phd= str2num(strtrim(line(13:end)));
                    if phd==0.
                        phd = 1.0;
                    end
                end
            case 3  % GISAXS
                if length(line) > 12 & line(1:12) =='% CenterBS :'
                    phd= str2num(strtrim(line(13:end)));
                    if phd==0.
                        phd = 1.0;
                    end
                end
        end

        if length(line) > 17 & line(1:16) =='% X-ray Energy :'
                eng= str2num(strtrim(line(17:end)));
        end    
        
        if length(line) > 20 & line(1:17) =='% Exposure time :'
                expt= str2num(strtrim(line(18:end)));                
        end    
        
        if length(line) > 23 & line(1:21) =='% Set Exposure time :'
                exptset= str2num(strtrim(line(22:end)));
        end           
    end 

catch
    phd = 1;
    ic1 = 1;
    eng = 14;
end	

    fclose(fid);
    
    if expt >0 & exptset >0
        phd = phd * exptset / expt;
        ic1 = ic1 * exptset / expt;
        %[num2str(exptset) ', ' num2str(expt)] 
    end
else
    phd = 1.;
    ic1 = 1.;
    eng = -1.;
end  
