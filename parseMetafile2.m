function [phd, ic1, eng, expt, metaInfo]=parseMetafile2(fileName, modeN)
% [phd, ic1, eng, expt, metaInfo]=parseMetafile(fileName, modeN)
% phd: It value
% ic1: Io value
% eng: x-ray energy
% expt: Exposure time
% 
%fid=fopen(fileName);
metaInfo.expttime='';
metaInfo.datadir='';
phd=1;
eng = -1;
ic1 = -1;
expt = -1.;
exptset = -1.;

fid=fopen(fileName);
if fid<0
    sprintf('%s not openned', fileName);
    return;
end

if (nargin == 1)
    modeN = -1; % photodiode number is correctly handled
end

if(fid ~= -1)
try
    while feof(fid)==0
        line=fgetl(fid);
        if (length(line) > 8 & line(1:7) =='% Date:')
                metaInfo.expttime= strtrim(line(9:end));
                continue;
        end    
        if (length(line) > 16 & line(1:13) =='% Data Direct')
                metaInfo.datadir= strtrim(line(20:end));
                continue;
        end           
        if (length(line) > 8 & line(1:7) =='% IC1 :')
                ic1= str2num(strtrim(line(9:end)));
                if ic1 == 0
                    ic1 =1.0;
                end    
                continue;
        end    

        if (length(line) > 8 & line(1:6) =='% I0 :')
                ic1= str2num(strtrim(line(8:end)));
                if ic1 == 0
                    ic1 =1.0;
                end    
                continue;
        end            
        
        switch modeN
            case -1
                if length(line) > 14 & lower(line(1:14)) =='% photodiode :'
                        phd= str2num(strtrim(line(15:end)));
                        if phd==0.
                            phd = 1.0;
                        end    
                        continue;
                end    
                
%                 if length(line) > 10 & line(1:10) =='% SAXSBS :'
%                         phd= str2num(strtrim(line(11:end)));
%                         if phd==0.
%                             phd = 1.0;
%                         end  
%                         continue;
%                 end
            case 0   % SAXS
                if length(line) > 10 & line(1:10) =='% SAXSBS :'
                        phd= str2num(strtrim(line(11:end)));
                        if phd==0.
                            phd = 1.0;
                        end  
                        continue;
                end
                
%                 if length(line) > 14 & lower(line(1:14)) =='% photodiode :'
%                         phd= str2num(strtrim(line(15:end)));
%                         if phd==0.
%                             phd = 1.0;
%                         end    
%                         continue;
%                 end    
                
            case 1  % GISAXS
                if length(line) > 12 & line(1:12) =='% GISAXSBS :'
                    phd= str2num(strtrim(line(13:end)));
                    if phd==0.
                        phd = 1.0;
                    end
                    continue;
                end
            case 2  % GISAXS
                if length(line) > 12 & line(1:12) =='% GIWAXSBS :'
                    phd= str2num(strtrim(line(13:end)));
                    if phd==0.
                        phd = 1.0;
                    end
                    continue;
                end
            case 3  % Central BS
                if length(line) > 12 & line(1:12) =='% CenterBS :'
                    phd= str2num(strtrim(line(13:end)));
                    if phd==0.
                        phd = 1.0;
                    end
                    continue;
                end
        end

        if length(line) > 17 & line(1:14) =='% X-ray Energy'
                %eng= str2num(strtrim(line(17:end)));
               [~, engStr]= strtok(line,':');
                eng = str2num(engStr(2:end));
                continue;
        end    
        
        if length(line) > 20 & lower(line(1:16)) =='% exposure time '
            [~, exptStr]= strtok(line,':');
            expt = str2num(exptStr(2:end));
            %expt= str2num(strtrim(line(18:end)));   
            continue;
        end    
        
        if length(line) > 23 & lower(line(1:20)) =='% set exposure time '
            [~, exptSetStr]= strtok(line,':');
            exptset = str2num(exptSetStr(2:end));    
            %exptset= str2num(strtrim(line(22:end)));
            continue;
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
