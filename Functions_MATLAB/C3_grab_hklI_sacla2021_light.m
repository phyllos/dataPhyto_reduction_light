function grab_hklI_sacla2021_light(binName)
% Run as: 
% grab_hklI_sacla2021_light('3ps')  % % the name of streamfile
% A. Hosseini, Feb-2022
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
  fidI = fopen([binName '.stream'],'r');
  directory = 'sacla2021_3ps_hklI';
  if ~exist(directory,'dir')
    system(['mkdir ' directory]);
  end
  fidO_snapshotInfo = fopen([binName '_snapshotInfo.dat'],'w');
  num_snapshot = 0;
  while ~feof(fidI)
    found_indexed_snapshot = false;
    while ~found_indexed_snapshot
      Image_filename_line = goto_tag(fidI,'Image filename:')
     %Event_line = goto_tag(fidI,'Event: //')
      Event_line = goto_tag(fidI,'Event: tag-')
      indexed_by_line = goto_tag(fidI,'indexed_by')
      if feof(fidI)
        fclose(fidO_snapshotInfo);
        fclose(fidI);
        return
      end
      found_indexed_snapshot = ~strcmp(indexed_by_line(14:17),'none')
    end
    eval([goto_tag(fidI,'num_reflections') ';'])
    goto_tag(fidI,'Reflections measured after indexing')
    fgetl(fidI);
    myData = nan(num_reflections,4);
    for jj=1:num_reflections
      A = sscanf(fgetl(fidI),'%4d%5d%5d%11f%11f%11f%11f%7f%7f%s');
      if isempty(A), break, end
%      myData(jj,:) = [asuP6_3(A(1:3)) A(4)];
       myData(jj,:) = [symP212121(A(1:3)) A(4)];
    end
    if any(isnan(myData(:)))|any(abs(myData(:))>1e10)|any(~isnumeric(myData(:)))
      continue
    end
    num_snapshot = num_snapshot+1;
    [runID,eventID] = extract_snapshot_id(Image_filename_line,Event_line);
%     runID = 1009987
%     eventID = 1498344208
%     num_reflections = 391
    fprintf(fidO_snapshotInfo,'%10d %16d %10d\n',runID,eventID,num_reflections);
%   fidO_hklI = fopen([directory '/cxig0915_' num2str(runID) '_' num2str(eventID) '_hklI.dat'],'w');
    fidO_hklI = fopen([directory '/sacla2021_run' num2str(runID) '_tag' num2str(eventID) '_hklI.dat'],'w');
    fprintf(fidO_hklI,'%4d%5d%5d%11.2f\n',myData');
    fclose(fidO_hklI);
    
  end

%%%%%%%%
function current_line=goto_tag(fid,tag)
  length_tag = length(tag);
  found = false;
  while (~found)&(~feof(fid))
    current_line = fgetl(fid);
    if (length(current_line)>length_tag-1)
      found = strcmp(current_line(1:length_tag),tag);
    end
  end
  if feof(fid)
    current_line = '';
  end
%EOFunc

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [runID,eventID]=extract_snapshot_id(Image_filename_line,Event_line)
  s0 = strfind(Image_filename_line,'/');
  s0 = s0(end);
  filename = Image_filename_line(s0+1:end);
  %s0 = strfind(filename,'-r');
  %s1 = strfind(filename,'.cxi');
  s0 = strfind(filename,'run')
  s1 = strfind(filename,'-light.h5');
  %runID = str2num(filename(s0+2:s1-1));
  runID = str2num(filename(s0+3:s1-1))

  s0 = strfind(Event_line,'tag-');
  %s0 = s0(end)
  s1 = strfind(Event_line,'//');
  %eventID = str2num(Event_line(s0+1:end))
  Event_line(s0:s1-1)
  eventID = str2num(Event_line(s0+4:s1-1));
%EOFunc

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function hkl_uwm = symP212121(hkl)
  h = hkl(1);
  k = hkl(2);
  l = hkl(3);

  if(h>=0 && k>=0 && l>=0)     %(h,k,l)
      hkl_uwm = [h,k,l];

  elseif(h>=0 && k<=0 && l<=0) %(h,-k,-l)
      hkl_uwm = [h,-k,-l];

  elseif(h<=0 && k>=0 && l<=0) %(-h,k,-l)
      hkl_uwm = [-h,k,-l];

  elseif(h<=0 && k<=0 && l>=0) %(-h,-k,l)
      hkl_uwm = [-h,-k,l];      
      
  elseif(h<=0 && k>=0 && l>=0) %(-h,k,l)
      hkl_uwm = [-h,k,l];

  elseif(h>=0 && k<=0 && l>=0) %(h,-k,l)
      hkl_uwm = [h,-k,l];

  elseif(h>=0 && k>=0 && l<=0) %(h,k,-l)
      hkl_uwm = [h,k,-l];

  elseif(h<0 && k<0 && l<0)    %(-h,-k,-l)
      hkl_uwm = [-h,-k,-l]; 
      
% All other cases    
  else
      hkl_uwm = [h,k,l];  
  end
%end function hkl_uwm
