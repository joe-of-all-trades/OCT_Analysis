% This script is written by Chao-Yuan Yeh. All copyrights reserved. 
function OCT_process_file_pub
root_folder = uigetdir('', 'select root folder containing all OCT data');
recur_proc(root_folder)
% Consider alternative implementation using built-in genpath function
matlab_path = strsplit(matlabpath, ';');
home_path = matlab_path{1};
cd(home_path)
disp('task completed')
end

function recur_proc(foldername)
% function that goes through folder structure resursively and process OCT
% files along the way

cd(foldername)
folderlist = dir;

sublist = folderlist([folderlist(:).isdir] &...
                      ~strcmp({folderlist(:).name}, '.') &...
                      ~strcmp({folderlist(:).name}, '..'));

if ~exist('filename_proccesed.log', 'file')
    fID = fopen('filename_proccesed.log','w');
    fclose(fID);
    
    filelist = dir('*.OCT');
    % Shorten file name(removing ID number) so filename can be handled by
    % downstream process. 
    for jj = 1 : length(filelist)
        movefile(filelist(jj).name,regexprep(filelist(jj).name,...
          ',\s[A-Z]\d{9}\s_\d*_','_'));
    end
    clear filelist
    
    % Convert Cross Line OCT files into tif
    filelist = dir('*Cross Line*.OCT');
    for jj = 1:length(filelist)
        tempstr = filelist(jj).name;
        fID = fopen(fullfile(pwd,tempstr));
        temp = fread(fID,'single');
        % Image size is specified by "OCT Window Height=768" and "XY Scan Length= 1019"
        temp = uint16(reshape(temp,768,1019,4)); 
        img1 = flip(temp(:,:,1),1);
        img2 = flip(temp(:,:,2),1);
        imwrite(img1,[tempstr(1:end-4),'_1.tif']);
        imwrite(img2,[tempstr(1:end-4),'_2.tif']);
        fclose(fID);
        clearvars tempstr fID temp img1 img2
    end
    clear filelist
    
    % Convert Line OCT files into tif
    filelist = dir('*_Line*.OCT');
    for jj = 1:length(filelist)
        tempstr = filelist(jj).name;
        fID = fopen(fullfile(pwd,tempstr));
        temp = fread(fID,'single');
        % Image size is specified by "OCT Window Height=956" and "XY Scan Length= 1019"
        temp = uint16(reshape(temp,956,1019,2)); 
        img1 = flip(temp(:,:,1),1);
        imwrite(img1,[tempstr(1:end-4),'.tif']);
        fclose(fID);
        clearvars tempstr fID temp img1 
    end
    clear filelist
end
    
if ~isempty(sublist)
    for ii = 1:length(sublist)
        sublist(ii).abspath = fullfile(pwd,sublist(ii).name);
    end
end

if ~isempty(sublist)
    for ii = 1:length(sublist)
        recur_proc(sublist(ii).abspath)
    end
end
end