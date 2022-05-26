%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%         Ruben Schulte         %%%%%
%%%%%   CARTESIUS DALES read data   %%%%%
%%%%%      started: 28-12-2020      %%%%%
%%%%%      changed: 07-04-2022      %%%%%
%%%%%       final: 23-05-2025       %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

close all
clear all
clc
direc = struct();
check = struct();


%% General settings
% These variables can be changed by the uses

% Time correction [hours] from UTC to local time (ECTS
check.corr_localtime = 2;      
% Define the relevant (data) directories 
direc.main_directory = '/scratch-shared/rschulte_temp';
% direc.main_directory = 'D:\DALES_data';
direc.run_name = 'nh3plume';
direc.exp_nr = '002'; 
% Difine directory separators for Linux & Windows
direc.dir_separator = '/';  % Linux (Snellius)
% direc.dir_separator = '\';  % Windows (laptop)

check.variable_background_indicator = 'b';
check.variable_plume_indicator = 'p';

%%%%%
% Additional options:
%%%%%

check.z_max = 2200;                     % Maximum z value [m] for cross-section output
check.center_sources = 'yes';           % Center source(s)?

% Combine (i.e. sum) variables  
check.combine = 'yes';              % Do variables need to be combined yes/no?
% Define variables to be combined (i.e. summed) 
    % sets are separated by ';', to be combined variables in a set separated by ','
check.combine_vars =    {'nh3_r1b', 'nh3_r1p' ...   % Reference set 1
                       }; 
% Define name of the new combined variable (separate sets by ';')
check.combine_newname = {'nh3r1' ...                                % Reference set
                         };      

% Define the surface.interactive.inp land_types which match emission for
% each nh3 scalar ({ "scalar name 1 ", [land_type nr 1]; "scalar name 2 ", [land_type nr 2]; ...} )
check.plume_receptor_source = {'nh3_r1b',   [1] ...
                             ; 'nh3_r1p',   [1] ...
                             ; 'nh3r1',     [1] ...
                               };
                           
% Process only specific cross-sections
check.process_cross = {'all'};          % Process all cross-sections
check.process_cross = {'xz' ...
                     ; 'xy02' ...
                     ; 'xy03' ...
                     ; 'xy04' ...
                     ; 'xy05' ...
                     ; 'xy06' ...
                     ; 'xy07' ...
                     ; 'xy08' ...
                     ; 'xy09' ...
                     ; 'xy10' ...
                     ; 'xy11' ...
                     ; 'xy12' ...
                     ; 'xy13' ...
                     ; 'xy14' ...
                     ; 'xy15' ...
                     ; 'xy16' ...
                     ; 'xy17' ...
                     ; 'xy18' ...
                     ; 'xy19' ...
                     ; 'xy20' ...
                     ; 'xy21' ...
                     ; 'xy22' ...
                     ; 'xy23' ...
                     ; 'xy24' ...
                      }; 
%                      ; 'xy02' ...
% Process only specific variables (options are 'w','thl','qt' and scalars)
check.process_var = {'all'};          % Process all variables
% check.process_var = {'w' ...
%                    ; 'thl' ...
%                    ; 'qt' ...
%                    ; 'nh3_r0b' ...
%                    ; 'nh3_r1b' ...
%                    ; 'nh3_r1p' ...
%                    ; 'nh3r1' ...
%                     }; 


%% End of settings



% !!!!!!!!!!!!!!!!!!!! Do not change the code below !!!!!!!!!!!!!!!!!!!!



%% Make sure the directories end with an '\' character

% Structure field names in direc structure
fields_direc = fields(direc);
% Loop through structure fields
for i = 1 : length(fields_direc) 
    fields_temp = direc.(char(fields_direc(i)));
    % If needed, add a '\' to the directory name string
    if strcmp(char(fields_temp(end)),direc.dir_separator) == 0
        direc.(char(fields_direc(i))) = strcat(fields_temp, direc.dir_separator);
    end
end

clear fields_* i


%% Start the timer

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Message: Data needs to be loaded
tload_start = datestr(now, 'dd-mm-yyyy HH:MM:SS');
fprintf(['f__DALESread: Loading the DALES output data for run ' direc.exp_nr(1:end-1) '\n\n' ... 
         '              --> Starting at: ' tload_start '\n\n']) 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Start loading the DALES output data


% Preallocate input structure
inp = struct();

%%%%%%%%%%
% Save the processed outputs in a separate folder
%%%%%%%%%%
name_list = dir([direc.main_directory direc.run_name direc.exp_nr ]);
name_list = {name_list.name}';
output_dir = 'output_v001';
dir_counter = 0;
% Determine the version of the output_old folder
while dir_counter == 0 
    if ismember({output_dir}, name_list) == 0
        dir_counter = 1;
    else
        olddir_name_nr = str2num(output_dir(end-2:end)) +1;
        if length(num2str(olddir_name_nr)) < 4
            output_dir = [output_dir( 1:end-length(num2str(olddir_name_nr)) ) num2str(olddir_name_nr)];
        end
    end
end
direc.output_dir = output_dir;
    
% Create the output folder
mkdir([ direc.main_directory direc.run_name direc.exp_nr ], direc.output_dir);

clear dir_counter name_list output_dir

%% Input files
% Relevant input files: 1) namoptions
%                       2) lscale.inp
%                       3) scalar.inp
%                       4) prof.inp

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf('              --> Loading relevant input files..\n') 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Get the names of all the files in the input directory
input_namelist = dir([direc.main_directory, direc.run_name, direc.exp_nr]);
input_namelist = {input_namelist.name}';

%%%%%%%%%%
%%%%% 1) namoptions
if ismember({['namoptions.' direc.exp_nr(1:end-1)]}, input_namelist)
    try
        % Define the location/name of the namoptions file
        namoptions_dir = [direc.main_directory, direc.run_name, direc.exp_nr 'namoptions.' direc.exp_nr(1:end-1)];
        % Call function to read the namoptions file and add to the inp
        % structure
        inp.namoptions = f__READnamoptions(namoptions_dir);
        %
        check.modelstart = str2num( inp.namoptions.value( ismember( inp.namoptions.option, {'xtime'} ) ) );
        % Clear temporary variables
        clear namoptions_dir
    catch
        fprintf(['              --> ERROR!\n' ...
                 '                  Something went wrong when reading namoptions.' direc.exp_nr(1:end-1) '.\n' ...
                 'f__DALESread: TERMINATED due to error\n'])
        return
    end
else
    fprintf(['              --> ERROR!\n' ...
             '                  No namoptions.' direc.exp_nr(1:end-1) ' exists in directory ' [direc.main_directory, direc.run_name, direc.exp_nr] '.\n' ...
             'f__DALESread: TERMINATED due to error\n'])
    return
end

%%%%%%%%%%
%%%%% 2) lscale.inp
if ismember({['lscale.inp.' direc.exp_nr(1:end-1)]}, input_namelist)
    try
        % Define the location/name of the lscale.inp file
        lscale_dir = [direc.main_directory, direc.run_name, direc.exp_nr 'lscale.inp.' direc.exp_nr(1:end-1)];
        % Call function to read the lscale.inp file and add to the inp
        % structure
        inp.lscale = f__READinput_textscan(lscale_dir, inp);
        % Clear temporary variables
        clear lscale_dir
    catch
        fprintf(['            --> ERROR!\n' ...
                 '                Something went wrong when reading lscale.inp.' direc.exp_nr(1:end-1) '.\n' ...
                 'f__DALESread: TERMINATED due to error\n'])
        return
    end
else
    fprintf(['              --> ERROR!\n' ...
             '                  No lscale.inp.' direc.exp_nr(1:end-1) ' exists in directory ' [direc.main_directory, direc.run_name, direc.exp_nr] '.\n' ...
             'f__DALESread: TERMINATED due to error\n'])
    return
end

%%%%%%%%%%
%%%%% 3) scalar.inp
if ismember({['scalar.inp.' direc.exp_nr(1:end-1)]}, input_namelist)
    try
        % Define the location/name of the scalar.inp file
        scalar_dir = [direc.main_directory, direc.run_name, direc.exp_nr 'scalar.inp.' direc.exp_nr(1:end-1)];
        % Call function to read the scalar.inp file and add to the inp
        % structure
        inp.scalar = f__READinput_textscan(scalar_dir, inp);
        % Clear temporary variables
        clear lscalar_dir
    catch
        fprintf(['              --> ERROR!\n' ...
                 '                  Something went wrong when reading scalar.inp.' direc.exp_nr(1:end-1) '.\n' ...
                 'f__DALESread: TERMINATED due to error\n'])
        return
    end
else
    fprintf(['              --> ERROR!\n' ...
             '                  No scalar.inp.' direc.exp_nr(1:end-1) ' exists in directory ' [direc.main_directory, direc.run_name, direc.exp_nr] '.\n' ...
             'f__DALESread: TERMINATED due to error\n'])
    return
end

%%%%%%%%%%
%%%%% 4) prof.inp
if ismember({['prof.inp.' direc.exp_nr(1:end-1)]}, input_namelist)
    try
        % Define the location/name of the scalar.inp file
        prof_dir = [direc.main_directory, direc.run_name, direc.exp_nr 'prof.inp.' direc.exp_nr(1:end-1)];
        % Call function to read the scalar.inp file and add to the inp
        % structure
        inp.prof = f__READinput_textscan(prof_dir, inp);
        % Clear temporary variables
        clear prof_dir
    catch
        fprintf(['              --> ERROR!\n' ...
                 '                  Something went wrong when reading prof.inp.' direc.exp_nr(1:end-1) '.\n' ...
                 'f__DALESread: TERMINATED due to error\n'])
        return
    end
else
    fprintf(['              --> ERROR!\n' ...
             '                  No prof.inp.' direc.exp_nr(1:end-1) ' exists in directory ' [direc.main_directory, direc.run_name, direc.exp_nr] '.\n' ...
             'f__DALESread: TERMINATED due to error\n'])
    return
end

% Save the inp data
save([direc.main_directory direc.run_name direc.exp_nr ...                  % Save directory
    direc.output_dir direc.dir_separator direc.exp_nr(1:end-1) 'inp.mat'], ...          % Save name
    '-struct', 'inp', '-v7.3');                                                        % Saved variables

clear scalar_dir input_namelist

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf('                  Relevant input files are loaded and saved!\n') 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% Define the vertical levels of the xy cross-secitons

if strcmp(inp.namoptions.value{ismember(inp.namoptions.option, {'lcross'})}, '.true.')
    if ismember({'crossheight'},inp.namoptions.option)
        cross_lvl_nr = str2num(inp.namoptions.value{ismember(inp.namoptions.option, {'crossheight'})});
        if cross_lvl_nr > 99
            cross_lvl_nr = [];
        end
    else
        cross_lvl_nr = 2;
    end
else
    cross_lvl_nr = [];
end

if strcmp(inp.namoptions.value{ismember(inp.namoptions.option, {'lfielddump'})}, '.true.')
    if ismember({'klow'},inp.namoptions.option)
        klow = str2num(inp.namoptions.value{ismember(inp.namoptions.option, {'klow'})});
    else
        klow = 1;
    end
    
    if ismember({'khigh'},inp.namoptions.option)
        khigh = str2num(inp.namoptions.value{ismember(inp.namoptions.option, {'khigh'})});
    else
        if ismember({'kmax'},inp.namoptions.option)
            khigh = str2num(inp.namoptions.value{ismember(inp.namoptions.option, {'kmax'})});
        else
            khigh = 96;
        end
    end
    if khigh > 100
        khigh = 99;
    end
    if ismember({'ncoarse'},inp.namoptions.option)
        ncoarse = str2num(inp.namoptions.value{ismember(inp.namoptions.option, {'ncoarse'})});
    else
        ncoarse = 1;
    end
    
    cross_lvl_nr = unique([cross_lvl_nr, klow : ncoarse : khigh]);
    
    clear ncoarse khigh klow
end

check.cross_lvl = cell(length(cross_lvl_nr),1);
for i = 1 : length(check.cross_lvl)
    if cross_lvl_nr(i) >= 10
        check.cross_lvl{i} = num2str(cross_lvl_nr(i));
    else
        check.cross_lvl{i} = ['0' num2str(cross_lvl_nr(i))];
    end
end

if isempty(check.z_max) || check.z_max > max(inp.prof.zc)
    check.z_max = max(inp.prof.zc);
end

clear cross_lvl_nr i

%% Define the cross-sections and variables that need to be processed

% See if all cross-sections need to be processed
if strcmp(check.process_cross{1}, 'all')
    
    check.process_cross = cell(length(check.cross_lvl)+2,1);
    check.process_cross{1} = 'xz';
    check.process_cross{2} = 'yz';
    for i = 1 : length(check.cross_lvl)
        check.process_cross{2+i} = ['xy' check.cross_lvl{i}];
    end
    
% Define the specific cross-sections that will be processed
else
    % Make sure the defined cross-sections actually exists and remove
    % incorrect ones
    for i = length(check.process_cross) : -1 : 1
        if strcmp(check.process_cross{i}(1:2), 'xy')
            if ismember({check.process_cross{i}(3:end)} , check.cross_lvl) == 0
                check.process_cross(i) = [];
            end
        elseif ismember(check.process_cross(i), {'xz','yz'}) == 0
            check.process_cross(i) = [];
        end
    end
end

% Define all available scalars
tmp_allvar = fields(inp.scalar);
tmp_allvar = tmp_allvar(ismember(tmp_allvar, {'zc','Properties','Row','Variables'}) == 0);
tmp_allvar = [tmp_allvar ; check.combine_newname;  {'w';'thl';'qt'} ];

% make sure that the predefined varialbes exist in the "all available variables"
if strcmp(check.process_var{1}, 'all')
    check.process_var = tmp_allvar;
else
    for i = length(check.process_var) :-1: 1
        if ismember(check.process_var(i), tmp_allvar) == 0
            check.process_var(i) = [];
        end
    end
end

clear i tmp_allvar


%% Load the genstat output
    % NOTE: This dataset is mandatory
    
% Check if the NAMGENSTAT output is .true. in namoptions 
if strcmp( inp.namoptions.value(find(ismember(inp.namoptions.option, {'&NAMGENSTAT'}))+1) ,'.true.') 
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   
    fprintf(['\n                  Loading output of the NAMGENSTAT routine.\n'])
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   
    tic

    % Call the function that reads the data of the NAMGENSTAT
    % outputs.
    f__NAMGENSTATread(direc, inp, check);
        % Output variables: structure with a separate field for 
        % each output file, preferably filled with tables.
        % Input variables:  direc, inp and the NAMGENSTAT field of
        % addons.

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
    fprintf(['                  - '])
    toc   
    fprintf(['                  Finished loading and saving NAMGENSTAT.\n'])
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   
    
else
    fprintf(['              --> ERROR!\n' ...
             '                  The mandatory namgenstat output is not available for this run.\n'])
end



%% Load the timestat output
    % NOTE: This dataset is mandatory
        
% Check if the NAMTIMESTAT outputs is .true. in namoptions 
if strcmp( inp.namoptions.value(find(ismember(inp.namoptions.option, {'&NAMTIMESTAT'}))+1) ,'.true.')
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   
    fprintf(['\n                  Loading output of the NAMTIMESTAT routine.\n'])
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   
    tic

    % Call the function that reads the data of the NAMTIMESTAT
    % outputs.
    f__NAMTIMESTATread(direc, inp, check);
        % Output variables: structure with a separate field for 
        % each output file, preferably filled with tables.
        % Input variables:  direc, inp and the NAMTIMESTAT field of
        % addons.

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
    fprintf(['                  - '])
    toc  
    fprintf(['                  Finished loading NAMTIMESTAT.\n'])
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  

else
    fprintf(['              --> ERROR!\n' ...
             '                  The mandatory namtimestat output is not available for this run.\n'])
end

% %% Load the cross-section
% 
% 
% 
% 
% % Make sure the cross-section settings are defined in namoptions
% if ismember({'lcross'}, inp.namoptions.option)
%     % Make sure the cross-section output is turned on and is set to be
%     % loaded & saved (check.process_cross)
%     if strcmp(inp.namoptions.value{ismember(inp.namoptions.option, {'lcross'})}, '.true.') 
%         
%         
%         
%         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
%         fprintf(['\n                  Loading output of the NAMCROSSSECTION routine.\n'])
%         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%         tic
%         
%         % Call the function that reads the basic variables of the 
%         % NAMCROSSSECTION outputs.
%         f__NAMCROSSSECTIONread_basic(direc, inp, check)
%         
%         
%         % XY cross-section
%         if ismember({['xy' check.cross_lvl{1}]}, check.process_cross)
%             % Call the function that reads the xy data of the NAMCROSSSECTION
%             % outputs.
%             f__NAMCROSSSECTIONread_xy(direc, inp, check);
%         else
%         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
%             fprintf(['\n                  - XY cross-section will not be loaded\n'])
%         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%         end
%         
%         % XZ cross-section
%         if ismember({'xz'}, check.process_cross)
%             % Call the function that reads the xz data of the NAMCROSSSECTION
%             % outputs.
%             f__NAMCROSSSECTIONread_xz(direc, inp, check);
%         else
%         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
%             fprintf(['\n                  - XZ cross-section will not be loaded\n'])
%         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%         end
%         
%         % YZ cross-section
%         if ismember({'yz'}, check.process_cross)
%             % Call the function that reads the yz data of the NAMCROSSSECTION
%             % outputs.
%             f__NAMCROSSSECTIONread_yz(direc, inp, check);
%         else
%         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
%             fprintf(['\n                  - YZ cross-section will not be loaded\n'])
%         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%         end
%         
%             
%         
%         clear JJ KK k kk cmb_* IN_switch outp_cross
% 
%         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%         fprintf(['                  - '])
%         toc  
%         fprintf(['                  Finished loading and saving NAMCROSSSECTION.\n'])
%         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
%         
%     end
% end
% 
% %% Load the fielddump xy cross-sections
% 
% % Make sure the cross-section settings are defined in namoptions
% if ismember({'lfielddump'}, inp.namoptions.option)
%     % Make sure the cross-section output is turned on and is set to be
%     % loaded & saved (check.process_cross)
%     if strcmp(inp.namoptions.value{ismember(inp.namoptions.option, {'lfielddump'})}, '.true.') 
%         
%         
%         for II = 2 : length(check.cross_lvl)
%             if ismember({['xy' check.cross_lvl{II}]}, check.process_cross)
%                 % Call the function that reads the xy data of the
%                 % NAMFIELDDUMP outputs.
%                 f__NAMFIELDDUMPread_xy(direc, inp, check, II)
%             end
%         end
%             
%         
%     end
% end
% 
% %% In case of warmstart, add the output of the coldstart run to the output of this warmstart run
% 
% % If the run is warmstarted with a different run, add the data from the 
% % warmstart run to the output
% if ismember({'lwarmstart'}, inp.namoptions.option) && ismember({'startfile'}, inp.namoptions.option)
%     if strcmp(inp.namoptions.value(ismember(inp.namoptions.option, {'lwarmstart'})), '.true.')
%         
%         % Define the experiment number of the warmstart file
%         cold_exp_nr = inp.namoptions.value(ismember(inp.namoptions.option, {'startfile'}));
%         if strcmp(cold_exp_nr{1}(end), '''')
%             cold_exp_nr = cold_exp_nr{1}(end-3:end-1);
%         else
%             cold_exp_nr = cold_exp_nr{1}(end-2:end);
%         end
%         cold_matdir = dir([direc.main_directory direc.run_name cold_exp_nr direc.dir_separator 'output_v*']);
%         cold_matdir = {cold_matdir.name}';
%         cold_matdir = cold_matdir{end};
%         cold_outp_files = dir([direc.main_directory direc.run_name cold_exp_nr direc.dir_separator cold_matdir]);
%         cold_outp_files = {cold_outp_files.name}';
%         
%         if strcmp(cold_exp_nr, direc.exp_nr(1:end-1)) == 0 && ...
%                 ismember({[cold_exp_nr 'cross.mat']}, cold_outp_files)
%             
%     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
%             fprintf(['\n                  Add the cross-section outputs of run ' ...
%                 cold_exp_nr ' to run ' direc.exp_nr(1:end-1) '\n']) 
%     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%             
%             % Add the cross-section output of the coldstart run to the
%             % cross-section output of this warmstart run
%             f__CROSScoldstartTOwarmstart(direc, inp)
%             
%             % Clear data
%             clear cold_*
%             
%         end     % <-- if strcmp(cold_exp_nr, direc.exp_nr(1:end-1))
%         
%     end     % <-- strcmp(inp.namoptions.value(ismember(inp.namoptions.option, {'lwarmstart'})), '.true.')
% end     % <-- ismember({'lwarmstart'}, inp.namoptions.option) && ismember({'startfile'}, inp.namoptions.option)



%% End the timer and display the runtime of the script

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf('  DATA IS LOADED AND SAVED!\n') 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Calculate time spent loading the data
tload_end = datestr(now, 'dd-mm-yyyy HH:MM:SS');
tload_temp = datenum(tload_end, 'dd-mm-yyyy HH:MM:SS')  - ...
    datenum(tload_start, 'dd-mm-yyyy HH:MM:SS');
tload_day   = tload_temp - mod(tload_temp,1);      
tload_temp  = tload_temp - tload_day;
tload_hour  = tload_temp - mod(tload_temp,1/24);  
tload_temp  = tload_temp - tload_hour;
tload_min   = tload_temp - mod(tload_temp,1/24/60);
tload_temp  = tload_temp - tload_min;
tload_sec   = tload_temp - mod(tload_temp,1/24/60/60);
tload_spent = [tload_day, tload_hour*24, tload_min*24*60, tload_sec*24*60*60];
% Message that loading the data is finished
fprintf('\n')
fprintf(['- DONE! Data is LOADED and SAVED.\n' ...
         '       function ended at: ' tload_end '\n' ...
         '       Time spent = '   num2str(tload_spent(1)) ' day(s)' ...
         '\n                    ' num2str(tload_spent(2)),' hour(s)' ...
         '\n                    ' num2str(tload_spent(2)),' hour(s)' ...
         '\n                    ' num2str(tload_spent(3)),' minutes(s) ' ...
         '\n                    ' num2str(tload_spent(4)) ' second(s)\n'])
        
% Remove temporary variables
clear tload_* data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%





























%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%                                                       %%%%%%%%%%
%%%%%%%%%%   Functions reading input files of the DALES addons   %%%%%%%%%%
%%%%%%%%%%                                                       %%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%









%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%% READ namoptions %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function func_out = f__READnamoptions(namoptions_dir)
    % The code underneath is generated by the matlab "Import Data"
    % functionality, based on the file of the aerosolrad_simple/001 run
    
    % Initialize variables.
    filename = namoptions_dir;
    delimiter = {'=','!'};

    % Read columns of data as text:
        % For more information, see the TEXTSCAN documentation.
    formatSpec = '%s%s%s%s%s%[^\n\r]';

    % Open the text file.
    fileID = fopen(filename,'r');
    
    % Read columns of data according to the format.
        % This call is based on the structure of the file used to generate this
        % code. If an error occurs for a different file, try regenerating the code
        % from the Import Tool.
    dataArray = textscan(fileID, formatSpec, 'Delimiter', delimiter, 'TextType', 'string',  'ReturnOnError', false);

    % Close the text file.
    fclose(fileID);

    % Convert the contents of columns containing numeric text to numbers.
        % Replace non-numeric text with NaN.
    raw = repmat({''},length(dataArray{1}),length(dataArray)-1);
    for col=1:length(dataArray)-1
        raw(1:length(dataArray{col}),col) = mat2cell(dataArray{col}, ones(length(dataArray{col}), 1));
    end
    numericData = NaN(size(dataArray{1},1),size(dataArray,2));

    % Converts text in the input cell array to numbers. Replaced non-numeric
    % text with NaN.
    rawData = dataArray{4};
    for row=1:size(rawData, 1)
        % Create a regular expression to detect and remove non-numeric prefixes and
        % suffixes.
        regexstr = '(?<prefix>.*?)(?<numbers>([-]*(\d+[\,]*)+[\.]{0,1}\d*[eEdD]{0,1}[-+]*\d*[i]{0,1})|([-]*(\d+[\,]*)*[\.]{1,1}\d+[eEdD]{0,1}[-+]*\d*[i]{0,1}))(?<suffix>.*)';
        try
            result = regexp(rawData(row), regexstr, 'names');
            numbers = result.numbers;

            % Detected commas in non-thousand locations.
            invalidThousandsSeparator = false;
            if numbers.contains(',')
                thousandsRegExp = '^\d+?(\,\d{3})*\.{0,1}\d*$';
                if isempty(regexp(numbers, thousandsRegExp, 'once'))
                    numbers = NaN;
                    invalidThousandsSeparator = true;
                end
            end
            % Convert numeric text to numbers.
            if ~invalidThousandsSeparator
                numbers = textscan(char(strrep(numbers, ',', '')), '%f');
                numericData(row, 4) = numbers{1};
                raw{row, 4} = numbers{1};
            end
        catch
            raw{row, 4} = rawData{row};
        end
    end

    % Split data into numeric and string columns.
    rawNumericColumns = raw(:, 4);
    rawStringColumns = string(raw(:, [1,2,3,5]));

    % Replace non-numeric cells with NaN
    R = cellfun(@(x) ~isnumeric(x) && ~islogical(x),rawNumericColumns); % Find non-numeric cells
    rawNumericColumns(R) = {NaN}; % Replace non-numeric cells
    
    % Create output variable
    func_out = table;
    func_out.option = rawStringColumns(:, 1);
    func_out.value = rawStringColumns(:, 2);
    func_out.description1 = rawStringColumns(:, 3);
    func_out.description2 = cell2mat(rawNumericColumns(:, 1));
    func_out.description3 = rawStringColumns(:, 4);

    % Clean up namoptions.option and namoption.value
    for i = 1 : length(func_out.option)
        
        % namoption.option: remove spaces
        IN_space = [1:1:length(func_out.option{i})];
        IN_space(ismember(func_out.option{i},' ')) = [];
        func_out.option{i} = func_out.option{i}(IN_space);
        % namoption.option: remove tabs
        IN_tab = [1:1:length(func_out.option{i})];
        IN_tab(ismember(func_out.option{i},'	')) = [];
        func_out.option{i} = func_out.option{i}(IN_tab);
        
        % namoption.value: remove spaces
        IN_space = [1:1:length(func_out.value{i})];
        IN_space(ismember(func_out.value{i},' ')) = [];
        func_out.value{i} = func_out.value{i}(IN_space);
        % namoption.value: remove tabs
        IN_tab = [1:1:length(func_out.value{i})];
        IN_tab(ismember(func_out.value{i},'	')) = [];
        func_out.value{i} = func_out.value{i}(IN_tab);
    end

    % Clear temporary variables
    clearvars filename delimiter formatSpec fileID dataArray ans raw col ...
        numericData rawData row regexstr result numbers ...
        invalidThousandsSeparator thousandsRegExp rawNumericColumns rawStringColumns R;
    
    
end


%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%% READ input with textscan function %%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function func_out = f__READinput_textscan(lscale_dir, inp)
   
    % Open the lscale.inp file
    fID = fopen(lscale_dir);
    
    % Get the Nz of the run
    Nz = str2num(inp.namoptions.value(ismember(inp.namoptions.option, {'kmax'})));
    
    % Get the headers of the lscale.inp file
        % Set the delimiters of the second line of the file
    header_delimiter = {' ',',','\t'};
        % read the second line, ignoring #
    header_row = textscan(fID, '%*1s %q %q %q %q %q %q %q %q %q %q %q %q %q %q %q %q %q %q %q %q %q %q %q %q %q %q %q %q %q %q %q %q %q %q %q %q %q %q %q %q %q %q', ...
        1, 'Delimiter', header_delimiter, 'TreatAsEmpty', {'#'}, 'HeaderLines',1);
        % Reorganize the variable
    empty_columns = NaN(size(header_row));
	for i = 1 : length(header_row)
        header_row{i} = header_row{i}{1};
        empty_columns(i) = isempty(header_row{i});
    end
    header_row = header_row(empty_columns == 0);
    % Read the actual data (line 3 to 3+Nz)
    data_format = [''];
    for i = 1 : length(header_row)
        data_format = [data_format '%f '];
    end
    file_data = textscan(fID, data_format, Nz, 'HeaderLines',1);
    
    % Close the lscale.inp  file
    fclose(fID);
    
    %%%%%%
    % Save the data as a table, with the headers as table variable names
    func_out = table();
    for i = 1 : length(header_row)
        % Matlab does not recognize "(" or ")" as variable names. Replace
        % them with "_"
        IN_replace = find(ismember(header_row{i},'('));
        IN_replace = [IN_replace find(ismember(header_row{i},')'))];
        header_row{i}(IN_replace) = '_';
        
        % Save the data in the table
        func_out.(header_row{i}) = file_data{i};
    end
    
end










%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%                                                       %%%%%%%%%%
%%%%%%%%%%  Functions reading output files of the DALES addons   %%%%%%%%%%
%%%%%%%%%%                                                       %%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%      READ NAMGENSTAT output       %%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function f__NAMGENSTATread(direc, inp, check)
    % This function reads the data of the NAMGENSTAT outputs.
        % Input variables:  direc, inp and the NAMGENSTAT field of addons
        % Output variables: structure with a separate field for each output 
        %                   file, preferably filled with tables.
    
    
    % Preallocation output structure
    outp_gen = struct();
    % Clarify how the output is saved
    outp_gen.README = {'Each variable is saved as an [height x time] function'};
    % Preallocate the info on the variables
    outp_gen.INFO = table({},{},{}, 'VariableNames',{'Name','LongName','Unit'});
    
    %% Read the NetCDF output
    
    % Make sure the "lnetcdf" switch is not set to " .false." 
    nc_exists = 1;  % 1 = NetCDF files exist, 0 = NetCDF files do not exist
    if ismember({'lnetcdf'}, inp.namoptions.option)
        if strcmp(inp.namoptions.value{...
                ismember(inp.namoptions.option,{'lnetcdf'})}, '.false.')
            nc_exists = 0;
        end
    end
    
    
    if nc_exists == 1
        % Define location of the NetCDF file
        nc_loc = [direc.main_directory, direc.run_name, direc.exp_nr 'profiles.' direc.exp_nr(1:end-1) '.nc'];
        % Read/save information on the NetCDF file
        nc_info = ncinfo(nc_loc);
        % Read the variable names that are in the NetCDF file
        nc_vars = {nc_info.Variables.Name}';
        
        % Preallocate temporary variable
        nc_units = table(cell(length(nc_vars),1), ...
                         cell(length(nc_vars),1), ...
                         cell(length(nc_vars),1), ...
                         'VariableNames',{'Name','LongName','Unit'});
        
        % Load/read all variables in the NetCDF files
        for i = 1 : length(nc_vars)
            % Read and save the variable data
            outp_gen.(nc_vars{i}) = ncread(nc_loc, nc_vars{i});
            % Read the variable information and store it in a temporary
            % variable
            temp_units = ncinfo(nc_loc, nc_vars{i});
            nc_units(i,1:3) = [{temp_units.Name}, ...
                               {temp_units.Attributes(1,1).Value}, ...
                               {temp_units.Attributes(1,2).Value}];
        end
        % Save the variable information in the temporary variable and store it
        % in the data sructure
        outp_gen.INFO = [outp_gen.INFO ; nc_units];
    else
    
        fprintf(['              --> There is no NetCDF output.\n' ...
                 '                - Please turn on the "lnetcdf" switch in Namoptions.\n'])
    end
    
    %% Rename the scalars in the output to match the names in scalar.inp
    
    % Get the input names of the scalars
    corr_inpscalars = inp.scalar.Properties.VariableNames(2:end);
    
    % Preallocation of matrix with names to correct
    corr_corrections = [{}, {}]; 
    % Get the LongNames as saved in f_output.INFO (will be used to search in)
    temp_longnames = (outp_gen.INFO.LongName);
    % Loop over the input names of the scalars to correct all of them
    temp_varfields = fields(outp_gen);
    temp_varfields(ismember(temp_varfields, {'README','INFO'})) = [];
    for i = 1 : length(corr_inpscalars)
        
        
        % Create the name that will be used to search in f_output.(cross_fields).INFO.LongNames
        corr_name = '000';
        corr_name(end-(length(num2str(i))-1) : end) = num2str(i);

        corr_name = [{['sv' corr_name]},{['Scalar ' corr_name]}]; 

        corr_IN = find(ismember(outp_gen.INFO.Name, corr_name(1)));
        
        if isempty(corr_IN) == 0
            outp_gen.INFO.Name{corr_IN} = corr_inpscalars{i};

            corr_longname = outp_gen.INFO.LongName{corr_IN};
            for j = 1 : length(corr_longname) - length(corr_name{2}) + 1
                if strcmp(corr_name{2}, corr_longname(j : j+length(corr_name{2})-1))
                    if j == 1
                        outp_gen.INFO.LongName{corr_IN} = ...
                            ['Scalar ', corr_inpscalars{i}, ...
                            corr_longname(j+length(corr_name{2}):end)];
                    else
                        outp_gen.INFO.LongName{corr_IN} = ...
                            [corr_longname(1:j-1), ' scalar ', corr_inpscalars{i}, ...
                            corr_longname(j+length(corr_name{2}):end)];
                    end
                end
            end
            corr_corrections = [corr_corrections ; ...
                temp_varfields(corr_IN), corr_inpscalars(i) ];
        end
        
        
    end
    
    if isempty(corr_corrections) == 0
        % Save the original field names of f_output
        outp_fields_old = fields(outp_gen);
        % Save the original f_output
        outp_gen_old = outp_gen;
        % Save the original field names of f_output

        % Base the new f_output field names on the old f_output field names
        outp_fields_new = outp_fields_old;
        % Apply the corrections to the f_output field names
        outp_fields_new(ismember(outp_fields_old, corr_corrections(:,1))) = ...
            corr_corrections(:,2);

        % Build a new and empty f_output structure
        outp_gen = struct();
        % Redefine f_output with the same data, but new names for scalars
        for i = 1 : length(outp_fields_new)
            outp_gen.(outp_fields_new{i}) = outp_gen_old.(outp_fields_old{i});
        end
    end
    
    
    %% Read the txt output
    
    % Define the location & name of the txt file you want to read
    txt_loc = [direc.main_directory, direc.run_name, direc.exp_nr 'field.' direc.exp_nr(1:end-1)];
%     txt_loc = 'C:\Users\schul068\PhD project\DALES\Runs\nh3plume\201\field.201'  
    %%%%% Get the headers of the txt file
    % Define the delimiters of header line
    txt_header_delimiter = {' '};
    % Define the line number where the header is located
    txt_header_linenr = 5;
    
    % Open the txt file
    txt_fID = fopen(txt_loc);
    % Read the line as cell array (it is made very long as we do not know
    % how many columns the file has and how many spaces are used between
    % headers
    txt_header_temp = textscan(txt_fID, '%q %q %q %q %q %q %q %q %q %q %q %q %q %q %q %q %q %q %q %q %q %q %q %q %q %q %q %q %q %q %q %q %q %q %q %q %q %q %q %q %q %q %q %q %q %q %q %q %q %q %q %q %q %q %q %q %q %q %q %q %q %q %q %q %q %q %q %q %q %q %q %q %q %q %q %q %q %q %q %q %q %q %q %q %q %q %q %q %q %q %q %q %q %q %q %q %q %q %q %q %q %q %q %q ', ...
        1, 'Delimiter', txt_header_delimiter, 'TreatAsEmpty', {'#'}, 'HeaderLines',txt_header_linenr);
    % Close the txt file
    fclose(txt_fID);
    
    %%%%% Restructure the header into a propper cell array
    % Preallocation
    txt_header = cell(1, length(txt_header_temp));
    txt_temp_prev = '';
    % Loop backwards over the temporary header cell array
    for i = length(txt_header_temp) : -1 : 1
        % Get the characters in the cell
        txt_temp = txt_header_temp{i};
        txt_temp = txt_temp{:};
        % Only keep the output when the cell contains an actual header
        if isempty(txt_temp)
            txt_header(i) = [];
        else
            % Ignore the "#" character
            txt_temp(ismember(txt_temp,'#')) = [];
            % Save the header
            txt_header{i} = txt_temp;
        end
        
        % The CLOUD FRACTION includes a space, which is a delimiter. Fix
        % this!
        if 1 < length(txt_header_temp)
            if strcmp(txt_temp, 'CLOUD') && strcmp(txt_temp_prev, 'FRACTION')
                txt_header{i} = 'CLOUD_FRACTION';
                txt_header(i+1) = [];
            end
        end
        % Save the content of the previous cell. Required for the above fix
        txt_temp_prev = txt_temp;
    end

    %%%%% Now read the actual data
    % Define the number of lines the headers take
    txt_ln_header = 7;
    % Define the number of lines the data takes (= kmax in NAMOPTIONS)
    txt_ln_data = str2num(inp.namoptions.value(ismember(inp.namoptions.option, {'kmax'})));
    % Define the number of times the data is written by DALES
    INgenstat = find(ismember(inp.namoptions.option, {'&NAMGENSTAT'}));
    INgen_tav = find(ismember(inp.namoptions.option(INgenstat:end), {'timeav'})) + INgenstat-1;
    txt_ln_time = str2num(inp.namoptions.value(ismember(inp.namoptions.option, {'runtime'}))) / ...
                  str2num(inp.namoptions.value{INgen_tav(1)});
    txt_ln_time = floor(txt_ln_time);
    
    % Define the delimiter of the data file
    txt_delimiter = {'\t'};
    
    % Preallocation
    txt_skip = 0;                                   % Lines to skip until data starts
    txt_data = NaN(txt_ln_data, 13, txt_ln_time);   % 3D matrix to write the data in
    for i = 1 :  txt_ln_time
        % Skip the header
        txt_skip = txt_skip + txt_ln_header;
        
        % Open the file
        txt_fID = fopen(txt_loc);
        % Read the data in a temporary variable
        temp_data = textscan(txt_fID, '%3q %f %f %f %f %f %f %f %f %f %f %f %f', ...
            txt_ln_data, 'Delimiter', txt_delimiter, 'HeaderLines',txt_skip);
        % Close the file
        fclose(txt_fID);
        temp_data{1,1} = [1 : 1 : txt_ln_data]';
        
        % Save the data from the temporary variable into the proper 3D matrix
        for j = 1 : 13
            txt_data(:,j,i) = temp_data{:,j};
        end 
        % Skip the already saved data
        txt_skip = txt_skip + txt_ln_data;
    end
    
    %% Specific things
    
    % Add the temperature data to f_output
    if sum(ismember(txt_header, 'TEMP')) == 1
        outp_temp = NaN(txt_ln_data, txt_ln_time);
        for i = 1 : txt_ln_time
            outp_temp(:,i) = txt_data(:,ismember(txt_header, 'TEMP'),i);
        end
        outp_gen.temp = outp_temp;
        outp_gen.INFO = [outp_gen.INFO; table({'temp'}, {'Temperature'}, {'K'}, ...
            'VariableNames',{'Name','LongName','Unit'})];
    end
    
    % Add the potential temperature data to f_output
    if sum(ismember(txt_header, 'THETA')) == 1
        outp_temp = NaN(txt_ln_data, txt_ln_time);
        for i = 1 : txt_ln_time
            outp_temp(:,i) = txt_data(:,ismember(txt_header, 'THETA'),i);
        end
        outp_gen.th = outp_temp;
        outp_gen.INFO = [outp_gen.INFO; table({'th'}, {'Potential temperature'}, {'K'}, ...
            'VariableNames',{'Name','LongName','Unit'})];
    end
    
    % Add the cloud fraction data to f_output
    if sum(ismember(txt_header, 'CLOUD_FRACTION')) == 1
        outp_temp = NaN(txt_ln_data, txt_ln_time);
        for i = 1 : txt_ln_time
            outp_temp(:,i) = txt_data(:,ismember(txt_header, 'CLOUD_FRACTION'),i);
        end
        outp_gen.cf = outp_temp;
        outp_gen.INFO = [outp_gen.INFO; table({'cf'}, {'Cloud fraction'}, {'-'}, ...
            'VariableNames',{'Name','LongName','Unit'})];
    end
    
    % Define the time array in hours local time
    outp_gen.time = outp_gen.time ./ 3600 + check.corr_localtime + ....
                str2num(inp.namoptions.value(ismember(inp.namoptions.option, {'xtime'})));
    
    
    
    %
    %%%%% TKE resolved
    %
    try
        outp_gen.tker = 0.5* ( outp_gen.u2r + ...
            outp_gen.v2r + outp_gen.w2r );
    catch
        outp_gen.tker = NaN(size(outp_gen.time));
    end

    outp_gen.INFO = [outp_gen.INFO ; ...
        table({'tker'},{'Resolved Turbulent Kinetic Energy'},{'m^2/s^2'}, ...
        'VariableNames',{'Name','LongName','Unit'})];

    %
    %%%%% TKE total
    %
    outp_gen.tke = outp_gen.tker + outp_gen.w2s;

    outp_gen.INFO = [outp_gen.INFO ; ...
        table({'tke'},{'Total Turbulent Kinetic Energy'},{'m^2/s^2'}, ...
        'VariableNames',{'Name','LongName','Unit'})];

    %
    %%%%% TKE resolved over total
    %
    outp_gen.tke_rot = abs(outp_gen.tker) ./ abs(outp_gen.tke);

    outp_gen.INFO = [outp_gen.INFO ; ...
        table({'tke_rot'},{'Resolved over total turbulent kinetic energy'},{'m^2/s^2'}, ...
        'VariableNames',{'Name','LongName','Unit'})];
    
    

    if strcmp(check.combine, 'yes')
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        fprintf('                  - Combine variables in the output\n') 
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        %%%%%
        % Combine variables in NAMGENSTAT data
        %%%%%
        
            
        for KK = 1 : size(check.combine_vars,1)

            if isempty(fields(outp_gen)) == 0 && ismember({check.combine_newname{KK}} , fields(outp_gen)) == 0


                cmb_fields = fields(outp_gen);
                
                if sum(ismember(cmb_fields, check.combine_vars(KK,:))) == length(check.combine_vars(KK,:))

                    % update the INFO table
                    cmb_infoIN = find( ismember(...
                        outp_gen.INFO.Name, check.combine_vars{KK,1}));
                    cmb_longname = outp_gen.INFO.LongName(cmb_infoIN);
                    cmb_longname = cmb_longname{:};

                    for JJ = 1 : length(cmb_longname) - length(check.combine_vars{KK,1}) + 1
                        if strcmp(check.combine_vars{KK,1}, cmb_longname(JJ : JJ+length(check.combine_vars{KK,1})-1))
                            if JJ == 1
                                cmb_tempname = ...
                                    [check.combine_newname{KK}, ...
                                    corr_longname(JJ+length(check.combine_vars{KK,1}):end)];
                            else
                                cmb_tempname = ...
                                    [cmb_longname(1:JJ-1), check.combine_newname{KK}, ...
                                    cmb_longname(JJ+length(check.combine_vars{KK,1}):end)];
                            end
                        end
                    end
                    outp_gen.INFO = [ outp_gen.INFO ; ...
                            table(check.combine_newname(KK), {cmb_tempname}, ...
                            outp_gen.INFO.Unit(cmb_infoIN), 'VariableNames',{'Name','LongName','Unit'})];

                    %  Add the new variable
                    outp_gen.(check.combine_newname{KK}) = outp_gen.(check.combine_vars{KK,1});
                    for k = 2 : length(check.combine_vars(KK,:))
                        outp_gen.(check.combine_newname{KK}) = outp_gen.(check.combine_newname{KK}) + ...
                            outp_gen.(check.combine_vars{KK,k});
                    end

                end     % sum(ismember(cmb_fields, check.combine_vars(KK,:))) == length(check.combine_vars(KK,:))

            end     % isempty(fields(outp_gen)) == 0 && ismember({check.combine_newname{KK}} , fields(outp_gen)) == 0

        end     % for KK = 1 : size(check.combine_vars,1)

    end
    
    % If the run is warmstarted with a different run, add the data from the 
    % warmstart run to the output
    if ismember({'lwarmstart'}, inp.namoptions.option) && ismember({'startfile'}, inp.namoptions.option)
        if strcmp(inp.namoptions.value(ismember(inp.namoptions.option, {'lwarmstart'})), '.true.')
            
            % Define the experiment number of the warmstart file
            ws_exp_nr = inp.namoptions.value(ismember(inp.namoptions.option, {'startfile'}));
            if strcmp(ws_exp_nr{1}(end), '''')
                ws_exp_nr = ws_exp_nr{1}(end-3:end-1);
            else
                ws_exp_nr = ws_exp_nr{1}(end-2:end);
            end
            
            % Define location of the .mat files of the warmstart experiment
            ws_matdir = dir([direc.main_directory direc.run_name ws_exp_nr direc.dir_separator 'output_v*']);
            ws_matdir = {ws_matdir.name}';
            ws_matdir = ws_matdir{end};
            % Define the names of the .mat files of the warmstart experiment
            ws_namelist = dir([direc.main_directory direc.run_name ws_exp_nr direc.dir_separator ...
                ws_matdir direc.dir_separator '*.mat']);
            ws_namelist = {ws_namelist.name}';
            % Check if the NAMGENSTAT output exists
            if ismember({[ws_exp_nr 'gen.mat']}, ws_namelist) && ...
                    strcmp(ws_exp_nr, direc.exp_nr(1:end-1)) == 0
                
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
                fprintf(['                  - Add the NAMGENSTAT output of run ' ...
                    ws_exp_nr ' to run ' direc.exp_nr(1:end-1) '\n']) 
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

                % Load the input data of the cold start run
                cs_outp = load([direc.main_directory direc.run_name ...
                    ws_exp_nr direc.dir_separator ws_matdir direc.dir_separator ...
                    ws_exp_nr 'gen.mat']);
                cs_fields = fields(cs_outp);
                cs_fields(ismember(cs_fields,{'README','INFO','zt','zm','zts'})) = [];
                
                % Define the fields of the warsm start run
                ws_fields = fields(outp_gen);
                ws_fields(ismember(ws_fields,{'README','INFO','zt','zm','zts'})) = [];
                ws_fields = sort(ws_fields);
                
                % Add the warmstart output to the output of this run
                for i = 1 : length(ws_fields)
                    if strcmp(ws_fields{i}, 'time')
                        outp_gen.(ws_fields{i}) = [cs_outp.(cs_fields{i}) ; outp_gen.(ws_fields{i})];
                    else
                        outp_gen.(ws_fields{i}) = [cs_outp.(cs_fields{i}) , outp_gen.(ws_fields{i})];
                    end
                end

            end     % <-- if ismember({[ws_exp_nr 'gen.mat']}, ws_namelist) && ...
            clear ws_*
        end     % <-- strcmp(inp.namoptions.value(ismember(inp.namoptions.option, {'lwarmstart'})), '.true.')
    end     % <-- ismember({'lwarmstart'}, inp.namoptions.option) && ismember({'startfile'}, inp.namoptions.option)
    
    
    
    % Save the outp_gen data
    save([direc.main_directory direc.run_name direc.exp_nr ...                  % Save directory
        direc.output_dir direc.dir_separator direc.exp_nr(1:end-1) 'gen.mat'], ...     % Save name
        '-struct', 'outp_gen', '-v7.3');                                                   % Saved variables
    
    
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%      READ NAMTIMESTAT output      %%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function f__NAMTIMESTATread(direc, inp, check)
    % This function reads the data of the NAMTIMESTAT outputs.
        % Input variables:  direc, inp and the NAMTIMESTAT field of addons
        % Output variables: structure with a separate field for each output 
        %                   file, preferably filled with tables.
    
    
    % Preallocation output structure
    outp_time = struct();
    % Clarify how the output is saved
    outp_time.README = {'Each variable is saved as an [time x 1] vector'};
    % Preallocate the info on the variables
    outp_time.INFO = table({},{},{}, 'VariableNames',{'Name','LongName','Unit'});
    
    %% Read the NetCDF output
    
    % Make sure the "lnetcdf" switch is not set to " .false." 
    nc_exists = 1;  % 1 = NetCDF files exist, 0 = NetCDF files do not exist
    if ismember({'lnetcdf'}, inp.namoptions.option)
        if strcmp(inp.namoptions.value{...
                ismember(inp.namoptions.option,{'lnetcdf'})}, '.false.')
            nc_exists = 0;
        end
    end
    
    if nc_exists == 1
        % Define location of the NetCDF file
        nc_loc = [direc.main_directory, direc.run_name, direc.exp_nr 'tmser.' direc.exp_nr(1:end-1) '.nc'];
        % Read/save information on the NetCDF file
        nc_info = ncinfo(nc_loc);
        % Read the variable names that are in the NetCDF file
        nc_vars = {nc_info.Variables.Name}';
        
        % Preallocate temporary variable
        nc_units = table(cell(length(nc_vars),1), ...
                         cell(length(nc_vars),1), ...
                         cell(length(nc_vars),1), ...
                         'VariableNames',{'Name','LongName','Unit'});
        
        % Load/read all variables in the NetCDF files
        for i = 1 : length(nc_vars)
            % Read and save the variable data
            outp_time.(nc_vars{i}) = ncread(nc_loc, nc_vars{i});
            % Read the variable information and store it in a temporary
            % variable
            temp_units = ncinfo(nc_loc, nc_vars{i});
            nc_units(i,1:3) = [{temp_units.Name}, ...
                               {temp_units.Attributes(1,1).Value}, ...
                               {temp_units.Attributes(1,2).Value}];
        end
        % Save the variable information in the temporary variable and store it
        % in the data sructure
        outp_time.INFO = [outp_time.INFO ; nc_units];
    else
    
        fprintf(['              --> There is no NetCDF output.\n' ...
                 '                - Please turn on the "lnetcdf" switch in Namoptions.\n'])
    end
    
    %% Read the txt output to add the missing variables
    
    % Define the location & name of the txt file you want to read
    txt_loc = [direc.main_directory, direc.run_name, direc.exp_nr 'tmlsm.' direc.exp_nr(1:end-1)];
    
    %%%%% Get the headers of the txt file
    % Open the txt file
    txt_fID = fopen(txt_loc);
    temp_vartype = repmat({'%q '},[999,1]); 
    temp_vartype = [temp_vartype{:}];
    % Read the header lines as cell array 
    txt_temp = textscan(txt_fID, '%q ', 8, 'Delimiter', '', 'TreatAsEmpty', {'#'}, 'HeaderLines',0);
    txt_temp = txt_temp{1};
    txt_IN_header = [];
    for i = 1 : length(txt_temp)
        if strcmp(txt_temp{i}(1), '#')
            txt_IN_header = [txt_IN_header; i];
        end
    end
    % Close the txt file
    fclose(txt_fID);
    
    if isempty(txt_IN_header)
        txt_oldexprnr = inp.namoptions.value(ismember(inp.namoptions.option, {'startfile'}));
        if strcmp(txt_oldexprnr{1}(end), '''')
            txt_oldexprnr = txt_oldexprnr{1}(end-3:end-1);
        else
            txt_oldexprnr = txt_oldexprnr{1}(end-2:end);
        end
        
        txt_fID2 = fopen([direc.main_directory, direc.run_name, ...
            txt_oldexprnr direc.dir_separator 'tmlsm.' txt_oldexprnr]);
        temp_vartype2 = repmat({'%q '},[999,1]); 
        temp_vartype2 = [temp_vartype2{:}];
        % Read the header lines as cell array 
        txt_temp2 = textscan(txt_fID2, '%q ', 8, 'Delimiter', '', 'TreatAsEmpty', {'#'}, 'HeaderLines',0);
        txt_temp2 = txt_temp2{1};
        txt_IN_header2 = [];
        txt_IN_datastart = [];
        for i = 1 : length(txt_temp2)
            if strcmp(txt_temp2{i}(1), '#')
                txt_IN_header2 = [txt_IN_header2; i];
            end
            if ismember({txt_temp2{i}(1)}, {'0','1','2','3','4','5','6','7','8','9'})
                txt_IN_datastart = [txt_IN_datastart; i];
            end
        end
        txt_IN_datastart = txt_IN_datastart(1);
        % Close the txt file
        fclose(txt_fID2);
        
        movefile(txt_loc, [txt_loc(1:end-4) '_noheader' txt_loc(end-3:end)]);
        
        txt_fID_new = fopen(txt_loc, 'w');
        for i = 1 : txt_IN_datastart-1 
            fprintf(txt_fID_new, [txt_temp2{i} '\n']);
        end
        txt_fID_old = fopen([txt_loc(1:end-4) '_noheader' txt_loc(end-3:end)]);
        txt_temp = textscan(txt_fID_old, '%q ', 'Delimiter', '', 'TreatAsEmpty', {'#'}, 'HeaderLines',0);
        txt_temp = txt_temp{1};
        fclose(txt_fID_old);
        
        for i = 1 : length(txt_temp)
            fprintf(txt_fID_new, [txt_temp{i} '\n']);
        end
        
        fclose(txt_fID_new);
        
        txt_IN_header = txt_IN_header2;
        
    end
    
    txt_headlength = txt_IN_header(2) - txt_IN_header(1);
    txt_header = repmat({''},[999,2]);
    for i = 1 : length(txt_IN_header)
        
        % Open the txt file
        txt_fID = fopen(txt_loc);
        txt_temp = textscan(txt_fID, '%q ', txt_headlength, ...
            'Delimiter', '', 'TreatAsEmpty', {'#'}, 'HeaderLines',txt_IN_header(i)-1);
        
        % Close the txt file
        fclose(txt_fID);
        
        
        txt_headertemp = [''];
        for j = 1 : length(txt_temp{:})
            txt_headertemp = [txt_headertemp, '   ', txt_temp{1,1}{j}];
        end
        txt_IN_nospace = find(ismember(txt_headertemp,' ')==0 & ismember(txt_headertemp,'#')==0);
        
        txt_counter = 1;
        for j = 1 : length(txt_IN_nospace)
            if j == 1
                txt_header{txt_counter,i} = [txt_header{txt_counter,i}, ...
                    txt_headertemp(txt_IN_nospace(j))];
            elseif (txt_IN_nospace(j) - txt_IN_nospace(j-1)) == 1
                txt_header{txt_counter,i} = [txt_header{txt_counter,i}, ...
                    txt_headertemp(txt_IN_nospace(j))];
            else
                txt_counter = txt_counter + 1;
                txt_header{txt_counter,i} = [txt_header{txt_counter,i}, ...
                    txt_headertemp(txt_IN_nospace(j))];
            end
        end
    end
    txt_header(ismember(txt_header(:,1),repmat({''},[1,1])),:) = [];
    
    % Manual unit corrections
    txt_header{ismember(txt_header(:,1),{'An'}),2} = '[mgCm2/s]';
    txt_header{ismember(txt_header(:,1),{'gcco2'}),2} = '[m/s]';
    
    
    %%%%% Now read the actual data
    % Define the number of times the data is written by DALES
    INtimestat = find(ismember(inp.namoptions.option, {'&NAMTIMESTAT'}));
    INtime_tav = find(ismember(inp.namoptions.option(INtimestat:end), {'dtav'})) + INtimestat-1;
    txt_ln_data = str2num(inp.namoptions.value(ismember(inp.namoptions.option, {'runtime'}))) / ...
                  str2num(inp.namoptions.value{INtime_tav(1)});
    
    % Define the delimiter of the data file
    txt_delimiter = {'\t'};
    % Define the number of lines the headers take
    txt_ln_header = txt_IN_header(end) + txt_headlength-1;
    
    % Preallocation
    txt_data = NaN(txt_ln_data, size(txt_header,1));   % 3D matrix to write the data in
   
    % Open the file
    txt_fID = fopen(txt_loc);
    %
    txt_varinfo = repmat({'%f '},[size(txt_header,1),1]); 
    txt_varinfo =  [txt_varinfo{:}];
    % Read the data in a temporary variable
    temp_data = textscan(txt_fID, txt_varinfo, ...
        txt_ln_data, 'Delimiter', txt_delimiter, 'HeaderLines',txt_ln_header);
    % Close the file
    fclose(txt_fID);

    % Save the data from the temporary variable into the proper 3D matrix
    for i = 1 : size(txt_header,1)
        if strcmp(txt_header{i}, 'time')
            txt_data(:,i) = temp_data{:,i} ./ 3600 + check.corr_localtime + ....
                str2num(inp.namoptions.value(ismember(inp.namoptions.option, {'xtime'})));
        else
            txt_data(:,i) = temp_data{:,i};
        end
    end
    
    
    txt_add  = {'tskin';'Resp';'wco2';'An';'gcco2'};
    txt_add_long = {'Skin temperature'; ...
                    'CO2 soil respiration'; ...
                    'Net CO2 flux'; ...
                    'Assimilation rate'; ...
                    'CO2 canopy conductance'};
    for i = 1 : length(txt_add)
        outp_time.(txt_add{i}) = ...
            txt_data(:, ismember(txt_header(:,1),txt_add(i)));
        outp_time.INFO = [outp_time.INFO ; ...
            table(txt_add(i),txt_add_long(i),txt_header(ismember(txt_header(:,1),txt_add(i)),2), ...
            'VariableNames',{'Name','LongName','Unit'})];
    end
    
    
    %% Rename the scalars in the output to match the names in scalar.inp
    
    % Get the input names of the scalars
    corr_inpscalars = inp.scalar.Properties.VariableNames(2:end);
    
    % Preallocation of matrix with names to correct
    corr_corrections = [{}, {}]; 
    % Get the LongNames as saved in f_output.INFO (will be used to search in)
    temp_longnames = (outp_time.INFO.LongName);
    % Loop over the input names of the scalars to correct all of them
    for i = 1 : length(corr_inpscalars)
        
        
        % Create the name that will be used to search in f_output.(cross_fields).INFO.LongNames
        corr_name = '000';
        corr_name(end-(length(num2str(i))-1) : end) = num2str(i);

        corr_name = [{['sv' corr_name]},{['Scalar ' corr_name]}]; 

        corr_IN = find(ismember(outp_time.INFO.Name, corr_name(1)));
        
        if isempty(corr_IN) == 0
            outp_time.INFO.Name{corr_IN} = corr_inpscalars{i};

            corr_longname = outp_time.INFO.LongName{corr_IN};
            for j = 1 : length(corr_longname) - length(corr_name{2}) + 1
                if strcmp(corr_name{2}, corr_longname(j : j+length(corr_name{2})-1))
                    if j == 1
                        outp_time.INFO.LongName{corr_IN} = ...
                            ['Scalar ', corr_inpscalars{i}, ...
                            corr_longname(j+length(corr_name{2}):end)];
                    else
                        outp_time.INFO.LongName{corr_IN} = ...
                            [corr_longname(1:j-1), ' scalar ', corr_inpscalars{i}, ...
                            corr_longname(j+length(corr_name{2}):end)];
                    end
                end
            end
        end
        
        
    end
    
    if isempty(corr_corrections) == 0
        % Save the original field names of f_output
        outp_fields_old = fields(outp_time);
        % Save the original f_output
        outp_time_old = outp_time;
        % Save the original field names of f_output

        % Base the new f_output field names on the old f_output field names
        outp_fields_new = outp_fields_old;
        % Apply the corrections to the f_output field names
        outp_fields_new(ismember(outp_fields_old, corr_corrections(:,1))) = ...
            corr_corrections(:,2);

        % Build a new and empty f_output structure
        outp_time = struct();
        % Redefine f_output with the same data, but new names for scalars
        for i = 1 : length(outp_fields_new)
            outp_time.(outp_fields_new{i}) = outp_time_old.(outp_fields_old{i});
        end
    end
    
    outp_time.time = outp_time.time ./ 3600 + check.corr_localtime + ....
                str2num(inp.namoptions.value(ismember(inp.namoptions.option, {'xtime'})));
    
    if strcmp(check.combine, 'yes')
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
        fprintf('                  - Combine variables in the NAMTIMESTAT output\n') 
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        %%%%%
        % Combine variables in NAMTIMESTAT data
        %%%%%
        for KK = 1 : size(check.combine_vars,1)

            if isempty(fields(outp_time)) == 0 && ismember({check.combine_newname{KK}} , fields(outp_time)) == 0


                cmb_fields = fields(outp_time);

                if sum(ismember(cmb_fields, check.combine_vars(KK,:))) == length(check.combine_vars(KK,:))

                    % update the INFO table
                    cmb_infoIN = find( ismember(...
                        outp_time.INFO.Name, check.combine_vars{KK,1}));
                    cmb_longname = outp_time.INFO.LongName(cmb_infoIN);
                    cmb_longname = cmb_longname{:};

                    for JJ = 1 : length(cmb_longname) - length(check.combine_vars{KK,1}) + 1
                        if strcmp(check.combine_vars{KK,1}, cmb_longname(JJ : JJ+length(check.combine_vars{KK,1})-1))
                            if JJ == 1
                                cmb_tempname = ...
                                    [check.combine_newname{KK}, ...
                                    corr_longname(JJ+length(check.combine_vars{KK,1}):end)];
                            else
                                cmb_tempname = ...
                                    [cmb_longname(1:JJ-1), check.combine_newname{KK}, ...
                                    cmb_longname(JJ+length(check.combine_vars{KK,1}):end)];
                            end
                        end
                    end
                    outp_time.INFO = [ outp_time.INFO ; ...
                            table(check.combine_newname(KK), {cmb_tempname}, ...
                            outp_time.INFO.Unit(cmb_infoIN), 'VariableNames',{'Name','LongName','Unit'})];

                    %  Add the new variable
                    outp_time.(check.combine_newname{KK}) = outp_time.(check.combine_vars{KK,1});
                    for k = 2 : length(check.combine_vars(KK,:))
                        outp_time.(check.combine_newname{KK}) = outp_time.(check.combine_newname{KK}) + ...
                            outp_time.(check.combine_vars{KK,k});
                    end

                end     % sum(ismember(cmb_fields, check.combine_vars(KK,:))) == length(check.combine_vars(KK,:))

            end     % isempty(fields(outp_time)) == 0 && ismember({check.combine_newname{KK}} , fields(outp_time)) == 0

        end     % for KK = 1 : size(check.combine_vars,1)

    end
    
    % If the run is warmstarted with a different run, add the data from the 
    % warmstart run to the output of this run
    if ismember({'lwarmstart'}, inp.namoptions.option) && ismember({'startfile'}, inp.namoptions.option)
        if strcmp(inp.namoptions.value(ismember(inp.namoptions.option, {'lwarmstart'})), '.true.')
            
            % Define the experiment number of the warmstart file
            ws_exp_nr = inp.namoptions.value(ismember(inp.namoptions.option, {'startfile'}));
            if strcmp(ws_exp_nr{1}(end), '''')
                ws_exp_nr = ws_exp_nr{1}(end-3:end-1);
            else
                ws_exp_nr = ws_exp_nr{1}(end-2:end);
            end
            
            % Define location of the .mat files of the warmstart experiment
            ws_matdir = dir([direc.main_directory direc.run_name ws_exp_nr direc.dir_separator 'output_v*']);
            ws_matdir = {ws_matdir.name}';
            ws_matdir = ws_matdir{end};
            % Define the names of the .mat files of the warmstart experiment
            ws_namelist = dir([direc.main_directory direc.run_name ws_exp_nr direc.dir_separator ...
                ws_matdir direc.dir_separator '*.mat']);
            ws_namelist = {ws_namelist.name}';
            % Check if the NAMGENSTAT output exists
            if ismember({[ws_exp_nr 'time.mat']}, ws_namelist) && ...
                    strcmp(ws_exp_nr, direc.exp_nr(1:end-1)) == 0
                
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
                fprintf(['                  - Add the NAMTIMESTAT output of run ' ...
                    ws_exp_nr ' to run ' direc.exp_nr(1:end-1) '\n']) 
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                
                % Load the input data of the warm start run
                ws_outp = load([direc.main_directory direc.run_name ws_exp_nr direc.dir_separator ...
                    ws_matdir direc.dir_separator ws_exp_nr 'time.mat']);

                ws_fields = fields(ws_outp);
                ws_fields(ismember(ws_fields,{'README','INFO'})) = [];
                % Add the warmstart output to the output o
                for i = 1 : length(ws_fields)
                    outp_time.(ws_fields{i}) = [ws_outp.(ws_fields{i}) ; outp_time.(ws_fields{i}) ];
                end

            end     % <-- if ismember({[ws_exp_nr 'time.mat']}, ws_namelist) && ...
            clear ws_*
        end     % <-- strcmp(inp.namoptions.value(ismember(inp.namoptions.option, {'lwarmstart'})), '.true.')
    end     % <-- ismember({'lwarmstart'}, inp.namoptions.option) && ismember({'startfile'}, inp.namoptions.option)
    
    

    % Save the outp_time data
    save([direc.main_directory direc.run_name direc.exp_nr ...              	% Save directory
        direc.output_dir direc.dir_separator direc.exp_nr(1:end-1) 'time.mat'], ...    % Save name
        '-struct', 'outp_time', '-v7.3');                                                  % Saved variables

    
    
end

   

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%% READ NAMCROSSSECTION basic variables  %%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function f__NAMCROSSSECTIONread_basic(direc, inp, check)
        
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
        fprintf('\n                  - Start loading basic data of the cross-section output\n') 
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    
    % Load the first xz file as a dummy
    base_files = dir([direc.main_directory, direc.run_name, direc.exp_nr, ...
        'crossxz.x*']);
    base_files = {base_files.name}';
    base_files = base_files{1};
    base_loc = [direc.main_directory, direc.run_name, direc.exp_nr, ...
        base_files];
    % Load the first xy file as a dummy
    base_filesxy = dir([direc.main_directory, direc.run_name, direc.exp_nr, ...
        'crossxy.*']);
    base_filesxy = {base_filesxy.name}';
    base_filesxy = base_filesxy{1};
    base_locxy = [direc.main_directory, direc.run_name, direc.exp_nr, ...
        base_filesxy];
%     base_varxy = {base_infoxy.Variables.Name}';
%     base_varxy = base_varxy(ismember(base_varxy, {}))
    base_infoxy = ncinfo(base_locxy);
    base_varxy = {base_infoxy.Variables}';
    base_varxy = base_varxy{1};
    base_varxy = base_varxy(ismember({base_varxy.Name}, {'yt','ym'}));
    % Save the infor from the files
    base_info = ncinfo(base_loc);
    % Define variables to be saved
    base_info = base_info.Variables;
    base_INinfo = find(ismember({base_info.Name},{'xm'}));
    base_info = [base_info(1:base_INinfo) , base_varxy , base_info(base_INinfo+1 : end)];
    base_var = {base_info.Name}';
    for i = 1 : length(base_var)
        if strcmp(base_var{i}(end-1:end), 'xz')
            base_var{i}(end-1:end) = [];
        end
    end
    % Define time & z vatiables
    base_time = ncread(base_loc, 'time') ./ 3600 + check.corr_localtime + ....
        str2num(inp.namoptions.value(ismember(inp.namoptions.option, {'xtime'})));
    base_zt = ncread(base_loc, 'zt');
    base_zm = ncread(base_loc, 'zm');
    base_zt = base_zt(base_zt <= check.z_max);
    base_zm = base_zm(1:length(base_zt));
    % Define the x variables
    base_dimx = str2num(inp.namoptions.value(ismember(inp.namoptions.option, {'itot'})));
    base_domx = str2num(inp.namoptions.value(ismember(inp.namoptions.option, {'xsize'})));
    base_dx = base_domx/base_dimx;
    base_xt = [base_dx/2 : base_dx : base_domx]';
    base_xm = [base_dx : base_dx : base_domx]';
    % Define the y variables
    base_dimy = str2num(inp.namoptions.value(ismember(inp.namoptions.option, {'jtot'})));
    base_domy = str2num(inp.namoptions.value(ismember(inp.namoptions.option, {'ysize'})));
    base_dy = base_domy/base_dimy;
    base_yt = [base_dy/2 : base_dy : base_domy]';
    base_ym = [base_dy : base_dy : base_domy]';
    % Define the positions of the cross-sections
    if ismember({'crossortho'}, inp.namoptions.option)
        base_xpos = base_xt(str2num(inp.namoptions.value(ismember(inp.namoptions.option, {'crossortho'}))));
    else
        base_xpos = base_xt(2);
    end
    if ismember({'crossortho'}, inp.namoptions.option)
        base_ypos = base_yt(str2num(inp.namoptions.value(ismember(inp.namoptions.option, {'crossplane'}))));
    else
        base_ypos = base_yt(2);
    end
    if ismember({'crossheight'}, inp.namoptions.option)
        base_zpos = base_zt(str2num(inp.namoptions.value(ismember(inp.namoptions.option, {'crossheight'}))));
    else
        base_zpos = base_zt(2);
    end
    base_zpos_fd = [];
    if ismember({'lfielddump'}, inp.namoptions.option)
        if ismember({'.true.'}, inp.namoptions.value(ismember(inp.namoptions.option,{'lfielddump'})))
            
            if ismember({'klow'}, inp.namoptions.option)
                base_zpos_fd = str2num(inp.namoptions.value(ismember(inp.namoptions.option, {'klow'})));
            else
                base_zpos_fd = 1;
            end
            if ismember({'khigh'}, inp.namoptions.option)
                base_zpos_fd = [base_zpos_fd, ...
                    str2num(inp.namoptions.value(ismember(inp.namoptions.option, {'khigh'})))];
            else
                base_zpos_fd = [base_zpos_fd, length(base_zt)];
            end
            
        end
    end
    base_zpos = [base_zpos ; base_zt(base_zpos_fd(1) : base_zpos_fd(2)) ];
    
    % Change the default scalar names into the true scalar names
    base_svnames_new = inp.scalar.Properties.VariableNames(2:end);
    base_svnames_old = cell(length(base_svnames_new),2);
    for i = 1 : length(base_svnames_new)
        % Create the name that will be used to search in f_output.(cross_fields).INFO.LongNames
        base_temp = '000';
        base_temp(end-(length(num2str(i))-1) : end) = num2str(i);

        base_svnames_old{i,1} = ['sv' base_temp];
        base_svnames_old{i,2} = ['scalar ' base_temp]; 
    end

    % Preallocate the "saving_struct" structure 
    saving_struct = struct();
    % Save an empty .mat file to be appended later
    save([direc.main_directory direc.run_name direc.exp_nr ...                  % Save directory
       direc.output_dir direc.dir_separator direc.exp_nr(1:end-1) 'cross.mat'], ...   % Save name
        '-struct', 'saving_struct', '-v7.3');                                                 % Saved variables

    % Preallocate the info on the variables
    saving_struct.INFO = table({},{},{}, 'VariableNames',{'Name','LongName','Unit'});
    % Fill the INFO field
    for i = 1 : size(base_info,2)

        if ismember(base_var(i), base_svnames_old(:,1))
            base_INsvname = find(ismember(base_svnames_old(:,1), base_var(i)));
            base_temp_longname = base_info(i).Attributes(1).Value;
            base_longname = base_temp_longname;
            for j = 1 : length(base_temp_longname) - length(base_svnames_old{base_INsvname,2}) + 1
                if strcmp(base_svnames_old{base_INsvname,2}, ...
                        base_temp_longname(j : j+length(base_svnames_old{base_INsvname,2})-1))

                    base_longname = ...
                        [base_temp_longname(1:j-1), 'scalar ', base_svnames_new{base_INsvname}, ...
                        base_temp_longname(j+length(base_svnames_old{base_INsvname,2}):end)];
                end
            end
            
            if length(base_longname) > 19
                if strcmp(base_longname(1:19), 'xz crosssection of ')
                    base_longname(1:19) = [];
                end
            end
            saving_struct.INFO = [saving_struct.INFO ; table(...
                base_svnames_new(base_INsvname), ...
                {base_longname}, {...
                'ppb'}, ...
                'VariableNames',{'Name','LongName','Unit'})];
        else
            base_longname = base_info(i).Attributes(1).Value;
            if length(base_longname) > 19
                if strcmp(base_longname(1:19), 'xz crosssection of ')
                    base_longname(1:19) = [];
                end
            end
            saving_struct.INFO = [saving_struct.INFO ; table(...
                base_var(i), {...
                base_longname}, {...
                base_info(i).Attributes(2).Value}, ...
                'VariableNames',{'Name','LongName','Unit'})];
        end

    end
    saving_struct.INFO = [saving_struct.INFO ; ...
        table({'xpos'}, {'Cross section x position'}, {'m'}, ...
            'VariableNames',{'Name','LongName','Unit'})];
    saving_struct.INFO = [saving_struct.INFO ; ...
        table({'ypos'}, {'Cross section y position'}, {'m'}, ...
            'VariableNames',{'Name','LongName','Unit'})];
    saving_struct.INFO = [saving_struct.INFO ; ...
        table({'zpos'}, {'Cross section z position'}, {'m'}, ...
            'VariableNames',{'Name','LongName','Unit'})];

    % If needed, combine the predefined variables
    if strcmp(check.combine, 'yes')
        for i = 1 : length(check.combine_newname)
            cmb_vars = check.combine_vars{i,1};
            cmb_longname = saving_struct.INFO.LongName{ismember( saving_struct.INFO.Name, cmb_vars)};
            
            for j = 1 : length(cmb_longname) - length(check.combine_vars{i,1}) + 1
                if strcmp(check.combine_vars{i,1}, cmb_longname(j : j+length(check.combine_vars{i,1})-1))
                    if j == 1
                        cmb_tempname = ...
                            [check.combine_newname{i}, ...
                            cmb_longname(j+length(check.combine_vars{i,1}):end)];
                    else
                        cmb_tempname = ...
                            [cmb_longname(1:j-1), check.combine_newname{i}, ...
                            cmb_longname(j+length(check.combine_vars{i,1}):end)];
                    end
                end
            end
            saving_struct.INFO = [ saving_struct.INFO ; ...
                table(check.combine_newname(i), {cmb_tempname}, ...
                saving_struct.INFO.Unit(ismember( saving_struct.INFO.Name, cmb_vars)), 'VariableNames',{'Name','LongName','Unit'})];

        end
        clear cmb_*
    end

    % Save the fields of the "saving" structure to the .mat file
    save([direc.main_directory direc.run_name direc.exp_nr ...                  % Save directory
        direc.output_dir direc.dir_separator direc.exp_nr(1:end-1) 'cross.mat'], ...   % Save name
        '-struct', 'saving_struct', '-append');    
    
    
    
    %%%%%
    % Recenter cross-section data around scalar emission sources 
    %%%%%

    % If needed, prepare recentering the cross-section data around the
    % scalar emission sources% Define the order in which the y-direction data will be saved 
    INnew_y = [1 : str2num( inp.namoptions.value(ismember(inp.namoptions.option,{'jtot'})) )];

    % Check if the source needs to be centered if lhetero is turned on
    if strcmp(check.center_sources, 'yes') && ismember({'lhetero'},inp.namoptions.option)   
        % Check if lhetero is turned on
        if strcmp(inp.namoptions.value(ismember(inp.namoptions.option,{'lhetero'})),'.true.')

            % Define the namoptions indices which could contain the
            % land_use(XX,YY) option
            src_INtemp = find(ismember(inp.namoptions.option,{'lhetero'}))+1;
            src_INtemp = [src_INtemp; src_INtemp-1 + find(ismember(inp.namoptions.option(src_INtemp:end),'/'))];
            src_INtemp = [src_INtemp(1) :1: src_INtemp(2)];

            % Save the YY in land_use(XX,YY), i.e. the y position of the
            % source(s)
            src_ypos = [];
            for i = 1 : length(src_INtemp)
                if length(inp.namoptions.option{src_INtemp(i)}) > 9
                    if strcmp(inp.namoptions.option{src_INtemp(i)}(1:9), 'land_use(')
                        if str2num(inp.namoptions.value{src_INtemp(i)}) > 0

                            src_IN_ypos = [find(ismember(inp.namoptions.option{src_INtemp(i)},','))+1 ...
                                       :1: ...
                                       find(ismember(inp.namoptions.option{src_INtemp(i)},')')) - 1];

                            src_ypos = [src_ypos ; str2num(inp.namoptions.option{src_INtemp(i)}(src_IN_ypos))];
                        end

                    end
                end
            end
            if isempty(src_ypos)
                src_ypos = 1;
            end

            % Define the domain gridpoints in y-direction
            src_ygrid = [str2num( inp.namoptions.value(ismember(inp.namoptions.option,{'jtot'})) )];
            % Define the number of patches in y-direction
            src_ypatches = [str2num( inp.namoptions.value(ismember(inp.namoptions.option,{'ypatches'})) )];

            if src_ygrid ~= src_ypatches
                sec_temp_ypos = src_ypos;
                src_ypos = [];
                for i = 1 : length(sec_temp_ypos)
                    src_ypos = [src_ypos, [(src_ygrid / src_ypatches * (sec_temp_ypos(i) -1) +1) ...
                                :1: ...
                                (src_ygrid / src_ypatches * sec_temp_ypos(i))] ];
                end
            end
            src_ypos = floor(mean([min(src_ypos), max(src_ypos)]));
            src_change = floor(src_ygrid/2) - src_ypos;
            if src_change > 0
                INnew_y = [INnew_y(end - src_change+1:end), INnew_y(1 : end-src_change)];
            elseif src_change < 0
                INnew_y = [INnew_y(-1*src_change+1:end), INnew_y(1 : -1*src_change)];
            end
        end
    end

    clear src_*
    
    % correct the y position of the cross-section after centering
    base_ypos = base_yt(round(base_yt(INnew_y),6) == round(base_ypos,6));
    
    %%%%%
    % Create a xy cross-section mask with the source position(s)
    %%%%%
    
    % Define the x boundaries of the surface hetero patches 
    src_xpatches = str2num(inp.namoptions.value( ...
        ismember(inp.namoptions.option,  {'xpatches'}))  );
    src_INx = [1 : ...
              length(base_xt) / src_xpatches : ...
              length(base_xt)];
    src_emis_xstart = base_xt(src_INx);
    src_INx = [length(base_xt) / src_xpatches : ...
              length(base_xt) / src_xpatches : ...
              length(base_xt)];
    src_emis_xend = base_xt(src_INx);
    % Define the y boundaries of the surface hetero patches 
    src_ypatches = str2num(inp.namoptions.value( ...
        ismember(inp.namoptions.option,  {'ypatches'}))  );
    src_INy = [1 : ...
              length(base_yt) / src_ypatches : ...
              length(base_yt)];
    src_emis_ystart = base_yt(src_INy);
    src_INy = [length(base_yt) / src_ypatches : ...
              length(base_yt) / src_ypatches : ...
              length(base_yt)];
    src_emis_yend = base_yt(src_INy);
        
    % Define the location of the emission sources
    src_INlanduse= [];   % Index on where in namoptions landuse can be found
    src_landusescheme = [];
    src_INemis = [];
    for j = 1 : length(inp.namoptions.option)
        src_temp_inp_opt = inp.namoptions.option{j};
        src_temp_inp_val = str2num(inp.namoptions.value{j});
        if length(src_temp_inp_opt) >= 13
            if strcmp(src_temp_inp_opt(1:9), 'land_use(') 
                src_INlanduse = [src_INlanduse; j];
                src_landusescheme = [src_landusescheme; src_temp_inp_val];


                src_INemis = [src_INemis; ...
                    str2num(src_temp_inp_opt(10 : find(ismember(src_temp_inp_opt, ','))-1)), ...
                    str2num(src_temp_inp_opt(find(ismember(src_temp_inp_opt, ','))+1 : end-1))];

            end
        end
    end



    src_INx_source = cell(length(src_INlanduse),1);
    src_INy_source = cell(length(src_INlanduse),1);
    for j = 1 : length(src_landusescheme)
        % Find indices of x-coordinates of the sources and the NaN band
        % around it
        src_INx_source{j} = find( ...
            base_xt >= src_emis_xstart(src_INemis(j,1)) & ...
            base_xt <= src_emis_xend(src_INemis(j,1)));
        % Find indices of y-coordinates of the sources and the NaN band
        % around it
        src_INy_source{j} = find( ...
            base_yt >= src_emis_ystart(src_INemis(j,2)) & ...
            base_yt <= src_emis_yend(src_INemis(j,2)) );
    end


    base_srcpos = struct();
    base_scalars = check.plume_receptor_source(:,1);
    for i = 1 : length(base_scalars)
        base_srcpos.(base_scalars{i}) = zeros( length(base_xt) , length(base_yt) ); 
        
        for j = 1 : length(src_landusescheme)
        
            if ismember(src_landusescheme(j), check.plume_receptor_source{i,2})
                base_srcpos.(base_scalars{i})(src_INx_source{j}, src_INy_source{j}) = 1;
            end

        end     % <-- for j = 1 : length(src_landusescheme)
        
        base_srcpos.(base_scalars{i}) = base_srcpos.(base_scalars{i})(:, INnew_y);
        
    end     % <-- for i = 1 : length(base_scalars)

    
    
    % reset the saving_struct structure
    saving_struct = struct();
    % Add the base variables to the saving_struct structure
    saving_struct.time = base_time;
    saving_struct.xt = base_xt;
    saving_struct.xm = base_xm;    
    saving_struct.yt = base_yt;
    saving_struct.ym = base_ym;
    saving_struct.zt = double(base_zt);
    saving_struct.zm = double(base_zm);
    saving_struct.xpos = base_xpos;
    saving_struct.ypos = base_ypos;
    saving_struct.zpos = double(base_zpos);
    saving_struct.cross_lvl = check.cross_lvl;
    saving_struct.srcpos = base_srcpos;
    % Save the fields of the "saving" structure to the .mat file
    save([direc.main_directory direc.run_name direc.exp_nr ...                  % Save directory
        direc.output_dir direc.dir_separator direc.exp_nr(1:end-1) 'cross.mat'], ...   % Save name
        '-struct', 'saving_struct', '-append');    
    
    % Clear the data and reset the saving_struct
    clear base_* saving_struct

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
        fprintf('                  - Basic data of the cross-section output is processed and saved!\n') 
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
end     % <-- End of function
        
        

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%  READ XY NAMCROSSSECTION output   %%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function f__NAMCROSSSECTIONread_xy(direc, inp, check)
    
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
        fprintf('\n                  - Start loading XY cross-section\n') 
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    
    
    
    % Find output files with an XY cross-section
    xy_files = dir([direc.main_directory, direc.run_name, direc.exp_nr, ...
        'crossxy.*']);
    xy_files = {xy_files.name}';
    
    if isempty(xy_files) == 0
        
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %    Recenter cross-section data around scalar emission sources   %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        % If needed, prepare recentering the cross-section data around the
        % scalar emission sources% Define the order in which the y-direction data will be saved 
        INnew_y = [1 : str2num( inp.namoptions.value(ismember(inp.namoptions.option,{'jtot'})) )];

        % Check if the source needs to be centered if lhetero is turned on
        if strcmp(check.center_sources, 'yes') && ismember({'lhetero'},inp.namoptions.option)   
            % Check if lhetero is turned on
            if strcmp(inp.namoptions.value(ismember(inp.namoptions.option,{'lhetero'})),'.true.')

                % Define the namoptions indices which could contain the
                % land_use(XX,YY) option
                src_INtemp = find(ismember(inp.namoptions.option,{'lhetero'}))+1;
                src_INtemp = [src_INtemp; src_INtemp-1 + find(ismember(inp.namoptions.option(src_INtemp:end),'/'))];
                src_INtemp = [src_INtemp(1) :1: src_INtemp(2)];

                % Save the YY in land_use(XX,YY), i.e. the y position of the
                % source(s)
                src_ypos = [];
                for i = 1 : length(src_INtemp)
                    if length(inp.namoptions.option{src_INtemp(i)}) > 9
                        if strcmp(inp.namoptions.option{src_INtemp(i)}(1:9), 'land_use(')
                            if str2num(inp.namoptions.value{src_INtemp(i)}) > 0

                                src_IN_ypos = [find(ismember(inp.namoptions.option{src_INtemp(i)},','))+1 ...
                                           :1: ...
                                           find(ismember(inp.namoptions.option{src_INtemp(i)},')')) - 1];

                                src_ypos = [src_ypos ; str2num(inp.namoptions.option{src_INtemp(i)}(src_IN_ypos))];
                            end
                            
                        end
                    end
                end
                if isempty(src_ypos)
                    src_ypos = 1;
                end

                % Define the domain gridpoints in y-direction
                src_ygrid = [str2num( inp.namoptions.value(ismember(inp.namoptions.option,{'jtot'})) )];
                % Define the number of patches in y-direction
                src_ypatches = [str2num( inp.namoptions.value(ismember(inp.namoptions.option,{'ypatches'})) )];

                if src_ygrid ~= src_ypatches
                    sec_temp_ypos = src_ypos;
                    src_ypos = [];
                    for i = 1 : length(sec_temp_ypos)
                        src_ypos = [src_ypos, [(src_ygrid / src_ypatches * (sec_temp_ypos(i) -1) +1) ...
                                    :1: ...
                                    (src_ygrid / src_ypatches * sec_temp_ypos(i))] ];
                    end
                end
                src_ypos = floor(mean([min(src_ypos), max(src_ypos)]));
                src_change = floor(src_ygrid/2) - src_ypos;
                if src_change > 0
                    INnew_y = [INnew_y(end - src_change+1:end), INnew_y(1 : end-src_change)];
                elseif src_change < 0
                    INnew_y = [INnew_y(-1*src_change+1:end), INnew_y(1 : -1*src_change)];
                end
            end
        end
        
        clear src_*
        
        
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %       Setup the information on the XY CROSSSECTION output       %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
       
        % Load the first file as a dummy
        xy_loc = [direc.main_directory, direc.run_name, direc.exp_nr, ...
            xy_files{1}];
        xy_info = ncinfo(xy_loc);
        % Define variables to be saved
        xy_var = {xy_info.Variables.Name}';
        % Preallocate time variable size
        xy_dimt = xy_info.Variables(ismember({xy_info.Variables.Name}, {'time'})).Size;
        % Preallocate x variable size
        xy_dimx = xy_info.Variables(ismember({xy_info.Variables.Name},{'xt'})).Size(1);
        xy_nprocx = str2num(inp.namoptions.value(ismember(inp.namoptions.option, 'nprocx')));
        % Preallocate y variable size
        xy_dimy = xy_info.Variables(ismember({xy_info.Variables.Name},{'yt'})).Size(1);
        xy_nprocy = str2num(inp.namoptions.value(ismember(inp.namoptions.option, 'nprocy')));
        % Determine the vertical index of the xy cross-section
        
        
        
        % Change the default scalar names into the true scalar names
        xy_svnames_new = inp.scalar.Properties.VariableNames(2:end);
        xy_svnames_old = cell(length(xy_svnames_new),2);
        for i = 1 : length(xy_svnames_new)
            
            % Add the cross-section indicator (xz) to the scalar name
            xy_svnames_new{i} = [xy_svnames_new{i}];
            
            % Create the name that will be used to search in f_output.(cross_fields).INFO.LongNames
            xy_temp = '000';
            xy_temp(end-(length(num2str(i))-1) : end) = num2str(i);
            
            xy_svnames_old{i,1} = ['sv' xy_temp 'xy'];
            xy_svnames_old{i,2} = ['scalar ' xy_temp]; 
        end
        
        % remove redundant variables (these are already saved in
        % cross.mat
        xy_var(ismember(xy_var, {'time','xt','xm','yt','ym','zt','zm'})) = [];
        tmp_processvar = check.process_var;
        for i = 1 : length(xy_svnames_new)
            if ismember(xy_svnames_new(i), tmp_processvar)
                tmp_processvar{ismember(tmp_processvar, xy_svnames_new(i))} = xy_svnames_old{i}(1:end-2);
            end
        end
        for i = 1 : length(tmp_processvar)
            tmp_processvar{i} = [tmp_processvar{i} 'xy'];
        end
        xy_var(ismember(xy_var, tmp_processvar)==0) = [];
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %                       Save each variable                        %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        for i =  1 : length(xy_var)
            
            % Reset the "saving_struct"
            saving_struct = struct();

            % Define the name of the variable that will be saved to the
            % .mat file. Here the scalar names in scalar.inp are taken
            % into account
            if ismember(xy_var(i), xy_svnames_old(:,1))
                xy_savevar = xy_svnames_new{ismember(xy_svnames_old(:,1), xy_var(i))};
            else
                xy_savevar = xy_var{i}(1:end-2);
            end

            saving_struct.(xy_savevar) = NaN(xy_dimx*xy_nprocx, xy_dimy*xy_nprocy, xy_dimt);
            for j = 1 : length(xy_files)
                % Determine the path of the to be loaded netCDF file 
                xy_loc = [direc.main_directory, direc.run_name, direc.exp_nr, xy_files{j}];
                % Define the core number that is read
                xy_locx = str2num(xy_files{j}(15:17));
                xy_locy = str2num(xy_files{j}(19:21));
                % Determine the indices that are filled by this netCDF
                % file
                xy_fillx = [1+xy_dimx*xy_locx : 1+xy_dimx*(xy_locx+1)-1];
                xy_filly = [1+xy_dimy*xy_locy : 1+xy_dimy*(xy_locy+1)-1];
                % Redefine the y-position of the data to match "INnew_y"
                xy_filly = find(ismember(INnew_y, xy_filly));
                
                % Load the data from the netCDF file
                saving_struct.(xy_savevar)(xy_fillx, xy_filly, :) =  ncread(xy_loc, xy_var{i});
            end
            
            
            % Save the variable to the .mat file
            save([direc.main_directory direc.run_name direc.exp_nr ...                  % Save directory
                direc.output_dir direc.dir_separator direc.exp_nr(1:end-1) ...
                'crossxy' check.cross_lvl{1} '_' xy_savevar '.mat'], ...   % Save name
                '-struct', 'saving_struct', '-v7.3'); 

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
            fprintf(['                    Variable "' xy_savevar '" is loaded and saved.\n']) 
            fprintf(['                    - XY cross-section ' check.cross_lvl{1} ' at ' num2str(i/length(xy_var)*100) '%%.\n']) 
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            % clear the data
            clear saving_struct xy_loc xy_filly_new xy_filly j xy_savevar
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %                  Combine variables (if needed)                  %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
         % If needed, combine the predefined variables
        if strcmp(check.combine, 'yes')
    
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
            fprintf('                  - Combine variables of the XY coss-section!\n') 
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            cmb_newname = check.combine_newname;
            cmb_newname = cmb_newname(ismember(cmb_newname, check.process_var));
            cmb_vars = check.combine_vars(find( ismember(check.combine_newname, cmb_newname) ),:);
            cmb_files = dir([direc.main_directory direc.run_name direc.exp_nr ... 
                        direc.output_dir direc.dir_separator]);
            cmb_files = {cmb_files.name}';
            
            % Combine the variables
            for i = 1 : length(cmb_newname)
                
                % Define the processed file names of the to-be-combined
                % variables
                cmb_varsfiles = cmb_vars(i,:);
                for j = 1 : length(cmb_varsfiles)
                    cmb_varsfiles{j} = [direc.exp_nr(1:end-1) ...
                        'crossxy' check.cross_lvl{1} '_' cmb_varsfiles{j} '.mat'];
                end
                
                % Make sure the to-be-combined vairables are already
                % processed
                if sum( ismember(cmb_varsfiles, cmb_files) ) == length(cmb_varsfiles)

                    % Reset saving_struct
                    saving_struct = struct();

                    % Add the first of the to be combined variables
                    cmb_matfile = matfile([direc.main_directory direc.run_name direc.exp_nr ... 
                        direc.output_dir direc.dir_separator direc.exp_nr(1:end-1) ...
                        'crossxy' check.cross_lvl{1} '_' cmb_vars{i,1} '.mat']);
                    saving_struct.(cmb_newname{i}) = cmb_matfile.(cmb_vars{i,1});

                    % Add the remaining to be combined variables
                    for j = 2 : size(cmb_vars,2)

                        cmb_matfile = matfile([direc.main_directory direc.run_name direc.exp_nr ... 
                            direc.output_dir direc.dir_separator ...
                            direc.exp_nr(1:end-1) 'crossxy' check.cross_lvl{1} '_' cmb_vars{i,j} '.mat']);

                        saving_struct.(cmb_newname{i}) = saving_struct.(cmb_newname{i}) + ...
                             cmb_matfile.(cmb_vars{i,j});
                    end

                    % Save the combined variable to the .mat file
                    save([direc.main_directory direc.run_name direc.exp_nr ...                  % Save directory
                        direc.output_dir direc.dir_separator direc.exp_nr(1:end-1) ...
                        'crossxy' check.cross_lvl{1} '_' cmb_newname{i} '.mat'], ...   % Save name
                        '-struct', 'saving_struct', '-v7.3'); 
                   % clear the data
                   clear saving_struct j cmb_matfile

                else
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    fprintf(['                    - !!! ERROR !!! \n'...
                        '                      The to-be-combined variables are not yet processed\n']) 
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                
                end         % <-- if sum( ismember(cmb_varsfiles, cmb_files) ) == length(cmb_varsfiles)
            
            end         % for i = 1 : length(cmb_newname)
                
        end         % <-- if strcmp(check.combine, 'yes')
        
    end         % <-- if isempty(xy_files) == 0
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
        fprintf('                  - XY coss-section output is processed and saved!\n') 
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    
end     % <-- End of function
        
        

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%  READ XZ NAMCROSSSECTION output   %%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function f__NAMCROSSSECTIONread_xz(direc, inp, check)
    
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
        fprintf('\n                  - Start loading XZ cross-section\n') 
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    
    
    
    % Find output files with an XZ cross-section
    xz_files = dir([direc.main_directory, direc.run_name, direc.exp_nr, ...
        'crossxz.x*']);
    xz_files = {xz_files.name}';
    
    if isempty(xz_files) == 0
        
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %    Recenter cross-section data around scalar emission sources   %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        % If needed, prepare recentering the cross-section data around the
        % scalar emission sources% Define the order in which the y-direction data will be saved 
        INnew_y = [1 : str2num( inp.namoptions.value(ismember(inp.namoptions.option,{'jtot'})) )];

        % Check if the source needs to be centered if lhetero is turned on
        if strcmp(check.center_sources, 'yes') && ismember({'lhetero'},inp.namoptions.option)   
            % Check if lhetero is turned on
            if strcmp(inp.namoptions.value(ismember(inp.namoptions.option,{'lhetero'})),'.true.')

                % Define the namoptions indices which could contain the
                % land_use(XX,YY) option
                src_INtemp = find(ismember(inp.namoptions.option,{'lhetero'}))+1;
                src_INtemp = [src_INtemp; src_INtemp-1 + find(ismember(inp.namoptions.option(src_INtemp:end),'/'))];
                src_INtemp = [src_INtemp(1) :1: src_INtemp(2)];

                % Save the YY in land_use(XX,YY), i.e. the y position of the
                % source(s)
                src_ypos = [];
                for i = 1 : length(src_INtemp)
                    if length(inp.namoptions.option{src_INtemp(i)}) > 9
                        if strcmp(inp.namoptions.option{src_INtemp(i)}(1:9), 'land_use(')
                            if str2num(inp.namoptions.value{src_INtemp(i)}) > 0

                                src_IN_ypos = [find(ismember(inp.namoptions.option{src_INtemp(i)},','))+1 ...
                                           :1: ...
                                           find(ismember(inp.namoptions.option{src_INtemp(i)},')')) - 1];

                                src_ypos = [src_ypos ; str2num(inp.namoptions.option{src_INtemp(i)}(src_IN_ypos))];
                            end
                            
                        end
                    end
                end
                if isempty(src_ypos)
                    src_ypos = 1;
                end

                % Define the domain gridpoints in y-direction
                src_ygrid = [str2num( inp.namoptions.value(ismember(inp.namoptions.option,{'jtot'})) )];
                % Define the number of patches in y-direction
                src_ypatches = [str2num( inp.namoptions.value(ismember(inp.namoptions.option,{'ypatches'})) )];

                if src_ygrid ~= src_ypatches
                    sec_temp_ypos = src_ypos;
                    src_ypos = [];
                    for i = 1 : length(sec_temp_ypos)
                        src_ypos = [src_ypos, [(src_ygrid / src_ypatches * (sec_temp_ypos(i) -1) +1) ...
                                    :1: ...
                                    (src_ygrid / src_ypatches * sec_temp_ypos(i))] ];
                    end
                end
                src_ypos = floor(mean([min(src_ypos), max(src_ypos)]));
                src_change = floor(src_ygrid/2) - src_ypos;
                if src_change > 0
                    INnew_y = [INnew_y(end - src_change+1:end), INnew_y(1 : end-src_change)];
                elseif src_change < 0
                    INnew_y = [INnew_y(-1*src_change+1:end), INnew_y(1 : -1*src_change)];
                end
            end
        end
        
        clear src_*
        
        
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %       Setup the information on the XZ CROSSSECTION output       %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
       
        % Load the first file as a dummy
        xz_loc = [direc.main_directory, direc.run_name, direc.exp_nr, ...
            xz_files{1}];
        xz_info = ncinfo(xz_loc);
        % Define variables to be saved
        xz_var = {xz_info.Variables.Name}';
        % Preallocate time variable size
        xz_dimt = xz_info.Variables(ismember({xz_info.Variables.Name}, {'time'})).Size;
        % Preallocate x variable size
        xz_dimx = xz_info.Variables(ismember({xz_info.Variables.Name},{'xt'})).Size(1);
        xz_nprocx = str2num(inp.namoptions.value(ismember(inp.namoptions.option, 'nprocx')));
        % Preallocate z variable size
        xz_zIN = find(inp.prof.zc <= check.z_max);
        xz_dimz = max(xz_zIN);
        
        % Change the default scalar names into the true scalar names
        xz_svnames_new = inp.scalar.Properties.VariableNames(2:end);
        xz_svnames_old = cell(length(xz_svnames_new),2);
        for i = 1 : length(xz_svnames_new)
            
            % Add the cross-section indicator (xz) to the scalar name
            xz_svnames_new{i} = [xz_svnames_new{i}];
            
            % Create the name that will be used to search in f_output.(cross_fields).INFO.LongNames
            xz_temp = '000';
            xz_temp(end-(length(num2str(i))-1) : end) = num2str(i);
            
            xz_svnames_old{i,1} = ['sv' xz_temp 'xz'];
            xz_svnames_old{i,2} = ['scalar ' xz_temp]; 
        end
        
        % remove redundant variables (these are already saved in
        % cross.mat
        xz_var(ismember(xz_var, {'time','xt','xm','yt','ym','zt','zm'})) = [];
        tmp_processvar = check.process_var;
        for i = 1 : length(xz_svnames_new)
            if ismember(xz_svnames_new(i), tmp_processvar)
                tmp_processvar{ismember(tmp_processvar, xz_svnames_new(i))} = xz_svnames_old{i}(1:end-2);
            end
        end
        for i = 1 : length(tmp_processvar)
            tmp_processvar{i} = [tmp_processvar{i} 'xz'];
        end
        xz_var(ismember(xz_var, tmp_processvar)==0)= [];
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %                       Save each variable                        %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        for i =  1 : length(xz_var)
            % Reset the "saving_struct"
            saving_struct = struct();

            % Define the name of the variable that will be saved to the
            % .mat file. Here the scalar names in scalar.inp are taken
            % into account
            if ismember(xz_var(i), xz_svnames_old(:,1))
                xz_savevar = xz_svnames_new{ismember(xz_svnames_old(:,1), xz_var(i))};
            else
                xz_savevar = xz_var{i}(1:end-2);
            end

            saving_struct.(xz_savevar) = NaN(xz_dimx*xz_nprocx, xz_dimz, xz_dimt);
            for j = 1 : length(xz_files)
                % Determine the path of the to be loaded netCDF file 
                xz_loc = [direc.main_directory, direc.run_name, direc.exp_nr, xz_files{j}];
                % Determine the indices that are filled by this netCDF
                % file
                xz_fillx = [1+xz_dimx*str2num(xz_files{j}(10:12)) :  1 : ...
                    1+xz_dimx*(str2num(xz_files{j}(10:12))+1)-1];
                
                saving_struct.(xz_savevar)(xz_fillx,: ,: ) = ncread(xz_loc, xz_var{i}, [1, 1, 1], [Inf, max(xz_zIN), Inf]);
            end

            % Save the variable to the .mat file
            save([direc.main_directory direc.run_name direc.exp_nr ...                  % Save directory
                direc.output_dir direc.dir_separator direc.exp_nr(1:end-1) 'crossxz_' xz_savevar '.mat'], ...   % Save name
                '-struct', 'saving_struct', '-v7.3'); 
            % clear the data

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
            fprintf(['                    Variable "' xz_savevar '" is loaded and saved.\n']) 
            fprintf(['                    - XZ cross-section at ' num2str(i/length(xz_var)*100) '%%.\n']) 
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            clear saving_struct xz_savevar xz_loc xz_fillx j 
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %                  Combine variables (if needed)                  %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
         % If needed, combine the predefined variables
        if strcmp(check.combine, 'yes')
             
            cmb_newname = check.combine_newname;
            cmb_newname = cmb_newname(ismember(cmb_newname, check.process_var));
            cmb_vars = check.combine_vars(find( ismember(check.combine_newname, cmb_newname) ),:);
            cmb_files = dir([direc.main_directory direc.run_name direc.exp_nr ... 
                        direc.output_dir direc.dir_separator]);
            cmb_files = {cmb_files.name}';
                
                
            % Combine the variables
            for i = 1 : length(cmb_newname)
                    
                % Define the processed file names of the to-be-combined
                % variables
                cmb_varsfiles = cmb_vars(i,:);
                for j = 1 : length(cmb_varsfiles)
                    cmb_varsfiles{j} = [direc.exp_nr(1:end-1) 'crossxz_' cmb_varsfiles{j} '.mat'];
                end
                
                % Make sure the to-be-combined vairables are already
                % processed
                if sum( ismember(cmb_varsfiles, cmb_files) ) == length(cmb_varsfiles)
                    
                    
                    
                    % Reset saving_struct
                    saving_struct = struct();
    
                    % Add the first of the to be combined variables
                    cmb_matfile = matfile([direc.main_directory direc.run_name direc.exp_nr ... 
                        direc.output_dir direc.dir_separator ...
                        direc.exp_nr(1:end-1) 'crossxz_' cmb_vars{i,1} '.mat']);
                    saving_struct.(cmb_newname{i}) = cmb_matfile.(cmb_vars{i,1});
                    
                    % Add the remaining to be combined variables
                    for j = 2 : size(cmb_vars,2)
                        
                        cmb_matfile = matfile([direc.main_directory direc.run_name direc.exp_nr ... 
                            direc.output_dir direc.dir_separator ...
                            direc.exp_nr(1:end-1) 'crossxz_' cmb_vars{i,j} '.mat']);
                        
                        saving_struct.(cmb_newname{i}) = saving_struct.(cmb_newname{i}) + ...
                             cmb_matfile.(cmb_vars{i,j});
                    end

                    % Save the combined variable to the .mat file
                    save([direc.main_directory direc.run_name direc.exp_nr ...                  % Save directory
                        direc.output_dir direc.dir_separator direc.exp_nr(1:end-1) 'crossxz_' cmb_newname{i} '.mat'], ...   % Save name
                        '-struct', 'saving_struct', '-v7.3'); 
                   % clear the data
                   clear saving_struct cmb_matfile j 
                
                else
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    fprintf(['                    - !!! ERROR !!! \n'...
                        '                      The to-be-combined variables are not yet processed\n']) 
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

                end         % <-- if sum( ismember(cmb_varsfiles, cmb_files) ) == length(cmb_varsfiles)
            
            end         % <-- for i = 1 : length(check.combine_newname)
            
        end         % <-- if strcmp(check.combine, 'yes')
        
    end         % <-- if isempty(xz_files) == 0
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
        fprintf('                  - XZ coss-section output is processed and saved!\n') 
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    
end     % <-- End of function
        
        

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%  READ YZ NAMCROSSSECTION output   %%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function f__NAMCROSSSECTIONread_yz(direc, inp, check)
    
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
        fprintf('\n                  - Start loading YZ cross-section\n') 
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    
    
    
    % Find output files with an YZ cross-section
    yz_files = dir([direc.main_directory, direc.run_name, direc.exp_nr, ...
        'crossyz.x*']);
    yz_files = {yz_files.name}';
    
    if isempty(yz_files) == 0
        
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %    Recenter cross-section data around scalar emission sources   %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        % If needed, prepare recentering the cross-section data around the
        % scalar emission sources% Define the order in which the y-direction data will be saved 
        INnew_y = [1 : str2num( inp.namoptions.value(ismember(inp.namoptions.option,{'jtot'})) )];

        % Check if the source needs to be centered if lhetero is turned on
        if strcmp(check.center_sources, 'yes') && ismember({'lhetero'},inp.namoptions.option)   
            % Check if lhetero is turned on
            if strcmp(inp.namoptions.value(ismember(inp.namoptions.option,{'lhetero'})),'.true.')

                % Define the namoptions indices which could contain the
                % land_use(XX,YY) option
                src_INtemp = find(ismember(inp.namoptions.option,{'lhetero'}))+1;
                src_INtemp = [src_INtemp; src_INtemp-1 + find(ismember(inp.namoptions.option(src_INtemp:end),'/'))];
                src_INtemp = [src_INtemp(1) :1: src_INtemp(2)];

                % Save the YY in land_use(XX,YY), i.e. the y position of the
                % source(s)
                src_ypos = [];
                for i = 1 : length(src_INtemp)
                    if length(inp.namoptions.option{src_INtemp(i)}) > 9
                        if strcmp(inp.namoptions.option{src_INtemp(i)}(1:9), 'land_use(')
                            if str2num(inp.namoptions.value{src_INtemp(i)}) > 0

                                src_IN_ypos = [find(ismember(inp.namoptions.option{src_INtemp(i)},','))+1 ...
                                           :1: ...
                                           find(ismember(inp.namoptions.option{src_INtemp(i)},')')) - 1];

                                src_ypos = [src_ypos ; str2num(inp.namoptions.option{src_INtemp(i)}(src_IN_ypos))];
                            end
                            
                        end
                    end
                end
                if isempty(src_ypos)
                    src_ypos = 1;
                end

                % Define the domain gridpoints in y-direction
                src_ygrid = [str2num( inp.namoptions.value(ismember(inp.namoptions.option,{'jtot'})) )];
                % Define the number of patches in y-direction
                src_ypatches = [str2num( inp.namoptions.value(ismember(inp.namoptions.option,{'ypatches'})) )];

                if src_ygrid ~= src_ypatches
                    sec_temp_ypos = src_ypos;
                    src_ypos = [];
                    for i = 1 : length(sec_temp_ypos)
                        src_ypos = [src_ypos, [(src_ygrid / src_ypatches * (sec_temp_ypos(i) -1) +1) ...
                                    :1: ...
                                    (src_ygrid / src_ypatches * sec_temp_ypos(i))] ];
                    end
                end
                src_ypos = floor(mean([min(src_ypos), max(src_ypos)]));
                src_change = floor(src_ygrid/2) - src_ypos;
                if src_change > 0
                    INnew_y = [INnew_y(end - src_change+1:end), INnew_y(1 : end-src_change)];
                elseif src_change < 0
                    INnew_y = [INnew_y(-1*src_change+1:end), INnew_y(1 : -1*src_change)];
                end
            end
        end
        
        clear src_*
        
        
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %       Setup the information on the YZ CROSSSECTION output       %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
       
        % Load the first file as a dummy
        yz_loc = [direc.main_directory, direc.run_name, direc.exp_nr, ...
            yz_files{1}];
        yz_info = ncinfo(yz_loc);
        % Define variables to be saved
        yz_var = {yz_info.Variables.Name}';
        % Preallocate time variable size
        yz_dimt = yz_info.Variables(ismember({yz_info.Variables.Name}, {'time'})).Size;
        % Preallocate y variable size
        yz_dimy = yz_info.Variables(ismember({yz_info.Variables.Name},{'yt'})).Size(1);
        yz_nprocy = str2num(inp.namoptions.value(ismember(inp.namoptions.option, 'nprocy')));
        % Preallocate z variable size
        yz_zIN = find(inp.prof.zc <= check.z_max);
        yz_dimz = max(yz_zIN);
        
        % Change the default scalar names into the true scalar names
        yz_svnames_new = inp.scalar.Properties.VariableNames(2:end);
        yz_svnames_old = cell(length(yz_svnames_new),2);
        for i = 1 : length(yz_svnames_new)
            
            % Add the cross-section indicator (xz) to the scalar name
            yz_svnames_new{i} = [yz_svnames_new{i}];
            
            % Create the name that will be used to search in f_output.(cross_fields).INFO.LongNames
            yz_temp = '000';
            yz_temp(end-(length(num2str(i))-1) : end) = num2str(i);
            
            yz_svnames_old{i,1} = ['sv' yz_temp 'yz'];
            yz_svnames_old{i,2} = ['scalar ' yz_temp]; 
        end
        
        % remove redundant variables (these are already saved in
        % cross.mat
        yz_var(ismember(yz_var, {'time','xt','xm','yt','ym','zt','zm'})) = [];
        tmp_processvar = check.process_var;
        for i = 1 : length(yz_svnames_new)
            if ismember(yz_svnames_new(i), tmp_processvar)
                tmp_processvar{ismember(tmp_processvar, yz_svnames_new(i))} = yz_svnames_old{i}(1:end-2);
            end
        end
        for i = 1 : length(tmp_processvar)
            tmp_processvar{i} = [tmp_processvar{i} 'yz'];
        end
        yz_var(ismember(yz_var, tmp_processvar)==0) = [];
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %                       Save each variable                        %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        for i =  1 : length(yz_var)
            % Reset the "saving_struct"
            saving_struct = struct();

            % Define the name of the variable that will be saved to the
            % .mat file. Here the scalar names in scalar.inp are taken
            % into account
            if ismember(yz_var(i), yz_svnames_old(:,1))
                yz_savevar = yz_svnames_new{ismember(yz_svnames_old(:,1), yz_var(i))};
            else
                yz_savevar = yz_var{i}(1:end-2);
            end

            saving_struct.(yz_savevar) = NaN(yz_dimy*yz_nprocy, yz_dimz, yz_dimt);
            for j = 1 : length(yz_files)
                % Determine the path of the to be loaded netCDF file 
                yz_loc = [direc.main_directory, direc.run_name, direc.exp_nr, yz_files{j}];
                % Define the core number that is read
                yz_locy = str2num(yz_files{j}(14:16));
                % Determine the indices that are filled by this netCDF
                % file
                yz_filly = [1+yz_dimy*yz_locy : 1+yz_dimy*(yz_locy+1)-1];
                % Redefine the y-position of the data to match "INnew_y"
                yz_filly_new = find(ismember(INnew_y, yz_filly));
                
                % Load the data from the netCDF file
                saving_struct.(yz_savevar)(yz_filly_new,: ,: ) = ...
                    ncread(yz_loc, yz_var{i}, [1, 1, 1], [Inf, max(yz_zIN), Inf]);
            end

            % Save the variable to the .mat file
            save([direc.main_directory direc.run_name direc.exp_nr ...                  % Save directory
                direc.output_dir direc.dir_separator direc.exp_nr(1:end-1) 'crossyz_' yz_savevar '.mat'], ...   % Save name
                '-struct', 'saving_struct', '-v7.3'); 
            % clear the data

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
            fprintf(['                    Variable "' yz_savevar '" is loaded and saved.\n']) 
            fprintf(['                    - YZ cross-section at ' num2str(i/length(yz_var)*100) '%%.\n']) 
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            clear saving_struct yz_savevar yz_loc yz_filly_new yz_filly j 
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %                  Combine variables (if needed)                  %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
         % If needed, combine the predefined variables
        if strcmp(check.combine, 'yes')
    
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
            fprintf('                  - Combine variables of the YZ coss-section!\n') 
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            cmb_newname = check.combine_newname;
            cmb_newname = cmb_newname(ismember(cmb_newname, check.process_var));
            cmb_vars = check.combine_vars(find( ismember(check.combine_newname, cmb_newname) ),:);
            cmb_files = dir([direc.main_directory direc.run_name direc.exp_nr ... 
                        direc.output_dir direc.dir_separator]);
            cmb_files = {cmb_files.name}';
                
                
            % Combine the variables
            for i = 1 : length(cmb_newname)
                    
                % Define the processed file names of the to-be-combined
                % variables
                cmb_varsfiles = cmb_vars(i,:);
                for j = 1 : length(cmb_varsfiles)
                    cmb_varsfiles{j} = [direc.exp_nr(1:end-1) 'crossyz_' cmb_varsfiles{j} '.mat'];
                end
                
                % Make sure the to-be-combined vairables are already
                % processed
                if sum( ismember(cmb_varsfiles, cmb_files) ) == length(cmb_varsfiles)

                    % Reset saving_struct
                    saving_struct = struct();
    
                    % Add the first of the to be combined variables
                    cmb_matfile = matfile([direc.main_directory direc.run_name direc.exp_nr ... 
                        direc.output_dir direc.dir_separator ...
                        direc.exp_nr(1:end-1) 'crossyz_' cmb_vars{i,1} '.mat']);
                    saving_struct.(cmb_newname{i}) = cmb_matfile.(cmb_vars{i,1});
                    
                    % Add the remaining to be combined variables
                    for j = 2 : size(cmb_vars,2)
                        
                        cmb_matfile = matfile([direc.main_directory direc.run_name direc.exp_nr ... 
                            direc.output_dir direc.dir_separator ...
                            direc.exp_nr(1:end-1) 'crossyz_' cmb_vars{i,j} '.mat']);
                        
                        saving_struct.(cmb_newname{i}) = saving_struct.(cmb_newname{i}) + ...
                             cmb_matfile.(cmb_vars{i,j});
                    end

                    % Save the combined variable to the .mat file
                    save([direc.main_directory direc.run_name direc.exp_nr ...                  % Save directory
                        direc.output_dir direc.dir_separator direc.exp_nr(1:end-1) 'crossyz_' cmb_newname{i} '.mat'], ...   % Save name
                        '-struct', 'saving_struct', '-v7.3'); 
                   % clear the data
                   clear saving_struct j cmb_matfile
                
                else
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    fprintf(['                    - !!! ERROR !!! \n'...
                        '                      The to-be-combined variables are not yet processed\n']) 
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    
                end         % <-- if sum( ismember(cmb_varsfiles, cmb_files) ) == length(cmb_varsfiles)
                
            end         % <-- for i = 1 : length(check.combine_newname)
            
        end         % <-- if strcmp(check.combine, 'yes')
        
    end         % <-- if isempty(yz_files) == 0
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
        fprintf('                  - YZ coss-section output is processed and saved!\n') 
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    
end     % <-- End of function
        
        

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%     READ XY FIELDDUMP output      %%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function f__NAMFIELDDUMPread_xy(direc, inp, check, II)
    
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
        fprintf(['\n                  - Start loading fielddump XY cross-section ' check.cross_lvl{II} '\n']) 
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    
    
    
    % Find output files with an XY cross-section
    fd_files = dir([direc.main_directory, direc.run_name, direc.exp_nr, ...
        'fielddump.*']);
    fd_files = {fd_files.name}';
    
    if isempty(fd_files) == 0
        
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %    Recenter cross-section data around scalar emission sources   %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        % If needed, prepare recentering the cross-section data around the
        % scalar emission sources% Define the order in which the y-direction data will be saved 
        INnew_y = [1 : str2num( inp.namoptions.value(ismember(inp.namoptions.option,{'jtot'})) )];

        % Check if the source needs to be centered if lhetero is turned on
        if strcmp(check.center_sources, 'yes') && ismember({'lhetero'},inp.namoptions.option)   
            % Check if lhetero is turned on
            if strcmp(inp.namoptions.value(ismember(inp.namoptions.option,{'lhetero'})),'.true.')

                % Define the namoptions indices which could contain the
                % land_use(XX,YY) option
                src_INtemp = find(ismember(inp.namoptions.option,{'lhetero'}))+1;
                src_INtemp = [src_INtemp; src_INtemp-1 + find(ismember(inp.namoptions.option(src_INtemp:end),'/'))];
                src_INtemp = [src_INtemp(1) :1: src_INtemp(2)];

                % Save the YY in land_use(XX,YY), i.e. the y position of the
                % source(s)
                src_ypos = [];
                for i = 1 : length(src_INtemp)
                    if length(inp.namoptions.option{src_INtemp(i)}) > 9
                        if strcmp(inp.namoptions.option{src_INtemp(i)}(1:9), 'land_use(')
                            if str2num(inp.namoptions.value{src_INtemp(i)}) > 0

                                src_IN_ypos = [find(ismember(inp.namoptions.option{src_INtemp(i)},','))+1 ...
                                           :1: ...
                                           find(ismember(inp.namoptions.option{src_INtemp(i)},')')) - 1];

                                src_ypos = [src_ypos ; str2num(inp.namoptions.option{src_INtemp(i)}(src_IN_ypos))];
                            end
                            
                        end
                    end
                end
                if isempty(src_ypos)
                    src_ypos = 1;
                end

                % Define the domain gridpoints in y-direction
                src_ygrid = [str2num( inp.namoptions.value(ismember(inp.namoptions.option,{'jtot'})) )];
                % Define the number of patches in y-direction
                src_ypatches = [str2num( inp.namoptions.value(ismember(inp.namoptions.option,{'ypatches'})) )];

                if src_ygrid ~= src_ypatches
                    sec_temp_ypos = src_ypos;
                    src_ypos = [];
                    for i = 1 : length(sec_temp_ypos)
                        src_ypos = [src_ypos, [(src_ygrid / src_ypatches * (sec_temp_ypos(i) -1) +1) ...
                                    :1: ...
                                    (src_ygrid / src_ypatches * sec_temp_ypos(i))] ];
                    end
                end
                src_ypos = floor(mean([min(src_ypos), max(src_ypos)]));
                src_change = floor(src_ygrid/2) - src_ypos;
                if src_change > 0
                    INnew_y = [INnew_y(end - src_change+1:end), INnew_y(1 : end-src_change)];
                elseif src_change < 0
                    INnew_y = [INnew_y(-1*src_change+1:end), INnew_y(1 : -1*src_change)];
                end
            end
        end
        
        clear src_*
        
        
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Setup the information on the fielddump  XY CROSSSECTION output  %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
       
        % Load the first file as a dummy
        fd_loc = [direc.main_directory, direc.run_name, direc.exp_nr, ...
            fd_files{1}];
        fd_info = ncinfo(fd_loc);
        % Define variables to be saved
        fd_var = {fd_info.Variables.Name}';
        % Preallocate time variable size
        fd_dimt = fd_info.Variables(ismember({fd_info.Variables.Name}, {'time'})).Size;
        % Preallocate x variable size
        fd_dimx = fd_info.Variables(ismember({fd_info.Variables.Name},{'xt'})).Size(1);
        fd_nprocx = str2num(inp.namoptions.value(ismember(inp.namoptions.option, 'nprocx')));
        % Preallocate y variable size
        fd_dimy = fd_info.Variables(ismember({fd_info.Variables.Name},{'yt'})).Size(1);
        fd_nprocy = str2num(inp.namoptions.value(ismember(inp.namoptions.option, 'nprocy')));
        
        
        
        % Change the default scalar names into the true scalar names
        fd_svnames_new = inp.scalar.Properties.VariableNames(2:end);
        fd_svnames_old = cell(length(fd_svnames_new),2);
        for i = 1 : length(fd_svnames_new)
            
            % Add the cross-section indicator (xz) to the scalar name
            fd_svnames_new{i} = [fd_svnames_new{i}];
            
            % Create the name that will be used to search in f_output.(cross_fields).INFO.LongNames
            fd_temp = '000';
            fd_temp(end-(length(num2str(i))-1) : end) = num2str(i);
            
            fd_svnames_old{i,1} = ['sv' fd_temp];
            fd_svnames_old{i,2} = ['scalar ' fd_temp]; 
        end
        
        % remove redundant variables (these are already saved in
        % cross.mat
        fd_var(ismember(fd_var, {'time','xt','xm','yt','ym','zt','zm'})) = [];
        tmp_processvar = check.process_var;
        for i = 1 : length(fd_svnames_new)
            if ismember(fd_svnames_new(i), tmp_processvar)
                tmp_processvar{ismember(tmp_processvar, fd_svnames_new(i))} = fd_svnames_old{i};
            end
        end
        fd_var(ismember(fd_var, tmp_processvar)==0) = [];
        
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %                       Save each variable                        %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        for i =  1 : length(fd_var)
            % Reset the "saving_struct"
            saving_struct = struct();

            % Define the name of the variable that will be saved to the
            % .mat file. Here the scalar names in scalar.inp are taken
            % into account
            if ismember(fd_var(i), fd_svnames_old(:,1))
                fd_savevar = fd_svnames_new{ismember(fd_svnames_old(:,1), fd_var(i))};
            else
                fd_savevar = fd_var{i};
            end

            saving_struct.(fd_savevar) = NaN(fd_dimx*fd_nprocx, fd_dimy*fd_nprocy, fd_dimt);
            for j = 1 : length(fd_files)
                % Determine the path of the to be loaded netCDF file 
                fd_loc = [direc.main_directory, direc.run_name, direc.exp_nr, fd_files{j}];
                % Define the core number that is read
                fd_locx = str2num(fd_files{j}(11:13));
                fd_locy = str2num(fd_files{j}(15:17));
                % Determine the indices that are filled by this netCDF
                % file
                fd_fillx = [1+fd_dimx*fd_locx : 1+fd_dimx*(fd_locx+1)-1];
                fd_filly = [1+fd_dimy*fd_locy : 1+fd_dimy*(fd_locy+1)-1];
                % Redefine the y-position of the data to match "INnew_y"
                fd_filly = find(ismember(INnew_y, fd_filly));
                
                % Load the data from the netCDF file
                saving_struct.(fd_savevar)(fd_fillx, fd_filly, :) =  ...
                    squeeze(ncread(fd_loc, fd_var{i}, [1 1 II-1 1], [Inf Inf 1 Inf]));
            end

            % Save the variable to the .mat file
            save([direc.main_directory direc.run_name direc.exp_nr ...                  % Save directory
                direc.output_dir direc.dir_separator direc.exp_nr(1:end-1) ...
                'crossxy' check.cross_lvl{II} '_' fd_savevar '.mat'], ...   % Save name
                '-struct', 'saving_struct', '-v7.3'); 
            % clear the data

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
            fprintf(['                    Variable "' fd_savevar '" is loaded and saved.\n']) 
            fprintf(['                    - XY fielddump cross-section ' check.cross_lvl{II} ' at ' num2str(i/length(fd_var)*100) '%%.\n']) 
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            clear saving_struct fd_savevar fd_loc fd_filly_new fd_filly j 
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %                  Combine variables (if needed)                  %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
         % If needed, combine the predefined variables
        if strcmp(check.combine, 'yes')
    
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
            fprintf('                  - Combine variables of the fielddump  XY coss-section!\n') 
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            cmb_newname = check.combine_newname;
            cmb_newname = cmb_newname(ismember(cmb_newname, check.process_var));
            cmb_vars = check.combine_vars(find( ismember(check.combine_newname, cmb_newname) ),:);
            cmb_files = dir([direc.main_directory direc.run_name direc.exp_nr ... 
                        direc.output_dir direc.dir_separator]);
            cmb_files = {cmb_files.name}';
            
            % Combine the variables
            for i = 1 : length(cmb_newname)
                
                % Define the processed file names of the to-be-combined
                % variables
                cmb_varsfiles = cmb_vars(i,:);
                for j = 1 : length(cmb_varsfiles)
                    cmb_varsfiles{j} = [direc.exp_nr(1:end-1) ...
                        'crossxy' check.cross_lvl{II} '_' cmb_varsfiles{j} '.mat'];
                end
                
                % Make sure the to-be-combined vairables are already
                % processed
                if sum( ismember(cmb_varsfiles, cmb_files) ) == length(cmb_varsfiles)
                    
                    % Reset saving_struct
                    saving_struct = struct();
    
                    % Add the first of the to be combined variables
                    cmb_matfile = matfile([direc.main_directory direc.run_name direc.exp_nr ... 
                        direc.output_dir direc.dir_separator direc.exp_nr(1:end-1) ...
                        'crossxy' check.cross_lvl{II} '_' cmb_vars{i,1} '.mat']);
                    saving_struct.(cmb_newname{i}) = cmb_matfile.(cmb_vars{i,1});
                    
                    % Add the remaining to be combined variables
                    for j = 2 : size(cmb_vars,2)
                        
                        cmb_matfile = matfile([direc.main_directory direc.run_name direc.exp_nr ... 
                            direc.output_dir direc.dir_separator ...
                            direc.exp_nr(1:end-1) 'crossxy' check.cross_lvl{II} '_' cmb_vars{i,j} '.mat']);
                        
                        saving_struct.(cmb_newname{i}) = saving_struct.(cmb_newname{i}) + ...
                             cmb_matfile.(cmb_vars{i,j});
                    end

                    % Save the combined variable to the .mat file
                    save([direc.main_directory direc.run_name direc.exp_nr ...                  % Save directory
                        direc.output_dir direc.dir_separator ...
                        direc.exp_nr(1:end-1) 'crossxy' check.cross_lvl{II} '_' cmb_newname{i} '.mat'], ...   % Save name
                        '-struct', 'saving_struct', '-v7.3'); 
                   % clear the data
                   clear saving_struct j cmb_matfile
                   
                else
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    fprintf(['                    - !!! ERROR !!! \n'...
                        '                      The to-be-combined variables are not yet processed\n']) 
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

                end         % <-- if sum( ismember(cmb_varsfiles, cmb_files) ) == length(cmb_varsfiles)
                
            end         % <-- for i = 1 : length(cmb_newname)
            
        end         % <-- if strcmp(check.combine, 'yes')
        
    end         % <-- if isempty(fd_files) == 0
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
        fprintf(['                  - The fielddump XY coss-section ' check.cross_lvl{II} ' output is processed and saved!\n']) 
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    
end     % <-- End of function
        
        

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%     READ XY FIELDDUMP output      %%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function f__CROSScoldstartTOwarmstart(direc, inp)
    	
    % Define the experiment number of the warmstart file
    cold_exp_nr = inp.namoptions.value(ismember(inp.namoptions.option, {'startfile'}));
    if strcmp(cold_exp_nr{1}(end), '''')
        cold_exp_nr = cold_exp_nr{1}(end-3:end-1);
    else
        cold_exp_nr = cold_exp_nr{1}(end-2:end);
    end
    % Define location of the .mat files of the warmstart experiment
    cold_matdir = dir([direc.main_directory direc.run_name cold_exp_nr direc.dir_separator 'output_v*']);
    cold_matdir = {cold_matdir.name}';
    cold_matdir = cold_matdir{end};
    % Define the names of the .mat files of the coldstart experiment
    cold_namelist = dir([direc.main_directory direc.run_name cold_exp_nr direc.dir_separator ...
        cold_matdir direc.dir_separator cold_exp_nr 'cross*']);
    cold_namelist = {cold_namelist.name}';
    % Remove the experiment ID from the namelist for comparison
    cold_namecompare = cell(length(cold_namelist),1);
    for i = 1 : length(cold_namelist)
        cold_namecompare{i} = cold_namelist{i}(4:end);
    end
    cold_namecompare(ismember(cold_namecompare, {'cross.mat'})) = [];
    % Define the names of the .mat files of the warmstart experiment
    warm_namelist = dir([direc.main_directory direc.run_name direc.exp_nr ... 
                    direc.output_dir direc.dir_separator ...
                    direc.exp_nr(1:end-1) 'cross*']);
    warm_namelist = {warm_namelist.name}';

    for i = 1 : length(warm_namelist)

        % Check if the NAMGENSTAT output exists
        if strcmp(warm_namelist{i}(4:end), 'cross.mat')

            % Load the input data of the coldstart run
            cold_matfile = matfile([direc.main_directory direc.run_name ...
                cold_exp_nr direc.dir_separator cold_matdir direc.dir_separator ...
                cold_exp_nr warm_namelist{i}(4:end)]);
            % Load the input data of the warmstart run
            warm_matfile = matfile([direc.main_directory direc.run_name ...
                direc.exp_nr direc.output_dir direc.dir_separator warm_namelist{i}]);
            % Allow changes to the warmstart output
            warm_matfile.Properties.Writable = true;

            % Add the coldstart output data to the warmstart output
            warm_matfile.time = cat(1, cold_matfile.time, warm_matfile.time);

        elseif ismember({warm_namelist{i}(4:end)}, cold_namecompare)

            % Load the input data of the coldstart run
            cold_matfile = matfile([direc.main_directory direc.run_name ...
                cold_exp_nr direc.dir_separator cold_matdir direc.dir_separator ...
                cold_exp_nr warm_namelist{i}(4:end)]);
            % Load the input data of the warmstart run
            warm_matfile = matfile([direc.main_directory direc.run_name ...
                direc.exp_nr direc.output_dir direc.dir_separator warm_namelist{i}]);
            % Allow changes to the warmstart output
            warm_matfile.Properties.Writable = true;
            % Define the variable name in the .mat file
            warm_var = whos(warm_matfile);
            warm_var = warm_var.name;

            % Add the coldstart output data to the warmstart output
            warm_matfile.(warm_var) = cat(3, cold_matfile.(warm_var), warm_matfile.(warm_var));

        end     % <-- if strcmp(warm_namelist{i}(4:end), 'cross.mat')

        clear warm_matfile warm_var cold_matfile

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
        fprintf(['                  - Adding run ' cold_exp_nr ' to run ' ...
            direc.exp_nr(1:end-1) ' is at ' num2str(i/length(warm_namelist)*100) '%%\n']) 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    end     % <-- for i = 1 : length(warm_namelist)
    
end     % End of function




