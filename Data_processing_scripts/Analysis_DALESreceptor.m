%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%         Ruben Schulte         %%%%%
%%%%%    DALES receptor analysis    %%%%%
%%%%%      started: 28-12-2020      %%%%%
%%%%%     restarted: 06-05-2021     %%%%%
%%%%%      changed: 22-07-2021      %%%%%
%%%%%       final: 23-05-2025       %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Note: This is a cleaned up version of Analysis_DALES_base2.m
% This script only covers the individual and bulk receptor analysis
%
% clear all

%% General settings
% These variables can be changed by the uses

% Define the relevant (data) directories 
direc.main_directory = 'D:\DALES_data\';
direc.run_name = 'nh3plume';
check.corr_loctime = 2;                     % Correct from UTC to local time
% Define indicators for background or plume concentration scalars
check.variable_background_indicator = 'b';  
check.variable_plume_indicator = 'p';


%% Analysis settings

% Receptors positions:
check.indiv_pos    = [ 750,  2390, 37.5 ...     % x-y-z positions [m] of the receptors {x1, y1, z1; x2, y2, z2; ....}
                     ];  
check.figure_width = 1000;
check.figure_font  = 18;    % Font size

%%%%%
% Receptor time series figure for the first individual receptor (check.indiv_pos)
%%%%%
check.switch_timeseries = 'no';        % Show result first INDIV receptor time series ('yes'/'no')
check.ts_var           = {...   % RunID, Scalar name, legend string, averaging time, color ID [1-4], LineStyle
                          '612', 'nh3b3' ,'10 s NH_{3,total}'  , 10/3600, 1, '-' ...
                        ; '612', 'nh3_b3b' ,'10 s NH_{3,bg}'     , 10/3600, 3, '-' ...
                        ; '612', 'nh3_b3p' ,'10 s NH_{3,plume}'  , 10/3600, 4, '-' ...
                        ; '612', 'nh3b3' ,'30 min NH_{3,total}'  , 1800/3600, 2, '-' ...
                        };
check.ts_txlim         = [12.5 17];         % Figure x-limits
check.ts_tanalysis     = [14.0, 17.0];      % start & end time analysis phase (color coded in figure)
check.ts_tbuffer       = [12.5, 14];        % start & end time buffer phase (color coded in figure)
check.ts_Clim          = [0 15];            % Figure C-limits (concentration limits, or y limits)
check.ts_ylabel        = 'NH_3 [ppb]';
check.ts_Cticks_manual = 'yes';        % Can only be used with linear axis
check.ts_Cticks_num    = unique([0:1:15]');     % y-ticks
check.ts_Cticks_str    = unique([0:5:15]');     % y-ticks labels

%%%%%
% Receptor statistics & time series figure for the first individual receptor (check.indiv_pos)
%%%%%
check.switch_recpstats = 'no';        % Show result first INDIV receptor time series ('yes'/'no')
check.rs_var           = {'612', 'nh3r1' ,'NH_{3,total}'  , 10/3600};
check.rs_txlim         = [14.0 17];
check.rs_tanalysis     = [14.0, 17.0];
check.rs_tbuffer       = [12.5, 14];
check.rs_Clim          = [0 16];
check.rs_ylabel        = 'NH_3 [ppb]';
check.rs_Cticks_num    = unique([0:1:15]');
check.rs_Cticks_str    = unique([0:5:15]');
check.rs_Ithreshhold   = 0.25;

% Receptor individual figures for the individual receptors
check.switch_INDIV     = 'no';        % Show statistics of individual receptors ('yes'/'no')
check.ind_runs         = {'012'}; 
check.ind_var          = 'nh3r1';
check.ind_avgt         = [10/3600];   
check.ind_twindow      = [14.0, 17.0];
check.ind_xylog        = 'no';
check.ind_Clim         = [0 15]; %[8 13]; %
check.ind_Ithreshhold  = 0.25;        % Threshold value [ppb] for Intermittency factor

% Receptor interactive figure for the first individual receptor
check.switch_interact   = 'no';        % Show interactive figure of first INDIV receptor ('yes'/'no')
check.int_var           = {'612', 'nh3r1'};
check.int_xy_xz         = 'xz';
check.int_avgt          = [10/3600, 0.5];
check.int_twindow       = [12.5 17];
check.int_logscale      = 'no';
% check.int_Clim          = [1e-5 1];
check.int_Clim          = [6 12];
% check.int_xlim          = [4 5250];
check.int_xlim          = [0 4500];
% check.int_ylim          = [0 4850];
check.int_ylim          = [0 2200];
check.int_Cticks_manual = 'yes';        % Can only be used with linear axis
% check.int_Cticks_num    = unique([6:0.5:13]');
% check.int_Cticks_str    = unique([6:1:13]');
% check.int_Cticks_num    = unique([1e-5 :1e-5: 1e-4, 1e-4 :1e-4: 0.001, 0.001 :0.001: 0.01, 0.01 : 0.01 : 0.1, 0.1 : 0.1 : 1]');
% check.int_Cticks_str    = unique([1e-5, 5e-5, 1e-4, 5e-4, 0.001, 0.005, 0.01, 0.05, 0.1, 0.5, 1]');
check.int_Cticks_num    = unique([6:0.5:12]');
check.int_Cticks_str    = unique([6:2:12]');


%%%%%%
% Receptor bulk analysis (distance Vs. statistics plots)
%%%%%%
check.switch_bulkstat  = 'no';      % Show receptor bulk analysis ('yes'/'no')
check.bs_xy_xz     = 'xy';          % type of cross-section ('xy'/'xz')
check.bs_xy_z      = 37.5;          % height [m] of the xy cross-seciton
check.bulk_pos_in  = { [600 : 20 : 9500] ...      % x-range of the in-plume receptor positions
                     , [2150 : 20 : 2640] ...     % y-range/z-range the in-plume of receptor positions  
                      };
check.bulk_pos_out = { [0 : 20 : 9500] ...     % x-range of out-of-plume receptor positions
                     , [0 : 20 : 200] ...     % y-range/z-range of out-of-plume receptor positions  
                      };
check.bs_var          = {'001', 'nh3r1' ...    % DALES run & variable
                        };
% check.bs_var          = {'801', 'nh3r1' ...    % DALES run & variable
%                        ; '801', 'nh3b1' ...
%                        ; '801', 'nh3b2' ...
%                        ; '801', 'nh3b3' ...
%                        ; '612', 'nh3r1' ...
%                        ; '804', 'nh3e1' ...
%                        ; '804', 'nh3e2' ...
%                        ; '804', 'nh3e3' ...
%                        ; '802', 'nh3c1' ...
%                        ; '612', 'nh3r1' ...
%                        ; '802', 'nh3c2' ...
%                        ; '802', 'nh3c3' ...
%                        ; '803', 'nh3d1' ...
%                        ; '803', 'nh3d2' ...
%                        ; '612', 'nh3r1' ...
%                        ; '803', 'nh3d3' ...
%                        ; '803', 'nh3d4' ...
%                        ; '311', 'nh3w1' ... 
%                        ; '312', 'nh3w2' ...
%                        ; '313', 'nh3w3' ...
%                        ; '612', 'nh3r1' ... 
%                        ; '314', 'nh3w4' ...
%                        ; '612', 'nh3r1' ...
%                         };
check.bs_stat         = {'fI',   [0 0.3] ...    % Statistic & fixed Ylim of figure...
                       ; 'F',   [0 0.3] ...    % Statistic & fixed Ylim of figure...
                       }; % Statistics, options are ['mean','std','fI','I','S','K', 'F']
check.bs_twindow      = [14.0 17.0];            % Analysis window
check.bs_avgt         = [10/3600];              % Averaging time [hours]
check.bs_Ithreshhold  = 0.25;                   % Threshold minimum value [ppb] for Intermittency factor
check.bs_dist_level   = [0.5 0.25, 0.10, 0.05]; % Threshold value [%] for blending-distance
check.bs_xlim         = [0 4500];               % x-limit of the figure
check.bs_show_xdist_fig = 'no';             % show figures for x-distance alongside the absolute distance figures


%%%%%%
% blending-distance comparison figures
%%%%%
check.switch_blenddist = 'no';          % Show blending-distance comparison ('yes'/'no')
% Statistics, options are ['mean','std','fI','I','S','K','F']
check.bd_stat         = 'F';           
% What to compare: 'variables', 'twindow', 'tavgt', 'zheight'. 
    % Multiple options in bd_compare create multiple figures
    % Can only compare 1 option for each figures, other options have first value as default
check.bd_compare      = { 'twindow' ...    % compare blending-distance for different twindows
                         };  
%                        ; 'zheight' ...    % compare blending-distance for different zheights
%                        ; 'variables' ...    % compare blending-distance for different variables
 
% Comparison options (first one is the default)
check.bd_zheight      = [[37.5, 7.5 : 5 :  117.5]'  ... 
                         ];
check.bd_var          = {'001', 'nh3r1' ...    % DALES run & variable
                        };
% check.bd_var          = {'801', 'nh3r1' ...    % DALES run & variable names
%                        ; '801', 'nh3b1' ...
%                        ; '801', 'nh3b2' ...
%                        ; '801', 'nh3b3' ...
%                        ; '612', 'nh3r1' ...
%                        ; '804', 'nh3e1' ...
%                        ; '804', 'nh3e2' ...
%                        ; '804', 'nh3e3' ...
%                        ; '802', 'nh3c1' ...
%                        ; '612', 'nh3r1' ...
%                        ; '802', 'nh3c2' ...
%                        ; '802', 'nh3c3' ...
%                        ; '803', 'nh3d1' ...
%                        ; '803', 'nh3d2' ...
%                        ; '612', 'nh3r1' ...
%                        ; '803', 'nh3d3' ...
%                        ; '803', 'nh3d4' ...
%                        ; '311', 'nh3w1' ... 
%                        ; '312', 'nh3w2' ...
%                        ; '313', 'nh3w3' ...
%                        ; '612', 'nh3r1' ... 
%                        ; '807', 'nh3w4' ...
%                        ; '806', 'nh3h1' ...
%                        ; '612', 'nh3r1' ...
%                        ; '111', 'nh3r1' ...
%                         };
check.bd_twindow      = [14.0, 17.0 ...     % Analysis phase start and end time [hours]
                       ; 13.0, 16.0 ...
                       ; 12.0, 15.0 ...
                       ; 11.0, 14.0 ...
                       ; 10.0, 13.0 ...
                       ; 09.0, 12.0 ...
                       ; 08.0, 11.0 ...
                        ];  
check.bd_avgt         = [10/3600 ...        % data averaging time [hours]
                       ; 60/3600 ...
                       ];  
check.bd_Ithreshhold  = 0.25;                   % Threshold value [ppb] for Intermittency factor
check.bd_dist_level   = [0.05 : 0.05: 0.5];     % Threshold value [ppb] for blending-distance
check.bd_peak_abs     = 'abs'; 
check.bd_max_centre   = 'max';      % Calculate bledning-distace using maximul ('max') or centreline ('centre') statistics
check.bd_savedata     = 'no';       % Save the resulting blending-distances ('yes'/'no')
% check.bd_savename     = 'BD_sense_F';    % Defaults: BD_sense_fI / BD_sense_F
check.bd_savename     = 'tempfilename';    % Save name for the blending-distance data




%% End of settings



% !!!!!!!!!!!!!!!!!!!! Do not change the code below !!!!!!!!!!!!!!!!!!!!



%% Make sure all runs are covered by "direc.exp_nr"

direc.exp_nr = {}; 
% check.ts_var: 
if strcmp(check.switch_timeseries, 'yes')
    for i = 1 : size(check.ts_var, 1)
        if ismember(check.ts_var(i,1), direc.exp_nr) == 0 
            direc.exp_nr(length(direc.exp_nr)+1) = check.ts_var(i,1);
        end
    end
end
% check.rs_var: 
if strcmp(check.switch_recpstats, 'yes')
    for i = 1 : size(check.rs_var, 1)
        if ismember(check.rs_var(i,1), direc.exp_nr) == 0 
            direc.exp_nr(length(direc.exp_nr)+1) = check.rs_var(i,1);
        end
    end
end
% check.ind_runs:
if strcmp(check.switch_INDIV, 'yes')
    for i = 1 : size(check.ind_var, 1)
        if ismember(check.ind_var(i,1), direc.exp_nr) == 0 
            direc.exp_nr(length(direc.exp_nr)+1) = check.ind_var(i,1);
        end
    end
end
% check.int_runs: 
if strcmp(check.switch_interact, 'yes')
    for i = 1 : size(check.int_var, 1)
        if ismember(check.int_var(i,1), direc.exp_nr) == 0 
            direc.exp_nr(length(direc.exp_nr)+1) = check.int_var(i,1);
        end
    end
end
% check.bs_var:
if strcmp(check.switch_bulkstat, 'yes')
    for i = 1 : size(check.bs_var, 1)
        if ismember(check.bs_var(i,1), direc.exp_nr) == 0 
            direc.exp_nr(length(direc.exp_nr)+1) = check.bs_var(i,1);
        end
    end
end
% check.bd_vars:
if strcmp(check.switch_blenddist, 'yes')
    for i = 1 : size(check.bd_var, 1)
        if ismember(check.bd_var(i,1), direc.exp_nr) == 0 
            direc.exp_nr(length(direc.exp_nr)+1) = check.bd_var(i,1);
        end
    end
end


if isempty(direc.exp_nr) == 1
    
    direc.exp_nr{1} = check.ts_var{1,1};
    
end

%% Make sure the directories end with an '\' character

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf('Prepare the DALES data\n')      
tload_start = datestr(now, 'dd-mm-yyyy HH:MM:SS');
fprintf([' --> Starting at: ' tload_start '\n']) 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Structure field names in direc structure
fields_direc = fields(direc);
% Loop through structure fields
for i = 1 : length(fields_direc) 
    fields_temp = direc.(char(fields_direc(i)));
    % If needed, add a '\' to the directory name string
    if strcmp(char(fields_temp(end)),'\') == 0
        direc.(char(fields_direc(i))) = strcat(fields_temp, '\');
    end
end

clear fields_* i

%% Load the DALES data

inp = struct();
for II = 1 : length(direc.exp_nr)
    
    runsavename = ['set' direc.exp_nr{II}(1:end-1)];
    if ismember(runsavename, fields(inp)) == 0
        
        dir_outp = dir([direc.main_directory direc.run_name direc.exp_nr{II} 'output_v*']);
        dir_outp = {dir_outp.name}';
        dir_outp = [dir_outp{end} '\'];
        
        % If it exists, load the cross.mat file as well
        check.(runsavename).loc = [direc.main_directory direc.run_name direc.exp_nr{II} dir_outp];
        check.(runsavename).namelist = dir(check.(runsavename).loc);    
        check.(runsavename).namelist = {check.(runsavename).namelist.name}';
        check.(runsavename).namelist(1:2) = [];
        
        if ismember({[direc.exp_nr{II}(1:end-1) 'cross.mat']}, ...
                check.(['set' direc.exp_nr{II}(1:end-1)]).namelist)
            
            inp.(runsavename)        = ...
                load([check.(runsavename).loc direc.exp_nr{II}(1:end-1) 'inp.mat']);
            outp_cross.(runsavename) = ...
                load([check.(runsavename).loc direc.exp_nr{II}(1:end-1) 'cross.mat']);
            
            outp_cross.(runsavename).loc = check.(runsavename).loc;
            outp_cross.(runsavename).namelist = ...
                dir([outp_cross.(runsavename).loc direc.exp_nr{II}(1:end-1) 'cross*']);    
            outp_cross.(runsavename).namelist = {outp_cross.(runsavename).namelist.name}';
        else
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            fprintf(['     !!!!! ERROR !!!!!\n'])
            fprintf(['     There is no cross-section output available for run ' direc.exp_nr{II}(1:end-1) '\n'])      
            fprintf(['     This matlab script is calceled. Try another DALES run.\n'])
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            return
        end

    else    
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        fprintf([' - ' runsavename ' is already loaded\n'])      
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    end
    clear runsavename
end
clear II


%% Define receptor positions in the grid and save the time-series

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    fprintf(' - Define the position of the receptors on the cross-section gird\n')      
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Translate the x-y-z ranges of the bulk-stat in-plume receptors to coordinates
if strcmp(check.bs_xy_xz, 'xz')
    check.bulk_pos_in{3} = check.bulk_pos_in{2};
    check.bulk_pos_in{2} = outp_cross.(['set' check.bs_var{1,1}]).ypos;
else
    check.bulk_pos_in{3} = check.bs_xy_z;
end
recp_pos = NaN(length(check.bulk_pos_in{1}) * ...
             length(check.bulk_pos_in{2}) * ...
             length(check.bulk_pos_in{3}), 3);
counter  = 1;
for i = 1 : length(check.bulk_pos_in{1})

    for j = 1 : length(check.bulk_pos_in{2})

        for k = 1 : length(check.bulk_pos_in{3})
            recp_pos(counter,:) = [check.bulk_pos_in{1}(i), ...
                                 check.bulk_pos_in{2}(j), ...
                                 check.bulk_pos_in{3}(k) ];

            counter = counter + 1;
        end
    end
end
check.bulk_pos_in = recp_pos;
clear counter recp_pos i j k

% Translate the x-y-z ranges of the bulk-stat out-of-plume receptors to coordinates
if strcmp(check.bs_xy_xz, 'xz')
    check.bulk_pos_out{3} = check.bulk_pos_out{2};
    check.bulk_pos_out{2} = outp_cross.(['set' check.bs_var{1,1}]).ypos;
else
    check.bulk_pos_out{3} = check.bs_xy_z;
end
recp_pos = NaN(length(check.bulk_pos_out{1}) * ...
             length(check.bulk_pos_out{2}) * ...
             length(check.bulk_pos_out{3}), 3);
counter  = 1;
for i = 1 : length(check.bulk_pos_out{1})

    for j = 1 : length(check.bulk_pos_out{2})

        for k = 1 : length(check.bulk_pos_out{3})
            recp_pos(counter,:) = [check.bulk_pos_out{1}(i), ...
                                 check.bulk_pos_out{2}(j), ...
                                 check.bulk_pos_out{3}(k) ];

            counter = counter + 1;
        end
    end
end
check.bulk_pos_out = recp_pos;
clear counter recp_pos i j k

recp_name = 'bulk_in';
outp_cross = f__receptor_position(inp, outp_cross, check.bulk_pos_in, recp_name);

recp_name = 'bulk_out';
outp_cross = f__receptor_position(inp, outp_cross, check.bulk_pos_out, recp_name);

recp_name = 'indiv';
outp_cross = f__receptor_position(inp, outp_cross, check.indiv_pos, recp_name);

clear recp_*

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
tload_spent = [tload_hour*24, tload_min*24*60, tload_sec*24*60*60];
% Message that loading the data is finished
fprintf([' --> DONE!\n' ...
         '       Time spent = '   num2str(tload_spent(1)),' hour(s)' ...
         '\n                    ' num2str(tload_spent(2)),' minutes(s) ' ...
         '\n                    ' num2str(tload_spent(3)) ' second(s)\n'])
 clear tload_start tload_end tload_spent tload_sec tload_min tload_hour tload_day tload_temp
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Receptor time series observation figures

if strcmp(check.switch_timeseries, 'yes')
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    fprintf('\nPlot receptor observation time series for the first individual receptor \n')  
    tload_start = datestr(now, 'dd-mm-yyyy HH:MM:SS');
    fprintf([' --> Starting at: ' tload_start '\n']) 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
    % Make a time series figure where multiple runs can be compared
    rc_position   = [10 230 check.figure_width 0.35*check.figure_width];  % Position of the figure
    rc_plotdx_pos = 50;                 % Movement towards the right of each subsequent figure [pixels]


    % Make the figure
    f__FIGtimeseries(outp_cross, check, rc_position)
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
    tload_spent = [tload_hour*24, tload_min*24*60, tload_sec*24*60*60];
    % Message that loading the data is finished
    fprintf([' --> DONE!\n' ...
             '       Time spent = '   num2str(tload_spent(1)),' hour(s)' ...
             '\n                    ' num2str(tload_spent(2)),' minutes(s) ' ...
             '\n                    ' num2str(tload_spent(3)) ' second(s)\n'])
 clear tload_start tload_end tload_spent tload_sec tload_min tload_hour tload_day tload_temp
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end

%% Receptor time series statistics figure

if strcmp(check.switch_recpstats, 'yes')
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    fprintf('\nPlot single receptor  time series with statistics\n')  
    tload_start = datestr(now, 'dd-mm-yyyy HH:MM:SS');
    fprintf([' --> Starting at: ' tload_start '\n']) 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
    % Make a time series figure where multiple runs can be compared
    rs_position   = [10 230 check.figure_width 0.55*check.figure_width];  % Position of the figure
    rs_plotdx_pos = 50;                 % Movement towards the right of each subsequent figure [pixels]

    % Make the figure
    f__FIGtimeseries_stats(outp_cross, check, rs_position)
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
    tload_spent = [tload_hour*24, tload_min*24*60, tload_sec*24*60*60];
    % Message that loading the data is finished
    fprintf([' --> DONE!\n' ...
             '       Time spent = '   num2str(tload_spent(1)),' hour(s)' ...
             '\n                    ' num2str(tload_spent(2)),' minutes(s) ' ...
             '\n                    ' num2str(tload_spent(3)) ' second(s)\n'])
 clear tload_start tload_end tload_spent tload_sec tload_min tload_hour tload_day tload_temp
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end


%% Receptor interactive figure

if strcmp(check.switch_interact, 'yes')
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    fprintf('- Plot receptor interactive figures\n')  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % For the first variable and run, make an interactive figure with the
    % cross section and the time series
    int_position   = [10 50 1500 550];  % Position of the figure
    f__FigReceptor_crosstime(outp_cross.(['set' check.int_var{1}]), check, check.int_var{2}, int_position)

end


%% Single receptor statistics


if strcmp(check.switch_INDIV,'yes') 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    fprintf('\nShowing individual receptor statistics\n')      
    tload_start = datestr(now, 'dd-mm-yyyy HH:MM:SS');
    fprintf([' --> Starting at: ' tload_start '\n']) 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if strcmp(check.ind_var(1:3), 'nh3')
        
        tic
        
        si_runs = cell(length(check.ind_runs),1);
        for i = 1 : length(check.ind_runs)
            si_runs{i} = ['set' check.ind_runs{i}];
        end

        si_figpos = [30 130 1500 670];
        for i = 1 : length(si_runs)

            % Show the position of the receptors
            f__FIGreceptor_xypos(outp_cross.(si_runs{i}), check, si_runs{i})


            %%%%%
            % Plot the individual receptor statistics
            %%%%%
            for j = size(outp_cross.(si_runs{i}).recp_indiv,1) :-1: 1
                si_figpos_new = si_figpos + [90*(i-1) -40*(size(outp_cross.(si_runs{i}).recp_indiv,1)-j) 0 0];


                f__FIGindividual_stats( outp_cross.(si_runs{i}), check, si_figpos_new, j )




            end     % <-- j = size(outp_cross.(si_runs{i}).recp_indiv,1) :-1: 1


        end     % <-- for i = 1 : length(si_runs)
        toc
        
    else
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        fprintf('  Individual receptor statistics are not shown.\n  Consider picking an NH3 variable for "indiv_var"\n')      
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
    end     % <-- if strcmp(check.ind_var(1:3), 'nh3')
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
    tload_spent = [tload_hour*24, tload_min*24*60, tload_sec*24*60*60];
    % Message that loading the data is finished
    fprintf([' --> DONE!\n' ...
             '       Time spent = '   num2str(tload_spent(1)),' hour(s)' ...
             '\n                    ' num2str(tload_spent(2)),' minutes(s) ' ...
             '\n                    ' num2str(tload_spent(3)) ' second(s)\n'])
 clear tload_start tload_end tload_spent tload_sec tload_min tload_hour tload_day tload_temp
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
end     % <-- if strcmp(check.switch_INDIV,'yes')



%% Time bulk statistics analysis

if strcmp(check.switch_bulkstat, 'yes')
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    fprintf('\nShowing bulk statistics with distance from the source\n')  
    tload_start = datestr(now, 'dd-mm-yyyy HH:MM:SS');
    fprintf([' --> Starting at: ' tload_start '\n']) 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
    
   
    nh3_blendthreshold = NaN(size(check.bs_var ,1), size(check.bs_stat,1),...
        length(check.bs_dist_level));
    for i = 1 : size(check.bs_var ,1)
        
        for j = 1 : size(check.bs_stat,1)
            
            
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                fprintf([' - Plotting ' check.bs_stat{j,1} ' Vs. distance for ' ...
                    check.bs_var{i,2} ' set ' check.bs_var{i,1} ' \n'])      
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            nh3_blendthreshold(i,j,:) = f__FIGbulk_stats(outp_cross, check, i, j);
            
        end
    end
    
    
    
    
    %%%%%
    % Prepare for the legend figure
    %%%%%
    nh3_colormap_Delta = length(check.bs_dist_level);
    if mod(nh3_colormap_Delta/2, 1) == 0
        nh3_colormap = [[ [094 : (255 - 094)/(nh3_colormap_Delta ) : 255 ]', ...
                         [060 : (255 - 060)/(nh3_colormap_Delta ) : 255 ]', ...
                         [153 : (255 - 153)/(nh3_colormap_Delta ) : 255 ]' ] ./255];
    else
        nh3_colormap = [[ [094 : (255 - 094)/ceil(nh3_colormap_Delta ) : 255 ]', ...
                         [060 : (255 - 060)/ceil(nh3_colormap_Delta ) : 255 ]', ...
                         [153 : (255 - 153)/ceil(nh3_colormap_Delta ) : 255 ]' ] ./255];
    end
    nh3_colorset = {[230,097,001]./255, [094,060,153]./255,...
                   [253,184,099]./255, [178,171,210]./255};
    
    nh3_thresholds = squeeze(nh3_blendthreshold(1,1,:));
    nh3_lgndstr = cell(1,length(nh3_thresholds));
    for i = 1 : length(nh3_thresholds)
        nh3_lgndstr{i} = [num2str(nh3_thresholds(i)*100) '% threshold'];
    end
    nh3_lgndstr = {'Individual receptors   _{ }','Maximum   _{ }', ...
        'Centreline _{ }','Blending-distance    _{ }', nh3_lgndstr{:}};
    if length(nh3_lgndstr) > 4
        nh3_lgnd_columns = 4;
    else
        nh3_lgnd_columns = length(nh3_lgndstr);
    end
    nh3_lgnd_rows = ceil(length(nh3_lgndstr) / nh3_lgnd_columns);
   
    
    %%%%%
    % Make the legend figure
    %%%%%
    figure('units','pixels', 'Color','w',...
        'innerposition', [10 600 check.figure_width 30*nh3_lgnd_rows], ...
        'Name', 'Receptor stats');
    subplot('Position',[-1 -1 0.8 0.5])
    hold on
    nh3_obj = gobjects(4 + length(check.bs_dist_level),1);
    nh3_obj(1) = scatter([0 1],[1 1], ...
        'Marker','o', 'LineWidth',3, 'SizeData',20, ...
        'MarkerFaceColor',nh3_colorset{3}, 'MarkerFaceAlpha',0.99, ...
        'MarkerEdgeColor','none', 'MarkerEdgeAlpha',0.99 );
    nh3_obj(2) = plot([0 1],[2 2],'LineStyle','-', 'LineWidth',3.5, 'Color',nh3_colorset{1});
    nh3_obj(3) = plot([0 1], [3 3], 'LineStyle',':', 'LineWidth',3.5, 'Color',nh3_colorset{1});
    nh3_obj(4) = plot([0 1],[4 4], '--k', 'LineWidth',3.5);
    for i = 1 : length(check.bs_dist_level)
        nh3_obj(4+i) = plot([0 1], [4+i 4+i], ...
            'LineStyle','-', 'LineWidth',3.0, 'Color',nh3_colormap(i,:));
    end
    set(gca, 'FontSize', check.figure_font)
    nh3_lgnd = legend(nh3_obj,nh3_lgndstr,'box','off', ...
        'NumColumns',nh3_lgnd_columns, 'Orientation','Horizontal');
    nh3_lgnd.Position = [0.01 0.04 0.25*nh3_lgnd_columns 0.96];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
    tload_spent = [tload_hour*24, tload_min*24*60, tload_sec*24*60*60];
    % Message that loading the data is finished
    fprintf([' --> DONE!\n' ...
             '       Time spent = '   num2str(tload_spent(1)),' hour(s)' ...
             '\n                    ' num2str(tload_spent(2)),' minutes(s) ' ...
             '\n                    ' num2str(tload_spent(3)) ' second(s)\n'])
 clear tload_start tload_end tload_spent tload_sec tload_min tload_hour tload_day tload_temp
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
end

%% blending-distance analysis

if strcmp(check.switch_blenddist, 'yes')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    fprintf('\nShowing the blending-distance at the plume centerline\n')  
    tload_start = datestr(now, 'dd-mm-yyyy HH:MM:SS');
    fprintf([' --> Starting at: ' tload_start '\n']) 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
    tic
    
    for i = 1 : length(check.bd_compare)
    
        bd_safedistance{i} = f__FIGsafedistance(outp_cross, check, check.bd_compare{i});
    end
    
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
    tload_spent = [tload_hour*24, tload_min*24*60, tload_sec*24*60*60];
    % Message that loading the data is finished
    fprintf([' --> DONE!\n' ...
             '       Time spent = '   num2str(tload_spent(1)),' hour(s)' ...
             '\n                    ' num2str(tload_spent(2)),' minutes(s) ' ...
             '\n                    ' num2str(tload_spent(3)) ' second(s)\n'])
	clear tload_start tload_end tload_spent tload_sec tload_min tload_hour tload_day tload_temp
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %% Prepare the data for the official blending distance figure

    if strcmp(check.bd_savedata,'yes')
        bd_Bd = bd_safedistance{1};
        save([direc.main_directory direc.run_name check.bd_savename '.mat'], 'bd_Bd', '-v7.3')
    end
end
    
%%
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf('\nThe script is FINISHED!\n')      
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%




'THE END'



















































function outp_cross = f__receptor_position(inp, outp_cross, recp_pos, recp_name)
    
    
    rc_runs = fields(inp);
    for i = 1 : length(rc_runs)

        % Position receptors on x-y-z- grid
        for j = 1 : size(recp_pos,1)
            % x
            [~, rc_sortIN] = sort(abs(outp_cross.(rc_runs{i}).xt - recp_pos(j,1)));
            recp_pos(j,1) = outp_cross.(rc_runs{i}).xt(rc_sortIN(1));
            % y
            [~, rc_sortIN] = sort(abs(outp_cross.(rc_runs{i}).yt - recp_pos(j,2)));
            recp_pos(j,2) = outp_cross.(rc_runs{i}).yt(rc_sortIN(1));
            % z
            [~, rc_sortIN] = sort(abs(outp_cross.(rc_runs{i}).zt - recp_pos(j,3)));
            recp_pos(j,3) = outp_cross.(rc_runs{i}).zt(rc_sortIN(1));
        end
        clear rc_sortIN
        
        for j = size(recp_pos,1) :-1: 1
            rc_IN = find(recp_pos(:,1) == recp_pos(j,1) & ...
                recp_pos(:,2) == recp_pos(j,2) & ...
                recp_pos(:,3) == recp_pos(j,3));
            if length(rc_IN) > 1
                recp_pos(j,:) = [];
            end
        end
        clear rc_IN j
        
        outp_cross.(rc_runs{i}).(['recp_' recp_name]) = NaN(size(recp_pos,1),3);
        outp_cross.(rc_runs{i}).(['recp_' recp_name '_file']) = cell(size(recp_pos,1),1);

        for j = 1 : size(recp_pos,1)
            
            if ismember(round(recp_pos(j,3),6), round(outp_cross.(rc_runs{i}).zpos,6))
                % Define the cross-section file in which the receptor can
                % be found
                INz = find( ismember(outp_cross.(rc_runs{i}).zpos, recp_pos(j,3)) );
                outp_cross.(rc_runs{i}).(['recp_' recp_name '_file']){j} = ...
                    [rc_runs{i}(4:end) 'crossxy' outp_cross.(rc_runs{i}).cross_lvl{INz} '_'];
                % Define the receptor position
                outp_cross.(rc_runs{i}).(['recp_' recp_name])(j,:) = recp_pos(j,:);
                
            elseif round(recp_pos(j,1),6) == round(outp_cross.(rc_runs{i}).xpos,6)
                % Define the cross-section file in which the receptor can
                % be found
                outp_cross.(rc_runs{i}).(['recp_' recp_name '_file']){j} = ...
                    [rc_runs{i}(4:end) 'crossyz_'];
                % Define the receptor position
                outp_cross.(rc_runs{i}).(['recp_' recp_name])(j,:) = recp_pos(j,:);
                
            elseif round(recp_pos(j,2),6) == round(outp_cross.(rc_runs{i}).ypos,6)
                % Define the cross-section file in which the receptor can
                % be found
                outp_cross.(rc_runs{i}).(['recp_' recp_name '_file']){j} = ...
                    [rc_runs{i}(4:end) 'crossxz_'];
                % Define the receptor position
                outp_cross.(rc_runs{i}).(['recp_' recp_name])(j,:) = recp_pos(j,:);
                
            else
                fprintf(['     !   WARNING   !\n'...
                    '     "' recp_name '" receptor nr ' num2str(j) ', at x-y-z = [' ...
                    num2str(recp_pos(j,1)) '-' num2str(recp_pos(j,2)) '-' num2str(recp_pos(j,3)) ...
                    '] is not located on one of the cross-sections.\n' ...
                    '     This receptor will be ignored\n'])
            end
        end

    end     % <-- for i = 1 : length(rc_runs)
    clear rc_* i j
    
end




%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%                                                       %%%%%%%%%%
%%%%%%%%%%        Functions for plotting profile figures         %%%%%%%%%%
%%%%%%%%%%                                                       %%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%










%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%  Plot first receptor time series  %%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function f__FIGtimeseries(outp, check, rc_position)
    
    
    
    rc_colorset = {[230,097,001]./255, [253,184,099]./255, ...
                   [094,060,153]./255, [178,171,210]./255};
    rc_linestyle = {'-',':','--','-.'};
    rc_fields = fields(outp);
    rc_xlabel = 'CEST';
    rc_ylabel = check.ts_ylabel;
    rc_xlim = check.ts_txlim;
    rc_title = ['Example of simulated NH_3 point measurement'];
    
    % Load the data
    rc_minmax = NaN(size(check.ts_var,1), 2);
    rc_data = cell(size(check.ts_var,1),1);
    for i = 1 : size(check.ts_var,1)
        
        rc_section = outp.(['set' check.ts_var{i,1}]).recp_indiv_file{1}(9:10);
        if strcmp(rc_section, 'xy')
            rc_INpos = [find( outp.(['set' check.ts_var{i,1}]).xt == outp.(['set' check.ts_var{i,1}]).recp_indiv(1,1) ), ...
                find( outp.(['set' check.ts_var{i,1}]).yt == outp.(['set' check.ts_var{i,1}]).recp_indiv(1,2) )];
        elseif strcmp(rc_section, 'xz')
            rc_INpos = [find( outp.(['set' check.ts_var{i,1}]).xt == outp.(['set' check.ts_var{i,1}]).recp_indiv(1,1) ), ...
                find( outp.(['set' check.ts_var{i,1}]).zt == outp.(['set' check.ts_var{i,1}]).recp_indiv(1,3) )];
        elseif strcmp(rc_section, 'yz')
            rc_INpos = [find( outp.(['set' check.ts_var{i,1}]).yt == outp.(['set' check.ts_var{i,1}]).recp_indiv(1,2) ), ...
                find( outp.(['set' check.ts_var{i,1}]).zt == outp.(['set' check.ts_var{i,1}]).recp_indiv(1,3) )];
        end
%         rc_INt = find( outp.(['set' check.ts_var{i,1}]).time >= check.ts_txlim(1) & ...
%             outp.(['set' check.ts_var{i,1}]).time <= check.ts_txlim(2) );
        rc_INt = [1: 1 : length(outp.(['set' check.ts_var{i,1}]).time )]';
        
        rc_time_temp = outp.(['set' check.ts_var{i,1}]).time(rc_INt);
        rc_matfile = matfile([check.(['set' check.ts_var{i,1}]).loc ...
             outp.(['set' check.ts_var{i,1}]).recp_indiv_file{1} check.ts_var{i,2}]);
        
        rc_time = rc_time_temp;
        if round(check.ts_var{i,4}, 6) > round(median( rc_time_temp(2:end) - rc_time_temp(1:end-1) ), 6)

            rc_time = [rc_time_temp(1) - median( rc_time_temp(2:end) - rc_time_temp(1:end-1) ) +  check.ts_var{i,4} ...
                : check.ts_var{i,4} : rc_time_temp(end)]';
            rc_C = NaN(length(rc_time),1);

            if length(rc_time) > 1
                for k = 1 : length(rc_time)
                    if k == 1 
                        rc_C(k) = mean( rc_matfile.(check.ts_var{i,2})(rc_INpos(1), rc_INpos(2), ...
                            rc_INt(rc_time_temp <= rc_time(k)))  );
                    else
                        rc_C(k) = mean( rc_matfile.(check.ts_var{i,2})(rc_INpos(1), rc_INpos(2), ...
                            rc_INt(rc_time_temp > rc_time(k-1) & rc_time_temp <= rc_time(k)))  );
                    end
                end
            end
            rc_data{i} = [rc_time, rc_C];
            
        elseif round(check.ts_var{i,4}, 6) == round(median( rc_time_temp(2:end) - rc_time_temp(1:end-1) ), 6)

            rc_data{i} = [rc_time, ...
                squeeze(rc_matfile.(check.ts_var{i,2})(rc_INpos(1), rc_INpos(2), rc_INt))];
        else
             rc_data{i} = [NaN, NaN; NaN, NaN];
        end
        
        rc_minmax(i,:) = [min(rc_data{i}(:,2)), max(rc_data{i}(:,2))];
        
    end
    
    rc_xtick = [0:0.25:24];
    rc_xticklabel = cell(size(rc_xtick));
    rc_xtick_dt = 0.5;
    for i = 1 : length(rc_xticklabel)
        if ismember(rc_xtick(i), [0 : rc_xtick_dt : 24])
            rc_hour = floor(rc_xtick(i));
            rc_min = round((rc_xtick(i) - rc_hour) * 60);
            
            rc_hour = num2str(rc_hour);
            if length(rc_hour) == 1
                rc_hour = ['0' rc_hour];
            end
            rc_min = num2str(rc_min);
            if length(rc_min) == 1
                rc_min = ['0' rc_min];
            end
            rc_xticklabel{i} = [rc_hour ':' rc_min];
        end
    end
    clear rc_hour rc_min
    
    %%
    
    % Build the figure
    figure('units','pixels', 'Color','w',...
        'innerposition', rc_position, ...
        'Name', ['t_window = ' num2str(check.ts_txlim(1)) '-' num2str(check.ts_txlim(2))]);
    rc_sub = subplot('Position',[0.085 0.215 0.7 0.76]);
    
    hold on         % show multiple lines
    grid on
    rc_ax1 = gca;   % axis variable
    
    for i = 1 : length(rc_xticklabel)
        if ismember(rc_xtick(i), [0 : rc_xtick_dt : 24])
            plot([rc_xtick(i) rc_xtick(i)], [-9999 9999], 'Color',[0 0 0 0.25], 'LineWidth', 1.2)
        end
    end
    for i = 1 : length(check.ts_Cticks_str)
        plot([-9999 9999], [check.ts_Cticks_str(i) check.ts_Cticks_str(i)], ...
            'Color',[0 0 0 0.20], 'LineWidth', 1)
    end
    % Show buffer phase
    patch([check.ts_tbuffer(1) check.ts_tbuffer(2) check.ts_tbuffer(2) check.ts_tbuffer(1)], ...
        [-9999 -9999 9999 9999], [1 1 1 1], 'FaceColor', 'k', 'FaceAlpha',0.1, 'EdgeColor','none')
    
    for i = 1 : size(check.ts_var,1)
            plot(rc_data{i}(:,1), rc_data{i}(:,2), 'Color',rc_colorset{check.ts_var{i,5}}, ...
                'LineStyle', check.ts_var{i,6}, 'LineWidth',4)
    end
        
    
    if sum(isnan(check.ts_Clim)) > 0
        rc_ylim = [min(rc_minmax(:,1)), max(rc_minmax(:,2))];
    else
        rc_ylim = check.ts_Clim;
    end
    % For 1 line, plot the average and the standard deviation range
    if size(check.ts_var,1) == 1
        plot([rc_xlim(1), rc_xlim(2)], [mean(rc_data{1}(2:end,2)), mean(rc_data{1}(2:end,2))], ...
            '-k', 'LineWidth',4)
        plot([rc_xlim(1), rc_xlim(2)], [mean(rc_data{1}(2:end,2)) - std(rc_data{1}(2:end,2)), mean(rc_data{1}(2:end,2)) - std(rc_data{1}(2:end,2))], ...
            ':k', 'LineWidth',4)
        plot([rc_xlim(1), rc_xlim(2)], [mean(rc_data{1}(2:end,2)) + std(rc_data{1}(2:end,2)), mean(rc_data{1}(2:end,2)) + std(rc_data{1}(2:end,2))], ...
            ':k', 'LineWidth',4)
        if strcmp(check.ts_var{1,2}(end),'p')
            plot([rc_xlim(1), rc_xlim(2)], [check.rs_Ithreshhold, check.rs_Ithreshhold], ...
                '--k', 'LineWidth',2)
        end
    end
    
    if strcmp(check.ts_Cticks_manual, 'yes')
        rc_yticks = cell(length(check.ts_Cticks_num),1);
        for i = 1 : length(check.ts_Cticks_num)
            if ismember( check.ts_Cticks_num(i), check.ts_Cticks_str )
                rc_yticks{i} = num2str(check.ts_Cticks_num(i));
            end
        end
            
        set(rc_ax1,'FontSize',check.figure_font, 'XLim',rc_xlim, 'YLim',rc_ylim, ...
            'YTick',check.ts_Cticks_num, 'YTickLabel',rc_yticks, ...
            'XTick',rc_xtick, 'XTickLabel',rc_xticklabel)
    else
        set(rc_ax1,'FontSize',check.figure_font, 'XLim',rc_xlim, 'YLim',rc_ylim, ...
            'XTick',rc_xtick, 'XTickLabel',rc_xticklabel)
    end
    xlabel(rc_xlabel);
    ylabel(rc_ylabel);
%     title(rc_title)    

%     plot([rc_data{1}(1,1), rc_data{1}(end,1)] , [mean(rc_data{1}(:,2)) mean(rc_data{1}(:,2))], ...
%         '-k','LineWidth', 3)
%     plot([rc_data{1}(1,1), rc_data{1}(end,1)] , ...
%         [mean(rc_data{1}(:,2)) + std(rc_data{1}(:,2)), mean(rc_data{1}(:,2)) + std(rc_data{1}(:,2))], ...
%         ':k','LineWidth', 3)
%     plot([rc_data{1}(1,1), rc_data{1}(end,1)] , [0.25 0.25], '--k','LineWidth', 3)
    
    % Make the figure legends
    rc_lgnd_pos  = [0.74     0.9   0.17    0.055];
    rc_lgnd_dpos = [0.00    -0.13  0.00    0.005];
    rc_lgndtavg = struct();
    rc_lgndrun = struct();    
    
    for j = 1 : size(check.ts_var,1)
        
        subplot('Position',[1+j 1+j 0.50 0.78])
        set(gca, 'FontSize',18);
        rc_lgndrunline = plot([0 1], [j-1 j], 'Color', rc_colorset{check.ts_var{j,5}}, ...
            'LineStyle',check.ts_var{j,6}, 'LineWidth', 4);
        hold on

        % Variable legend
        rc_lgndrun.(['line_' num2str(j)]) = legend(rc_lgndrunline, {''});
        % Set the legend position
        rc_lgndrun.(['line_' num2str(j)]).Position = ...
            rc_lgnd_pos + [0 0 0 0] + (j-1) .* rc_lgnd_dpos;
        % Turn the legend box off
        rc_lgndrun.(['line_' num2str(j)]).Box = 'off';
        % Add the run name to the legend
%         annotation('textbox', rc_lgnd_pos + [0.1 -0.4 0 0] + (j-1) .* rc_lgnd_dpos, ...
%                 'string', [check.ts_var{j,1} '-' check.ts_var{j,3}], ...
%                 'EdgeColor','none', 'FontSize',16, 'FontWeight','normal')
        annotation('textbox', rc_lgnd_pos + [0.1 0.0 0 0] + (j-1) .* rc_lgnd_dpos, ...
                'string', check.ts_var{j,3}, ...
                'EdgeColor','none', 'FontSize',16, 'FontWeight','normal')
    end   
    if size(check.ts_var,1) == 1
        subplot('Position',[1+j+2 1+j+2 0.50 0.78])
        set(gca, 'FontSize',18);
        rc_lgndrunline = plot([0 1], [j+2-1 j+2], '-k', 'LineWidth', 4);
        hold on
        % Variable legend
        rc_lgndrun.(['line_' num2str(j+2)]) = legend(rc_lgndrunline, {''});
        % Set the legend position
        rc_lgndrun.(['line_' num2str(j+2)]).Position = ...
            rc_lgnd_pos + [0 0 0 0] + (j+2-1) .* rc_lgnd_dpos;
        % Turn the legend box off
        rc_lgndrun.(['line_' num2str(j+2)]).Box = 'off';
        % Add the run name to the legend
    %         annotation('textbox', rc_lgnd_pos + [0.1 -0.4 0 0] + (j-1) .* rc_lgnd_dpos, ...
    %                 'string', [check.ts_var{j,1} '-' check.ts_var{j,3}], ...
    %                 'EdgeColor','none', 'FontSize',16, 'FontWeight','normal')
        annotation('textbox', rc_lgnd_pos + [0.1 0.0 0 0] + (j+2-1) .* rc_lgnd_dpos, ...
                'string', 'mean', ...
                'EdgeColor','none', 'FontSize',16, 'FontWeight','normal')
            
        subplot('Position',[1+j+3 1+j+3 0.50 0.78])
        set(gca, 'FontSize',18);
        rc_lgndrunline = plot([0 1], [j+3-1 j+3], ':k', 'LineWidth', 4);
        hold on
        % Variable legend
        rc_lgndrun.(['line_' num2str(j+3)]) = legend(rc_lgndrunline, {''});
        % Set the legend position
        rc_lgndrun.(['line_' num2str(j+3)]).Position = ...
            rc_lgnd_pos + [0 0 0 0] + (j+3-1) .* rc_lgnd_dpos;
        % Turn the legend box off
        rc_lgndrun.(['line_' num2str(j+3)]).Box = 'off';
        % Add the run name to the legend
    %         annotation('textbox', rc_lgnd_pos + [0.1 -0.4 0 0] + (j-1) .* rc_lgnd_dpos, ...
    %                 'string', [check.ts_var{j,1} '-' check.ts_var{j,3}], ...
    %                 'EdgeColor','none', 'FontSize',16, 'FontWeight','normal')
        annotation('textbox', rc_lgnd_pos + [0.1 0.0 0 0] + (j+3-1) .* rc_lgnd_dpos, ...
                'string', 'mean + \sigma', ...
                'EdgeColor','none', 'FontSize',16, 'FontWeight','normal')
            
        if strcmp(check.ts_var{1,2}(end),'p')
            subplot('Position',[1+j+4 1+j+4 0.50 0.78])
            set(gca, 'FontSize',18);
            rc_lgndrunline = plot([0 1], [j+4-1 j+4], '--k', 'LineWidth', 2);
            hold on
            % Variable legend
            rc_lgndrun.(['line_' num2str(j+4)]) = legend(rc_lgndrunline, {''});
            % Set the legend position
            rc_lgndrun.(['line_' num2str(j+4)]).Position = ...
                rc_lgnd_pos + [0 0 0 0] + (j+4-1) .* rc_lgnd_dpos;
            % Turn the legend box off
            rc_lgndrun.(['line_' num2str(j+4)]).Box = 'off';
            % Add the run name to the legend
        %         annotation('textbox', rc_lgnd_pos + [0.1 -0.4 0 0] + (j-1) .* rc_lgnd_dpos, ...
        %                 'string', [check.ts_var{j,1} '-' check.ts_var{j,3}], ...
        %                 'EdgeColor','none', 'FontSize',16, 'FontWeight','normal')
            annotation('textbox', rc_lgnd_pos + [0.1 0.0 0 0] + (j+4-1) .* rc_lgnd_dpos, ...
                    'string', 'NH_3 detection limit', ...
                    'EdgeColor','none', 'FontSize',16, 'FontWeight','normal')
        end
    end
    
    if check.ts_tbuffer(2) > check.ts_txlim(1) || ...
            check.ts_tbuffer(1) < check.ts_txlim(2)
        
        subplot('Position',[-1 -1 0.8 0.5])
        hold on
        rc_lgndrunline = patch( [0 1 1 0], [4 4 5 5], 1, ...
            'FaceColor','k', 'EdgeColor','none', 'FaceAlpha', 0.1);
        
        rc_lgndrun.(['line_' num2str(j+1)]) = legend(rc_lgndrunline, {''});
        rc_lgndrun.(['line_' num2str(j+1)]).Position = ...
            rc_lgnd_pos + [0 0 0 0] + (j+1-1) .* rc_lgnd_dpos;
        % Turn the legend box off
        rc_lgndrun.(['line_' num2str(j+1)]).Box = 'off';
        annotation('textbox', rc_lgnd_pos + [0.1 0.0 0 0] + (j+1-1) .* rc_lgnd_dpos, ...
                'string', 'Buffer phase', ...
                'EdgeColor','none', 'FontSize',16, 'FontWeight','normal')
    end
    ''
    
    
%     rc_temp = NaN(size(rc_data{1}(:,1)));
%     for i = 2 : length(rc_data{4}(:,2))
%         rc_INtemp = find(rc_data{1}(:,1) > rc_data{4}(i-1,1) & rc_data{1}(:,1) <= rc_data{4}(i,1));
%         rc_temp(rc_INtemp) = [rc_data{4}(i-1,2) : (rc_data{4}(i,2) - rc_data{4}(i-1,2)) / (length(rc_INtemp)-1) : rc_data{4}(i,2)];
%     end
%     rc_INtemp_t = find(rc_data{1}(:,1) > 14 & rc_data{1}(:,1) <= 17);
%     
%     figure
%     hold on; grid on
%     plot(rc_data{1}(rc_INtemp_t,1), rc_data{1}(rc_INtemp_t,2))
%     plot(rc_data{1}(rc_INtemp_t,1), rc_temp(rc_INtemp_t))
%     plot(rc_data{1}(rc_INtemp_t,1), rc_temp(rc_INtemp_t)-rc_data{1}(rc_INtemp_t,2))
%     plot(rc_data{1}(rc_INtemp_t,1), ones(size(rc_temp(rc_INtemp_t))) .* mean(rc_temp(rc_INtemp_t)-rc_data{1}(rc_INtemp_t,2)))
        
end



%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%  Plot first receptor statistics   %%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function f__FIGtimeseries_stats(outp, check, rs_position)
    
    
    
    
    rs_colorset = {[230,097,001]./255, [253,184,099]./255, ...
                   [094,060,153]./255, [178,171,210]./255};
    rs_xlabel = 'CEST';
    rs_ylabel = check.rs_ylabel;
    rs_xpos = outp.(['set' check.rs_var{1}]).recp_indiv(1,1) - 500;
    rs_title = ['Simulated NH_{3,total} at ' num2str(round(rs_xpos)) ' m distance from the emission source'];
    rs_xtick = [0:0.25:24];
    rs_xticklabel = cell(size(rs_xtick));
    rs_xtick_dt = 0.5;
    for i = 1 : length(rs_xticklabel)
        if ismember(rs_xtick(i), [0 : rs_xtick_dt : 24])
            rs_hour = floor(rs_xtick(i));
            rs_min = round((rs_xtick(i) - rs_hour) * 60);
            
            rs_hour = num2str(rs_hour);
            if length(rs_hour) == 1
                rs_hour = ['0' rs_hour];
            end
            rs_min = num2str(rs_min);
            if length(rs_min) == 1
                rs_min = ['0' rs_min];
            end
            rs_xticklabel{i} = [rs_hour ':' rs_min];
        end
    end
    clear rs_hour rs_min
    
    rs_section = outp.(['set' check.rs_var{1}]).recp_indiv_file{1}(9:10);
    if strcmp(rs_section, 'xy')
        rs_INpos = [find( outp.(['set' check.rs_var{1}]).xt == outp.(['set' check.rs_var{1}]).recp_indiv(1,1) ), ...
            find( outp.(['set' check.rs_var{1}]).yt == outp.(['set' check.rs_var{1}]).recp_indiv(1,2) )];
    elseif strcmp(rs_section, 'xz')
        rs_INpos = [find( outp.(['set' check.rs_var{1}]).xt == outp.(['set' check.rs_var{1}]).recp_indiv(1,1) ), ...
            find( outp.(['set' check.rs_var{1}]).zt == outp.(['set' check.rs_var{1}]).recp_indiv(1,3) )];
    elseif strcmp(rs_section, 'yz')
        rs_INpos = [find( outp.(['set' check.rs_var{1}]).yt == outp.(['set' check.rs_var{1}]).recp_indiv(1,2) ), ...
            find( outp.(['set' check.rs_var{1}]).zt == outp.(['set' check.rs_var{1}]).recp_indiv(1,3) )];
    end
    rs_INt = [1: 1 : length(outp.(['set' check.rs_var{1}]).time )]';

    rs_time_temp = outp.(['set' check.rs_var{1}]).time(rs_INt);
    rs_matfile = matfile([check.(['set' check.rs_var{1}]).loc ...
         outp.(['set' check.rs_var{1}]).recp_indiv_file{1} check.rs_var{2}]);
     

     % Make sure the variable is an nh3 variable
    rs_setname = check.rs_var{2};
    if strcmp(rs_setname(1:3), 'nh3') == 0
        rs_setname = outp.((['set' check.rs_var{1}])).INFO.Name{13};
    end
    % Find the name of the scalar set
    if strcmp(rs_setname(end - length(check.variable_plume_indicator)+1: end), ...
            check.variable_plume_indicator)

        rs_setname(end - length(check.variable_plume_indicator)+1: end) = '';
    elseif strcmp(rs_setname(end - length(check.variable_background_indicator)+1: end), ...
            check.variable_background_indicator)

        rs_setname(end - length(check.variable_background_indicator)+1: end) = '';
    end
    if strcmp(rs_setname(4), '_')
        rs_setname = rs_setname(5:end);
    else
        rs_setname = rs_setname(4:end);
    end
    % Make an exception for nhr_r0b, since set r0 does not have a plume
    if strcmp(check.rs_var{2}, 'nh3_r0b')
        rs_varplume = 'nh3_r1p';
    end
    rs_varplume = ['nh3_' rs_setname check.variable_plume_indicator];
    rs_matfile_plume = matfile([check.(['set' check.rs_var{1}]).loc ...
         outp.(['set' check.rs_var{1}]).recp_indiv_file{1} rs_varplume '.mat']);
     
     
     
    % Time averageing data
    rs_time = rs_time_temp;
    if round(check.rs_var{4}, 6) > round(median( rs_time_temp(2:end) - rs_time_temp(1:end-1) ), 6)

        rs_time = [rs_time_temp(1) - median( rs_time_temp(2:end) - rs_time_temp(1:end-1) ) +  check.rs_var{4} ...
            : check.rs_var{4} : rs_time_temp(end)]';
        rs_C_temp = NaN(length(rs_time),1);

        if length(rs_time) > 1
            for k = 1 : length(rs_time)
                if k == 1 
                    rs_C_temp(k) = mean( rs_matfile.(check.rs_var{2})(rs_INpos(1), rs_INpos(2), ...
                        rs_INt(rs_time_temp <= rs_time(k)))  );
                else
                    rs_C_temp(k) = mean( rs_matfile.(check.rs_var{2})(rs_INpos(1), rs_INpos(2), ...
                        rs_INt(rs_time_temp > rs_time(k-1) & rs_time_temp <= rs_time(k)))  );
                end
            end
        end
        rs_C = rs_C_temp;

    elseif round(check.rs_var{4}, 6) == round(median( rs_time_temp(2:end) - rs_time_temp(1:end-1) ), 6)

        rs_C = squeeze(rs_matfile.(check.rs_var{2})(rs_INpos(1), rs_INpos(2), rs_INt));
    else
        rs_C = NaN(size(rs_time));
    end
    rs_INt_analysis = find(rs_time > check.rs_tanalysis(1) & ...
        rs_time <= check.rs_tanalysis(2));
    
    % Calculate moving average concentration
    rs_averaging_time_h = 1;
    rs_Cma = movmean(rs_C, [rs_averaging_time_h / (rs_time(4) - rs_time(3)), 0]);
    
    % Time averageing PLUME data
    rs_time = rs_time_temp;
    if round(check.rs_var{4}, 6) > round(median( rs_time_temp(2:end) - rs_time_temp(1:end-1) ), 6)

        rs_time = [rs_time_temp(1) - median( rs_time_temp(2:end) - rs_time_temp(1:end-1) ) +  check.rs_var{4} ...
            : check.rs_var{4} : rs_time_temp(end)]';
        rs_Cplume = NaN(length(rs_time),1);

        if length(rs_time) > 1
            for k = 1 : length(rs_time)
                if k == 1 
                    rs_Cplume(k) = mean( rs_matfile_plume.(rs_varplume)(rs_INpos(1), rs_INpos(2), ...
                        rs_INt(rs_time_temp <= rs_time(k)))  );
                else
                    rs_Cplume(k) = mean( rs_matfile_plume.(rs_varplume)(rs_INpos(1), rs_INpos(2), ...
                        rs_INt(rs_time_temp > rs_time(k-1) & rs_time_temp <= rs_time(k)))  );
                end
            end
        end

    elseif round(check.rs_var{4}, 6) == round(median( rs_time_temp(2:end) - rs_time_temp(1:end-1) ), 6)

        rs_Cplume = [squeeze(rs_matfile_plume.(rs_varplume)(rs_INpos(1), rs_INpos(2), rs_INt))];
    else
        rs_Cplume = NaN(size(rs_time));
    end
    
    % Calculate C - Cma
    rs_Cstats = NaN(size(rs_C));
    rs_Cstats(rs_INt_analysis) = rs_C(rs_INt_analysis) - rs_Cma(rs_INt_analysis);
    
    % Calculate the statistics
    rs_mean = mean(rs_C(rs_INt_analysis));
    rs_mean_ma = mean(rs_Cstats(rs_INt_analysis));
    rs_std = std(rs_C(rs_INt_analysis));
    rs_std_ma = std(rs_Cstats(rs_INt_analysis));
    rs_fI = rs_std_ma / rs_mean;
    rs_I = length( find(rs_Cplume(rs_INt_analysis) >= check.rs_Ithreshhold) ) ...
        ./ length(rs_INt_analysis);
    
    rs_yticks = cell(length(check.rs_Cticks_num),1);
    for i = 1 : length(check.rs_Cticks_num)
        if ismember( check.rs_Cticks_num(i), check.rs_Cticks_str )
            rs_yticks{i} = num2str(check.rs_Cticks_num(i));
        end
    end
    
    
    %% Build the figure
    figure('units','pixels', 'Color','w',...
        'innerposition', rs_position, ...
        'Name', ['t_window = ' num2str(check.rs_txlim(1)) '-' num2str(check.rs_txlim(2))]);
    rs_sub = subplot('Position',[0.08 0.46 0.68 0.535]);
%     rs_toplgnd = gobjects(1,7);
    rs_toplgnd = gobjects(1,5);
    
    hold on         % show multiple lines
    grid on
    rs_ax1 = gca;   % axis variable
    
    for i = 1 : length(rs_xticklabel)
        if ismember(rs_xtick(i), [0 : rs_xtick_dt : 24])
            plot([rs_xtick(i) rs_xtick(i)], [-9999 9999], 'Color',[0 0 0 0.25], 'LineWidth', 1.2)
        end
    end
    for i = 1 : length(check.rs_Cticks_str)
        plot([-9999 9999], [check.rs_Cticks_str(i) check.rs_Cticks_str(i)], ...
            'Color',[0 0 0 0.20], 'LineWidth', 1)
    end
    
    % Plot the simulated measurements time series
    rs_toplgnd(1) = plot(rs_time, rs_C, 'Color',rs_colorset{1}, ...
        'LineStyle', '-', 'LineWidth',3);
%     % Plot the plume time series
%     rs_toplgnd(5) = plot(rs_time, rs_Cplume, 'Color',rs_colorset{2}, ...
%         'LineStyle', '-', 'LineWidth',3);
%     % Plot the intermittency threshold (NH3 detection limit)
%     rs_toplgnd(6) = plot([rs_time(rs_INt_analysis(1)) rs_time(rs_INt_analysis(end))], [check.rs_Ithreshhold check.rs_Ithreshhold], 'Color','k', ...
%         'LineStyle', ':', 'LineWidth',3);
%     % Plot the Moving average time series
%     rs_toplgnd(2) = plot(rs_time(rs_INt_analysis), rs_Cma(rs_INt_analysis), 'Color',rs_colorset{3}, ...
%         'LineStyle', '-', 'LineWidth',3);
%     % Plot the analysis window average
%     rs_toplgnd(3) = plot([rs_time(rs_INt_analysis(1)) rs_time(rs_INt_analysis(end))], ...
%         [mean(rs_C(rs_INt_analysis)), mean(rs_C(rs_INt_analysis))], ...
%         'Color',rs_colorset{3}, 'LineStyle', ':', 'LineWidth',3);
%     % Legend plot standard deviation
%     rs_toplgnd(4) = plot([-999 -990], [-999 -990], 'Color',rs_colorset{4}, ...
%         'LineStyle', '-', 'LineWidth',3);
%     % Show buffer phae
%     rs_toplgnd(7) = patch([check.rs_tbuffer(1) check.rs_tbuffer(2) check.rs_tbuffer(2) check.rs_tbuffer(1)], ...
%         [-9999 -9999 9999 9999], [1 1 1 1], 'FaceColor', 'k', 'FaceAlpha',0.1, 'EdgeColor','none');
    % Plot the plume time series
    rs_toplgnd(3) = plot(rs_time, rs_Cplume, 'Color',rs_colorset{2}, ...
        'LineStyle', '-', 'LineWidth',3);
    % Plot the intermittency threshold (NH3 detection limit)
    rs_toplgnd(5) = plot([rs_time(rs_INt_analysis(1)) rs_time(rs_INt_analysis(end))], [check.rs_Ithreshhold check.rs_Ithreshhold], 'Color','k', ...
        'LineStyle', ':', 'LineWidth',3);
    % Plot the Moving average time series
    rs_toplgnd(2) = plot(rs_time(rs_INt_analysis), rs_Cma(rs_INt_analysis), 'Color',rs_colorset{3}, ...
        'LineStyle', '-', 'LineWidth',3);
    % Plot the analysis window average
    % Legend plot standard deviation
    rs_toplgnd(4) = plot([-999 -990], [-999 -990], 'Color',rs_colorset{4}, ...
        'LineStyle', '-', 'LineWidth',3);
%     % Show buffer phae
%     rs_toplgnd(6) = patch([check.rs_tbuffer(1) check.rs_tbuffer(2) check.rs_tbuffer(2) check.rs_tbuffer(1)], ...
%         [-9999 -9999 9999 9999], [1 1 1 1], 'FaceColor', 'k', 'FaceAlpha',0.1, 'EdgeColor','none');
    
    
    rs_yticks = cell(length(check.rs_Cticks_num),1);
    for i = 1 : length(check.rs_Cticks_num)
        if ismember( check.rs_Cticks_num(i), check.rs_Cticks_str )
            rs_yticks{i} = num2str(check.rs_Cticks_num(i));
        end
    end

    set(rs_ax1,'FontSize',check.figure_font, 'XLim',check.rs_txlim, 'YLim',check.rs_Clim, ...
        'YTick',check.rs_Cticks_num, 'YTickLabel',rs_yticks, ...
        'XTick', rs_xtick, 'XTickLabel',cell(size(rs_xtick)))
    ylabel(rs_ylabel);
%     annotation('textbox',[0.08     0.96   0.68    0.05], 'String', rs_title, ...
%         'FontSize', 18, 'EdgeColor', 'none', 'FontWeight','Bold', 'HorizontalAlignment','center') 
%     rs_legend = legend(rs_toplgnd, ...
%         {check.rs_var{3}, [check.rs_var{3} '_{ MA}'],'NH_3', '\sigma' 'NH_{3,plume}','Detection limit','Buffer phase'});
    rs_legend = legend(rs_toplgnd, ...
        {check.rs_var{3}, [check.rs_var{3} '_{ MA}'], 'NH_{3,plume}', '\sigma_{detrend}', 'NH_3 Detection limit'});
    rs_legend.Box = 'off';
    rs_legend.Position =  [0.789     0.49   0.19    0.5];
%     temp_line = annotation('line','Position',[0.805, 0.40, 0.082 0.001],'LineWidth',1.2);
%     temp_line.Position = [0.829, 0.82, 0.038 0.001];
    
    
    rs_sub = subplot('Position',[0.08 0.135 0.68 0.31]);
    hold on         % show multiple lines
    grid on
    rs_ax2 = gca;   % axis variable
    set(rs_ax2,'FontSize',check.figure_font, 'XLim',check.rs_txlim, ...
        'Ylim', [floor(min(rs_Cstats(rs_INt_analysis))) ceil(max(rs_Cstats(rs_INt_analysis)))], ...
        'XTick', rs_xtick, 'XTickLabel', rs_xticklabel)
    for i = 1 : length(rs_xticklabel)
        if ismember(rs_xtick(i), [0 : rs_xtick_dt : 24])
            plot([rs_xtick(i) rs_xtick(i)], [-9999 9999], 'Color',[0 0 0 0.25], 'LineWidth', 1.2)
        end
    end
    plot([-9999 9999], [0 0], 'Color',[0 0 0 0.20], 'LineWidth', 1)
    set(rs_ax2, 'YTick', [-2 : 1 : 4], 'YTickLabel', {'-2','','0','','2','',''});
        
%     patch([check.rs_tbuffer(1) check.rs_tbuffer(2) check.rs_tbuffer(2) check.rs_tbuffer(1)], ...
%         [-9999 -9999 9999 9999], [1 1 1 1], 'FaceColor', 'k', 'FaceAlpha',0.1, 'EdgeColor','none');
    plot(rs_time, rs_Cstats, 'Color',rs_colorset{1}, 'LineStyle', '-', 'LineWidth',3);
    plot([rs_time(rs_INt_analysis(1)) rs_time(rs_INt_analysis(end))], [rs_mean_ma rs_mean_ma], ...
        'Color',rs_colorset{4}, 'LineStyle', ':', 'LineWidth',3);
    plot([rs_time(rs_INt_analysis(1)) rs_time(rs_INt_analysis(end))], [rs_mean_ma-rs_std_ma rs_mean_ma-rs_std_ma], ...
        'Color',rs_colorset{4}, 'LineStyle', '-', 'LineWidth',3);
    plot([rs_time(rs_INt_analysis(1)) rs_time(rs_INt_analysis(end))], [rs_mean_ma+rs_std_ma rs_mean_ma+rs_std_ma], ...
        'Color',rs_colorset{4}, 'LineStyle', '-', 'LineWidth',3);
    xlabel(rs_xlabel);
    ylabel('NH_{3,detrend} [ppb]')
    
    annotation('textbox',[0.78     0.125   0.210    0.286], 'String', ...
        '', 'FontSize', check.figure_font, 'Color','none', 'EdgeColor', 'k')
    annotation('textbox',[0.78     0.338   0.21    0.071], 'String', ...
        ['NH_{3,total} = ' num2str(round(rs_mean,1)) ' ppb'], 'FontSize', check.figure_font, ...
        'EdgeColor', 'none')
    annotation('line','Position',[0.787 0.398, 0.04 0.001],'LineWidth',1.2)
    annotation('textbox',[0.78     0.267   0.21    0.071], 'String', ...
        ['\sigma = ' num2str(round(rs_std_ma,3)) ' ppb'], 'FontSize', check.figure_font, ...
        'EdgeColor', 'none')
%     annotation('textbox',[0.78     0.196   0.21    0.071], 'String', ...
%         ['\sigma/NH_{3,total} = ' num2str(round(rs_fI,3))], 'FontSize', check.figure_font, ...
%         'EdgeColor', 'none')
%     annotation('line','Position',[0.85, 0.251, 0.043 0.001],'LineWidth',1.2)
    annotation('textbox',[0.78     0.196   0.21    0.071], 'String', ...
        ['fI = ' num2str(round(rs_fI,3))], 'FontSize', check.figure_font, 'EdgeColor', 'none')
    annotation('textbox',[0.78     0.125   0.21    0.071], 'String', ...
        ['I = ' num2str(round(rs_I,3))], 'FontSize', check.figure_font, 'EdgeColor', 'none')
    
    'Done'
    
    %% End of function
end


%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%      Show receptor positions      %%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function f__FIGreceptor_xypos(outp, check, si_setnr)
    
    st_mean = outp.time(4) - (outp.time(8)-outp.time(4));
    st_mean = [st_mean + 1, st_mean+2];
    if outp.time(end) < st_mean(2)
        st_mean = [outp.time(1), outp.time(end)];
    end
    

    %%%%%%
    % SHOW LOCATION OF ALL RECEPTORS + HIGHLIGHT check.INDIVreceptor_IN
    %%%%%%

    si_colorset = {[230,097,001]./255, [094,060,153]./255,...
                   [253,184,099]./255, [178,171,210]./255};

    si_colormap_Delta = 1000;
    si_colormap = [[ [094 : (255 - 094)/(si_colormap_Delta / 2) : 255 ]', ...
                     [060 : (255 - 060)/(si_colormap_Delta / 2) : 255 ]', ...
                     [153 : (255 - 153)/(si_colormap_Delta / 2) : 255 ]'] ...
                   ./255 ; ...
                   [ [255-(255 - 230)/(si_colormap_Delta / 2) : -(255 - 230)/(si_colormap_Delta / 2) : 230 ]', ...
                     [255-(255 - 097)/(si_colormap_Delta / 2) : -(255 - 097)/(si_colormap_Delta / 2) : 097 ]', ...
                     [255-(255 - 001)/(si_colormap_Delta / 2) : -(255 - 001)/(si_colormap_Delta / 2) : 001 ]'] ./255];

    % Define figure title
    si_title = [si_setnr ': The xy-mean ' check.ind_var ...
        ' at ' num2str(outp.recp_indiv(1,3)) ' m (' num2str(st_mean(1)) ...
        ' - ' num2str(st_mean(2)) ' CEST)'];
    % Define figure C label
    si_label = [check.ind_var ' [' outp.INFO.Unit{ ...
         ismember(outp.INFO.Name, {check.ind_var}) } ']'];

     %%%%% Define the variable to be plotted
    % Define the x and y axis data
    si_temp_x = outp.xm;
    si_temp_dx = si_temp_x(2) - si_temp_x(1);
    si_temp_y = outp.ym;
    si_temp_dy = si_temp_y(2) - si_temp_y(1);
    si_temp_source = outp.srcpos.(check.ind_var);
    si_temp_recep = zeros(length(outp.xt), length(outp.yt));
    si_temp_recep(ismember(outp.xt, outp.recp_indiv(:,1)) , ...
                  ismember(outp.yt, outp.recp_indiv(:,2)) ) = 1;
    
    si_matfile = matfile([check.(si_setnr).loc outp.recp_indiv_file{1} check.ind_var '.mat']);
    si_temp_c = squeeze( nanmean( si_matfile.(check.ind_var)(:,:, ...
        find(outp.time > st_mean(1) & outp.time <= st_mean(2)) ) ,3) ); 
    si_temp_c(si_temp_c<0) = 0;

    [si_crossX, si_crossY] = meshgrid(si_temp_x , si_temp_y);
    si_crossC = NaN(size(si_crossX));
    si_crossXgrid = NaN(4, length(si_temp_x) * length(si_temp_y));
    si_crossYgrid = NaN(4, length(si_temp_x) * length(si_temp_y));
    si_crossCgrid = NaN(length(si_temp_x) * length(si_temp_y),1);
    si_crossS = NaN(length(si_temp_x) * length(si_temp_y),1);
    si_crossR = NaN(length(si_temp_x) * length(si_temp_y),1);
    for j = 1 : length(si_temp_x)
        for k = 1 : length(si_temp_y)
            si_crossXgrid(:,k + (j-1)*length(si_temp_y)) = ...
                [si_temp_x(j); si_temp_x(j)+si_temp_dx; si_temp_x(j)+si_temp_dx; si_temp_x(j)];
            si_crossYgrid(:,k + (j-1)*length(si_temp_y)) = ...
                [si_temp_y(k); si_temp_y(k); si_temp_y(k)+si_temp_dy; si_temp_y(k)+si_temp_dy];


            si_crossC( (si_crossX == si_temp_x(j) & si_crossY == si_temp_y(k)) ) = si_temp_c(j,k);
            si_crossCgrid(k + (j-1)*length(si_temp_y)) = si_temp_c(j,k);

            si_crossS(k + (j-1)*length(si_temp_y)) = si_temp_source(j,k);
            si_crossR(k + (j-1)*length(si_temp_y)) = si_temp_recep(j,k);


        end
    end
    si_crossS(si_crossS == 0) = NaN;
    
    if isnan( check.ind_Clim(1) )
        check.ind_Clim(1) = min(si_crossC(:));
    end
    if isnan( check.ind_Clim(2) )
        check.ind_Clim(2) = max(si_crossC(:));
    end

    if strcmp(check.ind_xylog, 'yes') 
        if check.ind_Clim(1) == 0
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            fprintf(['\n  !!! ERROR !!!\nThe xy-positions of the receptors will not be shown.\n'...
                'With a log scale, "check.ind_Clim(1)" cannot be zero. \n\n'])      
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            return
        end
        % Figure color ticks
        if sum(isnan(check.ind_Clim)) == 0
            si_Cminmax = check.ind_Clim;
        else
            si_Cminmax = [min(si_crossC(:)), max(si_crossC(:))];
        end
        if min(si_crossC(:)) == max(si_crossC(:))
            si_Cminmax(2) = si_Cminmax(2)+0.1;
        end
        if min(si_crossC(:)) < check.ind_Clim(1)
            si_Cminmax(1) = check.ind_Clim(1);
        end
        if min(si_crossC(:)) > check.ind_Clim(2)
            si_Cminmax = [min(si_crossC(:)), max(si_crossC(:))];
        end
        si_ticklabelset = [ check.ind_Clim(1), 0 : check.ind_Clim(1)*5 : check.ind_Clim(1)*10];
        
        si_tickset = [0 : check.ind_Clim(1) : check.ind_Clim(1)*10];
        si_tickcounter = 0;
        while si_tickcounter == 0
            if si_tickset(end) >= 100 %check.ind_Clim(2)
                si_tickcounter = 1;
            else
                si_ticklabelset = [si_ticklabelset, 0 : si_tickset(end)*5 : si_tickset(end)*10];
                si_tickset = [si_tickset, 0 : si_tickset(end) : si_tickset(end)*10];
            end
        end
        si_temp_TickNum = unique(si_tickset);
        si_temp_TickLim = {si_temp_TickNum(si_temp_TickNum < si_Cminmax(1)) ; ...
                           si_temp_TickNum(si_temp_TickNum > si_Cminmax(2))};
        si_temp_TickLim = [si_temp_TickLim{1}(end), si_temp_TickLim{2}(1)];
        si_temp_TickNum = si_temp_TickNum(si_temp_TickNum >= si_temp_TickLim(1) & si_temp_TickNum <= si_temp_TickLim(2) );
        if length(si_temp_TickNum) <= 10
            si_ticklevels = NaN(((length(unique(si_tickset))-2)/9),1);
            for j = 1 : length(si_ticklevels)
                si_ticklevels(j) = check.ind_Clim(1)* 10^(j-1);
            end
            si_ticklevel = 1;
            si_counter = 0;
            while si_counter == 0
                if si_Cminmax(1) / si_ticklevels(si_ticklevel) >= 1 && ...
                        si_Cminmax(1) / si_ticklevels(si_ticklevel) <= 10
                    si_counter = 1;
                else
                    si_ticklevel = si_ticklevel + 1;
                end
            end
            si_temp_TickNum = [si_ticklevels(si_ticklevel) * floor(si_Cminmax(1) / si_ticklevels(si_ticklevel)) ...
                : si_ticklevels(si_ticklevel)/10 : ...
                si_ticklevels(si_ticklevel) * ceil(si_Cminmax(2) / si_ticklevels(si_ticklevel))];

            si_INmin = find( si_temp_TickNum - min(si_crossC(:)) > 0 );
            si_INmin = si_INmin(1) - 1;
            si_INmax = find( si_temp_TickNum - max(si_crossC(:)) > 0 );
            si_INmax = si_INmax(1);
            si_temp_TickNum = si_temp_TickNum(si_INmin:si_INmax);
            si_Cminmax = [si_temp_TickNum(1), si_temp_TickNum(end)];

            temp_TickLabelsStr = cell(length(si_temp_TickNum),1);
            temp_INlabelstr = [];
            if length(si_temp_TickNum) <= 11
                for j = 1 : length(si_temp_TickNum)
                    temp_TickLabelsStr{j} = num2str(si_temp_TickNum(j));
                    temp_INlabelstr = [temp_INlabelstr; j];
                end
            elseif length(si_temp_TickNum) > 11 && length(si_temp_TickNum) <= 21
                for j = 1 : 2 : length(si_temp_TickNum)
                    temp_TickLabelsStr{j} = num2str(si_temp_TickNum(j));
                    temp_INlabelstr = [temp_INlabelstr; j];
                end
            else
                for j = 1 : 5 : length(si_temp_TickNum)
                    temp_TickLabelsStr{j} = num2str(si_temp_TickNum(j));
                    temp_INlabelstr = [temp_INlabelstr; j];
                end
            end

            if length(temp_TickLabelsStr) <= 11
                si_temp_TickLabelsNum = si_temp_TickNum(...
                unique([1:1:floor(length(temp_TickLabelsStr)/2),  ...
                floor(length(temp_TickLabelsStr)/2) : 2 : length(si_temp_TickNum)]));
            elseif length(temp_TickLabelsStr) <= 21
                si_temp_TickLabelsNum = si_temp_TickNum(...
                    unique([1:1:temp_INlabelstr(3),  ...
                    temp_INlabelstr(3) : 4 : length(si_temp_TickNum)]));
            else
                si_temp_TickLabelsNum = si_temp_TickNum(...
                    unique([1:1:temp_INlabelstr(2),  ...
                    temp_INlabelstr(2) : 5 : length(si_temp_TickNum)]));
            end

        else
            if length(si_temp_TickNum) <= 15
                temp_TickLabelsStr = cell(length(si_temp_TickNum),1);
                for j = 1 : length(si_temp_TickNum)
                    temp_TickLabelsStr{j} = num2str(si_temp_TickNum(j));
                end
                si_temp_TickLabelsNum = si_temp_TickNum;
            else
                si_temp_TickLabelsNum = unique(si_ticklabelset);
                si_temp_TickLabelsNum = si_temp_TickLabelsNum(ismember(si_temp_TickLabelsNum, si_temp_TickNum));
                temp_TickLabelsStr = cell(length(si_temp_TickNum),1);
                for j = si_temp_TickLabelsNum
                    temp_TickLabelsStr{si_temp_TickNum == j} = num2str(j);
                end
            end
        end
    end
    


    subpos_ymax = 620;
    subpos_xmax = 1300;
    subpos_Csize = 100;  % Width required for the colorbar

    %%%%% Determine the subplot sizes
    % Define the model domain
    dom_x = outp.xm(end);
    dom_y = outp.ym(end);
    xy_dom_dx = subpos_xmax / dom_x;
    xy_dom_dy = subpos_ymax / dom_y;
    xy_posx = xy_dom_dy * dom_x;
    xy_posy = subpos_ymax;
    % Determine the size of the subplot 
    si_fig_subsize = [100 80 xy_posx xy_posy]; 

    % Definitive position and size of the figure
    si_fig_position = [0, 40, ...
        si_fig_subsize(3)+si_fig_subsize(1)+50+subpos_Csize, ...
        si_fig_subsize(4)+si_fig_subsize(2)+25+20];

    % Definitive position of the subplot in percentages (betweem 0 and 1)
    si_fig_subposition = [si_fig_subsize(1)/si_fig_position(3), ... 
                      si_fig_subsize(2)/si_fig_position(4), ...
                      si_fig_subsize(3)/si_fig_position(3), ...
                      si_fig_subsize(4)/si_fig_position(4)];


    si_fig = figure('units','pixels', 'Color','w',...
        'innerposition', si_fig_position);
    % Build the plotting space
    si_sub = subplot('Position',si_fig_subposition);
    hold on; 
    % Plot cross section
    si_gridplot = patch(si_crossXgrid, si_crossYgrid, si_crossCgrid);
    si_gridplot.LineWidth = 0.01;
    si_gridplot.EdgeAlpha = 0.2;
    % Colorbar
    si_cbar = colorbar;
    si_cbar.Position = [(si_fig_subsize(1) + si_fig_subsize(3) + 10) / si_fig_position(3), ...
                     (si_fig_subsize(2)) / si_fig_position(4), ...
                     (30) / si_fig_position(3), ...
                     (si_fig_subsize(4)) / si_fig_position(4)];
    si_cbar.Label.String = si_label;
    si_cbar.LineWidth = 1;
    si_cbar.TickLength = 0.02;
    si_gca = gca;
    if strcmp(check.ind_xylog, 'yes')
        set(si_gca, 'ColorScale','log', 'FontSize', 18)
        % Contour lines
        si_contour = contour(si_crossX, si_crossY, si_crossC, ...
            si_temp_TickLabelsNum, 'LineColor','k', 'LineWidth',1.7, ...
            'ShowText','on', 'TextListMode','manual', 'TextList', si_temp_TickLabelsNum);
    else
        set(si_gca, 'FontSize', 18)
        si_contour = contour(si_crossX, si_crossY, si_crossC, ...
            'LineColor','k', 'LineWidth',1.7, 'ShowText','on');
    end
    % Plot the location of the source and the receptors
    si_sourceplot = patch(si_crossXgrid(:,si_crossS == 1), si_crossYgrid(:,si_crossS == 1), si_colorset{1});
    si_sourceplot.EdgeColor = 'none';
    si_sourceplot.FaceColor = 'r';
    si_sourceplot.FaceAlpha = 1.00;
    if size(outp.recp_indiv,1) > 1

        si_recepplot = patch(si_crossXgrid(:,si_crossR == 1), si_crossYgrid(:,si_crossR == 1), si_colorset{1});
        si_recepplot.EdgeColor = 'none';
        si_recepplot.FaceColor = 'g';
        si_recepplot.FaceAlpha = 1.0;
    end
    % Settings
    si_gca = gca;
    xlabel('x [m]'); ylabel('y [m]')
    xlim([min(si_crossXgrid(:)) max(si_crossXgrid(:))])
    ylim([min(si_crossYgrid(:)) max(si_crossYgrid(:))])
    colormap(si_colormap);  
    if strcmp(check.ind_xylog, 'yes')
        
        si_cbar.Ticks = si_temp_TickNum;
        si_cbar.TickLabels = temp_TickLabelsStr;
        if si_temp_TickNum(1) == 0
            si_cbar.Limits = [si_temp_TickNum(2), si_temp_TickNum(end)];
            si_sub.CLim = [si_temp_TickNum(2), si_temp_TickNum(end)];
        else
            si_cbar.Limits = [si_temp_TickNum(1), si_temp_TickNum(end)];
            si_sub.CLim = [si_temp_TickNum(1), si_temp_TickNum(end)];
        end
        
    else
        
        si_cbar.Limits = si_sub.CLim;
    end
    annotation('textbox', [0.01 (si_fig_subsize(2)+si_fig_subsize(4))/si_fig_position(4), 0.98 0.01], ...
                'string', si_title, ...
                'EdgeColor','none', 'FontSize',18, 'FontWeight','normal', ...
                'HorizontalAlignment','center','VerticalAlignment','bottom', ...
                'FontWeight','Bold');
                
end 










%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%% Plot first receptor time series + cross-section %%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function f__FigReceptor_crosstime(outp, check, rc_var, rc_position)
    
    rc_colorset = {[230,097,001]./255, [094,060,153]./255,...
                   [253,184,099]./255, [178,171,210]./255};
    rc_linestyle = {'-',':','--','-.'};
    
    
    %% prepare data for the time series figure
    
    rc_plot = struct();
    rc_INt = find(outp.time > check.int_twindow(1) & outp.time <= check.int_twindow(2));
    rc_t = outp.time(rc_INt);
    % Preallocate the data
    rc_dt = median(rc_t(2:end) - rc_t(1:end-1));
    % Define matfile
    if strcmp(check.int_xy_xz, 'xy')
        rc_INx = find(round(outp.xt, 6) == round(outp.recp_indiv(1,1), 6));
        rc_INy = find(round(outp.yt, 6) == round(outp.recp_indiv(1,2), 6));
        
        rc_matfile = matfile([outp.loc check.int_var{1} 'crossxy' ...
            outp.cross_lvl{ round(outp.zpos, 6) == round(outp.recp_indiv(1,3), 6) } ...
            '_' check.int_var{1,2} '.mat']);
    elseif strcmp(check.int_xy_xz, 'xz')
        rc_INx = find(round(outp.xt, 6) == round(outp.recp_indiv(1,1), 6));
        rc_INy = find(round(outp.zt, 6) == round(outp.recp_indiv(1,3), 6));

        rc_matfile = matfile([outp.loc check.int_var{1} 'crossxz_' check.int_var{1,2} '.mat']);
    elseif strcmp(check.int_xy_xz, 'yz')
        rc_INx = find(round(outp.yt, 6) == round(outp.recp_indiv(1,2), 6));
        rc_INy = find(round(outp.zt, 6) == round(outp.recp_indiv(1,3), 6));

        rc_matfile = matfile([outp.loc check.int_var{1} 'crossyz_' check.int_var{1,2} '.mat']);
    end
    rc_plot.matfile = rc_matfile;
    
    
    rc_tnew = [rc_t(1) - rc_dt + min(check.int_avgt) : min(check.int_avgt) : rc_t(end)]';  
    rc_plot.var = NaN(length(rc_tnew), length(check.int_avgt));
    
    for j = 1 : length(check.int_avgt)
        
        if round(check.int_avgt(j),8) >= round(rc_dt,8)
            
            rc_t_temp = [rc_t(1) - rc_dt + check.int_avgt(j) : ...
                check.int_avgt(j) : rc_t(end)]';

            for JJ = 1 : length(rc_t_temp)
                if JJ == 1
                    
                    rc_plot.var(round(rc_tnew, 8) == round(rc_t_temp(JJ), 8),j) = ...
                        mean( rc_matfile.(check.int_var{1,2})(rc_INx, rc_INy, ...
                        rc_INt( round(rc_t, 8) <=  round(rc_t_temp(JJ), 8) ) ) );
                else
                    rc_plot.var(round(rc_tnew, 8) == round(rc_t_temp(JJ), 8),j) = ...
                        mean( rc_matfile.(check.int_var{1,2})(rc_INx, rc_INy, ...
                        rc_INt( round(rc_t, 8) >  round(rc_t_temp(JJ-1), 8) &  ...
                        round(rc_t, 8) <=  round(rc_t_temp(JJ), 8) ) ) );
                end

            end
        elseif round(check.int_avgt(j),8) < round(rc_dt,8)
            
            check.int_avgt(j) = rc_dt;
            rc_t_temp = [rc_t(1) - rc_dt + check.int_avgt(j) : ...
                check.int_avgt(j) : rc_t(end)]';
            rc_plot.var(round(rc_t, 8) == round(rc_t_temp, 8),j) = ...
                    rc_matfile.(check.int_var{1,2})(rc_INx, rc_INy, :);
        end

    end
    rc_plot.time = rc_tnew;
    
    %% Prepare the cross-section figure
    
    % Set the x and y limits
    if sum(isnan(check.int_xlim)) == 0
        rc_plot.cross_xlim = check.int_xlim;
        rc_plot.cross_ylim = check.int_ylim;
    else 
        if max(outp.xt) > max(outp.yt)
            rc_plot.cross_xlim = [0 max(outp.yt)];
            rc_plot.cross_ylim = [0 max(outp.yt)];
        else
            rc_plot.cross_xlim = [0 max(outp.xt)];
            rc_plot.cross_ylim = [0 max(outp.xt)];
        end
    end
    if sum(isnan(check.int_Clim)) > 0
        rc_clim = [min(rc_plot.var(:)) max(rc_plot.var(:))];
        if rc_clim(1) < rc_clim(2) / 1000 && rc_clim(1) > -0.0001
            rc_clim(1) = 0;
        end
    else
        rc_clim = check.int_Clim;
    end
    % Define basic parameters
    rc_plot.receppos = [outp.recp_indiv(1,1), outp.recp_indiv(1,2), outp.recp_indiv(1,2)];
    if strcmp(check.int_xy_xz, 'xy')
        temp_dx = outp.xt(2) - outp.xt(1);
        temp_IN = find(outp.xt >= rc_plot.cross_xlim(1) & ...
            outp.xt <= rc_plot.cross_xlim(2));
        rc_plot.xm = outp.xm( temp_IN );
        rc_plot.xt = outp.xt( temp_IN );
        
        temp_dy = outp.yt(2) - outp.yt(1);
        temp_IN = find(outp.yt >= rc_plot.cross_ylim(1) & ...
            outp.yt <= rc_plot.cross_ylim(2));
        rc_plot.ym = outp.ym( temp_IN );
        rc_plot.yt = outp.yt( temp_IN );

    elseif strcmp(check.int_xy_xz, 'xz')
        temp_dx = outp.xt(2) - outp.xt(1);
        temp_IN = find(outp.xt >= rc_plot.cross_xlim(1) & ...
            outp.xt <= rc_plot.cross_xlim(2));
        rc_plot.xm = outp.xm( temp_IN );
        rc_plot.xt = outp.xt( temp_IN );
        
        temp_dy = outp.zt(2) - outp.zt(1);
        temp_IN = find(outp.zt >= rc_plot.cross_ylim(1) & ...
            outp.zt <= rc_plot.cross_ylim(2));
        rc_plot.ym = outp.zm( temp_IN );
        rc_plot.yt = outp.zt( temp_IN );

    elseif  strcmp(check.int_xy_xz, 'yz')
        temp_dx = outp.yt(2) - outp.yt(1);
        temp_IN = find(outp.yt >= rc_plot.cross_xlim(1) & ...
            outp.yt <= rc_plot.cross_xlim(2));
        rc_plot.xm = outp.ym( temp_IN );
        rc_plot.xt = outp.yt( temp_IN );
        
        temp_dy = outp.zt(2) - outp.zt(1);
        temp_IN = find(outp.zt >= rc_plot.cross_ylim(1) & ...
            outp.zt <= rc_plot.cross_ylim(2));
        rc_plot.ym = outp.zm( temp_IN );
        rc_plot.yt = outp.zt( temp_IN );
    end
    clear temp_IN
    
    % Restructure x, y and c for the cross section figure
%     temp_c = outp.CROSSSECTION.([rc_var check.int_xy_xz]);
    rc_plot.X = NaN(4, length(rc_plot.xm) * length(rc_plot.ym));
    rc_plot.Y = NaN(4, length(rc_plot.xm) * length(rc_plot.ym));
    for j = 1 : length(rc_plot.xm)
        for k = 1 : length(rc_plot.ym)
            rc_plot.X(:,k + (j-1)*length(rc_plot.ym)) = ...
                [rc_plot.xm(j); rc_plot.xm(j)+temp_dx; rc_plot.xm(j)+temp_dx; rc_plot.xm(j)];
            rc_plot.Y(:,k + (j-1)*length(rc_plot.ym)) = ...
                [rc_plot.ym(k); rc_plot.ym(k); rc_plot.ym(k)+temp_dy; rc_plot.ym(k)+temp_dy];

        end
    end
    
    
    %% Build the figure
    
    global rc_fig
    
    % Preallocate the structure with all figure data
    rc_fig = struct();
    % Define the strings/labels needed for the figure
    rc_fig.ylabel = [outp.INFO.Name{ismember(outp.INFO.Name, rc_var)} ' [' ...
        outp.INFO.Unit{ismember(outp.INFO.Name, rc_var)} ']'];
    rc_fig.xlabel = 'Time [hours CEST]';
    rc_fig.Xlabel = [check.int_xy_xz(1) ' [m]'];
    rc_fig.Ylabel = [check.int_xy_xz(2) ' [m]'];
    rc_fig.Clabel = [rc_var ' ['  ...
        outp.INFO.Unit{ismember(outp.INFO.Name, rc_var)} ']'];
    rc_fig.title = [outp.INFO.LongName{ismember(outp.INFO.Name, rc_var)} ' at [' ...
        num2str(outp.recp_indiv(1,1)) ', ' ...
        num2str(outp.recp_indiv(1,2)) ', ' ...
        num2str(outp.recp_indiv(1,3))  ']'];
    % Define the figure settings
    rc_fig.Xlim = rc_plot.cross_xlim;
    rc_fig.Ylim = rc_plot.cross_ylim;
    rc_fig.xlim = [round(min(rc_plot.time)*10 )/10, round(max(rc_plot.time))];
    rc_fig.ylim = [min(rc_plot.var(:)) max(rc_plot.var(:))];
        if rc_fig.ylim(1) < rc_fig.ylim(2) / 1000 && rc_fig.ylim(1) > -0.0001
            rc_fig.ylim(1) = 0;
        end
    rc_fig.Clim = rc_clim;
    rc_fig.INcross_time = 1;   
%     rc_fig.INcross_time = find( round(outp.time, 8) == round( , 8));   
    if strcmp(check.int_xy_xz, 'xy')
        rc_fig.receppos = [rc_plot.receppos(1), rc_plot.receppos(2)];
    elseif strcmp(check.int_xy_xz, 'xz')
        rc_fig.receppos = [rc_plot.receppos(1), rc_plot.receppos(3)];
    elseif strcmp(check.int_xy_xz, 'yz')
        rc_fig.receppos = [rc_plot.receppos(2), rc_plot.receppos(3)];
    end
    
    
    
    %
    %%%%% Build the figure with the subplots
    %
    rc_fig.fig = figure('units','pixels', 'Color','w',...
        'innerposition', rc_position, ...
        'Name', rc_var);
    rc_crosssubsize = [100, 80, rc_position(4) - 80-20, rc_position(4) - 80-20-25];
    rc_fig.crosssub = subplot('Position', ...
        rc_crosssubsize ./ [rc_position(3) rc_position(4) rc_position(3) rc_position(4)]);
    rc_timesub_position = [rc_crosssubsize(1)+rc_crosssubsize(3)+160, rc_crosssubsize(2)+150, ...
        0, 0];
    rc_timesub_position = rc_timesub_position + ...
        [0, 0, rc_position(3)-rc_timesub_position(1)-60, rc_position(4)-rc_timesub_position(2)-50];
    rc_fig.timesub = subplot('Position',  ...
         rc_timesub_position ./ [rc_position(3) rc_position(4) rc_position(3) rc_position(4)]);
    %
    %%%%% fill the crossection subplot
    % 
    f__cr_CrossSectionFill(rc_plot, outp, check)
    
    
    %
    %%%%% fill the time series subplot
    % 
    figure(rc_fig.fig)
    subplot(rc_fig.timesub)
    hold on 
    grid on
    % Plot the time series
    for j = 1 : length(check.int_avgt)

        plot(rc_plot.time(isnan(rc_plot.var(:,j))==0), ...
             rc_plot.var(isnan(rc_plot.var(:,j))==0, j), ...
            'Color',rc_colorset{j}, 'LineStyle', rc_linestyle{1}, ...
            'LineWidth',4)
    end
    % Get the Y limits of the figure
    rc_ax1 = gca;
    % Plot the time at which the cross section is shown
    rc_fig.timeline = plot([rc_plot.time(rc_fig.INcross_time), ...
                            rc_plot.time(rc_fig.INcross_time)], [...
                            rc_fig.ylim(1), rc_fig.ylim(2)], '-k', 'LineWidth', 1.5);
    
    % Show the range of the cross section
    rc_fig.timeline_scale = scatter(...
        [rc_plot.time(rc_fig.INcross_time), ...
         rc_plot.time(rc_fig.INcross_time)], ...
         [rc_fig.crossLimits(1), rc_fig.crossLimits(2)] , ...
            '+k', 'LineWidth', 1.5);
    
    % Settings
    if strcmp(check.int_logscale, 'yes')
        [rc_x_tick, rc_x_tick_label, rc_y_tick, rc_y_tick_label] = ...
            f__CLASSfig_XYTicks(rc_fig.xlim, [0 rc_fig.ylim(2)]);
    else
        [rc_x_tick, rc_x_tick_label, rc_y_tick, rc_y_tick_label] = ...
            f__CLASSfig_XYTicks(rc_fig.xlim, rc_fig.ylim);
    end   
    set(rc_ax1,'FontSize',18, ...
        'XLim',rc_fig.xlim, 'XTick',rc_x_tick, 'XTickLabel', rc_x_tick_label, ...
        'YLim',rc_fig.ylim,'YTick',rc_y_tick, 'YTickLabel',rc_y_tick_label)
    xlabel(rc_fig.xlabel);
    ylabel(rc_fig.ylabel);
    
    rc_lgnd_pos  = [0.40     0.13   0.17    0.055];
    rc_lgnd_dpos = [0.00    -0.085  0.00    0.005];
    lgnd_count = 0;
    rc_lgndtavg = struct();
    % Make averaging time legend
    for j = 1 : length(check.int_avgt)
        subplot('Position',[-j -j 0.50 0.78])
        set(gca, 'FontSize',18);
        rc_line = plot([0 1], [j-1 j], 'LineStyle', rc_linestyle{1}, ...
            'Color',rc_colorset{j}, 'LineWidth', 4);
        hold on
        
        % Variable legend
        rc_lgndtavg.(['line_' num2str(j)]) = legend(rc_line, {''});
        % Set the legend position
        rc_lgndtavg.(['line_' num2str(j)]).Position = ...
            rc_lgnd_pos + [0 0.0 0 0] + (lgnd_count) .* rc_lgnd_dpos;
        % Turn the legend box off
        rc_lgndtavg.(['line_' num2str(j)]).Box = 'off';
        % Add the variable name to the legend
        if check.int_avgt(j)*60 >= 1
            annotation('textbox', rc_lgnd_pos + [0.1 0.01 0 0] + (lgnd_count) .* rc_lgnd_dpos, ...
                    'string', ['t_{avg} = ' num2str(check.int_avgt(j)*60) ' min'], ...
                    'EdgeColor','none', 'FontSize',18, 'FontWeight','normal')
        else
            annotation('textbox', rc_lgnd_pos + [0.1 0.01 0 0] + (lgnd_count) .* rc_lgnd_dpos, ...
                    'string', ['t_{avg} = ' num2str(check.int_avgt(j)*60+30) ' s'], ...
                    'EdgeColor','none', 'FontSize',18, 'FontWeight','normal')
        end
        
        lgnd_count = lgnd_count+1;
        if ismember(j, [2,4])
            rc_lgnd_pos(1) = rc_lgnd_pos(1)+0.17;
            lgnd_count = 0;
        end

    end 
    % Add the time of the cross-section
    temp_hours = floor( rc_plot.time(rc_fig.INcross_time));
    temp_minutes = floor((rc_plot.time(rc_fig.INcross_time) - temp_hours)*60);
    temp_seconds = floor((rc_plot.time(rc_fig.INcross_time) - temp_hours - temp_minutes/60)*3600);
    temp_time = {num2str(temp_hours), num2str(temp_minutes), num2str(temp_seconds)};
    for i = 2: 3
        if length(temp_time{i}) == 1
            temp_time{i} = ['0' temp_time{i}];
        end
    end
    rc_timestamp = [temp_time{1} ':' temp_time{2} ':' temp_time{3}];
    clear temp_*
    rc_fig.rc_timestamp_pos_main = rc_timesub_position ./ [rc_position(3) rc_position(4) rc_position(3) rc_position(4)];
    rc_timestamp_pos(1) = (rc_fig.rc_timestamp_pos_main(1) + (rc_fig.rc_timestamp_pos_main(3)-0.15) * ...
        rc_fig.INcross_time / length(rc_plot.time)-0.021);
    rc_timestamp_pos(2) = rc_fig.rc_timestamp_pos_main(2)-0.18;
    rc_timestamp_pos(3) = 0.1;
    rc_timestamp_pos(4) = 0.05;
    rc_fig.timestamp = annotation('textbox', rc_timestamp_pos, ...
        'string', rc_timestamp, 'HorizontalAlignment','center', ...
        'EdgeColor','none', 'FontSize',14);
    % Figure title
    annotation('textbox', [100, rc_position(4)-50, rc_position(3)-200, 50] ./ ...
                           [rc_position(3), rc_position(4), rc_position(3), rc_position(4)], ...
        'string', rc_fig.title, 'HorizontalAlignment','center', ...
        'EdgeColor','none', 'FontSize',18, 'FontWeight','bold')
    
    % Build the time slider
    rc_sliderstep = [(rc_plot.time(2) - rc_plot.time(1)), 0.5]./ ...
            (rc_plot.time(end) - rc_plot.time(1));
    rc_fig.timeslide = uicontrol(rc_fig.fig,...
        'Style','slider', 'String','Cross section time', ...
        'Position',[rc_fig.timesub.Position(1)*rc_position(3)-10, ...
                    (rc_fig.timesub.Position(2)-0.22)*rc_position(4) ,...
                    rc_fig.timesub.Position(3)*rc_position(3)+52, ...
                    0.04*rc_position(4)], ...
        'Min',rc_plot.time(1), 'Max',rc_plot.time(end), 'SliderStep',rc_sliderstep, ...
        'Value', rc_plot.time(1), 'Callback',{@f__cr_TimeSlider, rc_plot, outp, check});
    
    
    
end

function f__cr_CrossSectionFill(rc_plot, outp, check)
    
    global rc_fig
    
    rc_colormap_Delta = 1000;
    rc_colormap = [[ [094 : (255 - 094)/(rc_colormap_Delta / 2) : 255 ]', ...
                     [060 : (255 - 060)/(rc_colormap_Delta / 2) : 255 ]', ...
                     [153 : (255 - 153)/(rc_colormap_Delta / 2) : 255 ]'] ...
                   ./255 ; ...
                   [ [255-(255 - 230)/(rc_colormap_Delta / 2) : -(255 - 230)/(rc_colormap_Delta / 2) : 230 ]', ...
                     [255-(255 - 097)/(rc_colormap_Delta / 2) : -(255 - 097)/(rc_colormap_Delta / 2) : 097 ]', ...
                     [255-(255 - 001)/(rc_colormap_Delta / 2) : -(255 - 001)/(rc_colormap_Delta / 2) : 001 ]'] ./255];

    
    figure(rc_fig.fig)
    subplot(rc_fig.crosssub)
    
    % Define indices fo the to-be-loaded data
    if  strcmp(check.int_xy_xz, 'xy')
        rc_INxC = find(outp.xt >= rc_plot.cross_xlim(1) & outp.xt <= rc_plot.cross_xlim(2) );
        rc_INyC = find(outp.yt >= rc_plot.cross_ylim(1) & outp.yt <= rc_plot.cross_ylim(2) );
    elseif  strcmp(check.int_xy_xz, 'xz')
        rc_INxC = find(outp.xt >= rc_plot.cross_xlim(1) & outp.xt <= rc_plot.cross_xlim(2) );
        rc_INyC = find(outp.zt >= rc_plot.cross_ylim(1) & outp.zt <= rc_plot.cross_ylim(2) );
    elseif  strcmp(check.int_xy_xz, 'yz')
        rc_INxC = find(outp.yt >= rc_plot.cross_xlim(1) & outp.yt <= rc_plot.cross_xlim(2) );
        rc_INyC = find(outp.zt >= rc_plot.cross_ylim(1) & outp.zt <= rc_plot.cross_ylim(2) );
    end
    rc_INtC = find(round(outp.time,8) == round( rc_plot.time(rc_fig.INcross_time), 8));
    % Load the data from the matfile
    rc_data = rc_plot.matfile.(check.int_var{1,2})(rc_INxC, rc_INyC, rc_INtC);
    % Put the data in the correct format
    rc_C = NaN(1, length(rc_plot.xm) * length(rc_plot.ym));
    for j = 1 : length(rc_plot.xm)
        for k = 1 : length(rc_plot.ym)
            rc_C(k + (j-1)*length(rc_plot.ym)) = rc_data(j,k);
        end
    end
    
    
    
    
    grid on         % show grid
    hold on         % show multiple lines
    rc_ax1 = gca;   % axis variable
    colormap(rc_colormap); 
    
    % If needed, remove the color bar in order to update it
    if ismember({'colorbar'}, fields(rc_fig))
        delete(rc_fig.colorbar)
    end
    
    
    % Plot the cross section plot
    rc_fig.crossplot = patch(rc_plot.X, rc_plot.Y, rc_C );
    % Show the location of the receptor(s)
    rc_fig.receptor = scatter(rc_fig.receppos(1,1), rc_fig.receppos(1,2), ...
        'xk', 'SizeData', 100, 'LineWidth',2);
    
    if strcmp(check.int_logscale, 'yes')
        set(rc_ax1, 'ColorScale','log', 'FontSize', 18)
    else
        set(rc_ax1, 'FontSize', 18)
    end   
    rc_fig.crossplot.LineWidth = 0.01;
    rc_fig.crossplot.EdgeAlpha = 0.2;
    rc_fig.colorbar = colorbar;
    if strcmp(check.int_Cticks_manual, 'yes')
        rc_TicksStr = cell(length(check.int_Cticks_num), 1);
        for i = 1 : length(check.int_Cticks_num)
            if ismember(round(check.int_Cticks_num(i),8), round(check.int_Cticks_str,8))
                rc_TicksStr{i} = num2str(check.int_Cticks_num(i));
            else
                rc_TicksStr{i} = '';
            end
        end
        rc_fig.colorbar.Ticks = check.int_Cticks_num;
        rc_fig.colorbar.TickLabels = rc_TicksStr;
        rc_fig.colorbar.Limits = [min(check.int_Cticks_num) max(check.int_Cticks_num)];
        rc_fig.crosssub.CLim = [min(check.int_Cticks_num) max(check.int_Cticks_num)];
        rc_fig.crossLimits = [min(check.int_Cticks_num) max(check.int_Cticks_num)];
    else
%         rc_fig.colorbar.Limits = rc_fig.Clim;
%         rc_fig.crosssub.CLim = rc_fig.Clim;
        rc_fig.crossLimits = rc_fig.crosssub.CLim;
    end
    rc_fig.colorbar.Label.String = '';%rc_fig.Clabel; 
    xlabel(rc_fig.Xlabel)
    ylabel(rc_fig.Ylabel)
    xlim(rc_fig.Xlim)
    ylim(rc_fig.Ylim)
    
    
end

function f__cr_TimeSlider(~,~, rc_plot, outp, check)
    
    global rc_fig
    
    
    rc_fig.INcross_time = min(find(rc_plot.time >= rc_fig.timeslide.Value)); 
    
    
    % Update the cross section figure
    figure(rc_fig.fig)
    cla(rc_fig.crosssub)
    f__cr_CrossSectionFill(rc_plot, outp, check)
    
    % Update the time series figure
    figure(rc_fig.fig)
    subplot(rc_fig.timesub)
    delete(rc_fig.timeline)         % Clear the time line
    delete(rc_fig.timeline_scale)   % Clear the cross section range
    delete(rc_fig.timestamp)   % Clear the cross section range
    % Get the Y limits of the figure
    rc_ax1 = gca;
    rc_fig.rc_ylim = rc_ax1.YLim;
    % Plot the time at which the cross section is shown
    rc_fig.timeline = plot([rc_plot.time(rc_fig.INcross_time), ...
                            rc_plot.time(rc_fig.INcross_time)], [...
                            rc_fig.rc_ylim(1), rc_fig.rc_ylim(2)], '-k', 'LineWidth', 1.5);
    % Show the range of the cross section
    rc_fig.timeline_scale = scatter(...
        [rc_plot.time(rc_fig.INcross_time), rc_plot.time(rc_fig.INcross_time)], ...
         [rc_fig.crossLimits(1), rc_fig.crossLimits(2)] , '+k', 'LineWidth', 1.5);
    
    
      % Add the time of the cross-section
    temp_hours = floor( rc_plot.time(rc_fig.INcross_time));
    temp_minutes = floor((rc_plot.time(rc_fig.INcross_time) - temp_hours)*60);
    temp_seconds = floor((rc_plot.time(rc_fig.INcross_time) - temp_hours - temp_minutes/60)*3600);
    temp_time = {num2str(temp_hours), num2str(temp_minutes), num2str(temp_seconds)};
    for i = 2: 3
        if length(temp_time{i}) == 1
            temp_time{i} = ['0' temp_time{i}];
        end
    end
    rc_timestamp = [temp_time{1} ':' temp_time{2} ':' temp_time{3}];
    clear temp_*
    rc_timestamp_pos(1) = (rc_fig.rc_timestamp_pos_main(1) + (rc_fig.rc_timestamp_pos_main(3)-0.055) * ...
        rc_fig.INcross_time / length(rc_plot.time))-0.016;
    rc_timestamp_pos(2) = rc_fig.rc_timestamp_pos_main(2)-0.18;
    rc_timestamp_pos(3) = 0.1;
    rc_timestamp_pos(4) = 0.05;
    rc_fig.timestamp = annotation('textbox', rc_timestamp_pos, ...
        'string', rc_timestamp, 'HorizontalAlignment','center', ...
        'EdgeColor','none', 'FontSize',14);
     
end



%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%   Visualize individual receptor statistics    %%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function f__FIGindividual_stats(outp, check, si_position, INrecp)
    
    si_colorset = {[230,097,001]./255, [094,060,153]./255,...
                   [253,184,099]./255, [178,171,210]./255};
    si_INinfo = find(ismember(outp.INFO.Name, check.ind_var));
    si_statnames = {'mean','std','fI','I','S','K'; ...
                'Mean value','Standard deviation','Fluctuation Intencity', 'Intermittency','Skewness','Kurtosis' ; ...
                '\langleC\rangle', '\sigma', '\sigma/\langleC\rangle', 'I', 'S', 'K'};
    
    % Data preparations
    si_matfile = matfile([outp.loc outp.recp_indiv_file{INrecp} check.ind_var '.mat']);
    si_INt = find(outp.time >= check.ind_twindow(1)-1 & outp.time <= check.ind_twindow(2));
    if strcmp(outp.recp_indiv_file{INrecp}(9:10), 'xy')
        si_INpos = [find( round(outp.xt,6) == round(outp.recp_indiv(INrecp,1),6) ), ...
            find( round(outp.yt,6) == round(outp.recp_indiv(INrecp,2),6) )];
    elseif strcmp(outp.recp_indiv_file{INrecp}(9:10), 'xz')
        si_INpos = [find( round(outp.xt,6) == round(outp.recp_indiv(INrecp,1),6) ), ...
            find( round(outp.zt,6) == round(outp.recp_indiv(INrecp,3),6) )];
    elseif strcmp(outp.recp_indiv_file{INrecp}(9:10), 'yz')
        si_INpos = [find( round(outp.yt,6) == round(outp.recp_indiv(INrecp,2),6) ), ...
            find( round(outp.zt,6) == round(outp.recp_indiv(INrecp,3),6) )];
    end
    
    % "in-figure" positions
    si_subpos =    [0.06 0.34 0.48 0.60];
    si_subbarpos = [0.58 0.23 0.41 0.72];
    si_lgndpos =   [0.42, 0.08 0.12 0.13];
    si_annodpos =  [0.05, 0.175 0.05 0.045];
    % Preparations of figure strings
    si_title = [check.ind_var ...
        ' at x = ' num2str(outp.recp_indiv(INrecp,1)) ' m, '...
        ' at y = ' num2str(outp.recp_indiv(INrecp,2)) ' m, '...
        ' at z = ' num2str(outp.recp_indiv(INrecp,3)) ' m'];
    si_ylabel = [outp.INFO.Name{si_INinfo} ' ' outp.INFO.Unit{si_INinfo}];
    si_xlabel = 'Time CEST [h]';
    % Data statistics
    si_t    = outp.time(si_INt);
    si_dt   = median(si_t(2:end) - si_t(1:end-1));
    si_C    = squeeze(si_matfile.(check.ind_var)(si_INpos(1), si_INpos(2), si_INt));
    si_time_new = si_t;
    if round(check.ind_avgt, 6) > round(si_t(5) - si_t(4), 6)
        si_C_old = si_C;
        si_time_new = [si_t(1) - si_dt+check.bs_avgt : check.bs_avgt : si_t(end)]';
        si_C = NaN( length(si_time_new), 1 );
        for i = 1 : length(si_time_new)
            if i == 1 
                si_C(i) = nanmean(si_C_old(round(si_t, 8) <= round(si_time_new(i), 8) ) );
            else
                si_C(i) = nanmean(si_C_old(round(si_t, 8) >  round(si_time_new(i-1), 8) & ...
                    round(si_t, 8) <=  round(si_time_new(i), 8) ) );
            end

        end
    end
    
    % Calculate moving average concentration
    si_averaging_time_h = 1;
    si_Cma = movmean(si_C, [si_averaging_time_h / (si_time_new(4) - si_time_new(3)), 0]);
           
    INt_correct = find( si_time_new > check.ind_twindow(1) & si_time_new <= check.ind_twindow(2) );
    si_C = si_C(INt_correct);
    si_Cma = si_Cma(INt_correct);
    
    si_stat = NaN(6,1);
    si_stat(1) = mean(si_C);
    si_stat(2) = std( si_C - si_Cma );
    si_stat(3) = si_stat(2) / si_stat(1);
    si_stat(5) = skewness( si_C - si_Cma );
    si_stat(6) = kurtosis( si_C - si_Cma );
    % Check if the NH3 variable is a "plume" or "background" variable or not
    if strcmp(check.ind_var(end - length(check.variable_plume_indicator)+1: end), check.variable_plume_indicator) || ...
            strcmp(check.ind_var(end - length(check.variable_background_indicator)+1: end), check.variable_background_indicator)
        
        si_stat(4) = sum( si_C >= check.ind_Ithreshhold ) ./ length(si_C);
    else
        
        % Define plume concentration from the combined NH3 concentration
        si_var_p = check.ind_var(4:end);
        si_var_p = ['nh3_' si_var_p check.variable_plume_indicator];
        si_matfile_p = matfile([outp.loc outp.recp_indiv_file{INrecp} si_var_p '.mat']);
        si_C_p = squeeze(si_matfile_p.(si_var_p)(si_INpos(1), si_INpos(2), si_INt));
        
        if round(check.ind_avgt, 6) > round(si_t(5) - si_t(4), 6)
            si_C_p_old = si_C_p;
            si_C_p = NaN( length(si_time_new) );
            for i = 1 : length(si_time_new)
                if i == 1 
                    si_C_p(i) = nanmean(si_C_p_old(round(si_t, 8) <= round(si_time_new(i), 8) ) );
                else
                    si_C_p(i) = nanmean(si_C_p_old(round(si_t, 8) >  round(si_time_new(i-1), 8) & ...
                        round(si_t, 8) <=  round(si_time_new(i), 8) ) );
                end

            end
        end
        si_C_p = si_C_p(INt_correct);
        
        % Calculate the intermittency
        si_stat(4) = sum( si_C_p >= check.ind_Ithreshhold ) ./ length(si_C_p);
        
        clear si_C_p si_matfile_p si_N_60min si_N_30min
         
    end
    si_t = si_time_new(INt_correct);

    % Ylimits
    if isnan( check.ind_Clim(1) )
        check.ind_Clim(1) = min(si_C);
    end
    if isnan( check.ind_Clim(2) )
        check.ind_Clim(2) = max(si_C);
    end
    si_Clim = check.ind_Clim;
    if si_Clim(2) < min(si_C)
        si_Clim = [min(si_C), max(si_C)];
    elseif si_Clim(1) > max(si_C)
        si_Clim = [min(si_C), max(si_C)];
    end
    if si_Clim(1) == si_Clim(2)
        si_Clim(1) = si_Clim(1)-1;
        si_Clim(2) = si_Clim(2)+1;
    end
    if (strcmp(check.ind_var(end - length(check.variable_plume_indicator)+1: end), check.variable_plume_indicator) || ...
            strcmp(check.ind_var(end - length(check.variable_background_indicator)+1: end), check.variable_background_indicator)) ...
            == 0
        
        si_Clim(1) = min(si_C- si_Cma);
    end
    
    
    % Figure preparation
    figure('units','pixels', 'Color','w',...
        'innerposition', si_position);
    if strcmp(check.ind_var(end - length(check.variable_plume_indicator)+1: end), check.variable_plume_indicator) || ...
            strcmp(check.ind_var(end - length(check.variable_background_indicator)+1: end), check.variable_background_indicator)
        
        si_lgnd1_line = gobjects(1,4);
    else
        si_lgnd1_line = gobjects(1,6);
    end
    %%%%%
    % Subplot of the time series  
    %%%%%
    subplot('Position', si_subpos)
    hold on; grid on; 
    % Plot statistics
    si_lgnd1_line(4) = plot([-999 999], ...
        [check.ind_Ithreshhold, check.ind_Ithreshhold], ...
        'LineStyle',':', 'Color','k', 'LineWidth',1.5);
    si_lgnd1_line(2) = plot([-999 999], [si_stat(1) si_stat(1)], ...
        'Color',si_colorset{2}, 'LineWidth',2.5);
    si_lgnd1_line(3) = plot([-999 999], [si_stat(1)-si_stat(2) si_stat(1)-si_stat(2)], ...
        'Color',si_colorset{4}, 'LineWidth',2.5);
    plot([-999 999], [si_stat(1)+si_stat(2) si_stat(1)+si_stat(2)], ...
        'Color',si_colorset{4}, 'LineWidth',2.5);
    % Plot receptor observations
    si_lgnd1_line(1) = plot(si_t, si_C, 'Color',si_colorset{1}, 'LineWidth',3);
    % Show additional data in case of a combine NH3 variable
    if strcmp(check.ind_var(end - length(check.variable_plume_indicator)+1: end), check.variable_plume_indicator) || ...
            strcmp(check.ind_var(end - length(check.variable_background_indicator)+1: end), check.variable_background_indicator)
       
        si_lgnd2 = legend(si_lgnd1_line, {'C','\langleC\rangle','\sigma','threshhold'});
    else
         
        si_lgnd1_line(5) = plot(si_t, si_Cma, ...
            'Color',si_colorset{2}, 'LineStyle','--', 'LineWidth',1.5);
        si_lgnd1_line(6) = plot(si_t, (si_C - si_Cma), ...
            'Color',si_colorset{2}, 'LineStyle','-.', 'LineWidth',1.5);
        
        si_lgnd2 = legend(si_lgnd1_line, {'C','\langleC\rangle','\sigma','I threshhold', 'C_{MA}','(C - C_{MA})'});
    end
    si_lgnd2.Position =   si_lgndpos;
    si_lgnd2.Box = 'off';
    % Figure settings
    if (si_Clim(2) - si_Clim(1)) > 10
        si_Clim(1) = round(si_Clim(1)-0.5);
        si_Clim(2) = round(si_Clim(2)+0.5);
    end
    [si_x_tick, si_x_tick_label, si_y_tick, si_y_tick_label] = ...
        f__CLASSfig_XYTicks([si_t(1) si_t(end)], si_Clim);    
    set(gca, 'FontSize',16, ...
        'XLim',[si_t(1), si_t(end)], 'XTick',si_x_tick, 'XTickLabel',si_x_tick_label, ...
        'YLim',si_Clim,              'YTick',si_y_tick, 'YTickLabel',si_y_tick_label)
    xlabel(si_xlabel)
    ylabel(si_ylabel)
    title(si_title)
    % Show the statistics as annotations
    jjj = 1; 
    for j =  1 : 6
        if j <= 3
            si_anno_dx = 0;
        else
            if si_anno_dx == 0
                jjj = 1;
            end
            si_anno_dx = 0.18;
        end
        annotation('textbox',si_annodpos + [si_anno_dx, -0.045*(jjj-1) 0 0], ...
            'string',si_statnames{3,j}, ...
            'EdgeColor','none', 'FontSize',16, 'FontWeight','normal')
        annotation('textbox',si_annodpos + [si_anno_dx+0.05, -0.045*(jjj-1) 0.025 0], ...
            'string',['= ' num2str(round( si_stat(j) ,3))], ...
            'EdgeColor','none', 'FontSize',16, 'FontWeight','normal')
        jjj = jjj+1;
    end   
    %%%%%
    % Subpot of the bar plot
    %%%%%
    si_sub3 = subplot('Position',si_subbarpos);  
    hold on; grid on
    si_ax1 = gca;
    histogram(si_C, 20, 'FaceColor',si_colorset{1}, 'EdgeColor',si_colorset{1})
    plot([si_stat(1) si_stat(1)], si_ax1.YLim, 'Color',si_colorset{2}, 'LineWidth',2.5)
    plot([si_stat(1)-si_stat(2) si_stat(1)-si_stat(2)], si_ax1.YLim, 'Color',si_colorset{4}, 'LineWidth',2.5)
    plot([si_stat(1)+si_stat(2) si_stat(1)+si_stat(2)], si_ax1.YLim, 'Color',si_colorset{4}, 'LineWidth',2.5)
    ylabel('#')
    xlabel(si_ylabel)
    set(gca, 'FontSize',16, 'XLim',[min(si_C)-0.01*si_stat(1), max(si_C)+0.01*si_stat(1)])


end
%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%           Receptor bulk statistics            %%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function nh3_blendthreshold = f__FIGbulk_stats(outp, check, INi, INj)
    
    nh3_colorset = {[230,097,001]./255, [094,060,153]./255,...
                   [253,184,099]./255, [178,171,210]./255};
    
    nh3_set = check.bs_var{INi,1};
    nh3_var = check.bs_var{INi,2};
    nh3_statname = check.bs_stat{INj,1};
    nh3_ylim = check.bs_stat{INj,2};
    INt = find( outp.(['set' nh3_set]).time > check.bs_twindow(1)-1 & ...
        outp.(['set' nh3_set]).time <= check.bs_twindow(2) );
    nh3_t = outp.(['set' nh3_set]).time(INt);
    INtplume = find( outp.(['set' nh3_set]).time > check.bs_twindow(1) & ...
        outp.(['set' nh3_set]).time <= check.bs_twindow(2) );
    nh3_tplume = outp.(['set' nh3_set]).time(INtplume);
    nh3_dt = median(nh3_t(2:end) - nh3_t(1:end-1));
    nh3_x = outp.(['set' nh3_set]).xt;
    if strcmp(check.bs_xy_xz, 'xz')
        nh3_y = outp.(['set' nh3_set]).zt;
        
        nh3_matfile_base = [outp.(['set' nh3_set]).loc nh3_set 'crossxz_'];
    else
        nh3_y = outp.(['set' nh3_set]).yt;
        nh3_zpos = outp.(['set' nh3_set]).zpos;
        [~,nh3_INz] = sort(abs(nh3_zpos - outp.(['set' nh3_set]).recp_bulk_out(1,3)));
        
        nh3_cross_lvl = outp.(['set' nh3_set]).cross_lvl{ nh3_INz(1) };
        nh3_matfile_base = [outp.(['set' nh3_set]).loc nh3_set 'crossxy' nh3_cross_lvl '_'];
    end
    
    % Define the x-y position of the source
    nh3_srcpos_x = [];
    for i = 1 : length(nh3_y)
        if sum(outp.(['set' nh3_set]).srcpos.(nh3_var)(:,i) == 1) > 0
            nh3_srcpos_x = [nh3_srcpos_x ; nh3_x(...
                outp.(['set' nh3_set]).srcpos.(nh3_var)(:,i) == 1)];
        end
    end
    nh3_srcpos_x = nh3_srcpos_x + ...
        0.5*median(nh3_x(2:end) - nh3_x(1:end-1));
    nh3_srcpos_x = max(nh3_srcpos_x);
    nh3_srcpos_y = [];
    for i = 1 : length(nh3_x)
        if sum(outp.(['set' nh3_set]).srcpos.(nh3_var)(i,:) == 1) > 0
            nh3_srcpos_y = [nh3_srcpos_y ; nh3_y(...
                outp.(['set' nh3_set]).srcpos.(nh3_var)(i,:) == 1)];
        end
    end
    nh3_srcpos_y = median(nh3_srcpos_y);
    nh3_srcpos = [nh3_srcpos_x, nh3_srcpos_y];
    

    % Make sure the variable is an nh3 variable
    nh3_setname = nh3_var;
    if strcmp(nh3_setname(1:3), 'nh3') == 0
        nh3_setname = outp.((['set' nh3_set])).INFO.Name{13};
    end
    % Find the name of the scalar set
    if strcmp(nh3_setname(end - length(check.variable_plume_indicator)+1: end), ...
            check.variable_plume_indicator)

        nh3_setname(end - length(check.variable_plume_indicator)+1: end) = '';
    elseif strcmp(nh3_setname(end - length(check.variable_background_indicator)+1: end), ...
            check.variable_background_indicator)

        nh3_setname(end - length(check.variable_background_indicator)+1: end) = '';
    end
    if strcmp(nh3_setname(4), '_')
        nh3_setname = nh3_setname(5:end);
    else
        nh3_setname = nh3_setname(4:end);
    end
    % Make an exception for nhr_r0b, since set r0 does not have a plume
    if strcmp(nh3_var, 'nh3_r0b')
        nh3_varplume = 'nh3_r1p';
    end
    % Set the name of the in-plume and background variables of the set
    nh3_varplume = ['nh3_' nh3_setname check.variable_plume_indicator];
    nh3_varbackground = ['nh3_' nh3_setname check.variable_background_indicator];
    % Define the matfiles
    nh3_matfile = matfile([nh3_matfile_base nh3_var '.mat']);
    nh3_matfile_plume = matfile([nh3_matfile_base nh3_varplume '.mat']);
    nh3_matfile_background = matfile([nh3_matfile_base nh3_varbackground '.mat']);
    nh3_matfile_w = matfile([nh3_matfile_base 'w.mat']);
    
    if strcmp(nh3_statname, 'I')
        nh3_var = nh3_varplume;
        nh3_matfile = nh3_matfile_plume;
    end

    %% Load the relevant data
    
    %%%%%
    % Load the in-plume data
    %%%%%
    nh3_time_new = nh3_t;
    
    fprintf(['      - Loading ' check.bs_var{1} '-' nh3_var ' data: \n'])
    if round(check.bs_avgt, 8) > round(nh3_dt, 8)
%         if strcmp(nh3_set,'012') && round(check.bs_avgt, 8) == round(10/3600,8)
%             nh3_time_new = [nh3_t(1) - nh3_dt+check.bs_avgt : check.bs_avgt : nh3_t(end)];
%             nh3_in_C = NaN(length(nh3_x), length(nh3_y), length(nh3_time_new));
%             INt = find(ismember(round(nh3_t,8), round(nh3_time_new,8)));
%             
%             nh3_percentage = 0.20;
%             nh3_percentage_base = nh3_percentage;
%             for i = 1 : length(INt)
%                 nh3_in_C(:,:,i) = nh3_matfile.(nh3_var)(:,:, INt(i));
%                 
%                 if i/length(nh3_time_new) >= nh3_percentage_base
% 
%                     fprintf(['         ' check.bs_var{1} '-' nh3_var ' data at '...
%                         num2str(nh3_percentage_base * 100) ' %%\n'])
%                     nh3_percentage_base = nh3_percentage_base + nh3_percentage;
%                 end
%             end
%             
%         else
            % Define new time variable based on the new averaging time
            nh3_time_new = [nh3_t(1) - nh3_dt+check.bs_avgt : check.bs_avgt : nh3_t(end)];
            nh3_in_C = NaN(length(nh3_x), length(nh3_y), length(INt)/(check.bs_avgt*3600));

            nh3_percentage = 0.20;
            nh3_percentage_base = nh3_percentage;

            for i = 1 : length(nh3_time_new)
                if i == 1 
                    nh3_in_C(:,:,i) = mean(nh3_matfile.(nh3_var)(:,:, ...
                        INt(nh3_t <= nh3_time_new(i))), 3);
                else
                    nh3_in_C(:,:,i) = mean(nh3_matfile.(nh3_var)(:,:, ...
                        INt(nh3_t > nh3_time_new(i-1) & nh3_t <= nh3_time_new(i))), 3);
                end

                if i/length(nh3_time_new) >= nh3_percentage_base

                    fprintf(['         ' check.bs_var{1} '-' nh3_var ' data at '...
                        num2str(nh3_percentage_base * 100) ' %%\n'])
                    nh3_percentage_base = nh3_percentage_base + nh3_percentage;
                end
            end
%         end
    elseif round(check.bs_avgt, 8) == round(nh3_dt, 8)
        
        nh3_in_C = nh3_matfile.(nh3_var)(:,:, INt);
    end
    % Calculate moving average concentration, if needed
    if ismember({nh3_statname}, {'std', 'fI', 'S', 'K'})
        nh3_averaging_time_h = 1;
        nh3_in_Cma = movmean(nh3_in_C, [nh3_averaging_time_h / (nh3_time_new(4) - nh3_time_new(3)), 0],3);
    end
    % Only keep the time data within the analysis window
    INt_correct = find( nh3_time_new > check.bs_twindow(1) & nh3_time_new <= check.bs_twindow(2) );
    nh3_in_C = nh3_in_C(:,:,INt_correct);
    if ismember({nh3_statname}, {'std', 'fI', 'S', 'K'})
        nh3_in_Cma = nh3_in_Cma(:,:,INt_correct);
    end
    
    %%%%%
    % Load the background data
    %%%%%
    if strcmp(nh3_statname, 'I') == 0
        fprintf(['      - Loading ' check.bs_var{1} '-' nh3_varbackground ' data: \n'])
        nh3_time_new = nh3_t;
        if round(check.bs_avgt, 8) > round(nh3_dt, 8)
%             if strcmp(nh3_set,'012') && round(check.bs_avgt, 8) == round(10/3600,8)
%                 nh3_time_new = [nh3_t(1) - nh3_dt+check.bs_avgt : check.bs_avgt : nh3_t(end)];
%                 nh3_out_C = NaN(length(nh3_x), length(nh3_y), length(nh3_time_new));
%                 INt = find(ismember(round(nh3_t,8), round(nh3_time_new,8)));
% 
%                 nh3_percentage = 0.20;
%                 nh3_percentage_base = nh3_percentage;
%                 for i = 1 : length(INt)
%                     nh3_out_C(:,:,i) = nh3_matfile_background.(nh3_varbackground)(:,:, INt(i));
% 
%                     if i/length(nh3_time_new) >= nh3_percentage_base
% 
%                         fprintf(['         ' check.bs_var{1} '-' nh3_varbackground ' data at '...
%                             num2str(nh3_percentage_base * 100) ' %%\n'])
%                         nh3_percentage_base = nh3_percentage_base + nh3_percentage;
%                     end
%                 end
% 
%             else
                % Define new time variable based on the new averaging time
                nh3_time_new = [nh3_t(1) - nh3_dt+check.bs_avgt : check.bs_avgt : nh3_t(end)];
                nh3_out_C = NaN(length(nh3_x), length(nh3_y), length(INt)/(check.bs_avgt*3600));

                nh3_percentage = 0.20;
                nh3_percentage_base = nh3_percentage;

                for i = 1 : length(nh3_time_new)
                    if i == 1 
                        nh3_out_C(:,:,i) = mean(nh3_matfile_background.(nh3_varbackground)(:,:, ...
                            INt(nh3_t <= nh3_time_new(i))), 3);
                    else
                        nh3_out_C(:,:,i) = mean(nh3_matfile_background.(nh3_varbackground)(:,:, ...
                            INt(nh3_t > nh3_time_new(i-1) & nh3_t <= nh3_time_new(i))), 3);
                    end

                    if i/length(nh3_time_new) >= nh3_percentage_base
                        fprintf(['         ' check.bs_var{1} '-' nh3_varbackground ' data at '...
                            num2str(nh3_percentage_base * 100) ' %%\n'])
                        nh3_percentage_base = nh3_percentage_base + nh3_percentage;
                    end
                end
%             end
        elseif round(check.bs_avgt, 8) == round(nh3_dt, 8)
            
            nh3_out_C = nh3_matfile_background.(nh3_varbackground)(:,:, INt);
        end
        % Calculate moving average concentration, if needed
        if ismember({nh3_statname}, {'std', 'fI', 'S', 'K'})
            nh3_out_Cma = movmean(nh3_out_C, [nh3_averaging_time_h / (nh3_time_new(4) - nh3_time_new(3)), 0],3);
            clear nh3_averaging_time_h
        end
        % Now apply the correct time window (previous one was needed for moving
        % average calculation)
        nh3_out_C = nh3_out_C(:,:,INt_correct);
        if ismember({nh3_statname}, {'std', 'fI', 'S', 'K'})
            nh3_out_Cma = nh3_out_Cma(:,:,INt_correct);
        end
    end
    
    %%%%%
    % Load the vertical wind speed data
    %%%%%
    if strcmp(nh3_statname, 'F')
        
        fprintf(['      - Loading ' check.bs_var{1} '-w data: \n'])
        
        nh3_plume_time_new = nh3_tplume;
        if round(check.bs_avgt, 8) > round(nh3_dt, 8)
%             if strcmp(nh3_set,'012') && round(check.bs_avgt, 8) == round(10/3600,8)
%                 nh3_time_new = [nh3_t(1) - nh3_dt+check.bs_avgt : check.bs_avgt : nh3_t(end)];
%                 nh3_w = NaN(length(nh3_x), length(nh3_y), length(nh3_time_new));
%                 INt = find(ismember(round(nh3_t,8), round(nh3_time_new,8)));
% 
%                 nh3_percentage = 0.20;
%                 nh3_percentage_base = nh3_percentage;
%                 for i = 1 : length(INt)
%                     nh3_w(:,:,i) = nh3_matfile_w.w(:,:, INt(i));
% 
%                     if i/length(nh3_time_new) >= nh3_percentage_base
% 
%                         fprintf(['         ' check.bs_var{1} '-w data at '...
%                             num2str(nh3_percentage_base * 100) ' %%\n'])
%                         nh3_percentage_base = nh3_percentage_base + nh3_percentage;
%                     end
%                 end
% 
%             else
                % Define new time variable based on the new averaging time
                nh3_plume_time_new = [nh3_tplume(1) - nh3_dt+check.bs_avgt : check.bs_avgt : nh3_tplume(end)];
                nh3_w = NaN(length(nh3_x), length(nh3_y), length(nh3_tplume)/(check.bs_avgt*3600));

                nh3_percentage = 0.20;
                nh3_percentage_base = nh3_percentage;

                for i = 1 : length(nh3_plume_time_new)
                    if i == 1 
                        nh3_w(:,:,i) = mean(nh3_matfile_w.w(:,:, ...
                            INtplume(nh3_tplume <= nh3_plume_time_new(i))), 3);
                    else
                        nh3_w(:,:,i) = mean(nh3_matfile_w.w(:,:, ...
                            INtplume(nh3_tplume > nh3_plume_time_new(i-1) & nh3_tplume <= nh3_plume_time_new(i))), 3);
                    end

                    if i/length(nh3_plume_time_new) >= nh3_percentage_base
                    fprintf(['         ' check.bs_var{1} '-w data at '...
                        num2str(nh3_percentage_base * 100) ' %%\n'])
                        nh3_percentage_base = nh3_percentage_base + nh3_percentage;
                    end
                end
%             end
        elseif round(check.bs_avgt, 8) == round(nh3_dt, 8)

            nh3_w = nh3_matfile_w.w(:,:, INtplume);
        end
    end
    
    clear nh3_matfile* nh3_tplume i INt_correct
    
    %% Calculate statistics
    
    %%%%%
    % In-plume statistics
    %%%%% 
    fprintf(['      - Calculate in-plume statistics\n'])
    
    % Calculate the statistics
    if strcmp(nh3_statname, 'mean')
        nh3_in_stat =  mean( nh3_in_C ,3);
        nh3_out_stat =  mean( nh3_out_C ,3);
        nh3_ylabel = '[C] [ppb]';
        nh3_title = 'Mean';
        nh3_norm_stat = mean( (nh3_in_stat - nh3_out_stat) ./ nh3_out_stat ,3);
        
    elseif strcmp(nh3_statname, 'std')
        nh3_in_stat =  std( nh3_in_C-nh3_in_Cma , 0, 3);
        nh3_out_stat =  std( nh3_out_C-nh3_out_Cma , 0, 3);
        nh3_title = 'Standard deviation';
        nh3_ylabel = '\sigma [ppb]';
        nh3_norm_stat = abs( (nh3_in_stat - nh3_out_stat) ./ nh3_out_stat );
        
    elseif strcmp(nh3_statname, 'fI')
        nh3_in_stat = std( nh3_in_C-nh3_in_Cma ,0,3) ./ mean( nh3_in_C ,3);
        % NOTE: the background fluctuation intensity for a plume variable
        % can be extremely large, since the concentrations approach zero.
        % The background fI for a plume variable will therefore be set to
        % zero
        if strcmp(nh3_var(end - length(check.variable_plume_indicator)+1: end), check.variable_plume_indicator)
            nh3_out_stat = zeros(size(nh3_in_stat));
        else
            nh3_out_stat = std( nh3_out_C-nh3_out_Cma ,0,3) ./ mean( nh3_out_C ,3);
        end
        nh3_norm_stat = abs( (nh3_in_stat - nh3_out_stat) ./ nh3_out_stat );
        nh3_title = 'Fluctuation intensity';
        nh3_ylabel = '\sigma/[C] [-]';
        
    elseif strcmp(nh3_statname, 'I')
        
        % NOTE: Intermittency can only be calculated over plume scalars
        nh3_in_stat = sum( nh3_in_C >= check.bs_Ithreshhold ,3) ./ size(nh3_in_C,3);
        % The intermittency out of the plume is always zero
        nh3_out_stat =  zeros(size(nh3_in_stat));
        nh3_norm_stat = nh3_in_stat;
        nh3_ylabel = 'I [-]';
        nh3_title = 'Intermittency'; 
            
    elseif strcmp(nh3_statname, 'S')
        nh3_in_stat =  skewness( nh3_in_C-nh3_in_Cma ,1,3);
        nh3_out_stat =  skewness( nh3_out_C-nh3_out_Cma ,1,3);
        nh3_norm_stat = abs( (nh3_in_stat - nh3_out_stat) ./ nh3_out_stat );
        nh3_title = 'Skewness';
        nh3_ylabel = 'S [-]';
        
    elseif strcmp(nh3_statname, 'K')
        nh3_in_stat =  kurtosis( nh3_in_C-nh3_in_Cma ,1,3);
        nh3_out_stat =  kurtosis( nh3_out_C-nh3_out_Cma ,1,3);
        nh3_norm_stat = abs( (nh3_in_stat - nh3_out_stat) ./ nh3_out_stat );
        nh3_title = 'Kurtosis';
        nh3_ylabel = 'K [-]';
        
    elseif strcmp(nh3_statname, 'F')
        
        % Calculate the flux each 30 minutes
        nh3_fluxtime = [check.bs_twindow(1) : 0.5 : check.bs_twindow(2)];
        nh3_flux = NaN(size(nh3_in_C,1), size(nh3_in_C,2), length(nh3_fluxtime)-1);
        for i = 1 : length(nh3_fluxtime)-1
            if i == 1
                nh3_fluxIN = find(nh3_plume_time_new <= nh3_fluxtime(i+1));
            else
                nh3_fluxIN = find(nh3_plume_time_new > nh3_fluxtime(i) & ...
                    nh3_plume_time_new <= nh3_fluxtime(i+1));
            end
            nh3_flux(:,:,i) = nanmean(( nh3_w(:,:,nh3_fluxIN) - nanmean(nh3_w(:,:,nh3_fluxIN), 3) ) ...
                .* (nh3_in_C(:,:,nh3_fluxIN) - nanmean(nh3_in_C(:,:,nh3_fluxIN), 3)), 3);

        end
        % Calculate the averaged flux within the time window
        nh3_in_stat =  mean(nh3_flux, 3);
        clear nh3_fluxtime nh3_flux nh3_fluxIN i
        
        % Calculate the flux each 30 minutes
        nh3_fluxtime = [check.bs_twindow(1) : 0.5 : check.bs_twindow(2)];
        nh3_flux = NaN(size(nh3_out_C,1), size(nh3_out_C,2), length(nh3_fluxtime)-1);
        for i = 1 : length(nh3_fluxtime)-1
            if i == 1
                nh3_fluxIN = find(nh3_plume_time_new <= nh3_fluxtime(i+1));
            else
                nh3_fluxIN = find(nh3_plume_time_new > nh3_fluxtime(i) & ...
                    nh3_plume_time_new <= nh3_fluxtime(i+1));
            end
            nh3_flux(:,:,i) = nanmean(( nh3_w(:,:,nh3_fluxIN) - nanmean(nh3_w(:,:,nh3_fluxIN), 3) ) ...
                .* (nh3_out_C(:,:,nh3_fluxIN) - nanmean(nh3_out_C(:,:,nh3_fluxIN), 3)), 3);

        end
        % Calculate the averaged flux within the time window
        nh3_out_stat =  mean(nh3_flux, 3);
        nh3_norm_stat = abs( (nh3_in_stat - nh3_out_stat) ./ nh3_out_stat );
        
        clear nh3_fluxtime nh3_flux nh3_fluxIN i
        
        nh3_ylabel = 'F [ppb m/s]';
        nh3_title = 'Flux';
        
    end
    
    nh3_out_stat_mean = mean(nh3_out_stat(:));
    nh3_out_stat_std = std(nh3_out_stat(:));
    
    clear nh3_C nh3_Cma nh3_counter nh3_temp i j INx INy nh3_w nh3_w_matfile nh3_w_temp
    
    %% Prepare data for the figure
    
    fprintf(['      - Calculate blending-distances\n'])
    
    % Determine the threshold value for the blending-distance
    nh3_blendthreshold = NaN(length(check.bs_dist_level),1);
    for i = 1 : length(nh3_blendthreshold)
        if strcmp(nh3_statname, 'I') == 0

            if strcmp(check.bd_peak_abs, 'peak')
                nh3_blendthreshold(i) = check.bs_dist_level(i) * max(nh3_norm_stat(:));
            elseif strcmp(check.bd_peak_abs, 'abs')
                nh3_blendthreshold(i) = check.bs_dist_level(i);
            else
                nh3_blendthreshold(i) = 0;
            end
        else
            nh3_blendthreshold(i) = 0;
        end
    end
    
    % Calculate absolute distance
    nh3_pos = NaN(length(nh3_x), length(nh3_y));
    nh3_pos_x = NaN(length(nh3_x), length(nh3_y));
    nh3_pos_y = NaN(length(nh3_x), length(nh3_y));
    for i = 1 : length(nh3_x)
        nh3_pos(i,:) = sqrt( (nh3_x(i) - nh3_srcpos_x).^2 + (nh3_y - nh3_srcpos_y).^2 );
        nh3_pos_x(i,:) = nh3_x(i);
    end
    for i = 1 : length(nh3_y)
        nh3_pos_y(:,i) = nh3_y(i);
    end
    
    % Find the maximum values of the normalized statistics for absolute distance from the source          
    nh3_statmax = nh3_norm_stat(:);
    nh3_pos_statmax = nh3_pos(:);
    nh3_pos_statmax(isnan(nh3_statmax)) = [];
    nh3_statmax(isnan(nh3_statmax)) = [];

    temp_window = 50;       % Window of XX m within the maximum value is determined
    temp_maxplot = NaN(size(nh3_statmax));
    temp_maxdist = NaN(size(nh3_statmax));
    for i = 1 : length(nh3_statmax)
        INmaxdist = find( nh3_pos_statmax >= nh3_pos_statmax(i)-temp_window/2 & ...
            nh3_pos_statmax <= nh3_pos_statmax(i)+temp_window/2 );

        [temp_maxplot(i), temp_IN] = max( nh3_statmax(INmaxdist) );
        temp_maxdist(i) = nh3_pos_statmax(INmaxdist(temp_IN));
    end
    nh3_maxplot = unique(temp_maxplot);
    nh3_maxdist = NaN(size(nh3_maxplot));
    for i = 1 : length(nh3_maxplot)
        nh3_maxdist(i) = max(temp_maxdist(temp_maxplot == nh3_maxplot(i)));
    end
    [nh3_maxdist, INsort] = sort(nh3_maxdist);
    nh3_maxplot = nh3_maxplot(INsort);
    nh3_norm_maxstat = [nh3_maxdist, nh3_maxplot];
    if strcmp(nh3_statname, 'I')
        if nh3_norm_maxstat(end,2) == 0
            nh3_norm_maxstat(end,1) = nh3_norm_maxstat(end-1,1);
        end
    end
    clear nh3_statmax nh3_pos_statmax temp_* INmaxdist INsort
    
    
    % Find the values at the plume centerline
    nh3_in_centreline = [ abs(nh3_x - nh3_srcpos_x) , ...
        nh3_in_stat(:, round(nh3_y,6) == round(nh3_srcpos(2),6))];
    % Find the normalized values at the plume centerline
    nh3_norm_centreline = [ abs(nh3_x - nh3_srcpos_x) , ...
        nh3_norm_stat(:, round(nh3_y,6) == round(nh3_srcpos(2),6))];

    %%%%%
    % Calculate blending distance for absolute statistics, centreline
    %%%%%
    % Find the local maximum of the in-plume statistics within the first 500 m
    nh3_INmax = find( nh3_in_centreline(:,2) == max(nh3_in_centreline(nh3_in_centreline(:,1) <= 500, 2) ) );
    % Determine 'blending-distance'
%     nh3_INsafedist = find(nh3_in_statline <= check.bs_dist_level .* mean(nh3_out_statline));
    nh3_INsafedist = find(nh3_in_centreline(:,2) <= nh3_out_stat_mean + 3*nh3_out_stat_std+0.000001);
    % Remove values which are found before the maximum is reached.
    nh3_INsafedist(nh3_INsafedist < nh3_INmax) = [];
    if isempty(nh3_INsafedist)
        nh3_blend_centre_abs = NaN;
    else
        nh3_blend_centre_abs = nh3_in_centreline(nh3_INsafedist(1), 1);
    end

    %%%%%
    % Calculate blending distance for absolute statistics, Maximum
    %%%%%
    % Find the local maximum of the in-plume statistics within the first 500 m
    nh3_INmax = find( nh3_in_centreline(:,2) == max(nh3_in_centreline(nh3_in_centreline(:,1) <= 500, 2) ) );
    % Determine 'blending-distance'
%     nh3_INsafedist = find(nh3_in_statline <= check.bs_dist_level .* mean(nh3_out_statline));
    nh3_INsafedist = find(nh3_in_centreline(:,2) <= nh3_out_stat_mean + 3*nh3_out_stat_std+0.000001);
    if isempty(nh3_INsafedist)
        nh3_blend_centre_abs = NaN;
    else
        % Remove values which are found before the maximum is reached.
        nh3_INsafedist(nh3_INsafedist < nh3_INmax) = [];
        nh3_blend_centre_abs = nh3_in_centreline(nh3_INsafedist(1), 1);
    end

    nh3_blend_centre = NaN(length(check.bs_dist_level),1);
    nh3_blend_maxstat = NaN(length(check.bs_dist_level),1);
    for i = 1 : length(nh3_blendthreshold)
        %%%%%
        % Calculate blending distance for normalized statistics, centreline
        %%%%%
        if max(nh3_norm_centreline(:,2)) < nh3_blendthreshold(i)
            nh3_blend_centre(i) = 0;
        else
            % Find the local maximum of the in-plume statistics within the first 500 m
            nh3_INmax = find( nh3_norm_centreline(:,2) == max(nh3_norm_centreline(nh3_norm_centreline(:,1) <= 500, 2) ) );
            % Determine 'blending-distance'
            nh3_INsafedist = find(nh3_norm_centreline(:,2) <= nh3_blendthreshold(i));
            if isempty(nh3_INsafedist) == 0
                % Remove values which are found before the maximum is reached.
                nh3_INsafedist(nh3_INsafedist < nh3_INmax) = [];
                nh3_blend_centre(i) = nh3_norm_centreline(nh3_INsafedist(1), 1);
            end
        end

        %%%%%
        % Calculate blending distance for normalized statistics, Maximum
        %%%%%
        if max(nh3_norm_maxstat(:,2)) < nh3_blendthreshold(i)
            nh3_blend_maxstat(i) = 0;
        else
            % Find the local maximum of the in-plume statistics within the first 500 m
            nh3_INmax = find( nh3_norm_maxstat(:,2) == max(nh3_norm_maxstat(nh3_norm_maxstat(:,1) <= 500, 2) ) );
            % Determine 'blending-distance'
            nh3_INsafedist = find(nh3_norm_maxstat(:,2) <= nh3_blendthreshold(i));
            if isempty(nh3_INsafedist) == 0
                % Remove values which are found before the maximum is reached.
                nh3_INsafedist(nh3_INsafedist < nh3_INmax) = [];
                nh3_blend_maxstat(i) = nh3_norm_maxstat(nh3_INsafedist(1), 1);
            end
        end
    end
    
    % Relative statistics y-labels 
    if strcmp(nh3_statname, 'mean')
        nh3_rel_ylabel = 'PC_{NH_3} [%]';
    elseif strcmp(nh3_statname, 'std')
        nh3_rel_ylabel = 'PC_{\sigma} [%]';
    elseif strcmp(nh3_statname, 'I')
        nh3_rel_ylabel = 'I [-]';
    else
        nh3_rel_ylabel = ['PC_{' nh3_statname '} [%]'];
    end
    
    
    %% Make an interesting figure
    
%     fprintf(['      - Build 3 interesting figures\n'])
%     
%     % Absolute distance
%     figure
%     hold on; grid on
%     for i = 1 : size(nh3_in_stat, 2)
%         if i == 1
%             nh3_quickfig(1) = scatter(nh3_pos(:,i), nh3_out_stat(:,i), 'ok', 'SizeData', 30);
%         else
%             scatter(nh3_pos(:,i), nh3_out_stat(:,i), 'ok', 'SizeData', 30)
%         end
%     end
%     for i = 1 : size(nh3_in_stat, 2)
%         if i == 1
%             nh3_quickfig(2) = scatter(nh3_pos(:,i), nh3_in_stat(:,i), '.b', 'SizeData', 60);
%         else
%             scatter(nh3_pos(:,i), nh3_in_stat(:,i), '.b', 'SizeData', 60)
%         end
%     end
%     nh3_quickfig(3) = plot(nh3_in_centreline(:,1), nh3_in_centreline(:,2), '--g', 'LineWidth', 2.5);
% %     scatter(nh3_norm_maxstat(:,1), nh3_norm_maxstat(:,2), '.b', 'SizeData', 60);
%     legend(nh3_quickfig, {'background','in-plume','centreline'})
%     set(gca, 'FontSize', 16)
%     xlabel('|Distance| [m]')
%     ylabel(nh3_ylabel)
%     clear nh3_quickfig
%     
%     
%     % Distance in x-direciton
%     figure
%     hold on; grid on
%     for i = 1 : size(nh3_in_stat, 2)
%         if i == 1
%             nh3_quickfig(1) = scatter(nh3_pos_x(:,i) - nh3_srcpos_x, nh3_out_stat(:,i), 'ok', 'SizeData', 30);
%         else
%             scatter(nh3_pos_x(:,i) - nh3_srcpos_x, nh3_out_stat(:,i), 'ok', 'SizeData', 30)
%         end
%     end
%     for i = 1 : size(nh3_in_stat, 2)
%         if i == 1
%             nh3_quickfig(2) = scatter(nh3_pos_x(:,i) - nh3_srcpos_x, nh3_in_stat(:,i), '.b', 'SizeData', 60);
%         else
%             scatter(nh3_pos_x(:,i) - nh3_srcpos_x, nh3_in_stat(:,i), '.b', 'SizeData', 60)
%         end
%     end
%     legend(nh3_quickfig, {'background','in-plume'})
%     set(gca, 'FontSize', 16)
%     xlabel('x-distance [m]')
%     xlim([min(nh3_pos_x(:) - nh3_srcpos_x) max(nh3_pos_x(:) - nh3_srcpos_x)])
%     ylabel(nh3_ylabel)
%     clear nh3_quickfig
%     
%     
%     % Distance in y-direciton
%     figure
%     hold on; grid on
%     for i = 1 : size(nh3_in_stat, 2)
%         if i == 1
%             nh3_quickfig(1) = scatter(nh3_pos_y(:,i) - nh3_srcpos_y, nh3_out_stat(:,i), 'ok', 'SizeData', 30);
%         else
%             scatter(nh3_pos_y(:,i) - nh3_srcpos_y, nh3_out_stat(:,i), 'ok', 'SizeData', 30)
%         end
%     end
%     for i = 1 : size(nh3_in_stat, 2)
%         if i == 1
%             nh3_quickfig(2) = scatter(nh3_pos_y(:,i) - nh3_srcpos_y, nh3_in_stat(:,i), '.b', 'SizeData', 60);
%         else
%             scatter(nh3_pos_y(:,i) - nh3_srcpos_y, nh3_in_stat(:,i), '.b', 'SizeData', 60)
%         end
%     end
%     legend(nh3_quickfig, {'background','in-plume'})
%     set(gca, 'FontSize', 16)
%     xlabel('y-distance [m]')
%     xlim([min(nh3_pos_y(:) - nh3_srcpos_y) max(nh3_pos_y(:) - nh3_srcpos_y)])
%     ylabel(nh3_ylabel)
%     clear nh3_quickfig
%     
%     
%     clear nh3_C nh3_Cma nh3_counter nh3_temp i j INx INy  nh3_w nh3_w_matfile nh3_w_temp
    
    
    %% Build the x-distance absolute stats figure
    
    if strcmp(check.bs_show_xdist_fig, 'yes')
        fprintf(['      - Build the X-distrance Vs. statistics figure\n'])

        nh3_xdist = abs(nh3_pos_x - nh3_srcpos_x);

        figure('units','pixels', 'Color','w',...
            'innerposition', [10 200 1300 470], ...
            'Name', 'Receptor stats');
        subplot('Position',[0.10 0.16 0.67 0.71]);
        nh3_obj = gobjects(3,1);

        grid on      % show grid
        hold on         % show multiple lines
        nh3_ax1 = gca;   % axis variable
        % For the flux, show the zero line
        if strcmp(nh3_statname, 'F')
            if isnan(check.bs_xlim(1)) || isnan(check.bs_xlim(2))
                plot([min(nh3_x - nh3_srcpos(1)) max(nh3_x - nh3_srcpos(1))], [0 0], '-k', 'LineWidth', 1.5)
            else
                plot(check.bs_xlim, [0 0], '-k', 'LineWidth', 1.5)
            end
        end
        % In-plume Receptor statistics
        for i = 1 : size(nh3_in_stat,2)
            if i == 1
                nh3_obj(5) = scatter(nh3_pos_x(:,i), nh3_in_stat(:,i), ...
                        'Marker','o', 'LineWidth',3, 'SizeData',20, ...
                        'MarkerFaceColor',nh3_colorset{3}, 'MarkerFaceAlpha',0.99, ...
                        'MarkerEdgeColor','none', 'MarkerEdgeAlpha',0.99 );
            else
                scatter(nh3_xdist(:,i), nh3_in_stat(:,i), ...
                        'Marker','o', 'LineWidth',3, 'SizeData',20, ...
                        'MarkerFaceColor',nh3_colorset{3}, 'MarkerFaceAlpha',0.99, ...
                        'MarkerEdgeColor','none', 'MarkerEdgeAlpha',0.99 );
            end
        end
        % Out-of-plume receptor statistics
        for i = 1 : size(nh3_out_stat,2)
            if i == 1
                nh3_obj(6) = scatter(nh3_xdist(:,i), nh3_out_stat(:,i), ...
                            'Marker','o', 'LineWidth',3, 'SizeData',20, ...
                            'MarkerFaceColor',nh3_colorset{4}, 'MarkerFaceAlpha',0.99, ...
                            'MarkerEdgeColor','none', 'MarkerEdgeAlpha',0.99, 'LineWidth',2 );
            else
                scatter(nh3_xdist(:,i), nh3_out_stat(:,i), ...
                            'Marker','o', 'LineWidth',3, 'SizeData',20, ...
                            'MarkerFaceColor',nh3_colorset{4}, 'MarkerFaceAlpha',0.99, ...
                            'MarkerEdgeColor','none', 'MarkerEdgeAlpha',0.99, 'LineWidth',2 );
            end
        end
        % Plot averaged out of plume statistics statistics
        nh3_obj(3) = plot([min(nh3_xdist(:)) max(nh3_xdist(:))], [nh3_out_stat_mean nh3_out_stat_mean], ...
                'LineStyle','-', 'LineWidth',3.5, 'Color',nh3_colorset{2});
        nh3_obj(4) = plot([min(nh3_xdist(:)) max(nh3_xdist(:))], [nh3_out_stat_mean + nh3_out_stat_std nh3_out_stat_mean + nh3_out_stat_std], ...
                'LineStyle','--', 'LineWidth',3.5, 'Color',nh3_colorset{2});
        plot([min(nh3_xdist(:)) max(nh3_xdist(:))], [nh3_out_stat_mean - nh3_out_stat_std nh3_out_stat_mean - nh3_out_stat_std], ...
                'LineStyle','--', 'LineWidth',3.5, 'Color',nh3_colorset{2});


        nh3_INblend_centre = find(nh3_in_centreline(:,1) == nh3_blend_centre_abs);
        if isempty(nh3_INblend_centre)
            nh3_INblend_centre = length(nh3_in_centreline(:,1));
        else
            nh3_INblend_centre = nh3_INblend_centre(1);
        end
            % Blending distances
        nh3_obj(2) = plot([nh3_blend_centre_abs, nh3_blend_centre_abs, nh3_blend_centre_abs], ...
            [-9999 nh3_in_centreline(nh3_INblend_centre, 2) 9999], ...
            'LineStyle','-', 'LineWidth',3.5, 'Color','k');%nh3_colorset{1});

            % Centreline and maximum lines
        nh3_obj(1) = plot(nh3_in_centreline(:,1), nh3_in_centreline(:,2), ...
                    'LineStyle','-', 'LineWidth',3.5, 'Color',nh3_colorset{1});

        % Settings
        nh3_dylim = [min([min(nh3_out_stat(:)) min(nh3_in_stat(:))]), ...
            max([max(nh3_out_stat(:)) max(nh3_in_stat(:))])];
        set(nh3_ax1,'FontSize',18)
        xlabel('X-distance from the source [m]');
        ylabel(nh3_ylabel);
        ylim(nh3_ylim)
        ylim([nh3_dylim(1) - 0.05*(nh3_dylim(2) - nh3_dylim(1)), ...
            nh3_dylim(2) + 0.05*(nh3_dylim(2) - nh3_dylim(1))])
        if isnan(check.bs_xlim(1)) || isnan(check.bs_xlim(2))
            xlim([min(nh3_x - nh3_srcpos(1)) max(nh3_x - nh3_srcpos(1))])
        else
            xlim(check.bs_xlim)
        end
        if strcmp(check.bs_xy_xz, 'xz')
            title([ 'XZ, set' nh3_set ': ' nh3_title ' for ' nh3_var ' between ' ...
                num2str(check.bs_twindow(1)) '-' num2str(check.bs_twindow(2)) ' CEST at y = ' ...
                num2str(median(check.bulk_pos_in(:,2))) ' m'])
        else
            title([ 'XY, set' nh3_set ': ' nh3_title ' for ' nh3_var ' between ' ...
                num2str(check.bs_twindow(1)) '-' num2str(check.bs_twindow(2)) ' CEST at ' ...
                num2str(check.bs_xy_z) ' m height'])
        end
        nh3_lgnd = legend(nh3_obj,{'Centreline', '"Blending-distance"', ...
            'Mean_{background}','\sigma_{background}', ...
            'In-plumereceptor','Background receptor'},...
            'box','off');
        nh3_lgnd.Position = [0.78 0.3075 0.21 0.5625];
    end
    
    
     %% Build the NEW absolute distance figure
    
     
    fprintf(['      - Build the absolute distrance Vs. statistics figure\n'])
    if strcmp(nh3_statname, 'I') == 0

        nh3_colormap_Delta = length(nh3_blendthreshold);
        if mod(nh3_colormap_Delta/2, 1) == 0
            nh3_colormap = [[ [094 : (255 - 094)/(nh3_colormap_Delta ) : 255 ]', ...
                             [060 : (255 - 060)/(nh3_colormap_Delta ) : 255 ]', ...
                             [153 : (255 - 153)/(nh3_colormap_Delta ) : 255 ]' ] ./255];
        else
            nh3_colormap = [[ [094 : (255 - 094)/ceil(nh3_colormap_Delta ) : 255 ]', ...
                             [060 : (255 - 060)/ceil(nh3_colormap_Delta ) : 255 ]', ...
                             [153 : (255 - 153)/ceil(nh3_colormap_Delta ) : 255 ]' ] ./255];
        end

        % Define x- and y-lim of the main subplot
            % x-limit
        if isnan(check.bs_xlim(1)) || isnan(check.bs_xlim(2))
            nh3_xlim = [min(nh3_x - nh3_srcpos(1)) max(nh3_x - nh3_srcpos(1))];
        else

            nh3_xlim = check.bs_xlim;
        end
            % y-limit
        nh3_dylim = [min(nh3_norm_stat(:)), max(nh3_norm_stat(:))]; 
        if max( nh3_norm_stat(:) ) > max(nh3_blendthreshold(:)) 
            nh3_ylim = [nh3_dylim(1) - 0.05*(nh3_dylim(2) - nh3_dylim(1)), ...
                nh3_dylim(2) + 0.05*(nh3_dylim(2) - nh3_dylim(1))].*100;
        else
            nh3_ylim = [nh3_dylim(1) - 0.05*(nh3_dylim(2) - nh3_dylim(1)), ...
                max(nh3_blendthreshold(:)) + 0.05 * max(nh3_blendthreshold(:))] .*100;
        end
        % Define positions of subplot & zoom
        nh3_sub_position = [0.11 0.21 0.84 0.77];
        nh3_zoom_position = [nh3_sub_position(1) + nh3_sub_position(3) - 0.27, ...
            nh3_sub_position(2) + nh3_sub_position(4) - 0.50, 0.27, 0.50];
        % Define x- and y-boundarys for the zoom
        if sum(nh3_blend_maxstat == 0) ~= length(nh3_blend_maxstat)
            nh3_zoom = NaN(2, 2);
            nh3_zoom(:,1) = [min(nh3_blend_maxstat(nh3_blend_maxstat ~= 0)) - 100; ...
                max(nh3_blend_maxstat(nh3_blend_maxstat ~= 0)) + 100];
            if nh3_zoom(1,1) < 0
                nh3_zoom(1,1) = 0;
            end
            nh3_temp1 = nh3_pos(:);
            nh3_temp2 = nh3_norm_stat(:);
            nh3_zoom(1,2) = min(nh3_temp2( nh3_temp1 >= nh3_zoom(1,1) & nh3_temp1 <= nh3_zoom(2,1) ) )*100;
            nh3_zoom(2,2) = max(nh3_temp2( nh3_temp1 >= nh3_zoom(1,1) & nh3_temp1 <= nh3_zoom(2,1) ) )*100;
            nh3_Dzoom = nh3_zoom(2,2) - nh3_zoom(1,2);
            nh3_zoom(1,2) = nh3_zoom(1,2) - 0.04* nh3_Dzoom;
            nh3_zoom(2,2) = nh3_zoom(2,2) + 0.0* nh3_Dzoom;
            clear nh3_temp1 nh3_temp2
            % Define zoomplot position on the xy-grid
            nh3_zoomxy = NaN(2,2);
            nh3_zoomxy(1,:) = [nh3_zoom(1,1)+(nh3_zoom(2,1) - nh3_zoom(1,1))/2, nh3_zoom(2,2)];
            nh3_zoomxy(2,:) = [nh3_xlim(1) + (nh3_sub_position(3) - nh3_zoom_position(3)) / nh3_sub_position(3) * ...
                (nh3_xlim(2) - nh3_xlim(1)), ...
                nh3_ylim(1) + (nh3_sub_position(4) - nh3_zoom_position(4)/2) / nh3_sub_position(4) * ...
                (nh3_ylim(2) - nh3_ylim(1))];
        end
        


        figure('units','pixels', 'Color','w',...
            'innerposition', [10 200 check.figure_width 350], ...
            'Name', 'Receptor stats');
        nh3_subplot = subplot('Position',nh3_sub_position);

        grid on      % show grid
        hold on         % show multiple lines
        nh3_ax1 = gca;   % axis variable
        % Receptor statistics
        for i = 1 : size(nh3_norm_stat,2)
            scatter(nh3_pos(:,i)./1000, nh3_norm_stat(:,i).*100, ...
                    'Marker','o', 'LineWidth',3, 'SizeData',20, ...
                    'MarkerFaceColor',nh3_colorset{3}, 'MarkerFaceAlpha',1.0, ...
                    'MarkerEdgeColor','none', 'MarkerEdgeAlpha',1.0 );
        end


        nh3_blendthreshold_str = cell(1,length(nh3_blendthreshold));
        for i = 1 : length(nh3_blendthreshold)
            if strcmp(nh3_statname, 'I') == 0
                % Plot the threshold line 
                plot([min(nh3_pos(:)) max(nh3_pos(:))]./1000, [nh3_blendthreshold(i) nh3_blendthreshold(i)].*100, ...
                        'LineStyle','-', 'LineWidth',3.0, 'Color',nh3_colormap(i,:));
                nh3_blendthreshold_str{i} = ['Threshold = ' num2str(nh3_blendthreshold(i)*100) '%'];
            else
                plot([min(nh3_pos(:)) max(nh3_pos(:))]./1000, [0 0],'-k','LineWidth',1.5)
                nh3_blendthreshold_str{i} = [];
            end
        end

        for i = 1 : length(nh3_blendthreshold)
            nh3_INblend_max = find(nh3_norm_maxstat(:,1) == nh3_blend_maxstat(i));
            if isempty(nh3_INblend_max)
                nh3_INblend_max = length(nh3_norm_maxstat(:,1));
            else
                nh3_INblend_max = nh3_INblend_max(1);
            end
%             nh3_INblend_centre = find(nh3_norm_centreline(:,1) == nh3_blend_centre(i));
%             if isempty(nh3_INblend_centre)
%                 nh3_INblend_centre = length(nh3_norm_centreline(:,1));
%             else
%                 nh3_INblend_centre = nh3_INblend_centre(1);
%             end
                % Blending distances
            if nh3_blend_maxstat(i) ~= 0
                plot([nh3_blend_maxstat(i), nh3_blend_maxstat(i), nh3_blend_maxstat(i)]./1000, ...
                    [-9999, nh3_norm_maxstat(nh3_INblend_max,2).*100, 9999], ...
                    'LineStyle','--', 'LineWidth',3.5, 'Color',nh3_colormap(i,:));
            end
        end

        plot(nh3_norm_centreline(:,1)./1000, nh3_norm_centreline(:,2).*100, ...
                    'LineStyle',':', 'LineWidth',3.5, 'Color',nh3_colorset{1});
            % Centreline and maximum lines
        plot(nh3_norm_maxstat(:,1)./1000, nh3_norm_maxstat(:,2).*100, ...
                    'LineStyle','-', 'LineWidth',3.5, 'Color',nh3_colorset{1});
        % Dummy line for blending distance
        plot([-9990 -9999],[0 1], '--k', 'LineWidth',3.5);

        if sum(nh3_blend_maxstat == 0) ~= length(nh3_blend_maxstat)
            % Show zoomed in area
            plot([nh3_zoom(1,1) nh3_zoom(2,1)]./1000, [nh3_zoom(1,2) nh3_zoom(1,2)], ...
                '-k', 'LineWidth', 0.75)
            plot([nh3_zoom(1,1) nh3_zoom(2,1)]./1000, [nh3_zoom(2,2) nh3_zoom(2,2)], ...
                '-k', 'LineWidth', 0.75)
            plot([nh3_zoom(1,1) nh3_zoom(1,1)]./1000, [nh3_zoom(1,2) nh3_zoom(2,2)], ...
                '-k', 'LineWidth', 0.75)
            plot([nh3_zoom(2,1) nh3_zoom(2,1)]./1000, [nh3_zoom(1,2) nh3_zoom(2,2)], ...
                '-k', 'LineWidth', 0.75)
            plot([nh3_zoomxy(1,1), nh3_zoomxy(2,1)]./1000, [nh3_zoomxy(1,2), nh3_zoomxy(2,2)], ...
                '-k', 'LineWidth', 0.75)
        end

        % Settings
        set(nh3_ax1,'FontSize',check.figure_font)
        xlim(nh3_xlim./1000)
        ylim(nh3_ylim)
        xlabel('Absolute distance from the emission source [km]');
        ylabel(nh3_rel_ylabel);


        if sum(nh3_blend_maxstat == 0) ~= length(nh3_blend_maxstat)
            axes('Position',nh3_zoom_position)
            box on; grid on, hold on
            for i = 1 : size(nh3_norm_stat,2)
                scatter(nh3_pos(:,i)./1000, nh3_norm_stat(:,i).*100, ...
                        'Marker','o', 'LineWidth',3, 'SizeData',20, ...
                        'MarkerFaceColor',nh3_colorset{3}, 'MarkerFaceAlpha',1.0, ...
                        'MarkerEdgeColor','none', 'MarkerEdgeAlpha',1.0 );
            end
            nh3_ax2 = gca;   % axis variable
            nh3_ax2.YAxisLocation = 'right';

            for i = 1 : length(nh3_blendthreshold)
                nh3_INblend_max = find(nh3_norm_maxstat(:,1) == nh3_blend_maxstat(i));
                if isempty(nh3_INblend_max)
                    nh3_INblend_max = length(nh3_norm_maxstat(:,1));
                else
                    nh3_INblend_max = nh3_INblend_max(1);
                end
        %         nh3_INblend_centre = find(nh3_norm_centreline(:,1) == nh3_blend_centre(i));
        %         if isempty(nh3_INblend_centre)
        %             nh3_INblend_centre = length(nh3_norm_centreline(:,1));
        %         else
        %             nh3_INblend_centre = nh3_INblend_centre(1);
        %         end
                % Blending distances
                if nh3_blend_maxstat(i) ~= 0
                    plot([nh3_blend_maxstat(i), nh3_blend_maxstat(i), nh3_blend_maxstat(i)]./1000, ...
                            [-9999, nh3_norm_maxstat(nh3_INblend_max,2).*100, 9999], ...
                            'LineStyle','--', 'LineWidth',3.5, 'Color',nh3_colormap(i,:));
                end
        %         plot([nh3_blend_centre(i), nh3_blend_centre(i), nh3_blend_centre(i)], ...
        %             [-9999, nh3_norm_centreline(nh3_INblend_centre,2).*100, 9999], ...
        %             'LineStyle',':', 'LineWidth',3.5, 'Color',nh3_colormap(i,:));%nh3_colorset{1});
            end

            nh3_blendthreshold_str = cell(1,length(nh3_blendthreshold));
            for i = 1 : length(nh3_blendthreshold)
                if strcmp(nh3_statname, 'I') == 0
                    % Plot the threshold line 
                    plot([min(nh3_pos(:)) max(nh3_pos(:))]./1000, [nh3_blendthreshold(i) nh3_blendthreshold(i)].*100, ...
                            'LineStyle','-', 'LineWidth',3.0, 'Color',nh3_colormap(i,:));
                    nh3_blendthreshold_str{i} = ['Threshold = ' num2str(nh3_blendthreshold(i)*100) '%'];
                else
                    plot([min(nh3_pos(:)) max(nh3_pos(:))]./1000, [0 0],'-k','LineWidth',1.5)
                    nh3_blendthreshold_str{i} = [];
                end
            end

            plot(nh3_norm_centreline(:,1)./1000, nh3_norm_centreline(:,2).*100, ...
                        'LineStyle',':', 'LineWidth',3.5, 'Color',nh3_colorset{1});
                % Centreline and maximum lines
            plot(nh3_norm_maxstat(:,1)./1000, nh3_norm_maxstat(:,2).*100, ...
                        'LineStyle','-', 'LineWidth',3.5, 'Color',nh3_colorset{1});

            for i = 1 : length(nh3_blendthreshold)
                nh3_INblend_max = find(nh3_norm_maxstat(:,1) == nh3_blend_maxstat(i));
                if isempty(nh3_INblend_max)
                    nh3_INblend_max = length(nh3_norm_maxstat(:,1));
                else
                    nh3_INblend_max = nh3_INblend_max(1);
                end
                if nh3_blend_maxstat(i) ~= 0
                    scatter(nh3_blend_maxstat(i)./1000, nh3_norm_maxstat(nh3_INblend_max,2).*100, 'o', ...
                        'MarkerFaceColor',nh3_colormap(i,:), 'MarkerEdgeColor',nh3_colormap(i,:), ...
                        'Sizedata',45)
                end
            end

             % Settings
            set(nh3_ax2,'FontSize',check.figure_font)
            xlim(nh3_zoom(:,1)./1000)
            ylim(nh3_zoom(:,2))
            nh3_temp = nh3_ax2.YTickLabel;
            nh3_temp2 = cell(size(nh3_temp));
            for i = 1 : length(nh3_temp)
                if str2num(nh3_temp{i}) >= 0
                    nh3_temp2(i) = nh3_temp(i);
                else
                    nh3_temp2{i} = '';
                end
            end
            nh3_ax2.YTickLabel = nh3_temp2;
            clear nh3_temp nh3_temp2
            nh3_temp = nh3_ax2.XTickLabel;
            nh3_ax2.XTickLabel = nh3_temp;
            clear nh3_temp
        end
    end
    
    %% In case of intermittency, build a propper absolute distance figure
    
    if strcmp(nh3_statname, 'I')
        nh3_colormap_Delta = length(nh3_blendthreshold);
        if mod(nh3_colormap_Delta/2, 1) == 0
            nh3_colormap = [[ [094 : (255 - 094)/(nh3_colormap_Delta ) : 255 ]', ...
                             [060 : (255 - 060)/(nh3_colormap_Delta ) : 255 ]', ...
                             [153 : (255 - 153)/(nh3_colormap_Delta ) : 255 ]' ] ./255];
        else
            nh3_colormap = [[ [094 : (255 - 094)/ceil(nh3_colormap_Delta ) : 255 ]', ...
                             [060 : (255 - 060)/ceil(nh3_colormap_Delta ) : 255 ]', ...
                             [153 : (255 - 153)/ceil(nh3_colormap_Delta ) : 255 ]' ] ./255];
        end

        % Define x- and y-lim of the main subplot
        nh3_dylim = [min(nh3_norm_stat(:)), max(nh3_norm_stat(:))];
        nh3_ylim = [nh3_dylim(1) - 0.05*(nh3_dylim(2) - nh3_dylim(1)), ...
            nh3_dylim(2) + 0.05*(nh3_dylim(2) - nh3_dylim(1))];
        if isnan(check.bs_xlim(1)) || isnan(check.bs_xlim(2))
            nh3_xlim = [min(nh3_x - nh3_srcpos(1)) max(nh3_x - nh3_srcpos(1))];
        else
            nh3_xlim = check.bs_xlim;
        end
        % Define positions of subplot & zoom
        nh3_sub_position = [0.11 0.21 0.84 0.77];
        nh3_zoom_position = [nh3_sub_position(1) + nh3_sub_position(3) - 0.27, ...
            nh3_sub_position(2) + nh3_sub_position(4) - 0.50, 0.27, 0.50];
        % Define x- and y-boundarys for the zoom
        nh3_temp1 = nh3_pos(:);
        nh3_temp2 = nh3_norm_stat(:);
        nh3_zoom = NaN(2, 2);
        nh3_zoom(:,1) = [min(nh3_blend_maxstat) - 550; max(nh3_blend_maxstat) + 50];
        nh3_zoom(1,2) = -0.01;
        nh3_zoom(2,2) = max(nh3_temp2( nh3_temp1 >= nh3_zoom(1,1) & nh3_temp1 <= nh3_zoom(2,1) ) )+0.002;
        clear nh3_temp1 nh3_temp2
        % Define zoomplot position on the xy-grid
        nh3_zoomxy = NaN(2,2);
        nh3_zoomxy(1,:) = [nh3_zoom(1,1)+(nh3_zoom(2,1) - nh3_zoom(1,1))/2, nh3_zoom(2,2)];
        nh3_zoomxy(2,:) = [nh3_xlim(1) + (nh3_sub_position(3) - nh3_zoom_position(3)) / nh3_sub_position(3) * ...
            (nh3_xlim(2) - nh3_xlim(1)), ...
            nh3_ylim(1) + (nh3_sub_position(4) - nh3_zoom_position(4)/2) / nh3_sub_position(4) * ...
            (nh3_ylim(2) - nh3_ylim(1))];




        figure('units','pixels', 'Color','w',...
            'innerposition', [10 200 check.figure_width 350], ...
            'Name', 'Receptor stats');
        nh3_subplot = subplot('Position',nh3_sub_position);

        grid on      % show grid
        hold on         % show multiple lines
        nh3_ax1 = gca;   % axis variable
        % In-plume Receptor statistics
        for i = 1 : size(nh3_norm_stat,2)
            scatter(nh3_pos(:,i)./1000, nh3_norm_stat(:,i), ...
                    'Marker','o', 'LineWidth',3, 'SizeData',20, ...
                    'MarkerFaceColor',nh3_colorset{3}, 'MarkerFaceAlpha',1.0, ...
                    'MarkerEdgeColor','none', 'MarkerEdgeAlpha',1.0 );
        end


        nh3_blendthreshold_str = cell(1,length(nh3_blendthreshold));
        for i = 1 : length(nh3_blendthreshold)
            if strcmp(nh3_statname, 'I') == 0
                % Plot the threshold line 
                plot([min(nh3_pos(:)) max(nh3_pos(:))]./1000, [nh3_blendthreshold(i) nh3_blendthreshold(i)], ...
                        'LineStyle','-', 'LineWidth',3.0, 'Color',nh3_colormap(i,:));
                nh3_blendthreshold_str{i} = ['Threshold = ' num2str(nh3_blendthreshold(i)*100) '%'];
            else
                plot([min(nh3_pos(:)) max(nh3_pos(:))]./1000, [0 0],'-k','LineWidth',1.5)
                nh3_blendthreshold_str{i} = [];
            end
        end

        for i = 1 : length(nh3_blendthreshold)
            nh3_INblend_max = find(nh3_norm_maxstat(:,1) == nh3_blend_maxstat(i));
            if isempty(nh3_INblend_max)
                nh3_INblend_max = length(nh3_norm_maxstat(:,1));
            else
                nh3_INblend_max = nh3_INblend_max(1);
            end
        end
        nh3_INblend_max = nh3_INblend_max(1);
            % Blending distances
        plot([nh3_blend_maxstat(i), nh3_blend_maxstat(i), nh3_blend_maxstat(i)]./1000, ...
            [-9999, nh3_norm_maxstat(nh3_INblend_max,2), 9999], ...
            'LineStyle','--', 'LineWidth',3.5, 'Color','k');

        plot(nh3_norm_centreline(:,1)./1000, nh3_norm_centreline(:,2), ...
                    'LineStyle',':', 'LineWidth',3.5, 'Color',nh3_colorset{1});
            % Centreline and maximum lines
        plot(nh3_norm_maxstat(:,1)./1000, nh3_norm_maxstat(:,2), ...
                    'LineStyle','-', 'LineWidth',3.5, 'Color',nh3_colorset{1});
        % Dummy line for blending distance
        plot([-9990 -9999],[0 1], '--k', 'LineWidth',3.5);

        % Show zoomed in area
        plot([nh3_zoom(1,1) nh3_zoom(2,1)]./1000, [nh3_zoom(1,2) nh3_zoom(1,2)], ...
            '-k', 'LineWidth', 0.75)
        plot([nh3_zoom(1,1) nh3_zoom(2,1)]./1000, [nh3_zoom(2,2) nh3_zoom(2,2)], ...
            '-k', 'LineWidth', 0.75)
        plot([nh3_zoom(1,1) nh3_zoom(1,1)]./1000, [nh3_zoom(1,2) nh3_zoom(2,2)], ...
            '-k', 'LineWidth', 0.75)
        plot([nh3_zoom(2,1) nh3_zoom(2,1)]./1000, [nh3_zoom(1,2) nh3_zoom(2,2)], ...
            '-k', 'LineWidth', 0.75)
        plot([nh3_zoomxy(1,1), nh3_zoomxy(2,1)]./1000, [nh3_zoomxy(1,2), nh3_zoomxy(2,2)], ...
            '-k', 'LineWidth', 0.75)

        % Settings
        set(nh3_ax1,'FontSize',check.figure_font)
        xlabel('Absolute distance from the emission source [km]');
        ylabel(nh3_rel_ylabel);
        xlim(nh3_xlim./1000)
        ylim(nh3_ylim)


        axes('Position',nh3_zoom_position)
        box on; grid on, hold on
        for i = 1 : size(nh3_norm_stat,2)
            scatter(nh3_pos(:,i)./1000, nh3_norm_stat(:,i), ...
                    'Marker','o', 'LineWidth',3, 'SizeData',20, ...
                    'MarkerFaceColor',nh3_colorset{3}, 'MarkerFaceAlpha',1.0, ...
                    'MarkerEdgeColor','none', 'MarkerEdgeAlpha',1.0 );
        end
        nh3_ax2 = gca;   % axis variable
        nh3_ax2.YAxisLocation = 'right';

        for i = 1 : length(nh3_blendthreshold)
            nh3_INblend_max = find(nh3_norm_maxstat(:,1) == nh3_blend_maxstat(i));
            if isempty(nh3_INblend_max)
                nh3_INblend_max = length(nh3_norm_maxstat(:,1));
            else
                nh3_INblend_max = nh3_INblend_max(1);
            end
        end
        nh3_INblend_max = nh3_INblend_max(1);
            % Blending distances
        plot([nh3_blend_maxstat(i), nh3_blend_maxstat(i), nh3_blend_maxstat(i)]./1000, ...
            [-9999, nh3_norm_maxstat(nh3_INblend_max,2), 9999], ...
            'LineStyle','--', 'LineWidth',3.5, 'Color','k');

        nh3_blendthreshold_str = cell(1,length(nh3_blendthreshold));
        for i = 1 : length(nh3_blendthreshold)
            if strcmp(nh3_statname, 'I') == 0
                % Plot the threshold line 
                plot([min(nh3_pos(:)) max(nh3_pos(:))]./1000, [nh3_blendthreshold(i) nh3_blendthreshold(i)], ...
                        'LineStyle','-', 'LineWidth',3.0, 'Color',nh3_colormap(i,:));
                nh3_blendthreshold_str{i} = ['Threshold = ' num2str(nh3_blendthreshold(i)*100) '%'];
            else
                plot([min(nh3_pos(:)) max(nh3_pos(:))]./1000, [0 0],'-k','LineWidth',1.5)
                nh3_blendthreshold_str{i} = [];
            end
        end

        plot(nh3_norm_centreline(:,1)./1000, nh3_norm_centreline(:,2), ...
                    'LineStyle',':', 'LineWidth',3.5, 'Color',nh3_colorset{1});
            % Centreline and maximum lines
        plot(nh3_norm_maxstat(:,1)./1000, nh3_norm_maxstat(:,2), ...
                    'LineStyle','-', 'LineWidth',3.5, 'Color',nh3_colorset{1});

         % Settings
        set(nh3_ax2,'FontSize',check.figure_font)
        xlim(nh3_zoom(:,1)./1000)
        ylim(nh3_zoom(:,2))
        nh3_temp = nh3_ax2.YTickLabel;
        nh3_temp2 = cell(size(nh3_temp));
        for i = 1 : length(nh3_temp)
            if str2num(nh3_temp{i}) >= 0
                nh3_temp2(i) = nh3_temp(i);
            else
                nh3_temp2{i} = '';
            end
        end
        nh3_ax2.YTickLabel = nh3_temp2;
        clear nh3_temp nh3_temp2
        nh3_temp = nh3_ax2.XTickLabel;
        nh3_ax2.XTickLabel = nh3_temp;
        clear nh3_temp

    end
    
    
end
%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%          blending-distance comparisons            %%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function bd_blendingdist = f__FIGsafedistance(outp, check, bd_compare )
    
    bd_colorset = {[230,097,001]./255, [094,060,153]./255,...
                   [253,184,099]./255, [178,171,210]./255};
    
    % Define variables for each sensitivity study in the figure
    if strcmp(bd_compare, 'variables')
        
        % Runs that will be used
        bd_var      = cell(size(check.bd_var,1),1);
        bd_runs     = cell(size(check.bd_var,1),1);
        bd_zheight  = cell(size(check.bd_var,1),1);
        bd_twindow	= cell(size(check.bd_var,1),1);
        bd_avgt     = cell(size(check.bd_var,1),1);
        for i = 1 : size(check.bd_var,1)
            % Variables which change with each sensitivity study
            bd_var{i} = check.bd_var{i,2};
            bd_runs{i} = ['set' check.bd_var{i,1}];
            % Variables that are the same for each sensitivity study
            bd_zheight{i} = check.bd_zheight(1);
            bd_twindow{i} = check.bd_twindow(1,:);
            bd_avgt{i} = check.bd_avgt(1);
            
             if round(outp.(bd_runs{i}).time(5) - outp.(bd_runs{i}).time(4), 6) ...
                    > round(bd_avgt{i}, 6)
                % The averaging time is smaller then the time-steps. Make
                % them equal
                bd_avgt{i} = round(outp.(bd_runs{i}).time(5) - outp.(bd_runs{i}).time(4), 8);
            end
        end
        
        bd_yticks = [1 : size(check.bd_var,1)];
        bd_ytickslabel = cell(size(check.bd_var,1), 1);
        for i = 1 : length(bd_yticks)
            bd_ytickslabel{i} = [bd_runs{i} ' ' bd_var{i}];
        end
        bd_ylabel = 'Variables';
        
    elseif strcmp(bd_compare, 'twindow')
        
        % Runs that will be used
        bd_var      = cell(size(check.bd_twindow, 1) ,1);
        bd_runs     = cell(size(check.bd_twindow, 1) ,1);
        bd_zheight  = cell(size(check.bd_twindow, 1) ,1);
        bd_twindow  = cell(size(check.bd_twindow, 1) ,1);
        bd_avgt     = cell(size(check.bd_twindow, 1) ,1);
        for i = 1 : size(check.bd_twindow, 1)
            % Variables which change with each sensitivity study
            bd_twindow{i} = check.bd_twindow(i,:);
            % Variables that are the same for each sensitivity study
            bd_var{i} = check.bd_var{1,2};
            bd_runs{i} = ['set' check.bd_var{1,1}];
            bd_zheight{i} = check.bd_zheight(1);
            bd_avgt{i} = check.bd_avgt(1);
            
            % Make sure the averaging time is not smaller than the time steps
            if round(outp.(bd_runs{i}).time(5) - outp.(bd_runs{i}).time(4), 6) ...
                    > round(bd_avgt{i}, 6)
                % The averaging time is smaller then the time-steps. Make
                % them equal
                bd_avgt{i} = round(outp.(bd_runs{i}).time(5) - outp.(bd_runs{i}).time(4), 8);
            end
        end
        
        bd_yticks = [1 : size(check.bd_twindow,1)];
        bd_ytickslabel = cell(size(check.bd_twindow,1), 1);
        for i = 1 : length(bd_yticks)
            bd_ytickslabel{i} = [num2str(bd_twindow{i}(1)) '-' num2str(bd_twindow{i}(2)) ' CEST'];
        end
        bd_ylabel = 'Analysis window';
        
    elseif strcmp(bd_compare, 'tavgt')
        
        % Runs that will be used
        bd_var      = cell(size(check.bd_avgt, 1) ,1);
        bd_runs     = cell(size(check.bd_avgt, 1) ,1);
        bd_zheight  = cell(size(check.bd_avgt, 1) ,1);
        bd_twindow  = cell(size(check.bd_avgt, 1) ,1);
        bd_avgt     = cell(size(check.bd_avgt, 1) ,1);
        for i = 1 : size(check.bd_avgt, 1)
            % Variables which change with each sensitivity study
            bd_avgt{i} = check.bd_avgt(i);
            % Variables that are the same for each sensitivity study
            bd_var{i} = check.bd_var{1,2};
            bd_runs{i} = ['set' check.bd_var{1,1}];
            bd_twindow{i} = check.bd_twindow(1,:);
            bd_zheight{i} = check.bd_zheight(1);
            
            % Make sure the averaging time is not smaller than the time steps
            if round(outp.(bd_runs{i}).time(5) - outp.(bd_runs{i}).time(4), 6) ...
                    > round(bd_avgt{i}, 6)
                % The averaging time is smaller then the time-steps. Make
                % them equal
                bd_avgt{i} = round(outp.(bd_runs{i}).time(5) - outp.(bd_runs{i}).time(4), 8);
            end
        end
        
        bd_yticks = [1 : size(check.bd_avgt,1)];
        bd_ytickslabel = cell(size(check.bd_avgt,1), 1);
        for i = 1 : length(bd_yticks)
            bd_ytickslabel{i} = num2str(bd_avgt{i}*60);
        end
        bd_ylabel = 'Averaging time [min]';
        
    elseif strcmp(bd_compare, 'zheight')
        
        temp_zheight = NaN(size(check.bd_zheight, 1), 1);
        for i = 1 : length(temp_zheight)
            [~, temp_INz] = sort(abs( outp.(['set' check.bd_var{1,1}]).zpos - check.bd_zheight(i) ));
            temp_zheight(i) = outp.(['set' check.bd_var{1,1}]).zpos(temp_INz(1));
        end
        temp_zheight = unique(temp_zheight);
        
        % Runs that will be used
        bd_var      = cell(size(temp_zheight, 1) ,1);
        bd_runs     = cell(size(temp_zheight, 1) ,1);
        bd_zheight  = cell(size(temp_zheight, 1) ,1);
        bd_twindow  = cell(size(temp_zheight, 1) ,1);
        bd_avgt     = cell(size(temp_zheight, 1) ,1);
        for i = 1 : size(temp_zheight, 1)
            % Variables which change with each sensitivity study
            bd_zheight{i} = temp_zheight(i);
            % Variables that are the same for each sensitivity study
            bd_var{i} = check.bd_var{1,2};
            bd_runs{i} = ['set' check.bd_var{1,1}];
            bd_twindow{i} = check.bd_twindow(1,:);
            bd_avgt{i} = check.bd_avgt(1);
            
            % Make sure the averaging time is not smaller than the time steps
            if round(outp.(bd_runs{i}).time(5) - outp.(bd_runs{i}).time(4), 6) ...
                    > round(bd_avgt{i}, 6)
                % The averaging time is smaller then the time-steps. Make
                % them equal
                bd_avgt{i} = round(outp.(bd_runs{i}).time(5) - outp.(bd_runs{i}).time(4), 8);
            end
        end
        
        bd_yticks = [1 : size(temp_zheight,1)];
        bd_ytickslabel = cell(size(temp_zheight,1), 1);
        for i = 1 : length(bd_yticks)
            bd_ytickslabel{i} = num2str(bd_zheight{i});
        end
        bd_ylabel = 'Receptor height [m]';
        
    end
        
        
    % Preallocate the blending-distance for each sensitivity study
    bd_blendingdist = NaN(size(bd_var,1),length(check.bd_dist_level));
    % Calculate the blending-distance for each sensitivity study
    for II = 1 : length(bd_var)
        % Determine the source position
        bd_srcpos_x = [];
        for i = 1 : length(outp.(bd_runs{II}).yt)
            if sum(outp.(bd_runs{II}).srcpos.(bd_var{II})(:,i) == 1) > 0
                bd_srcpos_x = [bd_srcpos_x ; outp.(bd_runs{II}).xt(...
                    outp.(bd_runs{II}).srcpos.(bd_var{II})(:,i) == 1)];
            end
        end
        bd_srcpos_x = bd_srcpos_x + ...
            0.5*(outp.(bd_runs{II}).xt(5) - outp.(bd_runs{II}).xt(4));
        bd_srcpos_x = max(bd_srcpos_x);
        bd_srcpos_y = [];
        for i = 1 : length(outp.(bd_runs{II}).xt)
            if sum(outp.(bd_runs{II}).srcpos.(bd_var{II})(i,:) == 1) > 0
                bd_srcpos_y = [bd_srcpos_y ; outp.(bd_runs{II}).yt(...
                    outp.(bd_runs{II}).srcpos.(bd_var{II})(i,:) == 1)];
            end
        end
        [~, bd_tempy] = sort( abs(bd_srcpos_y - median(bd_srcpos_y)) );
        bd_srcpos_y = bd_srcpos_y(bd_tempy(1));
        bd_srcpos = [bd_srcpos_x, bd_srcpos_y];
        clear bd_srcpos_x bd_srcpos_y bd_tempy

        % Determine the to-be-loaded indices
        bd_INt = find( round( outp.(bd_runs{II}).time, 6)   >  round( bd_twindow{II}(1), 6)-1 & ...
            round( outp.(bd_runs{II}).time, 6) <= round( bd_twindow{II}(2), 6));
        bd_t = outp.(bd_runs{II}).time(bd_INt);
        bd_dt = median(bd_t(2:end) - bd_t(1:end-1));
        [~,bd_INtemp] = sort( abs(bd_zheight{II} - outp.(bd_runs{II}).zpos) );
        bd_crossheigt = outp.(bd_runs{II}).cross_lvl{bd_INtemp(1)};
        bd_posx = outp.(bd_runs{II}).xt;
        clear bd_INtemp
        
        % In case of intermittency, change the variable to the plume
        % variable
        if strcmp(check.bd_stat, 'I')
            bd_new_var = bd_var{II}(4:end);
            if strcmp(bd_new_var(1), '_')
                bd_new_var(1) = [];
                if strcmp(bd_new_var(end - length(check.variable_plume_indicator)+1: end), check.variable_plume_indicator)
                    bd_new_var(end - length(check.variable_plume_indicator)+1: end) = [];
                elseif strcmp(bd_new_var(end - length(check.variable_background_indicator)+1: end), check.variable_background_indicator)
                    bd_new_var(end - length(check.variable_background_indicator)+1: end) = [];
                end
            end
            bd_var{II} = ['nh3_' bd_new_var check.variable_plume_indicator];
        end
        
        
        % Get the data from the .mat files
        %%%
        % total data
        %%%
        
        % Get the raw data
        bd_matfile = matfile( [ ...
            outp.(bd_runs{II}).loc bd_runs{II}(4:end) 'crossxy' bd_crossheigt '_' bd_var{II} '.mat' ....
            ] );

        clear bd_C
        bd_time_new = bd_t;
        if round(bd_avgt{II}, 8) > round(bd_dt, 8)
            % Define new time variable based on the new averaging time
            bd_time_new = [bd_t(1) - bd_dt+bd_avgt{II} : bd_avgt{II} : bd_t(end)];
            bd_C = NaN(length(outp.(bd_runs{II}).xt), ...
                       length(outp.(bd_runs{II}).yt), ...
                       length(bd_INt)/(bd_avgt{II}*3600));

            bd_percentage = 0.20;
            bd_percentage_base = bd_percentage;

            for i = 1 : length(bd_time_new)
                if i == 1 
                    bd_C(:,:,i) = mean(bd_matfile.(bd_var{II})(:,:, ...
                        bd_INt(bd_t <= bd_time_new(i))), 3);
                else
                    bd_C(:,:,i) = mean(bd_matfile.(bd_var{II})(:,:, ...
                        bd_INt(bd_t > bd_time_new(i-1) & bd_t <= bd_time_new(i))), 3);
                end

                if i/length(bd_time_new) >= bd_percentage_base

                    fprintf(['         ' bd_runs{II} '-' bd_var{II} ' data at '...
                        num2str(bd_percentage_base * 100) ' %%\n'])
                    bd_percentage_base = bd_percentage_base + bd_percentage;
                end
            end
            clear bd_percentage_base bd_percentage i
        elseif round(bd_avgt{II}, 8) == round(bd_dt, 8)

            bd_C = bd_matfile.(bd_var{II})(:,:, bd_INt);
        end

        %%%
        % background data
        %%%
        if strcmp(check.bd_stat, 'I') == 0
            % Define the bavkground variable
            bd_bg_var = bd_var{II}(4:end);
            if strcmp(bd_bg_var(1), '_')
                bd_bg_var(1) = [];
                if strcmp(bd_bg_var(end - length(check.variable_plume_indicator)+1: end), check.variable_plume_indicator)
                    bd_bg_var(end - length(check.variable_plume_indicator)+1: end) = [];
                elseif strcmp(bd_bg_var(end - length(check.variable_background_indicator)+1: end), check.variable_background_indicator)
                    bd_bg_var(end - length(check.variable_background_indicator)+1: end) = [];
                end
            end
            bd_bg_var = ['nh3_' bd_bg_var check.variable_background_indicator];
            bd_bg_matfile = matfile( [ ...
                outp.(bd_runs{II}).loc bd_runs{II}(4:end) 'crossxy' bd_crossheigt '_' bd_bg_var '.mat' ....
                ] );            

            if round(bd_avgt{II}, 8) > round(bd_dt, 8)
                bd_bg_C = NaN(size(bd_C));

                bd_percentage = 0.20;
                bd_percentage_base = bd_percentage;

                for i = 1 : length(bd_time_new)
                    if i == 1 
                        bd_bg_C(:,:,i) = mean(bd_bg_matfile.(bd_bg_var)(:,:, ...
                            bd_INt(bd_t <= bd_time_new(i))), 3);
                    else
                        bd_bg_C(:,:,i) = mean(bd_bg_matfile.(bd_bg_var)(:,:, ...
                            bd_INt(bd_t > bd_time_new(i-1) & bd_t <= bd_time_new(i))), 3);
                    end

                    if i/length(bd_time_new) >= bd_percentage_base

                        fprintf(['         ' bd_runs{II} '-' bd_bg_var ' data at '...
                            num2str(bd_percentage_base * 100) ' %%\n'])
                        bd_percentage_base = bd_percentage_base + bd_percentage;
                    end
                end
                clear bd_percentage_base bd_percentage i
            elseif round(bd_avgt{II}, 8) == round(bd_dt, 8)

                bd_bg_C = bd_bg_matfile.(bd_bg_var)(:,:, bd_INt);
            end
        end
        
        %%%
        % w data
        %%%
        % Define centerline w data for flux
        if strcmp(check.bd_stat, 'F')

            bd_w_matfile = matfile( [ ...
                outp.(bd_runs{II}).loc bd_runs{II}(4:end) 'crossxy' bd_crossheigt '_w.mat' ....
                ] );

            if round(bd_avgt{II}, 8) > round(bd_dt, 8)
                bd_w = NaN(size(bd_C));

                bd_percentage = 0.20;
                bd_percentage_base = bd_percentage;

                for i = 1 : length(bd_time_new)
                    if i == 1 
                        bd_w(:,:,i) = mean(bd_w_matfile.w(:,:, ...
                            bd_INt(bd_t <= bd_time_new(i))), 3);
                    else
                        bd_w(:,:,i) = mean(bd_w_matfile.w(:,:, ...
                            bd_INt(bd_t > bd_time_new(i-1) & bd_t <= bd_time_new(i))), 3);
                    end

                    if i/length(bd_time_new) >= bd_percentage_base

                        fprintf(['         ' bd_runs{II} '-w data at '...
                            num2str(bd_percentage_base * 100) ' %%\n'])
                        bd_percentage_base = bd_percentage_base + bd_percentage;
                    end
                end
                clear bd_percentage_base bd_percentage i
            elseif round(bd_avgt{II}, 8) == round(bd_dt, 8)

                bd_w = bd_w_matfile.w(:,:, bd_INt);
            end
        end
        
        
        
        % Calculate moving average concentration, if needed
        if ismember({check.bd_stat}, {'std', 'fI', 'S', 'K'})
            bd_averaging_time_h = 1;
            bd_Cma = movmean(bd_C, [bd_averaging_time_h / (bd_time_new(4) - bd_time_new(3)), 0],3);
            bd_bg_Cma = movmean(bd_bg_C, [bd_averaging_time_h / (bd_time_new(4) - bd_time_new(3)), 0],3);
           clear nh3_averaging_time_h
        end
        % Now apply the correct time window (previous one was needed for moving
        % average calculation)
        bd_INt_correct = find( bd_time_new > bd_twindow{II}(1) & bd_time_new <= bd_twindow{II}(2) );
        bd_time_new = bd_time_new(bd_INt_correct);
        bd_C = bd_C(:,:,bd_INt_correct);
        bd_bg_C = bd_bg_C(:,:,bd_INt_correct);
        if ismember({check.bd_stat}, {'std', 'fI', 'S', 'K'})
            bd_Cma = bd_Cma(:,:,bd_INt_correct);
            bd_bg_Cma = bd_bg_Cma(:,:,bd_INt_correct);
        elseif ismember({check.bd_stat}, {'F'})
            bd_w = bd_w(:,:,bd_INt_correct);
        end
        
        % Calculate the centerline statistics
        if strcmp(check.bd_stat, 'mean')
            bd_total_stat =  mean( bd_C ,3);
            bd_bg_stat = mean( bd_bg_C, 3);

        elseif strcmp(check.bd_stat, 'std')
            bd_total_stat =  std( bd_C-bd_Cma , 0, 3);
            bd_bg_stat =  std( bd_bg_C-bd_bg_Cma , 0, 3);

        elseif strcmp(check.bd_stat, 'fI')
            bd_total_stat = std( bd_C-bd_Cma ,0,3) ./ mean( bd_C ,3);
            % In case of a plume variable, there are no fluctuations outside the plume
            if strcmp(bd_var{II}(end - length(check.variable_plume_indicator)+1: end), ...
                    check.variable_plume_indicator)
                bd_bg_stat = zeros(size(bd_total_stat));
            else
                bd_bg_stat = std( bd_bg_C-bd_bg_Cma ,0,3) ./ mean( bd_bg_C ,3);
            end

        elseif strcmp(check.bd_stat, 'I') 
            bd_total_stat = sum( bd_C >= check.bd_Ithreshhold ,3) ./ size(bd_C,3);
            bd_bg_stat = 0;

        elseif strcmp(check.bd_stat, 'S')
            bd_total_stat =  skewness( bd_C-bd_Cma ,1,3);
            bd_bg_stat =  skewness( bd_bg_C-bd_bg_Cma ,1,3);

        elseif strcmp(check.bd_stat, 'K') 
            bd_total_stat =  kurtosis( bd_C-bd_Cma ,1,3);
            bd_bg_stat =  kurtosis( bd_bg_C-bd_bg_Cma ,1,3);
        elseif strcmp(check.bd_stat, 'F')
            
              % Calculate the flux each 30 minutes
            bd_fluxtime = [bd_twindow{II}(1) : 0.5 : bd_twindow{II}(2)];
            bd_flux = NaN(size(bd_C,1), size(bd_C,2), length(bd_fluxtime)-1);
            bd_bg_flux = NaN(size(bd_C,1), size(bd_C,2), length(bd_fluxtime)-1);
            for i = 1 : length(bd_fluxtime)-1
                if i == 1
                    bd_fluxIN = find(bd_time_new <= bd_fluxtime(i+1));
                else
                    bd_fluxIN = find(bd_time_new > bd_fluxtime(i) & ...
                        bd_time_new <= bd_fluxtime(i+1));
                end
                bd_flux(:,:,i) = nanmean(( bd_w(:,:,bd_fluxIN) - nanmean(bd_w(:,:,bd_fluxIN), 3) ) ...
                    .* (bd_C(:,:,bd_fluxIN) - nanmean(bd_C(:,:,bd_fluxIN), 3)), 3);
                bd_bg_flux(:,:,i) = nanmean(( bd_w(:,:,bd_fluxIN) - nanmean(bd_w(:,:,bd_fluxIN), 3) ) ...
                    .* (bd_bg_C(:,:,bd_fluxIN) - nanmean(bd_bg_C(:,:,bd_fluxIN), 3)), 3);

            end
            % Calculate the averaged flux within the time window
            bd_total_stat =  mean(bd_flux, 3);
            bd_bg_stat = mean(bd_bg_flux, 3);
            

            clear bd_flux bd_bg_flux bd_fluxIN bd_fluxtime i
        end
        
        
        % Calculate the normalized statistics
        if strcmp(check.bd_stat, 'I')  == 0
            
            bd_norm_stat = abs( (bd_total_stat - bd_bg_stat) ./ bd_bg_stat);
        else
            bd_norm_stat = bd_total_stat;
        end
        
        
        % Determine the threshold value for the blending-distance
        bd_blendthreshold = NaN(length(check.bd_dist_level),1);
        for i = 1 : length(bd_blendthreshold)
        	if strcmp(check.bd_stat, 'I') == 0
                if strcmp(check.bd_peak_abs, 'peak')
                    bd_blendthreshold(i) = check.bd_dist_level(i) * max(bd_norm_stat(:));
                elseif strcmp(check.bd_peak_abs, 'abs')
                    bd_blendthreshold(i) = check.bd_dist_level(i);
                else
                    bd_blendthreshold(i) = 0;
                end
            else
                bd_blendthreshold(i) = 0;
            end
        end
        

        % Find the values at the plume centerline
        bd_centreline = [ abs(outp.(bd_runs{II}).xt - bd_srcpos(1)) , ...
            bd_norm_stat(:, round(outp.(bd_runs{II}).yt,6) == round(bd_srcpos(2),6))];
    
        % Calculate absolute distance
        bd_pos = NaN(length(outp.(bd_runs{II}).xt), length(outp.(bd_runs{II}).yt));
        for i = 1 : length(outp.(bd_runs{II}).xt)
            bd_pos(i,:) = sqrt( (outp.(bd_runs{II}).xt(i) - bd_srcpos(1)).^2 + ...
                (outp.(bd_runs{II}).yt - bd_srcpos(2)).^2 );
        end
        
        % Find the maximum values of the normalized statistics for absolute distance from the source          
        bd_statmax = bd_norm_stat(:);
        bd_pos_statmax = bd_pos(:);
        bd_pos_statmax(isnan(bd_statmax)) = [];
        bd_statmax(isnan(bd_statmax)) = [];
        
        temp_window = 50;       % Window of XX m within the maximum value is determined
        temp_maxplot = NaN(size(bd_statmax));
        temp_maxdist = NaN(size(bd_statmax));
        for i = 1 : length(bd_statmax)
            INmaxdist = find( bd_pos_statmax >= bd_pos_statmax(i)-temp_window/2 & ...
                bd_pos_statmax <= bd_pos_statmax(i)+temp_window/2 );

            [temp_maxplot(i), temp_IN] = max( bd_statmax(INmaxdist) );
            temp_maxdist(i) = bd_pos_statmax(INmaxdist(temp_IN));

        end
        bd_maxplot = unique(temp_maxplot);
        bd_maxdist = NaN(size(bd_maxplot));
        for i = 1 : length(bd_maxplot)
            bd_maxdist(i) = max(temp_maxdist(temp_maxplot == bd_maxplot(i)));
        end
        [bd_maxdist, INsort] = sort(bd_maxdist);
        bd_maxplot = bd_maxplot(INsort);
        bd_maxstat = [bd_maxdist, bd_maxplot];
        if strcmp(check.bd_stat, 'I') 
            if bd_maxstat(end,2) == 0
                bd_maxstat(end,1) = bd_maxstat(end-1,1);
            end
        end
        clear bd_statmax bd_pos_statmax temp_* INmaxdist INsort
        
        
        for i = 1 : length(bd_blendthreshold)
            %%%%%
            % Calculate blending distance for normalized statistics, Maximum
            %%%%%
            if max(bd_maxstat(:,2)) < bd_blendthreshold(i)
                bd_blend_maxstat = 0;
            else
                % Find the local maximum of the in-plume statistics within the first 1000 m
                bd_INmax = find( bd_maxstat(:,2) == max(bd_maxstat(bd_maxstat(:,1) <= 500, 2) ) );
                % Determine 'blending-distance'
                bd_INblendingdist = find(bd_maxstat(:,2) <= bd_blendthreshold(i));
                if isempty(bd_INblendingdist)
                    bd_blend_maxstat = NaN;
                else
                    % Remove values which are found before the maximum is reached.
                    bd_INblendingdist(bd_INblendingdist < bd_INmax) = [];
                    bd_blend_maxstat = bd_maxstat(bd_INblendingdist(1), 1);
                end
            end
            %%%%%
            % Calculate blending distance for normalized statistics, centreline
            %%%%%
            if max(bd_centreline(:,2)) < bd_blendthreshold(i)
                bd_blend_centre = 0;
            else
                % Find the local maximum of the in-plume statistics within the first 1000 m
                bd_INmax = find( bd_centreline(:,2) == max(bd_centreline(bd_centreline(:,1) <= 500, 2) ) );
                % Determine 'blending-distance'
                bd_INblendingdist = find(bd_centreline(:,2) <= bd_blendthreshold(i));
                if isempty(bd_INblendingdist)
                    bd_blend_centre = NaN;
                else
                    % Remove values which are found before the maximum is reached.
                    bd_INblendingdist(bd_INblendingdist < bd_INmax) = [];
                    bd_blend_centre = bd_centreline(bd_INblendingdist(1), 1);
                end
            end
            
            if strcmp(check.bd_max_centre, 'max')
                bd_blendingdist(II,i) = bd_blend_maxstat;
            elseif strcmp(check.bd_max_centre, 'centre')
                bd_blendingdist(II,i) = bd_blend_centre;
            end
        end
        
        
        clear bd_norm_stat bd_plume_stat bd_bg_stat bd_blend_maxstat ...
            bd_blend_centre bd_INblendingdist

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        fprintf(['     - Blending-distance calculations at ' num2str(II/length(bd_var)*100) ' %%\n'])      
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
    end         % <-- for II = 1 : length(bd_var)
    
    bd_colormap_Delta = length(check.bd_dist_level);
    bd_colormap_purple = [ [094 : (255 - 094)/(bd_colormap_Delta ) : 255 ]', ...
                     [060 : (255 - 060)/(bd_colormap_Delta ) : 255 ]', ...
                     [153 : (255 - 153)/(bd_colormap_Delta ) : 255 ]' ] ./255;
    bd_colormap_orange = [ [230 : (255 - 230)/(bd_colormap_Delta) : 255-(255 - 230)/(bd_colormap_Delta) ]', ...
                 [097 : (255 - 097)/(bd_colormap_Delta) : 255-(255 - 097)/(bd_colormap_Delta) ]', ...
                 [001 : (255 - 001)/(bd_colormap_Delta) : 255-(255 - 001)/(bd_colormap_Delta) ]'] ./255;
    
    figure('units','pixels', 'Color','w',...
        'innerposition', [10 150 1100 470], ...
        'Name', 'Receptor stats');
    hold on; grid on
    for i = 1 : length(bd_yticks)
        
        if mod(i/2,1) == 0
            bd_colormap = bd_colormap_purple;
        else
            bd_colormap = bd_colormap_orange;
        end 
        
        for j = 1 : length(check.bd_dist_level)
            scatter(bd_blendingdist(i,j), bd_yticks(i), 'SizeData', 80, ...
                'MarkerEdgeColor', bd_colormap(j,:), ...
                'MarkerFaceColor', bd_colormap(j,:))
        end
    end
    set(gca, 'FontSize', 16, 'YTick', bd_yticks, 'YTickLabel', bd_ytickslabel)
    ylabel(bd_ylabel)
    ylim([0 max(bd_yticks)+1])
    if max(bd_blendingdist(:)) > 2000
        xlim([0 max(bd_blendingdist(:))])
    else
        xlim([0 2000])
    end
    xlabel('"blending-distance" [m]')
    title([ '"blending-distance" sensitivity study for different ' bd_ylabel])
    
    
    
end
















%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%    Define x/y limits and ticks    %%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [x_tick, x_tick_label, y_tick, y_tick_label] = f__CLASSfig_XYTicks(x_lim, y_lim)
    
    % XLim
    dx_lim = (x_lim(2) - x_lim(1));

    II = 0;
    ii = 0.5;
    iter_limit = 0;
    dx_check = 0.0001;
    while II == 0
        if iter_limit > 100
            II = 1;
        elseif dx_lim > dx_check && ii == 0
            dx_check = dx_check *2;
            ii = mod(ii + mod(0.5,1),1);
            iter_limit = iter_limit+1;
        elseif dx_lim > dx_check && ii == 0.5

            dx_check = dx_check *5;
            ii = mod(ii + mod(0.5,1),1);
            iter_limit = iter_limit+1;
        else
            II = 1;
        end
    end

    if ii == 0
        dx = dx_check / 25;
    elseif ii == 0.5
        dx = dx_check / 10;
    end
    x_tick = [x_lim(1):dx:x_lim(2)];
    
    II = 0;
    iter_limit = 0;
    while II == 0
        if iter_limit > 100
            II = 1;
        elseif length(x_tick) < 15
            dx = dx / 2;
            x_tick = [x_lim(1):dx:x_lim(2)];
            iter_limit = iter_limit+1;
        elseif length(x_tick) > 25
            dx = dx * 2;
            x_tick = [x_lim(1):dx:x_lim(2)];
            iter_limit = iter_limit+1;
        else
            II = 1;
        end
    end

    if mod(x_tick(1),2) == 1
        if mod(x_tick(end),2) == 1
            IN_label = [2:2:length(x_tick)-1];
        else
            IN_label = [2:2:length(x_tick)];
        end
    else
        IN_label = [1:2:length(x_tick)];
    end
    x_tick_label = cell(length(x_tick),1);
    x_tick_label(IN_label) = num2cell(x_tick(IN_label));
        
    % YLim
    dy_lim = (y_lim(2) - y_lim(1));

    II = 0;
    ii = 0.5;
    iter_limit = 0;
    dy_check = 0.00001;
    while II == 0
        if iter_limit > 100
            II = 1;
        elseif dy_lim > dy_check && ii == 0
            dy_check = dy_check *2;
            ii = mod(ii + mod(0.5,1),1);
            iter_limit = iter_limit+1;
        elseif dy_lim > dy_check && ii == 0.5
            dy_check = dy_check *5;
            ii = mod(ii + mod(0.5,1),1);
            iter_limit = iter_limit+1;
        else
            II = 1;
        end
    end

    if ii == 0
        dy = dy_check / 25;
    elseif ii == 0.5
        dy = dy_check / 10;
    end
    y_tick = [y_lim(1):dy:y_lim(2)];
    
    II = 0;
    iter_limit = 0;
    while II == 0
        if iter_limit > 100
            II = 1;
        elseif length(y_tick) < 20
            dy = dy / 2;
            y_tick = [y_lim(1):dy:y_lim(2)];
            iter_limit = iter_limit+1;
        elseif length(y_tick) > 40
            dy = dy * 2;
            y_tick = [y_lim(1):dy:y_lim(2)];
            iter_limit = iter_limit+1;
        else
            II = 1;
        end
    end

    
    IN_label = [1:4:length(y_tick)];
    y_tick_label = cell(length(y_tick),1);
    for i = 1 : length(IN_label)
        y_tick_label{IN_label(i)} = round(y_tick(IN_label(i)),8);
    end
    
end




