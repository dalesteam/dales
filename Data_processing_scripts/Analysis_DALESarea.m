%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%         Ruben Schulte         %%%%%
%%%%%      DALES area analysis      %%%%%
%%%%%      started: 28-12-2020      %%%%%
%%%%%     restarted: 10-06-2021     %%%%%
%%%%%      changed: 14-06-2021      %%%%%
%%%%%       final: 23-05-2025       %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Note: This is a cleaned up version of Analysis_DALES_base2.m
% This script only covers the area statistics
%

%% General settings
% These variables can be changed by the uses

% Define the relevant (data) directories 
direc.main_directory = 'D:\DALES_data\';
direc.run_name = 'nh3plume';
check.corr_loctime = 2;               % Correct from UTC to local time
% Define indicators for background or plume concentration scalars
check.variable_background_indicator = 'b';
check.variable_plume_indicator = 'p';
check.figure_font  = 18;


%% Analysis settings

% Basic settings
check.st_xy_xz          = 'xy';                 % type of cross-section ('xy'/'xz')
check.st_var            = {'002', 'nh3_r1p', 'NH_{3, bg}'};
check.st_stats          = {'mean'};               % Statistics, options are ['mean','std','fI','I','S','K','flux'/'F']
check.st_PC             = 'no';                 % Calculate area percentage change ('yes'/'no')
check.st_PCthresholds   = [50; 25; 10; 5];      % Threshold [%] used to calculate the blending-distance
check.st_height         = 37.5;
check.st_twindow        = [14 17];
check.st_tavg           = [10/3600];
% Figure settings
check.st_fI_Cmin        = 0.00001;       % Minimum mean concnetration over which the fluctuation intesity is calculated
check.st_Ithreshhold    = 0.25;        % Threshold value [ppb] for Intermittency factor
check.st_logscale       = 'yes';
check.st_Clim           = [1e-5 0.5];
check.st_xlim           = [400 5200];
check.st_ylim           = [0 4800];
% Transect subplot
check.st_transects      = 'yes';    % Show transects ('yes'/'no')
check.st_transects_xpos = [750 1000 1250 1500 2000 2500 3000 3500 4000];    % x-distances [m] of transects
% manual tick settings
check.st_Cticks_manual  = 'yes';       % Manually set the Cticks of the figures
% Cticks
check.st_Cticks_num     = unique([1e-5 : 1e-5 : 9e-5, 1e-4 : 1e-4 : 9e-4, 1e-3 : 1e-3 : 9e-3, 1e-2 : 1e-2 : 9e-2, 0.1 : 0.1 : 0.5]');
% Ctick labels
check.st_Cticks_str     = unique([1e-5, 5e-5, 1e-4, 5e-4, 1e-3, 5e-3, 1e-2, 5e-2, 0.1, 0.5]');
% Other settings
check.st_skip_loading   = 'no';     % Skip loading (only possible when a figure is already made)
check.st_contourlines   = 'yes';    % Show contour lines ('yes'/'no')
check.st_save_stats     = 'no';     % Save the data in the workpace. This is an options for comparisson with another figure


%% End of settings



% !!!!!!!!!!!!!!!!!!!! Do not change the code below !!!!!!!!!!!!!!!!!!!!



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf('\nPrepare the DALES data\n')      
tload_start = datestr(now, 'dd-mm-yyyy HH:MM:SS');
fprintf([' --> Starting at: ' tload_start '\n']) 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Make sure all runs are covered by "direc.exp_nr"

direc.exp_nr = {}; 
% check.ts_var: 
for i = 1 : size(check.st_var, 1)
    if ismember(check.st_var(i,1), direc.exp_nr) == 0 
        direc.exp_nr(length(direc.exp_nr)+1) = check.st_var(i,1);
    end
end


%% Make sure the directories end with an '\' character

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


%%

if strcmp(check.st_skip_loading, 'yes') == 0
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    fprintf('\nPrepare the data for the statistics cross-section figure\n')      
    tload_start = datestr(now, 'dd-mm-yyyy HH:MM:SS');
    fprintf([' --> Starting at: ' tload_start '\n']) 
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % In case of percentage change, make sure the variable name is not the
    % NH3 plume scalar (which is not used in calculating the PC)
    if strcmp(check.st_PC, 'yes') && strcmp(check.st_stats, 'I') == 0
        if strcmp(check.st_var{2}(end-length(check.variable_plume_indicator)+1:end), ...
                check.variable_plume_indicator)
        
            check.st_var{2} = [ check.st_var{2}(1 : end-length(check.variable_plume_indicator)) ...
                check.variable_background_indicator];
        end
    end
    
    
    % Calculate the area statistics
    [st_stat, st_stat_title, st_stat_label, st_stat_unit, st_C] = ...
        f__calculate_area_stat(check, outp_cross); 
    
    
    
    
    % Calculate Percentage change (PC) if PC switch is turned on
    if strcmp(check.st_PC, 'yes') && strcmp(check.st_stats, 'I') == 0
        
       
       if strcmp(check.st_var{1,2}(1:3),'nh3')
           if strcmp(check.st_var{1,2}(4),'_')
               
               % The background NH3 statistics are already calculated
               st_statBG = st_stat;
               
               % prepare calculating the NH3 total statistics
               st_OR_st_var = check.st_var;
               check.st_var{1,2} = ['nh3' check.st_var{1,2}(5 : end - length(check.variable_background_indicator))];
               
               % Calculate the NH3 total statistics
               [st_statT, ~, ~, ~, ~] = f__calculate_area_stat(check, outp_cross);
               
               check.st_var = st_OR_st_var;
               
           else
               
               % The NH3 total statistics are already calculated
               st_statT = st_stat;
               
               % prepare calculating the background NH3 statistics
               st_OR_st_var = check.st_var;
               check.st_var{1,2} = ['nh3_' check.st_var{1,2}(4 : end) check.variable_background_indicator];
              
               % Calculate the background NH3 statistics
               [st_statBG, ~, ~, ~, ~] = f__calculate_area_stat(check, outp_cross);
                
               check.st_var = st_OR_st_var;
               
           end
           
            st_stat = abs( (st_statT - st_statBG) ./ st_statBG ) .* 100;
           
            st_stat_title = ['PC ' check.st_stats{1}];
            st_stat_label = ['PC_{' st_stat_label '}'];
            st_stat_unit = '[%]';
       end
               
       % Next, calculate the blending distance
       
    end
    
        fprintf(['\n- sigma = ' num2str(nanstd(st_stat(:))) '\n'])   
        fprintf(['\n- mean = ' num2str(nanmean(st_stat(:))) '\n'])   
        fprintf(['\n- median = ' num2str(nanmedian(st_stat(:))) '\n'])   
    
    
    % If switched on, save the st_stat data
    if strcmp(check.st_save_stats, 'yes')
        if exist('st_stat_SAVED') == 0
            st_stat_SAVED = st_stat;
            st_stat_SAVEDname = check.st_var(1,3);
        else
            st_temp = st_stat_SAVED;
            st_stat_SAVED = NaN(size(st_temp,1), size(st_temp,2), size(st_temp,3)+1);
            for i  = 1 : size(st_stat_SAVED,3)
                if i == size(st_stat_SAVED,3)
                    st_stat_SAVED(:,:,i) = st_stat;
                    st_stat_SAVEDname = [st_stat_SAVEDname ; check.st_var(1,3)];
                else
                    st_stat_SAVED(:,:,i) = st_temp(:,:,i);
                end
            end
            clear st_temp
        end
    end
%     st_stat = st_stat_SAVED(:,:,1) - st_stat_SAVED(:,:,2);

    % Preallocate variables for the patch plot
    st_run = ['set' check.st_var{1}];
    st_crossXgrid = NaN(4, size(st_stat,1) * size(st_stat,2));
    st_crossYgrid = NaN(4, size(st_stat,1) * size(st_stat,2));
    st_crossCgrid = NaN(size(st_stat,1) * size(st_stat,2), 1);
    st_crossSgrid = NaN(size(st_stat,1) * size(st_stat,2), 1);
    st_xt = outp_cross.(st_run).xt;
    st_dx = median(st_xt(2:101) - st_xt(1:100));
    st_xm = outp_cross.(st_run).xm - st_dx;
    if strcmp(check.st_xy_xz, 'xy')
        st_yt = outp_cross.(st_run).yt;
        st_dy = st_yt(2) - st_yt(1);
        st_ym = outp_cross.(st_run).ym - st_dy;
    elseif strcmp(check.st_xy_xz, 'xz')
        st_yt = outp_cross.(st_run).zt;
        st_ym = outp_cross.(st_run).zm;
        st_dy = st_yt(2) - st_yt(1);
        st_Sgrid = zeros(size(st_stat));
        st_Sxpos_IN = [];
        for i = 1 : length(outp_cross.(st_run).xt)
            for j = 1 : length(outp_cross.(st_run).yt)
                if outp_cross.(st_run).srcpos.(check.st_var{1,2})(i,j) == 1
                    st_Sxpos_IN = [st_Sxpos_IN; i];
                end
            end
        end
        st_Sgrid(st_Sxpos_IN,1) = 1;
    end
    % Preallocation variables for the contour plot
    [st_crossX, st_crossY] = meshgrid(st_xt , st_yt);
    st_crossC = NaN(size(st_crossX));

    % Quick fix for missing source position for scalar nh3_r0b
    if ismember(check.st_var(1,2), {'w','thl','qt','nh3_r0b'})
        outp_cross.(st_run).srcpos.(check.st_var{1,2}) = ...
            outp_cross.(st_run).srcpos.(outp_cross.(st_run).INFO.Name{13});
    end
    
    for j = 1 : size(st_stat,1)
        for k = 1 : size(st_stat,2)
            st_crossXgrid(:,k + (j-1)*size(st_stat,2)) = ...
                [st_xm(j); st_xm(j)+st_dx; st_xm(j)+st_dx; st_xm(j)]-st_dx;
            st_crossYgrid(:,k + (j-1)*size(st_stat,2)) = ...
                [st_ym(k); st_ym(k); st_ym(k)+st_dy; st_ym(k)+st_dy];

            st_crossC( (st_crossX == st_xt(j) & st_crossY == st_yt(k)) ) = st_stat(j,k);
            st_crossCgrid(k + (j-1)*size(st_stat,2)) = st_stat(j,k);

            if strcmp(check.st_xy_xz, 'xy')
                if outp_cross.(st_run).srcpos.(check.st_var{1,2})(j,k) == 1
                    st_crossSgrid(k + (j-1)*size(st_stat,2)) = ...
                        outp_cross.(st_run).srcpos.(check.st_var{1,2})(j,k);
                end
            elseif strcmp(check.st_xy_xz, 'xz')
                if st_Sgrid(j,k) == 1
                    st_crossSgrid(k + (j-1)*size(st_stat,2)) = st_Sgrid(j,k);
                end
            end
        end
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

end

%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf('\nBuild the cross-section of the statistics figure\n')      
tload_start = datestr(now, 'dd-mm-yyyy HH:MM:SS');
fprintf([' --> Starting at: ' tload_start '\n']) 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf(['      - The min/max in the data = [' ...
    num2str(round(min(st_crossCgrid(:)), 3)) ', ' num2str(round(max(st_crossCgrid(:)), 3)) ']\n'])      
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% Determine the indices of the data to be shown inthe figure
if isnan(check.st_xlim(1)) || isnan(check.st_xlim(2)) 
    st_INx = [1 : length(st_xt)]';
	st_xlim =  [min(st_crossXgrid(:)), max(st_crossXgrid(:))];
else
    st_INx = find(st_xt >= check.st_xlim(1) & st_xt <= check.st_xlim(2));
    st_INx = [st_INx(1) - 1; st_INx; st_INx(end)+1];
    if st_INx(1) == 0
        st_INx(1) = [];
    end
    if st_INx(end) > length(st_xt)
        st_INx(end) = [];
    end
    st_xlim = check.st_xlim;
end
if isnan(check.st_ylim(1)) || isnan(check.st_ylim(2)) 
    st_INy = [1 : length(st_yt)]';
	st_ylim =  [min(st_crossYgrid(:)), max(st_crossYgrid(:))];
else
    st_INy = find(st_yt >= check.st_ylim(1) & st_yt <= check.st_ylim(2));
    st_INy = [st_INy(1) - 1; st_INy; st_INy(end)+1];
    if st_INy(1) == 0
        st_INy(1) = [];
    end
    if st_INy(end) > length(st_yt)
        st_INy(end) = [];
    end
    st_ylim = check.st_ylim;
end




%%%%%
% Prepare the transect data
%%%%%
if strcmp(check.st_transects, 'yes')
    
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    fprintf('\n   - Prepare the data for the transect figure\n')      
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if ismember(check.st_stats, {'I'}) || ...
        strcmp(check.st_var{1,2}(end - length(check.variable_plume_indicator)+1: end), ...
        check.variable_plume_indicator)

        st_var_emit = check.st_var{1,2};

    elseif ismember(check.st_var(1,2), {'w','thl','qt','nh3_r0b'})
        st_varnames = {''};
        for i = length(outp_cross.(st_run).INFO.Name) : -1 : 1
            if length(outp_cross.(st_run).INFO.Name{i}) > 4
                if strcmp(outp_cross.(st_run).INFO.Name{i}(1:4), 'nh3_') && ...
                        strcmp(outp_cross.(st_run).INFO.Name{i}(end - length(check.variable_plume_indicator)+1: end), ...
                        check.variable_plume_indicator)
                    st_varnames = [st_varnames ; outp_cross.(st_run).INFO.Name{i}];
                end
            end
        end
        if length(st_varnames) > 1
            st_var_emit = st_varnames{2};
        else
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    fprintf('\n   !!!!! ERROR !!!!! There is no plume variable\n')      
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            return
        end
    else
        % Define the matfile
        st_var_emit = check.st_var{1,2}(4:end);
        if strcmp(st_var_emit(1), '_')
            st_var_emit(1) = [];
        end
        if strcmp(st_var_emit(end - length(check.variable_background_indicator)+1: end), ...
                check.variable_background_indicator)
            st_var_emit(end - length(check.variable_background_indicator)+1: end) = [];
        end
        st_var_emit = ['nh3_' st_var_emit check.variable_plume_indicator];
    end

    % Source x position
    st_srcpos_x = [0];
    for i = 1 : length(outp_cross.(st_run).yt)
        if sum(outp_cross.(st_run).srcpos.(st_var_emit)(:,i) == 1) > 0
            st_srcpos_x = [st_srcpos_x ; st_xt(...
                outp_cross.(st_run).srcpos.(st_var_emit)(:,i) == 1)];
        end
    end
    st_srcpos_x = max(st_srcpos_x) + (st_dx / 2);
    % Source y position
    st_srcpos_y = [];
    if strcmp(check.st_xy_xz, 'xy')
        for i = 1 : length(st_xt)
            if sum(outp_cross.(st_run).srcpos.(st_var_emit)(i,:) == 1) > 0
                st_srcpos_y = [st_srcpos_y ; st_yt(...
                    outp_cross.(st_run).srcpos.(st_var_emit)(i,:) == 1)];
            end
        end
    elseif strcmp(check.st_xy_xz, 'xz')
        st_srcpos_y = 0;
    end
    st_srcpos_y = mean([min(st_srcpos_y) max(st_srcpos_y)]);

    
    % Will it be a regular plume transect or are we building a PC figure?
    if strcmp(check.st_PC, 'yes') && strcmp(check.st_stats, 'I') == 0
        
        st_angles = [45, 0, -45];
        
        st_anglesX = NaN(length(st_angles), 2);
        st_anglesY = NaN(length(st_angles), 2);
        for j = 1 : length(st_angles)
            
            % For plume centreline 
            if st_angles(j) == 0
                st_anglesX(j,:) = [st_srcpos_x, check.st_xlim(2)];
                st_anglesY(j,:) = [st_srcpos_y, st_srcpos_y];
               
            % For all angles deviating from the plume centreline
            else
                
                temp_x = st_srcpos_x +  ...
                    (check.st_ylim(2) - st_srcpos_y) / tand( abs(st_angles(j)) );
                % For positive angles (towards the north)
                if st_angles(j) > 0 
                    temp_y = st_srcpos_y + ...
                        tand( abs(st_angles(j)) ) * (check.st_xlim(2) - st_srcpos_x);
                
                % For negative angles (towards the south)
                else 
                    temp_y = st_srcpos_y - ...
                        tand( abs(st_angles(j)) ) * (check.st_xlim(2) - st_srcpos_x);
                end         % <-- st_angle(i) > 0 
                
                % In case the line lenght is longer than the xlimit
                if temp_x > check.st_xlim(2)

                    st_anglesX(j,:) = [st_srcpos_x, check.st_xlim(2)];
                    st_anglesY(j,:) = [st_srcpos_y, temp_y];

                else 
                    st_anglesX(j,:) = [st_srcpos_x, temp_x];
                    if st_angles(j) > 0 
                        st_anglesY(j,:) = [st_srcpos_y, check.st_ylim(2)];
                    else
                        st_anglesY(j,:) = [st_srcpos_y, check.st_ylim(1)];
                    end

                end         % <-- temp_x > check.st_xlim(2)
            end         % <-- st_angle(j) == 0
            
        end         % <-- strcmp(check.st_PC, 'yes') && strcmp(check.st_stats, 'I') == 0
        
%         st_dangle
        
        
        
        
    else
        

        % Position the check.st_transects_xpos on the xt grid
        for i = 1 : length(check.st_transects_xpos)
            [~, st_sortIN] = sort(abs(st_xt - check.st_transects_xpos(i)));
            check.st_transects_xpos(i) = st_xt(st_sortIN(1));
        end
        % remove transect data outside xlim
        check.st_transects_xpos(check.st_transects_xpos < st_xlim(1) | ...
            check.st_transects_xpos > st_xlim(2) ) = [];
        % Determine the transect distance from the source
        st_transects_dist = check.st_transects_xpos - st_srcpos_x;

        % Load the emitted NH3 data
        if strcmp(check.st_skip_loading, 'yes') == 0 || exist('st_Cemit') == 0

            if ismember(check.st_stats, {'I'}) || ...
                    strcmp(check.st_var{1,2}(end - length(check.variable_plume_indicator)+1: end), ...
                    check.variable_plume_indicator)
                st_Cemit = st_C;
            else


                if strcmp(check.st_xy_xz, 'xy')
                    [~, st_INheight] = sort(abs(outp_cross.(st_run).zpos - check.st_height));

                    st_matfile_emit = matfile([outp_cross.(st_run).loc ...
                        check.st_var{1,1}  'crossxy' outp_cross.(st_run).cross_lvl{st_INheight(1)} ...
                        '_' st_var_emit '.mat']);

                elseif strcmp(check.st_xy_xz , 'xz')
                    st_matfile_emit = matfile([outp_cross.(st_run).loc ...
                        check.st_var{1,1}  'crossxz_' st_var_emit '.mat']);
                end

                % Load the emitted concentration data
                st_INt_emit = find(outp_cross.(st_run).time > check.st_twindow(1) & ...
                    outp_cross.(st_run).time <= check.st_twindow(2));

                st_t = outp_cross.(st_run).time(st_INt_emit);
                st_dt = median(st_t(5:105)-st_t(4:104));
                st_tnew = st_t;
                if round(check.st_tavg, 8) > round(st_dt, 8)

                    st_percentage = 0.20;
                    st_percentage_base = st_percentage;

                    st_tnew = [st_t(1) - st_dt+check.st_tavg : check.st_tavg : st_t(end)];
                    st_Cemit = NaN(length(st_xt), length(st_yt),length(st_tnew));
                    for i = 1 : length(st_tnew)
                        if i == 1 
                            st_Cemit(:,:,i) = mean(st_matfile_emit.(st_var_emit)(:,:, ...
                                st_INt_emit(st_t <= st_tnew(i))), 3);
                        else
                            st_Cemit(:,:,i) = mean(st_matfile_emit.(st_var_emit)(:,:, ...
                                st_INt_emit(st_t > st_tnew(i-1) & st_t <= st_tnew(i))), 3);
                        end

                        if i/length(st_tnew) >= st_percentage_base
                            fprintf(['      - Loading ' st_var_emit ' data: ' num2str(st_percentage_base * 100) ' %%\n'])
                            st_percentage_base = st_percentage_base + st_percentage;
                        end
                    end

                elseif round(check.st_tavg, 8) == round(st_dt, 8)

                    st_Cemit = st_matfile_emit.(st_var_emit)(:,:, st_INt_emit);

                end

            end
        end

        % Calculate in-plume average
        st_Cemit = mean(st_Cemit, 3);
        st_inplume  = zeros(size(st_Cemit));
        st_inplume(st_Cemit > check.st_fI_Cmin) = 1;
        % Determine plume width
        st_plumewidth = NaN(length(check.st_transects_xpos),1);
        for i = 1 : length(check.st_transects_xpos)
            if isempty( find(st_yt(st_inplume(st_xt == check.st_transects_xpos(i),:) == 1)) ) == 0
                
                st_plumewidth(i) = ...
                    [max(st_yt(st_inplume(st_xt == check.st_transects_xpos(i),:) == 1)) - ...
                    min(st_yt(st_inplume(st_xt == check.st_transects_xpos(i),:) == 1))];
            end
        end

        st_transIN = NaN(length(check.st_transects_xpos), 1);
        for i = 1 : length(st_transIN)
            st_transIN(i) = find(st_xt == check.st_transects_xpos(i));
        end
    end
end



    



%%%%%
% Prepare the color label ticks
%%%%%
if isnan(check.st_Clim(1)) || isnan(check.st_Clim(2))
    st_Cminmax = [min(st_crossCgrid(:)) max(st_crossCgrid(:))];
else
    st_Cminmax = [check.st_Clim(1), check.st_Clim(2)]; 
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    fprintf(['\n   - Area statistics min = ' num2str(min(st_crossCgrid(:))) ...
        ', and max = ' num2str(max(st_crossCgrid(:))) '\n'])      
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if strcmp(check.st_logscale, 'yes')
    if strcmp(check.st_Cticks_manual, 'yes')
        st_temp_TickNum = check.st_Cticks_num;
        st_temp_TickLabelsNum = check.st_Cticks_str;
        temp_TickLabelsStr = cell(length(st_temp_TickNum), 1);
        for i = 1 : length(st_temp_TickNum)
            if ismember(round(st_temp_TickNum(i),8), round(check.st_Cticks_str,8))
                temp_TickLabelsStr{i} = num2str(st_temp_TickNum(i));
            else
                temp_TickLabelsStr{i} = '';
            end
        end
    else
        if st_Cminmax(1) == 0
            st_Cminmax(1) = 0.00001;
        end
        st_ticklabelset = [st_Cminmax(1), 0 : st_Cminmax(1)*5 : st_Cminmax(1)*10];
        st_tickset = [0 : st_Cminmax(1) : st_Cminmax(1)*10];


        st_tickcounter = 0;
        while st_tickcounter == 0
            if st_tickset(end) >= st_Cminmax(2)
                st_tickcounter = 1;
            else
                st_ticklabelset = [st_ticklabelset, 0 : st_tickset(end)*5 : st_tickset(end)*10];
                st_tickset = [st_tickset, 0 : st_tickset(end) : st_tickset(end)*10];
            end
        end
        st_temp_TickNum = unique(st_tickset);
        st_temp_TickLim = {st_temp_TickNum(st_temp_TickNum <= st_Cminmax(1)) ; ...
                           st_temp_TickNum(st_temp_TickNum >= st_Cminmax(2))};
        st_temp_TickLim = [st_temp_TickLim{1}(end), st_temp_TickLim{2}(1)];
        if st_temp_TickLim(1) == 0
            st_temp_TickLim(1) = st_Cminmax(1);
        end
        st_temp_TickNum = st_temp_TickNum(st_temp_TickNum >= st_temp_TickLim(1) & st_temp_TickNum <= st_temp_TickLim(2) );
        if length(st_temp_TickNum) <= 10
            st_ticklevels = NaN(((length(unique(st_tickset))-2)/9),1);
            for j = 1 : length(st_ticklevels)
                st_ticklevels(j) = st_Cminmax(1)* 10^(j-1);
            end
            st_ticklevel = 1;
            st_counter = 0;
            while st_counter == 0
                if st_Cminmax(1) / st_ticklevels(st_ticklevel) >= 1 && ...
                        st_Cminmax(1) / st_ticklevels(st_ticklevel) <= 10
                    st_counter = 1;
                else
                    st_ticklevel = st_ticklevel + 1;
                end
            end
            st_temp_TickNum = [st_ticklevels(st_ticklevel) * floor(st_Cminmax(1) / st_ticklevels(st_ticklevel)) ...
                : st_ticklevels(st_ticklevel)/10 : ...
                st_ticklevels(st_ticklevel) * ceil(st_Cminmax(2) / st_ticklevels(st_ticklevel))];

            st_INmin = find( st_temp_TickNum - min(st_crossCgrid(:)) > 0 );
            st_INmin = st_INmin(1) - 1;
            if st_INmin == 0
                st_INmin = 1;
            end
            st_INmax = find( st_temp_TickNum - max(st_crossCgrid(:)) > 0 );
            st_INmax = st_INmax(1);
            st_temp_TickNum = st_temp_TickNum(st_INmin:st_INmax);
            st_Cminmax = [st_temp_TickNum(1), st_temp_TickNum(end)];

            temp_TickLabelsStr = cell(length(st_temp_TickNum),1);
            temp_INlabelstr = [];
            if length(st_temp_TickNum) <= 11
                for j = 1 : length(st_temp_TickNum)
                    temp_TickLabelsStr{j} = num2str(st_temp_TickNum(j));
                    temp_INlabelstr = [temp_INlabelstr; j];
                end
            elseif length(st_temp_TickNum) > 11 && length(st_temp_TickNum) <= 21
                for j = 1 : 2 : length(st_temp_TickNum)
                    temp_TickLabelsStr{j} = num2str(st_temp_TickNum(j));
                    temp_INlabelstr = [temp_INlabelstr; j];
                end
            else
                for j = 1 : 5 : length(st_temp_TickNum)
                    temp_TickLabelsStr{j} = num2str(st_temp_TickNum(j));
                    temp_INlabelstr = [temp_INlabelstr; j];
                end
            end

            if length(temp_TickLabelsStr) <= 11
                st_temp_TickLabelsNum = st_temp_TickNum(...
                unique([1:1:floor(length(temp_TickLabelsStr)/2),  ...
                floor(length(temp_TickLabelsStr)/2) : 2 : length(st_temp_TickNum)]));
            elseif length(temp_TickLabelsStr) <= 21
                st_temp_TickLabelsNum = st_temp_TickNum(...
                    unique([1:1:temp_INlabelstr(3),  ...
                    temp_INlabelstr(3) : 4 : length(st_temp_TickNum)]));
            else
                st_temp_TickLabelsNum = st_temp_TickNum(...
                    unique([1:1:temp_INlabelstr(2),  ...
                    temp_INlabelstr(2) : 5 : length(st_temp_TickNum)]));
            end

        else
            if length(st_temp_TickNum) <= 15
                temp_TickLabelsStr = cell(length(st_temp_TickNum),1);
                for j = 1 : length(st_temp_TickNum)
                    temp_TickLabelsStr{j} = num2str(st_temp_TickNum(j));
                end
                st_temp_TickLabelsNum = st_temp_TickNum;
            else
                st_temp_TickLabelsNum = unique(st_ticklabelset);
                st_temp_TickLabelsNum = st_temp_TickLabelsNum(ismember(st_temp_TickLabelsNum, st_temp_TickNum));
                temp_TickLabelsStr = cell(length(st_temp_TickNum),1);
                for j = st_temp_TickLabelsNum
                    temp_TickLabelsStr{st_temp_TickNum == j} = num2str(j);
                end
            end
        end
    end
else
    if strcmp(check.st_Cticks_manual, 'yes')
        st_temp_TickNum = check.st_Cticks_num;
        st_temp_TickLabelsNum = check.st_Cticks_str;
        temp_TickLabelsStr = cell(length(st_temp_TickNum), 1);
        for i = 1 : length(st_temp_TickNum)
            if ismember(round(st_temp_TickNum(i),8), round(check.st_Cticks_str,8))
                temp_TickLabelsStr{i} = num2str(st_temp_TickNum(i));
            else
                temp_TickLabelsStr{i} = '';
            end
        end
    else
        [st_temp_TickNum, temp_TickLabelsStr] = f__fig_CTicks(st_Cminmax);
        temp_TickLabelsStr_old = temp_TickLabelsStr;
        temp_TickLabelsStr = cell(length(temp_TickLabelsStr_old),1);
        for i = 1 : length(temp_TickLabelsStr)
            if isempty(temp_TickLabelsStr_old) == 0
                temp_TickLabelsStr{i} = num2str(temp_TickLabelsStr_old{i});
            else
                temp_TickLabelsStr{i} = '';
            end
        end
    end
end


%%
        
st_colormap_Delta = 1000;
st_colormap = [[ [094 : (255 - 094)/(st_colormap_Delta / 2) : 255 ]', ...
                 [060 : (255 - 060)/(st_colormap_Delta / 2) : 255 ]', ...
                 [153 : (255 - 153)/(st_colormap_Delta / 2) : 255 ]'] ...
               ./255 ; ...
               [ [255-(255 - 230)/(st_colormap_Delta / 2) : -(255 - 230)/(st_colormap_Delta / 2) : 230 ]', ...
                 [255-(255 - 097)/(st_colormap_Delta / 2) : -(255 - 097)/(st_colormap_Delta / 2) : 097 ]', ...
                 [255-(255 - 001)/(st_colormap_Delta / 2) : -(255 - 001)/(st_colormap_Delta / 2) : 001 ]'] ./255];

%%%%% Determine the figure & subplot sizes
st_subpos_xsize = 1000;         % Maximum width of the total figure [pixels]
st_subpostran_size = [300 300];    % x and y size of the transect figure [pixels]
st_subpos_Csize = 150;  % Width required for the colorbar [pixels]
% Define the x-y ratio of the cross-section data
st_dom_x = st_xm(st_INx(length(st_INx))) - st_xm(st_INx(1));
st_dom_y = st_ym(st_INy(length(st_INy))) - st_ym(st_INy(1));
st_xy_ratio = [st_dom_y / st_dom_x];

% Definitive position and size of the figure
% Determine the size of the subplot 
st_fig_subx = [ 95 st_subpos_xsize];
st_fig_subx(2) = st_fig_subx(2) - st_fig_subx(1) - st_subpos_Csize;
st_fig_suby = [75 st_fig_subx(2) * st_xy_ratio];
if strcmp(check.st_transects, 'yes')
    % Determine the size of the subplot 
    st_fig_subx(2) = st_fig_subx(2) - st_subpostran_size(1) - 50;
    st_fig_suby(2) = st_fig_subx(2) * st_xy_ratio;
    if st_fig_suby(2) < st_subpostran_size(2)
        st_fig_suby(2) = st_subpostran_size(2);
    end
end
st_fig_subsize = [st_fig_subx(1), st_fig_suby(1), ...
    st_fig_subx(2), st_fig_suby(2)]; 
st_fig_position = [20, 40, st_subpos_xsize, ...
    st_fig_subsize(4)+st_fig_subsize(2)+20];
% Definitive position of the subplot in percentages (betweem 0 and 1)
st_fig_subposition = [st_fig_subsize(1)/st_fig_position(3), ... 
                  st_fig_subsize(2)/st_fig_position(4), ...
                  st_fig_subsize(3)/st_fig_position(3), ...
                  st_fig_subsize(4)/st_fig_position(4)];
if strcmp(check.st_transects, 'yes')
    st_fig_transize = [st_fig_position(3) - st_subpostran_size(2) - 10, ... % 
        (st_fig_position(4) - st_subpostran_size(2) - 80)/2+80, ...
        st_subpostran_size(1), st_subpostran_size(2)];
    
    st_fig_tranposition = [st_fig_transize(1)/st_fig_position(3), ... 
                      st_fig_transize(2)/st_fig_position(4), ...
                      st_fig_transize(3)/st_fig_position(3), ...
                      st_fig_transize(4)/st_fig_position(4)];
end

%%

if strcmp(check.st_PC, 'yes')
%     st_fig_position_OR = [20 40 1000 500];
%     st_fig_position = st_fig_position_OR;
    st_move_transectplot = 50;
    st_fig_position = st_fig_position+[0 0 st_move_transectplot 0];
    st_fig = figure('units','pixels', 'Color','w',...
        'innerposition', st_fig_position);
    % Build the plotting space
    st_fig_subposition = [95/st_fig_position(3), st_fig_subposition(2), ...
        405/st_fig_position(3), st_fig_subposition(4)];
    st_subcross = subplot('Position', st_fig_subposition);
else
    st_fig = figure('units','pixels', 'Color','w',...
        'innerposition', st_fig_position);
    % Build the plotting space
    st_subcross = subplot('Position',st_fig_subposition);
end
hold on; 
% Plot cross section
st_gridplot = patch(st_crossXgrid, st_crossYgrid, st_crossCgrid);
st_gridplot.LineWidth = 0.01;
st_gridplot.EdgeAlpha = 0.2;
st_gridplot.EdgeColor = 'none';

st_gca = gca;
st_ax = axis;
if strcmp(check.st_contourlines, 'yes')
    % Contour lines
    st_contour_num = [];
    for i = 1 : length(temp_TickLabelsStr)
        if isempty(temp_TickLabelsStr{i}) == 0
            st_contour_num = [st_contour_num; str2num(temp_TickLabelsStr{i})];
        end
    end
%     st_contour_num = [1e-5 1e-4 1e-3 1e-2 1e-1 1 10 100];
    st_contour = contour(st_crossX, st_crossY, st_crossC, ...
        st_contour_num, 'LineColor','k', 'LineWidth',1.2, ...
        'ShowText','off', 'TextListMode','manual', 'TextList', st_contour_num);
end
if strcmp(check.st_PC,'yes')
        
    st_trans_cm_Delta = length(check.st_PCthresholds);
    if mod(st_trans_cm_Delta/2, 1) == 0
        st_trans_cm = [[ [094 : (255 - 094)/(st_trans_cm_Delta / 2) : 255 ]', ...
                         [060 : (255 - 060)/(st_trans_cm_Delta / 2) : 255 ]', ...
                         [153 : (255 - 153)/(st_trans_cm_Delta / 2) : 255 ]'] ...
                       ./255];
        st_trans_cm(end,:) = [];
        st_trans_cm = [st_trans_cm ; ...
                       [ [255-(255 - 230)/(st_trans_cm_Delta / 2) : -(255 - 230)/(st_trans_cm_Delta / 2) : 230 ]', ...
                         [255-(255 - 097)/(st_trans_cm_Delta / 2) : -(255 - 097)/(st_trans_cm_Delta / 2) : 097 ]', ...
                         [255-(255 - 001)/(st_trans_cm_Delta / 2) : -(255 - 001)/(st_trans_cm_Delta / 2) : 001 ]'] ./255];
    else
        st_trans_cm = [[ [094 : (255 - 094)/(st_trans_cm_Delta / 2) : 255 ]', ...
                     [060 : (255 - 060)/(st_trans_cm_Delta / 2) : 255 ]', ...
                     [153 : (255 - 153)/(st_trans_cm_Delta / 2) : 255 ]'] ...
                   ./255 ; ...
                   [ [255-(255 - 230)/(st_trans_cm_Delta / 2) : -(255 - 230)/(st_trans_cm_Delta / 2) : 230 ]', ...
                     [255-(255 - 097)/(st_trans_cm_Delta / 2) : -(255 - 097)/(st_trans_cm_Delta / 2) : 097 ]', ...
                     [255-(255 - 001)/(st_trans_cm_Delta / 2) : -(255 - 001)/(st_trans_cm_Delta / 2) : 001 ]'] ./255];
    end
    
    st_contour_pos = cell(length(check.st_PCthresholds), 1);
    for j = 1 : length(check.st_PCthresholds)
        
        
        [~, st_contour] = contour(st_crossX, st_crossY, st_crossC, ...
            [check.st_PCthresholds(j), check.st_PCthresholds(j)], 'ShowText','off', ...
            'LineColor',st_trans_cm(length(check.st_PCthresholds) -(j-1),:), 'LineWidth',3);
        
        st_INkeep = find(st_contour.ContourMatrix(1,:) > st_srcpos_x - 100);
        st_contour_pos{j} = st_contour.ContourMatrix(:,st_INkeep);
    end
end
% Plot the location of the sources
st_sourceplot = patch(st_crossXgrid(:,st_crossSgrid == 1), st_crossYgrid(:,st_crossSgrid == 1), 'r');
st_sourceplot.EdgeColor = 'none';
st_sourceplot.FaceAlpha = 1.00;
% plot the transect lines
if strcmp(check.st_transects, 'yes')
    if strcmp(check.st_PC,'yes')
        
        for i = 1 : size(st_anglesX, 1)
            
            plot(st_anglesX(i,:), st_anglesY(i,:), ...
                '--k', 'LineWidth', 1.5)
        end
        
    else
        
        plot(st_xt(st_xt > st_srcpos_x & st_xt <= max(check.st_transects_xpos)), ...
            repmat(st_srcpos_y, [length(st_xt(st_xt > st_srcpos_x & st_xt <= max(check.st_transects_xpos))),1]), ...
            '--k', 'LineWidth', 1.0)
        for i = 1 : length(check.st_transects_xpos)
            plot([check.st_transects_xpos(i) check.st_transects_xpos(i)], ...
                [min(st_yt(st_inplume(st_transIN(i),:) == 1)) max(st_yt(st_inplume(st_transIN(i),:) == 1))],...
                ':k', 'LineWidth', 1.5)
            if mod(i,2) == 0
                text(check.st_transects_xpos(i), min(st_yt(st_inplume(st_transIN(i),:) == 1)), ...
                    ['d_{x,' num2str(i) '}'], 'FontSize',check.figure_font, ...
                    'HorizontalAlignment','center', 'VerticalAlignment','top')
            else
                text(check.st_transects_xpos(i), max(st_yt(st_inplume(st_transIN(i),:) == 1)), ...
                    ['d_{x,' num2str(i) '}'], 'FontSize',check.figure_font, ...
                    'HorizontalAlignment','center', 'VerticalAlignment','bottom')
            end
        end
%       st_yt(st_inplume(st_transIN(i),:) == 1)
    end
    
    
end
% for i = 2 : 24
%     plot(check.st_xlim, [outp_cross.set012.zt(i) outp_cross.set012.zt(i)], ':k')
% end
% Settings
st_c = colorbar(st_subcross,'Ticks',st_temp_TickNum, 'Ticklabels', temp_TickLabelsStr);
% Colorbar Settings
if strcmp(check.st_logscale, 'yes')
    st_subcross.ColorScale = 'log';
    st_gca.ColorScale = 'log';
    st_gca.FontSize = check.figure_font;
else
    st_gca.FontSize = check.figure_font;
end
st_c.Position = [(st_fig_subsize(1) + st_fig_subsize(3) + 10) / st_fig_position(3), ...
                 (st_fig_subsize(2)) / st_fig_position(4), ...
                 (30) / st_fig_position(3), ...
                 (st_fig_subsize(4)) / st_fig_position(4)];
if strcmp(check.st_transects, 'yes') == 0 || strcmp(check.st_PC,'yes')
    st_c.Label.String = [st_stat_label ' ' st_stat_unit];
end
% % %     st_c.Label.String = ['NH_3 ' st_stat_unit];         % UNCOMMENT FOR MEAN PLOT
st_c.FontSize = check.figure_font;
st_c.LineWidth = 1.2;
st_c.TickLength = 0.03;
st_c.Ticks = st_temp_TickNum;
st_c.TickLabels = temp_TickLabelsStr;
st_c.Limits = [st_temp_TickNum(1) st_temp_TickNum(end)];   % Limits of the colorbar
if strcmp(check.st_PC,'yes')
    colormap(flipud(bone)); %
else
	colormap(st_colormap); %
end
st_subcross.CLim = [st_temp_TickNum(1) st_temp_TickNum(end)];   % Limits of the colors in the figure
xlabel([check.st_xy_xz(1) ' [m]']);
ylabel([check.st_xy_xz(2) ' [m]']); 
xlim(st_xlim)
ylim(st_ylim)
% Make a new subplot to add highlights to the colorbar
st_subbar = subplot('Position',st_c.Position);
hold on
for i = 1 : length(check.st_Cticks_str)
    if ismember(check.st_Cticks_str(i), check.st_Clim) == 0
        plot([0.3 1], [check.st_Cticks_str(i) check.st_Cticks_str(i)] ,'-k','LineWidth',1.7)
    end
end
set(st_subbar, 'XTick',[],'XLim',[0 1],'YLim',[check.st_Clim],...
    'YTick', [], 'XColor','none', 'YColor','none', 'Visible','off')
if strcmp(check.st_logscale, 'yes')
    set(st_subbar,'YScale','log')
end

%%%%%
% Transect subplot
%%%%%
if strcmp(check.st_transects, 'yes')
    
    % Wind direction line plot
    if strcmp(check.st_PC,'yes')
        
        st_trans_cm_Delta = length(check.st_PCthresholds);
        if mod(st_trans_cm_Delta/2, 1) == 0
            st_trans_cm = [[ [094 : (255 - 094)/(st_trans_cm_Delta / 2) : 255 ]', ...
                             [060 : (255 - 060)/(st_trans_cm_Delta / 2) : 255 ]', ...
                             [153 : (255 - 153)/(st_trans_cm_Delta / 2) : 255 ]'] ...
                           ./255];
            st_trans_cm(end,:) = [];
            st_trans_cm = [st_trans_cm ; ...
                           [ [255-(255 - 230)/(st_trans_cm_Delta / 2) : -(255 - 230)/(st_trans_cm_Delta / 2) : 230 ]', ...
                             [255-(255 - 097)/(st_trans_cm_Delta / 2) : -(255 - 097)/(st_trans_cm_Delta / 2) : 097 ]', ...
                             [255-(255 - 001)/(st_trans_cm_Delta / 2) : -(255 - 001)/(st_trans_cm_Delta / 2) : 001 ]'] ./255];
        else
            st_trans_cm = [[ [094 : (255 - 094)/(st_trans_cm_Delta / 2) : 255 ]', ...
                         [060 : (255 - 060)/(st_trans_cm_Delta / 2) : 255 ]', ...
                         [153 : (255 - 153)/(st_trans_cm_Delta / 2) : 255 ]'] ...
                       ./255 ; ...
                       [ [255-(255 - 230)/(st_trans_cm_Delta / 2) : -(255 - 230)/(st_trans_cm_Delta / 2) : 230 ]', ...
                         [255-(255 - 097)/(st_trans_cm_Delta / 2) : -(255 - 097)/(st_trans_cm_Delta / 2) : 097 ]', ...
                         [255-(255 - 001)/(st_trans_cm_Delta / 2) : -(255 - 001)/(st_trans_cm_Delta / 2) : 001 ]'] ./255];
        end
        
        
        st_subtran = subplot('Position',[(690+st_move_transectplot)/st_fig_position(3), ...
            st_fig_tranposition(2), 300/st_fig_position(3), st_fig_tranposition(4)]);
        hold on; grid on
        st_PClinestyle = {'-','--',':','-.'};
        st_lgndplot = gobjects(length(check.st_PCthresholds),1);
        st_contour_dist_max = 0;
        plot([-45 -45], [0 999999], 'Color','k', 'LineStyle', '--', 'LineWidth', 0.8)
        plot([0 0], [0 999999], 'Color','k', 'LineStyle', '--', 'LineWidth', 0.8)
        plot([45 45], [0 999999], 'Color','k', 'LineStyle', '--', 'LineWidth', 0.8)
        for i = 1 : length(check.st_PCthresholds)
            
            st_contour_angle_temp = (-1).* atand( (st_contour_pos{i}(2,:)' - st_srcpos_y) ./ ...
                (st_contour_pos{i}(1,:)' - st_srcpos_x) );
            st_contour_dist_temp = sqrt( ( st_contour_pos{i}(1,:)' - st_srcpos_x ).^2 + ...
                (st_contour_pos{i}(2,:)' - st_srcpos_y) .^2 );
            
            [st_contour_angle_temp, st_INsort] = sort(st_contour_angle_temp);
            st_contour_dist_temp = st_contour_dist_temp(st_INsort);
            st_contour_dist_max = max([st_contour_dist_max, max(st_contour_dist_temp)]);
            st_INmin = find(st_contour_dist_temp < 50);
            st_contour_angle_temp(st_INmin) = [];
%             st_contour_angle_temp = round(st_contour_angle_temp,2);
            st_contour_dist_temp(st_INmin) = [];
            
            st_contour_angle = unique(st_contour_angle_temp);
            st_contour_dist = NaN(size(st_contour_angle));
            for j = 1 : length(st_contour_angle)
                st_contour_dist(j) = max( ...
                    st_contour_dist_temp(st_contour_angle_temp == st_contour_angle(j)) );
            end
            
            
%             st_lgndplot(i) = plot( st_contour_angle, st_contour_dist./1000, ...
%                 'LineWidth',3, 'Color',st_trans_cm(length(check.st_PCthresholds)-i+1,:), ...
%                 'LineStyle', st_PClinestyle{1});
            st_lgndplot(i) = scatter( st_contour_angle, st_contour_dist./1000, ...
                'MarkerEdgeColor',st_trans_cm(length(check.st_PCthresholds)-i+1,:), ...
                'Marker', '.', 'SizeData', 150);
        end
        xlabel('Wind direction')
        ylabel(['BD_{' check.st_stats{1} '} [km]'])    
        xlim([-65 65])
        ylim([0, 100 * ceil(st_contour_dist_max/100)]./1000)
        st_xtick = sort(unique([-90 : 10 : 90, -45, 45]));
        st_xticklabel = repmat({''},[length(st_xtick), 1]);
        st_INxtick = find(ismember(st_xtick, [-45, 0, 45]));
        st_xticklabel(st_INxtick) = [{'NW'},{'W'},{'SW'}];
        st_ytick = sort(unique([0 : 0.25 : 10]));
        st_yticklabel = repmat({''},[length(st_ytick), 1]);
        st_yticklabel_num = sort(unique([0 : 0.5 : 10]));
        st_INytick = find(ismember(st_ytick, st_yticklabel_num ));
        for j = 1 : length(st_INytick)
            st_yticklabel{st_INytick(j)} = num2str(st_yticklabel_num(j));
        end
        set(gca, 'FontSize',check.figure_font, 'Xtick',st_xtick, 'XTickLabel',st_xticklabel, ...
            'Ytick',st_ytick, 'YTickLabel',st_yticklabel)
        
        %%%%%
        % Make the legend figure
        %%%%%
        st_lgnd_columns = 4;
        st_lgnd_rows = ceil(length(check.st_PCthresholds) / st_lgnd_columns);
        st_lgndstr = cell(length(check.st_PCthresholds), 1);
        for j = 1 : length(check.st_PCthresholds)
            if j == length(check.st_PCthresholds)
                st_lgndstr{j} = ['_{  }' num2str(check.st_PCthresholds(j)) ...
                    '% threshold'];
            else
                st_lgndstr{j} = ['_{  }' num2str(check.st_PCthresholds(j)) ...
                    '% threshold_{       }'];
            end
        end

        figure('units','pixels', 'Color','w',...
            'innerposition', [10 600 1000 36*st_lgnd_rows], ...
            'Name', 'Receptor stats');
        subplot('Position',[-1 -1 0.8 0.5])
        hold on
        st_obj = gobjects(length(check.st_PCthresholds),1);
        for i = 1 : length(check.st_PCthresholds)
            st_obj(i) = plot( [0 1], [i i], ...
                'LineWidth',3, 'Color',st_trans_cm(length(check.st_PCthresholds)-i+1,:));
        end
        set(gca, 'FontSize', check.figure_font)
        st_lgnd = legend(st_obj,st_lgndstr,'box','off', 'Orientation','vertical');
        st_lgnd.NumColumns = st_lgnd_columns;
        st_lgnd.Position = [0.022 0.06 0.98 0.91];
        
        
        
        
    % Regular transect plot
    else
        
        st_trans_cm_Delta = length(st_transIN);
        if mod(st_trans_cm_Delta/2, 1) == 0
            st_trans_cm = [[ [094 : (255 - 094)/(st_trans_cm_Delta / 2) : 255 ]', ...
                             [060 : (255 - 060)/(st_trans_cm_Delta / 2) : 255 ]', ...
                             [153 : (255 - 153)/(st_trans_cm_Delta / 2) : 255 ]'] ...
                           ./255];
            st_trans_cm(end,:) = [];
            st_trans_cm = [st_trans_cm ; ...
                           [ [255-(255 - 230)/(st_trans_cm_Delta / 2) : -(255 - 230)/(st_trans_cm_Delta / 2) : 230 ]', ...
                             [255-(255 - 097)/(st_trans_cm_Delta / 2) : -(255 - 097)/(st_trans_cm_Delta / 2) : 097 ]', ...
                             [255-(255 - 001)/(st_trans_cm_Delta / 2) : -(255 - 001)/(st_trans_cm_Delta / 2) : 001 ]'] ./255];
        else
            st_trans_cm = [[ [094 : (255 - 094)/(st_trans_cm_Delta / 2) : 255 ]', ...
                         [060 : (255 - 060)/(st_trans_cm_Delta / 2) : 255 ]', ...
                         [153 : (255 - 153)/(st_trans_cm_Delta / 2) : 255 ]'] ...
                       ./255 ; ...
                       [ [255-(255 - 230)/(st_trans_cm_Delta / 2) : -(255 - 230)/(st_trans_cm_Delta / 2) : 230 ]', ...
                         [255-(255 - 097)/(st_trans_cm_Delta / 2) : -(255 - 097)/(st_trans_cm_Delta / 2) : 097 ]', ...
                         [255-(255 - 001)/(st_trans_cm_Delta / 2) : -(255 - 001)/(st_trans_cm_Delta / 2) : 001 ]'] ./255];
        end

        st_lgndplot = gobjects(length(st_transIN),1);
        st_lgndstr = cell(length(st_transIN),1);


        st_subtran = subplot('Position',st_fig_tranposition);
        hold on; grid on
        st_trans_minmax = NaN(length(st_transIN),2);
        for i = 1 : length(st_transIN)
            st_lgndplot(i) = plot( (st_yt(st_inplume(st_transIN(i),:) == 1) - st_srcpos_y )./st_plumewidth(i), ...
                st_stat(st_transIN(i) ,st_inplume(st_transIN(i),:) == 1), ...
                'LineWidth',3, 'Color',st_trans_cm(length(st_transIN)-i+1,:));
            st_trans_minmax(i,1) = min((st_yt(st_inplume(st_transIN(i),:) == 1) - st_srcpos_y )./st_plumewidth(i));
            st_trans_minmax(i,2) = max((st_yt(st_inplume(st_transIN(i),:) == 1) - st_srcpos_y )./st_plumewidth(i));
            st_lgndstr{i} = ['d_{x,' num2str(i) '} = ' num2str(st_transects_dist(i)) ' m   _{ }'];
        end
        st_trans_minmax = [min(st_trans_minmax(:,1)) max(st_trans_minmax(:,2))];
        xlabel('y/y_{plume} [-]')
        ylabel([st_stat_label ' ' st_stat_unit])      % COMMENT FOR MEAN PLOT
        set(gca, 'FontSize',check.figure_font, 'Xtick',[-1 : 0.1 : 1], ...
            'XTickLabel',{'-1','','-0.8','','-0.6','','-0.4','','-0.2','','0','','0.2','','0.4','','0.6','','0.8','','1'})
        if st_temp_TickNum(1) < 0
            plot([-9999 9999], [0 0], '-k','LineWidth',0.5)
        end
        if strcmp(check.st_xy_xz, 'xy')
            plot([0 0], [-99999 99999],'-k', 'LineWidth',0.5)
            xlim([-1*max(abs(st_trans_minmax)), max(abs(st_trans_minmax))])
        elseif strcmp(check.st_xy_xz, 'xz')
            xlim([0 st_trans_minmax(2)])
        end
        if strcmp(check.st_logscale, 'yes')
            st_INylim = find(st_temp_TickNum >= 0.01 * 1*10^ceil( log10(st_temp_TickNum(end) )));
            set(gca, 'YTick',st_temp_TickNum(st_INylim), 'YTickLabel',temp_TickLabelsStr(st_INylim))
        else
            set(gca, 'YTick',st_temp_TickNum, 'YTickLabel',temp_TickLabelsStr)
        end
        ylim([st_temp_TickNum(1) st_temp_TickNum(end)])

        st_title = [check.st_var{3} ' ' st_stat_title];
        annotation('textbox', [st_fig_tranposition(1) - 0.1, 0.9, ...
            st_fig_tranposition(3)+0.1, 0.1], 'string', st_title, ...
                    'EdgeColor','none', 'FontSize',check.figure_font, 'FontWeight','normal', ...
                    'HorizontalAlignment','center','VerticalAlignment','bottom', ...
                    'FontWeight','Bold');

    % % %     % UNCOMMENT FOR MEAN PLOT
    % % %     st_mean_anno = annotation('line','Position',[0.787 0.398, 0.088 0.001],'LineWidth',1.2);
    % % %     st_mean_anno.Position = [0.63 0.449, 0.001 0.070];
    % % %     % UNCOMMENT FOR MEAN PLOT

        %%%%%
        % Make the legend figure
        %%%%%
        if length(check.st_transects_xpos) > 10
            st_lgnd_columns = 5;
            st_lgnd_rows = ceil(length(check.st_transects_xpos) / st_lgnd_columns);
        else
            st_lgnd_rows = 2;
            st_lgnd_rows = 10;
            st_lgnd_columns = ceil(length(check.st_transects_xpos)/st_lgnd_rows);
        end

        figure('units','pixels', 'Color','w',...
            'innerposition', [10 600 1000 36*st_lgnd_rows], ...
            'Name', 'Receptor stats');
        subplot('Position',[-1 -1 0.8 0.5])
        hold on
        st_obj = gobjects(length(check.st_transects_xpos),1);
        for i = 1 : length(check.st_transects_xpos)
            st_obj(i) = plot( [0 1], [i i], ...
                'LineWidth',3, 'Color',st_trans_cm(length(check.st_transects_xpos)-i+1,:));
        end
        set(gca, 'FontSize', check.figure_font)
        st_lgnd = legend(st_obj,st_lgndstr,'box','off', ...
            'NumColumns',st_lgnd_columns);
        st_lgnd.Position = [0.0001 0.06 0.9998 0.91];
        st_lgnd.Position = [0.0001 0.06 0.25 0.91];


    %     st_lgndstr_2 = st_lgndstr;
    %     for i = 1 : length(st_lgndstr_2)
    %         st_lgndstr_2{i} = [st_lgndstr_2{i} '_{     }'];
    %     end
    %     
    %     figure('units','pixels', 'Color','w',...
    %         'innerposition', [10 600 1000*2 36*st_lgnd_rows], ...
    %         'Name', 'Receptor stats');
    %     subplot('Position',[-1 -1 0.8 0.5])
    %     hold on
    %     st_obj = gobjects(length(check.st_transects_xpos),1);
    %     for i = 1 : length(check.st_transects_xpos)
    %         st_obj(i) = plot( [0 1], [i i], ...
    %             'LineWidth',3, 'Color',st_trans_cm(length(check.st_transects_xpos)-i+1,:));
    %     end
    %     set(gca, 'FontSize', check.figure_font)
    %     st_lgnd = legend(st_obj,st_lgndstr_2,'box','off', ...
    %         'NumColumns',st_lgnd_columns);
    %     st_lgnd.Position = [0.0001 0.06 0.9998 0.91];
    
    end
    
end
% temp_hours = {num2str(floor(check.st_twindow(1))), num2str(floor(check.st_twindow(2)))};
% temp_minutes = {num2str((check.st_twindow(1)- floor(check.st_twindow(1)))*60), ...
%     num2str((check.st_twindow(2)- floor(check.st_twindow(2)))*60)};
% for i = 1: length(temp_minutes)
%     if length(temp_minutes{i}) == 1
%         temp_minutes{i} = ['0' temp_minutes{i}];
%     end
% end

% figure 
% hold on
% for i = 1 : length(st_contour_pos)
%     scatter(st_contour_pos{i}(1,:), st_contour_pos{i}(2,:))
% end

%%


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





function [c_tick, c_tick_label] = f__fig_CTicks(c_lim)

    % CLim
    dc_lim = (c_lim(2) - c_lim(1));

    II = 0;
    ii = 0.5;
    iter_limit = 0;
    dc_check = 0.00001;
    while II == 0
        if iter_limit > 100
            II = 1;
        elseif dc_lim > dc_check && ii == 0
            dc_check = dc_check *2;
            ii = mod(ii + mod(0.5,1),1);
            iter_limit = iter_limit+1;
        elseif dc_lim > dc_check && ii == 0.5
            dc_check = dc_check *5;
            ii = mod(ii + mod(0.5,1),1);
            iter_limit = iter_limit+1;
        else
            II = 1;
        end
    end

    if ii == 0
        dy = dc_check / 25;
    elseif ii == 0.5
        dy = dc_check / 10;
    end
    c_tick = [c_lim(1):dy:c_lim(2)];
    
    II = 0;
    iter_limit = 0;
    while II == 0
        if iter_limit > 100
            II = 1;
        elseif length(c_tick) < 20
            dy = dy / 2;
            c_tick = [c_lim(1):dy:c_lim(2)];
            iter_limit = iter_limit+1;
        elseif length(c_tick) > 40
            dy = dy * 2;
            c_tick = [c_lim(1):dy:c_lim(2)];
            iter_limit = iter_limit+1;
        else
            II = 1;
        end
    end

    
    IN_label = [1:4:length(c_tick)];
    c_tick_label = cell(length(c_tick),1);
    for i = 1 : length(IN_label)
        c_tick_label{IN_label(i)} = round(c_tick(IN_label(i)),8);
    end
    
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%             Calculate the area statistics           %%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [st_stat, st_stat_title, st_stat_label, st_stat_unit, st_C] = ...
    f__calculate_area_stat(check, outp_cross)
    
    
    
    st_run = ['set' check.st_var{1}];

    st_INt = find(outp_cross.(st_run).time > check.st_twindow(1)-1 & outp_cross.(st_run).time <= check.st_twindow(2));
    st_t = outp_cross.(st_run).time(st_INt);
    
    

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    fprintf('   - Load the data from the .mat files\n')      
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


    if ismember(check.st_stats, {'I'}) && ...
            strcmp(check.st_var{1,2}(end - length(check.variable_plume_indicator)+1: end), ...
            check.variable_plume_indicator) == 0

        check.st_var{1,2} = check.st_var{1,2}(4:end);
        check.st_var{1,2} = ['nh3_' check.st_var{1,2} check.variable_plume_indicator];

    end

    st_tnew = st_t;
    st_dt = median(st_t(5:105)-st_t(4:104));
    if strcmp(check.st_xy_xz, 'xy')
        [~, st_INheight] = sort(abs(outp_cross.(st_run).zpos - check.st_height));
        st_height = num2str(round(outp_cross.(st_run).zpos(st_INheight(1)),2));
        st_matfile = matfile([outp_cross.(st_run).loc ...
            check.st_var{1,1}  'crossxy' outp_cross.(st_run).cross_lvl{st_INheight(1)} '_' check.st_var{1,2} '.mat']);
        
    elseif strcmp(check.st_xy_xz , 'xz')
        st_matfile = matfile([outp_cross.(st_run).loc ...
            check.st_var{1,1}  'crossxz_' check.st_var{1,2} '.mat']);
    end
        
        
    if round(check.st_tavg, 8) > round(st_dt, 8)

        st_tnew = [st_t(1) - st_dt+check.st_tavg : check.st_tavg : st_t(end)];
        if strcmp(check.st_xy_xz, 'xy')
            st_C = NaN(length(outp_cross.(st_run).xt), length(outp_cross.(st_run).yt), length(st_INt)/(check.st_tavg*3600));
        elseif strcmp(check.st_xy_xz , 'xz')
            st_C = NaN(length(outp_cross.(st_run).xt), length(outp_cross.(st_run).zt), length(st_INt)/(check.st_tavg*3600));
        end
        
        st_percentage = 0.20;
        st_percentage_base = st_percentage;

        for i = 1 : length(st_tnew)
            if i == 1 
                st_C(:,:,i) = mean(st_matfile.(check.st_var{1,2})(:,:, ...
                    st_INt(st_t <= st_tnew(i))), 3);
            else
                st_C(:,:,i) = mean(st_matfile.(check.st_var{1,2})(:,:, ...
                    st_INt(st_t > st_tnew(i-1) & st_t <= st_tnew(i))), 3);
            end

            if i/length(st_tnew) >= st_percentage_base
                fprintf(['      - Loading ' check.st_var{1} '-' check.st_var{2} ...
                    ' data : ' num2str(st_percentage_base * 100) ' %%\n'])
                st_percentage_base = st_percentage_base + st_percentage;
            end
        end

    elseif round(check.st_tavg, 8) == round(st_dt, 8)

        st_C = st_matfile.(check.st_var{1,2})(:,:, st_INt);

    end
    
    st_dt_new = median(st_tnew(5:105)-st_tnew(4:104));

    if ismember(check.st_stats, {'std','fI','S','K'})
        
        st_Cma = zeros(size(st_C));
    	if ismember(check.st_var(1,2), {'w','thl','qt'}) == 0
            st_Cma = movmean(st_C, [1 / (st_dt_new), 0],3);
        end
        
    elseif ismember(check.st_stats, {'flux', 'F'})
        
        
        if strcmp(check.st_xy_xz, 'xy')
            st_w_matfile = matfile([outp_cross.(st_run).loc ...
                check.st_var{1,1}  'crossxy' outp_cross.(st_run).cross_lvl{st_INheight(1)} '_w.mat']);

        elseif strcmp(check.st_xy_xz , 'xz')
            st_w_matfile = matfile([outp_cross.(st_run).loc ...
                check.st_var{1,1}  'crossxz_w.mat']);
        end
        

        st_tnew = st_t;
        if round(check.st_tavg, 8) > round(st_dt, 8)

            st_percentage = 0.20;
            st_percentage_base = st_percentage;

            st_tnew = [st_t(1) - st_dt+check.st_tavg : check.st_tavg : st_t(end)];
            if strcmp(check.st_xy_xz, 'xy')
                st_w = NaN(length(outp_cross.(st_run).xt), length(outp_cross.(st_run).yt), length(st_INt)/(check.st_tavg*3600));
            elseif strcmp(check.st_xy_xz , 'xz')
                st_w = NaN(length(outp_cross.(st_run).xt), length(outp_cross.(st_run).zt), length(st_INt)/(check.st_tavg*3600));
            end
            for i = 1 : length(st_tnew)
                if i == 1 
                    st_w(:,:,i) = mean(st_w_matfile.w(:,:, ...
                        st_INt(st_t <= st_tnew(i))), 3);
                else
                    st_w(:,:,i) = mean(st_w_matfile.w(:,:, ...
                        st_INt(st_t > st_tnew(i-1) & st_t <= st_tnew(i))), 3);
                end

                if i/length(st_tnew) >= st_percentage_base
                    fprintf(['      - Loading w data: ' num2str(st_percentage_base * 100) ' %%\n'])
                    st_percentage_base = st_percentage_base + st_percentage;
                end
            end

        elseif round(check.st_tavg, 8) == round(st_dt, 8)

            st_w = st_w_matfile.w(:,:, st_INt);

        end
        
    end
    st_INt_correct = find( st_tnew > check.st_twindow(1) & st_tnew <= check.st_twindow(2) );
    

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    fprintf('   - Processing the cross-section data\n')      
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    st_tnew = st_tnew(st_INt_correct);
    st_C = st_C(:,:,st_INt_correct);
    if ismember(check.st_stats, {'std','fI','S','K'})
        st_Cma = st_Cma(:,:,st_INt_correct);
    elseif ismember(check.st_stats, {'flux','F'})
        st_w = st_w(:,:,st_INt_correct);
    end
    


    if ismember(check.st_stats, {'mean'})

        st_stat = mean(st_C, 3);
        st_stat_title = 'mean';
        st_stat_label = '[NH_3]';
        st_stat_unit = '[ppb]';

    elseif ismember(check.st_stats, {'std'})

        st_stat = std( st_C-st_Cma , 0, 3);
        st_stat_title = 'standard deviation';
        st_stat_label = '\sigma';
        st_stat_unit = '[ppb]';

    elseif ismember(check.st_stats, {'fI'})

        st_temp = std( st_C-st_Cma ,0,3) ./ mean( st_C ,3);
        st_stat = NaN( size(st_C,1), size(st_C,2) );
        st_stat( mean(st_C, 3) > check.st_fI_Cmin) = st_temp( mean(st_C, 3) > check.st_fI_Cmin);
        clear st_temp
        st_stat_title = 'fluctuation intensity';
        st_stat_label = 'fI';
        st_stat_unit = '[-]';

    elseif ismember(check.st_stats, {'I'})

        st_stat = sum( st_C >= check.st_Ithreshhold , 3) ./ size(st_C,3);
        st_stat_title = 'intermittency';
        st_stat_label = 'I';
        st_stat_unit = '[-]';

    elseif ismember(check.st_stats, {'S'})

        st_stat =  skewness( st_C-st_Cma ,1,3);
        st_stat_title = 'skewness';
        st_stat_label = 'S';
        st_stat_unit = '[-]';

    elseif ismember(check.st_stats, {'K'})

        st_stat =  kurtosis( st_C-st_Cma ,1,3);
        st_stat_title = 'kutrosis';
        st_stat_label = 'K';
        st_stat_unit = '[-]';

    elseif ismember(check.st_stats, {'flux','F'})

        % Calculate the flux each 30 minutes
        st_fluxtime = [check.st_twindow(1) : 0.5 : check.st_twindow(2)];
        st_flux = NaN(size(st_C,1), size(st_C,2), length(st_fluxtime)-1);
        for i = 1 : length(st_fluxtime)-1
            if i == 1
                st_fluxIN = find(st_tnew <= st_fluxtime(i+1));
            else
                st_fluxIN = find(st_tnew > st_fluxtime(i) & ...
                    st_tnew <= st_fluxtime(i+1));
            end
            st_flux(:,:,i) = mean(( st_w(:,:,st_fluxIN) - mean(st_w(:,:,st_fluxIN), 3) ) ...
                .* (st_C(:,:,st_fluxIN) - mean(st_C(:,:,st_fluxIN), 3)), 3);

        end
        % Calculate the averaged flux within the time window
        st_stat =  mean(st_flux, 3);
        clear st_flux*
        st_stat_title = 'flux';
        st_stat_label = 'F';
        st_stat_unit = '[ppb m s^{-1}]';
          
        

    else
     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        fprintf(['\n     !!!!! ERROR !!!!!\n'])
        fprintf(['     "check.st_stats" input is incorrect\n'])      
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        return
    end
    
    
    
end

