%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%         Ruben Schulte         %%%%%
%%%%%     DALES analysis domain     %%%%%
%%%%%      started: 28-12-2020      %%%%%
%%%%%     restarted: 04-05-2021     %%%%%
%%%%%      changed: 05-05-2021      %%%%%
%%%%%       final: 23-05-2025       %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Note: This is a cleaned up version of Analysis_DALES_base2.m
% This script only covers the vertical profile, time series and general
% cross-section figures/analysis
%
clear all

%% General settings
% These variables can be changed by the uses

% Define the relevant (data) directories 
direc.main_directory = 'D:\DALES_data\';
% direc.main_directory = 'C:\Users\schul068\PhD project\DALES\Runs';
direc.run_name = 'nh3plume';
direc.exp_nr = {'612','613','615'};
direc.outp = 'save_output_dir'; 
check.corr_loctime = 2;

% Figure settings
check.customlgnd = 'yes';
check.exp_name = {'low-res','mid-res','high-res'}; 

%% Analysis settings

% Vertical profiles
check.vertical_profiles = 'no';    % Show vertical profile figures?
check.vp_var    = {'th','qt','u','wthlt','wthvt','nh3_r1b','wsv001t','tke'};    % Variables to plot the vertical profiles of         'wthlt','wthvt','wqtt','uwt','vwt',
% check.vp_var    = {'wsv003t','wsv001t','wsv005t','wsv007t'};    % Variables to plot the vertical profiles of         'wthlt','wthvt','wqtt','uwt','vwt',
check.vp_t      = [10 : 1 : 14];                         % Times to plot the vertical profiles (vp) [hours CEST]
check.vp_ylim   = [0 3000];


% Time series
check.time_series   = 'yes';    % Show time serie figures?
check.ts_twindow    = [8 17];
check.ts_tanalisis  = [8 12.5];
check.ts_var        = {{'ustar'},{'zi'},{'Qnet','H','LE','G0'},{'LE_frac'}};    % Variables to plot the time series off [{{'left y-axis','right y-axis'},{'left y-axis'}}]
check.ts_ylim       = [0 0.5 ...
                     ; 0 2200 ...
                     ; 0 600 ...
                     ; -0.2 1.5 ...
                     ; 0 1 ...
                       ];   % nr of y-limits has to match the number of variables

% Cross sections
check.crossection       = 'no';         % Show time serie figures?
check.cs_var            = {'nh3r1'};  % Variables to plot the cross sections
% check.cs_var            = {'thl', 'buoy'};         % Variables to plot the cross sections
check.cs_t              = [12+46/60];          % Times to plot the cros sections [hours CEST]
check.cs_sections       = {'xz'};      % Cross sections to be plotted (options: 'xy','xz','yz')
check.cs_xy_lvl         = '08';        % Vertical level of the xy cross-section
check.cs_colorscale     = 'linear';    % Options: 'logarithmic' / 'linear'
check.cs_xlim           = [0 5000];
check.cs_ylim           = [0 2100];
check.cs_Clim           = [6 12];
check.cs_Wcontour       = 'yes';
check.cs_WcontourLevel  = [0 : 0.5 : 2];   % Default = "[]", as the figure will find the levels itself


%% End of settings



% !!!!!!!!!!!!!!!!!!!! Do not change the code below !!!!!!!!!!!!!!!!!!!!



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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf('- Load the DALES data\n\n')      
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if exist('inp') == 0
    inp = struct();
    outp_gen = struct();
    outp_time = struct();
end

for II = 1 : length(direc.exp_nr)
    
    runsavename = ['set' direc.exp_nr{II}(1:end-1)];
    if ismember(runsavename, fields(inp)) == 0
        
        % Load the mandatory input, gen & time files
        try
            inp.(runsavename)           = ...
                load([direc.main_directory direc.run_name direc.exp_nr{II} direc.outp ...
                direc.exp_nr{II}(1:end-1) 'inp.mat']);
            outp_gen.(runsavename)      = ...
                load([direc.main_directory direc.run_name direc.exp_nr{II} direc.outp ...
                direc.exp_nr{II}(1:end-1) 'gen.mat']);
            outp_time.(runsavename)     = ...
                load([direc.main_directory direc.run_name direc.exp_nr{II} direc.outp ...
                direc.exp_nr{II}(1:end-1) 'time.mat']);
        catch
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            fprintf('- !!!!! ERROR !!!!! The inp, gen and/or time .mat files are missing.\n\n')      
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        end
        
        % If it exists, load the cross.mat file as well
        check.(['set' direc.exp_nr{II}(1:end-1)]).namelist = ...
            dir([direc.main_directory direc.run_name direc.exp_nr{II} direc.outp]);    
        check.(['set' direc.exp_nr{II}(1:end-1)]).namelist = ...
            {check.(['set' direc.exp_nr{II}(1:end-1)]).namelist.name}';
        check.(['set' direc.exp_nr{II}(1:end-1)]).namelist(1:2) = [];
        if ismember({[direc.exp_nr{II}(1:end-1) 'cross.mat']}, ...
                check.(['set' direc.exp_nr{II}(1:end-1)]).namelist)
            
            outp_cross.(runsavename)    = ...
                load([direc.main_directory direc.run_name direc.exp_nr{II} direc.outp ...
                direc.exp_nr{II}(1:end-1) 'cross.mat']);
            check.(['set' direc.exp_nr{II}(1:end-1)]).switch_cross = 1;
        else
            check.(['set' direc.exp_nr{II}(1:end-1)]).switch_cross = 0;
        end
        
        clear temp_*

    else    
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        fprintf(['     ' runsavename ' is already loaded\n'])      
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    end
    
    clear runsavename
end
clear old_inp II


%% Calculate tau, convective time scale

yes_no = 'no';

if strcmp(yes_no, 'yes')
    for i = 1 : length(direc.exp_nr)

        tau_g = 9.81;
        tau_tgen  = outp_gen.(['set' direc.exp_nr{i}(1:end-1)]).time;
        tau_ttime = outp_time.(['set' direc.exp_nr{i}(1:end-1)]).time;

        tau_INt = find( ismember( round(tau_ttime, 8), round(tau_tgen, 8) ) );

        tau_h = outp_time.(['set' direc.exp_nr{i}(1:end-1)]).zi;

        temp_th = NaN(size(tau_INt));
        for j = 1 : length(tau_INt)
            [tau_temp tau_INz] = sort(abs( tau_h(tau_INt(j))./2 - outp_gen.(['set' direc.exp_nr{i}(1:end-1)]).zt ));

            temp_th(j) = outp_gen.(['set' direc.exp_nr{i}(1:end-1)]).thl(tau_INz(1), j);
        end

        tau_th = NaN(size(tau_h));
        for j = 1 : length(temp_th)
            if j == 1
                tau_th(round(tau_ttime, 8) == round(tau_tgen(j), 8)) = temp_th(j);
            else
                tau_INtemp = find(round(tau_ttime, 8) >= round(tau_tgen(j-1), 8) & round(tau_ttime, 8) <= round(tau_tgen(j), 8));
                tau_temp = [temp_th(j-1) : (temp_th(j) - temp_th(j-1)) /(length(tau_INtemp)-1) : temp_th(j)];
                tau_th(tau_INtemp(2:end)) = tau_temp(2:end);
            end
        end

        tau_wtheta = outp_time.(['set' direc.exp_nr{i}(1:end-1)]).wtheta;
        tau_wstar = ( tau_g ./ tau_th .* tau_wtheta .* tau_h).^(1/3);
        tau = (tau_h ./ tau_wstar) ./60;

        outp_time.(['set' direc.exp_nr{i}(1:end-1)]).tau = tau;

        outp_time.(['set' direc.exp_nr{i}(1:end-1)]).INFO = [outp_time.(['set' direc.exp_nr{i}(1:end-1)]).INFO ; 
            table({'tau'}, {'Convective time scale'}, {'min'}, 'VariableNames',{'Name','LongName','Unit'})];

    end
end

%% Calculate average integrated TKE


yes_no = 'no';
yes_no = 'yes';

if strcmp(yes_no, 'yes')
    
    tke_twindow = [ 8,  9 ...
                 ;  9, 10 ...
                 ; 10, 11 ...
                 ; 11, 12 ...
                 ; 12, 13 ...
                 ; 13, 14 ...
                 ; 14, 15 ...
                 ; 15, 16 ...
                 ; 16, 17 ...
                 ; 14, 17 ...
                 ];
    
    tke_TKE = outp_gen.(['set' direc.exp_nr{1}(1:end-1)]).tke;
    tke_dz = median(outp_gen.(['set' direc.exp_nr{1}(1:end-1)]).zt(2:end) - ...
                    outp_gen.(['set' direc.exp_nr{1}(1:end-1)]).zt(1:end-1));
    tke_dz = double(tke_dz);
    
    tke_INTavg = NaN(size(tke_twindow, 1), 1);
    for i = 1 : size(tke_twindow, 1)
        
        TKE_INt = find(outp_gen.(['set' direc.exp_nr{1}(1:end-1)]).time > tke_twindow(i,1) & ...
            outp_gen.(['set' direc.exp_nr{1}(1:end-1)]).time <= tke_twindow(i,2) );
        
        tke_VPmean = mean(tke_TKE(:,TKE_INt),2);
        
        tke_INTavg(i,1) = sum(tke_VPmean .* tke_dz);
        
        
    end
    
    tke_INTavg
    
end



%% Add evaproation fraction to the data

evap_fields = fields(outp_time);
for i = 1 : length(evap_fields)
    
    if ismember({'LE_frac'}, outp_time.(evap_fields{i}).INFO.Name) == 0
        outp_time.(evap_fields{i}).LE_frac = outp_time.(evap_fields{i}).LE ./ ...
            (outp_time.(evap_fields{i}).LE + outp_time.(evap_fields{i}).H);

        outp_time.(evap_fields{i}).INFO = [outp_time.(evap_fields{i}).INFO ; ...
            table('LE_frac','Evaporation fraction, LvE / (LvE + H)','-', ...
            'VariableNames', {'Name','LongName','Unit'})];
    end
    
end


%% Add the G0 and An of the coldstart run to runs 612, 613, 614, 615












%%%%%
% Plot the figures
%%%%%












%% Plot vertical profile figures

if strcmp(check.vertical_profiles, 'yes')
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    fprintf('- Plot vertical profile figures\n\n')  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     close all
    
    try
        %%%%% Make the figures
        % Make a figure for each variable
        vp_plotx_pos  = 10;  % xposition of the figure [pixels]
        vp_plotdx_pos = 50;  % Movement towards the right of each subsequent figure [pixels]

        for i = 1 : length(check.vp_var)

            % Determine the position of the figure
            vp_position = [vp_plotx_pos + (i-1)*vp_plotdx_pos 50 750 700];
            % Make the figure
            f__FigProf(i, check.vp_var, check.vp_t, vp_position, outp_gen, inp, direc, check);

        end
        clear i vp_*
    catch
        fprintf('!!!!! ERROR !!!!! \n Something went wrong.\n Are you sure the NAMGENSTAT output is available?\n\n')  
        
    end
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    fprintf('- FINISHED plotting vertical profile figures\n\n')   
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
end
    

%% Plot time serie figures

if strcmp(check.time_series, 'yes')
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    fprintf('- Plot time serie figures\n\n')  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     close all
    
    try
        %%%%% Make the figures
        % Make a figure for each variable
        ts_plotx_pos  = 10;  % xposition of the figure [pixels]
        ts_plotdx_pos = 50;  % Movement towards the right of each subsequent figure [pixels]

        for i = 1 : size(check.ts_var,2)

            % Determine the position of the figure
            ts_position = [ts_plotx_pos + (i-1)*ts_plotdx_pos 230 1100 550];
            % Make the figure
            f__FigTimeSeries(i, check.ts_var, ts_position, outp_time, direc, check);

        end
        clear i ts_*
    catch
        
        fprintf('!!!!! ERROR !!!!! \n Something went wrong.\n Are you sure the NAMTIMESTAT output is available?\n\n')  
        
    end
    
    
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    fprintf('- FINISHED plotting time serie figures\n\n')   
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   
    
end

%% Cross section figures

if strcmp(check.crossection, 'yes')
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    fprintf('- Plot cross section figures\n\n')  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    
    %%%%% Make the figures
    % Make a figure for each variable
    cs_data = struct();
    cs_data.plotx_pos  = 00;  % xposition of the figure [pixels]
    cs_data.plotdx_pos = 50;  % Movement towards the right of each subsequent figure [pixels]

    cs_counter = 1;
    for i = 1 : size(check.cs_var,2)

        for II = 1 : length(direc.exp_nr)
            
            if str2num(check.cs_xy_lvl) > str2num(outp_cross.(['set' direc.exp_nr{II}(1:end-1)]).cross_lvl{end})
                check.cs_xy_lvl = data_cross.cross_lvl(length(outp_cross.(['set' direc.exp_nr{II}(1:end-1)]).cross_lvl));
            elseif str2num(check.cs_xy_lvl) < str2num(outp_cross.(['set' direc.exp_nr{II}(1:end-1)]).cross_lvl{1})
                check.cs_xy_lvl = outp_cross.(['set' direc.exp_nr{II}(1:end-1)]).cross_lvl(1);
            end

            for j = 1 : length(check.cs_t)

                % Determine the position of the figure
                cs_data.position_change = [cs_data.plotx_pos + (cs_counter-1)*cs_data.plotdx_pos ...
                    0 0 0];
                cs_data.var = check.cs_var{i};
                cs_data.tplot = check.cs_t(j);
                cs_data.exp_nr = direc.exp_nr{II};
                % Make the figure
                f__FigCrossSection(cs_data, check, direc, ...
                    outp_cross.(['set' direc.exp_nr{II}(1:end-1)]));

                cs_counter = cs_counter + 1;
            end

        end

    end

    clear cs_* i II j
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    fprintf('- FINISHED plotting cross section figures\n\n')   
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   
    
    
end

%%

'the end'






































%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%                                                       %%%%%%%%%%
%%%%%%%%%%                       Functions                       %%%%%%%%%%
%%%%%%%%%%                                                       %%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
















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
%%%%%%%%%%%%%%%%%%%%      Plot vertical profiles       %%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function f__FigProf(II, vp_var, vp_t, vp_position, outp_gen, inp, direc, check)
    
    vp_line = {'-',':','--','-.'};
    
    vp_fields = fields(inp);
    
    % Save the information of the to be plotted variables
    vp_info = table({}, {}, {}, 'VariableNames',{'Name','LongName','Unit'});
    for i = 1 : length(direc.exp_nr)
        try
            vp_info = [vp_info; outp_gen.(vp_fields{i}).INFO( ...
                ismember(outp_gen.(vp_fields{i}).INFO.Name, vp_var{II}) ,:)];
        catch
        end
    end
    % Make sure that vp_info is a table with 1 row
    if size(vp_info,1) == 0 
        vp_info = [vp_info; table({''}, {''}, {''}, 'VariableNames',{'Name','LongName','Unit'})];
    elseif size(vp_info,1) > 1 
        vp_info = vp_info(1,:);
    end
    % Sort the plotting times from low to high
    vp_t = sort(vp_t);
    
    % Find the run with the longest z vector
    vp_ylength = NaN(length(direc.exp_nr),1);
    vp_ylabel = [];
    for j = 1 : length(direc.exp_nr)
        try
            vp_ylength(j) = length(...
                outp_gen.(vp_fields{i}).zt);
            temp = [outp_gen.(vp_fields{i}).INFO{...
                ismember(outp_gen.(vp_fields{i}).INFO.Name,{'zt'}),1}{:} ...
                ' [' outp_gen.(vp_fields{i}).INFO{...
                ismember(outp_gen.(vp_fields{i}).INFO.Name,{'zt'}),3}{:} ']'];
            if length(temp) > length(vp_ylabel)
                vp_ylabel = temp;
            end
        catch
            vp_ylength(j) = 0;
        end
    end
    vp_ylength = max(vp_ylength);
    
    % Define the line colorsbd_colormap_Delta = size(bd_Bd, 2);
    vp_colormap_Delta = length(vp_t);
    if mod(vp_colormap_Delta/2, 1) == 0
        vp_colormap_temp = [[ [094 : (255 - 094)/(vp_colormap_Delta / 2) : 255 ]', ...
                         [060 : (255 - 060)/(vp_colormap_Delta / 2) : 255 ]', ...
                         [153 : (255 - 153)/(vp_colormap_Delta / 2) : 255 ]'] ...
                       ./255];
        vp_colormap_temp(end,:) = [];
        vp_colormap_temp = [vp_colormap_temp ; ...
                       [ [255-(255 - 230)/(vp_colormap_Delta / 2) : -(255 - 230)/(vp_colormap_Delta / 2) : 230 ]', ...
                         [255-(255 - 097)/(vp_colormap_Delta / 2) : -(255 - 097)/(vp_colormap_Delta / 2) : 097 ]', ...
                         [255-(255 - 001)/(vp_colormap_Delta / 2) : -(255 - 001)/(vp_colormap_Delta / 2) : 001 ]'] ./255];
    else
        vp_colormap_temp = [[ [094 : (255 - 094)/(vp_colormap_Delta / 2) : 255 ]', ...
                     [060 : (255 - 060)/(vp_colormap_Delta / 2) : 255 ]', ...
                     [153 : (255 - 153)/(vp_colormap_Delta / 2) : 255 ]'] ...
                   ./255 ; ...
                   [ [255-(255 - 230)/(vp_colormap_Delta / 2) : -(255 - 230)/(vp_colormap_Delta / 2) : 230 ]', ...
                     [255-(255 - 097)/(vp_colormap_Delta / 2) : -(255 - 097)/(vp_colormap_Delta / 2) : 097 ]', ...
                     [255-(255 - 001)/(vp_colormap_Delta / 2) : -(255 - 001)/(vp_colormap_Delta / 2) : 001 ]'] ./255];
    end
    vp_colormap_backwards = NaN(size(vp_colormap_temp));
    for j = 1 : size(vp_colormap_temp,1)
        vp_colormap_backwards(size(vp_colormap_temp,1)-j+1,:) = vp_colormap_temp(j,:);
    end

            
    
    % Build the figure
    vp_fig = figure('units','pixels', 'Color','w',...
        'innerposition', vp_position, ...
        'Name', vp_info{1,2}{:});
    % Build the plotting space
    vp_sub = subplot('Position',[0.15 0.15 0.50 0.78]);
%         vp_lgnd_line = gobjects(1,length(vp_t));

    
    % Preallocation figure variables
    vp_prof = gobjects(1,length(vp_t));
    vp_xlabel = [vp_info{1,2}{:} ' [' vp_info{1,3}{:} ']'];
    vp_ax1 = gca;   % axis variable

    % Figure settings part 1
    grid minor      % show grid
    hold on         % show multiple lines
    
    
    % Remove plotting times outside simulation time
%     vp_t(vp_t < vp_tstart | vp_t > vp_tend) = [];
    
    % Loop over all the selected times
%     vp_save = NaN(length( outp_gen.(vp_fields{1}).zt), length(vp_t));
    for i = 1 : length(vp_t)
        
        % Loop over the experiment numbers
        for j = 1 : length(direc.exp_nr)
            try
                % Organize the data of the experiment (number)
                data = outp_gen.(vp_fields{j})(:);
                data_inp = inp.(vp_fields{j})(:);
                
                % Define the starting time of the experiment
                vp_tstart = str2num(data_inp.namoptions.value{...
                    ismember(data_inp.namoptions.option, {'xtime'})...
                    }) + check.corr_loctime;
                % Define the time vector of the experiment
                vp_tday = data.time; 
                % Define the heights of the experiment
                vp_z = data.zt;

                % First profile when the first plotting time equals the 
                % initialization of the model
                if i == 1 && round(vp_t(i),6) == round(vp_tstart,6)
                    
                    if ismember(vp_var(II), fields(data_inp.prof)) % Try for a meteorological variable
                        % Set the data of the line
                        vp_plot = data_inp.prof.(vp_var{II})(:);
                    elseif ismember(vp_var(II), fields(data_inp.scalar))  % Try for a scalar
                        % Set the data of the line
                        vp_plot = data_inp.scalar.(vp_var{II})(:);
                    else
                        vp_plot = NaN(vp_ylength,1);
                    end
                        
                    
                    % Plot the line
                    plot(vp_plot, vp_z, 'LineStyle',vp_line{j}, 'Color','k', 'LineWidth', 4);
                    
                % Other profile when the first plotting time equals the 
                % initialization of the model
                elseif  ismember(round(vp_tstart,6), round(vp_t,6))
                    % Set the color of the line
%                     vp_ax1.ColorOrderIndex = i-1;
                    % Set the data of the line
                    vp_plot = data.(data.INFO.Name{...
                        ismember(data.INFO.Name, vp_var(II))})...
                        (:,ismember(round(data.time,5), round(vp_t(i),5)));
                    
                    
                    % Plot the line
                    plot(vp_plot, vp_z, 'LineStyle',vp_line{j}, 'LineWidth', 4, 'Color', vp_colormap_backwards(i-1,:));
                    
                % All profiles when the plotting times do not inlude the
                % initialization of the model
                else
                    % Set the color of the line
%                     vp_ax1.ColorOrderIndex = i;
                    % Set the data of the line
                    vp_plot = data.(vp_var{II})(:,round(data.time,6) == round(vp_t(i),6));
                    
                    % Plot the line
                    plot(vp_plot, vp_z, 'LineStyle',vp_line{j}, 'LineWidth', 4, 'Color', vp_colormap_backwards(i,:));
                end
%                 vp_save(:,i) = vp_plot;
            catch
            end
        end
            
    end
%     plot(mean(vp_save,2),data.zt, ':k','LineWidth', 3)
    % Figure settings part 2
    title(vp_info{1,2}{:})
    if sum(isnan(check.vp_ylim)) == 0
        ylim(check.vp_ylim)
    end
    xlabel(vp_xlabel);
    ylabel(vp_ylabel);
    set(vp_ax1, 'FontSize',18);

    %%%%% Set the profile timestamps
    % Preallocate time label
    vp_tlabel = repmat({''},[length(vp_t),1]);
    % Loop over the plotting times
    for j = 1 : length(vp_t)
        % Build the timestamps of the profiles
        if length(num2str(60 * mod(vp_t(j),1))) == 1
            vp_tlabel{j} = [num2str(vp_t(j)-mod(vp_t(j),1)), ...
                ':0' num2str(60 * mod(vp_t(j),1)) ' CEST'];
        else
            vp_tlabel{j} = [num2str(vp_t(j)-mod(vp_t(j),1)), ...
                ':' num2str(60 * mod(vp_t(j),1)) ' CEST'];
        end
    end

    % Positions of the legend and the timestamp annotations
    vp_lgnd_pos  = [0.65     0.90   0.14    0.055];
    vp_anno_pos  = [0.75     0.90   0.25    0.055];
    vp_lgnd_dpos = [0.00    -0.055  0.00    0.008];
    vp_anno_dpos = [0.00    -0.059  0.00    0.008];
    % Preallocate
    vp_lgnd = struct();
    % Loop over each plotting time 
    for j = 1 : length(vp_t)
        % Preallocate a dummy line
        vp_lgnd_line = gobjects(1,1);
        % Make a dummy subplot
        subplot('Position',[-j, -999, 1,1])
        hold on
        vp_temp_gca = gca;

        % Build the legend for the initialization profile if it is part
        % of the plotting times
        if vp_t(j) == vp_tstart
            % Make sure there is an initial profile defined
%             if sum(isnan(squeeze(vp_plot(i,:,j)))) ~= length(squeeze(vp_plot(i,:,j)))
                % Plot the dummy line in the dummy subplot
                vp_lgnd_line = plot([0 1],[0 1], '-k', 'LineWidth', 4);

                % Build the legend
                vp_lgnd.(['line_' num2str(j)]) = legend(vp_lgnd_line, {''});
                % Set the legend position
                vp_lgnd.(['line_' num2str(j)]).Position = vp_lgnd_pos + (j-1) .* vp_lgnd_dpos;
                % Turn the legend box off
                vp_lgnd.(['line_' num2str(j)]).Box = 'off';
                % Add the timestamp to the legend
                annotation('textbox', vp_anno_pos + (j-1) .* vp_anno_dpos, ...
                    'string', vp_tlabel{j}, ...
                    'EdgeColor','none', 'FontSize',18, 'FontWeight','normal')
                % Additionally mention that it is the initial profile
                % of the model
                annotation('textbox', vp_anno_pos + (j) .* vp_anno_dpos, ...
                    'string', '(initial profile)', ...
                    'EdgeColor','none', 'FontSize',18, 'FontWeight','normal')
%             end

        % Build the legend for the other plotting times, if the
        % initialization profile is plotted as well
        elseif  ismember(vp_tstart, vp_t) == 1
            % Set the color of the dummyline
%             vp_temp_gca.ColorOrderIndex = j-1;
            % Plot the dummy line in the dummy subplot
            vp_lgnd_line = plot([0 1],[0 1], '-', 'LineWidth', 4, 'Color', vp_colormap_backwards(j-1,:));

            % Build the legend
            vp_lgnd.(['line_' num2str(j)]) = legend(vp_lgnd_line, {''});
            % Set the legend position
            vp_lgnd.(['line_' num2str(j)]).Position = vp_lgnd_pos + (j) .* vp_lgnd_dpos;
            % Turn the legend box off
            vp_lgnd.(['line_' num2str(j)]).Box = 'off';
            % Add the timestamp to the legend
            annotation('textbox', vp_anno_pos + (j) .* vp_anno_dpos, ...
                'string', vp_tlabel{j}, ...
                'EdgeColor','none', 'FontSize',18, 'FontWeight','normal')

        % Build the legend for the other plotting times, if the
        % initialization profile is not plotted
        else
            % Set the color of the dummyline
            vp_temp_gca.ColorOrderIndex = j;
            % Plot the dummy line in the dummy subplot
            vp_lgnd_line = plot([0 1],[0 1], '-', 'LineWidth', 4, 'Color', vp_colormap_backwards(j,:));

            % Build the legend
            vp_lgnd.(['line_' num2str(j)]) = legend(vp_lgnd_line, {''});
            % Set the legend position
            vp_lgnd.(['line_' num2str(j)]).Position = vp_lgnd_pos + (j-1) .* vp_lgnd_dpos;
            % Turn the legend box off
            vp_lgnd.(['line_' num2str(j)]).Box = 'off';
            % Add the timestamp to the legend
            annotation('textbox', vp_anno_pos + (j-1) .* vp_anno_dpos, ...
                'string', vp_tlabel{j}, ...
                'EdgeColor','none', 'FontSize',18, 'FontWeight','normal')
        end
    end
    
    if length(direc.exp_nr) > 1

%         vp_runs = gobjects(1,length(direc.exp_nr));
        vp_lgndruns = struct();

        for j = 1 : length(direc.exp_nr)
            subplot('Position',[-j -j 0.50 0.78])
            set(gca, 'FontSize',18);
            vp_runs = plot([0 1], [j-1 j], 'Color', [150 150 150]./255, ...
                'LineStyle',vp_line{j}, 'LineWidth', 4);
            hold on
            
            if strcmp(check.customlgnd, 'yes')
                vp_runnames = check.exp_name(j);
            else
                vp_runnames = ['Run ' direc.exp_nr{j}(1:end-1)];
            end
            
            % Build the legend
            vp_lgndruns.(['line_' num2str(j)]) = legend(vp_runs, {''});
            % Set the legend position
            vp_lgndruns.(['line_' num2str(j)]).Position = ...
                vp_lgnd_pos - [0 0.5 0 0] + (j-1) .* vp_lgnd_dpos;
            % Turn the legend box off
            vp_lgndruns.(['line_' num2str(j)]).Box = 'off';
            % Add the timestamp to the legend
            annotation('textbox', vp_anno_pos - [0 0.5 0 0] + (j-1) .* vp_anno_dpos, ...
                    'string', vp_runnames, ...
                    'EdgeColor','none', 'FontSize',16, 'FontWeight','normal')
        end
    end
            

    
end




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%         Plot time series          %%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function f__FigTimeSeries(II, ts_var, ts_position, outp_time, direc, check)
    
    ts_colorset = {[230,097,001]./255, [094,060,153]./255,...
                   [253,184,099]./255, [178,171,210]./255};
    ts_linestyle = {'-',':','--','-.'};
    
    ts_fields = fields(outp_time);
    
    % Save the information of the to be plotted variables
    ts_datainfo = table({}, {}, {}, 'VariableNames',{'Name','LongName','Unit'});
    ts_INinfo = [];
    for i = 1 : length(direc.exp_nr)
        ts_datainfo = [ts_datainfo; ...
            outp_time.(ts_fields{i}).INFO];
    end
    ts_info_unique = unique(ts_datainfo.Name);
    for i = 1 : length(ts_info_unique)
         ts_temp = find(ismember(ts_datainfo.Name, ts_info_unique(i)));
         if isempty(ts_temp) == 0
            ts_INinfo = [ts_INinfo ; ts_temp(1)];
         end
    end
    ts_datainfo = ts_datainfo(ts_INinfo,:);
    % Define the info of the to be plotted variables
    ts_info = table({}, {}, {}, 'VariableNames',{'Name','LongName','Unit'});
    for i = 1 : length(ts_var{II})
        try
            ts_info = [ts_info; ts_datainfo( ismember( ts_datainfo.Name, ts_var{II}(i) ),: )];
        catch
            table({ts_var{II}{i}}, {'-'}, {'-'}, 'VariableNames',{'Name','LongName','Unit'})
        end
    end
    
    % Preallocation figure variables
    ts_variables = [];
    for i = 1 : length(ts_var{II})
        if i == 1
            ts_variables = ts_info.Name{i,1};
        else
            ts_variables = [ts_variables ' & ' ts_info.Name{i,1}];
        end
    end
    if length(ts_var{II}) == 1
        ts_ylabel = [ts_info.Name{1,1} ' [' ts_info.Unit{1,1} ']'];
    else
        for i = 1 : length(ts_var{II})
            if i == 1
                ts_ylabel = ts_info.Name{i,1};
            else
                ts_ylabel = [ts_ylabel ' & ' ts_info.Name{i,1}];
            end
        end
        ts_templabel = {};
        for i = 1 : length(ts_var{II})
            ts_templabel{i} = ts_info.Unit{i,1};
        end
        if length(unique(ts_templabel)) == 1
            ts_ylabel = [ ts_ylabel ' [' ts_info.Unit{1,1} ']'];
        else
            ts_ylabel = [ ts_ylabel ' [-]'];
        end
    end
    if length(ts_var{II}) == 1
        ts_title = ts_info.LongName{1,1};
    else
        for i = 1 : length(ts_var{II})
            if i == 1 
                ts_title = ts_info.LongName{i,1};
            else
                ts_title = [ts_title ' & ' ts_info.LongName{i,1}];
            end
        end
    end
    ts_xlabel = 'Time [hours CEST]';
    
    ts_lgnd_pos  = [0.73     0.85   0.17    0.055];
    ts_lgnd_dpos = [0.00    -0.085  0.00    0.005];
    ts_lgndvar = struct();
    ts_lgndrun = struct();
    
    % Build the figure
    ts_fig = figure('units','pixels', 'Color','w',...
        'innerposition', ts_position, ...
        'Name', ts_variables);
    % Build the plotting space
    ts_sub = subplot('Position',[0.095 0.15 0.68 0.78]);
    

    % Figure settings part 1
    grid on      % show grid
    hold on         % show multiple lines
    ts_ax1 = gca;   % axis variable
    
    patch([check.ts_tanalisis(1) check.ts_tanalisis(2) check.ts_tanalisis(2) check.ts_tanalisis(1)], ...
        [-9999 -9999 9999 9999], [1 1 1 1], 'FaceColor', 'k', 'FaceAlpha',0.1, 'EdgeColor','none')
    
    for i = 1 : length(direc.exp_nr)
        try
            
            data = outp_time.(ts_fields{i})(:);
            ts_INt = find(data.time >= check.ts_twindow(1) & data.time <= check.ts_twindow(2));
            data_t = data.time(ts_INt);
            
            ts_minmax = NaN(length(ts_var{II}), 2);
            for j = 1 : length(ts_var{II})
                ts_INvar =  find(ismember(data.INFO.Name, ts_var{II}(j)));
                ts_plot = data.( data.INFO.Name{ts_INvar(1)} )(ts_INt);
                
                ts_minmax(j,:) = [min(ts_plot) max(ts_plot)];
                
                plot(data_t(isnan(ts_plot)==0), ts_plot(isnan(ts_plot)==0), ...
                    'LineStyle', ts_linestyle{i}, 'Color',ts_colorset{j}, 'LineWidth', 4);
                
            end
        catch
        end
    end
%     plot([check.ts_twindow(1) check.ts_twindow(1)],[ts_ax1.YLim(1) ts_ax1.YLim(2)], 'LineWidth',2, 'LineStyle',':', 'Color','k')
%     plot([check.ts_twindow(2) check.ts_twindow(2)],[ts_ax1.YLim(1) ts_ax1.YLim(2)], 'LineWidth',2, 'LineStyle',':', 'Color','k')
    % Figure settings part 2
    ts_xlim = check.ts_twindow;
    try
        if sum( isnan(check.ts_ylim(II,:) ) ) == 0
            ts_ylim = check.ts_ylim(II,:);
        else
            ts_ylim = [min(ts_minmax(:,1)) max(ts_minmax(:,2))];
        end
    catch
        ts_ylim = [min(ts_minmax(:,1)) max(ts_minmax(:,2))];
    end
    [ts_x_tick, ts_x_tick_label, ts_y_tick, ts_y_tick_label] = ...
        f__CLASSfig_XYTicks(ts_xlim, ts_ylim);
    if ismember('obukh',ts_var{II})
        if ts_ylim(1) < -200
            ts_ylim(1) = -200;
        end
        if ts_ylim(2) > 200
            ts_ylim(2) = 200;
        end
    end
    set(ts_ax1,'FontSize',18, 'XLim',ts_xlim, 'XTick',ts_x_tick, ...
        'XTickLabel', ts_x_tick_label, ...
        'YLim',ts_ylim,'YTick',ts_y_tick, 'YTickLabel',ts_y_tick_label)
    xlabel(ts_xlabel);
    ylabel(ts_ylabel);
    title(ts_title)
    set(ts_ax1, 'FontSize',18);
    
    % Make variable legend if needed
    if length(ts_var{II}) > 1 
        for j = 1 : length(ts_var{II})
            subplot('Position',[-j -j 0.50 0.78])
            set(gca, 'FontSize',18);
            ts_line = plot([0 1], [j-1 j], 'LineStyle', ts_linestyle{1}, ...
                'Color',ts_colorset{j}, 'LineWidth', 4);
            hold on

            % Variable legend
            ts_lgndvar.(['line_' num2str(j)]) = legend(ts_line, {''});
            % Set the legend position
            ts_lgndvar.(['line_' num2str(j)]).Position = ...
                ts_lgnd_pos + [0 0.0 0 0] + (j-1) .* ts_lgnd_dpos;
            % Turn the legend box off
            ts_lgndvar.(['line_' num2str(j)]).Box = 'off';
            % Add the variable name to the legend
            annotation('textbox', ts_lgnd_pos + [0.1 0.01 0 0] + (j-1) .* ts_lgnd_dpos, ...
                    'string', ts_info.Name{j,1}, ...
                    'EdgeColor','none', 'FontSize',18, 'FontWeight','normal')
        end
    end
    % Make model run legend if needed
    if length(direc.exp_nr) > 1 
        for j = 1 : length(direc.exp_nr)
            subplot('Position',[1+j 1+j 0.50 0.78])
            set(gca, 'FontSize',18);
            ts_lgndrunline = plot([0 1], [j-1 j], 'Color', [150 150 150]./255, ...
                'LineStyle',ts_linestyle{j}, 'LineWidth', 4);
            hold on

            % Variable legend
            ts_lgndrun.(['line_' num2str(j)]) = legend(ts_lgndrunline, {''});
            % Set the legend position
            ts_lgndrun.(['line_' num2str(j)]).Position = ...
                ts_lgnd_pos + [0 -0.4 0 0] + (j-1) .* ts_lgnd_dpos;
            % Turn the legend box off
            ts_lgndrun.(['line_' num2str(j)]).Box = 'off';
            % Add the run name to the legend
            if strcmp(check.customlgnd, 'yes')
                ts_runnames = check.exp_name(j);
            else
                ts_runnames = ['Run ' direc.exp_nr{j}(1:end-1)];
            end
            annotation('textbox', ts_lgnd_pos + [0.1 -0.4 0 0] + (j-1) .* ts_lgnd_dpos, ...
                    'string', ts_runnames, ...
                    'EdgeColor','none', 'FontSize',16, 'FontWeight','normal')
        end
    end
    
end

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



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%        Plot cross sections        %%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function f__FigCrossSection(data_cs, check, direc, data_cross)
    
    cs_colorset = {[230,097,001]./255, [094,060,153]./255,...
                   [253,184,099]./255, [178,171,210]./255};
    
    cs_colormap_Delta = 1000;
    cs_colormap = [[ [094 : (255 - 094)/(cs_colormap_Delta / 2) : 255 ]', ...
                     [060 : (255 - 060)/(cs_colormap_Delta / 2) : 255 ]', ...
                     [153 : (255 - 153)/(cs_colormap_Delta / 2) : 255 ]'] ...
                   ./255 ; ...
                   [ [255-(255 - 230)/(cs_colormap_Delta / 2) : -(255 - 230)/(cs_colormap_Delta / 2) : 230 ]', ...
                     [255-(255 - 097)/(cs_colormap_Delta / 2) : -(255 - 097)/(cs_colormap_Delta / 2) : 097 ]', ...
                     [255-(255 - 001)/(cs_colormap_Delta / 2) : -(255 - 001)/(cs_colormap_Delta / 2) : 001 ]'] ./255];
    clear cs_colormap_Delta
    
    
    % Maximum size of the subplot in x and y direction
    subpos_ymax = 620;
    subpos_xmax = 750;
    subpos_Csize = 100;  % Width required for the colorbar
    
   
    
    % Define the time to be plotted
    if data_cs.tplot < data_cross.time(1)
        cs_INt = 1;
        fprintf(['!____ NOTE ____! ' ...
            '\nThe selected time is before the runtime bounds.\n'...
            'The selected time is set to ' num2str(data_cross.time(cs_INt)) ' CEST.\n\n']) 
    elseif data_cs.tplot > data_cross.time(end)
        cs_INt = length(data_cross.time);
        fprintf(['!____ NOTE ____! ' ...
            '\nThe selected time is past the runtime bounds.\n'...
            'The selected time is set to ' num2str(data_cross.time(cs_INt)) ' CEST.\n\n']) 
    else
        cs_INt = find( round(data_cross.time,6) == round(data_cs.tplot,6) );
    end
    

        
    
    
    %% Plot cross sections
    
    % Determine what cross section to be plotted
    for i = 1 : length(check.cs_sections)
        
        % Save the filename
        if strcmp('xy', check.cs_sections{i})
            cs_cross_select = [check.cs_sections{i} check.cs_xy_lvl];
            cs_filename_C = [direc.main_directory direc.run_name data_cs.exp_nr direc.outp ...
                data_cs.exp_nr(1:end-1) 'cross' cs_cross_select '_' data_cs.var '.mat'];
            cs_filename_w = [direc.main_directory direc.run_name data_cs.exp_nr direc.outp ...
                data_cs.exp_nr(1:end-1) 'cross' cs_cross_select '_w.mat'];
        elseif strcmp('xz', check.cs_sections{i})
            cs_cross_select = [check.cs_sections{i}];
            cs_filename_C = [direc.main_directory direc.run_name data_cs.exp_nr direc.outp ...
                data_cs.exp_nr(1:end-1) 'cross' check.cs_sections{i} '_' data_cs.var '.mat'];
            cs_filename_w = [direc.main_directory direc.run_name data_cs.exp_nr direc.outp ...
                data_cs.exp_nr(1:end-1) 'cross' check.cs_sections{i} '_w.mat'];
        elseif strcmp('yz', check.cs_sections{i})
            cs_cross_select = [check.cs_sections{i}];
            cs_filename_C = [direc.main_directory direc.run_name data_cs.exp_nr direc.outp ...
                data_cs.exp_nr(1:end-1) 'cross' check.cs_sections{i} '_' data_cs.var '.mat'];
            cs_filename_w = [direc.main_directory direc.run_name data_cs.exp_nr direc.outp ...
                data_cs.exp_nr(1:end-1) 'cross' check.cs_sections{i} '_w.mat'];
        end
        cs_dirlength = length([direc.main_directory direc.run_name data_cs.exp_nr direc.outp]);
        if ismember(cs_filename_C(cs_dirlength+1 : end), check.(['set' data_cs.exp_nr(1:end-1)]).namelist) == 0
            fprintf(['!____ ERROR ____! ' ...
                '\nThe selected combination of cross-section and variable does not exist in the output.' ...
                '\nSelected cross-section = ' cs_cross_select '\nSelected variable = ' data_cs.var '\n\n']) 
            return
        end
        % Determine the matfile
        cs_matfile_C = matfile(cs_filename_C);
        cs_matfile_w = matfile(cs_filename_w);
        clear cs_dirlength cs_cross_select cs_filename_*
        
        %%%%% Determine the subplot sizes
        % Make all figures on the same scale
        dxzy = [];
        % Determine the x and y pixel size of the subplots
        if strcmp('xy', check.cs_sections{i})
            if sum( isnan(check.cs_xlim) ) == 0
                cs_INx = find( data_cross.xt >= check.cs_xlim(1) - 100 & data_cross.xt <= check.cs_xlim(2) + 100 );
            else
                cs_INx = [1 : 1 : length(data_cross.xt)]';
            end
            if sum( isnan(check.cs_ylim) ) == 0
                cs_INy = find( data_cross.yt >= check.cs_ylim(1) - 100 & data_cross.yt <= check.cs_ylim(2) + 100 );
            else
                cs_INy = [1 : 1 : length(data_cross.yt)]';
            end
            
            dom_x = data_cross.xm(cs_INx(end));
            dom_y = data_cross.ym(cs_INy(end));
            dom_xy_dx = subpos_xmax / dom_x;
            dom_xy_dy = subpos_ymax / dom_y;
            dxzy = [dxzy; dom_xy_dx; dom_xy_dy];
        end
        if ismember('xz', check.cs_sections)
            if sum( isnan(check.cs_xlim) ) == 0
                cs_INx = find( data_cross.xt >= check.cs_xlim(1) - 100 & data_cross.xt <= check.cs_xlim(2) + 100 );
            else
                cs_INx = [1 : 1 : length(data_cross.xt)]';
            end
            if sum( isnan(check.cs_ylim) ) == 0
                cs_INy = find( data_cross.zt >= check.cs_ylim(1) - 100 & data_cross.zt <= check.cs_ylim(2) + 100 );
            else
                cs_INy = [1 : 1 : length(data_cross.zt)]';
            end
            
            dom_x = data_cross.xm(cs_INx(end));
            dom_z = data_cross.zm(cs_INy(end));
            dom_xz_dx = subpos_xmax / dom_x;
            dom_xz_dy = subpos_ymax / dom_z;
            dxzy = [dxzy; dom_xz_dx; dom_xz_dy];
        end
        if ismember('yz', check.cs_sections)
            if sum( isnan(check.cs_ylim) ) == 0
                cs_INx = find( data_cross.yt >= check.cs_xlim(1) - 100 & data_cross.yt <= check.cs_xlim(2) + 100 );
            else
                cs_INx = [1 : 1 : length(data_cross.yt)]';
            end
            if sum( isnan(check.cs_ylim) ) == 0
                cs_INy = find( data_cross.zt >= check.cs_ylim(1) - 100 & data_cross.zt <= check.cs_ylim(2) + 100 );
            else
                cs_INy = [1 : 1 : length(data_cross.zt)]';
            end
            
            dom_y = data_cross.ym(cs_INx(end));
            dom_z = data_cross.zm(cs_INy(end));
            dom_yz_dx = subpos_xmax / dom_y;
            dom_yz_dy = subpos_ymax / dom_z;
            dxzy = [dxzy; dom_yz_dx; dom_yz_dy];
        end
        % Find the smallest ratio to scale the figures to
        dxzy = min(dxzy);
        % scale the figures to the correct pixel/meter ratio
        if ismember('xy', check.cs_sections)
            cs_posx = dxzy * dom_x;
            cs_posy = dxzy * dom_y;
            clear dxzy dom_x dom_y dom_xy_dx dom_xy_dy
        end
        if ismember('xz', check.cs_sections)
            cs_posx = dxzy * dom_x;
            cs_posy = dxzy * dom_z;
            clear dxzy dom_x dom_z dom_xz_dx dom_xz_dy
        end
        if ismember('yz', check.cs_sections)
            cs_posx = dxzy * dom_y;
            cs_posy = dxzy * dom_z;
            clear dxzy dom_y dom_z dom_yz_dy dom_yz_dx
        end
        % Determine the size of the subplot 
        fig_subsize = [100 80 cs_posx+30 cs_posy]; 
        clear cs_posx cs_posy
        
        %%%%% Build the title string        
        % Build the timestamps
        if length(num2str(60 * mod(data_cross.time(cs_INt),1))) == 1
            cs_timestamp = [num2str(data_cross.time(cs_INt)-mod(data_cross.time(cs_INt),1)), ...
                ':0' num2str(60 * mod(data_cross.time(cs_INt),1)) ' CEST'];
        else
            cs_timestamp = [num2str(data_cross.time(cs_INt)-mod(data_cross.time(cs_INt),1)), ...
                ':' num2str(60 * mod(data_cross.time(cs_INt),1)) ' CEST'];
        end
        % Build the complete title string
        if strcmp(check.cs_sections{i}, 'xy')
            cs_title = [data_cs.exp_nr(1:end-1) ': The ' check.cs_sections{i} check.cs_xy_lvl ' cross section of ' data_cs.var ...
                ' at ' cs_timestamp ' at z = ' ...
                num2str(data_cross.zpos(ismember(data_cross.cross_lvl, {check.cs_xy_lvl}))) ' m'];
        elseif strcmp(check.cs_sections{i}, 'xz')
            cs_title = [data_cs.exp_nr(1:end-1) ': The ' check.cs_sections{i} ' cross section of ' data_cs.var ...
                ' at ' cs_timestamp ' at y = ' num2str(data_cross.ypos) ' m'];
        elseif strcmp(check.cs_sections{i}, 'yz')
            cs_title = [data_cs.exp_nr(1:end-1) ': The ' check.cs_sections{i} ' cross section of ' data_cs.var ...
                ' at ' cs_timestamp ' at x = ' num2str(data_cross.xpos) ' m'];
        end
        clear cs_timestamp
        % Build the Clabel string string
        cs_label_IN = find(ismember(data_cross.INFO.Name, {data_cs.var}));
        cs_label = [data_cs.var ' [' data_cross.INFO.Unit{cs_label_IN(1)} ']'];
        clear cs_label_IN
        
        
        
        % Definitive position and size of the figure
        fig_position = [0, 40, ...
            fig_subsize(3)+fig_subsize(1)+50+subpos_Csize-30, ...
            fig_subsize(4)+fig_subsize(2)+25]...%+20
            + data_cs.position_change;
        % Definitive position of the subplot in percentages (betweem 0 and 1)
        fig_subposition = [fig_subsize(1)/fig_position(3), ... 
                          fig_subsize(2)/fig_position(4), ...
                          fig_subsize(3)/fig_position(3), ...
                          fig_subsize(4)/fig_position(4)];
        
        %%%%% Define the variable to be plotted
        % Define the x and y axis data
        temp_x = data_cross.([check.cs_sections{i}(1) 'm'])(cs_INx);
        temp_dx = median(temp_x(2:end) - temp_x(1:end-1));
        temp_x = temp_x - temp_dx;
        temp_y = data_cross.([check.cs_sections{i}(2) 'm'])(cs_INy);
        temp_dy = median(temp_y(2:end) - temp_y(1:end-1));
        temp_c = cs_matfile_C.(data_cs.var)(cs_INx,cs_INy,cs_INt);
        
        crossX = NaN(4, length(temp_x) * length(temp_y));
        crossY = NaN(4, length(temp_x) * length(temp_y));
        crossC = NaN(length(temp_x) * length(temp_y),1);
        for j = 1 : length(temp_x)
            for k = 1 : length(temp_y)
                crossX(:,k + (j-1)*length(temp_y)) = ...
                    [temp_x(j); temp_x(j)+temp_dx; temp_x(j)+temp_dx; temp_x(j)];
                crossY(:,k + (j-1)*length(temp_y)) = ...
                    [temp_y(k); temp_y(k); temp_y(k)+temp_dy; temp_y(k)+temp_dy];
                crossC(k + (j-1)*length(temp_y)) = temp_c(j,k);
            end
        end
        clear temp_*


        
        
        %%
        
        %%%%%
        % Build the figure
        
        main_fig = figure('units','pixels', 'Color','w',...
            'innerposition', fig_position);
        % Build the plotting space
        main_sub = subplot('Position',fig_subposition);
        hold on; 
        % Plot cross section
        main_plot = patch(crossX, crossY, crossC);
        main_plot.LineWidth = 0.01;
        main_plot.EdgeAlpha = 0.2;
        
        % Plot vertical flow      
        if strcmp(check.cs_sections{i}, 'xy')
            [flowHor,flowVert] = meshgrid(data_cross.xm(cs_INx), data_cross.ym(cs_INy));
            flowHor = flowHor';  
            flowVert = flowVert';
            
        elseif strcmp(check.cs_sections{i}, 'xz')
            [flowHor,flowVert] = meshgrid(data_cross.xm(cs_INx), data_cross.zm(cs_INy));
            flowHor = flowHor';  
            flowVert = flowVert';
            
        elseif strcmp(check.cs_sections{i}, 'yz')
            [flowHor,flowVert] = meshgrid(data_cross.ym(cs_INx), data_cross.zm(cs_INy));
            flowHor = flowHor';  
            flowVert = flowVert';
            
        end
        flowW = cs_matfile_w.w(cs_INx,cs_INy,cs_INt);
        flowWup = NaN(size(flowW));
        flowWup((flowW >= 0)) = flowW(flowW >= 0);
        flowWdwn = NaN(size(flowW));
        flowWdwn((flowW < 0)) = flowW(flowW < 0);
        if strcmp(check.cs_Wcontour, 'yes')
            if isempty(check.cs_WcontourLevel) == 0
                [~,cs_windup]   = contour(flowHor, flowVert, flowWup, ...
                    'LineStyle', '-', 'Color','k', 'LineWidth',0.75, ...
                    'LevelList', check.cs_WcontourLevel);

                [~,cs_winddown] = contour(flowHor, flowVert, flowWdwn, ...
                    'LineStyle', '-', 'Color',[0.7 0.7 0.7], 'LineWidth',0.75, ...
                    'LevelList', check.cs_WcontourLevel.*(-1));

            else
                if sum(sum(isnan(flowWup))) ~= 0 
                    [~,cs_windup]   = contour(flowHor, flowVert, flowWup, ...
                        'LineStyle', '-', 'Color','k', 'LineWidth',0.5);
                end
                if sum(sum(isnan(flowWdwn))) ~= 0
                    contour(flowHor, flowVert, flowWdwn, ...
                        'LineStyle', '--', 'Color','k', 'LineWidth',0.5, ...
                        'LevelList', cs_windup.LevelList.*(-1));
                end
            end
        end
        
        % Optional logarithmic colorscale
        if strcmp(check.cs_colorscale, 'logarithmic')
        
            if sum( isnan(check.cs_Clim) ) == 0
                cr_Cminmax = check.cs_Clim;
            else
                cr_Cminmax = [min(crossC), max(crossC)];
            end
        %    rc_fig.crossLimits = [1.0*10^(0) 20]; 
%             cr_Cminmax = [cr_Cminmax(1), ceil(cr_Cminmax(2))];
%             if min(crossC ) < 0.1
%                 cr_Cminmax(1) = 0.1;
%             end

            rc_ax1 = gca;   % axis variable
            set(rc_ax1, 'ColorScale','log', 'FontSize', 18)
            
            yz_colorbar = colorbar;

            %%%%%
            % Logarithmic colors
            %%%%%
            yz_colorbar.Limits = cr_Cminmax;
            main_sub.CLim = cr_Cminmax;

            temp_TickNum = unique([1e-5:1e-5:1e-4, 1e-4:1e-4:1e-3, 1e-3:1e-3:1e-2, 0.01 : 0.01 : 0.1, 0 : 0.1 : 1, 0 : 1 : 10, 00 : 10 : 100, 0 : 100 : 1000, 0 : 1000 : 10000]);
            temp_TickNum = temp_TickNum(temp_TickNum >= cr_Cminmax(1) & temp_TickNum <= cr_Cminmax(2) );
            yz_colorbar.Ticks = temp_TickNum;
            temp_TickLabelsNum = unique([1e-5, 5e-5, 1e-4, 5e-4, 1e-3, 5e-3, 1e-2, 5e-2, 0.1, 0.5, 1, 5, 10, 50, 100, 500, 1000, 5000, 10000]);
            temp_TickLabelsNum = temp_TickLabelsNum(ismember(temp_TickLabelsNum, temp_TickNum));
            temp_TickLabelsStr = cell(length(yz_colorbar.Ticks),1);
            for j = temp_TickLabelsNum
                temp_TickLabelsStr{yz_colorbar.Ticks == j} = num2str(j);
            end
            yz_colorbar.TickLabels = temp_TickLabelsStr;

        else
            %%%%%
            % Regular colors
            %%%%%

            rc_ax1 = gca;   % axis variable
            set(rc_ax1, 'ColorScale','linear', 'FontSize', 18)

            yz_colorbar = colorbar;
            if sum( isnan(check.cs_Clim) ) == 0
                cr_Cminmax = check.cs_Clim;
            else
                cr_Cminmax = [min(crossC), max(crossC)];
            end
            yz_colorbar.Limits = cr_Cminmax;
            main_sub.CLim = cr_Cminmax;
        end
        yz_colorbar.Position = [(fig_subsize(1) + fig_subsize(3) + 10) / fig_position(3), ...
                         (fig_subsize(2)) / fig_position(4), ...
                         (30) / fig_position(3), ...
                         (fig_subsize(4)) / fig_position(4)];
        yz_colorbar.Label.String = cs_label;
        yz_colorbar.Label.String = 'NH_{3,total}';
           
        % Settings
        mainplot_gca = gca;
        set(mainplot_gca, 'FontSize', 18)
        xlabel([check.cs_sections{i}(1) ' [m]'])
        ylabel([check.cs_sections{i}(2) ' [m]'])
        if sum( isnan(check.cs_xlim) ) == 0
            xlim(check.cs_xlim)
        else
            xlim([min(crossX(:)) max(crossX(:))])
        end
        if sum( isnan(check.cs_ylim) ) == 0
            ylim(check.cs_ylim)
        else
            ylim([min(crossY(:)) max(crossY(:))])
        end
        colormap(cs_colormap);
%         annotation('textbox', [0.01 (fig_subsize(2)+fig_subsize(4)+0)/fig_position(4), 0.98 0.01], ...
%                     'string', cs_title, ...
%                     'EdgeColor','none', 'FontSize',18, 'FontWeight','normal', ...
%                     'HorizontalAlignment','center','VerticalAlignment','bottom', ...
%                     'FontWeight','Bold');
                
        
          
                
    end

end

%%





