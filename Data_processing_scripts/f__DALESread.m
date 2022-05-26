function [inp, outp_gen, outp_time, outp_crossxy, outp_crossxz, outp_crossyz, outp_field] = f__DALESread(direc, check, II)
    
    %% Initial definitions
    % Define names of relevant input files
    addons_options = {...
          'NAMPARTICLES' ...
         ;'NAMBUDGET' ...
         ;'NAMBULKMICROSTAT' ...
         ;'NAMCLOUDFIELD' ...
         ;'NAMCROSSSECTION' ...
         ;'NAMFIELDDUMP' ...
         ;'NAMGENSTAT' ...
         ;'NAMRADSTAT' ...
         ;'NAMSAMPLING' ...
         ;'NAMSTATTEND' ...
         ;'NAMTIMESTAT' ...
         };
    % Preallocate input structure
    inp = struct();
    % Preallocate output structures
    outp_gen = struct();
    outp_time = struct();
    outp_cross = struct();
    

     % [direc.main_directory, direc.run_name, direc.exp_nr{II}, direc.run_name(1:end-1), direc.exp_nr{II}(1:end-1) '.mat']

     %% Check for .mat file
        % If .mat file exists, load it
        % Else, load the DALES input and output and create a .mat file
     
    % Get the names of all the files in the directory
    name_list = dir([direc.main_directory, direc.run_name, direc.exp_nr{II}]);
    name_list = {name_list.name}';
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    fprintf(['f__DALESread: Loading the "' direc.run_name(1:end-1) '" experiment, run ' direc.exp_nr{II}(1:end-1) '.\n\n']) 
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % Check if the data is loaded before and saved as a .mat file
    if ismember({[direc.run_name(1:end-1), direc.exp_nr{II}(1:end-1) 'inp.mat']} , name_list  ) == 1 && ...
            ismember({[direc.run_name(1:end-1), direc.exp_nr{II}(1:end-1) 'outp_gen.mat']} , name_list  ) == 1 && ...
            ismember({[direc.run_name(1:end-1), direc.exp_nr{II}(1:end-1) 'outp_time.mat']} , name_list  ) == 1
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % !The data was indeed loaded before and saved as a .mat file!
        tload_start = datestr(now, 'dd-mm-yyyy HH:MM:SS');
        fprintf(['f__DALESread: Loading experiment '  ', run ' ]) 
        fprintf(['f__DALESread: The DALES run data has been processed and saved before.\n' ...
                 '              - Loading the existing file(s)\n' ...
                 '              - Starting at: ' tload_start '\n']) 
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        
        
        % Load the input data
        load([direc.main_directory, direc.run_name, direc.exp_nr{II}, ...
                     direc.run_name(1:end-1), direc.exp_nr{II}(1:end-1) 'inp.mat']);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        fprintf('                - Input files are loaded\n') 
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Load the NAMGENSTAT output data
        load([direc.main_directory, direc.run_name, direc.exp_nr{II}, ...
                     direc.run_name(1:end-1), direc.exp_nr{II}(1:end-1) 'outp_gen.mat']);
        
        gs_INwarmstart = find(ismember(inp.namoptions.option, {'lwarmstart'}));
        if isempty(gs_INwarmstart) == 0
            if strcmp(inp.namoptions.value{gs_INwarmstart}, '.true.')
                gs_warmstartfile = [inp.namoptions.value{ismember(inp.namoptions.option, {'startfile'})}(end-2:end) '\'];
                outp_gen_pre = load([direc.main_directory, direc.run_name, gs_warmstartfile, ...
                         direc.run_name(1:end-1), gs_warmstartfile(1:end-1) 'outp_gen.mat']);
                outp_gen_pre = outp_gen_pre.outp_gen;
                gs_fields = fields(outp_gen_pre);
                gs_INremove = find(ismember(gs_fields, {'README','INFO','zt','zm','zts'}));
                for j = 1 : length(gs_INremove)
                    outp_gen_pre = rmfield(outp_gen_pre, gs_fields{gs_INremove(j)});
                end
                
                gs_fields = fields(outp_gen_pre);
                for j = 1 : length(gs_fields)
                    if strcmp(gs_fields{j}, 'time')
                        outp_gen.(gs_fields{j}) = [outp_gen_pre.(gs_fields{j}); outp_gen.(gs_fields{j})];
                    else
                        outp_gen.(gs_fields{j}) = [outp_gen_pre.(gs_fields{j}), outp_gen.(gs_fields{j})];
                    end
                end
                
            end
        end
        clear gs_* outp_gen_pre
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        fprintf('                - Vertical profile output is loaded\n') 
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Load the NAMTIMESTAT output data
        load([direc.main_directory, direc.run_name, direc.exp_nr{II}, ...
                     direc.run_name(1:end-1), direc.exp_nr{II}(1:end-1) 'outp_time.mat']);
        
        
        ts_INwarmstart = find(ismember(inp.namoptions.option, {'lwarmstart'}));
        if isempty(ts_INwarmstart) == 0
            if strcmp(inp.namoptions.value{ts_INwarmstart}, '.true.')
                ts_warmstartfile = [inp.namoptions.value{ismember(inp.namoptions.option, {'startfile'})}(end-2:end) '\'];
                outp_time_pre = load([direc.main_directory, direc.run_name, ts_warmstartfile, ...
                         direc.run_name(1:end-1), ts_warmstartfile(1:end-1) 'outp_time.mat']);
                outp_time_pre = outp_time_pre.outp_time;
                ts_fields = fields(outp_time_pre);
                gs_INremove = find(ismember(ts_fields, {'README','INFO'}));
                for j = 1 : length(gs_INremove)
                    outp_time_pre = rmfield(outp_time_pre, ts_fields{gs_INremove(j)});
                end
                
                ts_fields = fields(outp_time_pre);
                for j = 1 : length(ts_fields)
                    outp_time.(ts_fields{j}) = [outp_time_pre.(ts_fields{j}); outp_time.(ts_fields{j})];
                end
                
            end
        end
        clear ts_* outp_time_pre
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        fprintf('                - Time series output is loaded\n\n') 
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        fprintf('                - Now loading the XY cross-section output\n') 
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Load the XY-crosssection data output data
        xy_temp = matfile([direc.main_directory, direc.run_name, direc.exp_nr{II}, ...
                     direc.run_name(1:end-1), direc.exp_nr{II}(1:end-1) 'outp_crossxy.mat']);
        xy_var = who(xy_temp);
        xy_var = xy_var(ismember(xy_var, [{'INFO','xt','xm','yt','ym','zc','time'}, check.xy_var_extra]));
        outp_crossxy.dir_matfile = [direc.main_directory, direc.run_name, direc.exp_nr{II}, ...
                     direc.run_name(1:end-1), direc.exp_nr{II}(1:end-1) 'outp_crossxy.mat'];
        outp_crossxy.matfile = matfile([direc.main_directory, direc.run_name, direc.exp_nr{II}, ...
                     direc.run_name(1:end-1), direc.exp_nr{II}(1:end-1) 'outp_crossxy.mat']);
        for i = 1 : length(xy_var)
            outp_crossxy.(xy_var{i}) = xy_temp.(xy_var{i});
        end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        fprintf('                --> XY cross-section output is loaded!!!\n') 
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        clear xt_temp xy_var
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        fprintf('                - Now loading the XZ cross-section output\n') 
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Load the XZ-crosssection data output data
        try
            xz_temp = matfile([direc.main_directory, direc.run_name, direc.exp_nr{II}, ...
                         direc.run_name(1:end-1), direc.exp_nr{II}(1:end-1) 'outp_crossxz.mat']);
            xz_var = who(xz_temp);
            xz_var = xz_var(ismember(xz_var, [{'INFO','xt','xm','zt','zm','zc','time'}, check.xy_var_extra]));
            outp_crossxz.dir_matfile = [direc.main_directory, direc.run_name, direc.exp_nr{II}, ...
                         direc.run_name(1:end-1), direc.exp_nr{II}(1:end-1) 'outp_crossxz.mat'];
            outp_crossxz.matfile = matfile([direc.main_directory, direc.run_name, direc.exp_nr{II}, ...
                         direc.run_name(1:end-1), direc.exp_nr{II}(1:end-1) 'outp_crossxz.mat']);
            for i = 1 : length(xz_var)
                outp_crossxz.(xz_var{i}) = xz_temp.(xz_var{i});

            end
        catch 
        end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        fprintf('                --> XZ cross-section output is loaded!!!\n') 
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        clear xz_temp xz_var
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        fprintf('                - Now loading the YZ cross-section output\n') 
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Load the YZ-crosssection data output data
%
%   NOTE: Not properly tested, since it was not used
%
%         yz_temp = matfile([direc.main_directory, direc.run_name, direc.exp_nr{II}, ...
%                      direc.run_name(1:end-1), direc.exp_nr{II}(1:end-1) 'outp_crossyz.mat']);
%         yz_var = who(yz_temp);
%         for i = 1 : length(yz_var)
%             outp_crossyz.(yz_var{i}) = yz_temp.(yz_var{i});
%             
% %             %%%%% Code in case there is too much data, so it needs to be
% %             %%%%% loaded 1 line at the time
% %             outp_crossyz.(yz_var{i}) = NaN(size(yz_temp.(yz_var{i})));
% %             for j = 1 : size(yz_temp.(yz_var{i}),1)
% %                 outp_crossyz.(yz_var{i})(j,:,:) = yz_temp.(yz_var{i})(j,:,:);
% %             end 
%         end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        fprintf('                --> YZ cross-section output is loaded!!!\n\n') 
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        clear yz_temp yz_var
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        fprintf('                - Now loading the field dump output\n') 
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Load the FIELDDUMP output data
        load([direc.main_directory, direc.run_name, direc.exp_nr{II}, ...
                     direc.run_name(1:end-1), direc.exp_nr{II}(1:end-1) 'outp_field.mat']);
        
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        fprintf('                --> Fieldd dump output is loaded!!!\n\n') 
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                 

                 
                 
               
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
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
        fprintf(['              --> DONE! Data is loaded.\n' ...
                 '                - Loading of data ended at: ' tload_end '\n' ...
                 '                - Time spent = '   num2str(tload_spent(1)) ' day(s)' ...
                 '\n                               ' num2str(tload_spent(2)),' hour(s)' ...
                 '\n                               ' num2str(tload_spent(2)),' hour(s)' ...
                 '\n                               ' num2str(tload_spent(3)),' minutes(s) ' ...
                 '\n                               ' num2str(tload_spent(4)) ' second(s)\n'])
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        % Remove temporary variables
        clear tload_* data

        fprintf(['f__DALESread: END OF FUNCTION \n'])

        return
    end
    % Remove temporary variables
    clear name_list
    
    
    
    %% In case the data was not yet processed, process the DALES data
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Message: Data needs to be loaded
    tload_start = datestr(now, 'dd-mm-yyyy HH:MM:SS');
    fprintf(['f__DALESread: No .mat files for this run exist.\n' ...
             '              - Loading the data and saving the .mat file...\n' ... 
             '              - Starting at: ' tload_start '\n\n']) 
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    
    

    %% Input files
    % Relevant input files: 1) namoptions
    %                       2) lscale.inp
    %                       3) scalar.inp
    %                       4) prof.inp
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    fprintf('              --> Loading relevant input files..\n') 
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % Get the names of all the files in the input directory
    input_namelist = dir([direc.main_directory, direc.run_name, direc.exp_nr{II}]);
    input_namelist = {input_namelist.name}';
    
    %%%%%%%%%%
    %%%%% 1) namoptions
    if ismember({['namoptions.' direc.exp_nr{II}(1:end-1)]}, input_namelist)
        try
            % Define the location/name of the namoptions file
            namoptions_dir = [direc.main_directory, direc.run_name, direc.exp_nr{II} 'namoptions.' direc.exp_nr{II}(1:end-1)];
            % Call function to read the namoptions file and add to the inp
            % structure
            inp.namoptions = f__READnamoptions(namoptions_dir);
            %
            check.modelstart = str2num( inp.namoptions.value( ismember( inp.namoptions.option, {'xtime'} ) ) );
            % Clear temporary variables
            clear namoptions_dir
        catch
            fprintf(['              --> ERROR!\n' ...
                     '                  Something went wrong when reading namoptions.' direc.exp_nr{II}(1:end-1) '.\n' ...
                     'f__DALESread: TERMINATED due to error\n'])
            return
        end
    else
        fprintf(['              --> ERROR!\n' ...
                 '                  No namoptions.' direc.exp_nr{II}(1:end-1) ' exists in the directory.\n' ...
                 'f__DALESread: TERMINATED due to error\n'])
        return
    end
    
    %%%%%%%%%%
    %%%%% 2) lscale.inp
    if ismember({['lscale.inp.' direc.exp_nr{II}(1:end-1)]}, input_namelist)
        try
            % Define the location/name of the lscale.inp file
            lscale_dir = [direc.main_directory, direc.run_name, direc.exp_nr{II} 'lscale.inp.' direc.exp_nr{II}(1:end-1)];
            % Call function to read the lscale.inp file and add to the inp
            % structure
            inp.lscale = f__READinput_textscan(lscale_dir);
            % Clear temporary variables
            clear lscale_dir
        catch
            fprintf(['            --> ERROR!\n' ...
                     '                Something went wrong when reading lscale.inp.' direc.exp_nr{II}(1:end-1) '.\n' ...
                     'f__DALESread: TERMINATED due to error\n'])
            return
        end
    else
        fprintf(['              --> ERROR!\n' ...
                 '                  No lscale.inp.' direc.exp_nr{II}(1:end-1) ' exists in the directory.\n' ...
                 'f__DALESread: TERMINATED due to error\n'])
        return
    end
    
    %%%%%%%%%%
    %%%%% 3) scalar.inp
    if ismember({['scalar.inp.' direc.exp_nr{II}(1:end-1)]}, input_namelist)
        try
            % Define the location/name of the scalar.inp file
            scalar_dir = [direc.main_directory, direc.run_name, direc.exp_nr{II} 'scalar.inp.' direc.exp_nr{II}(1:end-1)];
            % Call function to read the scalar.inp file and add to the inp
            % structure
            inp.scalar = f__READinput_textscan(scalar_dir);
            % Clear temporary variables
            clear lscalar_dir
        catch
            fprintf(['              --> ERROR!\n' ...
                     '                  Something went wrong when reading scalar.inp.' direc.exp_nr{II}(1:end-1) '.\n' ...
                     'f__DALESread: TERMINATED due to error\n'])
            return
        end
    else
        fprintf(['              --> ERROR!\n' ...
                 '                  No scalar.inp.' direc.exp_nr{II}(1:end-1) ' exists in the directory.\n' ...
                 'f__DALESread: TERMINATED due to error\n'])
        return
    end
    
    %%%%%%%%%%
    %%%%% 4) prof.inp
    if ismember({['prof.inp.' direc.exp_nr{II}(1:end-1)]}, input_namelist)
        try
            % Define the location/name of the scalar.inp file
            prof_dir = [direc.main_directory, direc.run_name, direc.exp_nr{II} 'prof.inp.' direc.exp_nr{II}(1:end-1)];
            % Call function to read the scalar.inp file and add to the inp
            % structure
            inp.prof = f__READinput_textscan(prof_dir);
            % Clear temporary variables
            clear prof_dir
        catch
            fprintf(['              --> ERROR!\n' ...
                     '                  Something went wrong when reading prof.inp.' direc.exp_nr{II}(1:end-1) '.\n' ...
                     'f__DALESread: TERMINATED due to error\n'])
            return
        end
    else
        fprintf(['              --> ERROR!\n' ...
                 '                  No prof.inp.' direc.exp_nr{II}(1:end-1) ' exists in the directory.\n' ...
                 'f__DALESread: TERMINATED due to error\n'])
        return
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    fprintf('                  Relevant input files are loaded!\n\n') 
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %% Process namoptions data
    
   
    %%%%% Save all routines defined by namoptions
        % Check if this is needed or not
    % Preallocation
    addons = struct();
    on_off = 0;
    IN_save = [];
    for i = 1 : length(inp.namoptions.option)
        if length(inp.namoptions.option{i}) > 0
            % Check if a new routine is defined in namoptions
            if strcmp(inp.namoptions.option{i}(1),'&')
                on_off = 1;
                temp_name = inp.namoptions.option{i}(2:end);

            % Check if it is the end of a routine in namoptions
            elseif strcmp(inp.namoptions.option{i}(1),'/')
                on_off = 0;
            end

            % Chack if this line states routine settings
            if on_off == 2
                % Save the index of the routine setting
                IN_save = [IN_save,i];
            % Chack if this line states the start of a routine definition
            elseif on_off == 1
                on_off = 2;
            % Chack if this line states the end of a routine definition
            elseif on_off == 0
                % Save all settings of the routine in its own table
                addons.(temp_name) = inp.namoptions(IN_save,1:2);
                % Clear index of routine settings
                IN_save = [];
            end
        end
    end
    
    addons_fields = fields(addons);
    
    %% Load output data for each routine
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    fprintf('              --> Loading relevant output files.\n') 
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    for i = 1 : length(addons_fields)
        if strcmp(addons_fields{i},'NAMPARTICLES')
            
            fprintf(['              --> The output of the NAMPARTICLES routine has not been implemented yet.\n' ...
                     '                - The script skips this routine output.\n'])
                 
        elseif strcmp(addons_fields{i},'NAMBUDGET')
            
            fprintf(['              --> The output of the NAMBUDGET routine has not been implemented yet.\n' ...
                     '                - The script skips this routine output.\n'])
                 
        elseif strcmp(addons_fields{i},'NAMBULKMICROSTAT')
            
            fprintf(['              --> The output of the NAMBULKMICROSTAT routine has not been implemented yet.\n' ...
                     '                - The script skips this routine output.\n'])
                 
        elseif strcmp(addons_fields{i},'NAMCLOUDFIELD')
            
            fprintf(['              --> The output of the NAMCLOUDFIELD routine has not been implemented yet.\n' ...
                     '                - The script skips this routine output.\n'])
                 
        elseif strcmp(addons_fields{i},'NAMCROSSSECTION')
            
            % Find the index for the on/of switch for the routine
            IN_switch = find(ismember(addons.(addons_fields{i}).option, {'lcross'}));
            
            % Check if the NAMCROSSSECTION routine switch is set to true
            if strcmp(addons.(addons_fields{i}).value(IN_switch),'.true.')
                
                fprintf(['              --> Loading output of the NAMCROSSSECTION routine.\n'])
                % Call the function that reads the data of the NAMCROSSSECTION
                % outputs.
                outp_cross = f__NAMCROSSSECTIONread(direc, inp, check, II);
                    % Output variables: structure with a separate field for 
                    % each output file, preferably filled with tables.
                    % Input variables:  direc, inp and the NAMCROSSSECTION field of
                    % addons.
                fprintf(['                  Finished loading NAMCROSSSECTION.\n'])
            end
                 
        elseif strcmp(addons_fields{i},'NAMFIELDDUMP')
            
             
            % Find the index for the on/of switch for the routine
            IN_switch = find(ismember(addons.(addons_fields{i}).option, {'lfielddump'}));
            
            % Check if the NAMCROSSSECTION routine switch is set to true
            if strcmp(addons.(addons_fields{i}).value(IN_switch),'.true.')
                
                fprintf(['              --> Loading output of the NAMFIELDDUMP routine.\n'])
                % Call the function that reads the data of the NAMCROSSSECTION
                % outputs.
                outp_field = f__NAMFIELDDUMPread(direc, inp, addons, II);
                    % Output variables: structure with a separate field for 
                    % each output file, preferably filled with tables.
                    % Input variables:  direc, inp and the NAMCROSSSECTION field of
                    % addons.
                fprintf(['                  Finished loading NAMFIELDDUMP.\n'])
            end
                 
        elseif strcmp(addons_fields{i},'NAMGENSTAT')
            
            % Find the index for the on/of switch for the routine
            IN_switch = find(ismember(addons.(addons_fields{i}).option, {'lstat'}));
            
            % Check if the NAMGENSTAT routine switch is set to true
            if strcmp(addons.(addons_fields{i}).value(IN_switch),'.true.')
                
                fprintf(['              --> Loading output of the NAMGENSTAT routine.\n'])
                % Call the function that reads the data of the NAMGENSTAT
                % outputs.
                outp_gen = f__NAMGENSTATread(direc, inp, addons, II);
                    % Output variables: structure with a separate field for 
                    % each output file, preferably filled with tables.
                    % Input variables:  direc, inp and the NAMGENSTAT field of
                    % addons.
                fprintf(['                  Finished loading NAMGENSTAT.\n'])
            end
                 
        elseif strcmp(addons_fields{i},'NAMRADSTAT')
            
            fprintf(['              --> The output of the NAMRADSTAT routine has not been implemented yet.\n' ...
                     '                - The script skips this routine output.\n'])
                 
        elseif strcmp(addons_fields{i},'NAMSAMPLING')
            
            fprintf(['              --> The output of the NAMSAMPLING routine has not been implemented yet.\n' ...
                     '                - The script skips this routine output.\n'])
                 
        elseif strcmp(addons_fields{i},'NAMSTATTEND')
            
            fprintf(['              --> The output of the NAMSTATTEND routine has not been implemented yet.\n' ...
                     '                - The script skips this routine output.\n'])
                 
        elseif strcmp(addons_fields{i},'NAMTIMESTAT')
            
            % Find the index for the on/of switch for the routine
            IN_switch = find(ismember(addons.(addons_fields{i}).option, {'ltimestat'}));
            
            % Check if the NAMGENSTAT routine switch is set to true
            if strcmp(addons.(addons_fields{i}).value(IN_switch),'.true.')
                
                fprintf(['              --> Loading output of the NAMTIMESTAT routine.\n'])
                % Call the function that reads the data of the NAMTIMESTAT
                % outputs.
                outp_time = f__NAMTIMESTATread(direc, inp, addons, II);
                    % Output variables: structure with a separate field for 
                    % each output file, preferably filled with tables.
                    % Input variables:  direc, inp and the NAMTIMESTAT field of
                    % addons.
                fprintf(['                  Finished loading NAMTIMESTAT.\n'])
            end
        end
    end
    
    

%% Combine variables

    
    if strcmp(check.combine, 'yes')
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        fprintf('\n')      
        fprintf('              --> Combine variables in the output\n') 
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
                
        
        %%%%%
        % Combine variables in NAMCROSSSECTION data
        %%%%%
        for KK = 1 : size(check.combine_vars,1)
            
            if isempty(fields(outp_cross)) == 0 && ismember({check.combine_newname{KK}} , fields(outp_cross)) == 0
                
                % Loop over all cross-sections 
                cmb_cross = fields(outp_cross);
                for kk = 1 : length(cmb_cross)
                    
                    % Define cross-section names
                    cmb_fields = fields(outp_cross.(cmb_cross{kk}));
                    % Add the crosssection indicator to the name
                    cmb_crossname = [check.combine_newname{KK} cmb_cross{kk}];

                    cmb_vars = [check.combine_vars{KK,1} cmb_cross{kk}];
                    cmb_infoIN = find( ismember(...
                        outp_cross.(cmb_cross{kk}).INFO.Name, ...
                        cmb_vars));
                    cmb_longname = outp_cross.(cmb_cross{kk}).INFO.LongName(cmb_infoIN);
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
                    outp_cross.(cmb_cross{kk}).INFO = [ outp_cross.(cmb_cross{kk}).INFO ; ...
                        table({[check.combine_newname{KK} cmb_cross{kk}]}, {cmb_tempname}, ...
                        outp_cross.(cmb_cross{kk}).INFO.Unit(cmb_infoIN), ...
                        'VariableNames',{'Name','LongName','Unit'})];
                    
                    %  Add the new variable
                    outp_cross.(cmb_cross{kk}).(cmb_crossname) = outp_cross.(cmb_cross{kk}).(cmb_vars);
                    for k = 2 : size(check.combine_vars,2)
                        if length(check.combine_vars{KK,k} > 0)
                            cmb_vars = [check.combine_vars{KK,k} cmb_cross{kk}];
                            outp_cross.(cmb_cross{kk}).(cmb_crossname) = ...
                                outp_cross.(cmb_cross{kk}).(cmb_crossname) + ...
                                outp_cross.(cmb_cross{kk}).(cmb_vars);
                        end
                    end

                end
                
            end     % isempty(fields(outp_time)) == 0 && ismember({check.combine_newname{KK}} , fields(outp_time)) == 0
        
        end     % for KK = 1 : size(check.combine_vars,1)
                
    end     % if strcmp(check.combine, 'yes')
    clear cmb_*
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    fprintf('                  Finished combining variables in the output\n\n')
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  


%% Add new variables to the NAMTIMESTAT data

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf('\n')      
fprintf('              --> Adding additional data to NAMTIMESTAT output\n\n')      
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    
    %
    %%%%% Bowen ratio
    %

    if isempty(fields(outp_time)) == 0

        try
            bowen = outp_time.H ./ outp_time.LE;
            
            outp_time.bow = NaN(size(bowen));
            outp_time.bow( abs(bowen) < 1) = bowen( abs(bowen) < 1 );

        catch
            outp_time.bow = NaN(size(outp_time.time));
        end

        outp_time.INFO = [outp_time.INFO ; ...
            table({'bow'},{'Bowen Ratio'},{'-'}, 'VariableNames',{'Name','LongName','Unit'})];
        
        clear bowen
    end
    
    
    %
    %%%%% NAMGENSTAT soil data, top layer
    %
    
    if isempty(fields(outp_time)) == 0 && isempty(fields(outp_gen)) == 0
        
        genfields = fields(outp_gen);
        var_list = {'tsoil','phiw','lambda','lambdas','gammas'};
        
        for j = 1 : length(var_list)
            try
                var_IN = find(ismember(genfields,var_list{j}));

                var_Name = outp_gen.INFO.Name{var_IN};
                var_LongName = [outp_gen.INFO.LongName{var_IN} ' (top layer)'];
                if strcmp(var_Name, 'tsoil')
                    var_Unit = 'K';
                elseif strcmp(var_Name, 'phiw')
                    var_Unit = 'm^3/m^3';
                else 
                    var_Unit = outp_gen.INFO.Unit{var_IN};
                end
                outp_time.INFO = [outp_time.INFO ; ...
                    table({var_Name},{var_LongName},{var_Unit}, 'VariableNames',{'Name','LongName','Unit'})];

                var_vptime = outp_gen.time;
                var_tstime = outp_time.time;

                outp_time.(var_Name) = NaN(size(var_tstime));

                outp_time.(var_Name)( ismember(var_tstime, var_vptime) ) = ...
                    outp_gen.(var_Name)(1,:)';
            catch
                
                var_tstime = outp_time.time;
                outp_time.(var_list{j}) = NaN(size(var_tstime));
                
                outp_time.INFO = [outp_time.INFO ; ...
                    table({var_list{j}},{'empty'},{'empty'}, 'VariableNames',{'Name','LongName','Unit'})];
            end
            
        end
        
        %
        %%%%% MXL variables
        %
        try

             mxl_vptime = outp_gen.time;
             mxl_tstime = outp_time.time;
             
             mxl_zt = outp_gen.zt;
             mxl_z = outp_time.zi;
             
             mxl_vp_t   = outp_gen.thl;
             mxl_vp_q   = outp_gen.qt;
             mxl_vp_nh3 = outp_gen.nh3;
             
             mxl_t   = NaN(size(mxl_tstime));
             mxl_q   = NaN(size(mxl_tstime));
             mxl_nh3 = NaN(size(mxl_tstime));
             
             mxlINt = find(ismember(mxl_tstime, mxl_vptime));
             
             for j = 1 : length(mxlINt)
                 mxl_t(mxlINt(j)) = nanmean( mxl_vp_t( mxl_zt <= mxl_z(mxlINt(j)), j ) );
                 mxl_q(mxlINt(j)) = nanmean( mxl_vp_q( mxl_zt <= mxl_z(mxlINt(j)), j ) );
                 mxl_nh3(mxlINt(j)) = nanmean( mxl_vp_nh3( mxl_zt <= mxl_z(mxlINt(j)), j ) );
             end
             
            outp_time.mxlthl   = mxl_t;
            outp_time.mxlqt   = mxl_q;
            outp_time.mxlnh3 = mxl_nh3;

            outp_time.INFO = [outp_time.INFO ; ...
                table({'mxlthl'},{'Mixed layer potential temperature'},{'K'}, 'VariableNames',{'Name','LongName','Unit'})];
            outp_time.INFO = [outp_time.INFO ; ...
                table({'mxlqt'},{'Mixed layer specific humidity'},{'kg/kg'}, 'VariableNames',{'Name','LongName','Unit'})];
            outp_time.INFO = [outp_time.INFO ; ...
                table({'mxlnh3'},{'Mixed layer nh3 specific mixing ration'},{'kg/kg'}, 'VariableNames',{'Name','LongName','Unit'})];
             
        catch
            mxl_tstime = outp_time.time;
            outp_time.mxlthl = NaN(size(mxl_tstime));
            outp_time.mxlqt = NaN(size(mxl_tstime));

            outp_time.INFO = [outp_time.INFO ; ...
                table({'mxlthl'},{'empty'},{'empty'}, 'VariableNames',{'Name','LongName','Unit'})];
            outp_time.INFO = [outp_time.INFO ; ...
                table({'mxlqt'},{'empty'},{'empty'}, 'VariableNames',{'Name','LongName','Unit'})];
            outp_time.INFO = [outp_time.INFO ; ...
                table({'mxlnh3'},{'empty'},{'empty'}, 'VariableNames',{'Name','LongName','Unit'})];
        end
        
        clear j genfields var_*
        
    elseif isempty(fields(outp_time)) == 0
        
        var_list = {'tsoil','phiw','lambda','lambdas','gammas'};
        
        for j = 1 : length(var_list)
        
            var_tstime = outp_time.time;
            outp_time.(var_list{j}) = NaN(size(var_tstime));

            outp_time.INFO = [outp_time.INFO ; ...
                table({var_list{j}},{'empty'},{'empty'}, 'VariableNames',{'Name','LongName','Unit'})];
        end
        
        mxl_tstime = outp_time.time;
        outp_time.mxlthl = NaN(size(mxl_tstime));
        outp_time.mxlqt =  NaN(size(mxl_tstime));
        outp_time.mxlnh3 = NaN(size(mxl_tstime));

        outp_time.INFO = [outp_time.INFO ; ...
            table({'mxlthl'},{'empty'},{'empty'}, 'VariableNames',{'Name','LongName','Unit'})];
        outp_time.INFO = [outp_time.INFO ; ...
            table({'mxlqt'},{'empty'},{'empty'}, 'VariableNames',{'Name','LongName','Unit'})];
        outp_time.INFO = [outp_time.INFO ; ...
            table({'mxlnh3'},{'empty'},{'empty'}, 'VariableNames',{'Name','LongName','Unit'})];
    end
    
    
    
    if isempty(fields(outp_gen)) == 0
        
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


        %
        %%%%% Buoyancy resolved over total
        %
        outp_gen.wthv_rot = abs(outp_gen.wthvr) ./ abs(outp_gen.wthvt);
        
        outp_gen.INFO = [outp_gen.INFO ; ...
            table({'wthv_rot'},{'Resolved over total buoyancy'},{'Km/s'}, ...
            'VariableNames',{'Name','LongName','Unit'})];
        
        %
        %%%%% heat flux resolved over total
        %
        outp_gen.wthl_rot = abs(outp_gen.wthlr) ./ abs(outp_gen.wthlt);
        
        outp_gen.INFO = [outp_gen.INFO ; ...
            table({'wthl_rot'},{'Resolved over total Theta_l flux'},{'Km/s'}, ...
            'VariableNames',{'Name','LongName','Unit'})];
        
    end
    
    
    
    
    %% Save the data
   
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    fprintf('              --> Saving the input and output data in a .mat file\n') 
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %%%%%
    % Save the data
    save([direc.main_directory direc.run_name direc.exp_nr{II} ...                  % Save directory
        direc.run_name(1:end-1), direc.exp_nr{II}(1:end-1) 'inp.mat'], ...          % Save name
        'inp', '-v7.3');                                                        % Saved variables
    save([direc.main_directory direc.run_name direc.exp_nr{II} ...                  % Save directory
        direc.run_name(1:end-1), direc.exp_nr{II}(1:end-1) 'outp_gen.mat'], ...     % Save name
        'outp_gen', '-v7.3');                                                   % Saved variables
    save([direc.main_directory direc.run_name direc.exp_nr{II} ...              	% Save directory
        direc.run_name(1:end-1), direc.exp_nr{II}(1:end-1) 'outp_time.mat'], ...    % Save name
        'outp_time', '-v7.3');                                                  % Saved variables
    save([direc.main_directory direc.run_name direc.exp_nr{II} ...                  % Save directory
        direc.run_name(1:end-1), direc.exp_nr{II}(1:end-1) 'outp_cross.mat'], ...   % Save name
        'outp_cross', '-v7.3');                                                 % Saved variables
    save([direc.main_directory direc.run_name direc.exp_nr{II} ...                  % Save directory
        direc.run_name(1:end-1), direc.exp_nr{II}(1:end-1) 'outp_field.mat'], ...   % Save name
        'outp_field', '-v7.3');                                                 % Saved variables
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    fprintf('                  DATA IS SAVED!\n') 
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    
    %% End of function
    
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
    fprintf(['              --> DONE! Data is LOADED and SAVED.\n' ...
             '                - function ended at: ' tload_end '\n' ...
             '                - Time spent = '   num2str(tload_spent(1)) ' day(s)' ...
             '\n                               ' num2str(tload_spent(2)),' hour(s)' ...
             '\n                               ' num2str(tload_spent(2)),' hour(s)' ...
             '\n                               ' num2str(tload_spent(3)),' minutes(s) ' ...
             '\n                               ' num2str(tload_spent(4)) ' second(s)\n' ...
             'f__DALESread: END OF FUNCTION\n\n'])
    % Remove temporary variables
    clear tload_* data
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
end


























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


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%% READ input with textscan function %%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function func_out = f__READinput_textscan(lscale_dir)
   
    % Open the lscale.inp file
    fID = fopen(lscale_dir);
    
    % Get the Nz of the run
        % First header delimiter is a tab
    firstline_delimiter = {'\t'};
        % Read the first line, ignoring spaces and # 
    firstline = textscan(fID, '%s %f %s %f', 1, 'Delimiter', firstline_delimiter, 'TreatAsEmpty', {' ','#'});
        % Save the Nz value
    Nz = firstline{2};
    
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

function f_output = f__NAMGENSTATread(direc, inp, addons, II)
    % This function reads the data of the NAMGENSTAT outputs.
        % Input variables:  direc, inp and the NAMGENSTAT field of addons
        % Output variables: structure with a separate field for each output 
        %                   file, preferably filled with tables.
    
    
    % Preallocation output structure
    f_output = struct();
    % Clarify how the output is saved
    f_output.README = {'Each variable is saved as an [height x time] function'};
    % Preallocate the info on the variables
    f_output.INFO = table({},{},{}, 'VariableNames',{'Name','LongName','Unit'});
    
    %% Read the NetCDF output
    
    % Make sure the "lnetcdf" switch is not set to " .false." 
    nc_exists = 1;  % 1 = NetCDF files exist, 0 = NetCDF files do not exist
    if ismember({'NAMNETCDFSTATS'}, fields(addons))
        if strcmp(addons.NAMNETCDFSTATS.value{...
                ismember(addons.NAMNETCDFSTATS.option,{'lnetcdf'})...
                },'.false.')
            nc_exists = 0;
        end
    end
    
    
    if nc_exists == 1
        % Define location of the NetCDF file
        nc_loc = [direc.main_directory, direc.run_name, direc.exp_nr{II} 'profiles.' direc.exp_nr{II}(1:end-1) '.nc'];
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
            f_output.(nc_vars{i}) = ncread(nc_loc, nc_vars{i});
            % Read the variable information and store it in a temporary
            % variable
            temp_units = ncinfo(nc_loc, nc_vars{i});
            nc_units(i,1:3) = [{temp_units.Name}, ...
                               {temp_units.Attributes(1,1).Value}, ...
                               {temp_units.Attributes(1,2).Value}];
        end
        % Save the variable information in the temporary variable and store it
        % in the data sructure
        f_output.INFO = [f_output.INFO ; nc_units];
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
    temp_longnames = (f_output.INFO.LongName);
    % Loop over the input names of the scalars to correct all of them
    temp_varfields = fields(f_output);
    temp_varfields(ismember(temp_varfields, {'README','INFO'})) = [];
    for i = 1 : length(corr_inpscalars)
        
        
        % Create the name that will be used to search in f_output.(cross_fields).INFO.LongNames
        corr_name = '000';
        corr_name(end-(length(num2str(i))-1) : end) = num2str(i);

        corr_name = [{['sv' corr_name]},{['Scalar ' corr_name]}]; 

        corr_IN = find(ismember(f_output.INFO.Name, corr_name(1)));
        
        if isempty(corr_IN) == 0
            f_output.INFO.Name{corr_IN} = corr_inpscalars{i};

            corr_longname = f_output.INFO.LongName{corr_IN};
            for j = 1 : length(corr_longname) - length(corr_name{2}) + 1
                if strcmp(corr_name{2}, corr_longname(j : j+length(corr_name{2})-1))
                    if j == 1
                        f_output.INFO.LongName{corr_IN} = ...
                            ['Scalar ', corr_inpscalars{i}, ...
                            corr_longname(j+length(corr_name{2}):end)];
                    else
                        f_output.INFO.LongName{corr_IN} = ...
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
        outp_fields_old = fields(f_output);
        % Save the original f_output
        f_output_old = f_output;
        % Save the original field names of f_output

        % Base the new f_output field names on the old f_output field names
        outp_fields_new = outp_fields_old;
        % Apply the corrections to the f_output field names
        outp_fields_new(ismember(outp_fields_old, corr_corrections(:,1))) = ...
            corr_corrections(:,2);

        % Build a new and empty f_output structure
        f_output = struct();
        % Redefine f_output with the same data, but new names for scalars
        for i = 1 : length(outp_fields_new)
            f_output.(outp_fields_new{i}) = f_output_old.(outp_fields_old{i});
        end
    end
    
    
    %% Read the txt output
    
    % Define the location & name of the txt file you want to read
    txt_loc = [direc.main_directory, direc.run_name, direc.exp_nr{II} 'field.' direc.exp_nr{II}(1:end-1)];
    
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
    txt_ln_data = str2num(addons.DOMAIN.value(ismember(addons.DOMAIN.option, {'kmax'})));
    % Define the number of times the data is written by DALES
    txt_ln_time = str2num(addons.RUN.value(ismember(addons.RUN.option, {'runtime'}))) / ...
                  str2num(addons.NAMGENSTAT.value(ismember(addons.NAMGENSTAT.option, {'timeav'})));
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
        temp_data = textscan(txt_fID, '%f %f %f %f %f %f %f %f %f %f %f %f %f', ...
            txt_ln_data, 'Delimiter', txt_delimiter, 'HeaderLines',txt_skip);
        % Close the file
        fclose(txt_fID);
        
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
        f_output.temp = outp_temp;
        f_output.INFO = [f_output.INFO; table({'temp'}, {'Temperature'}, {'K'}, ...
            'VariableNames',{'Name','LongName','Unit'})];
    end
    
    % Add the potential temperature data to f_output
    if sum(ismember(txt_header, 'THETA')) == 1
        outp_temp = NaN(txt_ln_data, txt_ln_time);
        for i = 1 : txt_ln_time
            outp_temp(:,i) = txt_data(:,ismember(txt_header, 'THETA'),i);
        end
        f_output.th = outp_temp;
        f_output.INFO = [f_output.INFO; table({'th'}, {'Potential temperature'}, {'K'}, ...
            'VariableNames',{'Name','LongName','Unit'})];
    end
    
    % Add the cloud fraction data to f_output
    if sum(ismember(txt_header, 'CLOUD_FRACTION')) == 1
        outp_temp = NaN(txt_ln_data, txt_ln_time);
        for i = 1 : txt_ln_time
            outp_temp(:,i) = txt_data(:,ismember(txt_header, 'CLOUD_FRACTION'),i);
        end
        f_output.cf = outp_temp;
        f_output.INFO = [f_output.INFO; table({'cf'}, {'Cloud fraction'}, {'-'}, ...
            'VariableNames',{'Name','LongName','Unit'})];
    end
    
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%      READ NAMTIMESTAT output      %%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function f_output = f__NAMTIMESTATread(direc, inp, addons, II)
    % This function reads the data of the NAMTIMESTAT outputs.
        % Input variables:  direc, inp and the NAMTIMESTAT field of addons
        % Output variables: structure with a separate field for each output 
        %                   file, preferably filled with tables.
    
    
    % Preallocation output structure
    f_output = struct();
    % Clarify how the output is saved
    f_output.README = {'Each variable is saved as an [time x 1] vector'};
    % Preallocate the info on the variables
    f_output.INFO = table({},{},{}, 'VariableNames',{'Name','LongName','Unit'});
    
    %% Read the NetCDF output
    
    % Make sure the "lnetcdf" switch is not set to " .false." 
    nc_exists = 1;  % 1 = NetCDF files exist, 0 = NetCDF files do not exist
    if ismember({'NAMNETCDFSTATS'}, fields(addons))
        if strcmp(addons.NAMNETCDFSTATS.value{...
                ismember(addons.NAMNETCDFSTATS.option,{'lnetcdf'})...
                },'.false.')
            nc_exists = 0;
        end
    end
    
    
    if nc_exists == 1
        % Define location of the NetCDF file
        nc_loc = [direc.main_directory, direc.run_name, direc.exp_nr{II} 'tmser.' direc.exp_nr{II}(1:end-1) '.nc'];
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
            f_output.(nc_vars{i}) = ncread(nc_loc, nc_vars{i});
            % Read the variable information and store it in a temporary
            % variable
            temp_units = ncinfo(nc_loc, nc_vars{i});
            nc_units(i,1:3) = [{temp_units.Name}, ...
                               {temp_units.Attributes(1,1).Value}, ...
                               {temp_units.Attributes(1,2).Value}];
        end
        % Save the variable information in the temporary variable and store it
        % in the data sructure
        f_output.INFO = [f_output.INFO ; nc_units];
    else
    
        fprintf(['              --> There is no NetCDF output.\n' ...
                 '                - Please turn on the "lnetcdf" switch in Namoptions.\n'])
    end
    
    %% Read the txt output to add the missing variables
    
    % Define the location & name of the txt file you want to read
    txt_loc = [direc.main_directory, direc.run_name, direc.exp_nr{II} 'tmlsm.' direc.exp_nr{II}(1:end-1)];
    
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
    txt_ln_data = str2num(addons.RUN.value(ismember(addons.RUN.option, {'runtime'}))) / ...
                  str2num(addons.NAMTIMESTAT.value(ismember(addons.NAMTIMESTAT.option, {'dtav'})));
    
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
        txt_data(:,i) = temp_data{:,i};
    end
    
    txt_add  = {'tskin';'Resp';'wco2';'An';'gcco2'};
    txt_add_long = {'Skin temperature'; ...
                    'CO2 soil respiration'; ...
                    'Net CO2 flux'; ...
                    'Assimilation rate'; ...
                    'CO2 canopy conductance'};
    for i = 1 : length(txt_add)
        f_output.(txt_add{i}) = ...
            txt_data(:, ismember(txt_header(:,1),txt_add(i)));
        f_output.INFO = [f_output.INFO ; ...
            table(txt_add(i),txt_add_long(i),txt_header(ismember(txt_header(:,1),txt_add(i)),2), ...
            'VariableNames',{'Name','LongName','Unit'})];
    end
    
    
    %% Rename the scalars in the output to match the names in scalar.inp
    
    % Get the input names of the scalars
    corr_inpscalars = inp.scalar.Properties.VariableNames(2:end);
    
    % Preallocation of matrix with names to correct
    corr_corrections = [{}, {}]; 
    % Get the LongNames as saved in f_output.INFO (will be used to search in)
    temp_longnames = (f_output.INFO.LongName);
    % Loop over the input names of the scalars to correct all of them
    for i = 1 : length(corr_inpscalars)
        
        
        % Create the name that will be used to search in f_output.(cross_fields).INFO.LongNames
        corr_name = '000';
        corr_name(end-(length(num2str(i))-1) : end) = num2str(i);

        corr_name = [{['sv' corr_name]},{['Scalar ' corr_name]}]; 

        corr_IN = find(ismember(f_output.INFO.Name, corr_name(1)));
        
        if isempty(corr_IN) == 0
            f_output.INFO.Name{corr_IN} = corr_inpscalars{i};

            corr_longname = f_output.INFO.LongName{corr_IN};
            for j = 1 : length(corr_longname) - length(corr_name{2}) + 1
                if strcmp(corr_name{2}, corr_longname(j : j+length(corr_name{2})-1))
                    if j == 1
                        f_output.INFO.LongName{corr_IN} = ...
                            ['Scalar ', corr_inpscalars{i}, ...
                            corr_longname(j+length(corr_name{2}):end)];
                    else
                        f_output.INFO.LongName{corr_IN} = ...
                            [corr_longname(1:j-1), ' scalar ', corr_inpscalars{i}, ...
                            corr_longname(j+length(corr_name{2}):end)];
                    end
                end
            end
        end
        
        
    end
    
    if isempty(corr_corrections) == 0
        % Save the original field names of f_output
        outp_fields_old = fields(f_output);
        % Save the original f_output
        f_output_old = f_output;
        % Save the original field names of f_output

        % Base the new f_output field names on the old f_output field names
        outp_fields_new = outp_fields_old;
        % Apply the corrections to the f_output field names
        outp_fields_new(ismember(outp_fields_old, corr_corrections(:,1))) = ...
            corr_corrections(:,2);

        % Build a new and empty f_output structure
        f_output = struct();
        % Redefine f_output with the same data, but new names for scalars
        for i = 1 : length(outp_fields_new)
            f_output.(outp_fields_new{i}) = f_output_old.(outp_fields_old{i});
        end
    end
    
    
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%    READ NAMCROSSSECTION output    %%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function f_output = f__NAMCROSSSECTIONread(direc, inp, check, II)
    
    
    loadvars = {'w';'thl';'qt'};
    basevars_sv = str2num(inp.namoptions.value(ismember(inp.namoptions.option,{'nsv'})));
    for i = 1 : basevars_sv
        if i >= 100
            loadvars = [loadvars; {['sv' num2str(i)]}];
        elseif i >= 10
            loadvars = [loadvars; {['sv0' num2str(i)]}];
        else
            loadvars = [loadvars; {['sv00' num2str(i)]}];
        end
    end
    loadvars = [loadvars ; check.addcrossvars];
    
    clear check
        
    %%%%%
    % Setup the information on the CROSSSECTION output
    
    % Preallocation output structure
    f_output = struct();

    %% X-Y cross section
    
    xy_loadvars = cell(size(loadvars));
    for i = 1 : length(loadvars)
        xy_loadvars{i} = [loadvars{i} 'xy'];
    end
    
    % Find output files with an X-Y cross section
    files_xy = dir([direc.main_directory, direc.run_name, direc.exp_nr{II}, ...
        'crossxy*']);
    if isempty(files_xy) == 0 
        files_xy = {files_xy.name}';
        
        % Preallocate the info on the variables
        f_output.xy.INFO = table({},{},{}, 'VariableNames',{'Name','LongName','Unit'});
        
        
        % Determine the height at which the cross section is taken
        xy_z = inp.prof.zc(str2num(files_xy{1}(9:12)));
        
       
        % Load the first file as a dummy
        xy_loc = [direc.main_directory, direc.run_name, direc.exp_nr{II}, ...
            files_xy{1}];
        xy_info = ncinfo(xy_loc);
        % Define variables to be saved
        xy_var = {xy_info.Variables.Name}';
        % Define time vatiable
        xy_time = ncread(xy_loc, 'time');
        % Preallocate variable sizes
        xy_dimx = xy_info.Variables(ismember({xy_info.Variables.Name},{'xt'})).Size(1);
        xy_nprocx = str2num(inp.namoptions.value(ismember(inp.namoptions.option, 'nprocx')));
        xy_dimy = xy_info.Variables(ismember({xy_info.Variables.Name},{'yt'})).Size(1);
        xy_nprocy = str2num(inp.namoptions.value(ismember(inp.namoptions.option, 'nprocy')));
        xy_dimt = length(xy_time);
        
        % Fill the INFO field
        for i = 1 : size(xy_info.Variables,2)
            
            f_output.xy.INFO = [f_output.xy.INFO ; table({...
                xy_var{i}}, {...
                xy_info.Variables(i).Attributes(1).Value}, {...
                xy_info.Variables(i).Attributes(2).Value}, ...
                'VariableNames',{'Name','LongName','Unit'})];
        end
        f_output.xy.INFO = [f_output.xy.INFO ; ...
            table({'zc'}, {'Cross section height'}, {'m'}, ...
                'VariableNames',{'Name','LongName','Unit'})];
        
        % Preallocation data
        for i = 1 : length(xy_var)
            if strcmp(xy_var{i},'time') == 0
                if ismember(xy_var(i),{'xt','xm'})
                    f_output.xy.(xy_var{i}) = NaN(xy_dimx*xy_nprocx, 1);
                elseif ismember(xy_var(i),{'yt','ym'})
                    f_output.xy.(xy_var{i}) = NaN(xy_dimy*xy_nprocy,1);                       
                elseif ismember(xy_var(i), xy_loadvars)
                    f_output.xy.(xy_var{i}) = ...
                        NaN(xy_dimx*xy_nprocx, xy_dimy*xy_nprocy, xy_dimt);
                end
            else
                f_output.xy.(xy_var{i}) = ncread(xy_loc, xy_var{i});
            end
        end
        
        
        % Determine the number of files
        xy_files_n1 = length(files_xy);
        xy_files_n2 = xy_nprocx * xy_nprocy;
        if xy_files_n1 == xy_files_n2
            for i = 1 : xy_files_n1
                % Directory location
                xy_loc = [direc.main_directory, direc.run_name, direc.exp_nr{II}, ...
                    files_xy{i}];
                
                xy_locx = str2num(files_xy{i}(15:17));
                xy_locy = str2num(files_xy{i}(19:21));
                
                xy_fillx = [1+xy_dimx*xy_locx : 1+xy_dimx*(xy_locx+1)-1];
                xy_filly = [1+xy_dimx*xy_locy : 1+xy_dimy*(xy_locy+1)-1];

%                 j = j+1
                for j = 1 : length(xy_var)
                    
                    if ismember(xy_var(j),{'xt','xm'}) 
                        f_output.xy.(xy_var{j})(xy_fillx) = ...
                            ncread(xy_loc, xy_var{j});
                    elseif ismember(xy_var(j),{'yt','ym'})
                        f_output.xy.(xy_var{j})(xy_filly) = ...
                            ncread(xy_loc, xy_var{j});
                    elseif ismember(xy_var(j),xy_loadvars)
                        f_output.xy.(xy_var{j})(xy_fillx, xy_filly, :) = ...
                            ncread(xy_loc, xy_var{j});
                    end
                    
                end
            end
        end
        f_output.xy.zc = xy_z;
        
    end
    clear xy_* i j files_xy
    
    
    %% X-Z cross section
    
    xz_loadvars = cell(size(loadvars));
    for i = 1 : length(loadvars)
        xz_loadvars{i} = [loadvars{i} 'xz'];
    end
    
    % Find output files with an X-Y cross section
    files_xz = dir([direc.main_directory, direc.run_name, direc.exp_nr{II}, ...
        'crossxz.x*']);
    if isempty(files_xz) == 0 
        files_xz = {files_xz.name}';
        
        %%%%%
        % Setup the information on the CROSSSECTION output
        
        % Preallocate the info on the variables
        f_output.xz.INFO = table({},{},{}, 'VariableNames',{'Name','LongName','Unit'});
        
        
        % Determine the y location at which the cross section is taken
        try
            xz_dy = str2num(inp.namoptions.value(ismember(inp.namoptions.option, {'ysize'}))) ...
                / str2num(inp.namoptions.value(ismember(inp.namoptions.option, {'jtot'})));
            xz_y = str2num(inp.namoptions.value(ismember(inp.namoptions.option, {'crossplane'}))) ...
                * xz_dy - 0.5 * xz_dy;
        catch
            xz_dy = str2num(inp.namoptions.value(ismember(inp.namoptions.option, {'ysize'}))) ...
                / str2num(inp.namoptions.value(ismember(inp.namoptions.option, {'jtot'})));
            xz_y = 2 * xz_dy - 0.5 * xz_dy;
        end
       
        % Load the first file as a dummy
        xz_loc = [direc.main_directory, direc.run_name, direc.exp_nr{II}, ...
            files_xz{1}];
        xz_info = ncinfo(xz_loc);
        % Define variables to be saved
        xz_var = {xz_info.Variables.Name}';
        % Define time vatiable
        xz_time = ncread(xz_loc, 'time');
        % Preallocate variable sizes
        xz_dimx = xz_info.Variables(ismember({xz_info.Variables.Name},{'xt'})).Size(1);
        xz_nprocx = str2num(inp.namoptions.value(ismember(inp.namoptions.option, 'nprocx')));
        xz_dimz = xz_info.Variables(ismember({xz_info.Variables.Name},{'zt'})).Size(1);
        xz_nprocz = str2num(inp.namoptions.value(ismember(inp.namoptions.option, 'kmax')));
        xz_dimt = length(xz_time);
        
        % Fill the INFO field
        for i = 1 : size(xz_info.Variables,2)
            
            f_output.xz.INFO = [f_output.xz.INFO ; table({...
                xz_var{i}}, {...
                xz_info.Variables(i).Attributes(1).Value}, {...
                xz_info.Variables(i).Attributes(2).Value}, ...
                'VariableNames',{'Name','LongName','Unit'})];
        end
        f_output.xz.INFO = [f_output.xz.INFO ; ...
            table({'yt'}, {'Cross section y location'}, {'m'}, ...
                'VariableNames',{'Name','LongName','Unit'})];
        
        % Preallocation data
        for i = 1 : length(xz_var)
            
            if strcmp(xz_var{i},'time')
                
                f_output.xz.(xz_var{i}) = ncread(xz_loc, xz_var{i});
                
            elseif ismember(xz_var(i),{'zt','zm'})
                
                f_output.xz.(xz_var{i}) = double(ncread(xz_loc, xz_var{i}));
                
            elseif ismember(xz_var(i),{'xt','xm'})
                
                f_output.xz.(xz_var{i}) = NaN(xz_dimx*xz_nprocx, 1);     
                
            elseif ismember(xz_var(i), xz_loadvars)
                
                f_output.xz.(xz_var{i}) = NaN(xz_dimx*xz_nprocx, xz_dimz, xz_dimt);
                
            end
        end
        
        
        % Determine the number of files
        xz_files_n1 = length(files_xz);
        xz_files_n2 = xz_nprocx;
        if xz_files_n1 == xz_files_n2
            for i = 1 : xz_files_n1
                % Directory location
                xz_loc = [direc.main_directory, direc.run_name, direc.exp_nr{II}, ...
                    files_xz{i}];
                
                xz_locx = str2num(files_xz{i}(10:12));
                
                xz_fillx = [1+xz_dimx*xz_locx : 1+xz_dimx*(xz_locx+1)-1];

%                 j = j+1
                for j = 1 : length(xz_var)
                    
                    if ismember(xz_var(j),{'xt','xm'}) 
                        f_output.xz.(xz_var{j})(xz_fillx) = ...
                            ncread(xz_loc, xz_var{j});
                    elseif ismember(xz_var(j), xz_loadvars)
                        f_output.xz.(xz_var{j})(xz_fillx, :, :) = ...
                            ncread(xz_loc, xz_var{j});
                    end
                    
                end
            end
        end
        f_output.xz.yt = xz_y;
        
    end
    clear xz_* i j files_xz
    
    
    %% Y-Z cross section
    
    yz_loadvars = cell(size(loadvars));
    for i = 1 : length(loadvars)
        yz_loadvars{i} = [loadvars{i} 'yz'];
    end
    
    % Find output files with an X-Y cross section
    files_yz = dir([direc.main_directory, direc.run_name, direc.exp_nr{II}, ...
        'crossyz.x*']);
    if isempty(files_yz) == 0 
        files_yz = {files_yz.name}';
        
        %%%%%
        % Setup the information on the CROSSSECTION output
        
        % Preallocate the info on the variables
        f_output.yz.INFO = table({},{},{}, 'VariableNames',{'Name','LongName','Unit'});
        
        
        % Determine the y location at which the cross section is taken
        try
            yz_dx = str2num(inp.namoptions.value(ismember(inp.namoptions.option, {'xsize'}))) ...
                / str2num(inp.namoptions.value(ismember(inp.namoptions.option, {'itot'})));
            yz_x = str2num(inp.namoptions.value(ismember(inp.namoptions.option, {'crossortho'}))) ...
                * yz_dx - 0.5 * yz_dx;
        catch
            yz_dx = str2num(inp.namoptions.value(ismember(inp.namoptions.option, {'xsize'}))) ...
                / str2num(inp.namoptions.value(ismember(inp.namoptions.option, {'itot'})));
            yz_x = 2 * yz_dx - 0.5 * yz_dx;
        end
       
        % Load the first file as a dummy
        yz_loc = [direc.main_directory, direc.run_name, direc.exp_nr{II}, ...
            files_yz{1}];
        yz_info = ncinfo(yz_loc);
        % Define variables to be saved
        yz_var = {yz_info.Variables.Name}';
        % Define time vatiable
        yz_time = ncread(yz_loc, 'time');
        % Preallocate variable sizes
        yz_dimy = yz_info.Variables(ismember({yz_info.Variables.Name},{'yt'})).Size(1);
        yz_nprocy = str2num(inp.namoptions.value(ismember(inp.namoptions.option, 'nprocy')));
        yz_dimz = yz_info.Variables(ismember({yz_info.Variables.Name},{'zt'})).Size(1);
        yz_nprocz = str2num(inp.namoptions.value(ismember(inp.namoptions.option, 'kmax')));
        yz_dimt = length(yz_time);
        
        % Fill the INFO field
        for i = 1 : size(yz_info.Variables,2)
            
            f_output.yz.INFO = [f_output.yz.INFO ; table({...
                yz_var{i}}, {...
                yz_info.Variables(i).Attributes(1).Value}, {...
                yz_info.Variables(i).Attributes(2).Value}, ...
                'VariableNames',{'Name','LongName','Unit'})];
        end
        f_output.yz.INFO = [f_output.yz.INFO ; ...
            table({'xt'}, {'Cross section x location'}, {'m'}, ...
                'VariableNames',{'Name','LongName','Unit'})];
        
        % Preallocation data
        for i = 1 : length(yz_var)
            if strcmp(yz_var{i},'time')
                
                f_output.yz.(yz_var{i}) = ncread(yz_loc, yz_var{i});
                
            elseif ismember(yz_var(i),{'zt','zm'})
                
                f_output.yz.(yz_var{i}) = double(ncread(yz_loc, yz_var{i}));
                
            elseif ismember(yz_var(i),{'yt','ym'})
                
                f_output.yz.(yz_var{i}) = NaN(yz_dimy*yz_nprocy, 1);     
                
            elseif ismember(yz_var(i), yz_loadvars)
                
                f_output.yz.(yz_var{i}) = NaN(yz_dimy*yz_nprocy, yz_dimz, yz_dimt);
                
            end
        end
        
        % Determine the number of files
        yz_files_n1 = length(files_yz);
        yz_files_n2 = yz_nprocy;
        if yz_files_n1 == yz_files_n2
            for i = 1 : yz_files_n1
                % Directory location
                yz_loc = [direc.main_directory, direc.run_name, direc.exp_nr{II}, ...
                    files_yz{i}];
                
                yz_locy = str2num(files_yz{i}(14:16));
                
                yz_filly = [1+yz_dimy*yz_locy : 1+yz_dimy*(yz_locy+1)-1];

                for j = 1 : length(yz_var)
                    
                    if ismember(yz_var(j),{'yt','ym'}) 
                        f_output.yz.(yz_var{j})(yz_filly) = ...
                            ncread(yz_loc, yz_var{j});
                    elseif ismember(yz_var(j), yz_loadvars)
                        f_output.yz.(yz_var{j})(yz_filly, :, :) = ...
                            ncread(yz_loc, yz_var{j});
                    end
                    
                end
            end
        end
        f_output.yz.xt = yz_x;
        
    end
    clear yz_* i j files_yz
    
    %% Rename the scalars in the output to match the names in scalar.inp
    
    % Get the input names of the scalars
    corr_inpscalars = inp.scalar.Properties.VariableNames(2:end);
    % Names of the three cross sectiond
    cross_fields = fields(f_output);
    
    % Loop over the three cross sections
    for II = 1 : length(cross_fields)

        for i = 1 : length(corr_inpscalars)

            
            % Create the name that will be used to search in f_output.(cross_fields).INFO.LongNames
            corr_name = '000';
            corr_name(end-(length(num2str(i))-1) : end) = num2str(i);
            
            corr_name = [{['sv' corr_name cross_fields{II}]},{['scalar ' corr_name]}]; 
            
            corr_IN = find(ismember(f_output.(cross_fields{II}).INFO.Name, corr_name(1)));
            
            f_output.(cross_fields{II}).INFO.Name{corr_IN} = ...
                [corr_inpscalars{i}, cross_fields{II}];
            
            corr_longname = f_output.(cross_fields{II}).INFO.LongName{corr_IN};
            for j = 1 : length(corr_longname) - length(corr_name{2}) + 1
                if strcmp(corr_name{2}, corr_longname(j : j+length(corr_name{2})-1))
                    
                    f_output.(cross_fields{II}).INFO.LongName{corr_IN} = ...
                        [corr_longname(1:j-1), ' scalar ', corr_inpscalars{i}, ...
                        corr_longname(j+length(corr_name{2}):end)];
                end
            end
            
            f_output.(cross_fields{II}).([corr_inpscalars{i}, cross_fields{II}]) = ...
                f_output.(cross_fields{II}).(corr_name{1});
            f_output.(cross_fields{II}) = rmfield(f_output.(cross_fields{II}), corr_name{1});
            
        end
    end
    
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%     READ NAMFIELDDUMP output      %%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function f_output = f__NAMFIELDDUMPread(direc, inp, addons, II)
    
    
        

    
    %%%%%
    % Setup the information on the FIELDDUMP output
    
    % Preallocation output structure
    f_output = struct();
    f_output.README = 'Original NetCDF files contain 4-D structure [x,y,z,t], which are too large';
    
    fd_files = dir([direc.main_directory, direc.run_name, direc.exp_nr{II}, ...
        'fielddump*']);
    fd_files = {fd_files.name}';
    % Load the first file as a dummy
    fd_fileloc = [direc.main_directory, direc.run_name, direc.exp_nr{II}, ...
        fd_files{1}];
    fd_info = ncinfo(fd_fileloc);
    % Define variables to be saved
    fd_var = {fd_info.Variables.Name}';
    % Define time vatiable
    fd_time = ncread(fd_fileloc, 'time');
    % Preallocate variable sizes
    fd_dimx = fd_info.Variables(ismember({fd_info.Variables.Name},{'xt'})).Size(1);
    fd_sizex = str2num(inp.namoptions.value(ismember(inp.namoptions.option, 'xsize')));
    fd_nprocx = str2num(inp.namoptions.value(ismember(inp.namoptions.option, 'nprocx')));
    fd_dimy = fd_info.Variables(ismember({fd_info.Variables.Name},{'yt'})).Size(1);
    fd_sizey = str2num(inp.namoptions.value(ismember(inp.namoptions.option, 'ysize')));
    fd_nprocy = str2num(inp.namoptions.value(ismember(inp.namoptions.option, 'nprocy')));
    fd_dimz = fd_info.Variables(ismember({fd_info.Variables.Name},{'zt'})).Size(1);
    fd_nprocz = str2num(inp.namoptions.value(ismember(inp.namoptions.option, 'kmax')));
    fd_dimt = length(fd_time);
    
    % Preallocate the info on the variables
    f_output.INFO = table({},{},{}, 'VariableNames',{'Name','LongName','Unit'});
    % Fill the INFO field
    for i = 1 : size(fd_info.Variables,2)

        f_output.INFO = [f_output.INFO ; table({...
            fd_var{i}}, {...
            fd_info.Variables(i).Attributes(1).Value}, {...
            fd_info.Variables(i).Attributes(2).Value}, ...
            'VariableNames',{'Name','LongName','Unit'})];
    end
    % Save the name of the directory
    f_output.direc = [direc.main_directory, direc.run_name, direc.exp_nr{II}];
    % Save an example file name
    f_output.ExampleName = fd_files{1};
    
    fd_IN = find(ismember(fd_var,{'time','xt','xm','yt','ym','zt','zm'}));
    
    % Preallocation data
    for i = 1 : length(fd_IN)
        if ismember(fd_var(fd_IN(i)),{'xt','xm'})
            temp = ncread(fd_fileloc, fd_var{fd_IN(i)});
            temp = double([temp(1) : temp(2)-temp(1): fd_sizex]');
            f_output.(fd_var{fd_IN(i)}) = temp;
        elseif ismember(fd_var(fd_IN(i)),{'yt','ym'})
            temp = ncread(fd_fileloc, fd_var{fd_IN(i)});
            temp = double([temp(1) : temp(2)-temp(1): fd_sizey]');
            f_output.(fd_var{fd_IN(i)}) = temp;
            
        else
            f_output.(fd_var{fd_IN(i)}) = double(ncread(fd_fileloc, fd_var{fd_IN(i)}));
        end
    end
    
    f_output.x_proc.filenames = cell(fd_nprocx,1);
    f_output.x_proc.xt = NaN(fd_nprocx,fd_dimx);
    f_output.x_proc.xm = NaN(fd_nprocx,fd_dimx);
    for i = 1 : fd_nprocx
        fd_nx = '000';
        fd_nx(end-(length(num2str(i-1))-1) : end) = num2str(i-1);
        f_output.x_proc.filenames{i} = fd_nx;
        
        f_output.x_proc.xt(i,:) = f_output.xt([1+((i-1)*fd_dimx) : (i*fd_dimx)]);
        f_output.x_proc.xm(i,:) = f_output.xm([1+((i-1)*fd_dimx) : (i*fd_dimx)]);
    end
    f_output.y_proc.filenames = cell(fd_nprocy,1);
    f_output.y_proc.yt = NaN(fd_nprocy,fd_dimy);
    f_output.y_proc.ym = NaN(fd_nprocy,fd_dimy);
    for i = 1 : fd_nprocy
        fd_ny = '000';
        fd_ny(end-(length(num2str(i-1))-1) : end) = num2str(i-1);
        f_output.y_proc.filenames{i} = fd_ny;
        
        f_output.y_proc.yt(i,:) = f_output.yt([1+((i-1)*fd_dimy) : (i*fd_dimy)]);
        f_output.y_proc.ym(i,:) = f_output.ym([1+((i-1)*fd_dimy) : (i*fd_dimy)]);
    end
    
    %
    %%%%% Change the scalar names
    %
    
    % Get the input names of the scalars
    corr_inpscalars = inp.scalar.Properties.VariableNames(2:end);
    

    for i = 1 : length(corr_inpscalars)

        % Create the name that will be used to search in f_output.(cross_fields).INFO.LongNames
        corr_name = '000';
        corr_name(end-(length(num2str(i))-1) : end) = num2str(i);

        corr_name = [{['sv' corr_name]},{['Scalar ' corr_name]}]; 

        corr_IN = find(ismember(f_output.INFO.Name, corr_name(1)));
        
        if isempty(corr_IN) == 0
            f_output.INFO.Name{corr_IN} = corr_inpscalars{i};

            corr_longname = f_output.INFO.LongName{corr_IN};
            for j = 1 : length(corr_longname) - length(corr_name{2}) + 1
                if strcmp(corr_name{2}, corr_longname(j : j+length(corr_name{2})-1))
                    if j == 1
                        f_output.INFO.LongName{corr_IN} = ...
                            ['Scalar ', corr_inpscalars{i}, ...
                            corr_longname(j+length(corr_name{2}):end)];
                    else
                        f_output.INFO.LongName{corr_IN} = ...
                            [corr_longname(1:j-1), ' scalar ', corr_inpscalars{i}, ...
                            corr_longname(j+length(corr_name{2}):end)];
                    end
                end
            end
        end
    end
    
    
    
end


        