%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%         Ruben Schulte         %%%%%
%%%%%    DALES modify input files   %%%%%
%%%%%      started: 08-01-2021      %%%%%
%%%%%      changed: 20-04-2021      %%%%%
%%%%%       final: 23-05-2025       %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all
close all
clc

%% General settings

% Define the relevant (data) directories 
direc = struct();
direc.main_directory = 'C:\Users\schul068\PhD project\DALES\Runs\';
direc.run_name = 'nh3plume';
direc.exp_nr = '801';



check = struct();

check.makefigures = 'no';
%%%%% Resolution
% 50x50 m horizontal grid
check.dz_min        = 5;    % highest vertical resolution [m]  (near surface)
check.dz_max        = 5;   % lowest vertical resolution [m] (higher atmosphere)
check.dz_trans      = 0;    % stepsize of resolution transition [m] (transition between high to low vert. resulution)
check.dz_trans_h    = 100;   % starting height of vertical resulution transition [m]
check.z_max         = 3000; % Maximum height of the domain

%%%%% Scalars
%
% Will the scalars be redefined by this script? yes/no
check.redefine_scalars = 'yes';
% Define the scalars as { 'name', [begin,end height], [begin,end concentration (ppb)]
%       Note that the script will interpolate the concentration between the
%       defined begin and end height
check.sv_vertprof   = {'nh3_r1b', [0,1400 ; 1400,3000], [13.1, 13.1; 1,1] ...     % Definition scalar #01 (reference)
                     ; 'nh3_r1p', [0,1400 ; 1400,3000], [ 0.0,  0.0; 0,0] ...     % Definition scalar #02 (reference)
                     ; 'nh3_b1b', [0,1400 ; 1400,3000], [6.5, 6.5; 1,1] ...     % Definition scalar #03 (background NH3 1)
                     ; 'nh3_b1p', [0,1400 ; 1400,3000], [ 0.0,  0.0; 0,0] ...     % Definition scalar #04 (background NH3 1)
                     ; 'nh3_b2b', [0,1400 ; 1400,3000], [19.7, 19.7; 1,1] ...     % Definition scalar #05 (background NH3 2)
                     ; 'nh3_b2p', [0,1400 ; 1400,3000], [ 0.0,  0.0; 0,0] ...     % Definition scalar #06 (background NH3 2)
                     ; 'nh3_b3b', [0,1400 ; 1400,3000], [33, 33; 1,1] ...     % Definition scalar #07 (background NH3 3)
                     ; 'nh3_b3p', [0,1400 ; 1400,3000], [ 0.0,  0.0; 0,0] ...     % Definition scalar #08 (background NH3 3)
                       };

% Non-periodic Boundary Condition (BC) settings in namoptions
check.sv_periodBCs = [ 2 ...               % Scalar indices with non-periodic BCs for the reference runs
                     ; 4 ...
                     ; 6 ...
                     ; 8 ...
                      ];  
                  
% Percentage-chemistry functionality settings in namoptions
check.sv_chemrate   = [ 1, -0.05 ...      % Scalar index (#01) with chemistry rate (loss/source) [%] (reference)
                      ; 2, -0.05 ...      % Scalar index (#02) with chemistry rate (loss/source) [%] (reference)
                      ; 3, -0.05 ...      % Scalar index (#01) with chemistry rate (loss/source) [%] (background NH3 1)
                      ; 4, -0.05 ...      % Scalar index (#02) with chemistry rate (loss/source) [%] (background NH3 1)
                      ; 5, -0.05 ...      % Scalar index (#03) with chemistry rate (loss/source) [%] (background NH3 2)
                      ; 6, -0.05 ...      % Scalar index (#04) with chemistry rate (loss/source) [%] (background NH3 2)
                      ; 7, -0.05 ...      % Scalar index (#05) with chemistry rate (loss/source) [%] (background NH3 3)
                      ; 8, -0.05 ...      % Scalar index (#06) with chemistry rate (loss/source) [%] (background NH3 3)
                       ];
                   
% Split-flux functionality settings  in namoptions
check.sv_splitflux  = [        ...      
                     ;  1,  1, 1 ...      % Scalar index (#01), set nr. & set order (reference)
                     ;  2,  1, 2 ...      % Scalar index (#02), set nr. & set order (reference)
                     ;  3,  2, 1 ...      % Scalar index (#03), set nr. & set order (background NH3 1)
                     ;  4,  2, 2 ...      % Scalar index (#04), set nr. & set order (background NH3 1)
                     ;  5,  3, 1 ...      % Scalar index (#05), set nr. & set order (background NH3 2)
                     ;  6,  3, 2 ...      % Scalar index (#06), set nr. & set order (background NH3 2)
                     ;  7,  4, 1 ...      % Scalar index (#05), set nr. & set order (background NH3 3)
                     ;  8,  4, 2 ...      % Scalar index (#06), set nr. & set order (background NH3 3)
                     ];

% Surface.interactive.inp settings
% --> Spin-up run
check.sv_flux       = [0,   0; ...      % For each typenr.: the surface flux for scalar #01 (reference)
                       0,   0; ...      % For each typenr.: the surface flux for scalar #02 (reference)
                       0,   0; ...      % For each typenr.: the surface flux for scalar #03 (background NH3 1)
                       0,   0; ...      % For each typenr.: the surface flux for scalar #04 (background NH3 1)
                       0,   0; ...      % For each typenr.: the surface flux for scalar #05 (background NH3 2)
                       0,   0; ...      % For each typenr.: the surface flux for scalar #06 (background NH3 2)
                       0,   0; ...      % For each typenr.: the surface flux for scalar #07 (background NH3 3)
                       0,   0; ...      % For each typenr.: the surface flux for scalar #08 (background NH3 3)
                       ];

% --> Analysis run
% check.sv_flux       = [-0.045,   0; ...      % For each typenr.: the surface flux for scalar #01 (background NH3 1)
%                        -0.045,  45; ...      % For each typenr.: the surface flux for scalar #02 (background NH3 1)
%                        -0.045,   0; ...      % For each typenr.: the surface flux for scalar #01 (background NH3 1)
%                        -0.045,  45; ...      % For each typenr.: the surface flux for scalar #02 (background NH3 1)
%                        -0.045,  0; ...      % For each typenr.: the surface flux for scalar #02 (background NH3 2)
%                        -0.045,  45; ...      % For each typenr.: the surface flux for scalar #02 (background NH3 2)
%                        -0.045,  0; ...      % For each typenr.: the surface flux for scalar #02 (background NH3 3)
%                        -0.045,  45; ...      % For each typenr.: the surface flux for scalar #02 (background NH3 3)
%                        ];

%% End of settings



% !!!!!!!!!!!!!!!!!!!! Do not change the code below !!!!!!!!!!!!!!!!!!!!

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    fprintf('Start the code: \n')  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


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

% Make a new directory for the old input files
name_list = dir([ direc.main_directory direc.run_name direc.exp_nr ]);
name_list = {name_list.name}';
olddir_name = 'input_old_v000';
dir_counter = 0;
while dir_counter == 0 
    if ismember({olddir_name}, name_list) == 0
        dir_counter = 1;
        mkdir([ direc.main_directory direc.run_name direc.exp_nr ], olddir_name);
    else
        olddir_name_nr = str2num(olddir_name(end-2:end)) +1;
        if length(num2str(olddir_name_nr)) < 4
            olddir_name = [olddir_name( 1:end-length(num2str(olddir_name_nr)) ) num2str(olddir_name_nr)];
        end
    end
end
% Create the input_old folder
% mkdir([ direc.main_directory direc.run_name direc.exp_nr ], olddir_name);

clear fields_* i


%% Make sure that the vertical resolution is not too small
    
if mod(check.dz_min/2, 0.01) ~= 0

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    fprintf('\n!!!!! ERROR !!!!!\n')  
    fprintf('- Your vertical resolution is so small that interpolation at 1 cm resolution is not enough\n\n')      
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    return
end



%% Load & adapt the baseprof input file


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    fprintf('\n- Update baseprof.inp\n')  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% Define the name of the file to open
bp_name = [ direc.main_directory direc.run_name direc.exp_nr 'baseprof.inp.' direc.exp_nr(1:end-1) ]; 
% Open the file
bp_fid = fopen(bp_name,'r');

% read the header (on the second line = 2x)
bp_firstline = fgetl(bp_fid);
bp_secondline = fgetl(bp_fid);
% Close the file
fclose(bp_fid);

% Save the header data separated by spaces and remove "#"
bp_INlinestarter = find(ismember(bp_secondline, '#'));
bp_secondline(1:bp_INlinestarter) = '';
bp_INheader = find(ismember(bp_secondline, ' ') | ismember(bp_secondline, '	') | ismember(bp_secondline, ','));
bp_header = cell(length(bp_INheader)+1,1);
for i = 1 : length(bp_INheader)+1
    if i == 1
        bp_header{i} = bp_secondline(1 : bp_INheader(i)-1);
    elseif i<= length(bp_INheader)
        bp_header{i} = bp_secondline(bp_INheader(i-1)+1 : bp_INheader(i)-1);
    else
        bp_header{i} = bp_secondline(bp_INheader(i-1)+1 : end);
    end
end

% Load the data in a matrix
bp_data = dlmread(bp_name, '',2,0);
% Calculate the vertical step sizes
bp_data_delta = NaN(size(bp_data));
for i = 1 : size(bp_data,1)-1
    bp_data_delta(i,:) = bp_data(i+1, :) - bp_data(i, :);
end





% Define the vertical domain of the model
bp_zdomain = [bp_data(1,1) - (bp_data(2,1)-bp_data(1,1))/2, ...
    bp_data(end,1) + (bp_data(2,1)-bp_data(1,1))/2];



% Define the old vertical resulution
zold = bp_data(:,1);
% Define first value of new zc vector
znew = [bp_zdomain(1) + check.dz_min/2 : ...
           check.dz_min : check.dz_trans_h] ;
% 
if znew(end) <  check.dz_trans_h
    znew = [znew, znew(end) + check.dz_min];
end
bp_tempdz = check.dz_min + check.dz_trans;
while bp_tempdz(end) < check.dz_max
    bp_tempdz = [bp_tempdz ; bp_tempdz(end)+check.dz_trans];
end
for i = 1 : length(bp_tempdz)
    znew = [znew, znew(end)+bp_tempdz(i)];
end
while znew(end) < bp_zdomain(2)
    znew = [znew, znew(end)+check.dz_max];
end
if znew(end) > bp_zdomain(2)
    znew = znew(1:end-1);
end
znew = znew';





%%%%%
% Interpolate & extrapolate the original data to a 1 cm vertical resolution
%%%%%
bp_int_data = NaN(length([bp_zdomain(1): 0.01 : bp_zdomain(2)]), length(bp_header));
bp_int_data(:,1) = [bp_zdomain(1): 0.01 : bp_zdomain(2)]';
for i = 1 : length(zold)-1

    bp_INz = find(bp_int_data(:,1) >= bp_data(i,1) & ...
                    bp_int_data(:,1) <= bp_data(i+1,1));

    for j = 2 : length(bp_header)

        if bp_data(i+1,j) == bp_data(i,j)
            bp_temp = bp_data(i,j) .* ones(length(bp_INz),1);
        else
            bp_temp = [bp_data(i,j) : ...
                      ( bp_data(i+1,j) - bp_data(i,j)) / (length(bp_INz)-1) : ...
                       bp_data(i+1,j)]';
        end
        bp_int_data(bp_INz,j) = bp_temp; 



        % Define data between 0m and z1
        if i == 1
            if  bp_data(i+1,j) == bp_data(i,j)
                bp_temp_sfc = bp_data(i,j) .* ones(length([1 : bp_INz(1)]),1);
            else
                bp_temp_sfc = ...
                    [bp_data(i,j) - (bp_temp(2) - bp_temp(1)) * (length([1 : bp_INz(1)])-1) : ...
                    bp_temp(2) - bp_temp(1) : ...
                    bp_data(i,j)]';
            end

            bp_int_data(1 : bp_INz(1), j) = bp_temp_sfc;

        % Define data between zend and the top of the domain
        elseif i == length(zold)-1 && bp_data(i+1,1) < bp_int_data(end,1)
            if  bp_data(i+1,j) == bp_data(i,j)

                bp_temp_top = bp_data(i,j) .* ones(length([bp_INz(end): size(bp_int_data,1)]),1);
            else

                bp_temp_top = ...
                    [ bp_data(i+1,j) : ...
                    bp_temp(2) - bp_temp(1) : ...
                     bp_data(i+1,j) + (bp_temp(2) - bp_temp(1)) * (length([bp_INz(end): size(bp_int_data,1)])-1)]';
            end
            bp_int_data(bp_INz(end) : size(bp_int_data,1), j) = bp_temp_top;

        end     % if i == 1

    end     % for j = 2 : length(bp_header)

end     % for i = 1 : length(zold)-1



bp_data_new = NaN(length(znew), length(bp_header));
bp_data_new(:,1) = znew;
for i = 1 : length(znew)

    for j = 2 : length(bp_header)

        bp_data_new(i,j) = bp_int_data( bp_int_data(:,1) == bp_data_new(i,1) ,j);

    end
end


if strcmp(check.makefigures, 'yes')
    for i = 2 : length(bp_header)
        figure
        scatter(bp_data(:,i), bp_data(:,1), 'xk')
        hold on
        scatter(bp_data_new(:,i), bp_data_new(:,1), '+r')
        plot(bp_int_data(:,i), bp_int_data(:,1))
        xlabel(bp_header{i})
        ylabel([bp_header{1} ' [m]'])
    end
end


%%%%%%%%%%
% Move original prof.inp to the input_old folder
%%%%%%%%%%
% Move the prof.inp to the new folder
movefile(bp_name, [ direc.main_directory direc.run_name direc.exp_nr olddir_name]);

% Create the new prof.inp file
bp_fid = fopen(bp_name, 'w');
% Write the first line in the prof.inp file
fprintf(bp_fid,'%s', bp_firstline);
% Add/change the second line (header) in the prof.inp file
for i = 1 : length(bp_header)
    if i == 1 
        fprintf(bp_fid,'\n%s', ['#' bp_header{i}]);
    else
        fprintf(bp_fid,'\t%s', bp_header{i});
    end
end
% Fill the prof.inp file
for i = 1 : size(bp_data_new,1)
    for j = 1 : size(bp_data_new,2)
        if j == 1
            fprintf(bp_fid,'\n%.2f', bp_data_new(i,j));
        else
            fprintf(bp_fid,'\t%.16f', bp_data_new(i,j));
        end
    end
end
% Close the file
fclose(bp_fid);
    
clear bp_*


%% Load & adapt the prof input file


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    fprintf('\n- Update prof.inp\n')  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    
if mod(check.dz_min/2, 0.01) ~= 0

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    fprintf('\n!!!!! ERROR !!!!!\n')  
    fprintf('- Your vertical resolution is so small that interpolation at 1 cm resolution is not enough\n\n')      
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    return
end

% Define the name of the file to open
pr_name = [ direc.main_directory direc.run_name direc.exp_nr 'prof.inp.' direc.exp_nr(1:end-1) ]; 
% Open the file
pr_fid = fopen(pr_name,'r');

% read the header (on the second line = 2x)
pr_firstline = fgetl(pr_fid);
pr_secondline = fgetl(pr_fid);
% Close the file
fclose(pr_fid);

% Save the header data separated by spaces and remove "#"
pr_secondline(ismember(pr_secondline, '#')) = '';
pr_INheader = find(ismember(pr_secondline, ' ') | ismember(pr_secondline, '	') | ismember(pr_secondline, ','));
pr_header = cell(length(pr_INheader)+1,1);
for i = 1 : length(pr_INheader)+1
    if i == 1
        pr_header{i} = pr_secondline(1 : pr_INheader(i)-1);
    elseif i<= length(pr_INheader)
        pr_header{i} = pr_secondline(pr_INheader(i-1)+1 : pr_INheader(i)-1);
    else
        pr_header{i} = pr_secondline(pr_INheader(i-1)+1 : end);
    end
end

% Load the data in a matrix
pr_data = dlmread(pr_name, '',2,0);
% Calculate the vertical step sizes
pr_data_delta = NaN(size(pr_data));
for i = 1 : size(pr_data,1)-1
    pr_data_delta(i,:) = pr_data(i+1, :) - pr_data(i, :);
end





% Define the vertical domain of the model
pr_zdomain = [pr_data(1,1) - (pr_data(2,1)-pr_data(1,1))/2, ...
    pr_data(end,1) + (pr_data(2,1)-pr_data(1,1))/2];


%%%%%
% Interpolate & extrapolate the original data to a 1 cm vertical resolution
%%%%%
pr_int_data = NaN(length([pr_zdomain(1): 0.01 : pr_zdomain(2)]), length(pr_header));
pr_int_data(:,1) = [pr_zdomain(1): 0.01 : pr_zdomain(2)]';
for i = 1 : length(zold)-1

    pr_INz = find(pr_int_data(:,1) >= pr_data(i,1) & ...
                    pr_int_data(:,1) <= pr_data(i+1,1));

    for j = 2 : length(pr_header)

        if pr_data(i+1,j) == pr_data(i,j)
            pr_temp = pr_data(i,j) .* ones(length(pr_INz),1);
        else
            pr_temp = [pr_data(i,j) : ...
                      ( pr_data(i+1,j) - pr_data(i,j)) / (length(pr_INz)-1) : ...
                       pr_data(i+1,j)]';
        end
        pr_int_data(pr_INz,j) = pr_temp; 



        % Define data between 0m and z1
        if i == 1
            if  pr_data(i+1,j) == pr_data(i,j)
                pr_temp_sfc = pr_data(i,j) .* ones(length([1 : pr_INz(1)]),1);
            else
                pr_temp_sfc = ...
                    [pr_data(i,j) - (pr_temp(2) - pr_temp(1)) * (length([1 : pr_INz(1)])-1) : ...
                    pr_temp(2) - pr_temp(1) : ...
                    pr_data(i,j)]';
            end

            pr_int_data(1 : pr_INz(1), j) = pr_temp_sfc;

        % Define data between zend and the top of the domain
        elseif i == length(zold)-1 && pr_data(i+1,1) < pr_int_data(end,1)
            if  pr_data(i+1,j) == pr_data(i,j)

                pr_temp_top = pr_data(i,j) .* ones(length([pr_INz(end): size(pr_int_data,1)]),1);
            else

                pr_temp_top = ...
                    [ pr_data(i+1,j) : ...
                    pr_temp(2) - pr_temp(1) : ...
                     pr_data(i+1,j) + (pr_temp(2) - pr_temp(1)) * (length([pr_INz(end): size(pr_int_data,1)])-1)]';
            end
            pr_int_data(pr_INz(end) : size(pr_int_data,1), j) = pr_temp_top;

        end     % if i == 1

    end     % for j = 2 : length(pf_header)

end     % for i = 1 : length(zold)-1



pr_data_new = NaN(length(znew), length(pr_header));
pr_data_new(:,1) = znew;
for i = 1 : length(znew)

    for j = 2 : length(pr_header)

        pr_data_new(i,j) = pr_int_data( pr_int_data(:,1) == pr_data_new(i,1) ,j);

    end
end


if strcmp(check.makefigures, 'yes')
    for i = 2 : length(pr_header)
        figure
        scatter(pr_data(:,i), pr_data(:,1), 'xk')
        hold on
        scatter(pr_data_new(:,i), pr_data_new(:,1), '+r')
        plot(pr_int_data(:,i), pr_int_data(:,1))
        xlabel(pr_header{i})
        ylabel([pr_header{1} ' [m]'])
    end
end


%%%%%%%%%%
% Move original prof.inp to the input_old folder
%%%%%%%%%%
% Move the prof.inp to the new folder
movefile(pr_name, [ direc.main_directory direc.run_name direc.exp_nr olddir_name]);

% Create the new prof.inp file
pr_fid = fopen(pr_name, 'w');
% Change the first line of the prof.inp file so it matches the new z array
pr_firstline_new = pr_firstline;
pr_temp = 'dz = ';
for i = 1 : length(pr_firstline) - length(pr_temp)
    if strcmp(pr_firstline(i : i+length(pr_temp)-1), pr_temp)
        pr_firstline_new = [pr_firstline(1:i-1) pr_temp num2str(znew(2) - znew(1))];
    end
end
pr_temp = ['Nz=' num2str(length(zold))];
for i = 1: length(pr_firstline_new) - length(pr_temp)
    if strcmp(pr_firstline_new(i : i+length(pr_temp)-1), pr_temp)
        pr_firstline_new = [pr_firstline_new(1:i-1) 'Nz= ' num2str(length(znew)) pr_firstline_new(i+length(pr_temp): end)];
    end
end
% Change the first line in the prof.inp file
fprintf(pr_fid,'%s', pr_firstline_new);
% Add/change the second line (header) in the prof.inp file
for i = 1 : length(pr_header)
    if i == 1 
        fprintf(pr_fid,'\n%s', ['#' pr_header{i}]);
    else
        fprintf(pr_fid,'\t%s', pr_header{i});
    end
end
% Fill the prof.inp file
for i = 1 : size(pr_data_new,1)
    for j = 1 : size(pr_data_new,2)
        if j == 1
            fprintf(pr_fid,'\n%.2f', pr_data_new(i,j));
        elseif strcmp(pr_header{j},'th')
            fprintf(pr_fid,'\t%.9f', pr_data_new(i,j));
        elseif strcmp(pr_header{j},'qt')
            fprintf(pr_fid,'\t%.14f', pr_data_new(i,j));
        else
            fprintf(pr_fid,'\t%f', pr_data_new(i,j));
        end
    end
end
% Close the file
fclose(pr_fid);

clear pr_*
    

%% Load & adapt the lscale input file


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    fprintf('\n- Update lscale.inp\n')  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
if mod(check.dz_min/2, 0.01) ~= 0

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    fprintf('\n!!!!! ERROR !!!!!\n')  
    fprintf('- Your vertical resolution is so small that interpolation at 1 cm resolution is not enough\n\n')      
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    return
end

% Define the name of the file to open
ls_name = [ direc.main_directory direc.run_name direc.exp_nr 'lscale.inp.' direc.exp_nr(1:end-1) ]; 
% Open the file
ls_fid = fopen(ls_name,'r');

% read the header (on the second line = 2x)
ls_firstline = fgetl(ls_fid);
ls_secondline = fgetl(ls_fid);
% Close the file
fclose(ls_fid);

% Save the header data separated by spaces and remove "#"
ls_secondline(ismember(ls_secondline, '#')) = '';
ls_INheader = find(ismember(ls_secondline, ' ') | ismember(ls_secondline, '	') | ismember(ls_secondline, ','));
ls_header = cell(length(ls_INheader)+1,1);
for i = 1 : length(ls_INheader)+1
    if i == 1
        ls_header{i} = ls_secondline(1 : ls_INheader(i)-1);
    elseif i<= length(ls_INheader)
        ls_header{i} = ls_secondline(ls_INheader(i-1)+1 : ls_INheader(i)-1);
    else
        ls_header{i} = ls_secondline(ls_INheader(i-1)+1 : end);
    end
end

% Load the data in a matrix
ls_data = dlmread(ls_name, '',2,0);
% Calculate the vertical step sizes
ls_data_delta = NaN(size(ls_data));
for i = 1 : size(ls_data,1)-1
    ls_data_delta(i,:) = ls_data(i+1, :) - ls_data(i, :);
end





% Define the vertical domain of the model
ls_zdomain = [ls_data(1,1) - (ls_data(2,1)-ls_data(1,1))/2, ...
    ls_data(end,1) + (ls_data(2,1)-ls_data(1,1))/2];




%%%%%
% Interpolate & extrapolate the original data to a 1 cm vertical resolution
%%%%%
ls_int_data = NaN(length([ls_zdomain(1): 0.01 : ls_zdomain(2)]), length(ls_header));
ls_int_data(:,1) = [ls_zdomain(1): 0.01 : ls_zdomain(2)]';
for i = 1 : length(zold)-1

    ls_INz = find(ls_int_data(:,1) >= ls_data(i,1) & ...
                    ls_int_data(:,1) <= ls_data(i+1,1));

    for j = 2 : length(ls_header)

        if ls_data(i+1,j) == ls_data(i,j)
            ls_temp = ls_data(i,j) .* ones(length(ls_INz),1);
        else
            ls_temp = [ls_data(i,j) : ...
                      ( ls_data(i+1,j) - ls_data(i,j)) / (length(ls_INz)-1) : ...
                       ls_data(i+1,j)]';
        end
        ls_int_data(ls_INz,j) = ls_temp; 



        % Define data between 0m and z1
        if i == 1
            if  ls_data(i+1,j) == ls_data(i,j)
                ls_temp_sfc = ls_data(i,j) .* ones(length([1 : ls_INz(1)]),1);
            else
                ls_temp_sfc = ...
                    [ls_data(i,j) - (ls_temp(2) - ls_temp(1)) * (length([1 : ls_INz(1)])-1) : ...
                    ls_temp(2) - ls_temp(1) : ...
                    ls_data(i,j)]';
            end

            ls_int_data(1 : ls_INz(1), j) = ls_temp_sfc;

        % Define data between zend and the top of the domain
        elseif i == length(zold)-1 && ls_data(i+1,1) < ls_int_data(end,1)
            if  ls_data(i+1,j) == ls_data(i,j)

                ls_temp_top = ls_data(i,j) .* ones(length([ls_INz(end): size(ls_int_data,1)]),1);
            else

                ls_temp_top = ...
                    [ ls_data(i+1,j) : ...
                    ls_temp(2) - ls_temp(1) : ...
                     ls_data(i+1,j) + (ls_temp(2) - ls_temp(1)) * (length([ls_INz(end): size(ls_int_data,1)])-1)]';
            end
            ls_int_data(ls_INz(end) : size(ls_int_data,1), j) = ls_temp_top;

        end     % if i == 1

    end     % for j = 2 : length(pf_header)

end     % for i = 1 : length(zold)-1



ls_data_new = NaN(length(znew), length(ls_header));
ls_data_new(:,1) = znew;
for i = 1 : length(znew)

    for j = 2 : length(ls_header)

        ls_data_new(i,j) = ls_int_data( ls_int_data(:,1) == ls_data_new(i,1) ,j);

    end
end



if strcmp(check.makefigures, 'yes')
    for i = 2 : length(ls_header)
        figure
        scatter(ls_data(:,i), ls_data(:,1), 'xk')
        hold on
        scatter(ls_data_new(:,i), ls_data_new(:,1), '+r')
        plot(ls_int_data(:,i), ls_int_data(:,1))
        xlabel(ls_header{i})
        ylabel([ls_header{1} ' [m]'])
    end
end


%%%%%%%%%%
% Move original lscale.inp to the input_old folder
%%%%%%%%%%
% Move the lscale.inp to the new folder
movefile(ls_name, [ direc.main_directory direc.run_name direc.exp_nr olddir_name]);

% Create the new lscale.inp file
ls_fid = fopen(ls_name, 'w');
% Change the first line of the lscale.inp file so it matches the new z array
ls_firstline_new = ls_firstline;
ls_temp = 'dz =';
for i = 1 : length(ls_firstline) - length(ls_temp)
    if strcmp(ls_firstline(i : i+length(ls_temp)-1), ls_temp)
        ls_firstline_new = [ls_firstline(1:i-1) ls_temp ' ' num2str(znew(2) - znew(1))];
    end
end
ls_temp1 = 'Nz=';
ls_temp2 = num2str(length(zold)); 
ls_temp_IN = [];
for i = 1: length(ls_firstline_new) - length(ls_temp1)
    if strcmp(ls_firstline_new(i : i+length(ls_temp1)-1), ls_temp1)
        ls_temp_IN = [ls_temp_IN; i];
    end
end
for i = 1: length(ls_firstline_new) - length(ls_temp2)
    if strcmp(ls_firstline_new(i : i+length(ls_temp2)-1), ls_temp2)
        ls_temp_IN = [ls_temp_IN; i];
    end
end
ls_firstline_new = [ls_firstline_new(1:ls_temp_IN(1)-1) 'Nz= ' num2str(length(znew)) ...
    ls_firstline_new(ls_temp_IN(2)+length(ls_temp2): end)];

% Change the first line in the lscale.inp file
fprintf(ls_fid,'%s', ls_firstline_new);
% Add/change the second line (header) in the lscale.inp file
for i = 1 : length(ls_header)
    if i == 1 
        fprintf(ls_fid,'\n%s', ['#' ls_header{i}]);
    else
        fprintf(ls_fid,'\t%s', ls_header{i});
    end
end
% Fill the lscale.inp file
for i = 1 : size(ls_data_new,1)
    for j = 1 : size(ls_data_new,2)
        if j == 1
            fprintf(ls_fid,'\n%.2f', ls_data_new(i,j));
        else
            fprintf(ls_fid,'\t%f', ls_data_new(i,j));
        end
    end
end
% Close the file
fclose(ls_fid);

clear ls_*


%% Load & adapt the scalar input file


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    fprintf('\n- Update scalar.inp\n')  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


if mod(check.dz_min/2, 0.01) ~= 0

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    fprintf('\n!!!!! ERROR !!!!!\n')  
    fprintf('- Your vertical resolution is so small that interpolation at 1 cm resolution is not enough\n\n')      
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    return
end

% Define the name of the file to open
sc_name = [ direc.main_directory direc.run_name direc.exp_nr 'scalar.inp.' direc.exp_nr(1:end-1) ]; 
% Open the file
sc_fid = fopen(sc_name,'r');

% read the header (on the second line = 2x)
sc_firstline = fgetl(sc_fid);
sc_secondline = fgetl(sc_fid);
% Close the file
fclose(sc_fid);



if strcmp(check.redefine_scalars, 'yes')
    % Define the new header line
    sc_header = [{'zc'}; check.sv_vertprof(:,1) ];

    % Define the vertical domain of the model
    sc_zdomain = [check.dz_min/2, check.z_max];


    % Define the old vertical resulution
    sc_data = dlmread(sc_name, '',2,0);
    % Define new input data
    sc_data_new = NaN(length(znew), length(sc_header));
    % Fill the zc vector
    sc_data_new(:,1) = znew;
    for i = 2 : length(sc_header)
        sc_def_z = check.sv_vertprof{i-1,2};     % The predefine height
        sc_def_c = check.sv_vertprof{i-1,3};     % The predefine concentration values
        for j = 1 : size(sc_def_z,1)

            % Find what z indices are within the predefines heights
            sc_INfill = find(sc_data_new(:,1) > sc_def_z(j,1) & ...
                             sc_data_new(:,1) <= sc_def_z(j,2));

            % Define the concentration data to fill the matrix with
            %       Note that it is linearly interpolated between the
            %       first and second predefined concentration
            if (sc_def_c(j,2)-sc_def_c(j,1)) == 0
                sc_sv_c = sc_def_c(j,1) .* ones(length(sc_INfill),1);
            else
                sc_sv_c = [sc_def_c(j,1) : ...  % starting point
                           (sc_def_c(j,2)-sc_def_c(j,1)) / (length(sc_INfill)-1) : ...  % Step size
                           sc_def_c(j,2)]';     % End point
            end

            % Fill the data matrix
            sc_data_new(sc_INfill,i) = sc_sv_c;
        end            
    end

    if strcmp(check.makefigures, 'yes')
        for i = 2 : length(sc_header)
            figure
            scatter(sc_data_new(:,i), sc_data_new(:,1), '+r')
            xlabel(sc_header{i})
            ylabel([sc_header{1} ' [m]'])
        end
    end

    %%%%%%%%%%
    % Move original scalar.inp to the input_old folder
    %%%%%%%%%%
    % Move the scalar.inp to the new folder
    movefile(sc_name, [ direc.main_directory direc.run_name direc.exp_nr olddir_name]);


    % Create the new scalar.inp file
    sc_fid = fopen(sc_name, 'w');
    % Change the first line of the scalar.inp file so it matches the new z array
%     sc_INfirstline = [];
%     for i = 2 : length(sc_firstline)
%         if strcmp(sc_firstline(i-1:i), 'Nz')
%             sc_INfirstline = [sc_INfirstline , i+2];
%         elseif strcmp(sc_firstline(i-1:i), 'dz')
%             sc_INfirstline = [sc_INfirstline , i-3];
%         end
%     end
%     sc_firstline_nr = sc_firstline(sc_INfirstline(1): sc_INfirstline(2));
%     sc_firstline_nr(ismember(sc_firstline_nr, ' ') | ismember(sc_firstline_nr, '	')  ) = '';
%     for i = 1 : length(sc_firstline) - length(sc_firstline_nr)
%         if strcmp(sc_firstline(i : i+length(sc_firstline_nr)-1), sc_firstline_nr)
%             sc_firstline_new = [sc_firstline(1:i-1) ...
%                                 num2str(length(znew)) ...
%                                 sc_firstline(i+length(sc_firstline_nr):end)];
%         end
%     end
%     sc_temp  = 'dz =	';
%     for i = 1 : length(sc_firstline) - length(sc_temp)
%         if strcmp(sc_firstline(i : i+length(sc_temp)-1), sc_temp)
%             if (znew(2) - znew(1)) == (znew(end) - znew(end-1))
%                 sc_firstline_new = [sc_firstline(1:i-1) sc_temp num2str(znew(2) - znew(1))];
%             else
%                 sc_firstline_new = [sc_firstline(1:i-1) sc_temp num2str(znew(2) - znew(1)) ...
%                     ' - ' num2str(znew(end) - znew(end-1))];
%             end
%         end
%     end
    sc_firstline_new = sc_firstline;
    sc_temp = 'dz =';
    for i = 1 : length(sc_firstline) - length(sc_temp)
        if strcmp(sc_firstline(i : i+length(sc_temp)-1), sc_temp)
            sc_firstline_new = [sc_firstline(1:i-1) sc_temp ' ' num2str(znew(2) - znew(1))];
        end
    end
    sc_temp1 = 'Nz=';
    sc_temp2 = num2str(length(zold)); 
    sc_temp_IN = [];
    for i = 1: length(sc_firstline_new) - length(sc_temp1)
        if strcmp(sc_firstline_new(i : i+length(sc_temp1)-1), sc_temp1)
            sc_temp_IN = [sc_temp_IN; i];
        end
    end
    for i = 1: length(sc_firstline_new) - length(sc_temp2)
        if strcmp(sc_firstline_new(i : i+length(sc_temp2)-1), sc_temp2)
            sc_temp_IN = [sc_temp_IN; i];
        end
    end
    sc_firstline_new = [sc_firstline_new(1:sc_temp_IN(1)-1) 'Nz= ' num2str(length(znew)) ...
        sc_firstline_new(sc_temp_IN(2)+length(sc_temp2): end)];
    % Change the first line in the scalar.inp file
    fprintf(sc_fid,'%s', sc_firstline_new);
    % Add/change the second line (header) in the scalar.inp file
    for i = 1 : length(sc_header)
        if i == 1 
            fprintf(sc_fid,'\n%s', ['#' sc_header{i}]);
        else
            fprintf(sc_fid,'\t%s', sc_header{i});
        end
    end
    % Fill the scalar.inp file
    for i = 1 : size(sc_data_new,1)
        for j = 1 : size(sc_data_new,2)
            if j == 1
                fprintf(sc_fid,'\n%.2f', sc_data_new(i,j));
            else
                fprintf(sc_fid,'\t%f', sc_data_new(i,j));
            end
        end
    end
    % Close the file
    fclose(sc_fid);
    
    clear sc_*

else
    % Save the header data separated by spaces and remove "#"
    sc_secondline(ismember(sc_secondline, '#')) = '';
    sc_INheader = find(ismember(sc_secondline, ' ') | ismember(sc_secondline, '	') | ismember(sc_secondline, ','));
    sc_header = cell(length(sc_INheader)+1,1);
    for i = 1 : length(sc_INheader)+1
        if i == 1
            sc_header{i} = sc_secondline(1 : sc_INheader(i)-1);
        elseif i<= length(sc_INheader)
            sc_header{i} = sc_secondline(sc_INheader(i-1)+1 : sc_INheader(i)-1);
        else
            sc_header{i} = sc_secondline(sc_INheader(i-1)+1 : end);
        end
    end

    % Load the data in a matrix
    sc_data = dlmread(sc_name, '',2,0);
    % Calculate the vertical step sizes
    sc_data_delta = NaN(size(sc_data));
    for i = 1 : size(sc_data,1)-1
        sc_data_delta(i,:) = sc_data(i+1, :) - sc_data(i, :);
    end





    % Define the vertical domain of the model
    sc_zdomain = [sc_data(1,1) - (sc_data(2,1)-sc_data(1,1))/2, check.z_max];
    
    %%%%%
    % Interpolate & extrapolate the original data to a 1 cm vertical resolution
    %%%%%
    sc_int_data = NaN(length([sc_zdomain(1): 0.01 : sc_zdomain(2)]), length(sc_header));
    sc_int_data(:,1) = [sc_zdomain(1): 0.01 : sc_zdomain(2)]';
    for i = 1 : length(zold)-1

        sc_INz = find(sc_int_data(:,1) >= sc_data(i,1) & ...
                        sc_int_data(:,1) <= sc_data(i+1,1));

        for j = 2 : length(sc_header)

            if sc_data(i+1,j) == sc_data(i,j)
                sc_temp = sc_data(i,j) .* ones(length(sc_INz),1);
            else
                sc_temp = [sc_data(i,j) : ...
                          ( sc_data(i+1,j) - sc_data(i,j)) / (length(sc_INz)-1) : ...
                           sc_data(i+1,j)]';
            end
            sc_int_data(sc_INz,j) = sc_temp; 



            % Define data between 0m and z1
            if i == 1
                if  sc_data(i+1,j) == sc_data(i,j)
                    sc_temp_sfc = sc_data(i,j) .* ones(length([1 : sc_INz(1)]),1);
                else
                    sc_temp_sfc = ...
                        [sc_data(i,j) - (sc_temp(2) - sc_temp(1)) * (length([1 : sc_INz(1)])-1) : ...
                        sc_temp(2) - sc_temp(1) : ...
                        sc_data(i,j)]';
                end

                sc_int_data(1 : sc_INz(1), j) = sc_temp_sfc;

            % Define data between zend and the top of the domain
            elseif i == length(zold)-1 && sc_data(i+1,1) < sc_int_data(end,1)
                if  sc_data(i+1,j) == sc_data(i,j)

                    sc_temp_top = sc_data(i,j) .* ones(length([sc_INz(end): size(sc_int_data,1)]),1);
                else

                    sc_temp_top = ...
                        [ sc_data(i+1,j) : ...
                        sc_temp(2) - sc_temp(1) : ...
                         sc_data(i+1,j) + (sc_temp(2) - sc_temp(1)) * (length([sc_INz(end): size(sc_int_data,1)])-1)]';
                end
                sc_int_data(sc_INz(end) : size(sc_int_data,1), j) = sc_temp_top;

            end     % if i == 1

        end     % for j = 2 : length(pf_header)

    end     % for i = 1 : length(zold)-1



    sc_data_new = NaN(length(znew), length(sc_header));
    sc_data_new(:,1) = znew;
    for i = 1 : length(znew)

        for j = 2 : length(sc_header)

            sc_data_new(i,j) = sc_int_data( sc_int_data(:,1) == sc_data_new(i,1) ,j);

        end
    end



    if strcmp(check.makefigures, 'yes')
        for i = 2 : length(sc_header)
            figure
            scatter(sc_data(:,i), sc_data(:,1), 'xk')
            hold on
            scatter(sc_data_new(:,i), sc_data_new(:,1), '+r')
            plot(sc_int_data(:,i), sc_int_data(:,1))
            xlabel(sc_header{i})
            ylabel([sc_header{1} ' [m]'])
        end
    end


    %%%%%%%%%%
    % Move original scalar.inp to the input_old folder
    %%%%%%%%%%
    % Move the scalar.inp to the new folder
    movefile(sc_name, [ direc.main_directory direc.run_name direc.exp_nr olddir_name]);

    % Create the new scalar.inp file
    sc_fid = fopen(sc_name, 'w');
    % Change the first line of the scalar.inp file so it matches the new z array
    sc_firstline_new = sc_firstline;
    sc_temp = 'dz =';
    for i = 1 : length(sc_firstline) - length(sc_temp)
        if strcmp(sc_firstline(i : i+length(sc_temp)-1), sc_temp)
            sc_firstline_new = [sc_firstline(1:i-1) sc_temp ' ' num2str(znew(2) - znew(1))];
        end
    end
    sc_temp1 = 'Nz=';
    sc_temp2 = num2str(length(zold)); 
    sc_temp_IN = [];
    for i = 1: length(sc_firstline_new) - length(sc_temp1)
        if strcmp(sc_firstline_new(i : i+length(sc_temp1)-1), sc_temp1)
            sc_temp_IN = [sc_temp_IN; i];
        end
    end
    for i = 1: length(sc_firstline_new) - length(sc_temp2)
        if strcmp(sc_firstline_new(i : i+length(sc_temp2)-1), sc_temp2)
            sc_temp_IN = [sc_temp_IN; i];
        end
    end
    sc_firstline_new = [sc_firstline_new(1:sc_temp_IN(1)-1) 'Nz= ' num2str(length(znew)) ...
        sc_firstline_new(sc_temp_IN(2)+length(sc_temp2): end)];   
    % Change the first line in the scalar.inp file
    fprintf(sc_fid,'%s', sc_firstline_new);
    % Add/change the second line (header) in the scalar.inp file
    for i = 1 : length(sc_header)
        if i == 1 
            fprintf(sc_fid,'\n%s', ['#' sc_header{i}]);
        else
            fprintf(sc_fid,'\t%s', sc_header{i});
        end
    end
    % Fill the scalar.inp file
    for i = 1 : size(sc_data_new,1)
        for j = 1 : size(sc_data_new,2)
            if j == 1
                fprintf(sc_fid,'\n%.2f', sc_data_new(i,j));
            else
                fprintf(sc_fid,'\t%f', sc_data_new(i,j));
            end
        end
    end
    % Close the file
    fclose(sc_fid);
    
    clear sc_*
end

%% Load & adapt the surface.interactive.inp input file


if strcmp(check.redefine_scalars, 'yes')
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    fprintf('\n- Update surface.interactive.inp\n')  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    si_name     = [ direc.main_directory direc.run_name direc.exp_nr ...
        'surface.interactive.inp.' direc.exp_nr(1:end-1) ]; 
    si_new_name = [ direc.main_directory direc.run_name direc.exp_nr ...
        'new_surface.interactive.inp.' direc.exp_nr(1:end-1) ]; 
    si_fid = fopen(si_name,'r');

    % Read the complete input file
    si_content = cell(100,1);
    si_counter = 1;
    si_stopwhile = 0;
    while si_stopwhile == 0

        si_readline = fgetl(si_fid);

        if si_readline == -1
            si_stopwhile = 1;
        else
            si_content{si_counter} = si_readline;
            si_counter = si_counter+1;
        end
    end
    fclose(si_fid);


    % Remove the empty cells in the content array
    si_INempty = [];
    for i = 1 : length(si_content)
        if isempty(si_content{i})
            si_INempty = [si_INempty; i];
        end
    end
    si_INempty = si_INempty(1);
    si_content = si_content(1:si_INempty-1);


    % find the starting index of the input data
    si_INspace = find(ismember(si_content{2}, ' '));
    si_INdata = [];
    si_INdata_end = si_INspace(1)-1;
    for i = 2 : length(si_INspace)
        if si_INspace(i) - si_INspace(i-1) > 1
            si_INdata = [si_INdata; si_INspace(i-1)+1];
            si_INdata_end = [si_INdata_end; si_INspace(i)-1];
        end

    end
    if si_INspace(end) ~= length(si_content{2})
        si_INdata = [si_INdata; si_INspace(end)+1];
    end
    si_INdata = [2 ; si_INdata];
    si_INdata_end = [si_INdata_end ; length(si_content{2})];

    % Define the headers
    si_header = cell(length(si_INdata),1);
    for i = 1 : length(si_INdata)
        si_header{i} = si_content{2}(si_INdata(i):si_INdata_end(i));
    end
    
    if ismember({'rssoilmin'}, si_header) 
        fprintf('\n!! WARNING !!\nThe surface.interactive.inp. still contains rssoilmin in the header.\nMake sur to remove rssoilmin from the input file when using v07 or up.\n')
    end

    % Define the positions of the scalar surface flux input
    si_varsize = NaN(length(si_header),1);
    si_INsvflux = [];
    for i = 1 : length(si_header)

        % define the stringlength that is available for the input data.
        if i ~= length(si_header)
            si_varsize(i) = si_INdata(i+1) - si_INdata(i) -1;
        end

        if length(si_header{i}) >= 7
            if strcmp(si_header{i}(1:7), 'wsvsurf')
                si_INsvflux = [si_INsvflux; i];
            end
        end
    end




    % Create a new "content" array for the update
    si_new_content = si_content;
    for i = 2 : length(si_new_content)
        si_new_content{i}(si_INdata(si_INsvflux(1))-1:end) = '';
    end
    si_new_header = si_header(1 : si_INsvflux(1)-1);

    % make sure the fluxes from settings are at the right size
    if size(check.sv_flux, 2) < length(si_new_content)-2
        % Too little typenrs are defined in check.sv_flux. The check.sv_flux is
        % completed by copying typenr 0 at the end of the array.
        for i = 1 : length(si_new_content)-2 - size(check.sv_flux, 2)
            check.sv_flux = [check.sv_flux, check.sv_flux(:,1)];
        end
    end


    % Update the new content array with the data from check.sv_flux
    for i = 2 : length(si_new_content)
        for j = 1 : size(check.sv_vertprof, 1)
            % Update the header
            if i == 2
                if j < 10
                    si_svnr = ['0' num2str(j)];
                elseif j > 99
                    fprintf('\n!! ERROR !!\nThere are over 99 scalars, which are too many.\nThe program is aborted.\n')
                else
                    si_svnr = num2str(j);
                end
                
                if j == 1
                    si_new_content{i} = [si_new_content{i} ' wsvsurf(' si_svnr ')  '];
                elseif j == size(check.sv_vertprof, 1)
                    si_new_content{i} = [si_new_content{i} 'wsvsurf(' si_svnr ')'];
                else
                    si_new_content{i} = [si_new_content{i} 'wsvsurf(' si_svnr ')  '];
                end

            % Update the input data
            else
                % Preallocate the
                sv_newinput = repmat(' ',1,11);

                sv_flux = num2str(check.sv_flux(j, i-2));
                if strcmp(sv_flux, '0')
                    sv_flux = '0.0';
                end
                if ismember('.',sv_flux) == 0
                    sv_flux = [sv_flux '.0'];
                end

                if length(sv_flux) > length(sv_newinput)
                    sv_temp = sv_flux(1 : length(sv_newinput));
                    if strcmp(sv_temp(end), '.')
                        sv_temp(end) = ' ';
                    end
                    sv_newinput = sv_temp;
                else
                    sv_newinput(1 : length(sv_flux)) = sv_flux;
                end

                if j == 1
                    si_new_content{i} = [si_new_content{i} ' ' sv_newinput '  '];
                elseif j == size(check.sv_vertprof, 1)
                    si_new_content{i} = [si_new_content{i} sv_newinput];
                else
                    si_new_content{i} = [si_new_content{i} sv_newinput '  '];
                end

            end
        end
    end

    % Write a new input file "new_surface.interactive.inp"  
    si_new_fid = fopen(si_new_name,'w');
    for i = 1 : length(si_new_content)
        if i == 1
            fprintf(si_new_fid,'%s', si_new_content{i});
        else
            fprintf(si_new_fid,'\n%s', si_new_content{i});
        end
    end
    fclose(si_new_fid);

    %%%%%%%%%%
    % Move original namoptions to the input_old folder
    %%%%%%%%%%
    % Move the prof.inp to the new folder
    movefile(si_name, [ direc.main_directory direc.run_name direc.exp_nr ...
        olddir_name '\' 'surface.interactive.inp.' direc.exp_nr(1:end-1)]);
    movefile(si_new_name, si_name);
    
    clear si_*
    
end

%% Update kmax & nsc in namoptions


if strcmp(check.redefine_scalars, 'yes')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    fprintf('\n- Update kmax and nsv in namoptions\n')  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\
else
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    fprintf('\n- Update kmax in namoptions\n')  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\
end


no_name     = [ direc.main_directory direc.run_name direc.exp_nr 'namoptions.' direc.exp_nr(1:end-1) ];  
no_new_name = [ direc.main_directory direc.run_name direc.exp_nr 'new_namoptions.' direc.exp_nr(1:end-1) ];  

no_fid      = fopen(no_name,'r');
no_new_fid  = fopen(no_new_name,'w');

% read the header (on the second line = 2x)
no_stopwhile = 0;
no_stopwhile_nsv = 0;
no_counter = 0;
no_changeline = 0;
no_changeline_nsv = 0;
while no_stopwhile < 2 && no_counter < 1000
    
    %%%%%
    % Change kmax
    %%%%%
    
    no_readline = fgetl(no_fid);
    
    if no_readline == -1
        no_stopwhile = 2;
    end
    
    if length(no_readline) >= 7
        if strcmp(no_readline(1:7),'&DOMAIN')
            no_stopwhile = 1;
        end
    end
    if no_stopwhile == 1 && length(no_readline) >= 4
        if strcmp(no_readline(1:4),'kmax')
            
            no_IN = find(ismember(no_readline, '=')) +3;
            no_INfollowup = find(ismember(no_readline(no_IN:end), ' ') | ismember(no_readline(no_IN:end), char(9))) - 1;
            no_INfollowup = no_IN + no_INfollowup(1) - 1;
            
            % Update the line kmax = .... with the new length of the
            % vertical grid
            no_readline_new = [no_readline(1:no_IN-1) num2str(length(znew)) no_readline(no_INfollowup+1:end)];
            
            no_changeline = 1;
        end
    end
    
    %%%%%
    % Change nsv
    %%%%%
    
    if strcmp(check.redefine_scalars, 'yes')
        
        if length(no_readline) >= 4
            if strcmp(no_readline(1:4),'&RUN')
                no_stopwhile_nsv = 1;
            end
        end
        if no_stopwhile_nsv == 1 && length(no_readline) >= 4
            if strcmp(no_readline(1:3),'nsv')

                no_IN = find(ismember(no_readline, '=')) +3;
                no_INfollowup = find(ismember(no_readline(no_IN:end), ' ') | ismember(no_readline(no_IN:end), char(9))) - 1;
                no_INfollowup = no_IN + no_INfollowup(1) - 1;

                % Update the line nsv  = .... with the new number of
                % scalars
                no_readline_nsv_new = [no_readline(1:no_IN-1) num2str(size(check.sv_vertprof,1)) no_readline(no_INfollowup+1:end)];

                no_changeline_nsv = 1;
            end
        end
    end
    
    
    %%%%%
    % Write the new namoptions file
    %%%%%
    
    if no_stopwhile ~= 2
        if no_changeline == 1
            fprintf(no_new_fid,'\n%s', no_readline_new);
            no_changeline = 0;
        elseif no_changeline_nsv == 1
                fprintf(no_new_fid,'\n%s', no_readline_nsv_new);
                no_changeline_nsv = 0;
        else
            if no_counter == 0
                fprintf(no_new_fid,'%s', no_readline);
            else
                fprintf(no_new_fid,'\n%s', no_readline);
            end
        end
    end
    
    no_counter = no_counter+1;
    
    
end
fclose(no_fid);
fclose(no_new_fid);

clear no_*


%% Update scalar BCs in namoptions

if strcmp(check.redefine_scalars, 'yes')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    fprintf('\n- Update scalar BCs in namoptions\n')  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    no_new_name = [ direc.main_directory direc.run_name direc.exp_nr 'new_namoptions.' direc.exp_nr(1:end-1) ];  
    
    no_content = cell(5000,1);
    no_stopwhile = 0;
    no_counter = 0;
    
    % save the content of new_namoptions
    no_new_fid  = fopen(no_new_name,'r');
    while no_stopwhile ~= 2
        
        no_counter = no_counter + 1;
        
        no_readline = fgetl(no_new_fid);
        if no_readline ~= -1
            no_content{no_counter} = no_readline;
        end
        
        % Check if you are cose to the end of the namoptions file
        if length(no_readline) >= 6
            if strcmp(no_readline(1:6), '&NAMDE')
                no_stopwhile = 1; 
            end
        end
        
        % Check if the end of the namoptions file is reached
        if no_readline == -1
            if no_stopwhile == 1
                no_stopwhile = 2; 
            end
        end
        
        % Emergency stop for the while function
        if no_counter == length(no_content)
            no_stopwhile = 2; 
        end
    end
    fclose(no_new_fid);
    % Clear the empty fields at the end of no_content
    no_content(no_counter:end) = [];
    
    % Find the location of the non-periodic BC settings
    no_INbc = [];
    for i = 1 : length(no_content)
        if length(no_content{i}) >= 15
            if strcmp(no_content{i}(1:15), 'lnonperiodbc_sv')
                no_INbc = [no_INbc; i];
            end
        end
    end
    % Save the indices before and after the settings (these will be replaced)
    no_INbc = [no_INbc(1)-1, no_INbc(end)+1];
    
    
    
    % Create a array with the new BC settings
    no_newcontent = cell(length(check.sv_periodBCs),1);
    no_fillertext =  no_content{no_INbc(1)+1};
    no_INspace = find(ismember(no_fillertext, ' '));
    no_fillertext = no_fillertext(1:no_INspace(end)-1);
    no_INperiod = find(ismember(no_fillertext, '.'));
    no_INperiod = no_INperiod(2)+1;
    for i = 1 : length(no_newcontent)
        if length(num2str(check.sv_periodBCs(i))) > 1 
            no_newcontent{i} = [no_fillertext(1:16)  num2str(check.sv_periodBCs(i))...
                ') = .true.' no_fillertext(no_INperiod:end) ' #' num2str(check.sv_periodBCs(i))];
        else
            no_newcontent{i} = [no_fillertext(1:16)  num2str(check.sv_periodBCs(i))...
                ')  = .true.' no_fillertext(no_INperiod:end) ' #' num2str(check.sv_periodBCs(i))];
        end
    end
    % Complete the newcontent array by adding the new BC settings in the
    % namoptions content
    no_newcontent = [no_content(1:no_INbc(1)); no_newcontent ; no_content(no_INbc(2):end)];
    
    % Rewrite the new namoptions file
    no_new_fid  = fopen(no_new_name,'w');
    for i = 1 : length(no_newcontent)
        if i == 1
            fprintf(no_new_fid,'%s', no_newcontent{i});
        else
            fprintf(no_new_fid,'\n%s', no_newcontent{i});
        end
    end
    fclose(no_new_fid);
    
    clear no_*
    
end


%% Update scalar chemrates in namoptions

if strcmp(check.redefine_scalars, 'yes')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    fprintf('\n- Update scalar chemrates in namoptions\n')  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    no_new_name = [ direc.main_directory direc.run_name direc.exp_nr 'new_namoptions.' direc.exp_nr(1:end-1) ];  
    
    no_content = cell(5000,1);
    no_stopwhile = 0;
    no_counter = 0;
    
    % save the content of new_namoptions
    no_new_fid  = fopen(no_new_name,'r');
    while no_stopwhile ~= 2
        
        no_counter = no_counter + 1;
        
        no_readline = fgetl(no_new_fid);
        if no_readline ~= -1
            no_content{no_counter} = no_readline;
        end
        
        % Check if you are cose to the end of the namoptions file
        if length(no_readline) >= 6
            if strcmp(no_readline(1:6), '&NAMDE')
                no_stopwhile = 1; 
            end
        end
        
        % Check if the end of the namoptions file is reached
        if no_readline == -1
            if no_stopwhile == 1
                no_stopwhile = 2; 
            end
        end
        
        % Emergency stop for the while function
        if no_counter == length(no_content)
            no_stopwhile = 2; 
        end
    end
    fclose(no_new_fid);
    % Clear the empty fields at the end of no_content
    no_content(no_counter:end) = [];
    
    % Find the location of the non-periodic BC settings
    no_INbc = [];
    for i = 1 : length(no_content)
        if length(no_content{i}) >= 11
            if strcmp(no_content{i}(1:11), 'pc_chemrate')
                no_INbc = [no_INbc; i];
            end
        end
    end
    % Save the indices before and after the settings (these will be replaced)
    no_INbc = [no_INbc(1)-1, no_INbc(end)+1];
    
     
    % Create a array with the new BC settings
    no_newcontent = cell(size(check.sv_chemrate,1),1);
    no_fillertext =  no_content{no_INbc(1)+1};
    no_INspace = find(ismember(no_fillertext, ' '));
    no_fillertext = no_fillertext(1:no_INspace(end-1)-1);
    no_INcomment = find(ismember(no_fillertext, '!'));
    for i = 1 : length(no_newcontent)
        
            
        if length(num2str(check.sv_chemrate(i,1))) > 1 
            no_newcontent{i} = [no_fillertext(1:12)  num2str(check.sv_chemrate(i,1))...
                ') = ' num2str(check.sv_chemrate(i,2)) '        ' no_fillertext(no_INcomment:end) ...
                ' #' num2str(check.sv_chemrate(i,1)) ' (' check.sv_vertprof{check.sv_chemrate(i,1), 1} ')'];
        else
            no_newcontent{i} = [no_fillertext(1:12)  num2str(check.sv_chemrate(i,1))...
                ')  = ' num2str(check.sv_chemrate(i,2)) '        ' no_fillertext(no_INcomment:end) ...
                ' #' num2str(check.sv_chemrate(i,1)) ' (' check.sv_vertprof{check.sv_chemrate(i,1), 1} ')'];
        end
    end
    % Complete the newcontent array by adding the new BC settings in the
    % namoptions content
    no_newcontent = [no_content(1:no_INbc(1)); no_newcontent ; no_content(no_INbc(2):end)];
    
    % Rewrite the new namoptions file
    no_new_fid  = fopen(no_new_name,'w');
    for i = 1 : length(no_newcontent)
        if i == 1
            fprintf(no_new_fid,'%s', no_newcontent{i});
        else
            fprintf(no_new_fid,'\n%s', no_newcontent{i});
        end
    end
    fclose(no_new_fid);
    
    clear no_*
    
end


%% Update split-flux settings in namoptions

if strcmp(check.redefine_scalars, 'yes')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    fprintf('\n- Update split-flux settings in namoptions\n')  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    no_new_name = [ direc.main_directory direc.run_name direc.exp_nr 'new_namoptions.' direc.exp_nr(1:end-1) ];  
    
    no_content = cell(5000,1);
    no_stopwhile = 0;
    no_counter = 0;
    
    % save the content of new_namoptions
    no_new_fid  = fopen(no_new_name,'r');
    while no_stopwhile ~= 2
        
        no_counter = no_counter + 1;
        
        no_readline = fgetl(no_new_fid);
        if no_readline ~= -1
            no_content{no_counter} = no_readline;
        end
        
        % Check if you are cose to the end of the namoptions file
        if length(no_readline) >= 6
            if strcmp(no_readline(1:6), '&NAMDE')
                no_stopwhile = 1; 
            end
        end
        
        % Check if the end of the namoptions file is reached
        if no_readline == -1
            if no_stopwhile == 1
                no_stopwhile = 2; 
            end
        end
        
        % Emergency stop for the while function
        if no_counter == length(no_content)
            no_stopwhile = 2; 
        end
    end
    fclose(no_new_fid);
    % Clear the empty fields at the end of no_content
    no_content(no_counter:end) = [];
    
    % Find the location of the non-periodic BC settings
    no_INbc = [];
    for i = 1 : length(no_content)
        if length(no_content{i}) >= 11
            if strcmp(no_content{i}(1:10), 'sf_scalars')
                no_INbc = [no_INbc; i];
            end
        end
    end
    % Save the indices before and after the settings (these will be replaced)
    no_INbc = [no_INbc(1)-1, no_INbc(end)+1];
    
     
    % Create a array with the new BC settings
    no_newcontent = cell(size(check.sv_splitflux,1),1);
    no_fillertext =  no_content{no_INbc(1)+1};
    for i = 1 : length(no_newcontent)
            
        no_newcontent{i} = [no_fillertext(1:11) ...
            num2str(check.sv_splitflux(i,2)) ',' num2str(check.sv_splitflux(i,3)) ') = ' ...
            num2str(check.sv_splitflux(i,1)) '        ! Scalar ' num2str(check.sv_splitflux(i,1)) ...
            ' (' check.sv_vertprof{check.sv_splitflux(i,1), 1} ') in split-flux set ' ...
            num2str(check.sv_splitflux(i,2)) ' is scalar nr ' num2str(check.sv_splitflux(i,3))];
    end
    % Complete the newcontent array by adding the new BC settings in the
    % namoptions content
    no_newcontent = [no_content(1:no_INbc(1)); no_newcontent ; no_content(no_INbc(2):end)];
    
    % Rewrite the new namoptions file
    no_new_fid  = fopen(no_new_name,'w');
    for i = 1 : length(no_newcontent)
        if i == 1
            fprintf(no_new_fid,'%s', no_newcontent{i});
        else
            fprintf(no_new_fid,'\n%s', no_newcontent{i});
        end
    end
    fclose(no_new_fid);
    
    clear no_*
    
    
end

        
%% Replace namoptions with new namoptions file

%%%%%%%%%%
% Move original namoptions to the input_old folder
%%%%%%%%%%
no_name     = [ direc.main_directory direc.run_name direc.exp_nr 'namoptions.' direc.exp_nr(1:end-1) ];  
no_new_name = [ direc.main_directory direc.run_name direc.exp_nr 'new_namoptions.' direc.exp_nr(1:end-1) ];  
% Move the prof.inp to the new folder
movefile(no_name, [ direc.main_directory direc.run_name direc.exp_nr olddir_name '\namoptions.' direc.exp_nr(1:end-1)]);
movefile(no_new_name, no_name);


%% The end

   
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    fprintf('\nFinished: The input files are updated\n\n')  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%












