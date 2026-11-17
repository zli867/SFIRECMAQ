% Improved J calculation and entrainment
%% Load Data
clearvars -except INPUT_MAT OUTPUT_MAT
clc;
if isempty(which('F0AM_ModelCore'))
    warning('F0AM not on path (standalone run) - adding default Model path.')
    addpath(genpath('/Users/zongrunli/Desktop/SFIRECMAQV2/F0AM/Model'))
end
if ~exist('INPUT_MAT','var')
    INPUT_MAT  = 'foam_input.mat'; 
end
load(INPUT_MAT)
%% Location Define
o = ones(size(P));
% location.latitude   = latitude * o;
% location.longitude  = longitude * o;
% location.altitude   = altitude * o;
cur_time = datetime(time);
dt = cur_time(2) - cur_time(1);
clear time
time.year = year(cur_time);
time.month          = month(cur_time);
time.day            = day(cur_time);
time.hour           = hour(cur_time);
time.min            = minute(cur_time);
time.sec            = second(cur_time);
time.UTC = 0;
% sun = sun_position(time, location);
%% Chemistry
ChemFiles = {...
   'CB6r2_K(Met)';
   'CB6r2_J(Met,3)';
   'CB6r2_AllRxns';
   };
%% DILUTION CONCENTRATIONS
% We are not diluting the chamber air, so this input is irrelevant (but still necessary).

BkgdConc = {...
%   names           values
    'DEFAULT'       0;   %0 for all zeros, 1 to use InitConc
    };

%% Model Options
ModelOptions.Verbose       = 1;
ModelOptions.EndPointsOnly = 1;
ModelOptions.LinkSteps     = 0;
ModelOptions.IntTime       = seconds(dt);
ModelOptions.SavePath      = 'SFIREOutput.mat';
ModelOptions.GoParallel    = 0;
%% Init Conc
fields = fieldnames(conc);
InitConc = cell(length(fields), 3);
for i = 1:length(fields)
    species = fields{i};
    InitConc{i, 1} = species;
    InitConc{i, 2} = 0;
    InitConc{i, 3} = 0;
end

InitMass = cell(length(fields), 2);
for i = 1:length(fields)
    species = fields{i};
    InitMass{i, 1} = species;
    InitMass{i, 2} = 0;
end
%% Main Logic
mass = cell(length(o), 2);  % Preallocate with 2 columns
res = cell(length(o), 2);  % Preallocate with 2 columns
res(1, :) = {cur_time(1),InitConc};
mass(1, :) = {cur_time(1),InitMass};
cur_time(end+1) = cur_time(end) + dt;

for t_step = 1:length(o)
    % Create current J
    cur_J = struct();
    % Get field names of J
    J_fields = fieldnames(J);

    % Iterate over each field and assign the first value to the new struct
    for i = 1:length(J_fields)
        fname = J_fields{i};        % Field name
        cur_J.(fname) = double(J.(fname)(t_step));  % First element of the array
    end

    Met = {...
        'P'          double(P(t_step));
        'T'          double(T(t_step));
        'RH'         double(RH(t_step));
        % 'SZA'        double(sun.zenith(t_step));
        'J'          cur_J
        };

    % Calculate this step's Concentration
    for i = 1:length(fields)
        species = fields{i};
        idx = strcmp(InitConc(:, 1), species);
        InitConc{idx, 2} = InitConc{idx, 2} + conc.(species)(t_step);
    end

    % Run F0AM Model
    S = F0AM_ModelCore(Met, InitConc, ChemFiles, BkgdConc, ModelOptions);

    % Calculate next step's Concentration
    if t_step < length(o)
        factor = (1/(1 + ratio(t_step))) * (V(t_step) / V(t_step + 1));
        cur_species_names = fieldnames(S.Conc);
        % Update InitConc with the current species' concentrations
        for i = 1:length(cur_species_names)
            if strcmp(cur_species_names{i}, 'ONE')
                continue;
            end
            idx = strcmp(InitConc(:, 1), cur_species_names{i});
            if any(idx)
                InitConc{idx, 2} = S.Conc.(cur_species_names{i}) * factor;  % Update concentration
            else
                InitConc{end + 1, 1} = cur_species_names{i};  % Add new species to InitConc
                InitConc{end, 2} = S.Conc.(cur_species_names{i}) * factor;
            end
            InitConc{end, 3} = 0;  % Default for the third column
        end
    end

    % calculate M_x^t, M_x^t = conc_x^t * V^t
    mass_factor = (1/(1 + ratio(t_step))) * V(t_step);
    disp(mass_factor)
    disp(factor)
    cur_species_names = fieldnames(S.Conc);
    % Update InitConc with the current species' concentrations
    for i = 1:length(cur_species_names)
        if strcmp(cur_species_names{i}, 'ONE')
            continue;
        end
        idx = strcmp(InitMass(:, 1), cur_species_names{i});
        if any(idx)
            InitMass{idx, 2} = S.Conc.(cur_species_names{i}) * mass_factor;
        else
            InitMass{end + 1, 1} = cur_species_names{i};
            InitMass{end, 2} = S.Conc.(cur_species_names{i}) * mass_factor;
        end
    end

    if t_step < length(o)
        res(t_step + 1, :) = {cur_time(t_step + 1), InitConc};
    end
    mass(t_step + 1, :) = {cur_time(t_step + 1),InitMass};
end

% % Convert Res to struct
% for i = 1: length(cur_species_names)
%     cur_species_name = cur_species_names(i);
%     if strcmp(cur_species_names{i}, 'ONE')
%         continue;
%     end
%     values = [0] * length(res);
%     times = [0] * length(res);
%     for j = 1: length(res)
%         times(j) = posixtime(res{j, 1});  % Convert to UNIX timestamp
%         cur_conc = res{j, 2};
%         idx = strcmp(cur_conc(:, 1), cur_species_name);
%         if any(idx)
%             value = cur_conc{idx, 2};
%         else
%             value = 0;
%         end
%         values(j) = value;
%     end
%     if i == length(cur_species_names)
%         eval('s.time=times;');
%     end
%     eval("s." + cur_species_name + "=values;");
% end
% save('sfire_res.mat','s')

% Convert Res to struct
for i = 1: length(cur_species_names)
    cur_species_name = cur_species_names(i);
    if strcmp(cur_species_names{i}, 'ONE')
        continue;
    end
    values = [0] * length(mass);
    times = [0] * length(mass);
    for j = 1: length(mass)
        times(j) = posixtime(mass{j, 1});  % Convert to UNIX timestamp
        cur_conc = mass{j, 2};
        idx = strcmp(cur_conc(:, 1), cur_species_name);
        if any(idx)
            value = cur_conc{idx, 2};
        else
            value = 0;
        end
        values(j) = value;
    end
    if i == length(cur_species_names)
        eval('s.time=times;');
    end
    eval("s." + cur_species_name + "=values;");
end
if ~exist('OUTPUT_MAT','var'); OUTPUT_MAT = 'sfire_res.mat'; end
save(OUTPUT_MAT, 's')