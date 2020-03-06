save_output = true;

% Simulations will be run for any study set to true.
runs = containers.Map;
runs('general')        = true;
runs('boundary_size')  = true;
runs('convergence')    = true;
runs('mesh_fineness')  = true;
runs('hiramatsu')      = true;


% GENERAL
if (runs('general'))
    close all; clear functions; clearvars -except save_output runs

    param_file = 'param_files/paper_param_files/earth_sun_model.m';
    earth_boxes = false;
    ic = 'sum_onebodies';
    fineness_level = 4;
    log2au_bound = 6;

    for sun_only = [true false]
        for earth_only = [true false]
            if (earth_only && sun_only)
                continue
            end

            close all; clear functions; clearvars -except save_output runs...
                param_file log2au_bound sun_only earth_only earth_boxes ic seed fineness_level
            try
                vain_findiff;
            except ME
                disp( getReport( ME, 'extended', 'hyperlinks', 'on') );
            end

        end
    end

    close all; clear functions; clearvars -except save_output runs

    param_file = 'param_files/paper_param_files/moon_earth_sun_model.m';
    earth_moon_boxes = false;
    ic = 'sum_onebodies';
    fineness_level = 4;
    log2au_bound = 6;

    for sun_only = [true false]
    for earth_only = [true false]
    for sun_earth_only = [true false]
    for earth_moon_only = [true false]
    for moon_only = [true false]
        if sun_only + earth_only + sun_earth_only + moon_only + earth_moon_only > 1
            continue
        end

        close all; clear functions; clearvars -except save_output runs...
            param_file log2au_bound sun_only sun_earth_only earth_only moon_only earth_moon_only earth_moon_boxes ic seed fineness_level
        %try
            vain_findiff_3d;
        %except ME
            %disp( getReport( ME, 'extended', 'hyperlinks', 'on') );
        %end

    end
    end
    end
    end
    end
end

% BOUNDARY SIZE
if (runs('boundary_size'))
    close all; clear functions; clearvars -except save_output runs

    param_file = 'param_files/paper_param_files/earth_sun_model.m';
    earth_boxes = false;
    earth_only = false;
    ic = 'sum_onebodies';
    fineness_level = 4;

    for sun_only = [true false]
        for log2au_bound = 1:9
            close all; clear functions; clearvars -except save_output runs...
                param_file log2au_bound sun_only earth_only earth_boxes ic seed fineness_level
            try
                vain_findiff;
            except ME
                disp( getReport( ME, 'extended', 'hyperlinks', 'on') );
            end

        end
    end



    close all; clear functions; clearvars -except save_output runs

    param_file = 'param_files/paper_param_files/moon_earth_sun_model.m';
    earth_moon_boxes = false;
    earth_only = false;
    moon_only = false;
    earth_moon_only = false;
    ic = 'sum_onebodies';
    fineness_level = 4;

    for sun_only = [true false]
        for sun_earth_only = [true false]
            if ~(sun_only && sun_earth_only)
                for log2au_bound = 1:9
                    close all; clear functions; clearvars -except save_output runs...
                        param_file log2au_bound sun_only earth_only moon_only sun_earth_only earth_moon_only earth_moon_boxes ic seed fineness_level
                    try
                        vain_findiff_3d;
                    except ME
                        disp( getReport( ME, 'extended', 'hyperlinks', 'on') );
                    end
                end
            end
        end
    end
end


% CONVERGENCE
if (runs('convergence'))
    close all; clear functions; clearvars -except save_output runs

    param_file = 'param_files/paper_param_files/earth_sun_model.m';
    earth_boxes = false;
    earth_only = false;
    fineness_level = 4;
    log2au_bound = 4;

    for sun_only = false; %[true false]
        ic = 'sum_onebodies';
        close all; clear functions; clearvars -except save_output runs...
            param_file log2au_bound sun_only earth_only earth_boxes ic seed fineness_level
        try
            vain_findiff;
        except ME
            disp( getReport( ME, 'extended', 'hyperlinks', 'on') );
        end

        ic = 'zero';
        close all; clear functions; clearvars -except save_output runs...
            param_file log2au_bound sun_only earth_only earth_boxes ic seed fineness_level
        try
            vain_findiff;
        except ME
            disp( getReport( ME, 'extended', 'hyperlinks', 'on') );
        end

        ic = 'random';
        for seed = 1:5
            close all; clear functions; clearvars -except save_output runs...
                param_file log2au_bound sun_only earth_only earth_boxes ic seed fineness_level
            try
                vain_findiff;
            except ME
                disp( getReport( ME, 'extended', 'hyperlinks', 'on') );
            end
        end
    end


    close all; clear functions; clearvars -except save_output runs

    param_file = 'param_files/paper_param_files/moon_earth_sun_model.m';
    earth_moon_boxes = false;
    sun_earth_only = false;
    earth_moon_only = false;
    earth_only = false;
    moon_only = false;
    fineness_level = 2;
    log2au_bound = 4;

    for sun_only = [true false]
        ic = 'sum_onebodies';
        close all; clear functions; clearvars -except save_output runs...
            param_file log2au_bound sun_only earth_only moon_only sun_earth_only earth_moon_only earth_moon_boxes ic seed fineness_level
        try
            vain_findiff_3d;
        except ME
            disp( getReport( ME, 'extended', 'hyperlinks', 'on') );
        end

        ic = 'zero';
        close all; clear functions; clearvars -except save_output runs...
            param_file log2au_bound sun_only earth_only moon_only sun_earth_only earth_moon_only earth_moon_boxes ic seed fineness_level
        try
            vain_findiff_3d;
        except ME
            disp( getReport( ME, 'extended', 'hyperlinks', 'on') );
        end

        ic = 'random';
        for seed = 1:5
            close all; clear functions; clearvars -except save_output runs...
                param_file log2au_bound sun_only earth_only moon_only sun_earth_only earth_moon_only earth_moon_boxes ic seed fineness_level
            try
                vain_findiff_3d;
            except ME
                disp( getReport( ME, 'extended', 'hyperlinks', 'on') );
            end
        end
    end
end


% MESH FINENESS
if (runs('mesh_fineness'))
    close all; clear functions; clearvars -except save_output runs

    param_file = 'param_files/paper_param_files/earth_sun_model.m';
    earth_boxes = false;
    earth_only = false;
    log2au_bound = 6;
    ic = 'sum_onebodies';

    for sun_only = [true false]
        for fineness_level = 1:6
            close all; clear functions; clearvars -except save_output runs...
                param_file log2au_bound sun_only earth_only earth_boxes ic seed fineness_level
            try
                vain_findiff;
            except ME
                disp( getReport( ME, 'extended', 'hyperlinks', 'on') );
            end
        end
    end


    close all; clear functions; clearvars -except save_output runs

    param_file = 'param_files/paper_param_files/moon_earth_sun_model.m';
    earth_moon_boxes = false;
    sun_earth_only = false;
    earth_moon_only = false;
    earth_only = false;
    moon_only = false;
    log2au_bound = 6;
    ic = 'sum_onebodies';

    for sun_only = [true false]
        for fineness_level = 1:6
            close all; clear functions; clearvars -except save_output runs...
                param_file log2au_bound sun_only earth_only moon_only sun_earth_only earth_moon_only earth_moon_boxes ic seed fineness_level
            try
                vain_findiff_3d;
            except ME
                disp( getReport( ME, 'extended', 'hyperlinks', 'on') );
            end
        end
    end

end


% HIRAMATSU BOUNDARY
if (runs('hiramatsu'))
    close all; clear functions; clearvars -except save_output runs

    param_file = 'param_files/paper_param_files/hiramatsu_model.m';
    fineness_level = 5;

    for log2au_bound = 1:10
        close all; clear functions; clearvars -except save_output runs...
            param_file log2au_bound fineness_level
        try
            vain_findiff;
        except ME
            disp( getReport( ME, 'extended', 'hyperlinks', 'on') );
        end

    end
end
