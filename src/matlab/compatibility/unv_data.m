classdef unv_data < handle
    properties
        dataset_label
        dataset_name
        dataset_location
        textual_id
        model_type
        analysis_type
        analysis_data
        data_characteristic
        result_type
        data_type
        nvaldc
        nodes
        values
    end
    
    properties (Constant, Access = private)
        dataset_locations = {
            1, 'Data at nodes'
            2, 'Data on elements'
            3, 'Data at nodes on elements'
            5, 'Data at points'
            };
        
        model_types = {
            0, 'Unknown'
            1, 'Structural'
            2, 'Heat transfer'
            3, 'Fluid flow'
            };
        
        analysis_types = {
            0,   'Unknown'
            1,   'Static'
            2,   'Normal modes'
            3,   'Complex eigenvalue first order'
            4,   'Transient'
            5,   'Frequency response'
            6,   'Buckling'
            7,   'Complex eigenvalue second order'
            9,   'Static non-linear'
            10,   'Craig-Bampton constraint modes'
            11,   'Equivalent attachment modes'
            12,   'Effective mass modes'
            13,   'Effective mass matrix'
            14,   'Effective mass matrix'
            };
        
        data_characteristics = {
            0,   'Unknown'
            1,   'Scalar'
            2,   '3 DOF global translation vector'
            3,   '6 DOF global translation & rotation vector'
            4,   'Symmetric global tensor'
            6,   'Stress resultants'
            };
        
        result_types = {
            8, 'Displacement'
            9, 'Reaction Force'
            11, 'Velocity'
            12, 'Acceleration'
            117, 'Pressure'
            301, 'Sound Pressure'
            302, 'Sound Power'
            303, 'Sound Intensity'
            304, 'Sound Energy'
            };
        
        data_types = {
            1, 'Integer'
            2, 'Single precision floating point'
            4, 'Double precision floating point'
            5, 'Single precision complex'
            6, 'Double precision complex'
            };
    end
    
    methods (Access = public)
        function read_from_cellstr(obj, input)
            
            obj.dataset_label = sscanf(input{1}, '%u');
            obj.dataset_name = input{2};
            
            rec3 = sscanf(input{3}, '%10u');
            obj.dataset_location = unv_data.id2text(rec3(1), unv_data.dataset_locations);
            
            obj.textual_id = input(4:8);
            rec9 = sscanf(input{9}, '%u', 6);
            obj.model_type = unv_data.id2text(rec9(1), unv_data.model_types);
            obj.analysis_type = unv_data.id2text(rec9(2), unv_data.analysis_types);
            obj.data_characteristic = unv_data.id2text(rec9(3), unv_data.data_characteristics);
            obj.result_type = unv_data.id2text(rec9(4), unv_data.result_types);
            obj.data_type = unv_data.id2text(rec9(5), unv_data.data_types);
            iscplx = ~isempty(strfind(obj.data_type, 'complex'));
            
            obj.nvaldc = rec9(6);
            
            integer_data = [sscanf(input{10}, '%10u', 8); sscanf(input{11}, '%10u', 2)];
            real_data = [sscanf(input{12}, '%g', 6); sscanf(input{13}, '%g', 6)];
            obj.get_analysis_data(integer_data, real_data);
            
            input = input(14:end);
            
            stream = strcat(input, {' '});
            stream = horzcat(stream{:});
            n = obj.nvaldc;
            if iscplx
                n = n * 2;
            end
            format = ['%u' repmat('%g', 1, n)];
            numbers = sscanf(stream, format, [n+1 Inf]);
            obj.nodes = numbers(1,:)';
            obj.values = numbers(2:end,:)';
            if iscplx
                obj.values = complex(obj.values(:,1:2:end), obj.values(:,2:2:end));
            end
        end
        
        function write_to_file(obj, fid)
            fprintf(fid, '%6d\n', -1);                      % block start
            fprintf(fid, '%6d\n', 2414);                    % data block
            fprintf(fid, '%10d\n', obj.dataset_label);
            fprintf(fid, '%s\n', obj.dataset_name);
            fprintf(fid, '%10d\n', unv_data.text2id(obj.dataset_location, unv_data.dataset_locations));
            for i = 1 : 5
                fprintf(fid, '%s\n', obj.textual_id{i});
            end
            
            fprintf(fid, '%10d%10d%10d%10d%10d%10d\n',...
                unv_data.text2id(obj.model_type, unv_data.model_types), ...
                unv_data.text2id(obj.analysis_type, unv_data.analysis_types), ...
                unv_data.text2id(obj.data_characteristic, unv_data.data_characteristics), ...
                unv_data.text2id(obj.result_type, unv_data.result_types), ...
                unv_data.text2id(obj.data_type, unv_data.data_types), ...
                obj.nvaldc);
            [idata, rdata] = obj.get_analysis_vectors();
            fprintf(fid, '%10d%10d%10d%10d%10d%10d%10d%10d\n%10d%10d\n', idata(1:10));
            fprintf(fid, [repmat('%13.5e', 1, 6) '\n' repmat('%13.5e', 1, 6) '\n'], rdata);
            for i = 1 : length(obj.nodes)
                fprintf(fid, '%10d\n', obj.nodes(i));
                fprintf(fid, [repmat('%13.5e', 1, 6), '\n'], [real(obj.values(i,:)); imag(obj.values(i,:))]);
            end
            fprintf(fid, '%6d\n', -1);
        end
    end
    
    methods (Access = private)
        function get_analysis_data(obj, integer_data, real_data)
            switch obj.analysis_type
                case 'Normal modes'
                    d.design_set_id = integer_data(1);
                    d.iteration_number = integer_data(2);
                    d.solution_set_id = integer_data(3);
                    d.boundary_condition = integer_data(4);
                    d.mode_number = integer_data(6);
                    d.creation_option = integer_data(9);
                    d.number_retained = integer_data(10);
                    d.frequency = real_data(2);
                    d.modal_mass = real_data(4);
                    d.viscous_damping_ratio = real_data(5);
                    d.hysteretic_damping_ratio = real_data(6);
                case 'Frequency response'
                    d.design_set_id = integer_data(1);
                    d.solution_set_id = integer_data(3);
                    d.boundary_condition = integer_data(4);
                    d.load_set = integer_data(5);
                    d.frequency_number = integer_data(8);
                    d.creation_option = integer_data(9);
                    d.number_retained = integer_data(10);
                    d.frequency = real_data(2);
                case 'Complex eigenvalue first order'
                    d.design_set_id = integer_data(1);
                    d.solution_set_id = integer_data(3);
                    d.boundary_condition = integer_data(4);
                    d.load_set = integer_data(5);
                    d.mode_number = integer_data(6);
                    d.creation_option = integer_data(9);
                    d.number_retained = integer_data(10);
                    d.eigenvalue = complex(real_data(7), real_data(8));
                    d.modal_A = complex(real_data(9), real_data(10));
                    d.modal_B = complex(real_data(11), real_data(12));
                otherwise
                    aaa = 3;
            end
            obj.analysis_data = d;
        end
        
        
        function [integer_data, real_data] = get_analysis_vectors(obj)
            d = obj.analysis_data;
            integer_data = zeros(16,1);
            real_data = zeros(12,1);
            switch obj.analysis_type
                case 'Normal modes'
                    integer_data(1) = d.design_set_id;
                    integer_data(2) = d.iteration_number;
                    integer_data(3) = d.solution_set_id;
                    integer_data(4) = d.boundary_condition;
                    integer_data(6) = d.mode_number;
                    integer_data(9) = d.creation_option;
                    integer_data(10) = d.number_retained;
                    real_data(2) = d.frequency;
                    real_data(4) = d.modal_mass;
                    real_data(5) = d.viscous_damping_ratio;
                    real_data(6) = d.hysteretic_damping_ratio;
                case 'Frequency response'
                    integer_data(1) = d.design_set_id;
                    integer_data(3) = d.solution_set_id;
                    integer_data(4) = d.boundary_condition;
                    integer_data(5) = d.load_set;
                    integer_data(8) = d.frequency_number;
                    integer_data(9) = d.creation_option;
                    integer_data(10) = d.number_retained;
                    real_data(2) = d.frequency;
                case 'Complex eigenvalue first order'
                    integer_data(1) = d.design_set_id;
                    integer_data(3) = d.solution_set_id;
                    integer_data(4) = d.boundary_condition;
                    integer_data(5) = d.load_set;
                    integer_data(6) = d.mode_number;
                    integer_data(9) = d.creation_option;
                    integer_data(10) = d.number_retained;
                    real_data(7) = real(d.eigenvalue);
                    real_data(8) = imag(d.eigenvalue);
                    real_data(9) = real(d.modal_A);
                    real_data(10) = imag(d.modal_A);
                    real_data(11) = real(d.modal_B);
                    real_data(12) = imag(d.modal_B);
                otherwise
                    aaa = 3;
            end
        end
        
    end
    
    methods (Static = true, Access = private)
        function txt = id2text(id, cell)
            sel = id == cell2mat(cell(:,1));
            txt = cell{sel,2};
        end
        
        function id = text2id(txt, cell)
            sel = strcmp(cell(:,2), txt);
            id = cell{sel,1};
        end
    end
end % of class
