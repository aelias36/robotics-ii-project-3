classdef robot_system < handle
    properties
        anchor_mat
        p_tag_front
        p_tag_back

        lidar
        binary_map
        utrue_prev
        LK % localization kalman

        lidar_history
        odom_history
        IMU_history
        dist_front_history
        dist_back_history
        cmd_history

        robot_state
        robot_state_history

        state_estimate

        Ts

        robot_plot
        robot_estimate_plot

        avoidance_strength
        avoidance_distance
        robot_radius
    end
    
    methods
        function obj = robot_system()
            obj.Ts = 1/15;
            obj.utrue_prev = [0;0];

            obj.anchor_mat = [ 0 10 10  0
                               0  0 10 10
                               3  3  3  3];
            obj.p_tag_front = [0.5; 0];
            obj.p_tag_back = [-0.5; 0];

            obj.LK = localization_kalman();
            obj.LK.T = obj.Ts;
            obj.LK.p_tag_front = obj.p_tag_front;
            obj.LK.p_tag_back  = obj.p_tag_back;
            obj.LK.anchor_mat = obj.anchor_mat;
            obj.LK.odom_cov_mat = eye(3)*0.1.^2;
            obj.LK.UWB_cov_mat = eye(8)*0.1.^2;
            obj.LK.IMU_cov_mat = eye(3)*0.1.^2;
            obj.LK.process_cov = diag([0 0 0 0.2 0.2 0.2 1 1]).^2;

            obj.odom_history = []; 
            obj.IMU_history = [];
            obj.utrue_prev = [0;0];
            obj.dist_front_history = [];
            obj.dist_back_history = [];
            
            obj.robot_state = [2; 1.5; 0];
            %obj.robot_state = [rand_pt_in_map(obj.binary_map)'; 0];
            obj.state_estimate = [2; 1.5; 0];
            obj.robot_state_history = [];

            obj.avoidance_strength = 1e-2;
            obj.avoidance_distance = 0.2;
            obj.robot_radius = 0.2;
        end

        function step_cl(obj, x_des, v_ff, n)
            v = [0;0;0];

            if ~isempty(x_des)
                err = obj.state_estimate - x_des;
                err(3) = wrapToPi(err(3));
                v_fb = -1 * err;
                v_fb(1:2) = rot2(-obj.state_estimate(3))*v_fb(1:2);
                v = v+v_fb;
            end

            if ~isempty(v_ff)
                v = v + v_ff;               
            end

            obj.step(v, n);
        end

        function step(obj,cmd, n)
            if norm(cmd(1:2)) > norm([1 1])
                cmd(1:2) = cmd(1:2) / norm(cmd(1:2)) * norm([1 1]);
            end

            if cmd(3)>2
                cmd(3)=2;
            elseif cmd(3)<-2
                cmd(1)=-2;
            end


            [obj.robot_state,utrue] = robot_model(obj.robot_state,cmd,obj.Ts,[0.2; 0.2; 0.2]);

            obj.robot_state_history(:,n) = obj.robot_state;
            %cmd_history(:,n) = cmd;
            obj.odom_history(:,n) = utrue + normrnd([0;0;0],[0.1; 0.1; 0.1]);
            obj.IMU_history(:,n) = [ (utrue(1)-obj.utrue_prev(1))/obj.Ts; (utrue(2)-obj.utrue_prev(2))/obj.Ts; utrue(3)] + normrnd([0;0;0],[0.1; 0.1;0.1]);
            obj.utrue_prev = utrue;

            %[ranges,angles] = obj.lidar(obj.robot_state',obj.binary_map);
            %lidar_history(:,n) = ranges;

            dist_vec_front = uwb_model(obj.robot_state, obj.p_tag_front, obj.anchor_mat)' + normrnd([0;0;0;0],[1;1;1;1]*10e-2);
            dist_vec_back = uwb_model(obj.robot_state, obj.p_tag_back, obj.anchor_mat)' + normrnd([0;0;0;0],[1;1;1;1]*10e-2);
        
            %dist_front_history(:,n) = dist_vec_front;
            %dist_back_history(:,n) = dist_vec_back;

            [state_c, ~] = obj.LK.multi_meas([dist_vec_front; dist_vec_back], obj.IMU_history(:,n), obj.odom_history(:,n));
            obj.state_estimate = state_c(1:3);
        end

        function init_draw_bot(obj)
            hold on
            obj.robot_plot = plot([0 0 0 0 0], [0 0 0 0 0], 'k');
            obj.robot_estimate_plot = plot([0 0 0 0 0], [0 0 0 0 0], 'r');
            hold off
            %obj.draw_bot()
        end

        function draw_bot(obj)
            update_robot_plot(obj.robot_plot, obj.robot_state);
            update_robot_plot(obj.robot_estimate_plot, obj.state_estimate);
        end

        function v = v_avoid_fun(obj, ranges)
            ranges = ranges - obj.robot_radius;
            ranges(ranges<0) = obj.avoidance_distance/100;
            
            v = -obj.avoidance_strength./(ranges.^2) .* (1./ranges - 1/obj.avoidance_distance);
            v(ranges > obj.avoidance_distance) = 0;
            v(isnan(v)) = 0;
        end

        function v_avoid = v_avoid(obj, robot_positions)
            [ranges,angles] = obj.lidar(obj.robot_state',obj.binary_map);
            
            % also add lidar readings for the other robots
            % make sure to rotate accordingly
            for i = 1:width(robot_positions)
                xy = robot_positions(:,i) - obj.robot_state(1:2);
                if norm(xy) < 1e-3
                    continue
                end
                ranges = [ranges; norm(xy)-obj.robot_radius];
                angles = [angles; atan2(xy(1), xy(2)) - obj.robot_state(3)];
            end

            avoidance_vecs = -[cos(angles) sin(angles)];
            avoidance_vels = obj.v_avoid_fun(ranges);
            v_avoid = -avoidance_vecs' * avoidance_vels;
            v_avoid = [v_avoid; 0];
        end
        
    end
end

function update_robot_plot(robot_plot, robot_state)
d_x = 0.4;
d_y = 0.2;

x_0 = [-d_x d_x d_x -d_x -d_x NaN d_x 2*d_x]/2;
y_0 = [-d_y -d_y d_y d_y -d_y NaN 0 0]/2;

xy = rot2(robot_state(3)) * [x_0; y_0] + robot_state(1:2);


set(robot_plot, "XData", xy(1,:), "YData", xy(2,:));

end
