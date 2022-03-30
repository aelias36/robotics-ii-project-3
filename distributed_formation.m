classdef distributed_formation < handle
    properties
        xy_mat_0
        theta_vec_0

        xy_formation
        theta_formation
        v_formation
        Ts
    end
    
    methods
        function obj = distributed_formation(xy_mat_0, theta_vec_0, Ts)
            obj.xy_mat_0 = xy_mat_0;
            obj.theta_vec_0 = theta_vec_0;
            obj.Ts = Ts;

            obj.xy_formation = [2; 1.5];
            obj.theta_formation = 0;
            obj.v_formation = [0;0;0];
        end
        
        function step(obj,v_formation)
            obj.xy_formation = obj.xy_formation + v_formation(1:2)*obj.Ts;
            obj.theta_formation = obj.theta_formation + v_formation(3)*obj.Ts;
            obj.v_formation = v_formation;
        end

        function step_pos(obj,pos_formation)
            delta = pos_formation - [obj.xy_formation; obj.theta_formation];
            delta(3) = wrapToPi(delta(3));
            obj.v_formation = delta / obj.Ts;

            obj.xy_formation = pos_formation(1:2);
            obj.theta_formation = pos_formation(3);
        end

        function v_i = v_ff(obj,i_rob)
        x_i = obj.x_des(i_rob);

        v_i = NaN([3 1]);    
        J = [0 -obj.v_formation(3); obj.v_formation(3) 0];
        v_i(1:2) = obj.v_formation(1:2) + J*rot2(x_i(3))*obj.xy_mat_0(:,i_rob);
        v_i(1:2) = rot2(-x_i(3))*v_i(1:2);
        v_i(3) = obj.v_formation(3);
        end

        function x_i = x_des(obj,i_rob)
            x_i = [obj.xy_formation + rot2(obj.theta_formation)*obj.xy_mat_0(:,i_rob)
                   obj.theta_formation+obj.theta_vec_0(i_rob)];
        end
    end
end
