function [K R1 R2] = getStereoCam(f, cx, cy, campos1, campos2)


K = [...
    f 0 cx; 
    0 f cy;
    0 0 1];

%%%%%%%%%%%%%%%%%%%
%% cam pos in world frame
R1 = eye(3);
R1 = [R1 campos1';0 0 0 1];

R2 = eye(3);
R2 = [R2 campos2';0 0 0 1];
%%%%%%%%%%%%%%%%%%%
%%make transformation to get camera frame
R1 = inv(R1);
R2 = inv(R2);
%%%%%%%%%%%%%%%%%%%

