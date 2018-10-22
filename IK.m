% Inverse Kinematics function that takes in a 4x4 homogeneous transform 
% matrix T and returns the joint variable values for a 5-DoF RRPRR 
% robotic arm. Refer documentation for DH parameters, link lengths, and
% joint variable ranges.
%
% Usage: j = IK(T)
% Input: T, 4x4 homogeneous transformation matrix, element of the group SE(3)
% Output: j, 1x5 row vector of joint variables
% If no solution is found, or the joint variables violate forward kinematic 
% invariance constraints, the element j(3) is set to -1.
%
% Function terminates if T is singular, or if the reconstruction of T from 
% the computed joint variables, TIK, is singular.
% 
% Copyright (c) 2018 Vighnesh Vatsal
% 
% Permission is hereby granted, free of charge, to any person obtaining a 
% copy of this software and associated documentation files (the "Software"), 
% to deal in the Software without restriction, including without limitation 
% the rights to use, copy, modify, merge, publish, distribute, sublicense, 
% and/or sell copies of the Software, and to permit persons to whom the 
% Software is furnished to do so, subject to the following conditions:
% 
% The above copyright notice and this permission notice shall be included in 
% all copies or substantial portions of the Software.
% 
% THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR 
% IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, 
% FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL 
% THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR 
% OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, 
% ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR 
% OTHER DEALINGS IN THE SOFTWARE.

function joint_vars = IK(T)

    if rank(T) ~= 4
        msg = 'The Forward Kinematics Matrix is Singular. Cannot solve for joint angles.';
        error(msg);    
    end

    % Link Lengths
    lengths = [-0.08, 0.045, 0.135];

    % Compute joint variables
    t_pred = zeros(1,5);

    v = reshape(T, [16,1]);
    thresh = 1e-8; % Threshold for R2z

    % Split into two cases:
    % Case 1: When R2z >> 0. R2z = v(7)
    if abs(v(7)) > thresh
        t_pred(1) = thet1(v,lengths(3));
        [c4, t_pred(2)] = thet2(t_pred(1), v);

        t_pred(4) = thet4(t_pred(2), c4, v);

        t_pred(5) = thet5(t_pred(2), c4, v);

        t_pred(3) = thet3(t_pred(1), t_pred(2), t_pred(4), t_pred(5), lengths, v);

    % Case 2: When R2z == 0
    else        
        t_pred(4) = 0;
        O0_O2 = [0; 0; lengths(1)];
        R1x = v(1);
        R1y = v(2);
        R1z = v(3);
        R3x = v(9);
        R3y = v(10);
        R3z = v(11);
        px = v(13);
        py = v(14);
        pz = v(15);           
        l3 = lengths(3);

        t_pred(1) = thet1(v,lengths(3));
        O0_O5 = [px - R1x*l3; py - R1y*l3; pz - R1z*l3];
        t1 = norm(O0_O2);
        t2 = norm(O0_O5);
        t3 = norm(O0_O2 - O0_O5);
        thet2_cand = pi - acos((t1^2 + t3^2 - t2^2)/(2*t1*t3));
        t_pred(2) = abs(thet2_cand);
        t_pred(5) = thet5(t_pred(2), 1, v);

    end
        t_pred(3) = thet3(t_pred(1), t_pred(2), t_pred(4), t_pred(5), lengths, v);


    if constraint_check(v, t_pred)
        t_pred(3) = -1;
    end

    % FK from IK predictions:
    TIK = eye(4);

    l1 = lengths(1);
    l2 = lengths(2);
    l3 = lengths(3);

    alphas = [pi/2 , pi/2 , 0 , pi/2 , pi/2, 0];
    ais = [0, 0, 0 , 0, 0, l3];
    di = [l1,0, t_pred(3) , l2, 0, 0];
    thetas = [t_pred, 0];
    thetas(3) = pi;

    for i = 1:size(thetas,2)
        T_new = [cos(thetas(i)), -cos(alphas(i))*sin(thetas(i)), sin(alphas(i))*sin(thetas(i)), ais(i)*cos(thetas(i));
                         sin(thetas(i)), cos(alphas(i))*cos(thetas(i)), -sin(alphas(i))*cos(thetas(i)), ais(i)*sin(thetas(i));
                         0, sin(alphas(i)), cos(alphas(i)), di(i);
                         0, 0, 0, 1];
        TIK = TIK*T_new;
    end

    if rank(TIK)==4
        joint_vars = t_pred;
    else
        msg = 'The reconstructed Forward Kinematics Matrix is Singular.';
        error(msg);    
    end

end


% Functions to Compute Inverse Kinematics terms based on analytical solution:

function t1 = thet1(v,l3)
    %Range: [-pi,pi]
    r1y = v(2);
    r1x = v(1);
    py = v(14);
    px = v(13);


    t1 = atan2(py - l3*r1y, px - l3*r1x);
end

function [c4,t2] = thet2(t1, v)
    %Range: [0,pi/2]

    r2x = v(5);
    r2y = v(6);
    r2z = v(7);

    A = [sin(t1), cos(t1)*r2z; -cos(t1), sin(t1)*r2z];
    avec = A\[r2x; r2y];

    %cos(t4) is also found
    c4 = avec(1);
    cand = atan2(1,avec(2));
    if cand >0 && cand<pi/2
        t2 = cand;
    elseif cand>pi/2 && cand<=pi
        t2 = -cand + pi;
    elseif cand<-pi/2 && cand>-pi
        t2 = pi + cand;
    else
        t2 = -cand;
    end
end

function d3 = thet3(t1, t2, t4, t5, l, v)
    %Range: [0.33,0.45]

    st = 0.33; %Start
    en = 0.45; %End

    thresh = 1e-8;
    pz = v(15);
    px = v(13);
    py = v(14);
    l1 = l(1);
    l2 = l(2);
    l3 = l(3);

    d3z = -pz/cos(t2) + l1/cos(t2) - l3*sin(t5) - l3*cos(t4)*cos(t5)*tan(t2) - l2;
    d3y = py/(sin(t1)*sin(t2)) - l2 - l3*sin(t5) - l3*(cos(t5)*cos(t1)*sin(t4) - cos(t5)*sin(t1)*cos(t2)*cos(t4))/(sin(t1)*sin(t2));
    d3x = px/(cos(t1)*sin(t2)) - l2 - l3*sin(t5) + l3*(cos(t5)*sin(t1)*sin(t4) + cos(t5)*cos(t1)*cos(t2)*cos(t4))/(cos(t1)*sin(t2));

    f1 = 0;
    f2 = 0;
    f3 = 0;

    if d3x <= en && d3x >= st
        f1 = 1;
    end
    if d3y <= en && d3y >= st
        f2 = 1;
    end
    if d3z <= en && d3z >= st
        f3 = 1;
    end

    if f1+f2+f3 == 3
        d3 = sum([d3x,d3y,d3z])/(3.0);
    elseif (f1+f2 == 2) && (isnan(d3z) || abs(cos(t2)) < thresh) 
         d3 = sum([d3x,d3y])/2.0;
    elseif (f1+f3 == 2) && (isnan(d3y) || abs(sin(t1)) < thresh || abs(sin(t2)) < thresh)
         d3 = sum([d3x,d3z])/2.0;
    elseif (f2+f3 == 2) && (isnan(d3x)  || abs(cos(t1)) < thresh || abs(sin(t2)) < thresh)
         d3 = sum([d3y,d3z])/2.0;
    elseif (f1) && (isnan(d3y) || abs(sin(t1)) < thresh || abs(sin(t2)) < thresh) && (isnan(d3z) || abs(cos(t2)) < thresh) 
         d3 = d3x;
    elseif (f2) && (isnan(d3x)  || abs(cos(t1)) < thresh || abs(sin(t2)) < thresh) && (isnan(d3z) || abs(cos(t2)) < thresh) 
        d3 = d3y;
    elseif (f3) && (isnan(d3x)  || abs(cos(t1)) < thresh || abs(sin(t2)) < thresh) && (isnan(d3y) || abs(sin(t1)) < thresh || abs(sin(t2)) < thresh)
        d3 = d3z;
    else
        d3 = -1; % No solution found
    end

end

function t4 = thet4(t2, c4, v)
    %Range: [-pi,pi]

    r2z = v(7);
    s4 = -r2z/sin(t2);

    t4 = atan2(s4,c4);
end

function t5 = thet5(t2, c4, v)
    %Range: [0,pi]

    r1z = v(3);
    r3z = v(11);

    A = [-cos(t2), -c4*sin(t2); -c4*sin(t2), cos(t2)];
    cvec = A\[r1z; r3z];

    t5 = atan2(cvec(1), cvec(2));
end

function result = constraint_check(v, t_pred)    
    result = 0;  
    
    r1z = v(3);
    r2z = v(7);
    r3z = v(11);
    r2x = v(5);
    r2y = v(6);
    thresh = 1e-10;
    
    %First constraint equation, c1 should be ~=0
    c1 = r1z*cos(t_pred(5))*sin(t_pred(4)) - r2z*cos(t_pred(4)) + ... 
        r3z*sin(t_pred(4))*sin(t_pred(5));
    if abs(c1) > thresh
        result = 1;
    end
    
    %Second constraint equation, c2 should be ~=0
    c2 = r2x*cos(t_pred(1))*sin(t_pred(2)) - r2z*cos(t_pred(2)) + ...
        r2y*sin(t_pred(1))*sin(t_pred(2));
    if abs(c2) > thresh
        result = 1;
    end
end