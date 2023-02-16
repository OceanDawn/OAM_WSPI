function p = W_angle(h)
%ANGLE  Phase angle.
%   ANGLE(H) returns the phase angles, in radians, of a matrix with
%   complex elements.  
%
%   Class support for input X:
%      float: double, single
%
%   See also ABS, UNWRAP.

%   Copyright 1984-2010 The MathWorks, Inc. 

p = atan2(imag(h), real(h));

% p =p+pi;

[r,c] = size(p);    % 读取行r、列c
for i = 1:r        % 建立for循环嵌套
    for k = 1:c
        if p(i,k) < 0
            p(i,k) = (pi-abs(p(i,k)))+pi;
        end
    end
end

