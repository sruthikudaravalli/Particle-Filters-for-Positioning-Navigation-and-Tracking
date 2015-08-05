%=====================================================================

% Purpose: To calculate the distance from sensors to target

% Arguments
% Input(s):
% target_x - horizontal location of target
% target_y - vertical location of target
% sensor_x - vector of horizontal locations of all sensors
% sensor_y - vector of vertical locations of all sensors
% Output(s):
% result - vector of distances from target to each sensor
%
% Assumptions:
% target_x - is a scalar
% target_y - is a scalar
% sensor_x - is a scalar or a vector with same length as sensor_y
% sensor_y - is a scalar or a vector with same length as sensor_x
%
% Function Declaration:
% function result = distance(target_x,target_y,sensor_x,sensor_y)
%=====================================================================
function result = distance(target_x,target_y,sensor_x,sensor_y)
no_of_sensors = length(sensor_x);
if no_of_sensors ~= length(sensor_y)
error('sensor_x and sensor_y should be same length vectors');
end
target_loc = target_x + i * target_y;
sensor_loc = sensor_x + i * sensor_y;
result = abs(sensor_loc - target_loc);
