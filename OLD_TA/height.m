function [value,isterminal,direction] = height(t,y)
% Locate the time when height passes through zero in a decreasing direction
% and stop integration.
value = y(2);     % value = 0
isterminal = 1;   % 1 to stop the integration; 0 to continue
direction = -1;   % -1 if can be approached for decreasing values / +1 if can be approached for increasing values  