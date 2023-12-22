function [angle] = qck(angle)

% qck.m - Reduce an angle between 0 and 2*pi
%
% PROTOTYPE:
%   [angle]=qck(angle)
%
% DESCRIPTION:
%   This function takes any angle and reduces it, if necessary,
% 	so that it lies in the range from 0 to 2 PI radians.
% 
% INPUTS:
%   ANGLE[1]    Angle to be reduced (in radians)
% 
% OUTPUTS:
%   QCK[1]      The angle reduced, if necessary, to the range
%               from 0 to 2 PI radians (in radians)
% 
% CALLED FUNCTIONS:
%   pi (from MATLAB)
%
% AUTHOR:
%   W.T. Fowler, July, 1978
%
% CHANGELOG:
%   8/20/90, REVISION: Darrel Monroe
%
% -------------------------------------------------------------------------

twopi = 2*pi;
 
diff = twopi * (fix(angle/twopi) + min([0,sign(angle)]));

angle = angle -diff;

return
