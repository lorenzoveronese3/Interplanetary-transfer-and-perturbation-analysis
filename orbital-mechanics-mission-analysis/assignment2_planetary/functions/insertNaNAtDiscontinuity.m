function [xOut, yOut] = insertNaNAtDiscontinuity(x, y, jumpThresholdDeg)
%INSERTNANATDISCONTINUITY Inserts NaN at discontinuities in x.
%
% PROTOTYPE:
%   [xOut, yOut] = insertNaNAtDiscontinuity(x, y, jumpThresholdDeg)
%
% PURPOSE:
%   When plotting a ground track, longitude is often wrapped
%   in the interval [-180, +180] (or [-pi, +pi]). When the trajectory crosses
%   the map edge, x may show a jump (e.g., +179 -> -179).
%   By inserting a NaN between the two points, MATLAB breaks the polyline and avoids
%   a "connecting" line that crosses the entire map.
%
% INPUT:
%   x                [Nx1]  Abscissas (typically longitude) [deg]
%   y                [Nx1]  Ordinates (typically latitude) [deg]
%   jumpThresholdDeg [1]    Threshold above which a jump in x is considered a discontinuity [deg]
%
% OUTPUT:
%   xOut, yOut        [Mx1] Vectors with NaNs inserted after each segment (M >= N)
%
% AUTHORS:
%   Simone Frescura
%   Nicola Pifano
%   Luca Puzzi
%   Lorenzo Veronese
%==========================================================================

x = x(:);
y = y(:);

% Indices where the jump between consecutive samples exceeds the threshold
dx = abs(diff(x));
breakIdx = find(dx > jumpThresholdDeg);

% Split into segments and add a NaN after each "closed" segment
xCells = cell(numel(breakIdx)+1, 1);
yCells = cell(numel(breakIdx)+1, 1);

iStart = 1;
for k = 1:numel(breakIdx)
    iEnd = breakIdx(k);
    xCells{k} = [x(iStart:iEnd); NaN];
    yCells{k} = [y(iStart:iEnd); NaN];
    iStart = iEnd + 1;
end

% Last segment (after the last jump, or the whole vector if there are no jumps)
xCells{end} = x(iStart:end);
yCells{end} = y(iStart:end);

xOut = vertcat(xCells{:});
yOut = vertcat(yCells{:});

end

