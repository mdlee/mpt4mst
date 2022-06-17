function [X, Y] = Raxes(A, xShift, yShift)
% function Raxes(A, xShift, yShift)
%
% make "R-style" axes
% axis handle A, and shift of xShift and yShift for the axes down and to
% the left, in the measurement units of the axes A
%
% may not work with subplots, or use of "axis square"
%
% Michael Lee (mdlee@uci.edu) 6-Apr-2016

% original properties needed for box on axes
boxOn = isequal(get(A, 'box'), 'on');
tickOut = isequal(get(A, 'tickdir'), 'out');
xColor = get(A, 'xcolor');
yColor = get(A, 'ycolor');
xLabelFontSize = get(get(A, 'xlabel'), 'fontsize'); 
yLabelFontSize = get(get(A, 'ylabel'), 'fontsize'); 

% make current axes invisible
set(A, 'xcolor', 'none', 'ycolor', 'none', 'box', 'off');

% copy axes and turn on relevant x or y
X = copyobj(A, gcf);
set(X, 'color', 'none', 'xcolor', xColor);
Y = copyobj(A, gcf);
set(Y, 'color', 'none', 'ycolor', yColor);

% remove anything on copied axes
set(X, 'ygrid', 'off');
set(X, 'xgrid', 'off');
C = get(X, 'children');
delete(C);
C = get(Y, 'children');
delete(C);

% move copies to make R-style
set(X, 'position', get(X, 'position') + [0 -xShift 0 0]);
set(Y, 'position', get(Y, 'position') + [-yShift 0 0 0]);
set(X, 'ygrid', 'off');
set(Y, 'xgrid', 'off');

% do it all again if upper and right axes if box is on
if boxOn
    
    % copy axes and turn on relevant x or y
    XU = copyobj(A, gcf);
    set(XU, 'color', 'none', 'xcolor', 'k', 'xticklabel', [], 'xlabel', [], 'tickdir', 'out');
    YR = copyobj(A, gcf);
    set(YR, 'color', 'none', 'ycolor', 'k', 'yticklabel', [], 'ylabel', [], 'tickdir', 'out');
    
    % reverse tick direction if needed
    if tickOut
        set(XU, 'tickdir', 'in');
        set(YR, 'tickdir', 'in');
    end
    
    % remove anything on copied axes
    C = get(XU, 'children');
    delete(C);
    C = get(YR, 'children');
    delete(C);
    
    % need position of axes
    position = get(A, 'position');
    
    % move copies to make R-style
    set(XU, 'position', [position(1) (position(2) + position(4) + xShift) position(3) position(4)]);
    set(YR, 'position', [(position(1) + position(3) + yShift) position(2) position(3) position(4)]);
    
end

% tidy
set(get(X, 'xlabel'), 'fontsize', xLabelFontSize);
set(get(Y, 'ylabel'), 'fontsize', yLabelFontSize);

