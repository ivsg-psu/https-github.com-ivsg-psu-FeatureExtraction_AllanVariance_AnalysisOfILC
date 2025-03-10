function avar_random_walk = fcn_AVAR_avarRandomWalk(random_walk_coefficient,...
                            list_of_correlation_time,varargin)
%% fcn_AVAR_avarRandomWalk
%   This function calculates AVAR of random walk from random walk
%   coefficient.
%
% FORMAT:
%   avar_random_walk = fcn_AVAR_avarRandomWalk(random_walk_coefficient,...
%                      list_of_correlation_time)
%
% INPUTS:
%   random_walk_coefficient: Noise coefficient for random walk [unit/sqrt(s)].
%   list_of_correlation_time: A M x 1 vector containing list of correlation time(s).
%   varargin: figure number for debugging.
%
% OUTPUTS:
%   avar_random_walk: A M x 1 vector of AVAR of random walk.
%
% EXAMPLES:
%   See the script:
%       script_test_fcn_AVAR_avarRandomWalk.m for a full test suite.
%
% This function was written on 2021_06_28 by Satya Prasad
% Questions or comments? szm888@psu.edu
% Updated: 2022/02/15

flag_do_debug = 0; % Flag to plot the results for debugging
flag_do_plot  = 0; % Flag to plot the results
flag_check_inputs = 1; % Flag to perform input checking

if flag_do_debug
    st = dbstack; %#ok<*UNRCH>
    fprintf(1, 'STARTING function: %s, in file: %s\n', st(1).name, st(1).file);
end

%% Check input arguments
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   _____                   _       
%  |_   _|                 | |      
%    | |  _ __  _ __  _   _| |_ ___ 
%    | | | '_ \| '_ \| | | | __/ __|
%   _| |_| | | | |_) | |_| | |_\__ \
%  |_____|_| |_| .__/ \__,_|\__|___/
%              | |                  
%              |_| 
% See: http://patorjk.com/software/taag/#p=display&f=Big&t=Inputs
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if flag_check_inputs
    % Are there the right number of inputs?
    if 2>nargin || 3<nargin
        error('Incorrect number of input arguments')
    end
    
    % Check input type and domain
    try
        fcn_AVAR_checkInputsToFunctions(random_walk_coefficient,'positive');
    catch ME
        assert(strcmp(ME.message,...
            'The random_walk_coefficient input must be a positive number'));
        fprintf(1, '%s\n\n', ME.message)
        return;
    end
    try
        fcn_AVAR_checkInputsToFunctions(list_of_correlation_time,'correlation time');
    catch ME
        assert(strcmp(ME.message,...
            'The list_of_correlation_time input must be a M x 1 vector of increasing positive numbers'));
        fprintf(1, '%s\n\n', ME.message)
        return;
    end
end

% Does the user want to make a plot at the end?
if 3 == nargin
    fig_num = varargin{end};
    flag_do_plot = 1;
else
    if flag_do_debug
        fig = figure;
        fig_for_debug = fig.Number;
        flag_do_plot = 1;
    end
end

%% Calculate AVAR of Random Walk
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   __  __       _       
%  |  \/  |     (_)      
%  | \  / | __ _ _ _ __  
%  | |\/| |/ _` | | '_ \ 
%  | |  | | (_| | | | | |
%  |_|  |_|\__,_|_|_| |_|
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
avar_random_walk = (random_walk_coefficient^2/3)*list_of_correlation_time;

%% Any debugging?
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   _____       _                 
%  |  __ \     | |                
%  | |  | | ___| |__  _   _  __ _ 
%  | |  | |/ _ \ '_ \| | | |/ _` |
%  | |__| |  __/ |_) | |_| | (_| |
%  |_____/ \___|_.__/ \__,_|\__, |
%                            __/ |
%                           |___/ 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if flag_do_plot
    figure(fig_num)
    clf
    width = 540; height = 400; right = 100; bottom = 200;
    set(gcf, 'position', [right, bottom, width, height])
    plot(list_of_correlation_time,avar_random_walk,'b','Linewidth',1.2)
    grid on
    set(gca,'XScale','log','YScale','log','FontSize',13)
    ylabel('Allan Variance $[Unit^2]$','Interpreter','latex','FontSize',18)
    xlabel('Correlation Time $[s]$','Interpreter','latex','FontSize',18)
end

if flag_do_debug
    fprintf(1, 'ENDING function: %s, in file: %s\n\n', st(1).name, st(1).file);
end

end