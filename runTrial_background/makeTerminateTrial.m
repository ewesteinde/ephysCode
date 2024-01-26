% Run NiDAQ acquisition in the background
% read in logged data & resaves it as a mat file

% create app to monitor whether I want to end the trial or not & print
% status to terminal
        
function button = makeTerminateTrial()

    % Create UIFigure and hide until all components are created
    button.UIFigure = uifigure('Visible', 'off');
    button.UIFigure.Position = [100 100 199 100];
    button.UIFigure.Name = 'MATLAB App';

    % Create EndTrialButton
    button.EndTrialButton = uibutton(button.UIFigure, 'state');
    button.EndTrialButton.Position = [41 33 121 36];
    button.EndTrialButton.Text = 'End Trial';
    button.EndTrialButton.ValueChangedFcn = @(src, event)CL_OL_ButtonPress;
    % Show the figure after all components are created
    button.UIFigure.Visible = 'on';
end

