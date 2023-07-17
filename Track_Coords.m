function [fullExcel, geom] = Track_Coords(fullExcel, ExcelArr, MatlabRows, readCoords, saveCoords, mainDir, tableName, onMac)

% Works with Track_Master to track a mouse in a 2D arena
% Allows you to set, save, and read coordinates from set Excel sheet
% Last updated 7/2023, Danny Lasky

%% Prepare variables for the coordinates of the arena border, mouse start, and toy locations
geom.Edges      = cell(length(MatlabRows),1);
geom.Ctr        = cell(length(MatlabRows),1);
geom.Perim      = cell(length(MatlabRows),1);
geom.DspEdges   = cell(length(MatlabRows),1);
geom.DspPerim   = cell(length(MatlabRows),1);
geom.MouseCtr   = cell(length(MatlabRows),1);
geom.ObjPos     = cell(length(MatlabRows),1);

%% Read in pre-set coordinates if toggle is on and coordinates are available
if readCoords == 1
    for n = 1:length(MatlabRows)
        geom.Edges{n}(1,1) = ExcelArr{n,'EdgesX1'};
        geom.Edges{n}(1,2) = ExcelArr{n,'EdgesY1'};
        geom.Edges{n}(2,1) = ExcelArr{n,'EdgesX2'};
        geom.Edges{n}(2,2) = ExcelArr{n,'EdgesY2'};
        geom.Edges{n}(3,1) = ExcelArr{n,'EdgesX3'};
        geom.Edges{n}(3,2) = ExcelArr{n,'EdgesY3'};
        geom.Edges{n}(4,1) = ExcelArr{n,'EdgesX4'};
        geom.Edges{n}(4,2) = ExcelArr{n,'EdgesY4'};
        
        geom.MouseCtr{n}(1,1) = ExcelArr{n,'MouseX'};
        geom.MouseCtr{n}(1,2) = ExcelArr{n,'MouseY'};
        
        if ExcelArr{n,'Phase'} >= 2
            geom.ObjPos{n}(1,1) = ExcelArr{n,'ObjAX'};
            geom.ObjPos{n}(1,2) = ExcelArr{n,'ObjAY'};
            geom.ObjPos{n}(2,1) = ExcelArr{n,'ObjBX'};
            geom.ObjPos{n}(2,2) = ExcelArr{n,'ObjBY'};
        end
        
        if ExcelArr{n,'Phase'} == 1
            errorSum = mean2(isnan(geom.Edges{n})) + mean2(isnan(geom.MouseCtr{n}));
        elseif ExcelArr{n,'Phase'} >= 2
            errorSum = mean2(isnan(geom.Edges{n})) + mean2(isnan(geom.MouseCtr{n})) + mean2(isnan(geom.ObjPos{n}));
        end
        
        if errorSum ~= 0
            error('Attempted to read-in unset coordinates. Rerun files with "readCoords = 0;"')
        end
        
%% Calculate other geometric parameters from those read in
        geom.Ctr{n}   = mean(geom.Edges{n});
        geom.Perim{n} = (0.8 .* (geom.Edges{n} - repmat(geom.Ctr{n}, size(geom.Edges{n},1),1))) ...
        	+ repmat(geom.Ctr{n}, size(geom.Edges{n},1),1);
      	geom.DspEdges{n} = [geom.Edges{n}; geom.Edges{n}(1,:)];
     	geom.DspPerim{n} = [geom.Perim{n}; geom.Perim{n}(1,:)];
    end
end

%% Set coordinates if readCoords toggle is off
if readCoords == 0
    for n = 1:length(MatlabRows)
        tempObj = char(ExcelArr{n,'VideoPath'});
        
        if onMac == 1
            tempObj = strrep(tempObj, 'P:\Jones_Maganti_Shared\Mouse Tracking', '/Volumes/mathewjones/Jones_Maganti_Shared/Mouse Tracking');
            tempObj = strrep(tempObj, '\', '/');
        end

        tempObj = strrep(tempObj,'.mp4','_Red.avi');
        vidObj = VideoReader(tempObj);
        
        startFrameNum = ExcelArr{n,'StartFrame'};
        startFrame = read(vidObj, startFrameNum);
        success = zeros(3,1);
        
        % Initialize drawing figure
        figure('units', 'pixels');
        set(gcf, 'Units', 'Normalized', 'OuterPosition', [0.025 0.025 0.475 0.95]);
        
        ax1 = subplot(2,1,1); 
    	image(startFrame); 
    	axis image
    	hold on
        figName = strcat("Animal ", char(ExcelArr{n,'AnimalID'}), " Wk ", num2str(ExcelArr{n,'Week'}), ...
                " Coh ", char(ExcelArr{n,'Cohort'}), " Phase ", num2str(ExcelArr{n,'Phase'}));
      	title(figName, 'Interpreter', 'none')
        
        subplot(2,1,2);
     	image(startFrame);
        axis image 
    
        % Draw boundaries
        axes(ax1)
        while ~success(1)
            uiwait(msgbox('Select four points defining the borders (double-click the last point)', ...
                'Select Maze Boundary', 'modal'));
            drawnow
            [x,y] = getpts(ax1);
            if length(x) == 4
               	success(1) = 1;
                geom.Edges{n} = [x,y];
                geom.Ctr{n}   = mean(geom.Edges{n});
                geom.Perim{n} = (0.8 .* (geom.Edges{n} - repmat(geom.Ctr{n}, size(geom.Edges{n},1),1))) ...
                    + repmat(geom.Ctr{n}, size(geom.Edges{n},1),1);
                geom.DspEdges{n} = [geom.Edges{n}; geom.Edges{n}(1,:)];
                geom.DspPerim{n} = [geom.Perim{n}; geom.Perim{n}(1,:)];
            	plot(geom.DspEdges{n}(:,1), geom.DspEdges{n}(:,2), 'y-');
                plot(geom.DspPerim{n}(:,1), geom.DspPerim{n}(:,2), 'c-');
                plot(geom.Ctr{n}(1), geom.Ctr{n}(2), 'r+');
            	drawnow
            end
        end
    
        % Select mouse center
        while ~success(2)
            uiwait(msgbox('Select the center of the mouse with a double-click', 'Select Mouse', 'modal'));
            drawnow
            [x,y] = getpts(ax1);
            if length(x) == 1
               	success(2) = 1;
                geom.MouseCtr{n} = [x,y];
            end
        end
    
        % Select toy locations if trial 2 or 3
        if ExcelArr{n,'Phase'} >= 2
            while ~success(3)
                uiwait(msgbox('First click on Object A (left/bottom-hand), then double-click on Object B (right/top-hand)', 'Select Objects', 'modal'));
                [x,y] = getpts(ax1);
                if ExcelArr.VideoOrientation == "Side"
                    if length(x) == 2 && y(1) > y(2)        % Since y-axis counts from top to bottom
                        success(3) = 1;
                        geom.ObjPos{n} = [x y];
                    end
                elseif ExcelArr.VideoOrientation == "Up"
                    if length(x) == 2 && x(2) > x(1)        % Since x-axis counts from left to right
                        success(3) = 1;
                        geom.ObjPos{n} = [x y];
                    end
                end
            end
        end

%% Saves off newly set coordinates to Excel if toggle is on
        if saveCoords == 1
            fullExcel{MatlabRows(n), 'EdgesX1'} = geom.Edges{n}(1,1);
            fullExcel{MatlabRows(n), 'EdgesY1'} = geom.Edges{n}(1,2);
            fullExcel{MatlabRows(n), 'EdgesX2'} = geom.Edges{n}(2,1);
            fullExcel{MatlabRows(n), 'EdgesY2'} = geom.Edges{n}(2,2);
            fullExcel{MatlabRows(n), 'EdgesX3'} = geom.Edges{n}(3,1);
            fullExcel{MatlabRows(n), 'EdgesY3'} = geom.Edges{n}(3,2);
            fullExcel{MatlabRows(n), 'EdgesX4'} = geom.Edges{n}(4,1);
            fullExcel{MatlabRows(n), 'EdgesY4'} = geom.Edges{n}(4,2);
            
            fullExcel{MatlabRows(n), 'MouseX'} = geom.MouseCtr{n}(1,1);
            fullExcel{MatlabRows(n), 'MouseY'} = geom.MouseCtr{n}(1,2);
            
            if ExcelArr{n,'Phase'} >= 2
                fullExcel{MatlabRows(n), 'ObjAX'} = geom.ObjPos{n}(1,1);
                fullExcel{MatlabRows(n), 'ObjAY'} = geom.ObjPos{n}(1,2);            
                fullExcel{MatlabRows(n), 'ObjBX'} = geom.ObjPos{n}(2,1);
                fullExcel{MatlabRows(n), 'ObjBY'} = geom.ObjPos{n}(2,2);
            end
        cd(mainDir)
        writetable(fullExcel, tableName)
        close all
        end
    end
end
