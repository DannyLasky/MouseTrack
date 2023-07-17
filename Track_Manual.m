%% Settings
onMac       = 0;
ExcelRows   = 507;
setFrames   = (100:3:1098)';     % Must be a vertical vector! (do every X if using every X, 1:X:2700
setFramesExp = (100:1098)';
everyX = 3;

tableName = "PatSep Trial Data Script 04-24-23.xlsx";

if onMac == 1
    mainDir     = "/Volumes/mathewjones/Jones_Maganti_Shared/Mouse Tracking";
    outputDir   = "/Volumes/mathewjones/Jones_Maganti_Shared/Mouse Tracking/Manual Mouse Coords";
else
    mainDir     = "P:\Jones_Maganti_Shared\Mouse Tracking";
    outputDir   = "P:\Jones_Maganti_Shared\Mouse Tracking\Manual Mouse Coords";
end

MatlabRows = ExcelRows - 1;     % Offset by 1 row since 1st row becomes the table header in Matlab

cd(mainDir)
fullExcel = readtable(tableName);
ExcelArr = fullExcel(MatlabRows,:);

fileID  = "Coh " + ExcelArr{1, 'Cohort'} + " Animal " + ExcelArr{1, 'AnimalID'} + " Wk " + ExcelArr{1, 'Week'} + " Phase " + ExcelArr{1, 'Phase'};

%% Prepare variables
if onMac == 1
    tempConvert = char(ExcelArr{1,'VideoPath'});
    tempConvert = strrep(tempConvert, 'P:\Jones_Maganti_Shared', '/Volumes/mathewjones/Jones_Maganti_Shared');
    videoPathMac = strrep(tempConvert, '\', '/');
    tempObj = videoPathMac;
else
    tempObj = char(ExcelArr.VideoPath);
end

tempObj = strrep(tempObj,'.mp4','_Red.avi');
vidObj  = VideoReader(tempObj);

startFrame      = ExcelArr.StartFrame;
setPath         = zeros(length(setFrames), 2);

%% Prepare frame by frame figure
figure;
set(gcf, 'Units', 'Normalized', 'OuterPosition', [0.025 0.025 0.95 0.95], 'Visible', 'on')



for n = 1:length(setFrames)
    frame = read(vidObj, setFrames(n) + startFrame - 1);
    frame = rgb2gray(frame);
    frameAdj = imadjust(frame);
    image(frameAdj);
    axis image
    axis on
    title(setFrames(n) + startFrame - 1)
    success = 0;
    while success == 0
        drawnow
        mousePoint = drawpoint;
        if length(mousePoint) == 1
            success = 1;
            setPath(n, 1) = mousePoint.Position(1);
            setPath(n, 2) = mousePoint.Position(2);
        end
    end
end

if ~isempty(everyX)
    xPath = cell(900, 1);
    yPath = cell(900, 1);
    
    for n = 1:length(setPath) - 1
        tempX = linspace(setPath(n, 1), setPath(n + 1, 1), everyX + 1);
        tempY = linspace(setPath(n, 2), setPath(n + 1, 2), everyX + 1);
        xPath{n}(1:3,1) = tempX(1:3);
        yPath{n}(1:3,1) = tempY(1:3);
    end
    
    xPathExp = cell2mat(xPath);
    yPathExp = cell2mat(yPath);
    setPath = [xPathExp, yPathExp];
end

manual.Path = [setFramesExp, setPath];

cd(outputDir)
save(fileID, 'manual')
close all


