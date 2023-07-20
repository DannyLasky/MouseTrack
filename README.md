# MouseTrack
## Mouse tracking code developed in the Jones and Maganti Labs

Hello! Welcome to the README for the mouse tracking scripts. This will explain the expectations for data being read into these scripts
and how to use them correcty. Feel free to reach me at djlasky@ucdavis.edu or dannyjlasky@gmail.com if you require necessary assistance.

Danny Lasky, 2023

**Matlab basics**
- I strongly recommend taking Matlab Onramp, a free online 2-hour course, before using the script at all
- I strongly recommend taking Matlab Fundamentals, a free online 16 hour course, before intending to make any significant changes to the script
- You need all related files in your Matlab path before running. On the "HOME" tab click "Set Path"
- Add folders/subfolders until all the code, input files, and Excel spreadsheet are in your Matlab path and save
- Path names will be different on a Mac compared to Windows, keep this in mind
- You will get errors! View these errors as challenges, not as impedences. Stack Overflow, Matlab answers, and Matlab documentation will be your best friends
- Google errors exactly as they pop up. These issues have been encountered before and are all online


**Expectations for tracking videos**
- Recorded originally in 1920 x 1080 quality and 30 fps (lower quality is fine)
- Have consistent lighting throughout video
- Use toys of a consistent size
- Keep arena directly under and all corners in view

**Pipeline for analyzing tracking videos**
- Enter the data in an Excel spreadsheet in a similar way as Vince PatSep Trial Data 06-14-23 (ask Matt Jones for this if you don't have it)
	- PatSep Trial Data 05-03-23 (Sophie's and Danny's experiments) should be almost exactly the same and can be used as a template
	- Will need to have the EXACT same variable names at the top of the Excel sheet

1. Need to downsample the video (speeds up analysis significantly and has no loss of accuracy)
  1. Open MovieResize.m
  2. Change rows 3-8 to match your requirements (between CHANGE ME)
     - These select the file path to the video, the Excel sheet, and the Excel rows to analyze in the spreadsheet
     - There is also a toggle to change the frame rate (1 is on, 0 is off) and what to set the frame rate to.
     - Shrink factor reduces each video dimension by that amount. If shrink factor = 4, then 1920 x 1080 becomes 480 x 270
  3. I recommend reducing videos to be 480 x 270 and 15 fps for faster analysis. Can view this by right clicking video > "Properties" > "Details"










**2. Select video start and end frames**
  1. In Matlab click the apps tab, open Video Viewer (under image processing and computer vision)
  2. Open your REDUCED video in Video Viewer
  3. Find the start frame to put in the Excel spreadsheet (arrow keys to go through frames, displayed in bottom right)
  4. End frame will be 2699 frames later (3 minutes at 15 fps = 2700 frames, but start and end frame both count, so subtract 1)
  5. Don't have any end frames where your hand/shadows enter into the arena. Shouldn't happen, but may have to pick an earlier end frame to avoid

**3. Preparing the main script**
  1. Open Track_Master. This script will call the functions: Track_Coords and Track_Master
  2. Set the Excel name, Excel rows to be run, and directory pathes
  3. Have readCoords = 0 if setting coordinates for the first time, otherwise have readCoords = 1 to read them in automatically
  4. You can enter directory pathes for both Mac and PC and toggle onMac accordingly (1 = on a Mac, 0 = on a PC)
  5. By default, the script saves the coordinates back to the Excel sheet, is expecting to work with a reduced video, has a object close
		threshold of 10 cm, and jump limiter of 50 (to keep the tracker on the mouse in bad situations), a wall length of 63 cm, and a
		version number of 7. Feel free to change these as necessary and update the version number if you make changes, it will saved back
		to the Excel sheet too, which helps you keep track of what version your data was run on. These variables actively affect others.
		If the wall length changes, it's important that it's changed in the script!

**4. Running the main script**
  1. CLOSE THE EXCEL SHEET (Otherwise output can not be saved there and you'll get an error)
  2. If readCoords = 0, you will be prompted to set the coordinates. If readCoords = 1 it will read them directly from the spreadsheet (assuming that they've been previously set, otherwise you'll get an error).
  3. An important note here is VideoOrientation. Look back through Sophie and Vince videos to get an idea of what the Up/Side orientations mean. If
		you have an orientation different than these you may have to rotate your videos or go into the TrackCoords script and change where it is
		used (lines 120-140). It limits where you can say the objects are by knowing what options make sense with the video orientation
	4. The script runs by creating a weighted background image (1/4 frames without a mouse, 3/4 frames with mouse throughout video and subtracting 
		that from each video frame to identify the biggest blob of difference, which should always be the mouse.
	5. You should get a heatmap showing where the mouse spent time, a trajectory path of the mouse, a double-speed, low resolution video of the tracking,
		a statistics dashboard, and a Matlab output with statistics.

**5. Creating group statisitics to compare treatments**
	1. Can do this directly through Excel or by using the Tracking_Stats_Function. This was not written by me, but works well for viewing different
		statisitcs like distance traveled and time near objects for your selected groups.

Good Luck,
Danny
