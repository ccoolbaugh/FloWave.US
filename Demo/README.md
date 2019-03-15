# Testing the FloWave.US Code  
*FloWave.US* was originally developed to analyze duplex ultrasound (GE Logiq Book e) recorded during muscle contraction experiments. Since its development, the code has been tested on other ultrasound equipment and applied to different experimental conditions including flow mediated dilation studies.    

We have provided examples of different ultrasound recordings to help you get started. If you are using *FloWave.US* for the first time, we recommend beginning with the muscle contraction video as you work through the [wiki](https://github.com/ccoolbaugh/FloWave.US/wiki) or [video](http://www.youtube.com/watch?feature=player_embedded&v=lehANYDmxTY "Youtube Tutorial") tutorials. All videos can also be used to practice creating ultrasound screen setting files with *UsSetup.m*.   

## Example Ultrasound Videos
The following video files are included in the "Demo" folder.

### Resting Blood Flow
**Experiment:** duplex ultrasound data were recorded in real-time in the brachial artery at rest.   
**Ultrasound Equipment:** General Electric (GE) LOGIQ Book e, 8-MHz linear array vascular probe (model 8L-RS).  
**Video Files:**   
  1. **DemoBMode.mov** - B-mode ultrasound data of the brachial artery  
    *  Note: This file can be used to test the *BMode* vessel diameter measurement code.     
  2. **GELogiqBookeDefault.mat** - ultrasound screen settings file used to calibrate *FloWave.US* to the on-screen markings.  
  
### Muscle Contraction Blood Flow  
**Experiment:** duplex ultrasound data were recorded in real-time in the anterior tibial artery during a 1-s, maximal effort dorsiflexion task.        
**Ultrasound Equipment:** General Electric (GE) LOGIQ Book e, 8-MHz linear array vascular probe (model 8L-RS)    
**Video Files:**   
  1. **DemoMuscleContraction.mov** - duplex ultrasound data with an automated calculation of TAMean velocity (teal line overlay)  
    *  Note: Analysis of this file is demonstrated in the wiki and video tutorials.     
  2. **DemoMuscleContraction_NoTAMean** - duplex ultrasound data without TAMean calculation.  
    *  Note: This video recording was obtained from a different research participant and will not be comparable to the TAMean data.   
  3. **GELogiqBookeDefault.mat** - ultrasound screen settings file used to calibrate *FloWave.US* to the on-screen markings.   
  
### Flow Mediated Dilation  
**Experiment:** duplex ultrasound data were recorded in the superficial femoral artery following 5 minutes of cuff occlusion.  
**Ultrasound Equipment:** General Electric (GE) LOGIQ Book e, 12L-RS wide-band linear array probe, 5-13 mHz  
**Credit:** Many thanks to Megan C Nelson at the University of Idaho for sharing this video.   
**Video Files:** 
  1. **Demo_FMD.mov** - duplex ultrasound data acquired during a flow mediated dilation (FMD) experiment. 
  2. **FMDtest_Settings.mat** - ultrasound screen settings file used to calibrate *FloWave.US* to the on-screen markings.

### Want to Help Grow the FloWave.US Community? 
If you've tested *FloWave.US* using different experimental conditions or ultrasound equipment, please let us know. We would love to hear how you are using the code and if you would be interested in sharing sample data with the community.   


## Getting Started with MATLAB
If you have not used MATLAB before, we recommend reviewing some of the links and short tutorial videos on the [Getting Started with MATLAB](http://www.mathworks.com/help/matlab/getting-started-with-matlab.html "MATLAB Help") to orient to the MATLAB environment.  

### Important Topics
- [Getting Started Video](https://www.mathworks.com/videos/getting-started-with-matlab-101684.html?s_tid=getstart_gs_vid "Video")  
- [Desktop Basics](https://www.mathworks.com/help/matlab/learn_matlab/desktop.html "Desktop")  
- [Workspace Variables](https://www.mathworks.com/help/matlab/learn_matlab/workspace.html "Variables")  
- [Entering Commands](https://www.mathworks.com/help/matlab/entering-commands.html "Commands")  

### Problems Opening FloWave.US?
MATLAB uses a search path to locate files. For FloWave.US to work correctly, all of the source code must be added to this search path. If a folder, subfolder, or file is not on the search path, MATLAB will return an error message ... "undefined function or variable XX".    

You can check to see if a folder is on the search path by looking at the "Current Folder" panel in the MATLAB desktop. Folders on the path will be in **bold** text. Folders *not* on the path will be in gray text. You can right-click on a folder in the "Current Folder" panel and use the menu option "Add to Path" to add the folder and subfolders to the search path, but note, these selected folders will resset when MATLAB is closed.  


## Need Help?
If you have problems running the code with these videos, please submit an [issue](https://github.com/ccoolbaugh/FloWave.US/issues "Bug Reports") on GitHub so we can help. 

Creating an issue will help us create a public forum of commonly asked questions that can be answered and viewed by all members of the *FloWave.US* community.  

