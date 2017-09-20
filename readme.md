The workflow is developed around Optovue OCT scanner. It begins with processing files and extract individual OCT images. File processing workflow might be unique to the particular machine that we develop the codes around, especially with regard to file structure.  

1. Run OCT_process_file.m to process file name and extract raw images and save in tiff format.
2. Run ImageAnalysis_for_CrossLine.m on individual cross lines scans. Detailed description of usage can be found [here](https://www.protocols.io/view/analysis-of-choroid-thickness-on-oct-image-using-m-imqcc5w)
3. After you've processed four crossline scan for both eyes, run OCT_Data_Comparison.m to generate interpolated thickness profile and comparison between the two eyes.  
4. Use OCT_Macular_Choroid_Area.m to process crossline scan through macula. 
