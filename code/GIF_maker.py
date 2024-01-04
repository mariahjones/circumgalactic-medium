import matplotlib.pyplot as plt
import os,imageio
import numpy as np

def gif_maker(gif_name,png_dir,gif_indx,num_gifs,dpi=90):
    # make png path if it doesn't exist already
    #if not os.path.exists(png_dir):
       # os.makedirs(png_dir)
	"""
	okay, i am trying to load in the files, save them with plt and append on the ith number
	(so that they are sorted in order), and then make the gif. but i cant figure out how to
	do this. probably like gif_indx.append(f)
		- i have since put this in, and will try this at some later date
	"""
	for f in range(6) 	
		png_dir = np.load(file_name) #file_name is wherever the downloaded .pngs are
# save each .png for GIF
    # lower dpi gives a smaller, grainier GIF; higher dpi gives larger, clearer GIF
    	plt.savefig(png_dir+'frame_'+str(gif_indx)+'_.png',dpi=dpi)
    #plt.close('all') # comment this out if you're just updating the x,y data
    	gif_indx = gif_indx.append(f)
	    if gif_indx==num_gifs-1:
	        # sort the .png files based on index used above
	        images,image_file_names = [],[]
	        for file_name in os.listdir(png_dir):
	            if file_name.endswith('.png'):
	                image_file_names.append(file_name)       
	        sorted_files = sorted(image_file_names, key=lambda y: int(y.split('_')[1]))

	   # define some GIF parameters
	        
	        frame_length = 0.5 # seconds between frames
	        end_pause = 1 # seconds to stay on last frame
	        # loop through files, join them to image array, and write to GIF called 'wind_turbine_dist.gif'
	        for ii in range(0,len(sorted_files)):       
	            file_path = os.path.join(png_dir, sorted_files[ii])
	            if ii==len(sorted_files)-1:
	                for jj in range(0,int(end_pause/frame_length)):
	                    images.append(imageio.imread(file_path))
	            else:
	                images.append(imageio.imread(file_path))
	        # the duration is the time spent on each image (1/duration is frame rate)
	        imageio.mimsave(gif_name, images,'GIF',duration=frame_length)