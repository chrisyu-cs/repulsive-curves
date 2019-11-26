ffmpeg -f image2 -framerate 25 -i build/frames/frame%04d.png -vcodec libx264 output.mp4
