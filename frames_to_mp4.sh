ffmpeg -f image2 -i build/frames/frame%04d.png -r 24 -vcodec libx264 output.mp4
