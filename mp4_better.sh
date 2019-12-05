ffmpeg -i build/frames/frame%04d.png -c:v libx264 -crf 18 -framerate 60 -profile:v main -pix_fmt yuv420p -c:a aac -ac 2 -b:a 128k -movflags faststart output.mp4
