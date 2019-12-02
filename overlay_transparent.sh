ffmpeg -n -i $1 -i $2 -filter_complex "[1:v]format=rgba,colorchannelmixer=aa=0.5[fg];[0][fg]overlay" -movflags +faststart overlay.mp4
