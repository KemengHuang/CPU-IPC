import cv2
import os
from tqdm import tqdm
fsp = 30
fourcc = cv2.VideoWriter_fourcc('D', 'I', 'V', 'X')
 
 
n = 0
video_path = 'video.mp4'
img_path = "saveScreen2"
img_path2 = "saveScreen1"
list_image = os.listdir(img_path)
list_image.sort()
list_image2 = os.listdir(img_path2)
list_image2.sort()
 
list_image = [os.path.join(img_path,x) for x in list_image]
list_image2 = [os.path.join(img_path2,x) for x in list_image2]
width = cv2.imread(list_image[0]).shape[1]
heighth = cv2.imread(list_image[0]).shape[0]
width = width * 2
video_out = cv2.VideoWriter(video_path, fourcc, fsp, (width,heighth))
print(len(list_image))
 
count = 0
for i in tqdm(range(len(list_image))):
     if i == 0:
          continue
     if os.path.exists(list_image[i]) and os.path.exists(list_image2[i]):
          img1 = cv2.imread(list_image[i])
          img2 = cv2.imread(list_image2[i])
          frame = cv2.hconcat([img1, img2])
          video_out.write(frame)
          count += 1
 
print('cout',count)
 
video_out.release()