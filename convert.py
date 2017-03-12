#!/usr/bin/python

import sys
from PIL import Image

def gray_to_rgb(input_file, output_file):
	
	N = 750
	M = 100

	img = Image.new('RGB', (int(N), int(M)), 'white')

	f = open(input_file, 'r')
	data = f.read()
	f.close()

	#print len(data)
	for i in range(0, M):
		for j in range(0, N):
	 		x = ord(data[i*N+j])
	 		#print x
	 		img.putpixel((j, i), (x, x, x))
	 		
	img.save(output_file, 'JPEG', quality = 100)

input_file = 'sigma_real_data.bin'
output_file = 'output.jpg'

if len(sys.argv) > 1 : 
	input_file =  sys.argv[1]
if len(sys.argv) > 2 : 
	output_file =  sys.argv[2]

gray_to_rgb(input_file, output_file)