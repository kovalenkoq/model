#include "CImg/CImg.h"
#include <iostream>

using namespace std;
using namespace cimg_library;
int main() {

	int H = 100, N = 750;

	unsigned char color[] = { 0,0,0 };
  	CImg<unsigned char> image(N,H,1,3,0);
  	image.fill(0);

	FILE *f_in = fopen("sigma_real_data.bin", "rb");
  	if (!f_in) return 0;

  	unsigned char ch = 0;
  	int i = 0, j = 0;
  	size_t err;
  	while(!feof(f_in))
  	{	
  		err = fread(&ch, sizeof(char), 1, f_in);
  		color[0] = ch;
  		color[1] = ch;
  		color[2] = ch;

  		image.draw_point(i++,j,color);
  		if (i == N)
  		{
  			j++;
  			i = 0;
  		}
  	}

  image.save_bmp("image");
  fclose(f_in);
  return 0;
}