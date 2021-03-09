#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <float.h>


/*
From : http://stackoverflow.com/questions/2654480/writing-bmp-image-in-pure-c-c-without-other-libraries
*/

int write_to_file_binary(int n, float *tmptab, int iter, float minval, float maxval){
	unsigned char *img = NULL;
	int i,j,x,y,w,h;
	float v,r,g,b;
	char *filename;
	char prefix[] = "out_";
	char suffix[] = ".bmp";
	char number[5];
	FILE *f;

	w = n; // width
	h = n; // height

	int filesize = 54 + 3*w*h;  //w is your image width, h is image height, both int
	if( img ) free( img );
	img = (unsigned char *)malloc(3*w*h);
	memset(img,0,sizeof(img));

	for(i=0; i<w; i++){
		for(j=0; j<h; j++){
			x=i; y=(w-1)-j;
			v = ((TMPTAB(j,i)-minval)/(maxval - minval));
			r = v*255; // Red channel
			g = v*255; // Green channel
			b = v*255; // Red channel
			if (r > 255) r=255;
			if (g > 255) g=255;
			if (b > 255) b=255;
			img[(x+y*w)*3+2] = (unsigned char)(r);
			img[(x+y*w)*3+1] = (unsigned char)(g);
			img[(x+y*w)*3+0] = (unsigned char)(b);
		}
	}

	unsigned char bmpfileheader[14] = {'B','M', 0,0,0,0, 0,0, 0,0, 54,0,0,0};
	unsigned char bmpinfoheader[40] = {40,0,0,0, 0,0,0,0, 0,0,0,0, 1,0, 24,0};
	unsigned char bmppad[3] = {0,0,0};

	bmpfileheader[ 2] = (unsigned char)(filesize    );
	bmpfileheader[ 3] = (unsigned char)(filesize>> 8);
	bmpfileheader[ 4] = (unsigned char)(filesize>>16);
	bmpfileheader[ 5] = (unsigned char)(filesize>>24);

	bmpinfoheader[ 4] = (unsigned char)(       w    );
	bmpinfoheader[ 5] = (unsigned char)(       w>> 8);
	bmpinfoheader[ 6] = (unsigned char)(       w>>16);
	bmpinfoheader[ 7] = (unsigned char)(       w>>24);
	bmpinfoheader[ 8] = (unsigned char)(       h    );
	bmpinfoheader[ 9] = (unsigned char)(       h>> 8);
	bmpinfoheader[10] = (unsigned char)(       h>>16);
	bmpinfoheader[11] = (unsigned char)(       h>>24);

	filename = (char *) malloc(sizeof(char)*(8+4+1));
	strcpy(filename,prefix);
	sprintf(number, "%d", iter);
	while(iter<1000){
		strcat(filename,"0");
		if(iter == 0) iter=1;
		iter *= 10;
	}	
	strcat(filename,number);
	strcat(filename,suffix);
	f = fopen(filename, "w");

	fwrite(bmpfileheader,1,14,f);
	fwrite(bmpinfoheader,1,40,f);
	for(i=0; i<h; i++){
		fwrite(img+(w*(h-i-1)*3),3,w,f);
		fwrite(bmppad,1,(4-(w*3)%4)%4,f);
	}
	fclose(f);

}


