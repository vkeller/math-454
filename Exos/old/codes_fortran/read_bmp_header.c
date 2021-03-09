#include <stdio.h>
#include <stdlib.h>
#include <string.h>
 
struct BitMap
{
  short Type;
  long Size;
  short Reserve1;
  short Reserve2;
  long OffBits;
  long biSize;
  long biWidth;
  long biHeight;
  short biPlanes;
  short biBitCount;
  long biCompression;
  long biSizeImage;
  long biXPelsPerMeter;
  long biYPelsPerMeter;
  long biClrUsed;
  long biClrImportant;
} Header;
 
int main( void )
{
  FILE *BMPFile = fopen ("out_110.bmp", "rb");
 
  if (BMPFile == NULL)
    return;
   
  memset(&Header, 0, sizeof(Header));
   
  fread(&Header.Type, 2, 1, BMPFile);
  fread(&Header.Size, 4, 1, BMPFile);
  fread(&Header.Reserve1, 2, 1, BMPFile);
  fread(&Header.Reserve2, 2, 1, BMPFile);
  fread(&Header.OffBits, 4, 1, BMPFile);
  fread(&Header.biSize, 4, 1, BMPFile);
  fread(&Header.biWidth, 4, 1, BMPFile);
  fread(&Header.biHeight, 4, 1, BMPFile);
  fread(&Header.biPlanes, 2, 1, BMPFile);
  fread(&Header.biBitCount, 2, 1, BMPFile);
  fread(&Header.biCompression, 4, 1, BMPFile);
  fread(&Header.biSizeImage, 4, 1, BMPFile);
  fread(&Header.biXPelsPerMeter, 4, 1, BMPFile);
  fread(&Header.biYPelsPerMeter, 4, 1, BMPFile);
  fread(&Header.biClrUsed, 4, 1, BMPFile);
  fread(&Header.biClrImportant, 4, 1, BMPFile);
               
  printf("\nType:%hd\n", Header.Type);
  printf("Size:%ld\n", Header.Size);  
  printf("Reserve1:%hd\n", Header.Reserve1);  
  printf("Reserve2:%hd\n", Header.Reserve2);
  printf("OffBits:%ld\n", Header.OffBits); 
  printf("biSize:%ld\n", Header.biSize);      
  printf("Width:%ld\n", Header.biWidth);    
  printf("Height:%ld\n", Header.biHeight);  
  printf("biPlanes:%hd\n", Header.biPlanes);  
  printf("biBitCount:%hd\n", Header.biBitCount);
  printf("biCompression:%ld\n", Header.biCompression);  
  printf("biSizeImage:%ld\n", Header.biSizeImage);  
  printf("biXPelsPerMeter:%ld\n", Header.biXPelsPerMeter);  
  printf("biYPelsPerMeter:%ld\n", Header.biYPelsPerMeter);  
  printf("biClrUsed:%ld\n", Header.biClrUsed);  
  printf("biClrImportant:%ld\n\n", Header.biClrImportant);  
   
  fclose(BMPFile);
   
  return 0;
}
