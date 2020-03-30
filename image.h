#pragma once

#include<stdio.h>
#include<stdlib.h>

typedef struct _Image {
  unsigned char*data;
  int W,H;
} Image;

#define IElem(_im,_x,_y,_c) (_im) -> data[(_y)*(_im)->W*3 +(_x)*3+(_c)]
#define isInsideImage(is,u,v) (0<=u && u<is->W && 0<=v && v<is->H)

Image*ImageAlloc(int W,int H);
Image*ImageRead(const char*name);
void ImageFree(Image*im);
void ImageWrite(const char*name,Image*im);
void ImageClear(Image*im);
void ImageImageProjectionAlpha(Image*id,Image*is,double a[3][3],double alpha);
void ImageDrawBox(Image*im,int x,int y);
void mult33(double a[3][3], double b[3][3], double c[3][3]);
