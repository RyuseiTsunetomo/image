/*

 gcc -O ransac2.c image.c -lm ; ./a.out 0.jpg 1.jpg
*/
#include<stdio.h>
#include<math.h>
#include<stdlib.h>
#include<string.h>
#include"image.h"

#define Elem(_a,_b,_c)  (_a)->data[(_a)->W*(_b)+(_c)]
#define Row(_a,_b)     ((_a)->data+(_a)->W*(_b))
#define Elem(_a,_b,_c)  (_a)->data[(_a)->W*(_b)+(_c)]
#define MAX 100
#define __rdtsc() ({ long long a,d; asm volatile ("rdtsc":"=a"(a),"=d"(d)); d<<32|a; })
#define DElem(_a,_b,_c)  (_a)->data[(_a)->W*(_c)+(_b)]

//第3回
double VP(double*a,double*b,int N){
  double s=0;
  int i;
  for(i=0;i<N;i++) s += a[i] * b[i] ;
  return s;
}

void VSS(double*d,double s,int N){
  int i;
  for(i=0;i<N;i++) d[i] *= s;
}

void VSA(double*d,double*a,double s,int N){
  int i;
  for(i=0;i<N;i++) d[i] += a[i] * s;
}
// matrix structure and functions
typedef struct {
  double *data;
  int W,H;
} Matrix;


void MatrixClear(Matrix*mt){
  memset(mt->data,0,mt->W*mt->H*sizeof(double));
}


void MatrixCopy(Matrix*mtD,Matrix*mt){
  memmove(mtD->data,mt->data,mt->W*mt->H*sizeof(double));
}

void MatrixCopyT(Matrix*mtD,Matrix*mt){
  int i,j;
  for(i=0;i<mtD->H;i++)
    for(j=0;j<mtD->W;j++)
      Elem(mtD,i,j) = Elem(mt,j,i);
}

void MatrixPrint(Matrix*mt){
  int i,j;
  for(i=0;i<mt->H;i++){
    for(j=0;j<mt->W;j++)
      printf("%f ",Elem(mt,i,j));
    printf("\n");
  }
}

void MatrixMultT(Matrix*mtD,Matrix*mtA,Matrix*mtB){
  // D = A B^T
  int i,j;
  for(i=0;i<mtA->H;i++)
    for(j=0;j<mtB->H;j++)
      Elem(mtD,i,j) = VP( Row(mtA,i), Row(mtB,j), mtA->W);
}


void MatrixQRDecompColMajor(Matrix*mtR,Matrix*mt){
  // Gram-Schmidt orthonormalization (R and Q)
    double t, *aT[] = { Row(mt,0), Row(mt,1), Row(mt,2), Row(mt,3), Row(mt,4),Row(mt,5), Row(mt,6), Row(mt,7)} ;
  int W = mt->W;
  MatrixClear(mtR);
  int i, j;

  for(i=0; i<8; i++){
      for(j=0; j<=i; j++){
          if(i > j){
              Elem(mtR,j,i) = t = VP(aT[j], aT[i], W);
              VSA(aT[i], aT[j], -t, W);
          }
          else{
              Elem(mtR,j,i) = t = sqrt(VP(aT[j], aT[i], W));
              VSS(aT[i], 1/t, W);
          }
      }
  }
}

void MatrixSimeqLr(Matrix*mtB,Matrix*mtR){
  // B = B L^{-1}
    long i, j;
    double * B = Row(mtB,0);
  
    for(i=7; i>=0; i--){
        for(j=7; i<j; j--){
            B[i] -= B[j]*Elem(mtR,i,j);
        }
        B[i] = B[i] / Elem(mtR,i,i);
    }
}

void Assign_Mat(double m[][3], Matrix*ans){
    m[0][0]=Elem(ans,0,0);
    m[0][1]=Elem(ans,0,1);
    m[0][2]=Elem(ans,0,2);
    m[1][0]=Elem(ans,0,3);
    m[1][1]=Elem(ans,0,4);
    m[1][2]=Elem(ans,0,5);
    m[2][0]=Elem(ans,0,6);
    m[2][1]=Elem(ans,0,7);
    m[2][2]=1;
}

void lsq(double m10[][3], int w[][4], 
         Matrix*cmA, Matrix*vt, Matrix*mtR, Matrix*tmp){
  int i;

  // create A (col-major)
  for(i=0;i<4;i++){
    Elem(cmA,0,i*2  )=w[i][0];
    Elem(cmA,1,i*2  )=w[i][1];
    Elem(cmA,2,i*2  )=1;
    Elem(cmA,3,i*2  )=0;
    Elem(cmA,4,i*2  )=0;
    Elem(cmA,5,i*2  )=0;
    Elem(cmA,6,i*2  )=-w[i][0]*w[i][2];
    Elem(cmA,7,i*2  )=-w[i][1]*w[i][2];
    Elem(cmA,0,i*2+1)=0;
    Elem(cmA,1,i*2+1)=0;
    Elem(cmA,2,i*2+1)=0;
    Elem(cmA,3,i*2+1)=w[i][0];
    Elem(cmA,4,i*2+1)=w[i][1];
    Elem(cmA,5,i*2+1)=1;
    Elem(cmA,6,i*2+1)=-w[i][0]*w[i][3];
    Elem(cmA,7,i*2+1)=-w[i][1]*w[i][3];
    Elem(vt ,0,i*2  )=w[i][2];
    Elem(vt ,0,i*2+1)=w[i][3];
  }  
  MatrixQRDecompColMajor(mtR,cmA);
  MatrixMultT(tmp,vt,cmA);
  MatrixSimeqLr(tmp,mtR);
  Assign_Mat(m10, tmp);
  
}

Matrix*MatrixAlloc(int _H,int _W){
    Matrix*mt=(Matrix*)malloc(sizeof(Matrix));;
    mt->W = _W,
        mt->H = _H;
    mt->data=(double*)malloc(mt->W*mt->H*sizeof(double));
    return mt;
}
//  第5回
void Imagefeature(Matrix*im2,Image*im){
    int x,y,u,v,W=7,ix,iy;
    double a;
    for(y=W+1;y<im->H-W-1;y++) for(x=W+1;x<im->W-W-1;x++){
            double ixx,ixy,iyy;
            ixx=iyy=ixy=0;
            for(v=-W;v<=W;v++) for(u=-W;u<=W;u++){
                    ix=IElem(im, x+u+1, y+v, 1) - IElem(im, x+u-1, y+v, 1);
                    iy=IElem(im, x+u, y+v+1, 1) - IElem(im, x+u, y+v-1, 1);
                    ixx+=ix*ix; 
                    ixy+=ix*iy;
                    iyy+=iy*iy;
                }
            a=((ixx+iyy)*(ixx+iyy)-4*(ixx*iyy-ixy*ixy));
            a=sqrt(a);
            DElem(im2, x, y) = (((ixx+iyy)-a)/2);
        }
}
int MatrixLocalMax(int w[][2], Matrix*im2, int M){
    int x,y,u,v,W=7,n=0,a, i;
    for(y=W+1;y<im2->H-W-1;y++) for(x=W+1;x<im2->W-W-1;x++){
            double max=-1;
            for(v=-W;v<=W;v++) for(u=-W;u<=W;u++){
                    if(max < DElem(im2, x+u, y+v)){
                        max = DElem(im2, x+u, y+v);
                    } 
                }
            if(max == DElem(im2, x, y)){
                a=n; if(n<M) n++;
                for(;a>0 && DElem(im2,w[a-1][0],w[a-1][1]) 
                        < DElem(im2,x,y); a--){
                    w[a][0]=w[a-1][0];
                    w[a][1]=w[a-1][1];
                } 
                w[a][0]=x; 
                w[a][1]=y;
            }
        }
    /*  for(i=0;i<M;i++){
          printf("%d:x=%d, y=%d, %lf\n", i, w[i][0], w[i][1], DElem(im2,w[i][0],w[i][1]));
        } */
    return n; // 記録した点の数
}

// 第6回
double ImageSSD(Image*im,int x1,int y1, Image*im2,int x2,int y2){
  int i,j,W=7;
  double sr=0,sg=0,sb=0,dr,dg,db;
  for(i=-W;i<=W;i++) for(j=-W;j<=W;j++){
    dr  = IElem(im, x1+j, y1+i, 0) - IElem(im2, x2+j , y2+i, 0);
    dg  = IElem(im, x1+j, y1+i, 1) - IElem(im2, x2+j , y2+i, 1);
    db  = IElem(im, x1+j, y1+i, 2) - IElem(im2, x2+j , y2+i, 2);
    sr += dr*dr;
    sg += dg*dg;
    sb += db*db;
  }
  return sr+sg+sb;
}


void calcSSDtable(Matrix*mt,
		  Image*im ,int x1[][2],int N1,
		  Image*im2,int x2[][2],int N2){
  int i,j;
  for(i=0;i<N1;i++)
    for(j=0;j<N2;j++)
      Elem(mt,i,j) = ImageSSD(im ,x1[i][0],x1[i][1],
			      im2,x2[j][0],x2[j][1]);
}

/*int matchMethod1(double w[][4],Matrix*mt,
                 Image*im ,int x1[][2],int N1,
                 Image*im2,int x2[][2],int N2){
  int i,j,k,ji,n=0;

  for(i=0;i<N1;i++){
    double sm=INFINITY,t;
    for(j=0;j<N2;j++){
      t=Elem(mt,i,j);
      if(sm>t) sm=t, ji=j;
    }
    printf("%d,%d,%d,%d,\n",
           x1[i][0],x1[i][1],
           x2[ji][0],x2[ji][1]);
    // 上の printf で表示されるものを w[n][] に格納.
    w[i][1] = x1[i][0];
    w[i][2] = x1[i][1];
    w[i][3] = x2[ji][0];
    w[i][4] = x2[ji][1];
    n++;
    for(k=0;k<N1;k++) Elem(mt,k,ji)=INFINITY;
  }

  return n;
}
*/

int matchMethod2(int w[][4],Matrix*mt,
                 Image*im ,int x1[][2],int N1,
                 Image*im2,int x2[][2],int N2){
    int i,j,k,ji,n=0,l,ii;

    for(l=0;l<N1;l++){
        double sm=INFINITY,t;
        for(i=0;i<N1;i++){
            for(j=0;j<N2;j++){
                t=Elem(mt,i,j);
                if(sm>t) sm=t, ji=j, ii=i;
            }
        }
        /*printf("%d,%d,%d,%d,\n",
               x1[ii][0],x1[ii][1],
               x2[ji][0],x2[ji][1]);*/
        // 上の printf で表示されるものを w[n][] に格納.
        w[l][0] = x1[ii][0];
        w[l][1] = x1[ii][1];
        w[l][2] = x2[ji][0];
        w[l][3] = x2[ji][1];
        n++;
        for(k=0;k<N1;k++) Elem(mt,k,ji)=INFINITY;
        for(k=0;k<N1;k++) Elem(mt,ii,k)=INFINITY;
    }

    return n;
}

//第7回
void initRndAry(int rndAry[]){
    int i;
    for(i=0;i<MAX;i++) rndAry[i]=i;
    srand48(__rdtsc()); // __rdtsc については第４回を参照．
}


void chooseFourNumbers(int rndAry[]){
    int i;
    for(i=0;i<4;i++){
        int j, t;
        j=(int)(drand48()*(MAX-i))+i; // 乱数関数は stdlib.h で宣言されている．
        t=rndAry[i]; rndAry[i]=rndAry[j]; rndAry[j]=t;
    }
}

void calcHomography(double H[3][3], int w[][4], int rndAry[4], Matrix*cmA, Matrix*vt, Matrix*mtR, Matrix*tmp){
    int a=rndAry[0], b=rndAry[1], c=rndAry[2], d=rndAry[3];

   
    int ww[][4]={
        w[a][0], w[a][1], w[a][2], w[a][3],
        w[b][0], w[b][1], w[b][2], w[b][3],    
        w[c][0], w[c][1], w[c][2], w[c][3],    
        w[d][0], w[d][1], w[d][2], w[d][3],    
    };
    lsq(H, ww, cmA, vt, mtR, tmp);
}

int calcScore(double H[3][3], int w[][4]){
    int i, score=0;
    for(i=0;i<MAX;i++){
        double x=w[i][0], y=w[i][1],
            u=w[i][2], v=w[i][3],
            x_prime, y_prime, z_prime, // 変換行列と (x,y,1)^T の積
            du,dv; // (x,y) を変換した座標(第2回の概要を参照)と (u,v) の差．
        x_prime = H[0][0]*x + H[0][1]*y + H[0][2];
        y_prime = H[1][0]*x + H[1][1]*y + H[1][2];
        z_prime = H[2][0]*x + H[2][1]*y + H[2][2];
        du = u - (x_prime/z_prime);
        dv = v - (y_prime/z_prime);
        if(du*du+dv*dv<3) {
            score++;
        }
    }
    return score;
}
void createPanorama(Image*id, Image*im, Image*im2, double bestH[3][3]){
    double m0d[][3]={
        1,0,-100,
        0,1,-100,
        0,0,1
    }; 
    int i;
    double m1d[3][3];
    ImageImageProjectionAlpha(id,im,m0d,.5);
    mult33(m1d, bestH, m0d);
    ImageImageProjectionAlpha(id,im2,m1d,.5);
/*    for(i=0; i<3;i++){
        printf("%f, %f, %f\n", m1d[i][0], m1d[i][1], m1d[i][2]);
    }
*/
    ImageWrite("out3.ppm", id);
}

main(int ac,char**av){
    Image *im,*im2, *id;
  Matrix*mt;
  Matrix*imf;
  int i;
  int w[MAX][4];
  int nm;
  int TRIAL=1000;
  
  int x1[MAX][2], N1;
  int x2[MAX][2], N2;
  

  if(ac<3) return 1;

  im =ImageRead(av[1]);
  imf=MatrixAlloc(im->H,im->W);
  Imagefeature(imf, im);
  N1=MatrixLocalMax(x1, imf, MAX);

  im2=ImageRead(av[2]);
  Imagefeature(imf, im2);
  N2=MatrixLocalMax(x2, imf, MAX);

  mt=MatrixAlloc(N1,N2);
  calcSSDtable(mt,im,x1,N1,im2,x2,N2);
  //nm=matchMethod1(w,mt,im,x1,N1,im2,x2,N2); // 特徴点の対応付け
  nm=matchMethod2(w,mt,im,x1,N1,im2,x2,N2); // 特徴点の対応付け
  /*for(i=0; i<nm; i++){
        printf("%d,%d,%d,%d,\n",
        w[i][0],w[i][1],
        w[i][2],w[i][3]);
        }
  */
  double bestH[3][3]; int trial, best_score=0 ; // 既に発見された最良の変換行列と得点
  int rndAry[MAX], score=0;
  Matrix*cmA,*vt,*mtR,*tmp; // 作業領域をここで確保して calcHomography で使う
  cmA=MatrixAlloc(8,8);
  vt=MatrixAlloc(1,8);
  mtR=MatrixAlloc(8,8);
  tmp=MatrixAlloc(1,8);
  initRndAry(rndAry);
  for(trial=0;trial<TRIAL;trial++){
      double H[3][3]; // 選んだ4点から計算される変換行列の置き場
      chooseFourNumbers(rndAry);
      calcHomography(H,w,rndAry, cmA,vt,mtR,tmp);
      score=calcScore(H,w);
      if(best_score < score){
              for(i=0; i<3; i++){
                  bestH[i][0] = H[i][0];
                  bestH[i][1] = H[i][1];
                  bestH[i][2] = H[i][2];
              }
              best_score=score;
      }
  }
  for(i=0; i<3;i++){
  printf("%f, %f, %f\n", bestH[i][0], bestH[i][1], bestH[i][2]);
  }
  printf("best_score=%d\n", best_score);
  id=ImageAlloc(1024,768);
  ImageClear(id);
  createPanorama(id, im,im2,bestH);
}
