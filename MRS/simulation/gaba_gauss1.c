#include <gamma.h>
#include <iostream>
#include <stdio.h>
#include <stdlib.h>

using namespace std;
int main(int argc, char *argv[])
{
  
  int i;
  int gaba_on;
  //  float j;
  int ss;
  FILE * fid;
  float * rho;
  //  float * phi;
  int rf_size;
  size_t result;
  int pw;
  int offset1;
  //  int offset2;
  char *outfile;
  if(argc!=3)
    exit(0);
  outfile=argv[1];
  offset1=atoi(argv[2]);


  cout<<argc<<endl;
  cout<<offset1<<endl;

  fid = fopen("gaba_gauss_gamma.wav","rb");

  if(fid==NULL) {fputs("File error",stderr); exit(1);}
  fseek(fid,0,SEEK_END);
  rf_size=ftell(fid);
  rewind(fid);
  cout<<"RF size is"<<rf_size<<endl;
  rho=(float *)malloc(rf_size);
  if(rho==NULL) {fputs("Memory error",stderr); exit(2);}
  result = fread(rho,sizeof(float),rf_size/sizeof(float),fid);
  cout<<"Result is"<< result<<endl;
  fclose(fid);

  pw=rf_size/sizeof(float);

  cout<<"pulse length is "<<pw<<endl;
  for(i=0;i<10;i++)
  cout<<rho[i]<<endl; 
 
  gaba_on=0;
  ss=0;

  spin_system sys;

  sys.read("gaba_gauss_on1.sys");
 
  cout<<sys<<endl;
  sys.offsetShifts(offset1);
  cout<<sys<<endl;


  float t_te=0.068;
  float t_dwell=0.000032;
  float t_12=0.006;
  float t_2g1=(t_12+t_te/2)/2-t_12-pw*t_dwell/2;
  float t_g13=t_te/2-t_2g1-pw*t_dwell;
  float t_3g2=(t_te/2-t_12)/2-pw*t_dwell/2;
  float t_g2r=(t_te/2-t_12)-t_3g2-pw*t_dwell;
 
  float angle_factor=0.24/0.24;
  cout<<"T_12 "<<t_12<<endl;
  cout<<"T_2g1 "<<t_2g1<<endl;
  cout<<"T_pulse "<<pw*t_dwell<<endl;
  cout<<"T_g13 "<<t_g13<<endl;
  cout<<"T_3g2 "<<t_3g2<<endl;
  cout<<"T_pulse "<<pw*t_dwell<<endl;

  cout<<"T_g2r "<<t_g2r<<endl;
  

  gen_op sigma0 = sigma_eq(sys);
  gen_op H = Ho(sys);
  gen_op detect = Fm(sys);


  gen_op udelay1=prop(H,t_12);
  gen_op udelay2g1=prop(H,t_2g1);
  gen_op udelayg13=prop(H,t_g13);
  gen_op udelay3g2=prop(H,t_3g2);
  gen_op udelayg2r=prop(H,t_g2r);



  gen_op udelay01=prop(H,0.005);
  gen_op udelay02=prop(H,0.034);
  gen_op udelay021=prop(H,0.017);
  gen_op udelay022=prop(H,0.017);
  gen_op udelay031=prop(H,0.017);
  gen_op udelay032=prop(H,0.012);

  gen_op udelay03=prop(H,0.029);

  gen_op udelay_dwell=prop(H,t_dwell);

 


  float angle=90;
 

	gen_op sigma1 = Iypuls(sys,sigma0,90.0);

    sigma1=evolve(sigma1,udelay1);
    sigma1=Iypuls(sys,sigma1,180.0);
    sigma1=evolve(sigma1,udelay2g1);

    /*****start mega press****/
  for(i=0;i<pw;i++)
    {
      
      angle=rho[i]*angle_factor*180/3.14159;

      sigma1=Ixypuls(sys,sigma1,0,angle);
      sigma1=evolve(sigma1,udelay_dwell);
    }
    sigma1=evolve(sigma1,udelayg13);
    sigma1=Iypuls(sys,sigma1,180.0);
    sigma1=evolve(sigma1,udelay3g2);

  for(i=0;i<pw;i++)
    {
      
      angle=rho[i]*angle_factor*180/3.14159;
      sigma1=Ixypuls(sys,sigma1,0,angle);
      sigma1=evolve(sigma1,udelay_dwell);

    }
  sigma1=evolve(sigma1,udelayg2r);


  block_1D data(2048);
  FID(sigma1,detect,H,0.0002,2048,data);

      MATLAB(outfile,"test_fid",data);

  for(i=0;i<10;i++)
    cout<<data(i)<<endl;
  

  free(rho); 

  return 0;

}
