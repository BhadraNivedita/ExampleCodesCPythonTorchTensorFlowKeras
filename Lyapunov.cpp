	////////Simulating Kapitza equation//////////

	//Compile :g++ <filename>.cpp

	//Run::./a.out

	#include <stdio.h>

	#include <iostream>

	#include <math.h>

	#include <fstream>
	
	#include<stdlib.h>

	#include<sstream>
	
	using namespace std;	

	#define DIMENSION	2

	#define ARRAY_SIZE	8

	#define Omega 		40			

	#define ITERATION	100000

	////////////////////////Defining class for pendulum/////////////////////////

	class Pendulum 

		{

		private:

		double pi,Real_Time,Time_Step_Size,t;

		public:

		void Derivatives(double,double*,double*);//This member function holds the equations to be integrated

        	Pendulum() { pi= 4.0*atan(1.0);}  //no argument constructor
 
		void Runge_Kutta4(double,double,double*,double*);/*Integrates those 8 equations by RK4 Algorithm*/ 

		void test(void); /*Initialises the eqn ,integrates and writes to the file*/

		};


	
	void Pendulum::Derivatives(double t,double* yy, double* dy)

	{

		double Amp,g,Length,beta; //Parameters for the system

		double j[DIMENSION][DIMENSION],e[DIMENSION][DIMENSION],elm[DIMENSION][DIMENSION],jem[DIMENSION][DIMENSION];

		double jlm[DIMENSION][DIMENSION],llm[DIMENSION][DIMENSION],tmp;
	
		g=9.8;		

		Length=1.0;

		Amp=0.110679718;

		beta=16.0;

	///////////////////////////Equations to be integrated///////////////////////////

		dy[0]  = yy[1];

		//dy[1] =sin(yy[0]);                                           

		dy[1]  = -(g*sin(yy[0])+Amp*Omega*Omega*sin(yy[0])*cos(Omega*t));

	//	dy[1]=-g*sin(yy[0]);

		//cout<<sin(pi/2.0);

	//	cout<<-g*sin(yy[0])-Amp*Omega*Omega*sin(yy[0])*cos(Omega*t)<<endl;

	//dy[1]=-yy[0];

		j[0][0] = 0.0;

		j[0][1] = 1.0;
	
		j[1][0] = -g*cos(yy[0])-Amp*Omega*Omega*cos(yy[0])*cos(Omega*t);

	//j[1][0]=-1.0;

		j[1][1] = 0.0;

		int k=DIMENSION;

        	for(int l=0;l<DIMENSION;l++){

           		for(int m=0;m<DIMENSION;m++){

             						 e[m][l] = yy[k];    

	//	cout<<k<<endl;
             							k=k+1;

           						}


		    
                              			}

	//cout<<yy[8]<<endl;	

	//cout<<e[1][0]*e[1][0]<<endl;

	tmp = sqrt(e[0][0]*e[0][0] + e[1][0]*e[1][0]);

	//cout<<tmp<<endl;

	//cout<<e[1][0]<<"\n";

	e[0][0] = e[0][0]/tmp;

	e[1][0] = e[1][0]/tmp;

	tmp=sqrt(e[0][0]*e[0][0]+e[1][0]*e[1][0]);/*cHECK IF THE VECTORS ARE ORTHONORMAL ~1.0 */ 

//	cout<<e[1][0]*e[1][0]<<"\n";

	//cout<<e[1][1]<<endl;

	tmp = sqrt(e[0][1]*e[0][1] + e[1][1]*e[1][1]);

	e[0][1] = e[0][1]/tmp;

	e[1][1] = e[1][1]/tmp;

	//tmp=sqrt(e[0][1]*e[0][1]+e[1][1]*e[1][1]);

	//cout<<tmp<<"\n";
	
        //cout<<(e[0][0]*e[0][1]+e[1][0]*e[1][1])<<endl;


	for( int l=0;l<DIMENSION;l++){

        	   for(int m=0;m<DIMENSION;m++){

		              elm[l][m]=0.0;

		              jem[l][m]=0.0;
              
				for( int k=0;k<DIMENSION;k++){

			                 elm[l][m] = elm[l][m] + e[k][l]*e[k][m];

			                 jem[l][m] = jem[l][m] + j[k][l]*e[k][m];
              
	
						}
           				}

        			}

					//cout<<elm[1][1]<<endl;
	
	for( int l=0;l<DIMENSION;l++){

   	     for (int m=0;m<DIMENSION;m++)	{	

             					 jlm[l][m]=0.0;

        	      for( int k=0;k<DIMENSION;k++	){

                	 			jlm[l][m] = jlm[l][m] + e[k][l]*jem[k][m]; 

 							}

                 				}
        		   		}


        for( int l=0;l<DIMENSION;l++){

        	for(int  m=0;m<DIMENSION;m++){

              	llm[l][m] = jlm[l][m] + jlm[m][l] + 2.0*beta*elm[l][m];

           					}

              llm[l][l] = jlm[l][l] + beta*(elm[l][l]-1.0);



        				}

		//cout<<llm[1][1]<<endl;


 	int i=DIMENSION;

        for( int m=0;m<DIMENSION;m++){

		for(int  k=0;k<DIMENSION;k++){

              	tmp=0.0;

              		for( int l=0;l<=m;l++){

                 	tmp = tmp + e[k][l]*llm[m][l];    

        //     cout<<e[0][0]<<endl;		

					    }

           				   dy[i] = jem[k][m] - tmp;

			              i=i+1;
           					}

        				}

	//			cout<<tmp<<endl;



			dy[6] = jlm[0][0];/*First Lyapunov exponent equation to be integrated*/

       			dy[7] = jlm[1][1];/*Second Lyapunov exponent equation to be integrated*/

	/*cout<<dy[7]<<endl;*/

	}
	
	//////////////RK4 Algorithm////////////////

	void Pendulum::Runge_Kutta4(double t, double Time_Step_Size,double *yin,double *yout)
	{	
	
	double k1[ARRAY_SIZE],k2[ARRAY_SIZE],k3[ARRAY_SIZE],k4[ARRAY_SIZE],y_k[ARRAY_SIZE];

	Derivatives(t,yin,yout);

	for (int s=0;s<ARRAY_SIZE;s++){

	k1[s]=yout[s]*Time_Step_Size;

	y_k[s]=yin[s]+k1[s]*0.5;
					}	

	Derivatives(t+Time_Step_Size*0.5,y_k,yout);

	for (int s=0;s<ARRAY_SIZE;s++){

	k2[s]=yout[s]*Time_Step_Size;
	
	y_k[s]=yin[s]+k2[s]*0.5;

					} 

	
	Derivatives(t+Time_Step_Size*0.5,y_k,yout);

	for(int s=0;s<ARRAY_SIZE;s++){

	k3[s]=yout[s]*Time_Step_Size;
	
	y_k[s]=yin[s]+k3[s];
			
			     		}

	Derivatives(t+Time_Step_Size,y_k,yout);

	for(int s=0;s<ARRAY_SIZE;s++){

	k4[s]=yout[s]*Time_Step_Size;
	
	yout[s]=yin[s]+(1.0/6.0)*(k1[s]+2*k2[s]+2*k3[s]+k4[s]);

	//cout<<s<<"\n";
					 }
	

	}

	

	void Pendulum::test(void)

	{

	//double yy[8], dy[nm];

	//cout<<nm<<"\n";

	t  = 0.0;/*Initial Time*/
	
	Time_Step_Size=0.001;

//	cout<<Time_Step_Size<<"\n";

	int Sampling_Interval=1;

	Real_Time=Time_Step_Size*Sampling_Interval;


	//////////////////////Initialization of the variables//////////////////////////

	double y_Initial[8]={pi,0.10,1.0/sqrt(2.0),1.0/sqrt(2.0),1.0/sqrt(2.0),-1.0/sqrt(2.0),0.0,0.0};


	ofstream Write_To_File("Lyapunov-threshold.dat");

	Write_To_File.setf(ios::scientific);

	Write_To_File.precision(4);

	double y_Final[ARRAY_SIZE];



	for(int iter=1; iter<ITERATION;iter++){

		for(int ii=0;ii<Sampling_Interval;ii++){

	Runge_Kutta4( t,Time_Step_Size,y_Initial,y_Final); /*Calling Runge_Kutta Algorithm*/              
		
	t+=Time_Step_Size;
	
	for( int jj=0;jj<ARRAY_SIZE;jj++){

				y_Initial[jj]=y_Final[jj];	  

        	   			  }

							}

	Write_To_File<<iter*Time_Step_Size<<"\t\t"<<y_Final[0]<<"\t\t"<<y_Final[1]<<"\t\t"<<y_Final[6]/(Time_Step_Size*iter*Sampling_Interval)

	<<"\t\t"<<y_Final[7]/(Time_Step_Size*iter*Sampling_Interval)<<endl;	
			
	  						}

	Write_To_File.close();
	
	}

	//////////////////////////////////////Main Function////////////////////////////////////////

		int main()
	
		{

		Pendulum Kapitza ;

		Kapitza.test();

		return 0;
        
		}



      
	
	



