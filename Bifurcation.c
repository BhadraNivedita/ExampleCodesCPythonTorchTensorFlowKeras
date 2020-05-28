
/*------------------------------------------------------------Logistic Map-------------------------------------------------------------------*/

	#include <stdio.h>

	#include <math.h>

	#define ITERATION 50

	#define CONDITION  100

	#define randomize() srand ((unsigned) time (NULL))


/*---------------------------------------CREATING RANDOM NUMBER FOR RANDOM VALUE OF R-PARAMETER----------------------------------------------*/

	FILE *output;  

	int random_condition(void)

	{
	     int ran;

	     ran =rand()%5;

	     return ran;
		
	} 

/*---------------------------------------------------------------------------------------------------------------------------------------------*/

	main()

	{

		int n,p;

		double r_PARA,x_INITIAL,x_FINAL,pi;

		pi=M_PI;


		output=fopen("logisticmap.dat", "w"); /*Writing to file*/

		for(p=0;p<CONDITION;p++){     /*Loop for different value of r*/

		x_INITIAL=0.10;

		//r_PARA=.15*random_condition();

		r_PARA=-1.0+(2.0/CONDITION)*p;

		 // printf("the initial value is %f\n",x0);

		//printf("the r value is %f\n",r_PARA);


					for(n=0;n<ITERATION;n++){    /*Loop for discrete sum*/

				        x_FINAL=x_INITIAL+r_PARA*x_INITIAL*(1.0-x_INITIAL); /*Updating value of x_INI*/

			                printf("the value of xn is %f\n",x_FINAL);

								//if(n>50){


									//fprintf(output,"%d\t%f\t%f\t%f\n",n,r_PARA,x_INITIAL,x_FINAL);

									//}

									x_INITIAL=x_FINAL;

						                }

									fprintf(output,"%d\t%f\t%f\t%f\n",n,r_PARA,x_INITIAL,x_FINAL);

									//fprintf(output,"\n");/*  Space to create blocks for several value of r*/
    
					}


				fclose(output);
	}
/*----------------------------------------------------------------------------------------------------------------------------------------------*/
