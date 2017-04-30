#include <stdio.h>
#define _USE_MATH_DEFINES
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <string>

#define MAXLINE 512
#define NA_MPMC 6.0221415e23
#define KB_MPMC 1.3806503e-23
#define ATM2PASCALS 101325.0
#define ATM2PSI 14.6959488

void output(string msg) {

    printf("%s", msg.c_str());
    fflush(stdout);

}


#define BACK_H2_ALPHA	1.033
#define BACK_H2_U0	38.488
#define BACK_H2_V00	9.746
#define BACK_H2_N	0.00
#define BACK_C		0.12

#define BACK_MAX_M	9
#define BACK_MAX_N	4

double h2_comp_back(double temperature, double pressure) {

	double alpha, y;				/* repulsive part of the compressibility factor */
	double V, V0, u, D[BACK_MAX_M][BACK_MAX_N];	/* attractive part */
	int n, m;					/* indices for double sum of attractive part */
	double comp_factor;
	double comp_factor_repulsive;
	double comp_factor_attractive;
	// double fugacity;  (unused variable)


	/* setup the BACK universal D constants */
	D[0][0] = -8.8043;	D[0][1] = 2.9396;	D[0][2] = -2.8225;	D[0][3] = 0.34;
	D[1][0] = 4.164627;	D[1][1] = -6.0865383;	D[1][2] = 4.7600148;	D[1][3] = -3.1875014;
	D[2][0] = -48.203555;	D[2][1] = 40.137956;	D[2][2] = 11.257177;	D[2][3] = 12.231796;
	D[3][0] = 140.4362;	D[3][1] = -76.230797;	D[3][2] = -66.382743;	D[3][3] = -12.110681;
	D[4][0] = -195.23339;	D[4][1] = -133.70055;	D[4][2] = 69.248785;	D[4][3] = 0.0;
	D[5][0] = 113.515;	D[5][1] = 860.25349;	D[5][2] = 0.0;		D[5][3] = 0.0;
	D[6][0] = 0.0;		D[6][1] = -1535.3224;	D[6][2] = 0.0;		D[6][3] = 0.0;
	D[7][0] = 0.0;		D[7][1] = 1221.4261;	D[7][2] = 0.0;		D[7][3] = 0.0;
	D[8][0] = 0.0;		D[8][1] = -409.10539;	D[8][2] = 0.0;		D[8][3] = 0.0;

	/* calculate attractive part */
	V0 = BACK_H2_V00*(1.0 - BACK_C*exp(-3.0*BACK_H2_U0/temperature));
	V = NA_MPMC*KB_MPMC*temperature/(pressure*ATM2PASCALS*1.0e-6);
	u = BACK_H2_U0*(1.0 + BACK_H2_N/temperature);

	comp_factor_attractive = 0;
	for(n = 0; n < BACK_MAX_N; n++)
		for(m = 0; m < BACK_MAX_M; m++)
			comp_factor_attractive += ((double)(m+1))*D[m][n]*pow(u/temperature, ((double)(n+1)))*pow(V0/V, ((double)(m+1)));

	/* calculate repulsive part */
	alpha = BACK_H2_ALPHA;
	y = (M_PI*sqrt(2.0)/6.0)*(pressure*ATM2PASCALS*1.0e-6)/(NA_MPMC*KB_MPMC*temperature)*V0;
	comp_factor_repulsive = 1.0 + (3.0*alpha - 2.0)*y;
	comp_factor_repulsive += (3.0*pow(alpha, 2) - 3.0*alpha + 1.0)*pow(y, 2);
	comp_factor_repulsive -= pow(alpha, 2)*pow(y, 3);
	comp_factor_repulsive /= pow((1.0 - y), 3);

	comp_factor = comp_factor_repulsive + comp_factor_attractive;
	return(comp_factor);

}



/* use the semi-empirical BACK equation of state */
/* Tomas Boublik, "The BACK equation of state for hydrogen and related compounds", Fluid Phase Equilibria, 240, 96-100 (2005) */

double h2_fugacity_back(double temperature, double pressure) {

	double fugacity_coefficient, fugacity;
	double comp_factor;
	double P, dP;
	char linebuf[MAXLINE];

	/* integrate (z-1)/P from 0 to P */
	fugacity_coefficient = 0;
	for(P = 0.001, dP = 0.001; P <= pressure; P += dP) {

		comp_factor = h2_comp_back(temperature, P);
		fugacity_coefficient += dP*(comp_factor - 1.0)/P;

	}
	fugacity_coefficient = exp(fugacity_coefficient);

	comp_factor = h2_comp_back(temperature, pressure);
	sprintf(linebuf, "INPUT: BACK compressibility factor at %.3f atm is %.3f\n", pressure, comp_factor);
	output(linebuf);

	fugacity = pressure*fugacity_coefficient;
	return(fugacity);

}



/* calculate the fugacity correction for H2 for 0 C and higher */
/* this empirical relation follows from: */
/* H.R. Shaw, D.F. Wones, American Journal of Science, 262, 918-929 (1964) */
double h2_fugacity_shaw(double temperature, double pressure) {

	double C1, C2, C3;
	double fugacity, fugacity_coefficient;

	C1 = -3.8402*pow(temperature, 1.0/8.0) + 0.5410;
	C1 = exp(C1);

	C2 = -0.1263*pow(temperature, 1.0/2.0) - 15.980;
	C2 = exp(C2);

	C3 = -0.11901*temperature - 5.941;
	C3 = exp(C3);
	C3 *= 300.0;

	fugacity_coefficient = C1*pressure - C2*pow(pressure, 2) + C3*exp(-pressure/300.0 - 1.0);
	fugacity_coefficient = exp(fugacity_coefficient);
	fugacity = fugacity_coefficient*pressure;

	return(fugacity);

}

/* fugacity for low temperature and up to 200 atm */
/* Zhou, Zhou, Int. J. Hydrogen Energy, 26, 597-601 (2001) */
double h2_fugacity_zhou(double temperature, double pressure) {

	double fugacity, fugacity_coefficient;

	pressure *= ATM2PSI;

	fugacity_coefficient = -1.38130e-4*pressure;
	fugacity_coefficient += 4.67096e-8*pow(pressure, 2)/2;
	fugacity_coefficient += 5.93690e-12*pow(pressure, 3)/3;
	fugacity_coefficient += -3.24527e-15*pow(pressure, 4)/4;
	fugacity_coefficient += 3.54211e-19*pow(pressure, 5)/5;

	pressure /= ATM2PSI;

	fugacity_coefficient = exp(fugacity_coefficient);
	fugacity = pressure*fugacity_coefficient;

	return(fugacity);

}



double h2_fugacity(double temperature, double pressure) {

	if((temperature == 77.0) && (pressure <= 200.0)) {

		output("INPUT: fugacity calculation using Zhou function\n");
		return(h2_fugacity_zhou(temperature, pressure));

	}  else if(temperature >= 273.15) {

		output("INPUT: fugacity calculation using Shaw function\n");
		return(h2_fugacity_shaw(temperature, pressure));

	} else {

		output("INPUT: fugacity calculation using BACK EoS\n");
		return(h2_fugacity_back(temperature, pressure));

	}

	//return(0); /* NOT REACHED */

}


#define MWCH4 16.043
#define BACK_CH4_ALPHA  1.000
#define BACK_CH4_U0     188.047
#define BACK_CH4_V00    21.532
#define BACK_CH4_N      2.40
#define BACK_C          0.12

#define BACK_MAX_M      9
#define BACK_MAX_N      4

double ch4_comp_back(double temperature, double pressure) {

	double alpha, y;                                /* repulsive part of the compressibility factor */
	double V, V0, u, D[BACK_MAX_M][BACK_MAX_N];     /* attractive part */
	int n, m;                                       /* indices for double sum of attractive part */
	double comp_factor;
	double comp_factor_repulsive;
	double comp_factor_attractive;
	// double fugacity;  (unused variable)
	
	/* setup the BACK universal D constants */
	D[0][0] = -8.8043;      D[0][1] = 2.9396;       D[0][2] = -2.8225;      D[0][3] = 0.34;
	D[1][0] = 4.164627;     D[1][1] = -6.0865383;   D[1][2] = 4.7600148;    D[1][3] = -3.1875014;
	D[2][0] = -48.203555;   D[2][1] = 40.137956;    D[2][2] = 11.257177;    D[2][3] = 12.231796;
	D[3][0] = 140.4362;     D[3][1] = -76.230797;   D[3][2] = -66.382743;   D[3][3] = -12.110681;
	D[4][0] = -195.23339;   D[4][1] = -133.70055;   D[4][2] = 69.248785;    D[4][3] = 0.0;
	D[5][0] = 113.515;      D[5][1] = 860.25349;    D[5][2] = 0.0;          D[5][3] = 0.0;
	D[6][0] = 0.0;          D[6][1] = -1535.3224;   D[6][2] = 0.0;          D[6][3] = 0.0;
	D[7][0] = 0.0;          D[7][1] = 1221.4261;    D[7][2] = 0.0;          D[7][3] = 0.0;
	D[8][0] = 0.0;          D[8][1] = -409.10539;   D[8][2] = 0.0;          D[8][3] = 0.0;
	
	/* calculate attractive part */
	V0 = BACK_CH4_V00*(1.0 - BACK_C*exp(-3.0*BACK_CH4_U0/temperature));
	V = NA_MPMC*KB_MPMC*temperature/(pressure*ATM2PASCALS*1.0e-6);
	u = BACK_CH4_U0*(1.0 + BACK_CH4_N/temperature);
	
	comp_factor_attractive = 0;
	for(n = 0; n < BACK_MAX_N; n++)
		for(m = 0; m < BACK_MAX_M; m++)
			comp_factor_attractive += ((double)(m+1))*D[m][n]*pow(u/temperature, ((double)(n+1)))*pow(V0/V, ((double)(m+1)));
	
	/* calculate repulsive part */
	alpha = BACK_CH4_ALPHA;
	y = (M_PI*sqrt(2.0)/6.0)*(pressure*ATM2PASCALS*1.0e-6)/(NA_MPMC*KB_MPMC*temperature)*V0;
	comp_factor_repulsive = 1.0 + (3.0*alpha - 2.0)*y;
	comp_factor_repulsive += (3.0*pow(alpha, 2) - 3.0*alpha + 1.0)*pow(y, 2);
	comp_factor_repulsive -= pow(alpha, 2)*pow(y, 3);
	comp_factor_repulsive /= pow((1.0 - y), 3);
	
	comp_factor = comp_factor_repulsive + comp_factor_attractive;
	return(comp_factor);

}




/* Incorporate BACK EOS */
double ch4_fugacity_back(double temperature, double pressure) {

        double fugacity_coefficient, fugacity;
        double comp_factor;
        double P, dP;
        char linebuf[MAXLINE];

        /* integrate (z-1)/P from 0 to P */
        fugacity_coefficient = 0;
        for(P = 0.001, dP = 0.001; P <= pressure; P += dP) {

                comp_factor = ch4_comp_back(temperature, P);
                fugacity_coefficient += dP*(comp_factor - 1.0)/P;

        }
        fugacity_coefficient = exp(fugacity_coefficient);

        comp_factor = ch4_comp_back(temperature, pressure);
        sprintf(linebuf, "INPUT: CH4 BACK compressibility factor at %.3f atm is %.3f\n", pressure, comp_factor);
        output(linebuf);

        fugacity = pressure*fugacity_coefficient;
        return(fugacity);
}
/* Apply the Peng-Robinson EoS to methane */
double ch4_fugacity_PR(double temperature, double pressure) {

	double Z, A, B, aa, bb, TcCH4, PcCH4, Tr;
	double alpha, alpha2, wCH4, R, Q, X, j, k, l;
	double theta, Q3;
	double uu, U, V, root1, root2, root3, stuff1, stuff2, stuff3; //, answer; (unused variable)
	double f1, f2, f3, f4, fugacity, lnfoverp;
	double pi=acos(-1.0);

	/*Peng Robinson variables and equations for CH4 units K,atm, L, mole*/
	TcCH4 = 190.564;    /* K */
	PcCH4 = 45.391;    /* atm */
	wCH4  = 0.01142;
	R     = 0.08206;   /* gas constant atmL/moleK */
	
	aa = (0.45724*R*R*TcCH4*TcCH4) / PcCH4;
	bb = (0.07780*R*TcCH4) / PcCH4;
	Tr = temperature / TcCH4;
	stuff1 = 0.37464 + 1.54226*wCH4 - 0.26992*wCH4*wCH4;
	stuff2=1.0-sqrt(Tr);
	alpha=1.0+stuff1*stuff2;
	alpha2=alpha*alpha;
	A=alpha2*aa*pressure/(R*R*temperature*temperature);
	B=bb*pressure/(R*temperature);
	
	/* solving a cubic equation part */
	j=-1.0*(1-B);
	k=A-3.0*B*B-2.0*B;
	l= -1*(A*B- B*B -B*B*B);
	Q=(j*j-3.0*k)/9.0;
	X=(2.0*j*j*j -9.0*j*k+27.0*l)/54.0;
	Q3=Q*Q*Q;

	/* Need to check X^2 < Q^3 */
	if((X*X)<(Q*Q*Q)){    /* THREE REAL ROOTS  */
		theta=acos((X/sqrt(Q3)));
		root1=-2.0*sqrt(Q)*cos(theta/3.0)-j/3.0;
		root2=-2.0*sqrt(Q)*cos((theta+2.0*pi)/3.0)-j/3.0;
		root3=-2.0*sqrt(Q)*cos((theta-2.0*pi)/3.0)-j/3.0;

		/*Choose the root closest to 1, which is "ideal gas law" */
		if((1.0-root1)<(1.0-root2) && (1.0-root1)<(1.0-root3))
			Z=root1;
		else if((1.0-root2)<(1.0-root3) && (1.0-root2)<(1.0-root1))
			Z=root2;
		else
			Z=root3;
	}
	else {   /* ONLY ONE real root */
		stuff3= X*X-Q*Q*Q;
		uu=X-sqrt(stuff3);
		/*Power function must have uu a positive number*/
		if(uu<0.0)
			uu=-1.0*uu;
		U=pow(uu,(1.0/3.0));
		V=Q/U;
		root1=U+V-j/3.0;
		Z=root1;
	}

	/* using Z calculate the fugacity */
	f1=(Z-1.0)-log(Z-B);
	f2=A/(2.0*sqrt(2.0)*B);
	f3=Z+(1.0+sqrt(2.0))*B;
	f4=Z+(1.0-sqrt(2.0))*B;
	lnfoverp=f1-f2*log(f3/f4);
	fugacity=exp(lnfoverp)*pressure;

	return(fugacity);
}


/* ***************************** CH4 EQUATION OF STATE *************************************** */
double ch4_fugacity(double temperature, double pressure) {

        if((temperature >= 298.0) && (temperature <= 300.0) && (pressure <= 500.0)) {

                output("INPUT: CH4 fugacity calculation using BACK EoS\n");
                return(ch4_fugacity_back(temperature, pressure));

        } else if((temperature == 150.0) && (pressure <= 200.0)) {

                output("INPUT: CH4 fugacity calculation using Peng-Robinson EoS\n");
                return(ch4_fugacity_PR(temperature, pressure));

        } else {

                output("INPUT: Unknown if CH4 fugacity will be correct at the requested temperature & pressure...defaulting to use the BACK EoS.\n");
                return(ch4_fugacity_back(temperature, pressure));

        }

        //return(0); /* NOT REACHED */

}



/* ******************************* END CH4 FUGACITY ****************************************** */





/* *************************** N2 BACK EQUATION OF STATE ************************************* */



#define BACK_N2_ALPHA   1.048
#define BACK_N2_U0      120.489
#define BACK_N2_V00     18.955
#define BACK_N2_N       10.81
#define BACK_C          0.12

#define BACK_MAX_M      9
#define BACK_MAX_N      4

double n2_comp_back(double temperature, double pressure) {

	double alpha, y;                                /* repulsive part of the compressibility factor */
	double V, V0, u, D[BACK_MAX_M][BACK_MAX_N];     /* attractive part */
	int n, m;                                       /* indices for double sum of attractive part */
	double comp_factor;
	double comp_factor_repulsive;
	double comp_factor_attractive;
	// double fugacity;  (unused variable)
	
	
	/* setup the BACK universal D constants */
	D[0][0] = -8.8043;      D[0][1] = 2.9396;       D[0][2] = -2.8225;      D[0][3] = 0.34;
	D[1][0] = 4.164627;     D[1][1] = -6.0865383;   D[1][2] = 4.7600148;    D[1][3] = -3.1875014;
	D[2][0] = -48.203555;   D[2][1] = 40.137956;    D[2][2] = 11.257177;    D[2][3] = 12.231796;
	D[3][0] = 140.4362;     D[3][1] = -76.230797;   D[3][2] = -66.382743;   D[3][3] = -12.110681;
	D[4][0] = -195.23339;   D[4][1] = -133.70055;   D[4][2] = 69.248785;    D[4][3] = 0.0;
	D[5][0] = 113.515;      D[5][1] = 860.25349;    D[5][2] = 0.0;          D[5][3] = 0.0;
	D[6][0] = 0.0;          D[6][1] = -1535.3224;   D[6][2] = 0.0;          D[6][3] = 0.0;
	D[7][0] = 0.0;          D[7][1] = 1221.4261;    D[7][2] = 0.0;          D[7][3] = 0.0;
	D[8][0] = 0.0;          D[8][1] = -409.10539;   D[8][2] = 0.0;          D[8][3] = 0.0;
	
	/* calculate attractive part */
	V0 = BACK_N2_V00*(1.0 - BACK_C*exp(-3.0*BACK_N2_U0/temperature));
	V = NA_MPMC*KB_MPMC*temperature/(pressure*ATM2PASCALS*1.0e-6);
	u = BACK_N2_U0*(1.0 + BACK_N2_N/temperature);

	comp_factor_attractive = 0;
	for(n = 0; n < BACK_MAX_N; n++)
		for(m = 0; m < BACK_MAX_M; m++)
			comp_factor_attractive += ((double)(m+1))*D[m][n]*pow(u/temperature, ((double)(n+1)))*pow(V0/V, ((double)(m+1)));
	
	/* calculate repulsive part */
	alpha = BACK_N2_ALPHA;
	y = (M_PI*sqrt(2.0)/6.0)*(pressure*ATM2PASCALS*1.0e-6)/(NA_MPMC*KB_MPMC*temperature)*V0;
	comp_factor_repulsive = 1.0 + (3.0*alpha - 2.0)*y;
	comp_factor_repulsive += (3.0*pow(alpha, 2) - 3.0*alpha + 1.0)*pow(y, 2);
	comp_factor_repulsive -= pow(alpha, 2)*pow(y, 3);
	comp_factor_repulsive /= pow((1.0 - y), 3);
	
	comp_factor = comp_factor_repulsive + comp_factor_attractive;
	return(comp_factor);
}




/* Incorporate BACK EOS */
double n2_fugacity_back(double temperature, double pressure) {

	double fugacity_coefficient, fugacity;
	double comp_factor;
	double P, dP;
	char linebuf[MAXLINE];
	
	/* integrate (z-1)/P from 0 to P */
	fugacity_coefficient = 0;
	for(P = 0.001, dP = 0.001; P <= pressure; P += dP) {
	
		comp_factor = n2_comp_back(temperature, P);
		fugacity_coefficient += dP*(comp_factor - 1.0)/P;
	
	}
	fugacity_coefficient = exp(fugacity_coefficient);
	
	comp_factor = n2_comp_back(temperature, pressure);
	sprintf(linebuf, "INPUT: BACK compressibility factor at %.3f atm is %.3f\n", pressure, comp_factor);
	output(linebuf);
	
	fugacity = pressure*fugacity_coefficient;
	return(fugacity);

}
/* Apply the Peng-Robinson EoS to N2 */
double n2_fugacity_PR(double temperature, double pressure) {

	double Z, A, B, aa, bb, TcN2, PcN2, Tr;
	double alpha, alpha2, wN2, R, Q, X, j, k, l;
	double theta, Q3;
	double uu, U, V, root1, root2, root3, stuff1, stuff2, stuff3; //, answer;  (unused variable)
	double f1, f2, f3, f4, fugacity, lnfoverp;
	double pi=acos(-1.0);
	
	/*Peng Robinson variables and equations for N2 units K,atm, L, mole*/
	TcN2 = 126.192;   /* K */
	PcN2 = 33.514;    /* atm */
	wN2  = 0.037;
	R    = 0.08206;   /* gas constant atmL/moleK */
	
	aa = (0.45724*R*R*TcN2*TcN2) / PcN2;
	bb = (0.07780*R*TcN2) / PcN2;
	Tr = temperature / TcN2;
	stuff1 = 0.37464 + 1.54226*wN2 - 0.26992*wN2*wN2;
	stuff2=1.0-sqrt(Tr);
	alpha=1.0+stuff1*stuff2;
	alpha2=alpha*alpha;
	A=alpha2*aa*pressure/(R*R*temperature*temperature);
	B=bb*pressure/(R*temperature);
	
	/* solving a cubic equation part */
	j=-1.0*(1-B);
	k=A-3.0*B*B-2.0*B;
	l= -1*(A*B- B*B -B*B*B);
	Q=(j*j-3.0*k)/9.0;
	X=(2.0*j*j*j -9.0*j*k+27.0*l)/54.0;
	Q3=Q*Q*Q;
	
	/* Need to check X^2 < Q^3 */
	if((X*X)<(Q*Q*Q)){    /* THREE REAL ROOTS  */
		theta=acos((X/sqrt(Q3)));
		root1=-2.0*sqrt(Q)*cos(theta/3.0)-j/3.0;
		root2=-2.0*sqrt(Q)*cos((theta+2.0*pi)/3.0)-j/3.0;
		root3=-2.0*sqrt(Q)*cos((theta-2.0*pi)/3.0)-j/3.0;

		/*Choose the root closest to 1, which is "ideal gas law" */
		if((1.0-root1)<(1.0-root2) && (1.0-root1)<(1.0-root3))
			Z=root1;
		else if((1.0-root2)<(1.0-root3) && (1.0-root2)<(1.0-root1))
			Z=root2;
		else
			Z=root3;
			
	} else {   /* ONLY ONE real root */
		
		stuff3= X*X-Q*Q*Q;
		uu=X-sqrt(stuff3);
		/*Power function must have uu a positive number*/
		if(uu<0.0)
			uu=-1.0*uu;
		U=pow(uu,(1.0/3.0));
		V=Q/U;
		root1=U+V-j/3.0;
		Z=root1;
	}

	/* using Z calculate the fugacity */
	f1=(Z-1.0)-log(Z-B);
	f2=A/(2.0*sqrt(2.0)*B);
	f3=Z+(1.0+sqrt(2.0))*B;
	f4=Z+(1.0-sqrt(2.0))*B;
	lnfoverp=f1-f2*log(f3/f4);
	fugacity=exp(lnfoverp)*pressure;

	return(fugacity);
}

/* Apply the Zhou function to N2 */
double n2_fugacity_zhou(double temperature, double pressure) {

	double fugacity_coefficient, fugacity;
	
	output("INPUT: N2 fugacity calculation using Zhou function\n");
	
	pressure *= ATM2PSI;
	
	fugacity_coefficient = -1.38130e-4*pressure;
	fugacity_coefficient += 4.67096e-8*pow(pressure, 2)/2;
	fugacity_coefficient += 5.93690e-12*pow(pressure, 3)/3;
	fugacity_coefficient += -3.24527e-15*pow(pressure, 4)/4;
	fugacity_coefficient += 3.54211e-19*pow(pressure, 5)/5;
	
	pressure /= ATM2PSI;
	
	fugacity_coefficient = exp(fugacity_coefficient);
	fugacity = pressure*fugacity_coefficient;
	
	return(fugacity);
}


double n2_fugacity(double temperature, double pressure) {

	if((temperature == 78.0) && (pressure <= 1.0)) {
	
		output("INPUT: N2 fugacity calculation using Zhou\n");
		return(n2_fugacity_zhou(temperature, pressure));
	
	} else if((temperature == 78.0) && (pressure >= 10.0) && (pressure <= 300.0)) {
	
		output("INPUT: N2 fugacity calculation using Peng-Robinson EoS\n");
		return(n2_fugacity_PR(temperature, pressure));
	
	} else if((temperature == 150.0) && (pressure < 175.0)) {
	
		output("INPUT: N2 fugacity calculation using Peng-Robinson EoS\n");
		return(n2_fugacity_PR(temperature, pressure));
	
	} else if((temperature == 150.0) && (pressure >= 175.0) && (pressure <= 325.0)) {
	
		output("INPUT: N2 fugacity calculation using BACK EoS\n");
		return(n2_fugacity_back(temperature, pressure));
	
	} else if((temperature >= 298.0) && (temperature <= 300.0) && (pressure <= 350.0)) {
	
		output("INPUT: N2 fugacity calculation using Peng-Robinson EoS\n");
		return(n2_fugacity_PR(temperature, pressure));
	
	} else {
		output("INPUT: Unknown if N2 fugacity will be correct at the requested temperature & pressure...defaulting to use the PR EoS.\n");
		return(n2_fugacity_PR(temperature, pressure));
	
	}
	
	//return(0); /* NOT REACHED */

}



/* ********************************** END N2 ****************************************************** */


/* reads in temperature in K, and pressure (of the ideal gas in the resevoir) in atm */
/* return the CO2 fugacity via the Peng-Robinson equation of state */
/* else return 0.0 on error - I don't have an error statement*/
/* units are  atm, K,  */

double co2_fugacity(double temperature, double pressure) {

	double Z, A, B, aa, bb, Tc, Pc, Tr;
	double alpha,alpha2, w, R, Q, X, j,k,l;
	double theta,  Q3;
	double uu,U,V,root1,root2,root3, stuff1, stuff2, stuff3; //, answer;  (unused variable)
	double f1, f2, f3, f4, fugacity, lnfoverp;
	double pi=acos(-1.0);

	/*Peng Robinson variables and equations for CO2 units K,atm, L, mole*/
	Tc=304.12;   /*K  */
	Pc=73.74/1.01325;    /*bar to atm*/
	w=0.225;
	R=0.08206;   /* gas constant atmL/moleK */
	
	aa=0.45724*R*R*Tc*Tc/Pc;
	bb=0.07780*R*Tc/Pc;
	Tr=temperature/Tc;
	stuff1=0.37464+1.54226*w -0.26992*w*w;
	stuff2=1.0-sqrt(Tr);
	alpha=1.0+stuff1*stuff2;
	alpha2=alpha*alpha;
	A=alpha2*aa*pressure/(R*R*temperature*temperature);
	B=bb*pressure/(R*temperature);

	/* solving a cubic equation part */
	j=-1.0*(1-B);
	k=A-3.0*B*B-2.0*B;
	l= -1*(A*B- B*B -B*B*B);
	Q=(j*j-3.0*k)/9.0;
	X=(2.0*j*j*j -9.0*j*k+27.0*l)/54.0;
	Q3=Q*Q*Q;

	/* Need to check X^2 < Q^3 */
	if((X*X)<(Q*Q*Q)){    /* THREE REAL ROOTS  */
		theta=acos((X/sqrt(Q3)));
		root1=-2.0*sqrt(Q)*cos(theta/3.0)-j/3.0;
		root2=-2.0*sqrt(Q)*cos((theta+2.0*pi)/3.0)-j/3.0;
		root3=-2.0*sqrt(Q)*cos((theta-2.0*pi)/3.0)-j/3.0;

		/*Choose the root closest to 1, which is "ideal gas law" */
		if((1.0-root1)<(1.0-root2) && (1.0-root1)<(1.0-root3))
			Z=root1;
		else if((1.0-root2)<(1.0-root3) && (1.0-root2)<(1.0-root1))
			Z=root2;
		else
			Z=root3;
			
	} else {   /* ONLY ONE real root */
		
		stuff3= X*X-Q*Q*Q;
		uu=X-sqrt(stuff3);
		/*Power function must have uu a positive number*/
		if(uu<0.0)
			uu=-1.0*uu;
		U=pow(uu,(1.0/3.0));
		V=Q/U;
		root1=U+V-j/3.0;
		Z=root1;
	}

	/* using Z calculate the fugacity */
	f1=(Z-1.0)-log(Z-B);
	f2=A/(2.0*sqrt(2.0)*B);
	f3=Z+(1.0+sqrt(2.0))*B;
	f4=Z+(1.0-sqrt(2.0))*B;
	lnfoverp=f1-f2*log(f3/f4);
	fugacity=exp(lnfoverp)*pressure;
	
	return(fugacity);

}


