

// Header file for input output functions
#include<iostream>
#include<stdlib.h>
#include<math.h>
#include<cmath>
#include<time.h>
#include <random>
#include <vector>
#include "nr3.h"
#include "mins.h"
#include "mins_ndim.h"
#include <fstream>
#include "QRdcmp.h"
#include "ludcmp.h"
#include "roots_multidim.h"
#include "sort.h"
#include <unistd.h>  // for sleep function

using namespace std;


double normalCDF(double x) // Phi(-∞, x) aka N(x)
{
    return std::erfc(-x/std::sqrt(2))/2;
}

double RationalApproximation(double t)
{
    // Abramowitz and Stegun formula 26.2.23.
    // The absolute value of the error should be less than 4.5 e-4.
    double c[] = {2.515517, 0.802853, 0.010328};
    double d[] = {1.432788, 0.189269, 0.001308};
    return t - ((c[2]*t + c[1])*t + c[0]) /
               (((d[2]*t + d[1])*t + d[0])*t + 1.0);
}

double NormalCDFInverse(double p)
{
    if (p <= 0.0 || p >= 1.0)
    {
        std::stringstream os;
        os << "Invalid input argument (" << p
           << "); must be larger than 0 but less than 1.";
        throw std::invalid_argument( os.str() );
    }

    // See article above for explanation of this section.
    if (p < 0.5)
    {
        // F^-1(p) = - G^-1(p)
        return -RationalApproximation( sqrt(-2.0*log(p)) );
    }
    else
    {
        // F^-1(p) = G^-1(1-p)
        return RationalApproximation( sqrt(-2.0*log(1-p)) );
    }
}


MatDoub generatedata(double coxlambda, double coxbeta1, double coxbeta2, double coxbeta3,
	    int workingn, double propdied) {

    // Rows are people,
    // Cols are bmkr value, trtrec, survtime, survstatus
    MatDoub myreturn(workingn,4);

	std::default_random_engine generator;

    for (int k=0;k<workingn;k++){
    	std::normal_distribution<double> distribution(0.5,sqrt(1.0/12.0));
    	double onenormMn0VarOneTwelfth = distribution(generator);
       	myreturn[k][0] = onenormMn0VarOneTwelfth;

    	double oneBernoulli  = rand() % 2;
    	myreturn[k][1] = oneBernoulli;

    	double thisrate = coxlambda * exp(coxbeta1*myreturn[k][0] + coxbeta2*myreturn[k][1] + coxbeta3*myreturn[k][0]*myreturn[k][1]);
    	std::exponential_distribution<double> distribution2(thisrate);
    	double oneExponential = distribution2(generator);
    	myreturn[k][2] = oneExponential;

    	std::bernoulli_distribution distribution3(propdied);
    	bool oneBernoulli4died = distribution3(generator);
    	myreturn[k][3] = oneBernoulli4died;
    }
    return myreturn;

}

double normalqfun(double p,double nmean,double nsd) {
	double myreturn = 0;
	double mycritval = NormalCDFInverse(p);
	myreturn = nmean + mycritval * nsd;
	return myreturn;
}

double coxbeta1fun(double k1, double k2) {
	double myreturn = 0;
	if (k1 == k2) {
		myreturn = 0;
	}
	if (k1 != k2) {
		double mynum = log(log(k1)/log(k2));
		double myden = normalqfun(0.25,0.5,sqrt(1.0/12.0)) - normalqfun(0.75,0.5,sqrt(1.0/12.0));
		myreturn = mynum/myden;
	}
	return myreturn;
}

double coxlambdafun(double k1, double k2, double t0) {
	double myreturn = 0;
	if (k1 == k2) {
		myreturn = -log(k1)/t0;
	}
	else {
		double mynum = log(k2/k1);
		double myden = t0*(exp(coxbeta1fun(k1, k2)*normalqfun(0.25,0.5,sqrt(1.0/12.0))) - exp(coxbeta1fun(k1, k2)*normalqfun(0.75,0.5,sqrt(1.0/12.0))));
		myreturn = mynum/myden;
	}
	return myreturn;
}

double coxbeta3fun(double k1, double k2, double k3, double k4, double t0) {
	double myreturn = 0;
	if (k1 == k2 && k3 == k4) {
		myreturn = 0;
	}
	else if (k1==k2 && k3 != k4) {
		myreturn = (k3 - k4)/(normalqfun(0.25,0.5,sqrt(1.0/12.0))-normalqfun(0.75,0.5,sqrt(1.0/12.0)));
	}
	else {
	 double mynum = log(log(k2)/log(k1)) - log(log(k4)/log(k3));
	 double myden = normalqfun(0.25,0.5,sqrt(1.0/12.0)) - normalqfun(0.75,0.5,sqrt(1.0/12.0)) ;
	 myreturn = mynum/myden;
	}
	return myreturn;
}

double coxbeta2fun(double k1, double k2, double k3, double k4, double t0) {
	double myreturn = 0;
	double k1t = log(-log(k1));
	double k2t = log(-log(k2));
	double k3t = log(-log(k4));
	double k4t = log(-log(k4));

  if (k1 == k2 && k3 == k4) {
    myreturn = log(log(k3) / log(k1));
  }
  else if (k1 == k2 && k3 != k4) {
    myreturn = k3t - k1t - normalqfun(0.25,0.5,sqrt(1.0/12.0)) * ((k4t - k3t) / (normalqfun(0.75,0.5,sqrt(1.0/12.0)) - normalqfun(0.25,0.5,sqrt(1.0/12.0))));
  }
  else {
    double prelognum = (-1)*log(k3*k4);
    double coxbetasum = coxbeta1fun(k1, k2) + coxbeta3fun(k1, k2, k3, k4, t0);
    double prelogden = t0*coxlambdafun(k1, k2, t0)*(exp(coxbetasum*normalqfun(0.25,0.5,sqrt(1.0/12.0))) + exp(coxbetasum*normalqfun(0.75,0.5,sqrt(1.0/12.0))));
    myreturn = log(prelognum / prelogden);
  }
  return myreturn;
}


VecDoub getcoxlambdasbetas(double k1, double k2, double k3, double k4, double t0) {
	VecDoub myreturn(4);
	myreturn[0] = coxlambdafun(k1,k2,t0);
	myreturn[1] = coxbeta1fun(k1,k2);
	myreturn[2] = coxbeta2fun(k1,k2,k3,k4,t0);
	myreturn[3] = coxbeta3fun(k1,k2,k3,k4,t0);
	return myreturn;
}


//  HERE IS THE LIKELIHOOD FOR THE COX REGRESSION
struct coxlikelihood{
	MatDoub simdata;
	coxlikelihood ();
	void Setcoxlikelihood( MatDoub &simdata );
	VecDoub operator()(const VecDoub_I x);

};

coxlikelihood::coxlikelihood() { }

void coxlikelihood::Setcoxlikelihood( MatDoub& simdata) {
	simdata = simdata;
}

VecDoub coxlikelihood::operator()(const VecDoub x) {
	int kkdn = simdata.nrows();
	VecDoub myreturn(4);
	int i;
	Doub thislambdai, thisbmkr, thistrt, thissurvtime, thissurvstat;
	for (i=0;i<4;i++) { myreturn[i]  = 0.0; }
	for (i=0;i<kkdn;i++) {
		thisbmkr = simdata[i][0];
		thistrt = simdata[i][1];
		thissurvtime = simdata[i][2];
		thissurvstat = simdata[i][3];
		thislambdai = exp(x[0] + x[1]*thisbmkr + x[2]*thistrt + x[3]*thisbmkr*thistrt);

		myreturn[0] = myreturn[0] - thissurvstat + thissurvtime * thislambdai;
		myreturn[1] = myreturn[1] - thissurvstat*thisbmkr + thissurvtime*thisbmkr*thislambdai;
		myreturn[2] = myreturn[2] - thissurvstat*thistrt + thissurvtime*thistrt*thislambdai;
		myreturn[3] = myreturn[3] - thissurvstat*thistrt*thisbmkr + thissurvtime*thistrt*thisbmkr*thislambdai;
	}
	return myreturn;
}

// HERE IS THE LIKELIHOOD FOR LOGISTIC REGRESSION

struct logisticlikelihood{
	MatDoub simdata;
	Doub t0;
	logisticlikelihood ();
	void Setlogisticlikelihood( MatDoub &simdata, Doub &t0 );
	VecDoub operator()(const VecDoub_I x);
};

logisticlikelihood::logisticlikelihood() { }

void logisticlikelihood::Setlogisticlikelihood( MatDoub& simdata, Doub &t0) {
	simdata = simdata;
	t0 = t0;
}

VecDoub logisticlikelihood::operator()(const VecDoub x) {
	int kkdn = simdata.nrows();
	VecDoub myreturn(4);
	int i;
	Doub thispii, thisli, thisbmkr, thistrt, thissurvtime, thissurvstat, thisBinaryResp;
	for (i=0;i<4;i++) { myreturn[i]  = 0.0; }
	for (i=0;i<kkdn;i++) {
		thisbmkr = simdata[i][0];
		thistrt = simdata[i][1];
		thissurvtime = simdata[i][2];
		thissurvstat = simdata[i][3];

		thisBinaryResp = 99;  // indicator of missing data
		if (thissurvtime < t0 && thissurvstat == 1) { thisBinaryResp = 1; } // bad outcome group
		if (thissurvtime >= t0) { thisBinaryResp = 0; } // good outcome group

		thisli = x[0] + x[1]*thisbmkr + x[2]*thistrt + x[3]*thisbmkr*thistrt;
		thispii = exp(thisli)/(1+exp(thisli));

		if (thisBinaryResp != 99) {
			myreturn[0] = myreturn[0] - thisBinaryResp + thispii;
			myreturn[1] = myreturn[1] - thisBinaryResp * thisbmkr + thisbmkr*thispii;
			myreturn[2] = myreturn[2] - thisBinaryResp * thistrt + thistrt*thispii;
			myreturn[3] = myreturn[3] - thisBinaryResp * thistrt * thisbmkr + thistrt*thisbmkr*thispii;
		}
	}
	return myreturn;
}




// HERE IS THE JANESTHETA STRUCT
struct JanesTheta {
public:
	MatDoub simdata;
	Doub t0;
	Doub beta0hatlogit, beta1hatlogit, beta2hatlogit, beta3hatlogit;
	Int kkdn;
	JanesTheta();
	void setJanesTheta ( MatDoub_I & , Doub & ,
	Doub & inputbeta0, Doub & inputbeta1, Doub & inputbeta2, Doub& inputbeta3);
	Doub Risk0 ( Doub thisbmkrval );
	Doub Risk1 ( Doub thisbmkrval );
	Doub Deltahat ( Doub thisbmkrval );
	Doub ProbDis1givenAis1andDeltaYover0 ( );
	Doub ProbDis1givenAis0andDeltaYover0 ( );
	Doub ProbDis1givenAis1andDeltaYunder0 ( );
	Doub ProbDis1givenAis0andDeltaYunder0 ( );
	Doub hatBsubneg( );
	Doub hatBsubpos( );
	Doub phatsubneg();
	Doub phatsubpos();
	Doub theta0hat ();
	Doub theta1hat ();
	MatDoub PercentileBootstrap (VecDoub_I percentiles, int nboots);
	MatDoub DoTheRefinedBootstrap (int nboots);
};

JanesTheta::JanesTheta() {}

void JanesTheta::setJanesTheta(MatDoub_I &inputdataset, Doub & inputt0,
	Doub & inputbeta0, Doub & inputbeta1, Doub & inputbeta2, Doub& inputbeta3)
{
    MatDoub simdata = inputdataset;
    kkdn = inputdataset.nrows();
	t0 = inputt0;
	beta0hatlogit = inputbeta0;
	beta1hatlogit = inputbeta1;
	beta2hatlogit = inputbeta2;
	beta3hatlogit = inputbeta3;
}

Doub JanesTheta::Risk0( Doub thisbmkrval ) {
	Doub thisreturn;
//	thisreturn = 1.0 - exp( -exp( beta0hatlogit + beta1hatlogit * thisbmkrval)); //Cox regression
	Doub myl1 = beta0hatlogit + beta1hatlogit * thisbmkrval;
	thisreturn = exp( myl1 )/(1+exp( myl1 )); // regression

	return thisreturn;
}

Doub JanesTheta::Risk1( Doub thisbmkrval ) {
	Doub thisreturn;

//	thisreturn = 1.0 - exp( -exp( beta0hatlogit + beta1hatlogit * thisbmkrval + beta2hatlogit + beta3hatlogit * thisbmkrval));
	Doub myl1 = beta0hatlogit + beta1hatlogit * thisbmkrval + beta2hatlogit + beta3hatlogit * thisbmkrval;
	thisreturn = exp( myl1 )/(1+exp( myl1 )); //logistic regression

	return thisreturn;
}

Doub JanesTheta::Deltahat( Doub thisbmkrval ) {
	Doub thisreturn;
	thisreturn = Risk0(thisbmkrval) - Risk1(thisbmkrval);
	return thisreturn;
}

Doub JanesTheta::ProbDis1givenAis1andDeltaYover0( ) {
	Doub myreturn;

	double mynumerator = 0.0, mydenominator = 0.0;
	Doub thisdeltahat, thisbmkrval, thistime, thistrt, thissurvstat, thisnotcensored;
	Doub thisnumsubset, thisdenomsubset;


	for (int i=0;i<kkdn;i++) {
		thisbmkrval = simdata[i][0];
		thistrt = simdata[i][1];
		thistime = simdata[i][2];
		thissurvstat = simdata[i][3];
		thisdeltahat = Deltahat(thisbmkrval);

		Bool condition1 = ((thistime > t0) || (thissurvstat == 1));
		Bool condition2 = (thisdeltahat > 0);
		Bool condition3 = (thistrt == 1);
		if (condition1 && condition2 && condition3) {
			mydenominator = mydenominator + 1.0;
		}
		Bool condition4 = (thistime < t0);
		if (condition1 && condition2 && condition3 && condition4) {
			mynumerator = mynumerator + 1.0;
		}
	}
	myreturn = mynumerator/mydenominator;

	return myreturn;
}

Doub JanesTheta::ProbDis1givenAis0andDeltaYover0( ) {
	Doub myreturn;

	double mynumerator = 0.0, mydenominator = 0.0;
	Doub thisdeltahat, thisbmkrval, thistime, thistrt, thissurvstat, thisnotcensored;
	Doub thisnumsubset, thisdenomsubset;


	for (int i=0;i<kkdn;i++) {
		thisbmkrval = simdata[i][0];
		thistrt = simdata[i][1];
		thistime = simdata[i][2];
		thissurvstat = simdata[i][3];
		thisdeltahat = Deltahat(thisbmkrval);

		Bool condition1 = ((thistime > t0) || (thissurvstat == 1));
		Bool condition2 = (thisdeltahat > 0);
		Bool condition3 = (thistrt == 0);
		if (condition1 && condition2 && condition3) {
			mydenominator = mydenominator + 1.0;
		}
		Bool condition4 = (thistime < t0);
		if (condition1 && condition2 && condition3 && condition4) {
			mynumerator = mynumerator + 1.0;
		}
	}
	myreturn = mynumerator/mydenominator;

	return myreturn;
}


Doub JanesTheta::hatBsubneg() {
	Doub myreturn;
	myreturn =  ProbDis1givenAis1andDeltaYunder0() - ProbDis1givenAis0andDeltaYunder0();
	return myreturn;
}

Doub JanesTheta::hatBsubpos() {
	Doub myreturn;
	myreturn =  ProbDis1givenAis0andDeltaYover0() - ProbDis1givenAis1andDeltaYover0();
	return myreturn;
}


Doub JanesTheta::phatsubneg( ) {
	Doub myreturn;
	Doub mynum= 0, mydenom = 0;
	Doub thistime, thissurvstat, thisbmkrval, thiscensored, thisnotcensored;

	for (int i=0;i<kkdn;i++) {
		thistime = simdata[i][2];
		thissurvstat = simdata[i][3];
		thisbmkrval = simdata[i][0];

		thiscensored = (thistime<5)*(thissurvstat==0);
		thisnotcensored = 1.-thiscensored;

		mydenom = thisnotcensored + mydenom;
		mynum = mynum + thisnotcensored*(Deltahat(thisbmkrval)<0);
	}
	myreturn = mynum/mydenom;

	return myreturn;
}

Doub JanesTheta::phatsubpos( ) {
	Doub myreturn;
	myreturn = 1. - phatsubneg();
	return myreturn;
}

Doub JanesTheta::theta1hat( ) {
	Doub myreturn;
	myreturn = hatBsubneg() * phatsubneg();
	// myreturn = max(myreturn,0.0);
	return myreturn;
}

Doub JanesTheta::theta0hat( ) {
	Doub myreturn;
	myreturn = hatBsubpos() * phatsubpos();
	// myreturn = max(myreturn,0.0);
	return myreturn;
}


Doub JanesTheta::ProbDis1givenAis1andDeltaYunder0( ) {
	Doub myreturn;

	double mynumerator = 0.0, mydenominator = 0.0;
	Doub thisdeltahat, thisbmkrval, thistime, thistrt, thissurvstat, thisnotcensored;
	Doub thisnumsubset, thisdenomsubset;
	kkdn = simdata.nrows();

	for (int i=0;i<kkdn;i++) {
		thisbmkrval = simdata[i][0];
		thistrt = simdata[i][1];
		thistime = simdata[i][2];
		thissurvstat = simdata[i][3];
		thisdeltahat = Deltahat(thisbmkrval);

		Bool condition1 = ((thistime > t0) || (thissurvstat == 1));
		Bool condition2 = (thisdeltahat < 0);
		Bool condition3 = (thistrt == 1);
		if (condition1 && condition2 && condition3) {
			mydenominator = mydenominator + 1.0;
		}
		Bool condition4 = (thistime < t0);
		if (condition1 && condition2 && condition3 && condition4) {
			mynumerator = mynumerator + 1.0;
		}
	}
	myreturn = mynumerator/mydenominator;

	return myreturn;
}

Doub JanesTheta::ProbDis1givenAis0andDeltaYunder0( ) {
	Doub myreturn;

	double mynumerator = 0.0, mydenominator = 0.0;
	Doub thisdeltahat, thisbmkrval, thistime, thistrt, thissurvstat, thisnotcensored;
	Doub thisnumsubset, thisdenomsubset;
	kkdn = simdata.nrows();

	for (int i=0;i<kkdn;i++) {
		thisbmkrval = simdata[i][0];
		thistrt = simdata[i][1];
		thistime = simdata[i][2];
		thissurvstat = simdata[i][3];
		thisdeltahat = Deltahat(thisbmkrval);

		Bool condition1 = ((thistime > t0) || (thissurvstat == 1));
		Bool condition2 = (thisdeltahat < 0);
		Bool condition3 = (thistrt == 0);
		if (condition1 && condition2 && condition3) {
			mydenominator = mydenominator + 1.0;
		}
		Bool condition4 = (thistime < t0);
		if (condition1 && condition2 && condition3 && condition4) {
			mynumerator = mynumerator + 1.0;
		}
	}
	myreturn = double(mynumerator) / double(max(mydenominator,0.01));

	return myreturn;
}



MatDoub JanesTheta::DoTheRefinedBootstrap (int nboots) {
	JanesTheta resubjt, bootjt, bootjtresamp;
	Doub A0, A1, B0, B1, C0, C1, D0, D1;
	MatDoub myreturn(2,nboots);
	// first row is theta0s, second row is theta1s
	coxlikelihood tempcoxlikelihood;
	MatDoub BootstrapSimData(kkdn,4);
	int i, thisboot, thisindex;
	VecDoub betaMLE(4);
	Bool tempbool;

	// RESUBSTITUTION ESTIMATES
	tempcoxlikelihood.simdata = simdata;
	betaMLE[0] = 0; betaMLE[1] = 0; betaMLE[2] = 0; betaMLE[3] = 0;
	newt(betaMLE,tempbool,tempcoxlikelihood);
	resubjt.simdata = simdata; // the full dataset is used
	resubjt.setJanesTheta(simdata,t0,betaMLE[0],betaMLE[1],betaMLE[2],betaMLE[3]);
	   // the beta estimates are taken from the bootstrap fit
	D0 = resubjt.theta0hat();
	D1 = resubjt.theta1hat();

    // END OF RESUBSTITUTEION ESTIMATES


	for (thisboot=0;thisboot<nboots;thisboot++) {

// 1) Draw a bootstrap sample and develop a model.
		// Create bootstrap dataset
		for (i=0;i<kkdn;i++) {
			thisindex = rand() % kkdn;
			for (int j=0;j<4;j++) {
				BootstrapSimData[i][j] = simdata[thisindex][j];
			}
		}

		// Maximum likelihood Estimation
		tempcoxlikelihood.simdata = BootstrapSimData;
		betaMLE[0] = 0; betaMLE[1] = 0; betaMLE[2] = 0; betaMLE[3] = 0;
		newt(betaMLE,tempbool,tempcoxlikelihood);

// 2) Apply the model to all n samples in the dataset to estimate theta0 and theta1,
		// call these A0 and A1.
		bootjt.simdata = simdata; // the full dataset is used
		bootjt.setJanesTheta(simdata,t0,betaMLE[0],betaMLE[1],betaMLE[2],betaMLE[3]);
		   // the beta estimates are taken from the bootstrap fit
		A0 = bootjt.theta0hat();
		A1 = bootjt.theta1hat();

// 3) Apply the model just to the samples used to train the model (resubstitution error),
		//  call it B0 and B1.
		bootjtresamp.simdata = BootstrapSimData;
		bootjtresamp.setJanesTheta(BootstrapSimData,t0,betaMLE[0],betaMLE[1],betaMLE[2],betaMLE[3]);
		B0 = bootjtresamp.theta0hat();
		B1 = bootjtresamp.theta1hat();

// 4) The overoptimism bias is C = A-B.
		C0 = A0 - B0;
		C1 = A1 - B1;

// 5) Now develop a model on the full dataset and estimate the error with resubstitution, D0 and D1.

		// THIS IS ABOVE

// 6) The final refined bootstrap estimate is D+C.

		myreturn[0][thisboot] = D0+C0;
		myreturn[1][thisboot] = D1+C1;
	}  // END OF BOOTSTRAP LOOP

	return myreturn;
}


MatDoub JanesTheta::PercentileBootstrap ( VecDoub_I percentiles, int nboots ) {
	MatDoub myreturn(2,2);
	MatDoub BSresults(2,nboots);
	BSresults = DoTheRefinedBootstrap(nboots);
	VecDoub tempvec0(nboots);
	VecDoub tempvec1(nboots);
	int i;

	double thisfullsampleTheta0 = theta0hat();
    double thisfullsampleTheta1 = theta1hat();

	for (i=0;i<nboots;i++) {
		tempvec0[i] = BSresults[0][i];
		tempvec1[i] = BSresults[1][i];
	}

    // percentile bootstrap code;
	int rank1 = round(0.025 * nboots);
	int rank2 = round(0.975 * nboots);

	myreturn[0][0] = select(rank1, tempvec0);
	myreturn[0][1] = select(rank2, tempvec0);

	myreturn[1][0] = select(rank1, tempvec1);
	myreturn[1][1] = select(rank2, tempvec1);

	return  myreturn;
}

VecDoub mySLR(VecDoub xs,VecDoub ys) {
	VecDoub myreturn(2);

	Doub myintercept, myslope;
	Doub xsSum = 0.0, ysSum  = 0.0;
	Doub xsSsq = 0.0, ysSsq  = 0.0, xySsq = 0;
	Int n = xs.size();

	for (int i=0;i<n;i++) {
		xsSum = xsSum + xs[i];
		ysSum = ysSum + ys[i];
		xsSsq = xsSsq + pow(xs[i],2);
		ysSsq = ysSsq + pow(ys[i],2);
		xySsq = xySsq + xs[i] * ys[i];
	}
	Doub xsMean = xsSum/double(n);
	Doub ysMean = ysSum/double(n);

	myslope = (xySsq - (double(n)*xsMean*ysMean))/(xsSsq - (n * pow(xsMean,2)));
	myintercept = ysMean - myslope * xsMean;

	myreturn[0] = myintercept;
	myreturn[1] = myslope;

	return myreturn;
}


// main function -
// where the execution of program begins
int main()
{
    // prints basic in and out
   double userk1, userk2, userk3, userk4, usert0, userpropdied, userwidth;
    cout<<"What is k1 (e.g. 0.25)?" << endl;
    cin >> userk1;
    while (userk1 > 1.0 || userk1 < 0.0) {
    	cout << "k1 must be between 0 and 1, please re-enter " << endl;
    	cin >> userk1;
    }
    cout<<"What is k2 (e.g. 0.75)?" << endl;
    cin >> userk2;
    while (userk2 > 1.0 || userk2 < 0.0) {
    	cout << "k2 must be between 0 and 1, please re-enter " << endl;
    	cin >> userk2;
    }
    cout<<"What is k3 (e.g. 0.75)?" << endl;
    cin >> userk3;
    while (userk3 > 1.0 || userk3 < 0.0) {
    	cout << "k3 must be between 0 and 1, please re-enter " << endl;
    	cin >> userk3;
    }
    cout<<"What is k4 (e.g. 0.25)?" << endl;
    cin >> userk4;
    while (userk4 > 1.0 || userk4 < 0.0) {
    	cout << "k4 must be between 0 and 1, please re-enter " << endl;
    	cin >> userk4;
    }
    cout<<"What is t0 (e.g. 5 years)?" << endl;
    cin >> usert0;
    cout << "What is the proportion of deaths expected? " << endl;
    cin >> userpropdied;
    while (userpropdied > 1.0 || userpropdied < 0.0) {
    	cout << "Proportion must be between 0 and 1, please re-enter " << endl;
    	cin >> userpropdied;
    }
    cout << "What is the targeted confidence interval width? " << endl;
    cin >> userwidth;
    while (userwidth > 1.0 || userwidth < 0.0) {
    	cout << "Targeted confidence interval width must be between 0 and 1, please re-enter " << endl;
    	cin >> userwidth;
    }

    Int thetaoneorzero;
    cout << "Are you interested in theta0 (enter 0) or theta1 (enter 1)? " << endl;
    cin >> thetaoneorzero;

   	cout << "Calculating ... " << endl;
   	cout << endl;

    VecDoub coxbeta(4);
	Doub coxlambda;
	Doub t0;
	Doub propdied;
    int kkdn = 600;
	Doub k1,k2,k3,k4;

    k1 = userk1;
    k2 = userk2;
    k3 = userk3;
    k4 = userk4;
    t0 = usert0;
    propdied = userpropdied;



    srand (time(NULL));
    // srand(1);

    coxbeta[1] = coxbeta1fun(k1,k2);
    // cout << "Beta1fun test is: "  << coxbeta[1] << endl;

    coxlambda = coxlambdafun(k1,k2,t0);
    // cout << "coxlambdafun test is: "  << coxlambda << endl;

    coxbeta[0] = log(coxlambda);
    // cout << "beta0hatlogit test is: "  << coxbeta[0] << endl;

    coxbeta[3] = coxbeta3fun(k1,k2,k3,k4,t0) ;
    // cout << "beta3fun test is: "  << coxbeta[3] << endl;

    coxbeta[2] = coxbeta2fun(k1,k2,k3,k4,t0) ;
    // cout << "beta2fun test is: "  << coxbeta[2] << endl;


    MatDoub simdata(kkdn,4);
    simdata = generatedata(coxlambda, coxbeta[1], coxbeta[2], coxbeta[3], kkdn, propdied);



    std::ofstream output("generateddataset.txt");
    for (int i=0;i<kkdn;i++)
    {
    	for (int j=0;j<4;j++){
			output << simdata[i][j] << " "; // behaves like cout - cout is also a stream
    	}
    	output << "\n";
	}
	output.close();


	Bool mybool;
    VecDoub myresult(4);


	int nboots = 200;

   	VecDoub mypercentiles(2);
   	mypercentiles[0] = 0.025;
   	mypercentiles[1] = 0.975;
   	MatDoub mypctlboot(2,2);


	int mcruns = 50;
   	VecInt kkdns(5);
   	VecDoub sumalltheta0widths(5);
   	for (int temp=0;temp<5;temp++) { sumalltheta0widths[temp] = 0; }
   	VecDoub sumalltheta1widths(5);
   	for (int temp=0;temp<5;temp++) { sumalltheta1widths[temp] = 0; }


   	for (int zip=0;zip<5;zip++) { kkdns[zip] = 200*(zip+1); }

   	Doub thistheta0width, thistheta1width;

   	time_t start, end;

   	time(&start);

   	MatDoub CImatrix(mcruns,5*4);

    std::ofstream theta0file("theta0hats.txt");
    std::ofstream theta1file("theta1hats.txt");

    VecInt n0s(5);
    VecInt n1s(5);
    for (int i=0;i<5;i++) {
    	n0s[i] = 0;
    	n1s[i] = 0;
    }

    cout << "Now running preliminary simulations." << endl;

   	for (int sampsizeind=0;sampsizeind<5;sampsizeind++) {
   		cout << "Now running simulations for sample size: " << kkdns[sampsizeind];
   		cout << "\n";
   		MatDoub simdata(kkdns[sampsizeind],4);
   		for (int final=0;final<mcruns;final++) {
   			simdata = generatedata(coxlambda, coxbeta[1], coxbeta[2], coxbeta[3], kkdns[sampsizeind], propdied);

		   	myresult[0]=0;
		   	myresult[1]=0;
		   	myresult[2]=0;
		   	myresult[3]=0;
   			logisticlikelihood logisticlikelihoodobject;
   			logisticlikelihoodobject.simdata = simdata;
   			logisticlikelihoodobject.t0 = t0;
   		    newt(myresult,mybool,logisticlikelihoodobject);

   			JanesTheta myjt;
   			myjt.simdata = simdata;
   			myjt.setJanesTheta(simdata,t0,myresult[0],myresult[1],myresult[2],myresult[3]);
   			theta0file << myjt.theta0hat() <<  "\n";
   			theta1file << myjt.theta1hat() <<  "\n";

   			mypctlboot = myjt.PercentileBootstrap(mypercentiles, nboots);
   			thistheta0width = mypctlboot[0][1] - mypctlboot[0][0];
   			thistheta1width = mypctlboot[1][1] - mypctlboot[1][0];

   			CImatrix[final][0+4*sampsizeind] = mypctlboot[0][0];
   			CImatrix[final][1+4*sampsizeind] = mypctlboot[0][1];
   			CImatrix[final][2+4*sampsizeind] = mypctlboot[1][1];
   			CImatrix[final][3+4*sampsizeind] = mypctlboot[1][1];

   			if (!isnan(thistheta0width)) {
   				n0s[sampsizeind] = n0s[sampsizeind] + 1;
   			   	sumalltheta0widths[sampsizeind] = thistheta0width + sumalltheta0widths[sampsizeind];
   			}
   			if (!isnan(thistheta1width)){
   			    n1s[sampsizeind] = n1s[sampsizeind] + 1;
   			   	sumalltheta1widths[sampsizeind] = thistheta1width + sumalltheta1widths[sampsizeind];
   			}
   		}

   	}
   	theta0file.close();
   	theta1file.close();

   	VecDoub finaltheta0widths(5);
   	for (int temp=0;temp<5;temp++) { finaltheta0widths[temp] = sumalltheta0widths[temp]/double(n0s[temp]); }
   	VecDoub finaltheta1widths(5);
   	for (int temp=0;temp<5;temp++) { finaltheta1widths[temp] = sumalltheta1widths[temp]/double(n1s[temp]); }


   VecInt validwidthYN(5);
   for (int i=0;i<5;i++) {
   	if ( (!isnan(finaltheta0widths[i])) && (!isnan(finaltheta1widths[i])) )
   		validwidthYN[i] = 1;
   	else validwidthYN[i] = 0; // should identify missing data as zero
   }
   Int howmanyvalid = 0;
   for (int i=0;i<5;i++) { howmanyvalid = howmanyvalid + validwidthYN[i]; }
   VecInt ValidIndices(howmanyvalid);
   Int mycurrent = 0;
   for (int j=0;j<5;j++) {
   	 	if(validwidthYN[j] == 1) {
   	 		ValidIndices[mycurrent] = j;
   	 		mycurrent = mycurrent + 1;
   	 	}
   	}

   	VecDoub myinvsqwidth0(howmanyvalid),myinvsqwidth1(howmanyvalid),revkkdns(howmanyvalid);
   	Int thisoneof5;
   	for (int i=0;i<howmanyvalid;i++) {
   		thisoneof5 = ValidIndices[i];
   		myinvsqwidth0[i] = 1/pow(finaltheta0widths[thisoneof5],2);
   		myinvsqwidth1[i] = 1/pow(finaltheta1widths[thisoneof5],2);
   		revkkdns[i] = kkdns[thisoneof5];
   	}

   	VecDoub SLRresult0(2), SLRresult1(2);

   	SLRresult0 = mySLR(myinvsqwidth0,revkkdns);
   	SLRresult1 = mySLR(myinvsqwidth1,revkkdns);

   	double slope0 = SLRresult0[1];
   	double intercept0 = SLRresult0[0];
   	double slope1 = SLRresult1[1];
   	double intercept1 = SLRresult1[0];
   	double mytarget = userwidth;
   	double samplesizeest0 = intercept0 + (slope0 * pow(1.0/mytarget,2));
   	double samplesizeest1 = intercept1 + (slope1 * pow(1.0/mytarget,2));

    // if (thetaoneorzero == 0) {
    // 	cout << "thistheta0width : " << finaltheta0widths[0] << " , " <<finaltheta0widths[1] << " , " <<finaltheta0widths[2] << " , " <<finaltheta0widths[3] << " , " <<finaltheta0widths[4] << endl ;
    // 	cout << "sample size estimate for theta0 is: " << round(samplesizeest0) << endl;
    // }
    // if (thetaoneorzero == 1) {
    //  	cout << "thistheta1width : " << finaltheta1widths[0] << " , " <<finaltheta1widths[1] << " , " <<finaltheta1widths[2] << " , " <<finaltheta1widths[3] << " , " <<finaltheta1widths[4] << endl ;
    // 	cout << "sample size estimate for theta1 is: " << round(samplesizeest1) << endl;
    // }


/////////////////////
// START FINE-TUNING
////////////////////

   	VecInt n4theta0localsearch(2);
   	Int theta0center;
   	theta0center = round(samplesizeest0);
   	Int Gap= round(theta0center/5); //40;

   	for (int i=0;i<2;i++) {
   		n4theta0localsearch[i] = theta0center - (Gap/2) + Gap*i;
 //  		cout << "Sample sizes for local interpolation is: " << n4theta0localsearch[i] << endl;
   	}
   	Int finetunemcruns = 400;
   	MatDoub localCImatrix(finetunemcruns,2*4);

   	VecDoub FTsumalltheta0widths(2);
   	for (int temp=0;temp<2;temp++) { FTsumalltheta0widths[temp] = 0; }
   	VecDoub FTsumalltheta1widths(2);
   	for (int temp=0;temp<2;temp++) { FTsumalltheta1widths[temp] = 0; }


   	VecInt FTn0s(2);
    VecInt FTn1s(2);
    for (int i=0;i<2;i++) {
    	FTn0s[i] = 0;
    	FTn1s[i] = 0;
    }


    cout << "Now running fine tuning simulations." << endl;

   	for (int sampsizeind=0;sampsizeind<2;sampsizeind++) {
   		cout << "Now running simulations for sample size: " << n4theta0localsearch[sampsizeind];
   		cout << "\n";
   		MatDoub simdata(n4theta0localsearch[sampsizeind],4);
   		for (int final=0;final<finetunemcruns;final++) {
   			simdata = generatedata(coxlambda, coxbeta[1], coxbeta[2], coxbeta[3], n4theta0localsearch[sampsizeind], propdied);

		   	myresult[0]=0;
		   	myresult[1]=0;
		   	myresult[2]=0;
		   	myresult[3]=0;
   			logisticlikelihood logisticlikelihoodobject;
   			logisticlikelihoodobject.simdata = simdata;
   			logisticlikelihoodobject.t0 = t0;
   		    newt(myresult,mybool,logisticlikelihoodobject);

   			JanesTheta myjt;
   			myjt.simdata = simdata;
   			myjt.setJanesTheta(simdata,t0,myresult[0],myresult[1],myresult[2],myresult[3]);

   			mypctlboot = myjt.PercentileBootstrap(mypercentiles, nboots);
   			thistheta0width = mypctlboot[0][1] - mypctlboot[0][0];
   			thistheta1width = mypctlboot[1][1] - mypctlboot[1][0];

   			localCImatrix[final][0+4*sampsizeind] = mypctlboot[0][0];
   			localCImatrix[final][1+4*sampsizeind] = mypctlboot[0][1];
   			localCImatrix[final][2+4*sampsizeind] = mypctlboot[1][1];
   			localCImatrix[final][3+4*sampsizeind] = mypctlboot[1][1];

   			if (!isnan(thistheta0width)) {
   				FTn0s[sampsizeind] = FTn0s[sampsizeind] + 1;
   			   	FTsumalltheta0widths[sampsizeind] = thistheta0width + FTsumalltheta0widths[sampsizeind];
   			}
   			if (!isnan(thistheta1width)){
   			    FTn1s[sampsizeind] = FTn1s[sampsizeind] + 1;
   			   	FTsumalltheta1widths[sampsizeind] = thistheta1width + FTsumalltheta1widths[sampsizeind];
   			}
   		}

   	}


   	VecDoub FTfinaltheta0widths(2);
   	for (int temp=0;temp<2;temp++) { FTfinaltheta0widths[temp] = FTsumalltheta0widths[temp]/double(FTn0s[temp]); }
   	VecDoub FTfinaltheta1widths(2);
   	for (int temp=0;temp<2;temp++) { FTfinaltheta1widths[temp] = FTsumalltheta1widths[temp]/double(n1s[temp]); }

   	for (int i=0;i<2;i++) {
//   		cout << "mean width for theta0 is " << FTfinaltheta0widths[i] << " when n is " << n4theta0localsearch[i] << endl;
   	}


   	time(&end);
   	cout << "timing is : " << double((end-start)/60.) << " minutes " << endl;

/////////////////////
// END FINE-TUNING
////////////////////



   	if (howmanyvalid < 5) {
   		cout << "Warning: only sample sizes " ;
   		for (int i=0;i<howmanyvalid;i++) {
   			cout << revkkdns[i] << ",  " ;
   		}
   		cout << "used." << endl;
   	}

   	cout << '\a';

   	Int SampleSizeEstimate;
   	SampleSizeEstimate = ceil(n4theta0localsearch[0] +
   		(n4theta0localsearch[1]-n4theta0localsearch[0]) *
   		(userwidth-FTfinaltheta0widths[0])/
   		(FTfinaltheta0widths[1]-FTfinaltheta0widths[0]));

   	cout << "Sample size estimate is: " << SampleSizeEstimate << "." << endl;


    return 0;

}
