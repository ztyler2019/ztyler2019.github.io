// conservation of string should be r+l=1
//This condition does not work in all cases, the string seems to strech in certain perameters and the coffee cup will push the string
//we could edit the program so that when l is decreasing, there should be a condition where r+l<=1
#include<algorithm>

void derivative(double x, double y[], double f[]) {
 
  double g=9.81;
  double b_wind = 0;
  double m=1;
  double M = 2;
  double r = y[0];
  double dr = y[1];
  double th = y[2];
  double dth = y[3];
 
  // Since the derivative of y[0] (i.e. f[0]) is v (i.e. y[1]), we 
  // simply set f[0] equal to y[1]
  f[0] = y[1]; 
  f[1] = (1/(M+m))*(m*r*pow(dth,2) - M*g + m*g*sin(th));
  f[2] = y[3];
  f[3] = -(2/r)*dr*dth + (g/r)*cos(th);
  return;
}

void rk4_step(double* y, double* ydiff, double* y_old, int nODE, double step, double x0) {
   
  // Setup of constans for the RK method.
  double const a2(1./5.), a3(3./10), a4(3./5.), a5(1.0), a6(7./8.);
  double const b21(1./5.);
  double const b31(3./40.),       b32(+9./40.);
  double const b41(3./10.),       b42(-9./10.),   b43(6./5.);
  double const b51(-11./54.),     b52(5./2.),     b53(-70./27),     b54(-35./27.);
  double const b61(1631./55296.), b62(175./512.), b63(575./13824.), b64(44275./110592.), b65(253./4096.);

  // coefficients for calculating y(nth-order) (i.e. c's in Eq.8.83)
  double const c[6] = {37./378., 0, 250./621., 125./594., 0, 512./1771.};
  // coefficients for calculating y(n-1-th order) (i.e. c*'s in Eq. 8.84)
  double const cs[6] = {2825./27648., 0, 18575./48384., 13525./55296., 277./14336., 1./4.};

  vector<double> f[6];
  for (int iFunc=0;iFunc<6;++iFunc) {
    for (int iN=0; iN<nODE; ++iN) {
      f[iFunc].push_back(0.0);
    }
    // at this point the size of each vector in f[iFunc] should be 2
    //cout << "Size of vector #" << iFunc << " is " << f[iFunc].size() << endl;
  }
  // 1st step
  // With x0, y_old, calculate derivatives f[0]
  derivative(x0,y_old,&f[0][0]);
  // With the f0 derivatives, calculate new y's.
  for (int i1=0;i1<nODE;++i1) {
    y[i1]=y_old[i1]+b21*step*f[0][i1];
  }
  //cout << "Step 1, y and y': " << y[0] << ", " << y[1] << endl;
  // Calculate new x.
  double x = x0 + a2*step;
  
  // 2nd step
  // With the x, and y's just calculated, calculate derivatives f[1]
  derivative(x,y,&f[1][0]);
  // With the f0 and f1, derivatives, calculate new y's.
  for (int i2=0;i2<nODE;++i2) {
    y[i2]=y_old[i2]+b31*step*f[0][i2]+b32*step*f[1][i2];    
  }
  //cout << "Step 2, y and y': " << y[0] << ", " << y[1] << endl;
  // Calculate new x.
  x = x0 + a3*step;
  
  //3d step
  // With the x, and y's just calculated, calculate derivatives f[2]
  derivative(x,y,&f[2][0]);
  // With the f0, f1, and f2 derivatives, calculate new y's.
  for (int i3=0;i3<nODE;++i3) {
    y[i3]=y_old[i3]+b41*step*f[0][i3]+b42*step*f[1][i3]+b43*step*f[2][i3];    
  }
  //cout << "Step 3, y and y': " << y[0] << ", " << y[1] << endl;
  // Calculate new x.
  x = x0 + a4*step;

  //4th step
  // With the x, and y's just calculated, calculate derivatives f[3]
  derivative(x,y,&f[3][0]);
  // With the f0, f1, ,f2 and f3 derivatives, calculate new y's.
  for (int i4=0;i4<nODE;++i4) {
    y[i4]=y_old[i4]+b51*step*f[0][i4]+b52*step*f[1][i4]+b53*step*f[2][i4]+b54*step*f[3][i4];    
  }
  //cout << "Step 4, y and y': " << y[0] << ", " << y[1] << endl;
  // Calculate new x.
  x = x0 + a5*step;

  //5th step
  // With the x, and y's just calculated, calculate derivatives f[4]
  derivative(x,y,&f[4][0]);
  // With the f0, f1, ,f2, f3 and f4 derivatives, calculate new y's.
  for (int i5=0;i5<nODE;++i5) {
    y[i5]=y_old[i5]+b61*step*f[0][i5]+b62*step*f[1][i5]+b63*step*f[2][i5]+b64*step*f[3][i5]+b65*step*f[4][i5];    
  }
  //cout << "Step 5, y and y': " << y[0] << ", " << y[1] << endl;
  // Calculate new x.
  x = x0 + a6*step;

  //6th step
  // With the x, and y's just calculated, calculate derivatives f[5]
  derivative(x,y,&f[5][0]);
  // With the f0, f1, ,f2, f3, f4 and f5 derivatives, calculate final y's. 
  // There will be two calculations, one using the c's and one using the c*'s.
  
  // create a vector to hold the 2nd solution
  vector<double> yp(nODE);
  for (int i6=0;i6<nODE;++i6) {
    // first add the old values to both solutions.
    y[i6]=y_old[i6];
    yp[i6]=y_old[i6];
    // next add the sum of the c[i]*f[i] for all 6 products
    // do the same for the solution using c*'s.
    for (int j1=0; j1<6; ++j1) {
      y[i6]+=step*c[j1]*f[j1][i6];
      yp[i6]+=step*cs[j1]*f[j1][i6];
    }
    // then calculate the difference between the solutions.
    ydiff[i6]=yp[i6]-y[i6];
  }
  //cout << "Step 6, y and y': " << y[0] << ", " << y[1] << endl;

  return;
}

void rk4_stepper(double* y_old, const int nODE, double xstart, double xmax, double hstep, double eps, int& nxstep) {

  cout << "nODE   : " << nODE <<  endl;
  cout << "xstart : " << xstart << endl;
  cout << "xmax   : " << xmax   << endl;
  cout << "hstep  : " << hstep << endl;
  cout << "eps    : " << eps   << endl;
  
  
  double heps=hstep*eps;
  double yerrmax=0.99; // max error in any given step
  double redPow = 1./5.; //power used in formula to reduce the step size
  double esmall = eps/100.; // lower limit of precision, if the difference is smaller, increase step size.

  vector<double> ydiff(nODE); // store the difference between the two R-K methods.
  vector<double> y(nODE);  // store the result of the calculation of the new dep. variable at each step.
  double hnew;
  double x0 = xstart; // value of the independent variable at each step. Initialize to xstart.

  double step_lolim = hstep/1000.; // lower limit for step size
  double step_hilim = hstep*10.; // upper limit for step size
  
  vector<double> xVals,yVals[nODE]; // vector of values of the indep and dep variables.  These are the values that will go into drawing the graphs of y(x) and y'(x), etc.
  
  // First, store the starting values.
  xVals.push_back(x0);
  for (int i=0; i<nODE; ++i) {
    yVals[i].push_back(y_old[i]); // store into final vector of results
    y[i] = y_old[i]; // store into array for current step, maybe not needed?
  }

  // precision.
  while (x0<xmax) {
    yerrmax=0.99; //reset max error for each step in the loop
    while (yerrmax<=0.99) { // when error is of this size, we found the necessary step size.
  
      rk4_step(&y[0], &ydiff[0], &y_old[0], nODE, hstep, x0);
      
       for (int j=0;j<nODE;++j) {
	yerrmax=max(yerrmax,fabs(heps/ydiff[j]));
	//cout << "Error for Eq# " << j << " : " << fabs(heps/ydiff[j]) << endl;
      }
  
      // or greater than 1 if ydiff is small (we're done)
      if (yerrmax==0.99) {
	// reduce the step size.
	hstep=0.5*hstep*pow(fabs(heps/ydiff[0]),redPow);
	//heps=hstep*eps;
	cout << "Reducing step size to  " << hstep << endl;
      }
      
      // We will repeat the calculation, unless the step size is too small.
      // Protect against it.
      if (hstep<step_lolim) {
	cout << "rk4_stepper: step size requested is smaller than lowest stepsize allowed" << endl;
	cout << "Fix by trying a lower step size of lowering the step_lolim" << endl;
	cout << "Exiting..." << endl;
	return;
      }
    }
    
    if (yerrmax>1./esmall) hstep*=2.;
    if (hstep>step_hilim) hstep=step_hilim;
         
    for (int k=0;k<nODE; ++k) {
      y_old[k]=y[k];
      yVals[k].push_back(y_old[k]); // store into final vector of results for dep. variable.
    }
    x0+=hstep; // propagate the indep. varialbe by one step.
    xVals.push_back(x0); // store into vector of indep variable values.
    nxstep++;
   
  }
  cout << "Total steps " << nxstep << endl;
  cout << "Size of x vector " << xVals.size() << endl;
  cout << "Size of y vector " << yVals[0].size() << endl;
  
  // Make TGraphs out of the output vectors.
  
  vector <double>  x, y, vx, vy;
 
  for (int i=0; i<yVals[0].size();i++) {
    x.push_back(-yVals[0][i]*cos(yVals[2][i]));
    y.push_back(-yVals[0][i]*sin(yVals[2][i]));
    vx.push_back(-yVals[1][i]*cos(yVals[2][i])+yVals[0][i]*yVals[3][i]*sin(yVals[2][i]));
    vy.push_back(-yVals[1][i]*sin(yVals[2][i])-yVals[0][i]*yVals[3][i]*cos(yVals[2][i]));
 }
  int end = 0;
  for (int i=0; i<yVals[0].size();i++) {
    if (yVals[1][i]> 0) {
      end = i;
      break;
    }
  }
  TGraph* yGraph = new TGraph(end +1,&(x[0]),&(y[0]));
  // Make the graphs look better, add color
  yGraph->SetLineColor(2);
  
  yGraph->SetMinimum(-1);
  yGraph->SetMaximum(+1);
  // Plot them in a canvas
  TCanvas* freeFallAdapStep = new TCanvas("freeFallAdapStep","Adaptive Step Size RK",500,500);
  yGraph->Draw("AL");
 
  yGraph->GetHistogram()->SetXTitle("x (m)");
  yGraph->GetHistogram()->SetYTitle("y (m)");
 
 
  gifmaker(x,y,xVals,vx,vy);

 return;

 
} 

void freeFallAdaptiveStep() {
  
  // Set up parameters for running the Runge-Kutta stepper
  // routine.
  
  //initial position and initial velocity.
  const int numberOfDiffEqs(4);

  double y_old[numberOfDiffEqs]; 
  y_old[0]=1; y_old[1]=0; y_old[2]=0; y_old[3]=0;
  
  
  // Initial and final values of the indep. variable
  double xstart = 0;
  double xmax = 2;
  
  // Step size and desired precision eps (given as percent of step size).
  // i.e. 0.1=10% of hstep will be used to compare against the difference
  // in the two solutions given by the 5th and 6th order r-k.

  double hstep = 0.005;
  double eps = 0.1; 
  int NumberOfSteps;
  rk4_stepper(y_old, numberOfDiffEqs, xstart, xmax, hstep, eps, NumberOfSteps);
  
  return;
}


void gifmaker(vector <double> x, vector <double> y, vector <double> t, vector <double> vx, vector <double> vy){
  TCanvas* c1 = new TCanvas("gif","Adaptive Step Size RK",500,500);
  for (int i = 0;i < x.size();i++) {
    TGraph* gif = new TGraph(i+1,&x[0],&y[0]);
    gif->Draw ("AL");
 gif->SetMinimum(-1);
  gif->SetMaximum(+1);
  gif->GetXaxis()->SetLimits(-1,1);
  gif->SetLineColor(kRed);
  TMarker m (x[i],y[i],20);
  m.Draw();
  TLine r (x[i],y[i],0,0);
  r.Draw();
  TMarker M (0,sqrt(pow(x[i],2)+pow(y[i],2))-1,21);
  M.Draw();
  TLine l (0,sqrt(pow(x[i],2)+pow(y[i],2))-1,0,0);
  l.Draw();
  TArrow v (x[i],y[i],x[i]+.05*vx[i],y[i]+.05*vy[i]); 
  v.Draw();
  v->SetLineColor(kBlue);
    c1->Update();
    c1->Modified();
    c1->SaveAs("pendulum.gif+");
  }
}
