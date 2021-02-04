int hw5(){

  TH2D * Histo2d = new TH2D("Histo2d","r for x",1000,1,4,1000,0,1);
    const int Niterations=200;
    double values[Niterations+1];
  values[0] = 0.2;

  for (double r=1;r<=4;r=r+3./1000.){
      
      	for (int i=0;i<Niterations;++i) {
    values[i+1]=recurrenceRelation(values[i],r);
	
    if (i>100){
      Histo2d->Fill(r,values[i]);
    }
	}
  }
  Histo2d->Draw("col");

return 0; 
}

double recurrenceRelation(double x,double r) {
return r*x*(1-x);
}
//r will bifurcate at a 2**n rate
