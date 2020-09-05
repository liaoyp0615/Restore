{
  TRandom3 ran;
  // if(gRandom) delete gRandom;
  // gRandon = new TRandom3(0);
  double m_rand;
  //define polya function
  
  //  char funname[6];
  TString funname;
  TF1* polya[3];
  double par_polya_0[3]={6257, 580.3, 443.2};
  double par_polya_1[3]={27.87,21.03,27.53};
  double par_polya_2[3]={0.4946,0.4088,0.4799};
 
  for(int i=0;i<3;i++){

    funname.Form("polya%d",i);
    // sprintf(funname,"polya%d",i);
    polya[i] = new TF1(funname, "[0]*pow(1+[2],1+[2])*pow(x/[1],[2])*exp(-(1+[2])*x/[1])/TMath::Gamma(1+[2])", 0, 500);
    
    polya[i]->SetParameter(0,par_polya_0[i]);
    polya[i]->SetParameter(1,par_polya_1[i]);
    polya[i]->SetParameter(2,par_polya_2[i]);
 
  }
  double trans1(0.867), trans2(0.61), trans3(0.613);
  //-----------------------------------------------
  
  //set distribution parameters
  double mu_x1(0.3838), sigma_x1(0.175);
  double mu_z1(0), sigma_z1(0.1689);
  double mu_t1(58.48), sigma_t1(2.147);

  double mu_x2(0.3774), sigma_x2(0.1757);
  double mu_z2(0), sigma_z2(0.1693);
  double mu_t2(58.31), sigma_t2(2.146);

  double mu_x12 = mu_x1+mu_x2;
  double sigma_x12 = sqrt(sigma_x1*sigma_x1+sigma_x2*sigma_x2);
  double mu_z12 = mu_z1+mu_z2;
  double sigma_z12 = sqrt(sigma_z1*sigma_z1+sigma_z2*sigma_z2);
  double mu_t12 = mu_t1+mu_t2;
  double sigma_t12 = sqrt(sigma_t1*sigma_t1+sigma_t2*sigma_t2);

  double Drift=2.5;
  double mu_x_par[2]={0.02627,0.4785};
  double sigma_x_par[4]={0.06928,0.09882,-0.01768,0.001678};
  double mu_z_par=0;
  double sigma_z_par[4]={0.06591,0.07845,-0.01192,0.0009879};
  double mu_t_par[2]={4.985,29.15};
  double sigma_t_par[4]={0.8703,1.092,-0.1897,0.01812};
  double mu_xd,sigma_xd,mu_zd,sigma_zd,mu_td,sigma_td;
  TF1* U_xd = new TF1("U_xd","0.02627+0.4785*x",0,5);
  TF1* Sigma_xd = new TF1("Sigma_xd","0.06928+0.09882*x-0.01768*x*x+0.001678*x*x*x",0,5);
  TF1* U_zd = new TF1("U_zd","0*x",0,5);
  TF1* Sigma_zd = new TF1("Sigma_zd","0.06591+0.07845*x-0.01192*x*x+0.0009879*x*x*x",0,5);
  TF1* U_td = new TF1("U_td","4.985+29.15*x",0,5);
  TF1* Sigma_td = new TF1("Sigma_td","0.8703+1.092*x-0.1897*x*x+0.01812*x*x*x",0,5);
  double mu_xd12 = mu_x1+mu_x2+U_xd->Eval(Drift);
  double sigma_xd12 = sqrt(sigma_x1*sigma_x1+sigma_x2*sigma_x2+(Sigma_xd->Eval(Drift))*(Sigma_xd->Eval(Drift)));
  double mu_zd12 = mu_z1+mu_z2+U_zd->Eval(Drift);
  double sigma_zd12 = sqrt(sigma_z1*sigma_z1+sigma_z2*sigma_z2+(Sigma_zd->Eval(Drift))*(Sigma_zd->Eval(Drift)));
  double mu_td12 = mu_t1+mu_t2+U_td->Eval(Drift);
  double sigma_td12 = sqrt(sigma_t1*sigma_t1+sigma_t2*sigma_t2+(Sigma_td->Eval(Drift))*(Sigma_td->Eval(Drift)));
  
  double x_dmin = U_xd->Eval(Drift) - 6*Sigma_xd->Eval(Drift);
  double x_dmax = U_xd->Eval(Drift) + 6*Sigma_xd->Eval(Drift);
  double z_dmin = U_zd->Eval(Drift) - 6*Sigma_zd->Eval(Drift);
  double z_dmax = U_zd->Eval(Drift) + 6*Sigma_zd->Eval(Drift);
  double t_dmin = U_td->Eval(Drift) - 6*Sigma_td->Eval(Drift);
  double t_dmax = U_td->Eval(Drift) + 6*Sigma_td->Eval(Drift);

  double X_up = mu_xd12+6*(sigma_x1+sigma_x2+Sigma_xd->Eval(Drift));
  double X_low = mu_xd12-6*(sigma_x1+sigma_x2+Sigma_xd->Eval(Drift));
  double Z_up = mu_zd12+6*(sigma_z1+sigma_z2+Sigma_zd->Eval(Drift));
  double Z_low = mu_zd12-6*(sigma_z1+sigma_z2+Sigma_zd->Eval(Drift));
  double T_up = mu_td12+6*(sigma_t1+sigma_t2+Sigma_td->Eval(Drift));
  double T_low = mu_td12-6*(sigma_t1+sigma_t2+Sigma_td->Eval(Drift));
      
				   				
  //  cout<<"T_up = "<<T_up<<endl;
  //  cout<<"T_low = "<<T_low<<endl;
  //----------------------------------------------------

  //define histogram for gain
  double gain_min(0), gain_max(500);
  int nbin=1000;
  TH1D* h_gain1 = new TH1D("h_gain1","Gain distribution for multiplied electrons",nbin,gain_min,gain_max);
  TH1D* h_gain2 = new TH1D("h_gain2","Gain distribution for multiplied electrons",nbin,gain_min,gain_max);
  TH1D* h_gain3 = new TH1D("h_gain3","Gain distribution for multiplied electrons",nbin,gain_min,gain_max);
  for(int i=1; i<nbin; i++){
    double g=gain_min+(gain_max-gain_min)/nbin*(i+0.5);
    h_gain1->SetBinContent(i+1,polya[0]->Eval(g));
    h_gain2->SetBinContent(i+1,polya[1]->Eval(g));
    h_gain3->SetBinContent(i+1,polya[2]->Eval(g));
  }
  
  double count1 = h_gain1->Integral();
  double count2 = h_gain2->Integral();
  double count3 = h_gain3->Integral();
  double zero1 = (1-trans1)/trans1*count1;
  double zero2 = (1-trans2)/trans2*count2;
  double zero3 = (1-trans3)/trans3*count3;
  h_gain1->SetBinContent(1,zero1);
  h_gain2->SetBinContent(1,zero2);
  h_gain3->SetBinContent(1,zero3);
  // cout<<h_gain1->GetRandom()<<endl;
  
  TH1D* h_gain23 = new TH1D("h_gain23","",nbin,0,3000);
  int nGain1, xGain1, nGain2, xGain2, nGain3, xGain3;
  int nGain3_count, Gain_total, nGain23;
  for(int i=0; i<100000; i++){
    m_rand = ran.Uniform(0.0,1.0);
    if(m_rand<trans2){
      xGain2=polya[1]->GetRandom();
      nGain2=int(xGain2);
    }
    else{
      nGain2=0;
    }  
    nGain3_count=0;
    //    cout<<"nGain2 = "<<nGain2<<endl;
    for(int j=0; j<nGain2; j++){     
      m_rand=ran.Uniform(0.0,1.0);
      if(m_rand<trans3){
	xGain3=polya[2]->GetRandom(1,500);
	nGain3=int(xGain3);
      }
      else{
	nGain3=0;
      }
      nGain3_count = nGain3_count + nGain3;
    }
    h_gain23->Fill(nGain3_count);
  }
  TH1D* h_gain_total=new TH1D("h_gain_total","",200,0,50000);
  for(int i=0; i<100000; i++){
    Gain_total=0;
    m_rand=ran.Uniform(0.0,1.0);
    if(m_rand<trans1){
      xGain1=polya[0]->GetRandom(1,500);
      nGain1=int(xGain1);
    }
    else{
      nGain1=0;
    }
    for(int m=0; m<nGain1; m++){
      m_rand = ran.Uniform(0.0,1.0);
      if(m_rand<trans2){
	xGain2=polya[1]->GetRandom(1,500);
	nGain2=int(xGain2);
      }
      else{
	nGain2=0;
      }  
      
      for(int j=0; j<nGain2; j++){     
	m_rand=ran.Uniform(0.0,1.0);
	if(m_rand<trans3){
	  xGain3=polya[2]->GetRandom(1,500);
	  nGain3=int(xGain3);
	}
	else{
	  nGain3=0;
	}
	Gain_total+=nGain3;
      }
    }
    h_gain_total->Fill(Gain_total);
  }

  //--------------------------------------------------------
  
  
  // define histogram
  double x1, x2, x3, z1, z2, z3, t1, t2, t3;
  double x, x_min(X_low), x_max(X_up);
  double t, t_min(150), t_max(250);
  double z, z_min(Z_low), z_max(Z_up);
  // x_dmin=x_min;x_dmax=d_max; z_dmin=z_min;z_dmax=z_max;
  //  t_dmin=0;t_dmax=150;
  int npos=200;
  int npos1=200;
  TH1D* hx = new TH1D("hx","",npos,x_min,x_max);
  TH1D* hz = new TH1D("hz","",npos,z_min,z_max);
  TH1D* ht = new TH1D("ht","",npos,t_min,t_max);
  
  TH1D* hx_ex = new TH1D("hx_ex","",npos1,x_dmin,x_dmax);
  TH1D* ht_ex = new TH1D("ht_ex","",npos1,t_dmin,t_dmax);

  TH1D* hx1 = new TH1D("hx1","",npos,x_min,x_max);
  TH1D* hz1 = new TH1D("hz1","",npos,z_min,z_max);
  TH1D* ht1 = new TH1D("ht1","",npos,t_min,t_max);
  TH1D* hx1_temp1 = new TH1D("hx1_temp1","",npos1,x_dmin,x_dmax);
  TH1D* hz1_temp1 = new TH1D("hz1_temp1","",npos1,z_dmin,z_dmax);
  TH1D* ht1_temp1 = new TH1D("ht1_temp1","",npos1,t_dmin,t_dmax);
  TH1D* hx1_temp2 = new TH1D("hx1_temp2","",npos,x_min,x_max);
  TH1D* hz1_temp2 = new TH1D("hz1_temp2","",npos,z_min,z_max);
  TH1D* ht1_temp2 = new TH1D("ht1_temp2","",npos,t_min,t_max);

  TH1D* hx2 = new TH1D("hx2","",npos,x_min,x_max);
  TH1D* hz2 = new TH1D("hz2","",npos,z_min,z_max);
  TH1D* ht2 = new TH1D("ht2","",npos,t_min,t_max);

  TH1D* nMiltiE = new TH1D("nMultiE","",200,0,50000);
  TH1D* nMiltiE1 = new TH1D("nMultiE1","",200,0,50000);
  TH1D* nMiltiE2 = new TH1D("nMultiE2","",200,0,50000);
  
  //-------------------------------------------------------

  //set GEM1 & GEM23 histogram array
  //Note that the length of these histo arrs should be npos!
  TH1D* hx23_arr[200];
  TH1D* hz23_arr[200];
  TH1D* ht23_arr[200];
  int npos_d = 200;// number of interval for drift d(0,5).
  double d, d_min(0), d_max(5);
  double delta_d=(d_max-d_min)/npos_d;
  //Note that the length of these histo arrs should be npos_d!
  TH1D* h_gaus_x[200];
  TH1D* h_gaus_z[200];
  TH1D* h_gaus_t[200];
  //set histogram arr for GEM1 
  for(int i=0; i<npos_d; i++){
    d=d_min+(d_max-d_min)/npos_d*(i+0.5);
    h_gaus_x[i]=new TH1D("","",npos1,x_dmin,x_dmax);
    h_gaus_z[i]=new TH1D("","",npos1,z_dmin,z_dmax);
    h_gaus_t[i]=new TH1D("","",npos1,t_dmin,t_dmax);
    mu_xd=U_xd->Eval(d);
    sigma_xd=Sigma_xd->Eval(d);
    mu_zd=U_zd->Eval(d);
    sigma_zd=Sigma_zd->Eval(d);
    mu_td=U_td->Eval(d);
    sigma_td=Sigma_td->Eval(d);
    // if(i==25) cout<<"mu_td = "<<mu_td<<endl; cout<<"sigma_td = "<<sigma_td<<endl;
    for(int j=0; j<npos1; j++){
      x=x_dmin+(x_dmax-x_dmin)/npos1*(j+0.5);
      z=z_dmin+(z_dmax-z_dmin)/npos1*(j+0.5);
      t=t_dmin+(t_dmax-t_dmin)/npos1*(j+0.5);
      h_gaus_x[i]->SetBinContent(j+1,TMath::Gaus(x,mu_xd,sigma_xd,1));
      h_gaus_z[i]->SetBinContent(j+1,TMath::Gaus(z,mu_zd,sigma_zd,1));
      h_gaus_t[i]->SetBinContent(j+1,TMath::Gaus(t,mu_td,sigma_td,1));
    }
  }
    
  for(int i=0; i<npos1; i++){
    x=x_dmin+(x_dmax-x_dmin)/npos1*(i+0.5);
    z=z_dmin+(z_dmax-z_dmin)/npos1*(i+0.5);
    t=t_dmin+(t_dmax-t_dmin)/npos1*(i+0.5);
    hx23_arr[i]=new TH1D("","",npos,x_min,x_max);
    hz23_arr[i]=new TH1D("","",npos,z_min,z_max);
    ht23_arr[i]=new TH1D("","",npos,t_min,t_max);
    for(int j=0;j<npos;j++){
      double x23=x_min+(x_max-x_min)/npos*(j+0.5);
      double z23=z_min+(z_max-z_min)/npos*(j+0.5);
      double t23=t_min+(t_max-t_min)/npos*(j+0.5);
      hx23_arr[i]->SetBinContent(j+1,TMath::Gaus(x23,x+mu_x12,sigma_x12,1));
      hz23_arr[i]->SetBinContent(j+1,TMath::Gaus(z23,z+mu_z12,sigma_z12,1));
      ht23_arr[i]->SetBinContent(j+1,TMath::Gaus(t23,t+mu_t12,sigma_t12,1));
    }
  }

  // time variable------------------------
  long start(0), end(0);
  long tx1(0), tz1(0), tt1(0);
  long tx2, tz2, tt2;
  long tx(0), tz(0), tt(0);
  //------------------------------------
  // histogram for mean and RMS
  int bin=200;
  TH1D* meanX = new TH1D("meanX","",bin,0,4);
  TH1D* meanZ = new TH1D("meanZ","",bin,-0.4,0.4);
  TH1D* meanT = new TH1D("meanT","",bin,0,220);
  TH1D* meanX1 = new TH1D("meanX1","",bin,0,4);
  TH1D* meanZ1 = new TH1D("meanZ1","",bin,-0.4,0.4);
  TH1D* meanT1 = new TH1D("meanT1","",bin,0,220);
  TH1D* meanX2 = new TH1D("meanX2","",bin,0,4);
  TH1D* meanZ2 = new TH1D("meanZ2","",bin,-0.4,0.4);
  TH1D* meanT2 = new TH1D("meanT2","",bin,0,220);

  TH1D* RMSX = new TH1D("RMSX","",bin,0,0.6);
  TH1D* RMSZ = new TH1D("RMSZ","",bin,0,0.6);
  TH1D* RMST = new TH1D("RMST","",bin,0,10);
  TH1D* RMSX1 = new TH1D("RMSX1","",bin,0,0.6);
  TH1D* RMSZ1 = new TH1D("RMSZ1","",bin,0,0.6);
  TH1D* RMST1 = new TH1D("RMST1","",bin,0,10);
  TH1D* RMSX2 = new TH1D("RMSX2","",bin,0,0.6);
  TH1D* RMSZ2 = new TH1D("RMSZ2","",bin,0,0.6);
  TH1D* RMST2 = new TH1D("RMST2","",bin,0,10);

  ofstream ofile;
  ofile.open("timeRatio.txt",std::ios::app);
  d=Drift;
  int id;
  double integ;
  int nEvt=8000;
  for(int dd=0; dd<2; dd++){
  for(int iEvt=0; iEvt<nEvt; iEvt++){

    //3 level sampling--------------------------
    // x sampling ------------------
    start=clock();
    mu_xd=U_xd->Eval(d);
    sigma_xd=Sigma_xd->Eval(d);
    Gain_total=0;
    m_rand=ran.Uniform(0.0,1.0);
    if(m_rand<trans1){
      xGain1=polya[0]->GetRandom(1,500);
      nGain1=int(xGain1);
    }
    else{
      nGain1=0;
    }
    for(int i=0; i<nGain1; i++){
      x1=ran.Gaus(mu_xd,sigma_xd);
      hx_ex->Fill(x1);
      m_rand=ran.Uniform(0.0,1.0);
      if(m_rand<trans2){
	xGain2=polya[1]->GetRandom(1,500);
	nGain2=int(xGain2);
      }
      else{
	nGain2=0;
      }
      for(int j=0; j<nGain2; j++){
	x2=ran.Gaus(x1+mu_x1,sigma_x1);
	m_rand=ran.Uniform(0.0,1.0);
	if(m_rand<trans3){
	  xGain3=polya[2]->GetRandom(1,500);
	  nGain3=int(xGain3);
	}
	else{
	  nGain3=0;
	}
	Gain_total+=nGain3;
	for(int k=0; k<nGain3; k++){
	  x3=ran.Gaus(x2+mu_x2,sigma_x2);
	  hx->Fill(x3);
	}
      }
    }
    end=clock();
    tx += end-start;
    //--------------------------------------
    //-- z sampling -----------------------
    start=clock();
    mu_zd=U_zd->Eval(d);
    sigma_zd=Sigma_zd->Eval(d);
    Gain_total=0;
    m_rand=ran.Uniform(0.0,1.0);
    if(m_rand<trans1){
      xGain1=polya[0]->GetRandom(1,500);
      nGain1=int(xGain1);
    }
    else{
      nGain1=0;
    }
    for(int i=0; i<nGain1; i++){
      z1=ran.Gaus(mu_zd,sigma_zd);
      m_rand=ran.Uniform(0.0,1.0);
      if(m_rand<trans2){
	xGain2=polya[1]->GetRandom(1,500);
	nGain2=int(xGain2);
      }
      else{
	nGain2=0;
      }
      for(int j=0; j<nGain2; j++){
	z2=ran.Gaus(z1+mu_z1,sigma_z1);
	m_rand=ran.Uniform(0.0,1.0);
	if(m_rand<trans3){
	  xGain3=polya[2]->GetRandom(1,500);
	  nGain3=int(xGain3);
	}
	else{
	  nGain3=0;
	}
	Gain_total+=nGain3;
	for(int k=0; k<nGain3; k++){
	  z3=ran.Gaus(z2+mu_z2,sigma_z2);
	  hz->Fill(z3);
	}
      }
    }
    end=clock();
    tz += end-start;
    //---------------------------------------
    //--- t sampling -----------------------
    start=clock();
    mu_td=U_td->Eval(d);
    sigma_td=Sigma_td->Eval(d);
    Gain_total=0;
    m_rand=ran.Uniform(0.0,1.0);
    if(m_rand<trans1){
      xGain1=polya[0]->GetRandom(1,500);
      nGain1=int(xGain1);
    }
    else{
      nGain1=0;
    }
    for(int i=0; i<nGain1; i++){
      t1=ran.Gaus(mu_td,sigma_td);
      ht_ex->Fill(t1);
      m_rand=ran.Uniform(0.0,1.0);
      if(m_rand<trans2){
	xGain2=polya[1]->GetRandom(1,500);
	nGain2=int(xGain2);
      }
      else{
	nGain2=0;
      }
      for(int j=0; j<nGain2; j++){
	t2=ran.Gaus(t1+mu_t1,sigma_t1);
	m_rand=ran.Uniform(0.0,1.0);
	if(m_rand<trans3){
	  xGain3=polya[2]->GetRandom(1,500);
	  nGain3=int(xGain3);
	}
	else{
	  nGain3=0;
	}
	Gain_total+=nGain3;
	for(int k=0; k<nGain3; k++){
	  t3=ran.Gaus(t2+mu_t2,sigma_t2);
	  ht->Fill(t3);
	}
      }
    }
    end=clock();
    tt += end-start;
    
    //---------------------------------------------------
    
    //-- simplified sampling ---------------------------
    
    //-- x sampling --------------------------
    start=clock();
    id=int((d-d_min)/delta_d);
    // cout<<id<<endl;
    Gain_total=0;
    xGain1=h_gain1->GetRandom();
    nGain1=int(xGain1);
    hx1_temp1->FillRandom(h_gaus_x[id],nGain1);
    double* arr_x=hx1_temp1->GetArray();
    for(int j=0; j<npos1; j++){
      if(arr_x[j+1]>0){
	nGain23=0;
	for(int k=0; k<arr_x[j+1]; k++){
	  nGain23+=int(h_gain23->GetRandom());
	}
	hx1_temp2->Add(hx23_arr[j],nGain23);
	Gain_total += nGain23;
      }
    }
    //cout<<Gain_total<<endl;
    integ=hx1_temp2->Integral();
    if(integ!=0) hx1->FillRandom(hx1_temp2,Gain_total);     
    hx1_temp1->Reset();
    end=clock();
    tx1 += end-start;
    //-- z sampling --------------------------
    start=clock();
    id=int((d-d_min)/delta_d);
    Gain_total=0;
    xGain1=h_gain1->GetRandom();
    nGain1=int(xGain1);
    hz1_temp1->FillRandom(h_gaus_z[id],nGain1);
    double* arr_z=hz1_temp1->GetArray();
    for(int j=0; j<npos1; j++){
      if(arr_z[j+1]>0){
	nGain23=0;
	for(int k=0; k<arr_z[j+1]; k++){
	  nGain23+=int(h_gain23->GetRandom());
	}
	hz1_temp2->Add(hz23_arr[j],nGain23);
	Gain_total += nGain23;
      }
    }
    integ=hz1_temp2->Integral();
    if(integ!=0) hz1->FillRandom(hz1_temp2,Gain_total);
    hz1_temp1->Reset();
    end=clock();
    tz1 += end-start;
    //-- t sampling --------------------------
    start=clock();
    id=int((d-d_min)/delta_d);
    Gain_total=0;
    xGain1=h_gain1->GetRandom();
    nGain1=int(xGain1);
    ht1_temp1->FillRandom(h_gaus_t[id],nGain1);
    double* arr_t=ht1_temp1->GetArray();
    for(int j=0; j<npos1; j++){
      if(arr_t[j+1]>0){
	nGain23=0;
	for(int k=0; k<arr_t[j+1]; k++){
	  nGain23+=int(h_gain23->GetRandom());
	}
	ht1_temp2->Add(ht23_arr[j],nGain23);
	Gain_total += nGain23;
      }
    }
    integ=ht1_temp2->Integral();
    if(integ!=0) ht1->FillRandom(ht1_temp2,Gain_total);
    ht1_temp1->Reset();
    end=clock();
    tt1+=end-start;
    /*
    // directly sampling
    Gain_total=int(h_gain_total->GetRandom());
    for(int i=0;i<Gain_total;i++){
      x=ran.Gaus(mu_xd12,sigma_xd12);
      hx2->Fill(x);
    }

    Gain_total=int(h_gain_total->GetRandom());
    for(int i=0;i<Gain_total;i++){
      z=ran.Gaus(mu_zd12,sigma_zd12);
      hz2->Fill(z);
    }
    Gain_total=int(h_gain_total->GetRandom());
    for(int i=0;i<Gain_total;i++){
      t=ran.Gaus(mu_td12,sigma_td12);
      ht2->Fill(t);
    */
    }
    
    hx1_temp2->Reset();
    hz1_temp2->Reset();
    ht1_temp2->Reset();
    
    // cout<<"iEvt="<<iEvt<<endl;
    /*
    meanX->Fill(hx->GetMean());
    meanZ->Fill(hz->GetMean());
    meanT->Fill(ht->GetMean());
    RMSX->Fill(hx->GetRMS());
    RMSZ->Fill(hz->GetRMS());
    RMST->Fill(ht->GetRMS());
    nMultiE->Fill(hx->GetEntries());
    
    meanX1->Fill(hx1->GetMean());
    meanZ1->Fill(hz1->GetMean());
    meanT1->Fill(ht1->GetMean());
    RMSX1->Fill(hx1->GetRMS());
    RMSZ1->Fill(hz1->GetRMS());
    RMST1->Fill(ht1->GetRMS());
    nMultiE1->Fill(hx1->GetEntries());
    
    meanX2->Fill(hx2->GetMean());
    meanZ2->Fill(hz2->GetMean());
    meanT2->Fill(ht2->GetMean());
    RMSX2->Fill(hx2->GetRMS());
    RMSZ2->Fill(hz2->GetRMS());
    RMST2->Fill(ht2->GetRMS());
    */
    hx->Reset();
    hz->Reset();
    ht->Reset();
    hx1->Reset();
    hz1->Reset();
    ht1->Reset();
    hx2->Reset();
    hz2->Reset();
    ht2->Reset();
    
  }
  
  ofile<<1.0*tx1/tx<<" "1.0*tz1/tz<<" "<<1.0*tt1/tt<<endl;
  }
  double scale;
  TLegend* lg;
  /*
  TCanvas c13;
  lg = new TLegend(0.6,0.7,0.9,0.9);
  nMultiE->SetLineColor(2);
  nMultiE->SetTitle("#it{Distribution of number of multiplied electrons for 8000 events}");
  nMultiE->SetStats(0);
  nMultiE->GetXaxis()->SetTitle("#it{number of multiplied electrons}");
  nMultiE->Draw();
  nMultiE1->Draw("same");
  lg->AddEntry("nMultiE","3-level sampling");
  lg->AddEntry("nMultiE1","Simplified sampling");
  lg->Draw();
  c13.Print("ne_2vs3.png");

  TCanvas c14;
  lg = new TLegend(0.6,0.7,0.9,0.9);
  h_gain23->SetTitle("#it{PDF of gain for GEM2-GEM3 multiplication}");
  h_gain23->SetStats(0);
  scale=1/h_gain23->Integral();
  h_gain23->Scale(scale,"nosw2");
  c14.SetLogy();
  h_gain23->GetXaxis()->SetTitle("#it{gain}");
  h_gain23->GetYaxis()->SetTitle("#it{probability}");
  h_gain23->Draw();
  c14.Print("gain23.png");
  
  TCanvas c15;
  lg = new TLegend(0.6,0.7,0.9,0.9);
  h_gain_total->SetTitle("#it{PDF of gain for GEM1-GEM2-GEM3 multiplication}");
  h_gain_total->SetStats(0);
  scale=1/h_gain_total->Integral();
  h_gain_total->Scale(scale,"nosw2");
  h_gain_total->GetXaxis()->SetTitle("#it{gain}");
  h_gain_total->GetYaxis()->SetTitle("#it{probability}");
  h_gain_total->Draw();
  c15.Print("gain_total.png");
  
  TCanvas c16;
  lg = new TLegend(0.5,0.7,0.9,0.9);
  h_gain1->SetTitle("#it{PDF of gain for GEM1, GEM2, GEM3 multiplication}");
  h_gain1->SetStats(0);
  h_gain1->GetXaxis()->SetTitle("#it{gain}");
  h_gain1->GetYaxis()->SetTitle("#it{probability}");
  scale=1/h_gain1->Integral();
  h_gain1->Scale(scale,"nosw2");
  scale=1/h_gain2->Integral();
  h_gain2->Scale(scale,"nosw2");
  scale=1/h_gain3->Integral();
  h_gain3->Scale(scale,"nosw2");
  h_gain1->SetLineColor(1);
  h_gain2->SetLineColor(2);
  h_gain3->SetLineColor(4);
  c16.SetLogy();
  h_gain1->Draw();
  h_gain2->Draw("same");
  h_gain3->Draw("same");
  lg->AddEntry("h_gain1","GEM1 multiplication");
  lg->AddEntry("h_gain2","GEM1 multiplication");
  lg->AddEntry("h_gain3","GEM1 multiplication");
  lg->Draw();
  c16.Print("gain_each.png");
  */
  
  /*  
  TCanvas c1;
  lg = new TLegend(0.6,0.7,0.9,0.9);
  meanX->SetLineColor(2);
  meanX2->SetTitle("#it{mean of #deltaX distribution for 8000 events}");
  meanX2->SetStats(0);
  meanX2->GetXaxis()->SetTitle("#it{mean (mm)}");
  meanX2->Draw();
  meanX->Draw("same");
  lg->AddEntry("meanX","3-level sampling");
  lg->AddEntry("meanX2","Direct sampling");
  lg->Draw();
  c1.Print("mean_xdvs3.png");

  TCanvas c2;
  lg = new TLegend(0.6,0.7,0.9,0.9);
  meanZ->SetLineColor(2);
  meanZ2->SetTitle("#it{mean of #deltaZ distribution for 8000 events}");
  meanZ2->SetStats(0);
  meanZ2->GetXaxis()->SetTitle("#it{mean (mm)}");
  meanZ2->Draw();
  meanZ->Draw("same");
  lg->AddEntry("meanZ","3-level sampling");
  lg->AddEntry("meanZ2","Direct sampling");
  lg->Draw();
  c2.Print("mean_zdvs3.png");

  TCanvas c3;
  lg = new TLegend(0.3,0.7,0.6,0.9);
  meanT->SetLineColor(2);
  meanT2->SetTitle("#it{mean of #deltaT distribution for 8000 events}");
  meanT2->SetStats(0);
  meanT2->GetXaxis()->SetTitle("#it{mean (ns)}");
  meanT2->Draw();
  meanT->Draw("same");
  lg->AddEntry("meanT","3-level sampling");
  lg->AddEntry("meanT2","Direct sampling");
  lg->Draw();
  c3.Print("mean_tdvs3.png");

  TCanvas c4;
  lg = new TLegend(0.6,0.7,0.9,0.9);
  RMSX->SetLineColor(2);
  RMSX2->SetTitle("#it{RMS of #deltaX distribution for 8000 events}");
  RMSX2->SetStats(0);
  RMSX2->GetXaxis()->SetTitle("#it{RMS (mm)}");
  RMSX2->Draw();
  RMSX->Draw("same");
  lg->AddEntry("RMSX","3-level sampling");
  lg->AddEntry("RMSX2","Direct sampling");
  lg->Draw();
  c4.Print("RMS_xdvs3.png");

  TCanvas c5;
  lg = new TLegend(0.6,0.7,0.9,0.9);
  RMSZ->SetLineColor(2);
  RMSZ2->SetTitle("#it{RMS of #deltaZ distribution for 8000 events}");
  RMSZ2->SetStats(0);
  RMSZ2->GetXaxis()->SetTitle("#it{RMS (mm)}");
  RMSZ2->Draw();
  RMSZ->Draw("same");
  lg->AddEntry("RMSZ","3-level sampling");
  lg->AddEntry("RMSZ2","Direct sampling");
  lg->Draw();
  c5.Print("RMS_zdvs3.png");

  TCanvas c6;
  lg = new TLegend(0.6,0.7,0.9,0.9);
  RMST->SetLineColor(2);
  RMST2->SetTitle("#it{RMS of #deltaT distribution for 8000 events}");
  RMST2->SetStats(0);
  RMST2->GetXaxis()->SetTitle("#it{RMS (ns)}");
  RMST2->Draw();
  RMST->Draw("same");
  lg->AddEntry("RMST","3-level sampling");
  lg->AddEntry("RMST2","Direct sampling");
  lg->Draw();
  c6.Print("RMS_tdvs3.png");

  TCanvas c7;
  lg = new TLegend(0.6,0.7,0.9,0.9);
  meanX->SetLineColor(2);
  meanX1->SetTitle("#it{mean of #deltaX distribution for 8000 events}");
  meanX1->SetStats(0);
  meanX1->GetXaxis()->SetTitle("#it{mean (mm)}");
  meanX1->Draw();
  meanX->Draw("same");
  lg->AddEntry("meanX","3-level sampling");
  lg->AddEntry("meanX1","Simplified sampling");
  lg->Draw();
  c7.Print("mean_x2vs3.png");

  TCanvas c8;
  lg = new TLegend(0.6,0.7,0.9,0.9);
  meanZ->SetLineColor(2);
  meanZ1->SetTitle("#it{mean of #deltaZ distribution for 8000 events}");
  meanZ1->SetStats(0);
  meanZ1->GetXaxis()->SetTitle("#it{mean (mm)}");
  meanZ1->Draw();
  meanZ->Draw("same");
  lg->AddEntry("meanZ","3-level sampling");
  lg->AddEntry("meanZ1","Simplified sampling");
  lg->Draw();
  c8.Print("mean_z2vs3.png");

  TCanvas c9;
  lg = new TLegend(0.3,0.7,0.6,0.9);
  meanT->SetLineColor(2);
  meanT1->SetTitle("#it{mean of #deltaT distribution for 8000 events}");
  meanT1->SetStats(0);
  meanT1->GetXaxis()->SetTitle("#it{mean (mm)}");
  meanT1->Draw();
  meanT->Draw("same");
  lg->AddEntry("meanT","3-level sampling");
  lg->AddEntry("meanT1","Simplified sampling");
  lg->Draw();
  c9.Print("mean_t2vs3.png");

  TCanvas c10;
  lg = new TLegend(0.6,0.7,0.9,0.9);
  RMSX->SetLineColor(2);
  RMSX1->SetTitle("#it{RMS of #deltaX distribution for 8000 events}");
  RMSX1->SetStats(0);
  RMSX1->GetXaxis()->SetTitle("#it{RMS (mm)}");
  RMSX1->Draw();
  RMSX->Draw("same");
  lg->AddEntry("RMSX","3-level sampling");
  lg->AddEntry("RMSX1","Simplified sampling");
  lg->Draw();
  c10.Print("RMS_x2vs3.png");

  TCanvas c11;
  lg = new TLegend(0.6,0.7,0.9,0.9);
  RMSZ->SetLineColor(2);
  RMSZ1->SetTitle("#it{RMS of #deltaZ distribution for 8000 events}");
  RMSZ1->SetStats(0);
  RMSZ1->GetXaxis()->SetTitle("#it{RMS (mm)}");
  RMSZ1->Draw();
  RMSZ->Draw("same");
  lg->AddEntry("RMSZ","3-level sampling");
  lg->AddEntry("RMSZ1","Simplified sampling");
  lg->Draw();
  c11.Print("RMS_z2vs3.png");

  TCanvas c12;
  lg = new TLegend(0.6,0.7,0.9,0.9);
  RMST->SetLineColor(2);
  RMST1->SetTitle("#it{RMS of #deltaT distribution for 8000 events}");
  RMST1->SetStats(0);
  RMST1->GetXaxis()->SetTitle("#it{RMS (ns)}");
  RMST1->Draw();
  RMST->Draw("same");
  lg->AddEntry("RMST","3-level sampling");
  lg->AddEntry("RMST1","Simplified sampling");
  lg->Draw();
  c12.Print("RMS_t2vs3.png");
  */
  
  /*
  TCanvas c4;
  lg = new TLegend(0.6,0.7,0.9,0.9);
  hx->SetTitle("#it{#deltaX distribution comparison (non-scaled)}");
  hx->GetXaxis()->SetTitle("#it{#deltaX (mm)}");
  hx->SetStats(0);
  hx->SetLineColor(2);
  // scale=hx->Integral()/hx2->Integral();
  // hx2->Scale(scale,"nosw2");
  hx->Draw();
  hx2->Draw("same");
  // gPad->Legend();
  lg->AddEntry("hx","3 level sampling");
  lg->AddEntry("hx2","Direct sampling");
  lg->Draw();
  c4.Print("x_dvs3ns.png");
  
  TCanvas c5;
  lg = new TLegend(0.6,0.7,0.9,0.9);
  hz->SetTitle("#it{#deltaZ distribution comparison (non-scaled)}");
  hz->GetXaxis()->SetTitle("#it{#deltaZ (mm)}");
  hz->SetStats(0);
  hz->SetLineColor(2);
  //  scale=hz->Integral()/hz2->Integral();
  //  hz2->Scale(scale,"nosw2");
  hz->Draw();
  hz2->Draw("same");
  lg->AddEntry("hz","3 level sampling");
  lg->AddEntry("hz2","Direct sampling");
  lg->Draw();
  c5.Print("z_dvs3ns.png");
  
  TCanvas c6;
  lg = new TLegend(0.6,0.7,0.9,0.9);
  ht->SetTitle("#it{#deltaT distribution comparison (non-scaled)}");
  ht->GetXaxis()->SetTitle("#it{#deltaT (ns)}");
  ht->SetStats(0);
  ht->SetLineColor(2);
  //  scale=ht->Integral()/ht2->Integral();
  //  ht2->Scale(scale,"nosw2");
  ht->Draw();
  ht2->Draw("same");
  lg->AddEntry("ht","3 level sampling");
  lg->AddEntry("ht2","Direct sampling");
  lg->Draw();
  c6.Print("t_dvs3ns.png");
  
  TCanvas c1;
  lg = new TLegend(0.6,0.7,0.9,0.9);
  hx->SetTitle("#it{#deltaX distribution comparison (scaled)}");
  hx->GetXaxis()->SetTitle("#it{#deltaX (mm)}");
  hx->SetStats(0);
  hx->SetLineColor(2);
  scale=hx->Integral()/hx2->Integral();
  hx2->Scale(scale,"nosw2");
  hx->Draw();
  hx2->Draw("same");
  // gPad->Legend();
  lg->AddEntry("hx","3 level sampling");
  lg->AddEntry("hx2","Direct sampling");
  lg->Draw();
  c1.Print("x_dvs3s.png");
  
  TCanvas c2;
  lg = new TLegend(0.6,0.7,0.9,0.9);
  hz->SetTitle("#it{#deltaZ distribution comparison (scaled)}");
  hz->GetXaxis()->SetTitle("#it{#deltaZ (mm)}");
  hz->SetStats(0);
  hz->SetLineColor(2);
  scale=hz->Integral()/hz2->Integral();
  hz2->Scale(scale,"nosw2");
  hz->Draw();
  hz2->Draw("same");
  lg->AddEntry("hz","3 level sampling");
  lg->AddEntry("hz2","Direct sampling");
  lg->Draw();
  c2.Print("z_dvs3s.png");
  
  TCanvas c3;
  lg = new TLegend(0.6,0.7,0.9,0.9);
  ht->SetTitle("#it{#deltaT distribution comparison (scaled)}");
  ht->GetXaxis()->SetTitle("#it{#deltaT (ns)}");
  ht->SetStats(0);
  ht->SetLineColor(2);
  scale=ht->Integral()/ht2->Integral();
  ht2->Scale(scale,"nosw2");
  ht->Draw();
  ht2->Draw("same");
  lg->AddEntry("ht","3 level sampling");
  lg->AddEntry("ht2","Direct sampling");
  lg->Draw();
  c3.Print("t_dvs3s.png");
  
  TCanvas c10;
  lg = new TLegend(0.6,0.7,0.9,0.9);
  hx->SetTitle("#it{#deltaX distribution comparison (non-scaled)}");
  hx->GetXaxis()->SetTitle("#it{#deltaX (mm)}");
  hx->SetStats(0);
  hx->SetLineColor(2);
  //  scale=hx->Integral()/hx1->Integral();
  //  hx1->Scale(scale,"nosw2");
  hx->Draw();
  hx1->Draw("same");
  // gPad->Legend();
  lg->AddEntry("hx","3-level sampling");
  lg->AddEntry("hx1","Simplified sampling");
  lg->Draw();
  c10.Print("x_2vs3ns.png");
  
  TCanvas c11;
  lg = new TLegend(0.6,0.7,0.9,0.9);
  hz->SetTitle("#it{#deltaZ distribution comparison (non-scaled)}");
  hz->GetXaxis()->SetTitle("#it{#deltaZ (mm)}");
  hz->SetStats(0);
  hz->SetLineColor(2);
  //  scale=hz->Integral()/hz1->Integral();
  //  hz1->Scale(scale,"nosw2");
  hz->Draw();
  hz1->Draw("same");
  lg->AddEntry("hz","3-level sampling");
  lg->AddEntry("hz1","Simplified sampling");
  lg->Draw();
  c11.Print("z_2vs3ns.png");
  
  TCanvas c12;
  lg = new TLegend(0.6,0.7,0.9,0.9);
  ht->SetTitle("#it{#deltaT distribution comparison (non-scaled)}");
  ht->GetXaxis()->SetTitle("#it{#deltaT (ns)}");
  ht->SetStats(0);
  ht->SetLineColor(2);
  //  scale=ht->Integral()/ht1->Integral();
  //  ht1->Scale(scale,"nosw2");
  ht->Draw();
  ht1->Draw("same");
  lg->AddEntry("ht","3-level sampling");
  lg->AddEntry("ht1","Simplified sampling");
  lg->Draw();
  c12.Print("t_2vs3ns.png");

  
  TCanvas c7;
  lg = new TLegend(0.6,0.7,0.9,0.9);
  hx->SetTitle("#it{#deltaX distribution comparison (scaled)}");
  hx->GetXaxis()->SetTitle("#it{#deltaX (mm)}");
  hx->SetStats(0);
  hx->SetLineColor(2);
  scale=hx->Integral()/hx1->Integral();
  hx1->Scale(scale,"nosw2");
  hx->Draw();
  hx1->Draw("same");
  // gPad->Legend();
  lg->AddEntry("hx","3-level sampling");
  lg->AddEntry("hx1","Simplified sampling");
  lg->Draw();
  c7.Print("x_2vs3s.png");
  
  TCanvas c8;
  lg = new TLegend(0.6,0.7,0.9,0.9);
  hz->SetTitle("#it{#deltaZ distribution comparison (scaled)}");
  hz->GetXaxis()->SetTitle("#it{#deltaZ (mm)}");
  hz->SetStats(0);
  hz->SetLineColor(2);
  scale=hz->Integral()/hz1->Integral();
  hz1->Scale(scale,"nosw2");
  hz->Draw();
  hz1->Draw("same");
  lg->AddEntry("hz","3-level sampling");
  lg->AddEntry("hz1","Simplified sampling");
  lg->Draw();
  c8.Print("z_2vs3s.png");
  
  TCanvas c9;
  lg = new TLegend(0.6,0.7,0.9,0.9);
  ht->SetTitle("#it{#deltaT distribution comparison (scaled)}");
  ht->GetXaxis()->SetTitle("#it{#deltaT (ns)}");
  ht->SetStats(0);
  ht->SetLineColor(2);
  scale=ht->Integral()/ht1->Integral();
  ht1->Scale(scale,"nosw2");
  ht->Draw();
  ht1->Draw("same");
  lg->AddEntry("ht","3-level sampling");
  lg->AddEntry("ht1","Simplified sampling");
  lg->Draw();
  c9.Print("t_2vs3s.png");
  
  */

  
  
  /*
  auto legend = new TLegend(0.3,0.7,0.5,0.8);
  TCanvas tcan;
  tcan.SetLogy();
  h_gain1->SetLineColor(1);
  h_gain1->Draw();
  h_gain2->SetLineColor(2);
  h_gain3->SetLineColor(4);
  h_gain2->Draw("same");
  h_gain3->Draw("same");
  legend->AddEntry("h_gain1","Gain for GEM1");
  legend->AddEntry("h_gain2","Gain for GEM2");
  legend->AddEntry("h_gain3","Gain for GEM3");
  legend->SetTextFont(42);
  legend->Draw();
  */
  /*
  double scale;
  TCanvas c1;
  hx->Draw();
  hx1->Draw("same");
  //  hz->Draw();
  // ht->Draw();
  // h_gaus_x[id]->Draw();

  TCanvas c3;
  hz->Draw();
  hz1->Draw("same");

  TCanvas c4;
  ht->SetLineColor(2);
  ht->Draw();
  ht1->Draw("same");

  TCanvas c6;
  hx_ex->SetLineColor(2);
  scale = hx_ex->Integral()/h_gaus_x[id]->Integral();
  hx_ex->Draw();
  h_gaus_x[id]->Scale(scale,"nosw2");
  h_gaus_x[id]->Draw("same");
  

  TCanvas c5;
  ht_ex->SetLineColor(2);
  scale = ht_ex->Integral()/h_gaus_t[id]->Integral();
  ht_ex->Draw();
  h_gaus_t[id]->Scale(scale,"nosw2");
  h_gaus_t[id]->Draw("same");
  
  TCanvas c2;
  c2.SetLogy();
  h_gain23->Draw();
  */
  /*
  TCanvas c7;
  meanX->SetLineColor(2);
  meanX->Draw();
  meanX1->Draw("same");

  TCanvas c8;
  meanZ->SetLineColor(2);
  meanZ->Draw();
  meanZ1->Draw("same");

  TCanvas c9;
  meanT->SetLineColor(2);
  meanT->Draw();
  meanT1->Draw("same");

  TCanvas c10;
  RMSX->SetLineColor(2);
  RMSX->Draw();
  RMSX1->Draw("same");

  TCanvas c11;
  RMSZ->SetLineColor(2);
  RMSZ->Draw();
  RMSZ1->Draw("same");

  TCanvas c12;
  RMST->SetLineColor(2);
  RMST->Draw();
  RMST1->Draw("same");

  cout<<"time for x: "<<1.0*tx1/tx<<endl;
  cout<<"time for z: "<<1.0*tz1/tz<<endl;
  cout<<"time for t: "<<1.0*tt1/tt<<endl;
  */
}
    

    
  
