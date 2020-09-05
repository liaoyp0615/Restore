#include <iostream>
#include <fstream>
#include "TFile.h"
#include "TString.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TH1.h"
#include "TRandom.h"
#include "TRandom3.h"
#include "TF1.h"
using namespace std;

int main(){
	//PreDefine--------------------------
	double x1, x2, x3;
	//double mu_x1, sigma_x1, mu_x2, sigma_x2, mu_x12, sigma_x12;
	TRandom3 ran;
	double m_rand;
	TString funname;
	// time variable--------------------
	long start(0), end(0);
	//long tx1(0), tz1(0), tt1(0);
	//long tx2, tz2, tt2;
	long tx(0);
	
	// Define polya function------------------
	TF1* polya[3];
	double par_polya_0[3]={6257, 580.3, 443.2};
	double par_polya_1[3]={27.87,21.03,27.53};
	double par_polya_2[3]={0.4946,0.4088,0.4799};
	for(int i=0;i<3;i++){
		funname.Form("polya%d",i);
		polya[i] = new TF1(funname, "[0]*pow(1+[2],1+[2])*pow(x/[1],[2])*exp(-(1+[2])*x/[1])/TMath::Gamma(1+[2])", 0, 500);
		polya[i]->SetParameter(0,par_polya_0[i]);
		polya[i]->SetParameter(1,par_polya_1[i]);
		polya[i]->SetParameter(2,par_polya_2[i]);
	}
 	
	double trans1(0.867), trans2(0.61), trans3(0.613);
	int Gain_total, xGain1, nGain1, xGain2, nGain2, xGain3, nGain3; 
	int nGain3_count, nGain23;
	
	// Set distribution parameters-------------------------
	double mu_x1=0.3838;
	double sigma_x1 = 0.175;
	double mu_x2 = 0.3774;
	double sigma_x2 = 0.1757;
	double mu_x12 = mu_x1 + mu_x2;  // Define Convolution of two parameters 
	double sigma_x12 = sqrt(sigma_x1*sigma_x1+sigma_x2*sigma_x2);
	
	double mu_xd,sigma_xd; //	???
	TF1* U_xd = new TF1("U_xd","0.02627+0.4785*x",0,5);
	TF1* Sigma_xd = new TF1("Sigma_xd","0.06928+0.09882*x-0.01768*x*x+0.001678*x*x*x",0,5);
	
	double Drift=2.5; // ! set drifting distance
	double mu_x_par[2]={0.02627,0.4785};
	double sigma_x_par[4]={0.06928,0.09882,-0.01768,0.001678};
	double mu_xd12 = mu_x1+mu_x2+U_xd->Eval(Drift);
	double sigma_xd12 = sqrt(sigma_x1*sigma_x1+sigma_x2*sigma_x2+(Sigma_xd->Eval(Drift))*(Sigma_xd->Eval(Drift)));
	
	double x_max = mu_xd12+6*(sigma_x1+sigma_x2+Sigma_xd->Eval(Drift)); //set range of histograms
	double x_min = mu_xd12-6*(sigma_x1+sigma_x2+Sigma_xd->Eval(Drift));
	double x_dmin = U_xd->Eval(Drift) - 6*Sigma_xd->Eval(Drift);
	double x_dmax = U_xd->Eval(Drift) + 6*Sigma_xd->Eval(Drift);
	
	
	// Define histogram for gain
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

	// Sampling number of gaining electrons-------------------------

	TH1D* h_gain23 = new TH1D("h_gain23","",nbin,0,3000);
	
	for(int i=0; i<100000; i++){  //why do it test for 100000 times?   ????
		m_rand = ran.Uniform(0.0,1.0);
		if(m_rand<trans2){
			xGain2=polya[1]->GetRandom();
			nGain2=int(xGain2);
		}
		else{
			nGain2=0;
		}  
		nGain3_count=0;
	
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

	
	
	// Define histogram ----------------------------
	int npos=200;
	int npos1=200;
	TH1D* hx = new TH1D("hx","",npos,x_min,x_max);
	TH1D* hx_ex = new TH1D("hx_ex","",npos1,x_dmin,x_dmax);
	
	
	//3 level sampling--------------------------
    // x sampling ------------------
    start=clock();  //set clock
    
	mu_xd=U_xd->Eval(Drift);  // 让 U_xd 取d并赋值给mu_xd     but,what is d  ???
	sigma_xd=Sigma_xd->Eval(Drift);
    
	Gain_total=0;
	m_rand=ran.Uniform(0.0,1.0);
	if(m_rand<trans1){  //首先，需要考虑透过率 
		xGain1=polya[0]->GetRandom(1,500);
		nGain1=int(xGain1);  //利用polya分布确定第一次倍增电子数 
	}
	else{
		nGain1=0;
	}
	for(int i=0; i<nGain1; i++){  //对第一次倍增的电子进行抽样（只考虑x） 
		x1 = ran.Gaus(mu_xd,sigma_xd);  //对x进行抽样 **** 
		hx_ex->Fill(x1);  // histogram hx_ex 
		m_rand=ran.Uniform(0.0,1.0);  //第二次倍增 
		if(m_rand<trans2){
			xGain2=polya[1]->GetRandom(1,500);
			nGain2=int(xGain2);
		}
		else{
			nGain2=0;
		}
		for(int j=0; j<nGain2; j++){
			x2 = ran.Gaus(x1+mu_x1,sigma_x1); //注意这里，利用第一次抽样后的一个位移x1 
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
				x3 = ran.Gaus(x2+mu_x2,sigma_x2);
				hx->Fill(x3); //histogram hx  最终结果 
			}
		}
	}
	
	// Draw histogram 
	TCanvas c1;
	TLegend* lg;
	lg = new TLegend(0.6,0.7,0.9,0.9);
	hx->SetTitle("#it{#deltaX distribution comparison (scaled)}");
	hx->GetXaxis()->SetTitle("#it{#deltaX (mm)}");
	hx->SetStats(0);
	hx->SetLineColor(2);
	hx->Draw();
	// gPad->Legend();
	lg->AddEntry("hx","Direct sampling");
	lg->Draw();
	c1.Print("x_dSampling.png");
		
    end=clock();
    tx += end-start;
    printf("tx=%d",&tx);
    
}
