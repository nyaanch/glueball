#include <iostream>
#include <vector>
#include <map>
#include <fstream>
#include <string>
#include <sstream>
#include <iomanip>
#include "math.h"
#include "math.h"

using namespace std;
#include "stdlib.h"

#include "TMath.h"
#include "TTree.h"
#include "TF1.h"
#include "TH1F.h"
#include "TGraph.h"
#include "TFile.h"
#include "TH1D.h"
#include "TLegend.h"
#include "TCanvas.h"
#include "TAxis.h"
#include "TGraphErrors.h"
#include "TLatex.h"
#include "TLine.h"
#include "TMinuit.h"
#include "TStyle.h"
#include "TProfile.h"
#include "TVector.h"
#include "TMatrix.h"
#include "TMarker.h"
#include "TRandom3.h"
#include "TMultiGraph.h"
TRandom3 random_generator ;

void theta(TCanvas *c1, TH2F* th_px, TH2F* th_py, TH2F* th_px1, TH2F* th_py1)
{
	th_px->Draw("colz");
	th_py->Draw("colz");
	
	th_px1->Draw("colz");
	//th_py1->Draw("colz");
}


void de_dx(TCanvas *c1, TH2F* dedx)
{
	gStyle->SetOptTitle(0);
	gStyle->SetOptStat(0);
	
	dedx->GetYaxis()->SetTitle("dE/dX (MeV/cm)");
	dedx->GetXaxis()->SetTitle("track momentum (GeV)");
	
	dedx->Draw("colz");	
}

void prot_an(TCanvas *c1, TH2F* prot_an_x, TH2F* prot_an_y)
{
	gStyle->SetOptTitle(0);
	gStyle->SetOptStat(0);
	
	prot_an_x->GetYaxis()->SetTitle("#theta_{x}L");
	prot_an_x->GetXaxis()->SetTitle("#theta_{x}R");
	
	prot_an_x->Draw("colz");

	prot_an_y->GetYaxis()->SetTitle("M1 (GeV)");
	prot_an_y->GetXaxis()->SetTitle("M2 (GeV)");
	
	prot_an_y->Draw("colz");
}

void pt_pro(TCanvas *c1, TH2F* pt_cut)
{	
	gStyle->SetOptTitle(0);
	
	TH1D* pt_proX = pt_cut->ProjectionX();
	TH1D* pt_proY = pt_cut->ProjectionY();

	pt_proY->GetYaxis()->SetTitle("M1 (GeV)");
	pt_proY->GetXaxis()->SetTitle("M2 (GeV)");

	pt_proX->Add(pt_proY);

	pt_proX->Draw();
}

void m_pro(TCanvas *c1, TH2F* m, Bool_t separate)
{
	gStyle->SetOptTitle(0);
	gStyle->SetLineScalePS(.3);

	TH1D* m_proX = m->ProjectionX();
	TH1D* m_proY = m->ProjectionY();

	m_proX->GetYaxis()->SetTitle("Events");
	m_proX->GetXaxis()->SetTitle("M (GeV)");
	
	m_proY->GetYaxis()->SetTitle("Events");
	m_proY->GetXaxis()->SetTitle("M (GeV)");
	
	m_proX->Sumw2();
	m_proX->SetMarkerStyle(4);

	m_proX->Add(m_proY);
	
	if(!separate)
	{
		///FIT FUNCTION
	
		TF1* fit = new TF1("fit", "[0]*exp(-((x-[1])*(x-[1]))/(2*[2]*[2])) + [3]*exp(-((x-[4])*(x-[4]))/(2*[5]*[5])) + [6]*pow((x-[7]),[8])*exp(([9]*x)+([10]*x*x)+([11]*x*x*x)) + [12]*exp(-((x-[13])*(x-[13]))/(2*[14]*[14]))", 0, 2);
	
		fit->SetParameters(1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1); //11
		fit->SetParameter(11, 1);
		fit->SetParameter(12, 10000);
		fit->FixParameter(13, 1.3);
		fit->SetParameter(14, 0.4);

		fit->SetParameters(1.92941e+03, 7.43308e-01, -5.43028e-02, 1.38772e+03, 4.90000e-01, 1.85130e-02, 5.26645e+01 ,   -1.05882e+00,    2.97467e+01 ,  -1.78941e+01);

	/*
		fit->FixParameter(0, 1.92941e+03);
		fit->FixParameter(1, 7.43308e-01);           //rho
		fit->FixParameter(2, -5.43028e-02);
	*/	
	
		fit->FixParameter(3, 2.04548e+03); 
		fit->FixParameter(4, 4.97263e-01);           //K_0
		fit->FixParameter(5, 1.08482e-02);
	/*
		fit->FixParameter(6, -4.04679e+02);
		fit->FixParameter(7, 2.09713e+00);
		fit->FixParameter(8, 1.00000e+00);            //background
		fit->FixParameter(9, 3.89325e-01); 
	*/	
		Double_t left, right;
	
		left = 0.33;
		right = 1.475;
	
		m_proX->Fit("fit", "", "", left, right);
	
		Int_t bin1, bin2;
	
		bin1 = m_proX->FindBin(left);
		bin2 = m_proX->FindBin(right);
	
		cout<<left<<"->"<<bin1<<endl;
		cout<<right<<"->"<<bin2<<endl;

		fit->SetLineColor(kRed);
	
		m_proX->Draw("");
	
	}
	else
	{
		//separate fit functions
		TF1* rho = new TF1("kaon", "gaus", 0.5, 1);
		TF1* kaon = new TF1("pion", "gaus", 0.3, 0.8);
	
		rho->SetParameter(0, 1.92941e+03);
		rho->SetParameter(1, 7.43308e-01);
		rho->SetParameter(2, -5.43028e-02);
	
		kaon->SetParameter(0, 2.04548e+03);
		kaon->SetParameter(1, 4.97263e-01);
		kaon->SetParameter(2, 1.08482e-02);
	
		kaon->SetLineColor(kBlue);
		rho->SetLineColor(kMagenta);
	
		m_proX->Draw("");
		rho->Draw("SAME");
		kaon->Draw("SAME");
	}
}

void m_pt_cut(TCanvas *c1, TH2F* m, TH2F* pt_cut)
{
	gStyle->SetOptTitle(0);
	gStyle->SetLineScalePS(.3);	

	m->GetYaxis()->SetTitle("M1 (GeV)");
	m->GetXaxis()->SetTitle("M2 (GeV)");
	
	m->Draw("colz");

	pt_cut->GetYaxis()->SetTitle("M1 (GeV)");
	pt_cut->GetXaxis()->SetTitle("M2 (GeV)");
	
	pt_cut->Draw("colz");
}

void data()
{
	TCanvas *c1 = new TCanvas("Figure","Figure",1000,800);

	TFile *input = new TFile("MinBias.root", "read");	

	TTree *tree = (TTree*)input->Get("tree;1"); //skim
	
	//root file variables
	
	Int_t ntrk, trk_id[1000], trk_q[1000];
	Float_t trk_pt[1000], trk_eta[1000], trk_phi[1000], trk_dedx[1000], ThxR, ThxL, ThyR, ThyL; //TRK variables

        Int_t entries = (Int_t)tree->GetEntries();

	//TREE BRANCH ADDRESS
	tree->SetBranchAddress("trk_pt",trk_pt);
	tree->SetBranchAddress("ntrk", &ntrk);
	tree->SetBranchAddress("trk_phi",trk_phi);
	tree->SetBranchAddress("trk_q",trk_q);
	tree->SetBranchAddress("trk_eta",trk_eta);
   	tree->SetBranchAddress("trk_dedx", trk_dedx);
	tree->SetBranchAddress("ThxR",&ThxR);
	tree->SetBranchAddress("ThxL",&ThxL);
	tree->SetBranchAddress("ThyR",&ThyR);
	tree->SetBranchAddress("ThyL",&ThyL);

	//arrays and 'local' variables
	
	Double_t pion_p[2][4], pion_m[2][4];
	//0a - 1 particle, 1a - 2 particle etc
	//0-energy, 1-px, 2-py, 3-pz, 4-charge
	
	Int_t n_p, n_m;
	
	Double_t p;
	
	Double_t px, py, pz, E;
	Double_t px_s[1000], py_s[1000], pz_s[1000];
	Double_t px_reco, py_reco, E_reco;
	Double_t p11[5], p22[5], p12[5], p21[5];
	Double_t x_angle, y_angle, p_x, p_y;
	
	Double_t pion_mass = 0.13957039; //[GeV]
	
	Int_t n_4trk = 0;
	
	Bool_t pt_low;
	
	//histograms
	
	//two dimensional histograms
	TH2F* m = new TH2F("m", "m", /* X-dimension */ 400, 0, 2, /* Y-dimension */ 400, 0, 2);
	TH2F* pt_cut = new TH2F("pt_cut", "pt_cut", /* X-dimension */ 40, 0, 2, /* Y-dimension */ 40, 0, 2);
	
	TH2F* px_2d = new TH2F("px_2d", "px_2d", /* X-dimension */ 100, -2, 2, /* Y-dimension */ 100, -2, 2);
	TH2F* py_2d = new TH2F("py_2d", "py_2d", /* X-dimension */ 100, -2, 2, /* Y-dimension */ 100, -2, 2);
	
	TH2F* prot_an_x = new TH2F("prot_an_x", "prot_an_x", 100, -2e-4, 2e-4, 100, -2e-4, 2e-4);
	TH2F* prot_an_y = new TH2F("prot_an_y", "prot_an_y", 100, -2e-4, 2e-4, 100, -2e-4, 2e-4);
	
	TH2F* dedx = new TH2F("dedx", "dedx", 300, 0, 2, 300, 0, 20);
	
	TH2F* th_px = new TH2F("th_px", "th_px", 300, -5, 5, 300, -2e-3, 2e-3);
	TH2F* th_py = new TH2F("th_py", "th_py", 300, -5, 5, 300, -2e-3, 2e-3);
	
	TH2F* th_px1 = new TH2F("th_px1", "th_px1", 300, -5, 5, 300, -5, 5);
	TH2F* th_py1 = new TH2F("th_py1", "th_py1", 300, -5, 5, 300, -2e-3, 2e-3);
	
	Bool_t mass_projection, mass_ptrans_cut, pt_projection, twoD, dedx_, proton_angle, th;
	
	//for functions
	mass_projection = false;
	mass_ptrans_cut = false;
	pt_projection = false;
	twoD = false;
	dedx_ = false;
	proton_angle = false;
	th = false;
	
	//separate fit or general (for mass projection)
	
	Bool_t separate = false;
	
	cout<<"entries = "<<entries<<endl;

	for(Int_t e=0; e<entries; ++e)   
	{
		tree->GetEntry(e);
		   
		n_p=0;
		n_m=0;
					
		px_reco = 0;
		py_reco = 0;
		
		x_angle = 0;
		y_angle = 0;
			
		if(ntrk == 4) //selecting 4-tracks
		{
			pt_low = false;
			
			n_4trk++; //4-track counter
			
			for(Int_t j=0; j < ntrk; ++j)	
			{	
				if(trk_pt[j] < 0.25) pt_low = true;
				
				prot_an_x->Fill(ThxL, ThxR);
				prot_an_y->Fill(ThyL, ThyR);
				
				x_angle = ThxL + ThxR;
				y_angle = ThyL + ThyR;
				
				p = trk_pt[j]*cosh(trk_eta[j]);
				
				dedx->Fill(p, trk_dedx[j]);

				px = trk_pt[j]*cos(trk_phi[j]);
				py = trk_pt[j]*sin(trk_phi[j]);
				pz = trk_pt[j]*sinh(trk_eta[j]);
				
				px_reco += trk_pt[j]*cos(trk_phi[j]);
				py_reco += trk_pt[j]*sin(trk_phi[j]);
			
				E = sqrt(pion_mass*pion_mass + trk_pt[j]*trk_pt[j] + pz*pz);
				
				if(((trk_q[j] == 1) && (n_p<2)) || ((n_m<2) && (trk_q[j] == -1))) //selecting by charge
				{
				
					//cout<<"j =="<<j<<endl;
					
					if((trk_q[j] == 1) && (n_p<2)) 
					{
						pion_p[n_p][0] = E;
						pion_p[n_p][1] = px;
						pion_p[n_p][2] = py;
						pion_p[n_p][3] = pz;
						n_p++;
						
						
						//cout<<"n_p = "<<n_p<<endl;
					}
					if((trk_q[j] == -1) && (n_m<2))
					{
						pion_m[n_m][0] = E;
						pion_m[n_m][1] = px;
						pion_m[n_m][2] = py;
						pion_m[n_m][3] = pz;
						n_m++;
					}	
				}		
			}
			
			for(Int_t a=0; a < 4; ++a)	
			{
				p11[a] = pion_p[0][a] + pion_m[0][a];
				p22[a] = pion_p[1][a] + pion_m[1][a];
				p12[a] = pion_p[0][a] + pion_m[1][a];
				p21[a] = pion_p[1][a] + pion_m[0][a];
			}

			p11[4] =  sqrt(p11[0]*p11[0] - p11[1]*p11[1] - p11[2]*p11[2] - p11[3]*p11[3]);
			p22[4] =  sqrt(p22[0]*p22[0] - p22[1]*p22[1] - p22[2]*p22[2] - p22[3]*p22[3]);
			p12[4] =  sqrt(p12[0]*p12[0] - p12[1]*p12[1] - p12[2]*p12[2] - p12[3]*p12[3]);
			p21[4] =  sqrt(p21[0]*p21[0] - p21[1]*p21[1] - p21[2]*p21[2] - p21[3]*p21[3]);

			m->Fill(p11[4], p22[4]);
			m->Fill(p12[4], p21[4]);

			if(!pt_low) 
			{
				pt_cut->Fill(p11[4], p22[4]);
				pt_cut->Fill(p12[4], p21[4]);
			}
			
			p_x = x_angle*6500;
			p_y = y_angle*6500;
			
			th_px->Fill(px_reco, x_angle);
			th_py->Fill(py_reco, y_angle);
			
			th_px1->Fill(px_reco, p_x);
			th_py1->Fill(py_reco, p_y);
		}
				
	}
	
	cout<<"n_4trk = "<<n_4trk<<endl;

	if(mass_projection) m_pro(c1, m, separate); //mass projection
	
	if(mass_ptrans_cut) m_pt_cut(c1, m, pt_cut); //2d invariant mass and pt histogram
	
	if(pt_projection) pt_pro(c1, pt_cut); //pt projection
	
	if(dedx_) de_dx(c1, dedx);
	
	if(proton_angle) prot_an(c1, prot_an_x, prot_an_y);
	
	if(th) theta(c1, th_px, th_py, th_px1, th_py1);
	
}

int main()
{
	data();
	return 0;
}
