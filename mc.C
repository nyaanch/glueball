#include <iostream>
#include <vector>
#include <map>
#include <fstream>
#include <string>
#include <sstream>
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

void de_dx(TCanvas *c1, TH2F* dedx)
{
	gStyle->SetOptTitle(0);
	gStyle->SetOptStat(0);
	
	dedx->GetYaxis()->SetTitle("dE/dX (MeV/cm)");
	dedx->GetXaxis()->SetTitle("track momentum (GeV)");
	
	dedx->Draw("colz");	
}

void px_py_2d(TCanvas *c1, TH2F* px_2d, TH2F* py_2d, TH2F* px_2d_recoProt, TH2F* py_2d_recoProt)
{
	gStyle->SetOptTitle(0);
	gStyle->SetOptStat(0);
	
	px_2d->GetYaxis()->SetTitle("#Sigma p_{x, gen p} (GeV)");
	px_2d->GetXaxis()->SetTitle("#Sigma p_{x, gen #pi} (GeV)");
	
	py_2d->GetYaxis()->SetTitle("#Sigma p_{y, gen p } (GeV)");
	py_2d->GetXaxis()->SetTitle("#Sigma p_{y, gen #pi} (GeV)");
	
	
	px_2d_recoProt->GetYaxis()->SetTitle("#Sigma p_{x, gen p} (GeV)");
	px_2d_recoProt->GetXaxis()->SetTitle("#Sigma p_{x, reco #pi} (GeV)");
	
	py_2d_recoProt->GetYaxis()->SetTitle("#Sigma p_{y, gen p} (GeV)");
	py_2d_recoProt->GetXaxis()->SetTitle("#Sigma p_{y, reco #pi} (GeV)");
	
	px_2d->Draw("colz");

	py_2d->Draw("colz");

	px_2d_recoProt->Draw("colz");
	
	py_2d_recoProt->Draw("colz");
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

void mc()
{

	TCanvas *c1 = new TCanvas("Figure","Figure",1000,800);

	TFile *input = new TFile("MinBias.root", "read");	

	TTree *tree = (TTree*)input->Get("analysis/tree"); //skim
	
	//root file variables
	
	Int_t ntrk, trk_id[1000], trk_q[1000];
	Float_t trk_pt[1000], trk_eta[1000], trk_phi[1000], trk_dedx[1000]; //TRK variables
	
	Int_t  ngentrk, gentrk_id[1000];
	Float_t gentrk_pt[1000], gentrk_m[1000], gentrk_eta[1000], gentrk_phi[1000]; //GENTRK variables
	
        Int_t entries = (Int_t)tree->GetEntries();
	
	//branch addresses
	
	//TREE BRANCH ADDRESS
	tree->SetBranchAddress("trk_pt",trk_pt);
	tree->SetBranchAddress("ntrk", &ntrk);
	tree->SetBranchAddress("trk_phi",trk_phi);
	tree->SetBranchAddress("trk_q",trk_q);
	tree->SetBranchAddress("trk_eta",trk_eta);
	
   	tree->SetBranchAddress("trk_dedx", trk_dedx);

	//GENTRK BRANCH ADDRESS
	tree->SetBranchAddress("ngentrk", &ngentrk);
	tree->SetBranchAddress("gentrk_pt",gentrk_pt);
	tree->SetBranchAddress("gentrk_phi",gentrk_phi);
	tree->SetBranchAddress("gentrk_id",gentrk_id);
	tree->SetBranchAddress("gentrk_m",gentrk_m);
	tree->SetBranchAddress("gentrk_eta",gentrk_eta);
	
	
	//arrays and 'local' variables
	
	Double_t pion_p[2][4], pion_m[2][4];
	//0a - 1 particle, 1a - 2 particle etc
	//0-energy, 1-px, 2-py, 3-pz, 4-charge
	
	Int_t n_p, n_m;
	
	Double_t p;
	
	Double_t px, py, pz, E;
	Double_t px_s[1000], py_s[1000], pz_s[1000];
	Double_t px_prot, py_prot, E_prot;
	Double_t px_pion, py_pion, E_pion;
	Double_t px_reco, py_reco, E_reco;
	Double_t p11[5], p22[5], p12[5], p21[5];
	
	Double_t pion_mass = 0.13957039; //[GeV]

	Int_t n_4trk = 0;
	
	Bool_t pt_low;
	
	//histograms
	
	//two dimensional histograms
	TH2F* m = new TH2F("m", "m", /* X-dimension */ 400, 0, 2, /* Y-dimension */ 400, 0, 2);
	TH2F* pt_cut = new TH2F("pt_cut", "pt_cut", /* X-dimension */ 40, 0, 2, /* Y-dimension */ 40, 0, 2);
	
	TH2F* px_2d = new TH2F("px_2d", "px_2d", /* X-dimension */ 100, -2, 2, /* Y-dimension */ 100, -2, 2);
	TH2F* py_2d = new TH2F("py_2d", "py_2d", /* X-dimension */ 100, -2, 2, /* Y-dimension */ 100, -2, 2);
	
	TH2F* px_2d_recoProt = new TH2F("px_2d_recoProt", "px_2d_recoProt", /* X-dimension */ 100, -2, 2, /* Y-dimension */ 100, -2, 2);
	TH2F* py_2d_recoProt = new TH2F("py_2d_recoProt", "py_2d_recoProt", /* X-dimension */ 100, -2, 2, /* Y-dimension */ 100, -2, 2);
	
	TH2F* dedx = new TH2F("dedx", "dedx", 300, 0, 2, 300, 0, 20);
	
	cout<<"entries = "<<entries<<endl;
	
	Bool_t mass_projection, mass_ptrans_cut, pt_projection, twoD, dedx_, proton_angle, th;
	
	//for functions
	mass_projection = false;
	mass_ptrans_cut = false;
	pt_projection = false;
	twoD = false;
	dedx_ = false;
		
	//separate fit or general (for mass projection)
	
	Bool_t separate = false;

	for(Int_t e=0; e<entries; ++e)   
	{
		tree->GetEntry(e);   
		
		n_p=0;
		n_m=0;
					
		px_pion = 0;
		py_pion = 0;
		
		px_prot = 0;
		py_prot = 0;
				
		for(Int_t i=0; i < ngentrk; ++i)
		{	
			if(gentrk_id[i] == 2212) 
			{	
				px_prot += gentrk_pt[i]*cos(gentrk_phi[i]);
				py_prot += gentrk_pt[i]*sin(gentrk_phi[i]);
				
			}
			else
			{
				px_pion += gentrk_pt[i]*cos(gentrk_phi[i]);
				py_pion += gentrk_pt[i]*sin(gentrk_phi[i]);
			}
		}
		
		
		px_2d->Fill(px_prot, px_pion);
		py_2d->Fill(py_prot, py_pion); 
		
		px_reco = 0;
		py_reco = 0;

		if(ntrk == 4) //selecting 4-tracks
		{
			pt_low = false;
			
			n_4trk++;
			
			for(Int_t j=0; j < ntrk; ++j)	
			{	
				if(trk_pt[j] < 0.25) pt_low = true;
				
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
				
					if((trk_q[j] == 1) && (n_p<2))
					{
						pion_p[n_p][0] = E;
						pion_p[n_p][1] = px;
						pion_p[n_p][2] = py;
						pion_p[n_p][3] = pz;
						n_p++;
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

			px_2d_recoProt->Fill(px_prot, px_reco);
			py_2d_recoProt->Fill(py_prot, py_reco);
		}
				
	}
	
	cout<<"n_4trk = "<<n_4trk<<endl;
	
	if(mass_projection) m_pro(c1, m, separate); //mass projection
	
	if(mass_ptrans_cut) m_pt_cut(c1, m, pt_cut); //2d invariant mass and pt histogram
	
	if(pt_projection) pt_pro(c1, pt_cut); //pt projection
	
	if(twoD) px_py_2d(c1, px_2d, py_2d, px_2d_recoProt, py_2d_recoProt);
	
	if(dedx_) de_dx(c1, dedx);
}

int main()
{
	mc();
	return 0;
}
