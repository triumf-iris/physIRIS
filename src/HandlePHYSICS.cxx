// HandlePhysics.cxx
// calculating Physics variables requiring data from more than one detector

#include <iostream>
#include <assert.h>
#include <fstream>
#include <string>
#include <vector>

#include <TFile.h>
#include <TChain.h>
#include <TMath.h>
#include "TCutG.h"

#include "HandlePHYSICS.h"
#include "eloss.h"
#include "nucleus.h"
#include "runDepPar.h"
#include "CalibPHYSICS.h"
#include "Graphsdedx.h"
#include "geometry.h"

//#define pd // (p,d) reaction analysis
extern TEvent *IrisEvent;
extern TFile* treeFile;
extern TTree* tree;

TChain *input_chain;

CalibPHYSICS calPhys;
geometry geoP;
Graphsdedx dedx_i, dedx_l, dedx_h;
IDet detec; // calibrated variables from detectors, to be passed to HandlePhysics
IDet *det = &detec;
ITdc timeArray;
ITdc *tdc = &timeArray;
IScaler scal;
IScaler *pscaler = &scal;
PTrack lP;
PTrack *plP = &lP;
PTrack hP;
PTrack *phP = &hP;

Int_t Run, Event;
Int_t prevRun =0;
double Qgs = 0;

TChain* createChain(std::vector<Int_t> runs, std::string Directory)
{
    printf("***Creating New TChain***\n");
  	
	TChain *ch = new TChain("Iris");
	Int_t nRuns = runs.size();
	printf("Adding %d runs!\n",nRuns);
  	for(int i=0; i< nRuns; ++i) {
		std::string filename;
		std::string runNo;
		if(runs.at(i)==0) runNo = "";
		else if(runs.at(i)>99999) runNo = Form("%d",runs.at(i));
		else runNo = Form("0%04d",runs.at(i));
    	filename = Directory + runNo + ".root";
    	
		int addReturn = ch->AddFile(filename.data());
    	TFile tFile(filename.data());
    	if(addReturn == 1 && !tFile.IsZombie()) {
      		printf("Adding run: %s\n",filename.data());
    	} 
		else{
      		printf("WARNING: Tried to add an invalid file: %s ==> Skipping.\n",filename.data());
    	}
  	}
  	printf("----------\n\n");
  	return ch;
}

const int  YdMul= 128;

float ICELossCorr=1.; // has to be implemented!

const int Nchannels = 24;
const int binlimit = 1900;
 		
const Double_t ICLength=22.9*0.062; //cm*mg/cm^3 at 19.5 Torr
const Double_t ICWindow1=0.05*3.44*0.1; //mu*g/cm^3*0.1
const Double_t ICWindow2=0.05*3.44*0.1; //mu*g/cm^3*0.1


Double_t eATgt[100], eAIso[100], eAWndw[100], eAAg[100];	
Double_t dedxATgt[100], dedxAIso[100], dedxAWndw[100], dedxAAg[100];	
Double_t eBSi[100], eBTgt[100], eBSiO2[100], eBB[100], eBP[100], eBAl[100], eBMy[100], eBCsI[100], eBAg[100];	
Double_t dedxBSi[100], dedxBTgt[100], dedxBSiO2[100], dedxBB[100], dedxBP[100], dedxBAl[100], dedxBMy[100], dedxBCsI[100], dedxBAg[100];	
Double_t ebTgt[100], ebSi[100], ebAl[100], ebB[100], ebMy[100], ebP[100], ebCsI[100], ebSiO2[100], ebAg[100];	
Double_t dedxbTgt[100], dedxbSi[100], dedxbAl[100], dedxbB[100], dedxbMy[100], dedxbP[100], dedxbCsI[100], dedxbSiO2[100], dedxbAg[100];	

double EBAC = 0.; //Beam energy from accelerator
runDep runDepPar; // run dependant parameters
nucleus beam; // beam particle properties
nucleus target; // target properties
nucleus lej; // light ejetcile properties
nucleus hej; // heavy ejectile properties  
Int_t useYCalc = 0;//Use calculated YY1 energy  : do we need to ? :Jaspreet
Double_t PResid; //Momentum of residue
Double_t PBeam; // Calculated beam momentum after scattering off Ag
Double_t PA; //Beam momentum before reaction



Double_t A,B,C; //Used for quadratic equations
Double_t MBeam = 0.; // Beam mass
Double_t kBF = 108.904752/(MBeam/931.494013); //Ratio of Beam particle mass and 109-Ag foil nucleus mass
const Double_t MFoil = 931.494013*108.; //Ag foil mass AS
//Double_t geoP.FoilThickness = 5.7;                      /// ?????????? confrim if 5.7 or 5.44 --Jaspreet
Double_t energy = 0;
Double_t cosTheta = 1.;// S3 angle cosine
Double_t EBeam;
Double_t betaCM, gammaCM; //CM velocity

Double_t mA= 0.; //11.;//Beam mass //Reassigned in HandleBOR
Double_t ma = 0.; //target mass
Double_t mb= 0.; //1.;//Light ejectile mass
Double_t mB= 0.; // Heavy ejectile mass

Double_t Pb1,Pb2,PbU,PbUSd; //Light ejectile momentum
Double_t thetaR=sqrt(-1.),thetaR2=sqrt(-1.);
Double_t thetaD=sqrt(-1.),thetaD2=sqrt(-1.);
Double_t thetaCM=sqrt(-1.),thetaCM2=sqrt(-1.);
Double_t ECsI1=sqrt(-1.),ECsI2=sqrt(-1.);
Double_t EYY1=sqrt(-1.),EYY12=sqrt(-1.);
Double_t Eb1=sqrt(-1.), Eb2=sqrt(-1.), EbU=sqrt(-1.), EbUSd=sqrt(-1.);
Double_t EB1=sqrt(-1.), EB2=sqrt(-1.), EBU=sqrt(-1.), EBUSd=sqrt(-1.);
Double_t PB1=sqrt(-1.), PB2=sqrt(-1.), PBU=sqrt(-1.), PBUSd=sqrt(-1.);
Double_t Q1=sqrt(-1.), Q2=sqrt(-1.), QU=sqrt(-1.), QUSd=sqrt(-1.);
Double_t Ex1 = sqrt(-1.), Ex2 = sqrt(-1.), ExU = sqrt(-1.) ;
Double_t Pb1y=sqrt(-1.), Pb2y=sqrt(-1.), PbUy=sqrt(-1.), PbUSdy=sqrt(-1.);
Double_t Pb1xcm=sqrt(-1.), Pb2xcm=sqrt(-1.), PbUxcm=sqrt(-1.), PbUSdxcm=sqrt(-1.);

Int_t CCsI1=-1;
Int_t CCsI2=-1;
	
TCutG *YdCsIGate = NULL;
TCutG *SdGate = NULL;
TCutG *YuGate = NULL;
TCutG *SuGate = NULL;




void getRunPar(Int_t runNo)
{
	FILE *pFile;
	int i = 0;
	double a,b;
	int run_for_corr = 0;
	
	pFile = fopen(runDepPar.runPar.data(), "r");	
	if (pFile == NULL) {
		printf("No run dependent energy and target thickness. Skipping correction.\n");
	}
	else  {
		printf("Reading config file '%s'\n",runDepPar.runPar.data());
	
		while (!feof(pFile)){
			fscanf(pFile,"%d\t%lf\t%lf\n",&i,&a,&b);
			//printf("%d\t%lf\t%lf\n",i,a,b);
			if(i==runNo){ 
				run_for_corr = i;
				runDepPar.energy = a;
			   	geoP.TargetThickness = b;	
			}
		}
		fclose (pFile);	
	
		if(run_for_corr==0){ 
			printf("Run %d not in list. No correction applied!\n",runNo);
		}
		else{
			printf("Run: %d\tBeam energy: %f\tTarget thickness: %f\n",run_for_corr,runDepPar.energy,geoP.TargetThickness);
		}
	}
}

void calculateBeamEnergy(Double_t E)
{
	runDepPar.EBAC = E;
	printf("New Beam Energy: %f\n" ,E);
	Double_t temp_E = E;
	Double_t HalfTargE;
	if(calPhys.boolIC==kTRUE){
		E = E-eloss(E,ICWindow1,eAWndw,dedxAWndw);  
		E = E-eloss(E,ICLength,eAIso,dedxAIso)/ICELossCorr;  
		E = E-eloss(E,ICWindow2,eAWndw,dedxAWndw);  
		printf("Energy loss in IC (including windows): %.3f MeV\n" ,temp_E-E);

		temp_E = E;
	}
	if(geoP.TargetOrientation==kTRUE){//foil facing downstream of target
		E = E-eloss(E,geoP.TargetThickness/2.,eATgt,dedxATgt);  
		printf("Energy loss in half target: %.3f MeV\n" ,temp_E-E);
		printf("Beam energy in center of target: %.3f MeV\n" ,E);
		HalfTargE = E;
		E = E-eloss(E,geoP.TargetThickness/2.,eATgt,dedxATgt);// eloss in remainder of target  

		temp_E = E;

		if(geoP.FoilThickness>0.) E = E-eloss(E,geoP.FoilThickness,eAAg,dedxAAg);  
		else E = temp_E;
		printf("Energy loss in silver foil: %.3f MeV\n" ,temp_E-E);
		printf("Energy after silver foil: %.3f MeV\n",E);
	}
	else{ //foil facing upstream of target
	        if(geoP.FoilThickness>0.) E = E-eloss(E,geoP.FoilThickness,eAAg,dedxAAg);  
		else E = temp_E;
		printf("Energy loss in silver foil: %.3f MeV\n" ,temp_E-E);
		printf("Energy after silver foil: %.3f MeV\n",E);

		temp_E = E;
		E = E-eloss(E,geoP.TargetThickness/2.,eATgt,dedxATgt);  
		printf("Energy loss in half target: %.3f MeV\n" ,temp_E-E);
		printf("Beam energy in center of target: %.3f MeV\n" ,E);
		HalfTargE = E;
	}
	
	runDepPar.energy = HalfTargE;;
	runDepPar.momentum = sqrt(runDepPar.energy*runDepPar.energy+2.*runDepPar.energy*beam.mass);//beam momentum
	runDepPar.beta = runDepPar.momentum/(runDepPar.energy + beam.mass + target.mass);
	runDepPar.gamma = 1./sqrt(1.-runDepPar.beta*runDepPar.beta);

	EBAC = runDepPar.EBAC;
	Qgs = runDepPar.Q ;
	EBeam= runDepPar.energy;
	PA = runDepPar.momentum;//beam momentum
	betaCM = runDepPar.beta;
	gammaCM = runDepPar.gamma;
	printf("Beam momentum in center of target: %.3f MeV\n",PA);
	printf("Beta: %f\tGamma: %f\n",betaCM,gammaCM);
}

//---------------------------------------------------------------------------------
void HandleBOR_PHYSICS(std::string BinPath, std::string Directory, std::string CalibFile, std::string OutputFile)
{
  	printf("In HandleBOR_PHYSICS...\n");	
	if(CalibFile=="") printf("No calibration file specified!\n\n");
	calPhys.Load(CalibFile);
	calPhys.Print();

	Int_t runTmp;
	std::vector<Int_t> runs;
	FILE *rFile;
	rFile = fopen(calPhys.fileRunList.data(), "r");	
	char rLine[256]; 
	if (rFile == NULL || calPhys.boolRunList==false) {
		printf("No list of runs.\n");
   	}  
 	else  {
		printf("Reading run List '%s'\n",calPhys.fileRunList.data());

		while (!feof(rFile)){
			fgets(rLine,256,rFile);
			printf("%s",rLine);
       		sscanf(rLine,"%*s %*s %*s %d",&runTmp);
			if(runs.size()==0 || runTmp>runs.at(runs.size()-1)) runs.push_back(runTmp);
     	}
		fclose(rFile);	
		printf("\n");
 	}

	input_chain = createChain(runs,Directory);
	
	// deactivate some branches, only relevant for simulated data from simIris
	if(input_chain->GetListOfBranches()->FindObject("Evnt")) input_chain->SetBranchStatus("Evnt",0); 
	if(input_chain->GetListOfBranches()->FindObject("reacX")) input_chain->SetBranchStatus("reacX",0); 
	if(input_chain->GetListOfBranches()->FindObject("reacY")) input_chain->SetBranchStatus("reacY",0); 
	if(input_chain->GetListOfBranches()->FindObject("reacZ")) input_chain->SetBranchStatus("reacZ",0); 
	if(input_chain->GetListOfBranches()->FindObject("wght")) input_chain->SetBranchStatus("wght",0); 
	if(input_chain->GetListOfBranches()->FindObject("Qgen")) input_chain->SetBranchStatus("Qgen",0); 
	if(input_chain->GetListOfBranches()->FindObject("qdet")) input_chain->SetBranchStatus("Qdet",0); 
	if(input_chain->GetListOfBranches()->FindObject("ICdE")) input_chain->SetBranchStatus("ICdE",0);
	if(input_chain->GetListOfBranches()->FindObject("SSBdE")) input_chain->SetBranchStatus("SSBdE",0);
	if(input_chain->GetListOfBranches()->FindObject("yd")) input_chain->SetBranchStatus("yd*",0);
	if(input_chain->GetListOfBranches()->FindObject("csi")) input_chain->SetBranchStatus("csi*",0);
	if(input_chain->GetListOfBranches()->FindObject("sd1")) input_chain->SetBranchStatus("sd1*",0);
	if(input_chain->GetListOfBranches()->FindObject("sd2")) input_chain->SetBranchStatus("sd2*",0);
	//***************************************************************************
	
	if(input_chain->GetListOfBranches()->FindObject("det")){
	   	input_chain->SetBranchAddress("det",&det);
	}
	if(input_chain->GetListOfBranches()->FindObject("tdc")){ 
		input_chain->SetBranchAddress("tdc",&tdc);
	}
	if(input_chain->GetListOfBranches()->FindObject("scaler")){
	   	input_chain->SetBranchAddress("scaler",&pscaler);
	}
	if(input_chain->GetListOfBranches()->FindObject("lP")){
	   	input_chain->SetBranchAddress("lP",&plP);
	}
	if(input_chain->GetListOfBranches()->FindObject("hP")){
	   	input_chain->SetBranchAddress("hP",&phP);
	}		
	if(input_chain->GetListOfBranches()->FindObject("Run")){
	   	input_chain->SetBranchAddress("Run",&Run);
	}
	if(input_chain->GetListOfBranches()->FindObject("Event")){
	   	input_chain->SetBranchAddress("Event",&Event);
	}
	treeFile = new TFile(OutputFile.data(),"RECREATE");
	tree=input_chain->CloneTree(0);	
	IrisEvent = new TEvent();
 	tree->Branch("IrisEvent","TEvent",&IrisEvent,32000,99);

	if(calPhys.boolFGate==kTRUE){
		TFile *fYdCsIGates = new TFile(calPhys.fileGate.data());
   		printf("opened file %s\n",calPhys.fileGate.data());
   
		if(calPhys.boolNGate==kTRUE){
			YdCsIGate = (TCutG*)fYdCsIGates->FindObjectAny(calPhys.nameGate.data());
  			if(!YdCsIGate) printf("No Yd/CsI gate.\n");  
			else printf("Grabbed Yd/CsI gate %s.\n",calPhys.nameGate.data());
		}
		else printf("No Yd/CsI gate.\n"); 
	}
	else printf("No Yd/CsI gate.\n");  
	
	if(calPhys.boolFSdGate==kTRUE){
		TFile *fSdGates = new TFile(calPhys.fileSdGate.data());
   		printf("opened file %s\n",calPhys.fileSdGate.data());
   
		if(calPhys.boolNSdGate==kTRUE){
			SdGate = (TCutG*)fSdGates->FindObjectAny(calPhys.nameSdGate.data());
  			if(!SdGate) printf("No S3 gate.\n");  
			else printf("Grabbed S3 gate %s.\n",calPhys.nameSdGate.data());
		}
		else printf("No S3 gate.\n"); 
	}
	else printf("No S3 gate.\n");  

	if(calPhys.boolFYuGate==kTRUE){
		TFile *fYuGates = new TFile(calPhys.fileYuGate.data());
   		printf("opened file %s\n",calPhys.fileYuGate.data());
   
		if(calPhys.boolNYuGate==kTRUE){
			YuGate = (TCutG*)fYuGates->FindObjectAny(calPhys.nameYuGate.data());
  			if(!YuGate) printf("No Yu gate.\n");  
			else printf("Grabbed Yu gate %s.\n",calPhys.nameYuGate.data());
		}
		else printf("No Yu gate.\n"); 
	}
	else printf("No Yu gate.\n");  
	
	if(calPhys.boolFSuGate==kTRUE){
		TFile *fSuGates = new TFile(calPhys.fileSuGate.data());
   		printf("opened file %s\n",calPhys.fileSuGate.data());
   
		if(calPhys.boolNSuGate==kTRUE){
			SuGate = (TCutG*)fSuGates->FindObjectAny(calPhys.nameSuGate.data());
  			if(!SuGate) printf("No upstream S3 gate.\n");  
			else printf("Grabbed upstream S3 gate %s.\n",calPhys.nameSuGate.data());
		}
		else printf("No upstream S3 gate.\n"); 
	}
	else printf("No upstream S3 gate.\n");  
	
	geoP.ReadGeometry(calPhys.fileGeometry.data());
	geoP.Print();

//	ICELossCorr = 1.;	


	//	// Time dependent correction of IC energy loss 
//	FILE *pFile;
//	int Chan = 0;
//	double a,b;
//	int run_for_corr = 0;
//	
//	if (pFile == NULL || calPhys.boolTCorrIC==false) {
//		//fprintf(logFile,"No time dependent correction for IC energy loss. Skipping correction.\n");
//		printf("No time dependent correction for IC energy loss. Skipping correction.\n");
//		ICELossCorr =1.;
//	}
//	else  {
//		printf("Reading config file '%s'\n",calPhys.fileTCorrIC.data());
//
//		while (!feof(pFile)){
//    		fscanf(pFile,"%d%lf%lf",&Chan,&a,&b);
//			if(Chan==run){ 
//				run_for_corr = Chan;
//				ICELossCorr = b; 
//			}
//    	}
//    	fclose (pFile);	
//
//		if(run_for_corr==0){ 
//			printf("Run %d not in list. No correction applied!\n",run);
//		}
//		else{
//			printf("Run: %d\tIC Gain correction: %f\n\n",Chan,ICELossCorr);
//		}
//  	}

	if(calPhys.boolRunDepPar){
		runDepPar.setRunDepPar(calPhys.fileRunDepPar);// setting run dependent parameters.
		runDepPar.Print();

		beam.getInfo(BinPath,runDepPar.nA);
		target.getInfo(BinPath,runDepPar.na);
		hej.getInfo(BinPath,runDepPar.nB);
		lej.getInfo(BinPath,runDepPar.nb);

		mA = beam.mass; //Beam mass //Reassigned in HandleBOR
		ma = target.mass;
		mb = lej.mass; //Light ejectile mass
		mB = hej.mass;
		printf("mA=%f, ma=%f, mb=%f, mB=%f\n",mA,ma,mb,mB); 
		kBF = MFoil/mA;
	
		printf("Beam energy: %f\n", runDepPar.energy);
		// printf("Target thickness: %f\n",geoP.TargetThickness);
	
		if(calPhys.boolIdedx==kTRUE){
			printf("\n\nLoading dedx Graphs for incoming %s ...\n",runDepPar.nA.data());
			dedx_i.Load(calPhys.fileIdedx);
			dedx_i.Print();
			if(dedx_i.boolAg==kTRUE) loadELoss(dedx_i.Ag,eAAg,dedxAAg,mA);	
			if(dedx_i.boolIso==kTRUE) loadELoss(dedx_i.Iso,eAIso,dedxAIso,mA);	
			if(dedx_i.boolWndw==kTRUE) loadELoss(dedx_i.Wndw,eAWndw,dedxAWndw,mA);	
			if(dedx_i.boolTgt==kTRUE) loadELoss(dedx_i.Tgt,eATgt,dedxATgt,mA);	
		}

		if(calPhys.boolLdedx==kTRUE){
			printf("\n\nLoading dedx Graphs for target-like %s...\n",runDepPar.nb.data());
			dedx_l.Load(calPhys.fileLdedx);
			dedx_l.Print();
			if(dedx_l.boolSi==kTRUE) loadELoss(dedx_l.Si,ebSi,dedxbSi,mb);	
			if(dedx_l.boolAl==kTRUE) loadELoss(dedx_l.Al,ebAl,dedxbAl,mb);	
			if(dedx_l.boolB==kTRUE) loadELoss(dedx_l.B,ebB,dedxbB,mb);	
			if(dedx_l.boolP==kTRUE) loadELoss(dedx_l.P,ebP,dedxbP,mb);	
			if(dedx_l.boolMy==kTRUE) loadELoss(dedx_l.My,ebMy,dedxbMy,mb);	
			if(dedx_l.boolTgt==kTRUE) loadELoss(dedx_l.Tgt,ebTgt,dedxbTgt,mb);	
			if(dedx_l.boolAg==kTRUE) loadELoss(dedx_l.Ag,ebAg,dedxbAg,mb);
			if(dedx_l.boolSiO2==kTRUE) loadELoss(dedx_l.SiO2,ebSiO2,dedxbSiO2,mb);	
		}

		if(calPhys.boolHdedx==kTRUE){
			printf("\n\nLoading dedx Graphs for beam-like %s...\n",runDepPar.nB.data());
			dedx_h.Load(calPhys.fileHdedx);
			dedx_h.Print();
			if(dedx_h.boolTgt==kTRUE) loadELoss(dedx_h.Tgt,eBTgt,dedxBTgt,mB);	
			if(dedx_h.boolSi==kTRUE) loadELoss(dedx_h.Si,eBSi,dedxBSi,mB);	
			if(dedx_h.boolAl==kTRUE) loadELoss(dedx_h.Al,eBAl,dedxBAl,mB);	
			if(dedx_h.boolB==kTRUE) loadELoss(dedx_h.B, eBB,dedxBB,mB);	
			if(dedx_h.boolP==kTRUE) loadELoss(dedx_h.P, eBP,dedxBP,mB);	
			if(dedx_h.boolSiO2==kTRUE) loadELoss(dedx_h.SiO2,eBSiO2,dedxBSiO2,mB);
			if(dedx_h.boolAg==kTRUE) loadELoss(dedx_h.Ag, eBAg,dedxBAg,mB);
			if(dedx_h.boolMy==kTRUE) loadELoss(dedx_h.My,eBMy,dedxBMy,mB);			
		}

		// Initialize runPar with values from first run in chain	
		if(runDepPar.bool_runPar==kTRUE) getRunPar(runs.at(0));
		prevRun=runs.at(0);
		calculateBeamEnergy(runDepPar.energy);
		
		printf("MBeam: %f\t MFoil: %f\t kBF: %f\n",beam.mass,MFoil,kBF);
		printf("beam mass: %f MeV (%f)\ttarget mass: %f MeV (%f)\n",mA,mA/931.494061,ma,ma/931.494061);
		printf("heavy ejectile mass: %f MeV (%f)\tlight ejectile mass: %f MeV (%f)\n",mB,mB/931.494061,mb,mb/931.494061);
	}
 	if (calPhys.boolICGates==kFALSE){
		runDepPar.ICmin=0;
		runDepPar.ICmax=4096;
	}
	
//--------------------------------------------------------------------------------
		printf("End of HandleBOR_Physics\n");
} //HandleBOR_Physics

void HandlePHYSICS()
{
  	Int_t nEntries = input_chain->GetEntries();
  	printf("%d entries in total.\n",nEntries);
  	for(Int_t i=0; i<nEntries; i++)
    {
		IrisEvent->Clear();
        Long64_t check_entry = input_chain->LoadTree(i);
		if(check_entry<0) break;
		input_chain->GetEntry(i);
		if((i%100)==0) printf("Processing event %d\r",i);
		
		if (det->TICEnergy.size()==0) continue; 
		if (det->TICEnergy.at(0)<runDepPar.ICmin || det->TICEnergy.at(0)>runDepPar.ICmax) continue;
		if(YdCsIGate!=NULL||SdGate!=NULL||YuGate!=NULL||SuGate!=NULL) //Skip if NO Gates
		{Bool_t conditionsmet = 0;
		if (YdCsIGate!=NULL&&(det->TYdMul==0)) continue; // event has YY1
		for(int z=0;z<det->TYdMul;z++){
		if (YdCsIGate!=NULL&&calPhys.numGate!=2&&det->TCsI1Mul>0&&YdCsIGate->IsInside(det->TCsI1Energy.at(z),det->TYdEnergy.at(z)*cos(TMath::DegToRad()*det->TYdTheta.at(z)))==1) conditionsmet = 1; // event in proton/deuteron/etc YdCsIGate?
		break;}
		for(int z=0;z<det->TYdMul;z++){
		if (YdCsIGate!=NULL&&calPhys.numGate==2&&det->TCsI2Mul>0&&YdCsIGate->IsInside(det->TCsI2Energy.at(z),det->TYdEnergy.at(z)*cos(TMath::DegToRad()*det->TYdTheta.at(z)))==1) conditionsmet = 1; // event in proton/deuteron/etc YdCsIGate?
		break;}
		if (SdGate!=NULL&&(det->TSd1rEnergy.size()==0||det->TSd1sEnergy.size()==0)) continue; // event has S3 hit?
		for(int z=0;z<det->TSd1rMul;z++){
		if (SdGate!=NULL&&det->TSd2rMul==det->TSd1rMul&&SdGate->IsInside(det->TSd2rEnergy.at(z),det->TSd1rEnergy.at(z)*cos(TMath::DegToRad()*det->TSd1Theta.at(z)))==1) conditionsmet = 1; // event in proton/deuteron/etc SdGate?
		break;}
		if (YuGate!=NULL&&det->TYuEnergy.size()==0) continue; // event has Yu hit?
		if (YuGate!=NULL&&YuGate->IsInside(det->TYuTheta.at(0),det->TYuEnergy.at(0))==0) continue; // event in proton/deuteron/etc YuGate?
		if (SuGate!=NULL&&(det->TSurEnergy.size()==0||det->TSusEnergy.size()==0)) continue; // event has Su hit?
		if (SuGate!=NULL&&SuGate->IsInside(det->TSuTheta.at(0),det->TSurEnergy.at(0))==0) continue; // event in proton/deuteron/etc SuGate?
		if (conditionsmet == 0) continue;}
		if(runDepPar.bool_runPar == kTRUE && Run != prevRun){
		  	getRunPar(Run);
		  	calculateBeamEnergy(runDepPar.energy);
		  	prevRun = Run;
		}
		
		IrisEvent->fEBAC = EBAC;
		IrisEvent->fmA = mA;
		IrisEvent->fma = ma;
		IrisEvent->fmB = mB;
		IrisEvent->fmb = mb;
		IrisEvent->fkBF = kBF;
		IrisEvent->fEBeam = EBeam;
		IrisEvent->fbetaCM = betaCM;
		IrisEvent->fgammaCM = gammaCM;
		IrisEvent->fPA = PA;

		//adding dead layer energy losses for Target-like Particle
		if(det->TSd1rEnergy.size()>0 && det->TSd1sEnergy.size()>0){
				for(int z=0;z<det->TSd1rMul;z++){
				Double_t sd1theta=TMath::DegToRad()* (det->TSd1Theta.at(z));
				cosTheta = cos(sd1theta);
	  			//Sd2 ring side
	  			energy = 0;
				if(det->TSd1rMul==det->TSd2rMul && det->TSd2rEnergy.at(z)>0. && det->TSd2sEnergy.at(z)>0){
				energy = det->TSd2rEnergy.at(z);
                energy = energy+elossFi(energy,0.1*1.822*0.5/cosTheta,ebP,dedxbP); //phosphorus implant
                energy = energy+elossFi(energy,0.1*2.7*0.3/cosTheta,ebAl,dedxbAl); //first metal
                energy = energy+elossFi(energy,0.1*2.7*0.3/cosTheta,ebAl,dedxbAl); //first metal
                energy = energy+elossFi(energy,0.1*1.822*0.5/cosTheta,ebP,dedxbP); //phosphorus implant
				}                
                energy = energy + det->TSd1rEnergy.at(z);// energy lost and measured in Sd1
    
                energy = energy+elossFi(energy,0.1*2.35*0.5/cosTheta,ebB,dedxbB); //boron junction implant
	  			energy = energy+elossFi(energy,0.1*2.7*0.3/cosTheta,ebAl,dedxbAl); //first metal
	  			energy = energy+elossFi(energy,0.1*2.65*3.5/cosTheta,ebSiO2,dedxbSiO2); //SiO2
	  			energy = energy+elossFi(energy,0.1*2.7*1.5/cosTheta,ebAl,dedxbAl); //second metal
			
				//target/foil Orientation
                if(geoP.TargetOrientation==kTRUE){ //Foil downstream of target
                  energy = energy+elossFi(energy,geoP.FoilThickness/cosTheta,eBAg,dedxBAg); // Ag foil energyloss
                }
				energy = energy+elossFi(energy,geoP.TargetThickness/2./cosTheta,eBTgt,dedxBTgt); // SHT energy loss
				
			PResid = sqrt(2.*det->TSdETot*mA);     //Beam momentum in MeV/c
			A = kBF-1.;                              //Quadratic equation parameters
			B = 2.0*PResid* cos(sd1theta);
			C = -1.*(kBF+1)*PResid*PResid; 
			if (A!=0)    PBeam = (sqrt(B*B-4.*A*C)-B)/(2*A);
			IrisEvent->fPBeam_tlP.push_back(PBeam);
			IrisEvent->fPResid_tlP.push_back(PResid);
			IrisEvent->fA_tlP.push_back(A);
			IrisEvent->fB_tlP.push_back(B);
			IrisEvent->fC_tlP.push_back(C);
			IrisEvent->TSdETot_tlP.push_back(energy);
			//to calculate residue energy from beam
			IrisEvent->fEB_tlP.push_back(energy);
			//IrisEvent->fEB=  IrisEvent->fEB + elossFi(det->TSdETot,geoP.FoilThickness/2.,eAAg,dedxAAg); //energy loss from the end of H2 to the center of Ag.
			
		        Double_t SdThetaCM = TMath::RadToDeg()*atan(tan(sd1theta)/sqrt(gammaCM-gammaCM*betaCM*(mA+energy)/(PBeam*cos(sd1theta))));// check if this is still correct for H2 target tk
				IrisEvent->TSdThetaCM_tlP.push_back(SdThetaCM);
			}
			}		
		
			// Qvalue calculation for Target like particle. CsI1Energy
		if (det->TYdEnergy.size()>0&&det->TYdRing.size()>0&&det->TYdMul>0) {    //check if in the proton/deuteron YdCsIGate
			
			for(int z=0;z<det->TYdMul;z++){

			thetaR2=det->TYdTheta.at(z)*TMath::Pi()/180.;//randomized angles from treeIris (or simIris)
			thetaD2 = thetaR2*TMath::RadToDeg();
			IrisEvent->fTheta1_tlP.push_back(thetaD2);
			EYY1 = det->TYdEnergy.at(z);
			ECsI1 = 0;
			
			if(det->TCsI1Energy.at(z)>0.01){
			CCsI1= det->TCsI1Channel.at(z);
			ECsI1= det->TCsI1Energy.at(z);
			
			ECsI1= ECsI1+elossFi(ECsI1,0.1*1.4*6./cos(thetaR2),ebMy,dedxbMy); //Mylar
			ECsI1= ECsI1+elossFi(ECsI1,0.1*2.702*0.3/cos(thetaR2),ebAl,dedxbAl); //0.3 u Al
			ECsI1= ECsI1+elossFi(ECsI1,0.1*1.822*0.1/cos(thetaR2),ebP,dedxbP);} // 0.1Phosphorus			
			Eb2= ECsI1+EYY1; //use measured Yd // change june28;

			Eb2= Eb2+elossFi(Eb2,0.1*2.35*0.05/cos(thetaR2),ebB,dedxbB); //0.05 u B 
			Eb2= Eb2+elossFi(Eb2,0.1*2.702*0.1/cos(thetaR2),ebAl,dedxbAl); //0.1 u Al
			IrisEvent->fEYd1_tlP.push_back(Eb2-ECsI1);
			Eb2= Eb2+elossFi(Eb2,geoP.TargetThickness/2./cos(thetaR2),ebTgt,dedxbTgt); //deuteron energy midtarget
			
			if(geoP.TargetOrientation==kTRUE){ //Foil downstream of target
			  Eb2 = Eb2+elossFi(Eb2,geoP.FoilThickness/cos(thetaR2),ebAg,dedxbAg); //
			}

			IrisEvent->YdCsI1ETot_tlP.push_back(Eb2);
			Pb2 = sqrt(Eb2*Eb2+2.*Eb2*mb);
			Pb2y = Pb2*sin(thetaR2);
			Pb2xcm = gammaCM*betaCM*(Eb2+mb)- gammaCM*Pb2*cos(thetaR2);
			EB2 = EBeam+mA+ma-Eb2-mb;
			PB2 = sqrt(PA*PA+Pb2*Pb2-2.*PA*Pb2*cos(thetaR2));
			//Q2 = mA+ma-mb- sqrt(mA*mA+mb*mb-ma*ma-2.*(mA+EBeam)*(mb+Eb2)+2.*PA*Pb2*cos(thetaR2)+2.*(EBeam+mA+ma-Eb2-mb)*ma);  //Alisher's equation 
			Q2 = mA+ma-mb-sqrt(EB2*EB2-PB2*PB2); //Equivalent to the previous equation
			Ex2 = Qgs - Q2; 
			IrisEvent->fCCsI1_tlP.push_back(CCsI1);
			IrisEvent->fECsI1_tlP.push_back(ECsI1);
			IrisEvent->fEb1_tlP.push_back(Eb2);
			IrisEvent->fPb1_tlP.push_back(Pb2);
			IrisEvent->EB1_det_tlP.push_back(EB2);
			IrisEvent->PB1_det_tlP.push_back(PB2);
			IrisEvent->fPby1_tlP.push_back(Pb2y);
			IrisEvent->fPbxcm1_tlP.push_back(Pb2xcm);
			IrisEvent->Qdet1_tlP.push_back(Q2);
			IrisEvent->Ex1_tlP.push_back(Ex2);
			thetaCM = TMath::RadToDeg()*atan(Pb2y/Pb2xcm);
			thetaCM = (thetaCM<0) ? thetaCM+180. : thetaCM;
			IrisEvent->fThetacm1_tlP.push_back(thetaCM);
			
		}
	}

		// Qvalue calculation for Target like particle. CsI2Energy
		if (det->TYdEnergy.size()>0&&det->TYdRing.size()>0&&det->TYdMul>0) {    //check if in the proton/deuteron YdCsIGate

			for(int z=0;z<det->TYdMul;z++){

			thetaR2=det->TYdTheta.at(z)*TMath::Pi()/180.;//randomized angles from treeIris (or simIris)
			thetaD2 = thetaR2*TMath::RadToDeg();
			IrisEvent->fTheta2_tlP.push_back(thetaD2);
			EYY1 = det->TYdEnergy.at(z);
			ECsI2 = 0;
     
			if(det->TCsI2Energy.at(z)>0.01){
			CCsI2= det->TCsI2Channel.at(z);
			ECsI2= det->TCsI2Energy.at(z);
			
			ECsI2= ECsI2+elossFi(ECsI2,0.1*1.4*6./cos(thetaR2),ebMy,dedxbMy); //Mylar
			ECsI2= ECsI2+elossFi(ECsI2,0.1*2.702*0.3/cos(thetaR2),ebAl,dedxbAl); //0.3 u Al
			ECsI2= ECsI2+elossFi(ECsI2,0.1*1.822*0.1/cos(thetaR2),ebP,dedxbP);} // 0.1Phosphorus			
			Eb2= ECsI2+EYY1; //use measured Yd // change june28
			
			Eb2= Eb2+elossFi(Eb2,0.1*2.35*0.05/cos(thetaR2),ebB,dedxbB); //0.05 u B 
			Eb2= Eb2+elossFi(Eb2,0.1*2.702*0.1/cos(thetaR2),ebAl,dedxbAl); //0.1 u Al
			IrisEvent->fEYd2_tlP.push_back(Eb2-ECsI2);
			Eb2= Eb2+elossFi(Eb2,geoP.TargetThickness/2./cos(thetaR2),ebTgt,dedxbTgt); //deuteron energy midtarget
			
			if(geoP.TargetOrientation==kTRUE){ //Foil downstream of target
			  Eb2 = Eb2+elossFi(Eb2,geoP.FoilThickness/cos(thetaR2),ebAg,dedxbAg); //
			}

			IrisEvent->YdCsI2ETot_tlP.push_back(Eb2);
			Pb2 = sqrt(Eb2*Eb2+2.*Eb2*mb);
			Pb2y = Pb2*sin(thetaR2);
			Pb2xcm = gammaCM*betaCM*(Eb2+mb)- gammaCM*Pb2*cos(thetaR2);
			EB2 = EBeam+mA+ma-Eb2-mb;
			PB2 = sqrt(PA*PA+Pb2*Pb2-2.*PA*Pb2*cos(thetaR2));
			//Q2 = mA+ma-mb- sqrt(mA*mA+mb*mb-ma*ma-2.*(mA+EBeam)*(mb+Eb2)+2.*PA*Pb2*cos(thetaR2)+2.*(EBeam+mA+ma-Eb2-mb)*ma);  //Alisher's equation 
			Q2 = mA+ma-mb-sqrt(EB2*EB2-PB2*PB2); //Equivalent to the previous equation
			Ex2 = Qgs- Q2; 
			IrisEvent->fCCsI2_tlP.push_back(CCsI2);
			IrisEvent->fECsI2_tlP.push_back(ECsI2);
			IrisEvent->fEb2_tlP.push_back(Eb2);
			IrisEvent->fPb2_tlP.push_back(Pb2);
			IrisEvent->EB2_det_tlP.push_back(EB2);
			IrisEvent->PB2_det_tlP.push_back(PB2);
			IrisEvent->fPby2_tlP.push_back(Pb2y);
			IrisEvent->fPbxcm2_tlP.push_back(Pb2xcm);
			IrisEvent->Qdet2_tlP.push_back(Q2);
			IrisEvent->Ex2_tlP.push_back(Ex2);
			thetaCM2 = TMath::RadToDeg()*atan(Pb2y/Pb2xcm);
			thetaCM2 = (thetaCM2<0) ? thetaCM2+180. : thetaCM2;
			IrisEvent->fThetacm2_tlP.push_back(thetaCM2);
		}
	}
	
		//adding dead layer energy losses for Beam-like Particle.
		if(det->TSd1rEnergy.size()>0 && det->TSd1sEnergy.size()>0 ){
				for(int z=0;z<det->TSd1rMul;z++){
				Double_t sd1theta=TMath::DegToRad()* (det->TSd1Theta.at(z));
				cosTheta = cos(sd1theta);
	  			//Sd2 ring side
	  			energy = 0;
				if(det->TSd1rMul==det->TSd2rMul && det->TSd2rEnergy.at(z)>0. && det->TSd2sEnergy.at(z)>0.){
				energy = det->TSd2rEnergy.at(z);
                energy = energy+elossFi(energy,0.1*1.822*0.5/cosTheta,eBP,dedxBP); //phosphorus implant
                energy = energy+elossFi(energy,0.1*2.7*0.3/cosTheta,eBAl,dedxBAl); //first metal
                energy = energy+elossFi(energy,0.1*2.7*0.3/cosTheta,eBAl,dedxBAl); //first metal
                energy = energy+elossFi(energy,0.1*1.822*0.5/cosTheta,eBP,dedxBP); //phosphorus implant
				}                
                energy = energy + det->TSd1rEnergy.at(z);// energy lost and measured in Sd1
    
                energy = energy+elossFi(energy,0.1*2.35*0.5/cosTheta,eBB,dedxBB); //boron junction implant
	  			energy = energy+elossFi(energy,0.1*2.7*0.3/cosTheta,eBAl,dedxBAl); //first metal
	  			energy = energy+elossFi(energy,0.1*2.65*3.5/cosTheta,eBSiO2,dedxBSiO2); //SiO2
	  			energy = energy+elossFi(energy,0.1*2.7*1.5/cosTheta,eBAl,dedxBAl); //second metal
			
				//target/foil Orientation
                if(geoP.TargetOrientation==kTRUE){ //Foil downstream of target
                  energy = energy+elossFi(energy,geoP.FoilThickness/cosTheta,eBAg,dedxBAg); // Ag foil energyloss
                }
				energy = energy+elossFi(energy,geoP.TargetThickness/2./cosTheta,eBTgt,dedxBTgt); // SHT energy loss
				
			PResid = sqrt(2.*det->TSdETot*mA);     //Beam momentum in MeV/c
			A = kBF-1.;                              //Quadratic equation parameters
			B = 2.0*PResid* cos(sd1theta);
			C = -1.*(kBF+1)*PResid*PResid; 
			if (A!=0)    PBeam = (sqrt(B*B-4.*A*C)-B)/(2*A);
			IrisEvent->fPBeam_blP.push_back(PBeam);
			IrisEvent->fPResid_blP.push_back(PResid);
			IrisEvent->fA_blP.push_back(A);
			IrisEvent->fB_blP.push_back(B);
			IrisEvent->fC_blP.push_back(C);
			IrisEvent->TSdETot_blP.push_back(energy);
			//to calculate residue energy from beam
			IrisEvent->fEB_blP.push_back(energy);
			//IrisEvent->fEB=  IrisEvent->fEB + elossFi(det->TSdETot,geoP.FoilThickness/2.,eAAg,dedxAAg); //energy loss from the end of H2 to the center of Ag.
			
		        Double_t SdThetaCM = TMath::RadToDeg()*atan(tan(sd1theta)/sqrt(gammaCM-gammaCM*betaCM*(mA+energy)/(PBeam*cos(sd1theta))));// check if this is still correct for H2 target tk
				IrisEvent->TSdThetaCM_blP.push_back(SdThetaCM);
			}
			}

		// Qvalue calculation for Beam like particle. CsI1Energy
		if (det->TYdEnergy.size()>0&&det->TYdRing.size()>0&&det->TYdMul>0) {    //check if in the proton/deuteron YdCsIGate
			for(int z=0;z<det->TYdMul;z++){

			thetaR=det->TYdTheta.at(z)*TMath::Pi()/180.;//randomized angles from treeIris (or simIris)
			thetaD = thetaR*TMath::RadToDeg();
			IrisEvent->fTheta1_blP.push_back(thetaD);
			EYY1 = det->TYdEnergy.at(z);
			ECsI1 = 0;
     
		   //check if in the proton/deuteron YdCsIGate
			if(det->TCsI1Energy.at(z)>0.01){
			CCsI1= det->TCsI1Channel.at(z);
			ECsI1= det->TCsI1Energy.at(z);

			ECsI1= ECsI1+elossFi(ECsI1,0.1*1.397*6./cos(thetaR),eBMy,dedxBMy); //Mylar
			ECsI1= ECsI1+elossFi(ECsI1,0.1*2.702*0.3/cos(thetaR),eBAl,dedxBAl); //0.3 u Al
			ECsI1= ECsI1+elossFi(ECsI1,0.1*1.8219*0.1/cos(thetaR),eBP,dedxBP);} // 0.1Phosphorus			
			Eb1= ECsI1+EYY1; //use measured Yd // change june28

			Eb1= Eb1+elossFi(Eb1,0.1*2.3502*0.05/cos(thetaR),eBB,dedxBB); //0.05 u B 
			Eb1= Eb1+elossFi(Eb1,0.1*2.702*0.1/cos(thetaR),eBAl,dedxBAl); //0.1 u Al
			IrisEvent->fEYd1_blP.push_back(Eb1-ECsI1);
			Eb1= Eb1+elossFi(Eb1,geoP.TargetThickness/2./cos(thetaR),eBTgt,dedxBTgt); //deuteron energy midtarget

			if(geoP.TargetOrientation==kTRUE){ //Foil downstream of target
			  Eb1 = Eb1+elossFi(Eb1,geoP.FoilThickness/cos(thetaR),eBAg,dedxBAg); //
			}

			IrisEvent->YdCsI1ETot_blP.push_back(Eb1);
			Pb1 = sqrt(Eb1*Eb1+2.*Eb1*mB);
			Pb1y = Pb1*sin(thetaR);
			Pb1xcm = gammaCM*betaCM*(Eb1+mB)- gammaCM*Pb1*cos(thetaR);
			EB1 = EBeam+mA+ma-Eb1-mB;
			PB1 = sqrt(PA*PA+Pb1*Pb1-2.*PA*Pb1*cos(thetaR));
			//Q1 = mA+ma-mb- sqrt(mA*mA+mb*mb-ma*ma-2.*(mA+EBeam)*(mb+Eb1)+2.*PA*Pb1*cos(thetaR)+2.*(EBeam+mA+ma-Eb1-mb)*ma);  //Alisher's equation 
			Q1 = mA+ma-mB-sqrt(EB1*EB1-PB1*PB1); //Equivalent to the previous equation
			//Ex1 = 12 - Q1; 
			IrisEvent->fCCsI1_blP.push_back(CCsI1);
			IrisEvent->fECsI1_blP.push_back(ECsI1);
			IrisEvent->fEb1_blP.push_back(Eb1);
			IrisEvent->fPb1_blP.push_back(Pb1);
			IrisEvent->EB1_det_blP.push_back(EB1);
			IrisEvent->PB1_det_blP.push_back(PB1);
			IrisEvent->fPby1_blP.push_back(Pb1y);
			IrisEvent->fPbxcm1_blP.push_back(Pb1xcm);
			IrisEvent->Qdet1_blP.push_back(Q1);
			IrisEvent->Ex1_blP.push_back(Ex1);
			thetaCM = TMath::RadToDeg()*atan(Pb1y/Pb1xcm);
			thetaCM = (thetaCM<0) ? thetaCM+180. : thetaCM;
			IrisEvent->fThetacm1_blP.push_back(thetaCM);
		}
	}

	// Qvalue calculation for Beam like particle. CsI2Energy
	if (det->TYdEnergy.size()>0&&det->TYdRing.size()>0&&det->TYdMul>0) {    //check if in the proton/deuteron YdCsIGate
			for(int z=0;z<det->TYdMul;z++){

			thetaR=det->TYdTheta.at(z)*TMath::Pi()/180.;//randomized angles from treeIris (or simIris)
			thetaD = thetaR*TMath::RadToDeg();
			IrisEvent->fTheta2_blP.push_back(thetaD);
			EYY1 = det->TYdEnergy.at(z);
			ECsI1 = 0;
     
		   //check if in the proton/deuteron YdCsIGate
			if(det->TCsI2Energy.at(z)>0.01){
			CCsI1= det->TCsI2Channel.at(z);
			ECsI1= det->TCsI2Energy.at(z);

			ECsI1= ECsI1+elossFi(ECsI1,0.1*1.397*6./cos(thetaR),eBMy,dedxBMy); //Mylar
			ECsI1= ECsI1+elossFi(ECsI1,0.1*2.702*0.3/cos(thetaR),eBAl,dedxBAl); //0.3 u Al
			ECsI1= ECsI1+elossFi(ECsI1,0.1*1.8219*0.1/cos(thetaR),eBP,dedxBP);} // 0.1Phosphorus			
			Eb1= ECsI1+EYY1; //use measured Yd // change june28

			Eb1= Eb1+elossFi(Eb1,0.1*2.3502*0.05/cos(thetaR),eBB,dedxBB); //0.05 u B 
			Eb1= Eb1+elossFi(Eb1,0.1*2.702*0.1/cos(thetaR),eBAl,dedxBAl); //0.1 u Al
			IrisEvent->fEYd2_blP.push_back(Eb1-ECsI1);
			Eb1= Eb1+elossFi(Eb1,geoP.TargetThickness/2./cos(thetaR),eBTgt,dedxBTgt); //deuteron energy midtarget

			if(geoP.TargetOrientation==kTRUE){ //Foil downstream of target
			  Eb1 = Eb1+elossFi(Eb1,geoP.FoilThickness/cos(thetaR),eBAg,dedxBAg); //
			}

			IrisEvent->YdCsI2ETot_blP.push_back(Eb1);
			Pb1 = sqrt(Eb1*Eb1+2.*Eb1*mB);
			Pb1y = Pb1*sin(thetaR);
			Pb1xcm = gammaCM*betaCM*(Eb1+mB)- gammaCM*Pb1*cos(thetaR);
			EB1 = EBeam+mA+ma-Eb1-mB;
			PB1 = sqrt(PA*PA+Pb1*Pb1-2.*PA*Pb1*cos(thetaR));
			//Q1 = mA+ma-mb- sqrt(mA*mA+mb*mb-ma*ma-2.*(mA+EBeam)*(mb+Eb1)+2.*PA*Pb1*cos(thetaR)+2.*(EBeam+mA+ma-Eb1-mb)*ma);  //Alisher's equation 
			Q1 = mA+ma-mB-sqrt(EB1*EB1-PB1*PB1); //Equivalent to the previous equation
			//Ex1 = 12 - Q1; 
			IrisEvent->fCCsI2_blP.push_back(CCsI1);
			IrisEvent->fECsI2_blP.push_back(ECsI1);
			IrisEvent->fEb2_blP.push_back(Eb1);
			IrisEvent->fPb2_blP.push_back(Pb1);
			IrisEvent->EB2_det_blP.push_back(EB1);
			IrisEvent->PB2_det_blP.push_back(PB1);
			IrisEvent->fPby2_blP.push_back(Pb1y);
			IrisEvent->fPbxcm2_blP.push_back(Pb1xcm);
			IrisEvent->Qdet2_blP.push_back(Q1);
			IrisEvent->Ex2_blP.push_back(Ex1);
			thetaCM = TMath::RadToDeg()*atan(Pb1y/Pb1xcm);
			thetaCM = (thetaCM<0) ? thetaCM+180. : thetaCM;
			IrisEvent->fThetacm2_blP.push_back(thetaCM);
		}
	} 

	// Upstream
		// Upstream
		// Upstream
		
		if (det->TYuEnergy.size()>0&&det->TYuRing.size()>0) {    //check if in the proton/deuteron YuGate
		        //thetaR = atan2(geoP.YdInnerRadius+((det->TYuRing.at(0)+0.5)*(geoP.YdOuterRadius-geoP.YdInnerRadius)/16),geoP.YuDistance);
			thetaR=det->TYuTheta.at(0)*TMath::Pi()/180.;//randomized angles from treeIris (or simIris)
			thetaD = thetaR*TMath::RadToDeg();
			IrisEvent->fThetaDU = thetaD;
			cosTheta = cos(TMath::Pi()-thetaR);
			
			EbU=  det->TYuEnergy.at(0); //use measured Yd // change june28
			
			EbU= EbU+elossFi(EbU,0.1*2.35*0.05/cosTheta,ebB,dedxbB); //0.05 u B 
			EbU= EbU+elossFi(EbU,0.1*2.70*0.1/cosTheta,ebAl,dedxbAl); //0.1 u Al
			EbU= EbU+elossFi(EbU,geoP.TargetThickness/2./cosTheta,ebTgt,dedxbTgt); //deuteron energy midtarget
			
			if(geoP.TargetOrientation==kFALSE){ //Foil upstream of target
			  EbU = EbU+elossFi(EbU,geoP.FoilThickness/cosTheta,ebAg,dedxbAg); //
			}

			PbU = sqrt(EbU*EbU+2.*EbU*mb);
			PbUy = PbU*sin(thetaR);
			PbUxcm = gammaCM*betaCM*(EbU+mb)- gammaCM*PbU*cos(thetaR);
			EBU = EBeam+mA+ma-EbU-mb;
			PBU = sqrt(PA*PA+PbU*PbU-2.*PA*PbU*cos(thetaR));
			QU = mA+ma-mb-sqrt(EBU*EBU-PBU*PBU);
			//ExU = Qgs - QU ; 
			IrisEvent->fEbU = EbU;
			IrisEvent->fPbU = PbU;
			IrisEvent->fEBU = EBU;
			IrisEvent->fPBU = PBU;
			IrisEvent->fPbUy = PbUy;
			IrisEvent->fPbUxcm = PbUxcm;
			IrisEvent->fQvU = QU;
			IrisEvent->fExU = ExU ;
			thetaCM = TMath::RadToDeg()*atan(PbUy/PbUxcm);
			thetaCM = (thetaCM<0) ? thetaCM+180. : thetaCM;
			IrisEvent->fThetacmU = thetaCM;
		}

		if(det->TSurEnergy.size()>0){
		        //thetaR = atan2(geoP.SdInnerRadius+((geoP.SdOuterRadius-geoP.SdInnerRadius)/24.)*(det->TSurChannel.at(0)+1.)-0.5,geoP.SuDistance);
			thetaR=det->TSuTheta.at(0)*TMath::Pi()/180.;//randomized angles from treeIris (or simIris)
			thetaD = TMath::RadToDeg()*thetaR;
			IrisEvent->fThetaDUSd = thetaD;
			cosTheta = cos(TMath::Pi()-thetaR);

			EbUSd = det->TSurEnergy.at(0);
			//sector side
			EbUSd = EbUSd+elossFi(EbUSd,0.1*1.822*0.5/cosTheta,ebP,dedxbP); //phosphorus implant
		        EbUSd = EbUSd+elossFi(EbUSd,0.1*2.7*0.3/cosTheta,ebAl,dedxbAl); //metal
			EbUSd = EbUSd+elossFi(EbUSd,geoP.TargetThickness/2./cosTheta,ebTgt,dedxbTgt); //deuteron energy midtarget

			if(geoP.TargetOrientation==kFALSE){ //Foil upstream of target
			  EbUSd = EbUSd+elossFi(EbUSd,geoP.FoilThickness/cosTheta,ebAg,dedxbAg); //
			}
					  
			PbUSd = sqrt(EbUSd*EbUSd+2.*EbUSd*mb);
			PbUSdy = PbUSd*sin(thetaR);
			PbUSdxcm = gammaCM*betaCM*(EbUSd+mb)- gammaCM*PbUSd*cos(thetaR);
			EBUSd = EBeam+mA+ma-EbUSd-mb;
			PBUSd = sqrt(PA*PA+PbUSd*PbUSd-2.*PA*PbUSd*cos(thetaR));
			QUSd = mA+ma-mb-sqrt(EBUSd*EBUSd-PBUSd*PBUSd);
			IrisEvent->fEbUSd = EbUSd;
			IrisEvent->fPbUSd = PbUSd;
			IrisEvent->fEBUSd = EBUSd;
			IrisEvent->fPBUSd = PBUSd;
			IrisEvent->fPbUSdy = PbUSdy;
			IrisEvent->fPbUSdxcm = PbUSdxcm;
			IrisEvent->fQvUSd = QUSd;
			thetaCM = TMath::RadToDeg()*atan(PbUSdy/PbUSdxcm);
			thetaCM = (thetaCM<0) ? thetaCM+180. : thetaCM;
			IrisEvent->fThetacmUSd = thetaCM;
		}
   
		tree->Fill();
	} // end event loop
	
}

void HandleEOR_PHYSICS()
{
  	printf(" in Physics EOR\n");
	tree->AutoSave();
	treeFile->Close();
}
 
