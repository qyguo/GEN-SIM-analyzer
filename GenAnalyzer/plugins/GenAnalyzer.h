/*
 * =====================================================================================
 *
 *       Filename:  GenAnalyzer.h
 *
 *    Description:
 *
 *        Version:  1.0
 *        Created:  12.10.2016 09:03:18
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Ramkrishna Sharma (Ram), ramkrishna.sharma71@gmail.com
 *   Organization:  University Of Delhi, Delhi, India
 *
 * =====================================================================================
 */

#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <cmath>
#include "TTree.h"
#include "TClonesArray.h"
#include "TLorentzVector.h"

#include "DataFormats/Math/interface/deltaR.h"
#include "Math/GenVector/VectorUtil.h"
#include "DataFormats/GeometryVector/interface/VectorUtil.h"

#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Common/interface/TriggerNames.h"
#include "HLTrigger/HLTcore/interface/HLTConfigProvider.h"

#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include "SimDataFormats/GeneratorProducts/interface/LHEEventProduct.h"
#include "DataFormats/JetReco/interface/GenJet.h"
#include "DataFormats/JetReco/interface/GenJetCollection.h"
#include "DataFormats/METReco/interface/GenMETCollection.h"
#include "DataFormats/METReco/interface/GenMET.h"

#include "TLorentzVector.h"

using namespace edm;
using namespace std;
using namespace reco;
using namespace ROOT::Math::VectorUtil;

class GenAnalyzer : public edm::one::EDAnalyzer<edm::one::SharedResources>  {                                                                  
public:
  explicit GenAnalyzer(const edm::ParameterSet&);
  ~GenAnalyzer();
  
  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);
  static bool reorder(const TLorentzVector &a, const TLorentzVector &b);
  static bool reorder_M(const TLorentzVector &a, const TLorentzVector &b);
  static bool mZ1_mZ2(unsigned int& L1, unsigned int& L2, unsigned int& L3, unsigned int& L4, bool makeCuts, vector<float> GENlep_pt, vector<float> GENlep_eta, vector<float> GENlep_phi, vector<float> GENlep_mass, vector<float> GENlep_id, vector<float> GENlep_RelIso);
  //https://stackoverflow.com/a/26230635/2302094
  static bool jetCleaning(const reco::GenJet  * genAK8jet,const vector<reco::GenJet>* firstJetCollection, const double r_seperation=0.8);
  static void SortedCleanedJetVector(const std::vector<reco::GenJet>* jetCollection1, const std::vector<reco::GenJet>* jetCollection2, const std::vector<TLorentzVector> &photons, std::vector<TLorentzVector> &outputLorentzVector, const double r_seperation=0.8);
  static void SortedCleanedJetVector_(const std::vector<reco::GenJet>* jetCollection1, const std::vector<reco::GenJet>* jetCollection2, const std::vector<TLorentzVector> &leps__, std::vector<TLorentzVector> &outputLorentzVector, const double r_seperation=0.8);
  static void indexOfSelectedJet(const std::vector<TLorentzVector> &inputLorentzVector, int &index1, int &index2);
  static void indexOfSelectedJet(const std::vector<TLorentzVector> &inputLorentzVector, double massComp, int &index1, int &index2, int in_index1=-1, int in_index2=-1);
  static TLorentzVector maxPtLorentzVector(const std::vector<TLorentzVector> &inputLorentzVector);
  static TLorentzVector minMassLorentzVector(const std::vector<TLorentzVector> &inputLorentzVector, const double mass);
  static TLorentzVector minMassLorentzVector(const std::vector<TLorentzVector> &inputLorentzVector, const double mass, int &position, bool skip);
  static void minMassLorentzVector(const std::vector<TLorentzVector> &inputLorentzVector, const double mass, TLorentzVector &leadingJet, TLorentzVector &subleadingJet);
  static std::vector<TLorentzVector> minMassLorentzVector(const std::vector<TLorentzVector> &input_AK4LorentzVector, const std::vector<TLorentzVector> &input_AK8LorentzVector);

  void SetBranches();
  void Clear();
  
private:
  virtual void beginJob() override;
  virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
  virtual void endJob() override;
  
  // ----------member data ---------------------------
  virtual void AddBranch(std::vector<std::string>*, std::string name);
  virtual void AddBranch(std::vector<double>*, std::string name);
  virtual void AddBranch(std::vector<int>*, std::string name);
  virtual void AddBranch(int* vec, std::string name);
  virtual void AddBranch(double* vec, std::string name);
  virtual void computeAngles(TLorentzVector thep4H, TLorentzVector thep4Z1, TLorentzVector thep4M11, TLorentzVector thep4M12, TLorentzVector thep4Z2, TLorentzVector thep4M21, TLorentzVector thep4M22, double& costheta1, double& costheta2, double& Phi, double& costhetastar, double& Phi1);
  
  edm::EDGetTokenT<reco::GenParticleCollection> genParticlesToken;
  edm::EDGetTokenT<LHEEventProduct> LHEEventToken;
  edm::EDGetTokenT<reco::GenJetCollection> genAK4jetToken;
  edm::EDGetTokenT<reco::GenJetCollection> genAK8jetToken;
  edm::EDGetTokenT<reco::GenMETCollection> genMetCaloToken;
  
  bool	Verbose_;
  TString OutPutFileName_;
  
  edm::Service<TFileService> fs;
  TFile * outputFile_;
  TTree* tree;
  std::ofstream file1;
  
  //std::vector<int> pdgID_;
  
  //  std::vector<std::string> LHEWeightIDs_;
  //  std::vector<double> LHEWeights_;
  
  int nEVENT=-999;
  
  //  int		isMuMinus_ = -999;
  //  double 	LHELeptPt_ = -999.0;
  //  double 	LHELeptEta_ = -999.0;
  //  double 	LHELeptPhi_ = -999.0;
  //  double 	LHELeptM_ = -999.0;
  //  double 	LHELeptE_ = -999.0;

  vector<float> gen_Lep_id;
  vector<float> gen_Muon_Pt;
  int gen_Lep_Hindex[4];
  vector<float> gen_Lep_RelIso;
  vector<float> gen_Lep_pt;
  vector<float> gen_Lep_eta;
  vector<float> gen_Lep_phi;
  vector<float> gen_Lep_mass;
  vector<float> gen_Zboson_mass;
  vector<float> gen_Zboson_pt;
  vector<float> gen_Zboson_eta;
  vector<float> gen_Zboson_phi;
  vector<float> gen_Higgs_mass;
  vector<float> genJet_pt;
  vector<float> genJet_eta;
  vector<float> genJet_phi;
  vector<float> genJet_mass;
  vector<float> genMuon_pt;
  vector<float> genMuon_eta;
  vector<float> genMuon_phi;
  vector<float> genMuon_mass;
  float gen_Higgs_mass_FourLep;
  float gen_Higgs_pt_FourLep;
  float gen_Higgs_eta_FourLep;
  float gen_Higgs_phi_FourLep;
  float gen_Higgs_mass_TwoZboson;
  float gen_Higgs_pt_TwoZboson;
  float gen_Higgs_eta_TwoZboson;
  float gen_Higgs_phi_TwoZboson;
  float gen_Zboson1_pt;
  float gen_Zboson1_eta;
  float gen_Zboson1_phi;
  float gen_Zboson1_mass;
  float gen_Zboson2_pt;
  float gen_Zboson2_eta;
  float gen_Zboson2_phi;
  float gen_Zboson2_mass;
  int num_Jet;
  int num_Muon;
  int genLep_nLep;
  
//  
  double genJetAK4_njets_ = -999.0;
  double AK8Gen_HiggsJet_njets_ = -999.0;

};

//
//
// constants, enums and typedefs
//

//
// static data member definitions
//

//
// constructors and destructor
//
GenAnalyzer::GenAnalyzer(const edm::ParameterSet& iConfig)
{
  //now do what ever initialization is needed
  usesResource("TFileService");
  LHEEventToken = consumes<LHEEventProduct>(iConfig.getParameter<edm::InputTag>("LHEEventInputTag"));
  genParticlesToken = consumes<reco::GenParticleCollection>(iConfig.getParameter<edm::InputTag>("genParticlesInputTag"));
  //genAK4jetToken = consumes<reco::GenJetCollection>(iConfig.getParameter<edm::InputTag>("ak4GenJetsInputTag"));
  //genAK8jetToken = consumes<reco::GenJetCollection>(iConfig.getParameter<edm::InputTag>("ak8GenJetsInputTag"));
  Verbose_ = iConfig.getParameter<bool>("Verbose");
  OutPutFileName_ = iConfig.getParameter<std::string>("OutPutFileName");
}


GenAnalyzer::~GenAnalyzer()
{
  // do anything here that needs to be done at desctruction time
  // (e.g. close files, deallocate resources etc.)
  delete outputFile_;
}

void GenAnalyzer::AddBranch(std::vector<std::string>* vec, std::string name){
  tree->Branch(name.c_str(),vec);
}

void GenAnalyzer::AddBranch(std::vector<double>* vec, std::string name){
  tree->Branch(name.c_str(),vec);
}

void GenAnalyzer::AddBranch(std::vector<int>* vec, std::string name){
  tree->Branch(name.c_str(),vec);
}

void GenAnalyzer::AddBranch(int* vec, std::string name){
  tree->Branch(name.c_str(),vec,(name+"/I").c_str());
}

void GenAnalyzer::AddBranch(double* vec, std::string name){
  tree->Branch(name.c_str(),vec,(name+"/D").c_str());
}


//////////////////////////////////
//// P A P E R   4 - V E C T O R   D E F I N I T I O N   O F   P H I   A N D   P H I 1
//////////////////////////////////
void GenAnalyzer::computeAngles(TLorentzVector thep4H, TLorentzVector thep4Z1, TLorentzVector thep4M11, TLorentzVector thep4M12, TLorentzVector thep4Z2, TLorentzVector thep4M21, TLorentzVector thep4M22, double& costheta1, double& costheta2, double& Phi, double& costhetastar, double& Phi1){
  ///////////////////////////////////////////////
  // check for z1/z2 convention, redefine all 4 vectors with convention
  ///////////////////////////////////////////////
  TLorentzVector p4H, p4Z1, p4M11, p4M12, p4Z2, p4M21, p4M22;
  p4H = thep4H;
  
  p4Z1 = thep4Z1; p4M11 = thep4M11; p4M12 = thep4M12;
  p4Z2 = thep4Z2; p4M21 = thep4M21; p4M22 = thep4M22;
  //// costhetastar
  TVector3 boostX = -(thep4H.BoostVector());
  TLorentzVector thep4Z1inXFrame( p4Z1 );
  TLorentzVector thep4Z2inXFrame( p4Z2 );
  thep4Z1inXFrame.Boost( boostX );
  thep4Z2inXFrame.Boost( boostX );
  TVector3 theZ1X_p3 = TVector3( thep4Z1inXFrame.X(), thep4Z1inXFrame.Y(), thep4Z1inXFrame.Z() );
  TVector3 theZ2X_p3 = TVector3( thep4Z2inXFrame.X(), thep4Z2inXFrame.Y(), thep4Z2inXFrame.Z() );
  costhetastar = theZ1X_p3.CosTheta();
  
  //// --------------------------- costheta1
  TVector3 boostV1 = -(thep4Z1.BoostVector());
  TLorentzVector p4M11_BV1( p4M11 );
  TLorentzVector p4M12_BV1( p4M12 );
  TLorentzVector p4M21_BV1( p4M21 );
  TLorentzVector p4M22_BV1( p4M22 );
  p4M11_BV1.Boost( boostV1 );
  p4M12_BV1.Boost( boostV1 );
  p4M21_BV1.Boost( boostV1 );
  p4M22_BV1.Boost( boostV1 );
  
  TLorentzVector p4V2_BV1 = p4M21_BV1 + p4M22_BV1;
  //// costheta1
  costheta1 = -p4V2_BV1.Vect().Dot( p4M11_BV1.Vect() )/p4V2_BV1.Vect().Mag()/p4M11_BV1.Vect().Mag();
  
  //// --------------------------- costheta2
  TVector3 boostV2 = -(thep4Z2.BoostVector());
  TLorentzVector p4M11_BV2( p4M11 );
  TLorentzVector p4M12_BV2( p4M12 );
  TLorentzVector p4M21_BV2( p4M21 );
  TLorentzVector p4M22_BV2( p4M22 );
  p4M11_BV2.Boost( boostV2 );
  p4M12_BV2.Boost( boostV2 );
  p4M21_BV2.Boost( boostV2 );
  p4M22_BV2.Boost( boostV2 );
  
  TLorentzVector p4V1_BV2 = p4M11_BV2 + p4M12_BV2;
  //// costheta2
  costheta2 = -p4V1_BV2.Vect().Dot( p4M21_BV2.Vect() )/p4V1_BV2.Vect().Mag()/p4M21_BV2.Vect().Mag();
  
  //// --------------------------- Phi and Phi1
  //    TVector3 boostX = -(thep4H.BoostVector());
  TLorentzVector p4M11_BX( p4M11 );
  TLorentzVector p4M12_BX( p4M12 );
  TLorentzVector p4M21_BX( p4M21 );
  TLorentzVector p4M22_BX( p4M22 );
  
  p4M11_BX.Boost( boostX );
  p4M12_BX.Boost( boostX );
  p4M21_BX.Boost( boostX );
  p4M22_BX.Boost( boostX );
  
  TVector3 tmp1 = p4M11_BX.Vect().Cross( p4M12_BX.Vect() );
  TVector3 tmp2 = p4M21_BX.Vect().Cross( p4M22_BX.Vect() );
  
  TVector3 normal1_BX( tmp1.X()/tmp1.Mag(), tmp1.Y()/tmp1.Mag(), tmp1.Z()/tmp1.Mag() );
  TVector3 normal2_BX( tmp2.X()/tmp2.Mag(), tmp2.Y()/tmp2.Mag(), tmp2.Z()/tmp2.Mag() );
  
  //// Phi
  TLorentzVector p4Z1_BX = p4M11_BX + p4M12_BX;
  double tmpSgnPhi = p4Z1_BX.Vect().Dot( normal1_BX.Cross( normal2_BX) );
  double sgnPhi = tmpSgnPhi/fabs(tmpSgnPhi);
  Phi = sgnPhi * acos( -1.*normal1_BX.Dot( normal2_BX) );
  
  //////////////
  TVector3 beamAxis(0,0,1);
  TVector3 tmp3 = (p4M11_BX + p4M12_BX).Vect();
  
  TVector3 p3V1_BX( tmp3.X()/tmp3.Mag(), tmp3.Y()/tmp3.Mag(), tmp3.Z()/tmp3.Mag() );
  TVector3 tmp4 = beamAxis.Cross( p3V1_BX );
  TVector3 normalSC_BX( tmp4.X()/tmp4.Mag(), tmp4.Y()/tmp4.Mag(), tmp4.Z()/tmp4.Mag() );
  
  //// Phi1
  double tmpSgnPhi1 = p4Z1_BX.Vect().Dot( normal1_BX.Cross( normalSC_BX) );
  double sgnPhi1 = tmpSgnPhi1/fabs(tmpSgnPhi1);
  Phi1 = sgnPhi1 * acos( normal1_BX.Dot( normalSC_BX) );
  
  //    std::cout << "extractAngles: " << std::endl;
  //    std::cout << "costhetastar = " << costhetastar << ", costheta1 = " << costheta1 << ", costheta2 = " << costheta2 << std::endl;
  //    std::cout << "Phi = " << Phi << ", Phi1 = " << Phi1 << std::endl;
}


void GenAnalyzer::SetBranches(){
  //AddBranch(&pdgID_,	"pdgID");
  //  AddBranch(&isMuMinus_ , "isMuMinus");
  //  AddBranch(&LHELeptPt_ ,	"LHELeptPt");
  //  AddBranch(&LHELeptEta_ ,	"LHELeptEta");
  //  AddBranch(&LHELeptPhi_ ,	"LHELeptPhi");
  //  AddBranch(&LHELeptM_ ,	"LHELeptM");
  //  AddBranch(&LHELeptE_ ,	"LHELeptE");
  //
  //  AddBranch(&LHEWeightIDs_, "LHEWeightIDs");
  //  AddBranch(&LHEWeights_, "LHEWeights");

  //AddBranch(&gen_Lep_pt, "gen_Lep_pt");
  //AddBranch(&gen_Zboson_mass, "gen_Zboson_mass");
  tree->Branch("gen_Lep_id", &gen_Lep_id);
  tree->Branch("gen_Muon_Pt", &gen_Muon_Pt);
  tree->Branch("gen_Lep_Hindex", &gen_Lep_Hindex, "gen_Lep_Hindex[4]/I");
  tree->Branch("gen_Lep_RelIso", &gen_Lep_RelIso);
  tree->Branch("gen_Lep_pt", &gen_Lep_pt);
  tree->Branch("gen_Lep_eta", &gen_Lep_eta);
  tree->Branch("gen_Lep_phi", &gen_Lep_phi);
  tree->Branch("gen_Lep_mass", &gen_Lep_mass);
  tree->Branch("gen_Zboson_mass", &gen_Zboson_mass);
  tree->Branch("gen_Zboson_pt", &gen_Zboson_pt);
  tree->Branch("gen_Zboson_eta", &gen_Zboson_eta);
  tree->Branch("gen_Zboson_phi", &gen_Zboson_phi);
  tree->Branch("gen_Higgs_mass", &gen_Higgs_mass);
  tree->Branch("genJet_pt", &genJet_pt);
  tree->Branch("genJet_eta", &genJet_eta);
  tree->Branch("genJet_phi", &genJet_phi);
  tree->Branch("genJet_mass", &genJet_mass);
  tree->Branch("genMuon_pt", &genMuon_pt);
  tree->Branch("genMuon_eta", &genMuon_eta);
  tree->Branch("genMuon_phi", &genMuon_phi);
  tree->Branch("genMuon_mass", &genMuon_mass);
  tree->Branch("gen_Higgs_mass_FourLep", &gen_Higgs_mass_FourLep);
  tree->Branch("gen_Higgs_pt_FourLep", &gen_Higgs_pt_FourLep);
  tree->Branch("gen_Higgs_eta_FourLep", &gen_Higgs_eta_FourLep);
  tree->Branch("gen_Higgs_phi_FourLep", &gen_Higgs_phi_FourLep);
  tree->Branch("gen_Higgs_mass_TwoZboson", &gen_Higgs_mass_TwoZboson);
  tree->Branch("gen_Higgs_pt_TwoZboson", &gen_Higgs_pt_TwoZboson);
  tree->Branch("gen_Higgs_eta_TwoZboson", &gen_Higgs_eta_TwoZboson);
  tree->Branch("gen_Higgs_phi_TwoZboson", &gen_Higgs_phi_TwoZboson);
  tree->Branch("gen_Zboson1_pt", &gen_Zboson1_pt);
  tree->Branch("gen_Zboson1_eta", &gen_Zboson1_eta);
  tree->Branch("gen_Zboson1_phi", &gen_Zboson1_phi);
  tree->Branch("gen_Zboson1_mass", &gen_Zboson1_mass);
  tree->Branch("gen_Zboson2_pt", &gen_Zboson2_pt);
  tree->Branch("gen_Zboson2_eta", &gen_Zboson2_eta);
  tree->Branch("gen_Zboson2_phi", &gen_Zboson2_phi);
  tree->Branch("gen_Zboson2_mass", &gen_Zboson2_mass);
  tree->Branch("num_Jet", &num_Jet);
  tree->Branch("num_Muon", &num_Muon);
  tree->Branch("genLep_nLep", &genLep_nLep);
  
  AddBranch(&genJetAK4_njets_, "genJetAK4_njets");
  AddBranch(&AK8Gen_HiggsJet_njets_, "genJetAK8_njets");

}

void GenAnalyzer::Clear(){
  //pdgID_.clear();
  //  LHEWeightIDs_.clear();
  //  LHEWeights_.clear();
  gen_Lep_id.clear();
  gen_Muon_Pt.clear();
  //gen_Lep_Hindex[0]=-1;
  //gen_Lep_Hindex[1]=-1;
  //gen_Lep_Hindex[2]=-1;
  //gen_Lep_Hindex[3]=-1;
  gen_Lep_RelIso.clear();
  gen_Lep_pt.clear();
  gen_Lep_eta.clear();
  gen_Lep_phi.clear();
  gen_Lep_mass.clear();
  gen_Zboson_mass.clear();
  gen_Zboson_pt.clear();
  gen_Zboson_eta.clear();
  gen_Zboson_phi.clear();
  gen_Higgs_mass.clear();
  genJet_pt.clear();
  genJet_eta.clear();
  genJet_phi.clear();
  genJet_mass.clear();
  genMuon_pt.clear();
  genMuon_eta.clear();
  genMuon_phi.clear();
  genMuon_mass.clear();
  gen_Higgs_mass_FourLep = 0;
  gen_Higgs_pt_FourLep = 0;
  gen_Higgs_eta_FourLep = 0;
  gen_Higgs_phi_FourLep = 0;
  gen_Higgs_mass_TwoZboson = 0;
  gen_Higgs_pt_TwoZboson = 0;
  gen_Higgs_eta_TwoZboson = 0;
  gen_Higgs_phi_TwoZboson = 0;
  gen_Zboson1_pt = 0;
  gen_Zboson1_eta = 0;
  gen_Zboson1_phi = 0;
  gen_Zboson1_mass = 0;
  gen_Zboson2_pt = 0;
  gen_Zboson2_eta = 0;
  gen_Zboson2_phi = 0;
  gen_Zboson2_mass = 0;
  num_Jet = 0;
  num_Muon = 0;
  genLep_nLep = 0;
  
  genJetAK4_njets_ = -999.0;
  AK8Gen_HiggsJet_njets_ = -999.0;

}

/**
 * This helps to identify the higher pT lorentzVector.
 * @param  a first lorentz vector
 * @param  b second lorentz vector
 * @return   [bool] True if first LV has higher pT than second otherwise False.
 */
bool GenAnalyzer::reorder(const TLorentzVector &a, const TLorentzVector &b)
{
  return a.Pt() > b.Pt();
}

bool GenAnalyzer::reorder_M(const TLorentzVector &a, const TLorentzVector &b)
{
  return a.M() > b.M();
}

bool GenAnalyzer::mZ1_mZ2(unsigned int& L1, unsigned int& L2, unsigned int& L3, unsigned int& L4, bool makeCuts, vector<float> GENlep_pt, vector<float> GENlep_eta, vector<float> GENlep_phi, vector<float> GENlep_mass, vector<float> GENlep_id, vector<float> GENlep_RelIso)
{

    float genIsoCutEl=0.35;
    float genIsoCutMu=0.35;
    float Zmass=91.1876;
    double offshell = 999.0; bool findZ1 = false; bool passZ1 = false;

    L1 = 0; L2 = 0;

    unsigned int N = GENlep_pt.size();

    for(unsigned int i=0; i<N; i++){
        for(unsigned int j=i+1; j<N; j++){


            if((GENlep_id[i]+GENlep_id[j])!=0) continue;

            TLorentzVector li, lj;
            li.SetPtEtaPhiM(GENlep_pt[i],GENlep_eta[i],GENlep_phi[i],GENlep_mass[i]);
            lj.SetPtEtaPhiM(GENlep_pt[j],GENlep_eta[j],GENlep_phi[j],GENlep_mass[j]);

            //if (verbose) cout<<"gen lep i id: "<<GENlep_id[i]<<" pt: "<<li.Pt()<<" lep j id: "<<GENlep_id[j]<<" pt: "<<lj.Pt()<<endl;

            if (makeCuts) {
                if ( abs(GENlep_id[i]) == 13 && (li.Pt() < 5.0 || abs(li.Eta()) > 2.4)) continue;
                if ( abs(GENlep_id[i]) == 11 && (li.Pt() < 7.0 || abs(li.Eta()) > 2.5)) continue;
                if ( GENlep_RelIso[i]>((abs(GENlep_id[i])==11)?genIsoCutEl:genIsoCutMu)) continue;
                
                if ( abs(GENlep_id[j]) == 13 && (lj.Pt() < 5.0 || abs(lj.Eta()) > 2.4)) continue;
                if ( abs(GENlep_id[j]) == 11 && (lj.Pt() < 7.0 || abs(lj.Eta()) > 2.5)) continue;
                if ( GENlep_RelIso[j]>((abs(GENlep_id[i])==11)?genIsoCutEl:genIsoCutMu)) continue;                
            }

            TLorentzVector mll = li+lj;
            //if (verbose) cout<<"gen mass ij: "<<mll.M()<<endl;

            if(abs(mll.M()-Zmass)<offshell){
                double mZ1 = mll.M();
                //if (verbose) cout<<"foundZ1"<<endl;
                L1 = i; L2 = j; findZ1 = true; offshell = abs(mZ1-Zmass);          
            }
        }    
    }

    TLorentzVector l1, l2;
    l1.SetPtEtaPhiM(GENlep_pt[L1],GENlep_eta[L1],GENlep_phi[L1],GENlep_mass[L1]);
    l2.SetPtEtaPhiM(GENlep_pt[L2],GENlep_eta[L2],GENlep_phi[L2],GENlep_mass[L2]);
    TLorentzVector ml1l2 = l1+l2;

    if(ml1l2.M()>40 && ml1l2.M()<120 && findZ1) passZ1 = true;
    if (!makeCuts) passZ1 = true;

    double pTL34 = 0.0; bool findZ2 = false; 
    //bool m4lwindow=false; double window_lo=70.0; double window_hi=140.0;
   
    //cout<<"findZ2"<<endl;
    for(unsigned int i=0; i<N; i++){
        if(i==L1 || i==L2) continue; // can not be the lep from Z1
        for(unsigned int j=i+1; j<N; j++){
            if(j==L1 || j==L2) continue; // can not be the lep from Z1
            if((GENlep_id[i]+GENlep_id[j])!=0) continue;            

            TLorentzVector li, lj;
            li.SetPtEtaPhiM(GENlep_pt[i],GENlep_eta[i],GENlep_phi[i],GENlep_mass[i]);
            lj.SetPtEtaPhiM(GENlep_pt[j],GENlep_eta[j],GENlep_phi[j],GENlep_mass[j]);
            TLorentzVector Z2 = li+lj;

            if (makeCuts) {
                if ( abs(GENlep_id[i]) == 13 && (li.Pt() < 5.0 || abs(li.Eta()) > 2.4)) continue;
                if ( abs(GENlep_id[i]) == 11 && (li.Pt() < 7.0 || abs(li.Eta()) > 2.5)) continue;
                if ( GENlep_RelIso[i]>((abs(GENlep_id[i])==11)?genIsoCutEl:genIsoCutMu)) continue;
                
                if ( abs(GENlep_id[j]) == 13 && (lj.Pt() < 5.0 || abs(lj.Eta()) > 2.4)) continue;
                if ( abs(GENlep_id[j]) == 11 && (lj.Pt() < 7.0 || abs(lj.Eta()) > 2.5)) continue;
                if ( GENlep_RelIso[j]>((abs(GENlep_id[i])==11)?genIsoCutEl:genIsoCutMu)) continue;
            }

            if ( (li.Pt()+lj.Pt())>=pTL34 ) {
                double mZ2 = Z2.M();
                //if (verbose) cout<<"GEN mZ2: "<<mZ2<<endl;
                if( (mZ2>12 && mZ2<120) || (!makeCuts) ) {
                    L3 = i; L4 = j; findZ2 = true; 
                    pTL34 = li.Pt()+lj.Pt();
                    //if (verbose) cout<<"is the new GEN cand"<<endl;
                    //if (m4l>window_lo && m4l<window_hi) m4lwindow=true;
                } else {
                    // still assign L3 and L4 to this pair if we don't have a passing Z2 yet
                    if (findZ2 == false) {L3 = i; L4 = j;}
                    //cout<<"is not new GEN cand"<<endl;
                }
            }
            
        } // lj
    } // li

    if(passZ1 && findZ2) return true;
    else return false;
    
}


/**
 * This takes one AK4 (AK8) jet and checks from AK8 (AK4) jet collection if its cleaned or not.
 * @param  genAK8jet          [reco::GenJet] one gen jet
 * @param  firstJetCollection [const vector<reco::GenJet>] vector of genJet collection 
 * @param  r_seperation       [double] this is the value of delta R     
 * @return                    [bool] returns true if jet is cleaned else false.
 */
bool GenAnalyzer::jetCleaning(const reco::GenJet  * genAK8jet,const vector<reco::GenJet>* firstJetCollection, const double r_seperation)
{
  for(vector<reco::GenJet>::const_iterator genjet = firstJetCollection->begin(); genjet != firstJetCollection->end(); genjet++) {
    if (deltaR(genAK8jet->pt(), genAK8jet->eta(), genjet->pt(), genjet->eta()) < r_seperation) return false;
  }
  return true;
}

/**
 * This function returns the sorted TLorentzVector that contains information from first
 * passed jetCollection.
 * @param firstJetCollection         first genjetCollection whose information will be passed to 
 *                            TLorentzVector local_Vec_genJetAK5
 * @param secondJetCollection         Another genjetCollection which from which the first jet
 *                            collection is going to be cleaned.
 * @param Vec_Photons         Sorted TLorentzVector that contains the photon collection.
 * @param local_Vec_genJetAK4 vector of TLorentzVector that we need
 * @param r_seperation        This is the delta-R seperation distance beetween the two 
 *                            jets.
 */
void GenAnalyzer::SortedCleanedJetVector(const std::vector<reco::GenJet>* firstJetCollection, const std::vector<reco::GenJet>* secondJetCollection, const std::vector<TLorentzVector> &Vec_Photons, std::vector<TLorentzVector> &local_Vec_genJetAK4, const double r_seperation)
{
  // int nAK4jets = 0;
  TLorentzVector temp_genJetAK4;
  for(vector<reco::GenJet>::const_iterator genjet = firstJetCollection->begin(); genjet != firstJetCollection->end(); genjet++) {
    if (deltaR(genjet->eta(),genjet->phi(), Vec_Photons[0].Eta(),Vec_Photons[0].Phi())>0.4 && deltaR(genjet->eta(),genjet->phi(), Vec_Photons[1].Eta(),Vec_Photons[1].Phi())>0.4) {
      if ( GenAnalyzer::jetCleaning(&(*genjet), secondJetCollection, r_seperation ) ) {
        temp_genJetAK4.SetPtEtaPhiE(genjet->pt(), genjet->eta(), genjet->phi(), genjet->energy());
        local_Vec_genJetAK4.push_back(temp_genJetAK4);
        // std::cout << "**> " << genjet->pt() << "\t" << genjet->eta() << "\t" << genjet->phi() << "\t" << genjet->energy()<< std::endl;
      }
    }
  }
  std::sort(local_Vec_genJetAK4.begin(), local_Vec_genJetAK4.end(), GenAnalyzer::reorder);
}

void GenAnalyzer::SortedCleanedJetVector_(const std::vector<reco::GenJet>* firstJetCollection, const std::vector<reco::GenJet>* secondJetCollection, const std::vector<TLorentzVector> &Vec_Lep, std::vector<TLorentzVector> &local_Vec_genJetAK4, const double r_seperation)
{
  // int nAK4jets = 0;
  TLorentzVector temp_genJetAK4;
  for(vector<reco::GenJet>::const_iterator genjet = firstJetCollection->begin(); genjet != firstJetCollection->end(); genjet++) {
    double min_deltaR=999;
    for(int ii=0; (unsigned)ii<Vec_Lep.size(); ii++)
    {
      double temp_deltaR = deltaR(genjet->eta(),genjet->phi(), Vec_Lep[ii].Eta(), Vec_Lep[ii].Phi());
        if (temp_deltaR<min_deltaR)
          min_deltaR = temp_deltaR;

    }
    if (min_deltaR>0.4){
    //if (deltaR(genjet->eta(),genjet->phi(), Vec_Photons[0].Eta(),Vec_Photons[0].Phi())>0.4 && deltaR(genjet->eta(),genjet->phi(), Vec_Photons[1].Eta(),Vec_Photons[1].Phi())>0.4) {
      if ( GenAnalyzer::jetCleaning(&(*genjet), secondJetCollection, r_seperation ) ) {
        temp_genJetAK4.SetPtEtaPhiE(genjet->pt(), genjet->eta(), genjet->phi(), genjet->energy());
        local_Vec_genJetAK4.push_back(temp_genJetAK4);
        // std::cout << "**> " << genjet->pt() << "\t" << genjet->eta() << "\t" << genjet->phi() << "\t" << genjet->energy()<< std::endl;
      }
    }
  }
  std::sort(local_Vec_genJetAK4.begin(), local_Vec_genJetAK4.end(), GenAnalyzer::reorder);
}

/**
 * This takes a vector of TLorentzVector and give a TLorentzVector having highest pT
 * @param  inputLorentzVector Vector of TLorentzVector
 * @return                    TLorentzVector having highest pT
 */
TLorentzVector GenAnalyzer::maxPtLorentzVector(const std::vector<TLorentzVector> &inputLorentzVector) {
  double temp_pt = -999.0;
  TLorentzVector AK8Gen_HiggsJet_MaxPt;
  for (std::vector<TLorentzVector>::const_iterator jet = inputLorentzVector.begin(); jet != inputLorentzVector.end(); ++jet)
  {
    if (jet->Pt()>temp_pt)
    {
          AK8Gen_HiggsJet_MaxPt.SetPtEtaPhiE(jet->Pt(), jet->Eta(), jet->Phi(), jet->Energy());
          temp_pt = jet->Pt();      
    }
  }
  return AK8Gen_HiggsJet_MaxPt;
}

/**
 * This takes a vector of TLorentzVector and gives a TLorentzVector having minimum mass difference from the specified mass.
 * @param  inputLorentzVector [TLorentzVector] Input vector of TLorentzVector
 * @param  mass               [double] mass from which it will compare.
 * @return                    TLorentzVector having minimum mass difference w.r.t. specified mass.
 */
TLorentzVector GenAnalyzer::minMassLorentzVector(const std::vector<TLorentzVector> &inputLorentzVector, const double mass) {
  double temp_AK8jet_deltaM = 9999.0;
  TLorentzVector AK8Gen_HiggsJet_MaxPt;
  for (std::vector<TLorentzVector>::const_iterator jet = inputLorentzVector.begin(); jet != inputLorentzVector.end(); ++jet)
  {
    if (abs(jet->M()-mass)<temp_AK8jet_deltaM)
    {
          AK8Gen_HiggsJet_MaxPt.SetPtEtaPhiE(jet->Pt(), jet->Eta(), jet->Phi(), jet->Energy());
          temp_AK8jet_deltaM = abs(jet->M()-mass);
    }
  }
  return AK8Gen_HiggsJet_MaxPt;
}

/**
 * This takes a vector of TLorentzVector and returns a TLorentzVector having minimum delta mass w.r.t. the 
 * specified mass. Also, it takes a input positional arguments. It just skip that while checking.
 * @param  inputLorentzVector Input vector of TLorentzVector.
 * @param  mass               mass from which it will compare
 * @param  position           position of the vector that need to be skip from calculation
 * @param  skip               [description]
 * @return                    TLorentzVector having minimum mass difference w.r.t. specified mass.
 */
TLorentzVector GenAnalyzer::minMassLorentzVector(const std::vector<TLorentzVector> &inputLorentzVector, const double mass, int &position, bool skip) {
  double temp_AK8jet_deltaM = 9999.0;
  TLorentzVector AK8Gen_HiggsJet_MaxPt;
  int counter = -1;
  // std::cout << "===================" << position << std::endl;
  for (std::vector<TLorentzVector>::const_iterator jet = inputLorentzVector.begin(); jet != inputLorentzVector.end(); ++jet)
  {
    counter++;
    // std::cout << "==> " << counter << "\t" << position << std::endl;
    if (skip) 
    {
      if (counter == position) continue;
    }
    // std::cout << counter << "\t" << position << "\t" << jet->Pt() <<  "\t" << abs(jet->M()-mass) << "\t" << temp_AK8jet_deltaM<< std::endl;
    if (abs(jet->M()-mass)<temp_AK8jet_deltaM)
    {
      AK8Gen_HiggsJet_MaxPt.SetPtEtaPhiE(jet->Pt(), jet->Eta(), jet->Phi(), jet->Energy());
      temp_AK8jet_deltaM = abs(jet->M()-mass);
      if (!skip) 
      {
        position = counter;
        // std::cout << "position : " << position << std::endl;
      }
    }
  }
  return AK8Gen_HiggsJet_MaxPt;
}

/**
 * [GenAnalyzer::minMassLorentzVector description]
 * @param inputLorentzVector [description]
 * @param mass               [description]
 * @param leadingJet         [description]
 * @param subleadingJet      [description]
 */
void GenAnalyzer::minMassLorentzVector(const std::vector<TLorentzVector> &inputLorentzVector, const double mass, TLorentzVector &leadingJet, TLorentzVector &subleadingJet) {
  double temp_AK8jet_deltaM = 9999.0;
  for (std::vector<TLorentzVector>::const_iterator jet = inputLorentzVector.begin(); jet != inputLorentzVector.end()-1; ++jet)
    for (std::vector<TLorentzVector>::const_iterator jet1 = jet+1; jet1 != inputLorentzVector.end(); ++jet1)
    {
      if (abs((*jet+*jet1).M()-mass)<temp_AK8jet_deltaM)
      {
        temp_AK8jet_deltaM = abs((*jet+*jet1).M() - mass);
        if (jet->Pt()>jet1->Pt()) {
          leadingJet.SetPtEtaPhiE((jet)->Pt(), (jet)->Eta(), (jet)->Phi(), (jet)->Energy());          
          subleadingJet.SetPtEtaPhiE((jet1)->Pt(), (jet1)->Eta(), (jet1)->Phi(), (jet1)->Energy());          
        } else {
          leadingJet.SetPtEtaPhiE((jet1)->Pt(), (jet1)->Eta(), (jet1)->Phi(), (jet1)->Energy());   
          subleadingJet.SetPtEtaPhiE((jet)->Pt(), (jet)->Eta(), (jet)->Phi(), (jet)->Energy());          
        }
      }
    }
}


/**
 * Returns the two index from jet collection that satisfies condition abs(mass-80) minimum.
 * @param inputLorentzVector Input TLorentzVector that contains AK8 jet information
 * @param index1             Index of first on-shell selected jet
 * @param index2             Index of second on-shell selected jet
 */
void GenAnalyzer::indexOfSelectedJet(const std::vector<TLorentzVector> &inputLorentzVector, int &index1, int &index2) {
  double tempMass1 = 9999.0;
  for ( int ak4_jet1 = 0; ak4_jet1 < int(inputLorentzVector.size())-1; ++ak4_jet1)
  {
    for ( int ak4_jet2 = ak4_jet1+1; ak4_jet2 < int(inputLorentzVector.size()); ++ak4_jet2)
    {
      double mass = (inputLorentzVector[ak4_jet1] + inputLorentzVector[ak4_jet2]).M();      
      if (abs(mass - 80.0) < tempMass1)
      {
        tempMass1 = abs(mass - 80.0);
        index1 = ak4_jet1;
        index2 = ak4_jet2;
      }
    }
  }
}

/**
 * Return the two index from the jet collection vector of type TLorentzVector.
 * @param inputLorentzVector Input TLorentzVector that contains jet information
 * @param massComp           Value of mass for which the minimum delta mass is going to be calculated.
 * @param index1             Index of first selected jet.
 * @param index2             Index of second selected jet.
 * @param in_index1          Index of fist jet to be ignored.
 * @param in_index2          Index of second jet to be ignored.
 */
void GenAnalyzer::indexOfSelectedJet(const std::vector<TLorentzVector> &inputLorentzVector, double massComp, int &index1, int &index2, int in_index1, int in_index2) {
  double tempMass1 = 9999.0;
  for ( int ak4_jet1 = 0; ak4_jet1 < int(inputLorentzVector.size())-1; ++ak4_jet1)
  {
    double mass = 0.0;
    if (ak4_jet1 == in_index1) continue;
    if (ak4_jet1 == in_index2) continue;
    for ( int ak4_jet2 = ak4_jet1+1; ak4_jet2 < int(inputLorentzVector.size()); ++ak4_jet2)
    {
      if (ak4_jet2 == in_index1) continue;
      if (ak4_jet2 == in_index2) continue;
      if (in_index1 == -1 && in_index2 == -1) {
          mass = (inputLorentzVector[ak4_jet1] + inputLorentzVector[ak4_jet2]).M();      
      } else {
          mass = (inputLorentzVector[ak4_jet1] + inputLorentzVector[ak4_jet2] + 
          inputLorentzVector[in_index1] + inputLorentzVector[in_index2]).M();
      }
      if (abs(mass - massComp) < tempMass1)
      {
        tempMass1 = abs(mass - massComp);
        index1 = ak4_jet1;
        index2 = ak4_jet2;
      }
    }
  }
}


/**
 * @brief      Selection of two AK4 jet and one AK8 jet. This selection is done
 *             by searching the combination of two AK4 and one AK8 which gives 
 *             us the minimum delta mass w.r.t. the Higgs mass.
 *
 * @param      input_AK4LorentzVector  The input vector of AK4 lorentz vector
 * @param      input_AK8LorentzVector  The input vector of AK8 lorentz vector
 *
 * @return     Returns the vector of TLorentzVector having size 3. Whose first two elements
 *             contains information about the AK4 jets and the third element contains 
 *             information of AK8 jet.
 */
std::vector<TLorentzVector> GenAnalyzer::minMassLorentzVector(const std::vector<TLorentzVector> &input_AK4LorentzVector, const std::vector<TLorentzVector> &input_AK8LorentzVector)
{
  std::vector<TLorentzVector> outputLorentzVector;
  double tempMass1 = 9999.0;

  for (std::vector<TLorentzVector>::const_iterator i = input_AK4LorentzVector.begin(); i != input_AK4LorentzVector.end()-1; ++i)
  {
    for (std::vector<TLorentzVector>::const_iterator j = i+1; j != input_AK4LorentzVector.end(); ++j)
    {
      for (std::vector<TLorentzVector>::const_iterator k = input_AK8LorentzVector.begin(); k != input_AK8LorentzVector.end(); ++k)
      {
        double mass = (*i + *j + *k).M();
        if (abs(mass - 125.0) < tempMass1)
        {
          outputLorentzVector.clear();
          outputLorentzVector.push_back(*i);
          outputLorentzVector.push_back(*j);
          outputLorentzVector.push_back(*k);
        }
      }      
    }
  }
  if (outputLorentzVector.size()>3) 
  {
    std::cout << "size of output vector seems more than three... please check the code." << std::endl;
    exit(EXIT_FAILURE);
  }
  return outputLorentzVector;
}
