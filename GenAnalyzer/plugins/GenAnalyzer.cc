// -*- C++ -*-
//
// Package:    GenAnalyzerNew/GenAnalyzer
// Class:      GenAnalyzer
// 
/**\class GenAnalyzer GenAnalyzer.cc GenAnalyzerNew/GenAnalyzer/plugins/GenAnalyzer.cc
 
 Description: Read the GENSIM file
 
 Implementation:
 [Notes on implementation]
 */
//
// Original Author:  Ram Krishna Sharma
//         Created:  Wed, 21 Sep 2016 22:13:15 GMT
//
//


// system include files
#include <memory>
#include <typeinfo>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/one/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "SimDataFormats/GeneratorProducts/interface/LHEEventProduct.h"
#include "SimDataFormats/GeneratorProducts/interface/GenEventInfoProduct.h"
#include "GenAnalyzer.h"

#include "DataFormats/JetReco/interface/GenJet.h"
#include "DataFormats/JetReco/interface/GenJetCollection.h"
#include "DataFormats/METReco/interface/GenMETCollection.h"
#include "DataFormats/METReco/interface/GenMET.h"

#include "DataFormats/Math/interface/deltaR.h"

#include <vector>
#include "TMath.h"
#include "TLorentzVector.h"
//
// class declaration
//

// If the analyzer does not use TFileService, please remove
// the template argument to the base class so the class inherits
// from  edm::one::EDAnalyzer<> and also remove the line from
// constructor "usesResource("TFileService");"
// This will improve performance in multithreaded jobs.

using namespace std;
using namespace reco;
using namespace edm;
//

// ------------ method called for each event  ------------
void
GenAnalyzer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  using namespace edm;
  
  edm::Handle<reco::GenParticleCollection> genpsHandle;
  iEvent.getByToken(genParticlesToken, genpsHandle);
  
  //edm::Handle<reco::GenJetCollection> genAK4jetHandle;
  //iEvent.getByToken(genAK4jetToken, genAK4jetHandle);
  //
  //edm::Handle<reco::GenJetCollection> genAK8jetHandle;
  //iEvent.getByToken(genAK8jetToken, genAK8jetHandle);
  
  const vector<reco::GenParticle>* genps_coll = genpsHandle.product();
  //const vector<reco::GenJet>* genAK4_coll = genAK4jetHandle.product();
  //const vector<reco::GenJet>* genAK8_coll = genAK8jetHandle.product();

  edm::Handle<LHEEventProduct> LHEEventHandle;
  iEvent.getByToken(LHEEventToken, LHEEventHandle);
  //need for other
  const LHEEventProduct* LHE = 0;
  
  //if (Verbose_)  std::cout<<"===================\n\n"<<endl;
  
  int nGenParticle=0;
  int nGenLep=0;
  int nGenZ=0;
  int nGenH=0;
  int nGenJet=0;
  
  TLorentzVector PHO;
  TLorentzVector Wpquarks;
  TLorentzVector Wmquarks;
  TLorentzVector Wboson1;
  TLorentzVector Wboson2;
  TLorentzVector Higgs1;
  TLorentzVector Higgs2;
  TLorentzVector Lep;
  TLorentzVector Zboson1;
  TLorentzVector Zboson2;
  TLorentzVector Jet;
  TLorentzVector Muon;
  
  std::vector<TLorentzVector> Vec_Photons;
  std::vector<TLorentzVector> Vec_wpJET;
  std::vector<TLorentzVector> Vec_wmJET;
  std::vector<TLorentzVector> Vec_Wboson;
  std::vector<TLorentzVector> Vec_Higgs;
  std::vector<TLorentzVector> Vec_Zboson;
  std::vector<TLorentzVector> Vec_Lep;
  std::vector<TLorentzVector> Vec_Jet;
  std::vector<TLorentzVector> Vec_Muon;
  std::vector<float> Vec_Lep_id;
  //std::vector<float> Vec_Muon_Pt;
  
  if (Verbose_) std::cout<<"Size = "<<genps_coll->size()<<std::endl;
 
  //std::cout<<"<event>"<<std::endl; 
  for(vector<reco::GenParticle>::const_iterator genps_it = genps_coll->begin(); genps_it != genps_coll->end(); genps_it++)
  {
    nGenParticle++;
    if (Verbose_) std::cout<<"nGenParticle: "<<nGenParticle<<std::endl;
    //if ( ( abs(genps_it->pdgId())==13 ) && abs(genps_it->mother()->pdgId())==23 && genps_it->isHardProcess() && (genps_it->status()==23 || genps_it->status()==1) )
    //if ( ( abs(genps_it->pdgId())==13 ) && abs(genps_it->mother()->pdgId())==23 && genps_it->isHardProcess() && (genps_it->status()==23 || genps_it->status()==1) )
    if ( ( abs(genps_it->pdgId())==13 ) && (genps_it->status()==23 || genps_it->status()==1) )
    {
        gen_Muon_Pt.push_back(genps_it->pt());
        Muon.SetPtEtaPhiE(genps_it->pt(), genps_it->eta(), genps_it->phi(), genps_it->energy());
        Vec_Muon.push_back(Muon);

    }
    //if ( ( abs(genps_it->pdgId())==15 ) && abs(genps_it->mother()->pdgId())==23 && genps_it->isHardProcess() && (genps_it->status()==23 || genps_it->status()==1) )
    //{
    //    if (Verbose_) std::cout<<"15 Tau and its mother is Z"<<std::endl;
    ////if ( ( abs(genps_it->pdgId())==11 || abs(genps_it->pdgId())==13 ) && abs(genps_it->mother()->pdgId())==15 && (genps_it->mother()->mother()->pdgId()==23 && genps_it->mother()->mother()->mother()->pdgId()==25 ) && genps_it->isHardProcess() && (genps_it->status()==23 || genps_it->status()==1) )///false isHardProcess
    ////    std::cout<<"e or muon -> Tau -> Z -> Higgs"<<std::endl;
    //}
    if (Verbose_) std::cout<<"aaaaaaaaa"<<std::endl;
//////    if ( ( abs(genps_it->pdgId())==11 || abs(genps_it->pdgId())==13 ) && genps_it->mother()->mother()->pdgId()==23 && genps_it->mother()->mother()->mother()->pdgId()==25 && abs(genps_it->mother()->pdgId())==15 && (genps_it->status()==23 || genps_it->status()==1))
//////    {
//////        if(Verbose_) std::cout<<"e or muon -> Tau -> Z -> Higgs"<<std::endl;
//////        if (genps_it->isHardProcess()) std::cout<<"isHardProcess"<<std::endl;///false isHardProcess
//////        //if (genps_it->mother()->isHardProcess()) std::cout<<"e or muon -> Tau isHardProcess"<<std::endl;///no function of mother
////////        nGenLep++;
////////        Lep.SetPtEtaPhiE(genps_it->mother()->pt(), genps_it->mother()->eta(), genps_it->mother()->phi(), genps_it->mother()->energy());
////////        Vec_Lep.push_back(Lep);
//////    }
    if (Verbose_) std::cout<<"bbbbbbb"<<std::endl;

    //if ( ( abs(genps_it->pdgId())==11 || abs(genps_it->pdgId())==13 ) && abs(genps_it->mother()->pdgId())==23 && genps_it->isHardProcess() && (genps_it->status()==23 || genps_it->status()==1) ) 
    ////if ( ( abs(genps_it->pdgId())==13 ) && abs(genps_it->mother()->pdgId())==23 && genps_it->isHardProcess() && (genps_it->status()==23 || genps_it->status()==1) ) 
    //{
    //    if (Verbose_) std::cout<<"in ele or mu and its mother is Z"<<std::endl;
    //    nGenLep++;
    //    Lep.SetPtEtaPhiE(genps_it->pt(), genps_it->eta(), genps_it->phi(), genps_it->energy());
    //    Vec_Lep.push_back(Lep);
    //    Vec_Lep_id.push_back(genps_it->pdgId());

    //}
    /////////new
    //if ( ( abs(genps_it->pdgId())==11 || abs(genps_it->pdgId())==13 ) && abs(genps_it->mother()->pdgId())==23 && genps_it->isHardProcess() && (genps_it->status()==23 || genps_it->status()==1) ) 
    //{
    //    //std::cout<<"    "<<genps_it->pdgId()<<"    "<<genps_it->px()<<"    "<<genps_it->py()<<"    "<<genps_it->pz()<<"    "<<genps_it->energy()<<"    "<<std::endl;
    //}
    ///////////new
    //if ( ( abs(genps_it->pdgId())==23 ) && abs(genps_it->mother()->pdgId())==25 && genps_it->isHardProcess() ) 
    //{
    //    //std::cout<<"    "<<genps_it->pdgId()<<"    "<<genps_it->px()<<"    "<<genps_it->py()<<"    "<<genps_it->pz()<<"    "<<genps_it->energy()<<"    "<<genps_it->mass()<<"    "<<std::endl;
    //    //std::cout<<"    "<<genps_it->pdgId()<<"    "<<genps_it->px()<<"    "<<genps_it->py()<<"    "<<genps_it->pz()<<"    "<<genps_it->energy()<<"    "<<std::endl;
    //}
    /////////new
    //if ( (abs(genps_it->pdgId())<7 || abs(genps_it->pdgId())==21) && (genps_it->status()==23 || genps_it->status()==1 ) && genps_it->isHardProcess() )
    //if ( (abs(genps_it->pdgId())<7 || abs(genps_it->pdgId())==21) && (genps_it->status()==23 || genps_it->status()==1 ))
    //{
    //    //std::cout<<"    "<<genps_it->pdgId()<<"    "<<genps_it->px()<<"    "<<genps_it->py()<<"    "<<genps_it->pz()<<"    "<<genps_it->energy()<<"    "<<std::endl;
    //}
    /////////end new
    if (Verbose_) std::cout<<"cccccc"<<std::endl;
    if (Verbose_) std::cout<<"num_lep: "<<nGenLep<<std::endl;
    //if ( ( abs(genps_it->pdgId())==23 ) && (abs(genps_it->mother()->pdgId())==25) && genps_it->isHardProcess() )
    //{
    //    if (Verbose_) std::cout<<"in Z and its mother is H"<<std::endl;
    //    nGenZ++;
    //    Zboson1.SetPtEtaPhiE(genps_it->pt(), genps_it->eta(), genps_it->phi(), genps_it->energy());
    //    Vec_Zboson.push_back(Zboson1);
    //}
    if (Verbose_) std::cout<<"dddddd"<<std::endl;
    if (Verbose_) std::cout<<"num_Z: "<<nGenZ<<std::endl;
    //if ( genps_it->pdgId()==25 && genps_it->isHardProcess() )
    //if ( genps_it->pdgId()==25 && genps_it->status()==22 )
    //{
    //    nGenH++;
    //    Higgs1.SetPtEtaPhiE(genps_it->pt(), genps_it->eta(), genps_it->phi(), genps_it->energy());
    //    Vec_Higgs.push_back(Higgs1);
    //    //std::cout<<"    25    "<<Higgs1.Px()<<"    "<<Higgs1.Py()<<"    "<<Higgs1.Pz()<<"    "<<Higgs1.E()<<std::endl;
    //    //std::cout<<"    25    "<<Higgs1.Pt()<<std::endl;
    //}
    if (Verbose_) std::cout<<"eeeee"<<std::endl;
    if (Verbose_) std::cout<<"num_H: "<<nGenH<<std::endl;
    if ( (abs(genps_it->pdgId())<7 || abs(genps_it->pdgId())==21) && (genps_it->status()==23 || genps_it->status()==1 ) && genps_it->isHardProcess() )
    {
        nGenJet++;
        Jet.SetPtEtaPhiE(genps_it->pt(), genps_it->eta(), genps_it->phi(), genps_it->energy());
        Vec_Jet.push_back(Jet);
    }
    if (Verbose_) std::cout<<"ffffffff"<<std::endl;
    if (Verbose_) std::cout<<"num_Jet: "<<nGenJet<<std::endl;
    //if ( abs(genps_it->pdgId())==21 && (genps_it->status()==23 ) && genps_it->isHardProcess() )
    //    if (Verbose_) std::cout<<"gluon "<<std::endl;

    if (Verbose_) std::cout<<"ffffffff"<<std::endl;

  }
  genLep_nLep = Vec_Lep.size();
  num_Jet=Vec_Jet.size();
  if (Verbose_) cout<<"count again: "<<genLep_nLep<<"; :"<<num_Jet<<std::endl;
  num_Muon=Vec_Muon.size();
  if (Verbose_) cout<<"count again: "<<num_Jet<<std::endl;

  if (Verbose_) std::cout<<"hhhhhhhh"<<std::endl;
  if (Verbose_) std::cout<<"num_GenParticle: "<<nGenParticle<<std::endl;


  //if (Vec_Lep.size()==4 && Vec_Zboson.size()==2 && Vec_Higgs.size()==1)
  //{
  //    TLorentzVector FourLep, TwoZboson;
  //    //std::sort(Vec_Lep.begin(), Vec_Lep.end(), GenAnalyzer::reorder);
  //    std::sort(Vec_Zboson.begin(), Vec_Zboson.end(), GenAnalyzer::reorder_M);
  //    for (int i=0; (unsigned)i<Vec_Lep.size() ; i++)
  //    {
  //        gen_Lep_pt.push_back(Vec_Lep[i].Pt());
  //        gen_Lep_eta.push_back(Vec_Lep[i].Eta());
  //        gen_Lep_phi.push_back(Vec_Lep[i].Phi());
  //        gen_Lep_mass.push_back(Vec_Lep[i].M());
  //        gen_Lep_id.push_back(Vec_Lep_id[i]);
  //        FourLep+=Vec_Lep[i];

  //        TLorentzVector thisLep;
  //        thisLep.SetPtEtaPhiM(Vec_Lep[i].Pt(),Vec_Lep[i].Eta(),Vec_Lep[i].Phi(),Vec_Lep[i].M());
  //        // GEN iso calculation
  //        double this_GENiso=0.0;
  //        double this_GENneutraliso=0.0;
  //        double this_GENchargediso=0.0;
  //        if (Verbose_) cout<<"gen iso calculation"<<endl;
  //        for(size_t k=0; k<genps_coll->size();k++){
  //        //for(vector<reco::GenParticle>::const_iterator genps_coll = genps_coll->begin(); genps_coll != genps_coll->end(); genps_coll++)
  //            if( (*genps_coll)[k].status() != 1 ) continue; // stable particles only         
  //            if (abs((*genps_coll)[k].pdgId())==12 || abs((*genps_coll)[k].pdgId())==14 || abs((*genps_coll)[k].pdgId())==16) continue; // exclude neutrinos
  //            if ((abs((*genps_coll)[k].pdgId())==11 || abs((*genps_coll)[k].pdgId())==13)) continue; // exclude leptons
  //            //if (gen_fsrset.find(k)!=gen_fsrset.end()) continue; // exclude particles which were selected as fsr photons
  //            double this_dRvL = deltaR(thisLep.Eta(), thisLep.Phi(), (*genps_coll)[k].eta(), (*genps_coll)[k].phi());
  //            if(this_dRvL<0.3) {
  //                if (Verbose_) cout<<"adding to geniso id: "<<(*genps_coll)[k].pdgId()<<" status: "<<(*genps_coll)[k].status()<<" pt: "<<(*genps_coll)[k].pt()<<" dR: "<<this_dRvL<<endl;
  //                this_GENiso = this_GENiso + (*genps_coll)[k].pt();
  //                if ((*genps_coll)[k].charge()==0) this_GENneutraliso = this_GENneutraliso + (*genps_coll)[k].pt();
  //                if ((*genps_coll)[k].charge()!=0) this_GENchargediso = this_GENchargediso + (*genps_coll)[k].pt();
  //            }
  //        } // GEN iso loop
  //        this_GENiso = this_GENiso/thisLep.Pt();
  //        if (Verbose_) cout<<"gen lep pt: "<<thisLep.Pt()<<" rel iso: "<<this_GENiso<<endl;
  //        gen_Lep_RelIso.push_back(this_GENiso);
  //        // END GEN iso calculation

  //    }
  //    
  //    unsigned int L1=99; unsigned int L2=99; unsigned int L3=99; unsigned int L4=99;
  //    GenAnalyzer::mZ1_mZ2(L1, L2, L3, L4, false, gen_Lep_pt, gen_Lep_eta, gen_Lep_phi, gen_Lep_mass, gen_Lep_id, gen_Lep_RelIso);
  //    gen_Lep_Hindex[0]=L1;
  //    gen_Lep_Hindex[1]=L2;
  //    gen_Lep_Hindex[2]=L3;
  //    gen_Lep_Hindex[3]=L4;
  //    gen_Higgs_mass_FourLep=FourLep.M();
  //    gen_Higgs_pt_FourLep=FourLep.Pt();
  //    gen_Higgs_eta_FourLep=FourLep.Eta();
  //    gen_Higgs_phi_FourLep=FourLep.Phi();
  //    for (Int_t i = 0; (unsigned)i<Vec_Zboson.size(); i++)
  //    {
  //        gen_Zboson_mass.push_back(Vec_Zboson[i].M());
  //        gen_Zboson_pt.push_back(Vec_Zboson[i].Pt());
  //        gen_Zboson_eta.push_back(Vec_Zboson[i].Eta());
  //        gen_Zboson_phi.push_back(Vec_Zboson[i].Phi());
  //        TwoZboson+=Vec_Zboson[i];
  //    }
  //    gen_Higgs_mass_TwoZboson=TwoZboson.M();
  //    gen_Higgs_pt_TwoZboson=TwoZboson.Pt();
  //    gen_Higgs_eta_TwoZboson=TwoZboson.Eta();
  //    gen_Higgs_phi_TwoZboson=TwoZboson.Phi();
  //    gen_Zboson1_pt=Vec_Zboson[0].Pt();
  //    gen_Zboson1_eta=Vec_Zboson[0].Eta();
  //    gen_Zboson1_phi=Vec_Zboson[0].Phi();
  //    gen_Zboson1_mass=Vec_Zboson[0].M();
  //    gen_Zboson2_pt=Vec_Zboson[1].Pt();
  //    gen_Zboson2_eta=Vec_Zboson[1].Eta();
  //    gen_Zboson2_phi=Vec_Zboson[1].Phi();
  //    gen_Zboson2_mass=Vec_Zboson[1].M();
  //    //for (Int_t i = 0; (unsigned)i<Vec_Higgs.size(); i++)
  //    //{
  //    //    //std::cout<<Vec_Higgs[i].M()<<endl;
  //    //    gen_Higgs_mass.push_back(Vec_Higgs[i].M());
  //    //}
  //    
  //    std::sort(Vec_Jet.begin(), Vec_Jet.end(), GenAnalyzer::reorder);
  //    for(int i=0; (unsigned)i<Vec_Jet.size(); i++)
  //    {
  //        genJet_pt.push_back(Vec_Jet[i].Pt());
  //        genJet_eta.push_back(Vec_Jet[i].Eta());
  //        genJet_phi.push_back(Vec_Jet[i].Phi());
  //        genJet_mass.push_back(Vec_Jet[i].M());
  //    }

  //}
  std::sort(Vec_Jet.begin(), Vec_Jet.end(), GenAnalyzer::reorder);
  for(int i=0; (unsigned)i<Vec_Jet.size(); i++)
  {
      genJet_pt.push_back(Vec_Jet[i].Pt());
      genJet_eta.push_back(Vec_Jet[i].Eta());
      genJet_phi.push_back(Vec_Jet[i].Phi());
      genJet_mass.push_back(Vec_Jet[i].M());
  }

  std::sort(Vec_Muon.begin(), Vec_Muon.end(), GenAnalyzer::reorder);
  for(int i=0; (unsigned)i<Vec_Muon.size(); i++)
  {
      genMuon_pt.push_back(Vec_Muon[i].Pt());
      genMuon_eta.push_back(Vec_Muon[i].Eta());
      genMuon_phi.push_back(Vec_Muon[i].Phi());
      genMuon_mass.push_back(Vec_Muon[i].M());
  }


  
  tree->Fill();
  //std::cout<<"<rwgt>"<<std::endl;
  //std::cout<< "<wgt id='SM'> 1 </wgt>" << std::endl;  
  //std::cout<< "<wgt id='kc10'> 1 </wgt>" << std::endl;  
  //std::cout<<"</rwgt>"<<std::endl;
  ///need to use good
  if(LHEEventHandle.isValid())
  {
     // clear the defined vectors before start

     LHE = LHEEventHandle.product();

     //std::cout<<"<rwgt>"<<std::endl;
     ////std::cout<< "<wgt id='SM'> 1 </wgt>" << std::endl;  
     ////std::cout<< "<wgt id='kc1'> 1 </wgt>" << std::endl;  
     ////need for other
     for(const auto& weight : LHE->weights()) {
     //   //std::cout<< "weight.id: " << weight.id << "\t weight.wgt: " << weight.wgt << std::endl;
     //   std::cout<< "<wgt id='" << weight.id << "'> " << weight.wgt << " </wgt>" << std::endl;
     }
     //std::cout<<"</rwgt>"<<std::endl;
  }
  //std::cout<<"</event>"<<std::endl;
  
  Vec_Wboson.clear();
  Vec_Higgs.clear();
  Vec_wmJET.clear();
  Vec_wpJET.clear();
  Vec_Photons.clear();
  Vec_Zboson.clear();
  Vec_Lep.clear();
  Vec_Lep_id.clear();
  Vec_Jet.clear();
  Vec_Muon.clear();
  //Vec_Muon_Pt.clear();
  for(int i=0; i<4; i++)
  {
      gen_Lep_Hindex[i]=-1;
  }
  
  Clear();
}

// ------------ method called once each job just before starting event loop  ------------
void 
GenAnalyzer::beginJob()
{
  std::cout<<"Inside beginJob()"<<std::endl;
  outputFile_ = new TFile(OutPutFileName_,"RECREATE");
  outputFile_->SetCompressionLevel(2);
  tree = new TTree("otree","GenParticles Basic Info");
  
  SetBranches();
}

// ------------ method called once each job just after ending the event loop  ------------
void 
GenAnalyzer::endJob() 
{
  outputFile_->Write();
  outputFile_->Close();
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
GenAnalyzer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}


//define this as a plug-in
DEFINE_FWK_MODULE(GenAnalyzer);
