// Test sorting speed using different approaches 
// To use in a bash loop (and print time in seconds for each file) do
// for i in run???.root; do t1=$SECONDS; root -b -q -l '../code/SpeedTest.C("'${i}'")'; t2=$SECONDS; echo $(($t2-$t1)); done

void SpeedTest(TString filename){

  // Open the raw data ROOT file and get the TTree
  TFile *infile = TFile::Open(filename);
  TTree *rawdata = (TTree*)infile->Get("rawdata");
  
  // Create an output file
  TString outfilename = filename;
  outfilename.ReplaceAll(".root","_hists.root");
  TFile *outfile = new TFile(outfilename, "RECREATE");
  
  // Use LoadBaskets to load TTree into memory and drastically speed up random access, e.g. when using an index to loop through the events
  rawdata->LoadBaskets();

  // Set alias for x and y position
  rawdata->SetAlias("x","(xpos + 1024*(mcpdID-1))/8.");
  rawdata->SetAlias("y","ypos/8.");


  
  // ----------- Rate calc --------------//
  
  // Get the number of entries
  // For the boundary use a window of 2 wire pitches either side of the centre
  int entries = rawdata->GetEntries();
  int entriesS1 = rawdata->GetEntries("mcpdID==1");
  int entriesS2 = rawdata->GetEntries("mcpdID==2");
  int entriesB = rawdata->GetEntries("abs(xpos-127.5)<2");
  
  // Get the total run time
  ULong64_t t_min = rawdata->GetMinimum("time");
  ULong64_t t_max = rawdata->GetMaximum("time");
  ULong64_t t_tot = t_max-t_min; // in units of 100ns
  double t_secs = t_tot*100.e-9;
  // -----------------------------------//
  
  // ----------- Position spectra ------ //
  TH1I *hx = new TH1I("hx","hx",512,0,256);
  TH1I *hy = new TH1I("hy","hy",256,0,128);
  rawdata->Draw("x>>hx(512,0,256)","","goff");
  rawdata->Draw("y>>hy(256,0,128)","","goff");

  // Use the chi squared of a pol0 fit to quantify the fine structure in the position spectra
  TF1 *fx = new TF1("fx","[0]");
  // Boundary region
  hx->Fit(fx,"","Q",117.5,137.5);
  double xposB_mean = fx->GetParameter(0);
  double xposB_chi = fx->GetChisquare();
  double xposB_NDF = fx->GetNDF();
  // Seg1
  hx->Fit(fx,"","Q",97.5,117.5);
  double xposS1_mean = fx->GetParameter(0);
  double xposS1_chi = fx->GetChisquare();
  double xposS1_NDF = fx->GetNDF();
  // Seg2
  hx->Fit(fx,"","Q",137.5,157.5);
  double xposS2_mean = fx->GetParameter(0);
  double xposS2_chi = fx->GetChisquare();
  double xposS2_NDF = fx->GetNDF();

  // ----------------------------------- //
  
  // ---------------- ToT spectra -------------------- //
  TH1I *hToT = new TH1I("hToT","hToT",256,0,256);
  TH1I *hToTS1 = new TH1I("hToTS1","hToTS1",256,0,256);
  TH1I *hToTS2 = new TH1I("hToTS2","hToTS2",256,0,256);
  TH1I *hToTB = new TH1I("hToTB","hToTB",256,0,256);
  rawdata->Draw("amp>>hToT(256,0,256)","","goff");
  rawdata->Draw("amp>>hToTS1(256,0,256)","x<125","goff");
  rawdata->Draw("amp>>hToTS2(256,0,256)","x>130","goff");
  rawdata->Draw("amp>>hToTB(256,0,256)","abs(x-127.5)<2","goff");
  
  // Skewed Gaussian for fitting ToT spectra
  // [0] amplitude, [1] centroid, [2] width, [3] left hand skew, [4] right hand skew
  TF1 *f = new TF1("sgf","[0]*exp(-pow((x-[1])/([2]+ ( (x<[1])*[3] + (x>[1])*[4] ) *(x-[1])),2))");
  // Set initial parameters
  f->SetParameters(hToT->GetMaximum(), hToT->GetMaximumBin(), 10.);
  hToT->Fit(f,"","Q");
  double ToT_max = f->GetParameter(1);
  double ToT_width = f->GetParameter(2);
  f->SetParameters(hToTS1->GetMaximum(), hToTS1->GetMaximumBin(), 10.);
  hToTS1->Fit(f,"","Q");
  double ToTS1_max = f->GetParameter(1);
  double ToTS1_width = f->GetParameter(2);
  f->SetParameters(hToTS2->GetMaximum(), hToTS2->GetMaximumBin(), 10.);
  hToTS2->Fit(f,"","Q");
  double ToTS2_max = f->GetParameter(1);
  double ToTS2_width = f->GetParameter(2);
  f->SetParameters(hToTB->GetMaximum(), hToTB->GetMaximumBin(), 10.);
  hToTB->Fit(f,"","Q");
  double ToTB_max = f->GetParameter(1);
  double ToTB_width = f->GetParameter(2);  
  
  // ------------------------------------------------- //
  
  // -------- Waiting time --------- //
  // Histogram for storing the waiting time between events 
  TH1D *timediff = new TH1D("timediff","timediff",1000000,0,1e6);

  // Set time variable for use in loop 
  ULong64_t time;
  rawdata->SetBranchAddress("time", &time);
  
  // Set variables for diagnosing events with small waiting times
  UShort_t xpos, ypos, amp;
  rawdata->SetBranchAddress("xpos", &xpos);
  rawdata->SetBranchAddress("ypos", &ypos);
  rawdata->SetBranchAddress("amp", &amp);
  UChar_t mcpdID;
  rawdata->SetBranchAddress("mcpdID", &mcpdID);

  // x and y position of events with short waiting times
  TH1I *hxdt = new TH1I("hxdt","hxdt",512,0,256); 
  TH1I *hydt = new TH1I("hydt","hydt",256,0,128);
  TH2I *hxydt = new TH2I("hxydt","hxydt",256,0,256,128,0,128);
  TH1I *hampdt = new TH1I("hampdt","hampdt",256,0,256);
  TH1I *hdxdt = new TH1I("hdxdt","hdxdt",512,0,256);

  // Build an index of the TTree using event time
  // Use this to loop through events in chronological order
  rawdata->BuildIndex("time");
  TTreeIndex *index = (TTreeIndex*)rawdata->GetTreeIndex();

  // Fill timediff using all events and file sorted in time order 
  ULong64_t ti = 0;
  ULong64_t dt, tf;
  ULong64_t j;
  
  // Loop over all events using time index for order
  // This part is vastly sped up by the TTree->LoadBaskets() method
  for (int i = 1; i < entries; i++){
    j = index->GetIndex()[i]; 
    rawdata->GetEntry(j);
    tf = time;
    dt = tf-ti;
    timediff->Fill(dt);
    ti=tf;
    // Select events based on time since previous event
    /* 
    if (dt<2){
      int segID = 1;
      if (mcpdID==0x02) segID = 2;
      double xdt = (xpos + 1024 * (segID-1))/8.;
      double ydt = ypos / 8.;
      hxdt->Fill(xdt);
      hydt->Fill(ydt);
      hxydt->Fill(xdt,ydt);
      hampdt->Fill(amp);
    }
   */ 
  }

  // ---------------------------------//


  // ------------- Check time ordering of events ---------------- //
  // Track the last recorded time in each segment separately
  ULong64_t tS1 = 0;
  ULong64_t tS2 = 0;
      
  // Loop over all events and print a message whenever the time regresses
  
  for (int i = 1; i < entries; i++){
    rawdata->GetEntry(i);
    if(mcpdID==0x01){
      if (time < tS1) cout << "entry " << i << " " << tS1 << " " << time << endl;
      tS1 = time;
    }
    else{
      if (time < tS2) cout << "entry " << i << " " << tS2 << " " << time << endl;
      tS2 = time;
    }
  }
   
 
  // ------------------------------------------------------------ //

  // Output results
  // entries, entriesS1, entriesS2, entriesB
  // timediff

  cout << t_secs << "\t" << entries << "\t" << entriesS1 << "\t" << entriesS2 << "\t" << entriesB << endl;
  cout << xposS1_mean << "\t" << xposS1_chi << "\t" << xposS1_NDF << endl;
  cout << xposS2_mean << "\t" << xposS2_chi << "\t" << xposS2_NDF << endl;
  cout << xposB_mean << "\t" << xposB_chi << "\t" << xposB_NDF << endl;
  cout << ToT_max << "\t" << ToT_width << endl;
  cout << ToTS1_max << "\t" << ToTS1_width << endl; 
  cout << ToTS2_max << "\t" << ToTS2_width << endl; 
  cout << ToTB_max << "\t" << ToTB_width << endl; 
  

  timediff->Write();
  hx->Write();
  hy->Write();
  hToT->Write();
  hToTS1->Write();
  hToTS2->Write();
  hToTB->Write();

  hxdt->Write();
  hydt->Write();
  hxydt->Write();
  hampdt->Write();

  outfile->Close();
  
  //delete outfile;
  //delete timediff;
}
