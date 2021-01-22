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
  double t_secs = t_tot/100.e-9;
  // -----------------------------------//
  
  // ----------- Position spectra ------ //
  TH1I *hx = new TH1I("hx","hx",512,0,256);
  TH1I *hy = new TH1I("hy","hy",256,0,128);
  rawdata->Draw("x>>hx(512,0,256)","","goff");
  rawdata->Draw("y>>hy(256,0,128)","","goff");
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
  hToT->Fit(f,"","");
  double ToT_max = f->GetParameter(1);
  double ToT_width = f->GetParameter(2);
  f->SetParameters(hToTS1->GetMaximum(), hToTS1->GetMaximumBin(), 10.);
  hToTS1->Fit(f,"","");
  double ToTS1_max = f->GetParameter(1);
  double ToTS1_width = f->GetParameter(2);
  f->SetParameters(hToTS2->GetMaximum(), hToTS2->GetMaximumBin(), 10.);
  hToTS2->Fit(f,"","");
  double ToTS2_max = f->GetParameter(1);
  double ToTS2_width = f->GetParameter(2);
  f->SetParameters(hToTB->GetMaximum(), hToTB->GetMaximumBin(), 10.);
  hToTB->Fit(f,"","");
  double ToTB_max = f->GetParameter(1);
  double ToTB_width = f->GetParameter(2);  
  
  // ------------------------------------------------- //
  
  // -------- Waiting time --------- //
  // Histogram for storing the waiting time between events 
  TH1D *timediff = new TH1D("timediff","timediff",100000,0,1e6);

  // Set time variable for use in loop 
  ULong64_t time;
  rawdata->SetBranchAddress("time", &time);
  
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
  }
  // ---------------------------------//


  // Output results
  // entries, entriesS1, entriesS2, entriesB
  // timediff
  cout << entries << "\t" << entriesS1 << "\t" << entriesS2 << "\t" << entriesB << endl;
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
  outfile->Close();
  
  //delete outfile;
  //delete timediff;
}
