// Add a Gaussian fitting routine for finding events centroids
// Still a work in progress

// ------------------------------- //
// ------- Data structures ------- //
// ------------------------------- // 
  
// Raw data entries
struct Entry{
  uint16_t xpos;           // wire number (=200 if stripe event)
  uint16_t ypos;           // stripe number (=200 if wire event)
  uint16_t amp;            // ToT in clock cycles (12.5 ns)
  unsigned long long time;     // The full time stamp in clocks (12.5 ns)
  uint8_t eventID;        // 0 for real events, 1 for self triggers 
  uint32_t eventTS;        // The 19 bit time stamp within the buffer
  uint8_t mcpdID;         // 1 for segment 1, 2 for segment 2
  uint8_t status;         // 
  unsigned long long  param0;         // unused
  unsigned long long  param1;         // unused
  unsigned long long  param2;         // unused
  unsigned long long  param3;         // unused
} entry;

// Events - entries in the sorted TTree
struct Event{
  float xpos=0;                // calculated centroid in x
  float ypos=0;                // calculated centroid in y
  int ToTx=0;                  // summed ToT in x
  int ToTy=0;                  // summed ToT in y
  int xToTx=0;                 // summed product of ToT and position for wires
  int yToTy=0;                 // summed product of ToT and position for stripes
  int multx=0;                 // multiplicity in x
  int multy=0;                 // multiplicity in y
  uint64_t time=0;            // time of first entry in event
  int rawevtnum=0;             // location of first entry in raw data file
  int minx=1000; int maxx=0;
  int miny=1000; int maxy=0;
  int widthx=0; int widthy=0;  // Difference in max and min channels in a single event
  int seg=0;                   // Segment number
  long long dtime=0;           // Time since last event
  float xfit = 0;		// Fitted x position
  float yfit = 0;		// Fitted y position
  //TGraph gx;
  //TGraph gy;
} event;


// Time window for events
const int time_window = 30;
const Event emptyevent;

// Event buffer for both segments
Event evtbuff[2];

// Define global TGraphs for x and y ToT distribution
// Consider adding these as Struct members later
//TGraph *gx = new TGraph();
//TGraph *gy = new TGraph();
TF1 *gfit = new TF1("gfit","gaus");
TGraph gx[2];
TGraph gy[2];

bool DrawEvent = 0;  // Switch on to draw ToT plots
TCanvas *can;
TPad *histpad;
TPad *histtitle;
TText *titletext;

// ------------------------- //
// ------- Functions ------- //
// ------------------------- // 

void AddEntry(Entry entry, int row);
void CalculateEvent(int seg);


// ----------------------------------------- // 
// --- Add an entry to the current event --- //cd 
// ----------------------------------------- // 

void AddEntry(Entry entry, int row, TTree *data){

  // Get the segment number (0 or 1 to match array)
  // This determines which event buffer to work with
  int seg = entry.mcpdID - 1;
  
  //cout << "Time: " << entry.time << endl;
  //cout << "Time: " << entry.time << endl;
  //cout << "TimeBuf: " << evtbuff[seg].time << endl;
  
  // If the time to the last event is greater than the time window then event is over
  if (entry.time - evtbuff[seg].time > time_window){
  
    //cout << "\n### Event over ###" << endl;
    //cout << "Row number: " << evtbuff[seg].rawevtnum << endl;
    //cout << "Segment: " << seg+1 << endl;
    //cout << "Timestamp: " << evtbuff[seg].time << endl;
    //cout << "Multx: " << evtbuff[seg].multx << endl;
    //cout << "Multy: " << evtbuff[seg].multy << endl;
       
    
    // Unless the x and y multiplicity is at least 1, ignore event
    if(evtbuff[seg].multx>0 && evtbuff[seg].multy>0){
      CalculateEvent(seg);
      event = evtbuff[seg];
      
      // Only fill if there is at least multiplicity 1 in x and y
      //if (evtbuff[seg].multx>0 && evtbuff[seg].multy>0){
      data->Fill();
   
      if(gx[seg].GetN()>0 && gy[seg].GetN()>0 &&DrawEvent){
        histpad->cd(1);
        gx[seg].Draw("AP");
        histpad->cd(2);
        gy[seg].Draw("AP");
        cout << "x: " << evtbuff[seg].xfit << endl;
        gx[seg].Print();
        cout << "y: " << evtbuff[seg].yfit << endl;
        gy[seg].Print();
        can->WaitPrimitive();
      }
    }
    
    else{
      //cout << "Invalid event" << endl;
    }
    
    //cout << "### Event over ###\n" << endl;
    
    
    // Still need to do this step - dtime calculation needs adjustment
    
    // Store the time difference
    long long dtime = entry.time - evtbuff[seg].time;
    // overwrite the buffer with an empty event
    evtbuff[seg] = emptyevent;  
    evtbuff[seg].time = entry.time;
    evtbuff[seg].rawevtnum = row;
    evtbuff[seg].seg = seg;
    evtbuff[seg].dtime = dtime;
    // Reset the TGraphs
    gx[seg].Set(0);
    gy[seg].Set(0);	
  }  
  
  // Fill if wire
  if(entry.ypos==0){
    //cout << "wire: " << entry.xpos << " segment: " << seg+1 << endl;
    evtbuff[seg].ToTx += entry.amp;
    evtbuff[seg].xToTx += (entry.xpos * entry.amp);
    evtbuff[seg].multx ++;
    if(entry.xpos>evtbuff[seg].maxx) evtbuff[seg].maxx = entry.xpos;
    if(entry.xpos<evtbuff[seg].minx) evtbuff[seg].minx = entry.xpos;
    // Add point to x coord TGraph
    int gxN = gx[seg].GetN(); // Number of indices in current TGraph
    gx[seg].SetPoint(gxN, entry.xpos, entry.amp);
    
  }
  // Fill if stripe (and remove 512 channel offset)
  else{
    //cout << "stripe: " << entry.ypos << " segment: " << seg+1 << endl;
    evtbuff[seg].ToTy += entry.amp;
    evtbuff[seg].yToTy += ((entry.ypos - 512) * entry.amp);
    evtbuff[seg].multy ++;
    if(entry.ypos>evtbuff[seg].maxy) evtbuff[seg].maxy = entry.ypos;
    if(entry.ypos<evtbuff[seg].miny) evtbuff[seg].miny = entry.ypos;
    int gyN = gy[seg].GetN(); // Number of indices in current TGraph
    gy[seg].SetPoint(gyN, entry.ypos - 512, entry.amp);
  }  
}


// ------------------------------------------ //
// ------- Calculate event parameters ------- //
// ------------------------------------------ //

void CalculateEvent(int seg){

  // x position - only if there are wire signals	
  if (evtbuff[seg].multx > 0){
    evtbuff[seg].xpos = float(evtbuff[seg].xToTx)/evtbuff[seg].ToTx;
    evtbuff[seg].xpos += (seg*128);
    evtbuff[seg].widthx = evtbuff[seg].maxx - evtbuff[seg].minx;
    
    // Fit the ToT data with a Gaussian if the multiplicity is at least 2
    gx[seg].Sort();
    if (evtbuff[seg].multx > 1){
      gfit->SetParameter(1, evtbuff[seg].xpos - (seg*128));
      gx[seg].Fit(gfit, "Q");
      evtbuff[seg].xfit = gfit->GetParameter(1) + (seg*128);
    }
    else evtbuff[seg].xfit = evtbuff[seg].xpos;
  }
  // If no wire signals set xpos=-10
  else{
    evtbuff[seg].xpos = -10;
    evtbuff[seg].xfit = -10;
  }
  
  // y position - only if there are stripe signals
  if (evtbuff[seg].multy > 0){
    evtbuff[seg].ypos = float(evtbuff[seg].yToTy)/evtbuff[seg].ToTy;
    evtbuff[seg].widthy = evtbuff[seg].maxy - evtbuff[seg].miny;
    
    // Fit the ToT data with a Gaussian if the multiplicity is at least 2
    gy[seg].Sort();
    if (evtbuff[seg].multy > 1){
      gfit->SetParameter(1, evtbuff[seg].ypos);
      gy[seg].Fit(gfit, "Q");
      evtbuff[seg].yfit = gfit->GetParameter(1);
    }
    else evtbuff[seg].yfit = evtbuff[seg].ypos;
  }
  // If no stripe signals set ypos=-10
  else{
    evtbuff[seg].ypos = -10;
    evtbuff[seg].yfit = -10;
  }
}


// -------------------- //
// ------- Main ------- //
// -------------------- // 


void GaussSort(TString filename){
  
  // Get the names of the input and output ROOT files from the argument passed
  TString outfilename = filename;
  outfilename.ReplaceAll(".root","_sorted.root");
    
  
  // ---------------------------------- //
  // ------- Get raw data TTree ------- //
  // ---------------------------------- // 

  // Open the raw data ROOT file and get the TTree
  TFile *infile = TFile::Open(filename);
  TTree *rawdata = (TTree*)infile->Get("rawdata"); 
  
  rawdata->SetBranchAddress("xpos",&entry.xpos);  
  rawdata->SetBranchAddress("ypos",&entry.ypos);
  rawdata->SetBranchAddress("amp",&entry.amp);
  rawdata->SetBranchAddress("time",&entry.time);
  rawdata->SetBranchAddress("eventID",&entry.eventID);
  rawdata->SetBranchAddress("eventTS",&entry.eventTS);
  rawdata->SetBranchAddress("mcpdID",&entry.mcpdID);
  rawdata->SetBranchAddress("status",&entry.status);
  rawdata->SetBranchAddress("param0",&entry.param0);
  rawdata->SetBranchAddress("param1",&entry.param1);
  rawdata->SetBranchAddress("param2",&entry.param2);
  rawdata->SetBranchAddress("param3",&entry.param3);
    
  
  // Create a new TFile and TTree for the sorted data
  TFile *outfile = new TFile(outfilename,"RECREATE");
  TTree *data = new TTree("data","Sorted data");
    
  
  // ----------------------------------- //
  // ------- Set up sorted TTree ------- //
  // ----------------------------------- //
  
  
  data->Branch("xpos", &event.xpos, "xpos/F");
  data->Branch("ypos", &event.ypos, "ypos/F");
  data->Branch("ToTx", &event.ToTx, "ToTx/I");
  data->Branch("ToTy", &event.ToTy, "ToTy/I");
  data->Branch("multx", &event.multx, "multx/I");
  data->Branch("multy", &event.multy, "multy/I");
  data->Branch("time", &event.time, "time/L");
  data->Branch("rawevtnum", &event.rawevtnum, "rawevtnum/I");
  data->Branch("seg", &event.seg, "seg/I");
  data->Branch("widthx", &event.widthx, "widthx/I");
  data->Branch("widthy", &event.widthy, "widthy/I");
  data->Branch("dtime", &event.dtime, "dtime/L");
  data->Branch("xfit", &event.xfit, "xfit/F");
  data->Branch("yfit", &event.yfit, "yfit/F");
  //data->Branch("gx", &event.gx); 
  //data->Branch("gy", &event.gy); 
  // New branch for x and y position from Gaussian fits
  // Also store the TGraphs themselves
  
  
  // Add new entries for Gaussian centroid in x and y
  // Also add a TGraph in x and y for each entry?
  
  
  // Create the canvas if drawing is enabled
  if (DrawEvent==1){
    can = new TCanvas("can", "can", 1000,600);
    histpad = new TPad("Plots", "Plots", 0.01, 0.05, 0.95, 0.95);
    histpad->Draw();
    histpad->cd();
    histpad->Divide(2,1);
    
    histtitle = new TPad("title", "title", 0.1,0.96,0.9,0.99);
    histtitle->Draw();
    titletext = new TText(0.5,0.5,"Title goes here");
    titletext->SetTextSize(0.8);
    titletext->SetTextAlign(22);
    
    gx[0].SetMarkerStyle(3);
    gx[1].SetMarkerStyle(3);
    gy[0].SetMarkerStyle(3);
    gy[1].SetMarkerStyle(3);

  }
  
  // ------------------------------------- //
  // ------- Loop over all entries ------- //
  // ------------------------------------- //  
  
  // Get number of rows and start loop
  int num_rows = rawdata->GetEntries();
  
  for (int row=0; row<num_rows; row++){

    //cout << "Entry " << row << endl;
    
    // Fill 'entry' with raw data from the current row
    rawdata->GetEntry(row); 
    
    // Print status
    if(row%10000 == 0) cout << "Entry " << row << " of " << num_rows << "\r" << flush; 
    
    // Check the eventID - only read if 0
    if(entry.eventID !=0) continue;
 
    // Pass the entry and the current row to be processed
    AddEntry(entry, row, data);
    
  }
  
  
  // ------------------------------------- //
  // ------- Read out final events ------- //
  // ------------------------------------- //
  
  // Only if there are both x and y data in the buffers
  if(evtbuff[0].multx>0 && evtbuff[0].multy>0){
    CalculateEvent(0);
    event = evtbuff[0];
    data->Fill();
  }
  if(evtbuff[1].multx>0 && evtbuff[1].multy>0){
    CalculateEvent(1);
    event = evtbuff[1];
    data->Fill();
  }
  
  
  // ----------------------- //
  // ------- Tidy up ------- //
  // ----------------------- //    
  
  data->Write();
  outfile->Close();
  infile->Close();

}
