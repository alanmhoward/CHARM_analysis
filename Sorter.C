// Macro for analysing root TTree containing raw data from digital readout

// Try to create events and hold them in a buffer
// Check that new events do not belong to other events in the buffer
// Eventually write out events and make space for new events in the buffer

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
} event;


// Time window for events
const int time_window = 30;
const Event emptyevent;

// Event buffer up to nine segments
Event evtbuff[9];


// ------------------------- //
// ------- Functions ------- //
// ------------------------- // 

void AddEntry(Entry entry, int row);
void CalculateEvent(int seg);


// ----------------------------------------- // 
// --- Add an entry to the current event --- //
// ----------------------------------------- // 

void AddEntry(Entry entry, int row, TTree *data){

  // Get the segment number (0 or 1 to match array)
  // This determines which event buffer to work with
  int seg = entry.mcpdID - 1;
  
  // For first event need to avoid writing
  
  // If the time to the last event is greater than the time window then event is over
  if (entry.time - evtbuff[seg].time > time_window){
  
    CalculateEvent(seg);
    event = evtbuff[seg];
    
    data->Fill();
    
    // Store the time difference
    long long dtime = entry.time - evtbuff[seg].time;
    // overwrite the buffer with an empty event
    evtbuff[seg] = emptyevent;  
    evtbuff[seg].time = entry.time;
    evtbuff[seg].rawevtnum = row;
    evtbuff[seg].seg = seg;
    evtbuff[seg].dtime = dtime;
  }  
  
  // Fill if wire
  if(entry.ypos==0){
    evtbuff[seg].ToTx += entry.amp;
    evtbuff[seg].xToTx += (entry.xpos * entry.amp);
    evtbuff[seg].multx ++;
    if(entry.xpos>evtbuff[seg].maxx) evtbuff[seg].maxx = entry.xpos;
    if(entry.xpos<evtbuff[seg].minx) evtbuff[seg].minx = entry.xpos;
    
  }
  // Fill if stripe (and remove 512 channel offset)
  else{
    evtbuff[seg].ToTy += entry.amp;
    evtbuff[seg].yToTy += ((entry.ypos - 512) * entry.amp);
    evtbuff[seg].multy ++;
    if(entry.ypos>evtbuff[seg].maxy) evtbuff[seg].maxy = entry.ypos;
    if(entry.ypos<evtbuff[seg].miny) evtbuff[seg].miny = entry.ypos;
  }  
}


// ------------------------------------------ //
// ------- Calculate event parameters ------- //
// ------------------------------------------ //

void CalculateEvent(int seg){
  if (evtbuff[seg].multx > 0){
    evtbuff[seg].xpos = float(evtbuff[seg].xToTx)/evtbuff[seg].ToTx;
    evtbuff[seg].xpos += (seg*128);
    evtbuff[seg].widthx = evtbuff[seg].maxx - evtbuff[seg].minx;
  }
  else evtbuff[seg].xpos = -10;
  if (evtbuff[seg].multy > 0){
    evtbuff[seg].ypos = float(evtbuff[seg].yToTy)/evtbuff[seg].ToTy;
    evtbuff[seg].widthy = evtbuff[seg].maxy - evtbuff[seg].miny;
  }
  else evtbuff[seg].ypos = -10;
}


// -------------------- //
// ------- Main ------- //
// -------------------- // 


void Sorter(TString filename){
  
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
  
  
  // ------------------------------------- //
  // ------- Loop over all entries ------- //
  // ------------------------------------- //  
  
  // Get number of rows and start loop
  int num_rows = rawdata->GetEntries();
  
  for (int row=0; row<num_rows; row++){
    
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
 

  for (int i=0; i<9; i++){
    CalculateEvent(i);
    event = evtbuff[i];
    data->Fill();
  }
  
  // ----------------------- //
  // ------- Tidy up ------- //
  // ----------------------- //    
  
  data->Write();
  outfile->Close();
  infile->Close();

}
