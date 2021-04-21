// Macro for analysing root TTree containing sorted data from digital readout

// Try to build the events split between segments 

// ------------------------------- //
// ------- Data structures ------- //
// ------------------------------- // 
  

// Events - entries in the sorted TTree
struct Event{
  float xpos=0;                // calculated centroid in x
  float ypos=0;                // calculated centroid in y
  int ToTx=0;                  // summed ToT in x
  int ToTy=0;                  // summed ToT in y
  int multx=0;                 // multiplicity in x
  int multy=0;                 // multiplicity in y
  long long time=0;            // time of first entry in event
  int rawevtnum=0;             // location of first entry in raw data file
  int minx=1000; int maxx=0;
  int miny=1000; int maxy=0;
  int widthx=0; int widthy=0;  // Difference in max and min channels in a single event
  int seg=0;                   // Segment number
  long long dtime=0;           // Time since last event
} event;

// Debug on/off - if on then extra parameters will be added to the final TTree
bool debug = 1;

// Debug parameters - only used if debug == 1
int boundtime = 0;

// The time window correlated events must lie within
int time_window = 30;

// x window to set around the boundary - i.e. 127.5 +- x_window
float x_window = 3.0;

// Boundary events must have matching y position within the y_window
float y_window = 3.0;


// -------------------- //
// ------- Main ------- //
// -------------------- // 

void Correlator(TString filename, int event_start=0, int num_events=0){
  
  // Get the names of the input and output ROOT files from the argument passed
  TString outfilename = filename;
  outfilename.ReplaceAll("_sorted.root","_final.root");
    
  
  // ---------------------------------- //
  // ------- Get raw data TTree ------- //
  // ---------------------------------- // 

  // Open the raw data ROOT file and get the TTree
  TFile *infile = TFile::Open(filename);
  TTree *data = (TTree*)infile->Get("data"); 
  
  
  data->SetBranchAddress("xpos",&event.xpos);
  data->SetBranchAddress("ypos",&event.ypos);
  data->SetBranchAddress("ToTx", &event.ToTx);
  data->SetBranchAddress("ToTy", &event.ToTy);
  data->SetBranchAddress("multx", &event.multx);
  data->SetBranchAddress("multy", &event.multy);
  data->SetBranchAddress("time", &event.time);
  data->SetBranchAddress("rawevtnum", &event.rawevtnum);
  data->SetBranchAddress("seg", &event.seg);
  data->SetBranchAddress("widthx", &event.widthx);
  data->SetBranchAddress("widthy", &event.widthy);
  data->SetBranchAddress("dtime", &event.dtime);	// time since previous event
  
  
  // --------------------------------------- //
  // -------- Clone original TTree --------- //
  // --------------------------------------- // 
  
  // After fetching an entry, the values can be altered before writing the new TTree
  // Entries which are not required can be skipped without writing
  TFile *outfile = new TFile(outfilename,"RECREATE");
  TTree *d = data->CloneTree(0);
  
  // Use the seg entry to record boudary events (set seg = -1)
  // widthx - sum entries
  // mult - sum entries for x, largest of the two individual values for y
  // ToT - sum entries
  if (debug){
    d->Branch("boundtime", &boundtime, "boundtime/s");
  }
  
  
  // ------------------------------------- //
  // ------- Loop over all entries ------- //
  // ------------------------------------- //  
  
  // Vector to track entries which are combined
  vector <int> vec;
  Event evtbuff;
  
  // Get number of rows and start loop
  int num_rows = data->GetEntries();
  
  for (int row=0; row<num_rows; row++){
    
    // Get the sorted data entry
    data->GetEntry(row);
    
    if (row%10000 == 0) cout << "Reading entry " << row << " of " << num_rows << "\r" << flush;
    
    // Skip if the multiplicity in either x or y is zero
    if (event.multx == 0 || event.multy == 0) continue;
    
    // Skip if the event is in the vector of already matched events (and remove entry from vector)
    if (std::find(vec.begin(), vec.end(), row) != vec.end()){
      vec.erase(std::remove(vec.begin(), vec.end(), row), vec.end());
      continue;
    }
    
    // Look for edge events
    // Valid only in the two-segment case
    //if (event.xpos > 125 && event.xpos < 130){
    if (abs(event.xpos - 127.5) < x_window){
    
      // Record the current time and y position since these need to match
      evtbuff = event;

      // look for a corresponding event in the other segment
      // Loop over subsequent rows until a later event is found in the other segment
      // Need to then keep track of the matched event so that it can be skipped later
      int check_row = row;
      while (true && check_row < num_rows){
        
        // Get the next entry to check against
        check_row++;
        data->GetEntry(check_row);

		  // Check if also a boundary event
		  if (abs(event.xpos - 127.5) > x_window) continue;

        // Skip if the multiplicity is low
        if (event.multx == 0 || event.multy == 0) continue;
        
        // If the entry is from the same segment skip forward
        if (event.seg == evtbuff.seg) continue;
        
        // If entries are within the time and y windows correlate them
        if (abs(int64_t(event.time - evtbuff.time)) < time_window){
        
          if ( abs( event.ypos - evtbuff.ypos ) > y_window) continue; 
          
          // Build the new event parameters
          int ToTx = event.ToTx + evtbuff.ToTx;
          int ToTy = event.ToTy + evtbuff.ToTy;
          float xpos = ((event.ToTx * event.xpos) + (evtbuff.ToTx * evtbuff.xpos))/float(ToTx);
          float ypos = ((event.ToTy * event.ypos) + (evtbuff.ToTy * evtbuff.ypos))/float(ToTy);
          
          evtbuff.ToTx = ToTx;
          evtbuff.ToTy = ToTy;
          evtbuff.xpos = xpos;
          evtbuff.ypos = ypos;
          evtbuff.seg = -1; // Set seg=-1 for reconstructed events
          evtbuff.multx += event.multx;
          if (event.multy > evtbuff.multy) evtbuff.multy = event.multy;
          if (event.dtime > evtbuff.dtime) evtbuff.dtime = event.dtime;
          
          boundtime = abs(event.time - evtbuff.time);
          
          event = evtbuff;
          d->Fill();
          
          // Reset debug parameters to zero
          boundtime = 0;
          
          vec.push_back(check_row);
          break;

        }
        
        // If the timestamp on the checked event is earlier, keep looping
        else if(event.time < evtbuff.time) continue;
        
        // No matching events found - go back to original entry and write
        else{ 
          data->GetEntry(row);
          d->Fill();  
          break;
        }
        
      } // End of while loop
      
    }// End of edge event reconstruction
    
    // if the loop was skipped, event is still that given by row
    
    // otherwise the current event is from sometime later
    
    // write the event if the if loop was skipped
    else d->Fill();
    
  }
  
  
  // ----------------------- //
  // ------- Tidy up ------- //
  // ----------------------- //    
  
  d->AutoSave();
  delete infile;
  delete outfile;
  //outfile->Close();
  //infile->Close();

}
