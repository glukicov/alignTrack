//file loopdir.C

// root > TFile f("myfile with subdirs.root");
//  root > .L loopdir.C
//  root > loopdir(&f);

void loopdir (TDirectory *dir) {
   //loop on all keys of dir including possible subdirs
   //print a message for all keys with a class TH1
   TIter next (dir->GetListOfKeys());
   TKey *key;
   while ((key = (TKey*)next())) {
      // if (strstr(key->GetClassName(),"TH1F")) {
      //    printf (" key : %s is a %s in %s\n",
      //            key->GetName(),key->GetClassName(),dir->GetPath());
      // }
      if (strstr(key->GetClassName(),"TRatioPlot")) {
         printf (" key : %s is a %s in %s\n",
                 key->GetName(),key->GetClassName(),dir->GetPath());
      }
      // if (strstr(key->GetClassName(),"TCanvas")) {
      //    printf (" key : %s is a %s in %s\n",
      //            key->GetName(),key->GetClassName(),dir->GetPath());
        
      //    cout << "gROOT->FindObject(x): " << gROOT->FindObject(key->GetName()); 
      // }
      if (!strcmp(key->GetClassName(),"TDirectory")) {
         dir->cd(key->GetName());
         TDirectory *subdir = gDirectory;
         loopdir(subdir);
         dir->cd();
      }
   }
}
