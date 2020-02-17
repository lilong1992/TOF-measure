Double_t triggertimefastcal(TDoubleWaveform& rawwaveform,Bool_t subtract=true,Double_t halfbaseline=450,Double_t F=0.5,Int_t D=3){
  if (subtract){
    Double_t baseline=rawwaveform.Sum(0,100)/100;
    Double_t cft0=0,cft1=1;
    Int_t index=100;
    while (cft0<20&&index<200) {
      cft0=F*rawwaveform.At(index+D)-rawwaveform.At(index)+(1-F)*baseline;
      index++;
    }
    if(index>=200) return -1;
    cft1=F*rawwaveform.At(index+D)-rawwaveform.At(index)+(1-F)*baseline;
    while (cft1>0&&index<200) {
      index++;
      cft0=cft1;
      cft1=F*rawwaveform.At(index+D)-rawwaveform.At(index)+(1-F)*baseline;
    }
    if(index>=200) return -1.5;
    return index+cft0/(cft0-cft1); //zero crossing point
  }
  else {
    Double_t cft0=0,cft1=1;
    Int_t index=0;
    while (cft0<900&&index<200) {
      cft0=0.5*rawwaveform.At(index+D)-rawwaveform.At(index)+halfbaseline;
      index++;
    }
    if(index>=200) return -100;
    cft1=0.5*rawwaveform.At(index+D)-rawwaveform.At(index)+halfbaseline;
    while (cft1>0&&index<200) {
      index++;
      cft0=cft1;
      cft1=0.5*rawwaveform.At(index+D)-rawwaveform.At(index)+halfbaseline;
    }
    if(index>=200) return -150;
    return index+cft0/(cft0-cft1);
  }
}

void beamenergytof(TString filename){
  TChain* tree = new TChain( "sis3316tree" );
  TString filenameBase = filename( 0, filename.Last( '_' ) + 1 );
  filenameBase.Append( "*.root" );
  Int_t nFilesAdded = tree->Add( filenameBase );


  UInt_t nSamples;
  UShort_t waveform[65536];
  UShort_t peakHighValue;
  UShort_t channelID;
  ULong64_t timestamp;
  UInt_t mawMaximumValue;
  UInt_t mawValueAfterTrigger;
  UInt_t mawValueBeforeTrigger;

  tree->SetBranchAddress( "nSamples", &nSamples );
  tree->SetBranchAddress( "waveform", waveform );
  tree->SetBranchAddress( "peakHighValue", &peakHighValue );
  tree->SetBranchAddress( "channelID", &channelID );
  tree->SetBranchAddress( "timestamp", &timestamp );
  tree->SetBranchAddress( "mawMaximumValue", &mawMaximumValue );
  tree->SetBranchAddress( "mawValueAfterTrigger", &mawValueAfterTrigger );
  tree->SetBranchAddress( "mawValueBeforeTrigger", &mawValueBeforeTrigger );

  Long64_t nEntries = tree->GetEntries();
  Long64_t nevents=nEntries/2;
  cout<<"total events:"<<nevents<<endl;
  auto h1=new TH2F("h1","PSD vs TOF;TOF (ns);PSD",8000,0,800,300,0,3);
  auto h2=new TH1F("h2","0_degree_signals;index;counts",3020,-2,300);
  auto h3=new TH1F("h3","TOF histogram;TOF (ns);counts",8000,0,800);
  auto h4=new TH2F("h4","PSD vs Integral;Integral (ADCunit);PSD",10000,0,100000,300,0,3);
  auto h5=new TH1F("h5","TOF histogram;TOF (ns);counts",8000,0,800);
  auto h6=new TH2F("h6","Integral vs TOF;TOF (ns);Integral (ADCunit)",8000,0,800,10000,0,100000);
  auto h7=new TH1F("h7","n MAW-CFD;ns;counts",6000,0,60);
  auto h8=new TH1F("h8","g MAW-CFD;ns;counts",6000,0,60);
  auto h9=new TH1F("h9","energyspectrum tof>600&&tof<645;integral;counts",1000,0,100000);
  auto h10=new TH1F("h10","0_degree_signals;index;counts",3020,-2,300);
  auto h11=new TH1F("h11","TOF histogram all;TOF (ns);counts",8000,0,800);
  auto h13=new TH1F("h13","TOF histogram LE;TOF (ns);counts",8000,0,800);
  auto h12=new TH1F("h12","neutron energy;keV;counts",1000,0,200);
  TDoubleWaveform rawwaveform;
  ULong64_t timestampcheck=1;
  Double_t bd0,bpm,mawtrigger,tof=0,tgam=230.9,dist=169.6;
  Double_t PSD=0,Integral=0;
  Bool_t notsaturated=true;
  Int_t peaktime=4,gaptime=1;

  for(Long64_t entryNumber=0;entryNumber<nEntries;entryNumber++){
    if (entryNumber%30000==0) cout<<entryNumber/2<<"events done!"<<endl;

      tree->GetEntry(entryNumber);
      rawwaveform.SetData(waveform,nSamples);
      if(channelID==0){//calculate PSD and Integral
        if (peakHighValue>16380) notsaturated=false;
        Double_t baseline=rawwaveform.Sum(0,50)/50;
        bd0=triggertimefastcal(rawwaveform);
        Int_t j=(int) (bd0+0.5);
        Double_t SG=rawwaveform.Sum(j-5,j+7)-baseline*12;
        Integral=rawwaveform.Sum(j-5,j+170)-baseline*175;
        PSD=Integral/SG;

        //the following calculation is for debugging DAQ, no need for TOF
        TDoubleWaveform tmpwaveform;
        tmpwaveform.SetLength(nSamples - peaktime);
        tmpwaveform[0]=rawwaveform.Sum(1,1+peaktime);
        for (Int_t i=1;i<tmpwaveform.GetLength();i++) tmpwaveform[i]=tmpwaveform[i-1]-rawwaveform[i]+rawwaveform[i+peaktime];
        TDoubleWaveform TFilteredwavefrom;
        TFilteredwavefrom.SetLength(nSamples-2*peaktime-gaptime);
        TFilteredwavefrom=tmpwaveform.SubWaveform(peaktime+gaptime,tmpwaveform.GetLength());
        tmpwaveform.MakeSimilarTo(TFilteredwavefrom);
        TFilteredwavefrom-=tmpwaveform;

        //the following calculation is for debugging DAQ, no need for TOF
        Int_t maxindex=100;
        Double_t halfmawM=0,mawA,mawB;
        for(Int_t i=100;i<200;i++){
          if(TFilteredwavefrom.At(i)>halfmawM){
            halfmawM=TFilteredwavefrom.At(i);
            maxindex=i;
          }
        }
        halfmawM/=2;
        mawB=TFilteredwavefrom.At(maxindex);
        mawA=TFilteredwavefrom.At(maxindex+1);
        while (mawA>halfmawM) {
          maxindex++;
          mawB=mawA;
          mawA=TFilteredwavefrom.At(maxindex);
        }
        mawtrigger=maxindex+(halfmawM-mawB)/(mawA-mawB);


      }
      if(channelID==4){//find two bpm time use the nearist but come before bd0
        bpm=triggertimefastcal(rawwaveform,false);
      }
    //end of an event
    if (timestamp==timestampcheck){
      if (notsaturated && bd0>0){
        tof=4*(bd0-bpm);
        Double_t triggerdiff=4*(mawtrigger+13-bd0);
        Double_t ntof=tof+400;
        if(tof<200) tof+=400;
        h1->Fill(tof,PSD);
        h2->Fill(bd0);
        h10->Fill(mawtrigger);
        h6->Fill(tof,Integral);
        h4->Fill(Integral,PSD);
        h11->Fill(tof);
        if (Integral<1000) h13->Fill(tof);
        if(PSD>1.28) {h3->Fill(tof);
          h7->Fill(triggerdiff);
          //if(tof>600&&tof<645) h9->Fill(Integral);
          Double_t en=1000*0.5*939.56/(1+29.98*(ntof-tgam)/dist)/(1+29.98*(ntof-tgam)/dist);
          h12->Fill(en);
        }
        else {if (Integral>2000) h5->Fill(tof);
          h8->Fill(triggerdiff);}
      }
      notsaturated=true; //clean variables
  }
  else timestampcheck=timestamp;
  }
  auto c=new TCanvas("c","tofplots");
  c->Divide(2,2);
  c->cd(1);
  h1->Draw("COLZ");
  c->cd(2);
  h6->Draw("COLZ");
  //h7->Draw();
  c->cd(3);
  h3->SetLineColor(1);
  h3->Draw();
  h5->SetLineColor(2);
  h5->Draw("same");
  leg=new TLegend(0.7,0.8,0.8,0.89);
  leg->AddEntry(h3,"neutron");
  leg->AddEntry(h5,"gamma");
  leg->Draw();
  //h8->Draw();
  //gStyle->SetOptStat(0);
  c->cd(4);
  h4->Draw("COLZ");
  //h9->Draw();
  c->Update();
  auto c2=new TCanvas("c2","tofplots");
  c2->Divide(1,2);
  c2->cd(1);
  h11->Draw();
  c2->cd(2);
  h12->Draw();
  //h13->Draw();
  c2->Update();

}
