# EcalStudies

See the cfi file in python/ and the example cfg in test/ for detailed information. Developed in CMSSW_8_0_19.

## Example commands:

For test pulse runs:

cmsRun runPulseTree.py runtype=TP inputFiles=file:run284191_ls0001_streamDQM_mrg-c2f12-31-01.dat outputFile=test.root

For laser runs:

cmsRun runPulseTree.py runtype=Laser inputFiles=file:run284984_ls0003_streamDQM_mrg-c2f12-31-01.dat inputFiles=file:run284984_ls0004_streamDQM_mrg-c2f12-31-01.dat outputFile=test.root

For PhiSym stream:

cmsRun runPulseTree.py runtype=PhiSym inputFiles=/store/data/Run2016H/AlCaPhiSym/RAW/v1/000/282/814/00000/0A8037BA-7F8F-E611-B6D1-02163E0145FE.root outputFile=test.root
