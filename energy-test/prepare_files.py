import ROOT as R
import argparse

parser = argparse.ArgumentParser(description="Prepare the data.")
parser.add_argument('filename')
parser.add_argument('--nevents',default = -1, type = int)
args = parser.parse_args()
filename = args.filename
nevents = args.nevents

rootfile = R.TFile(filename+'.root', 'read')
roottree = rootfile.Get("DalitzEventList")
if nevents == -1:
    nevents = roottree.GetEntries()
output = open(filename+"_"+str(nevents)+".txt",'w')
print('Considering',nevents,'candidates')
count = 0
for entry in roottree:
    count +=1
    if count > nevents:
        break
    #create 4 vectors
    kaon = R.TLorentzVector()
    kaon.SetPxPyPzE(entry.Kp_Px,entry.Kp_Py,entry.Kp_Pz,entry.Kp_E)
    piplus = R.TLorentzVector()
    piplus.SetPxPyPzE(entry.pip_Px,entry.pip_Py,entry.pip_Pz,entry.pip_E)
    piminusa = R.TLorentzVector()
    piminusa.SetPxPyPzE(entry.pim0_Px,entry.pim0_Py,entry.pim0_Pz,entry.pim0_E)
    piminusb = R.TLorentzVector()
    piminusb.SetPxPyPzE(entry.pim1_Px,entry.pim1_Py,entry.pim1_Pz,entry.pim1_E)

    #create the relevant variables
    #we will want both pi+ pi- and kaon pi-, this will set 4 variables
    #fifth variable? highest mass K+,pi-,pi+ combination?
    kpipi_a = (kaon+piplus+piminusa).M()
    kpipi_b = (kaon+piplus+piminusb).M()

    piminus0 = R.TLorentzVector()
    piminus1 = R.TLorentzVector()

    if(kpipi_a > kpipi_b):
        piminus0 = piminusa
        piminus1 = piminusb
    else:
        piminus0 = piminusb
        piminus1 = piminusa

    output.write(str((piplus+piminus0).M2())+" "+str((piplus+piminus1).M2())+" "+str((kaon+piminus0).M2())+" "+str((kaon+piminus1).M2())+" "+str((kaon+piplus+piminus0).M2())+"\n")
    
    
