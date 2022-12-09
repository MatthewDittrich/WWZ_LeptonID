import ROOT as r
import math
from ROOT import TFile
from ROOT import TCanvas
from ROOT import TH1F
from ROOT import TH2F
from ROOT import TLegend
import csv
from array import array
import sys
import numpy as np
#-----------------------------------------------------Defining Some Histograms--------------------------------------------------------------------
e_pt_bin = [10,15,20,30,40,50,200]
m_pt_bin = [10,15,20,30,40,50,200]
e_eta_bin = [-2.5,-2.0,-1.566,-1.4442,-0.8,0,0.8,1.4442,1.566,2.0,2.5]
m_eta_bin = [-2.5,-1.1,0,1.1,2.5]

den_e_fine = TH2F("den_e_fine","Gen Level Electrons Fine",180,-2.5,2.5,180,0,200)
num_e_fine = TH2F("num_e_fine","Passed Electrons Fine",180,-2.5,2.5,180,0,200)
den_e = TH2F("den_e","Gen Level Electrons",10,array('d',e_eta_bin),6,array('d',e_pt_bin))
num_e = TH2F("num_e","Passed Electrons",10,array('d',e_eta_bin),6,array('d',e_pt_bin))

den_m_fine = TH2F("den_m_fine","Gen Level Muons Fine",180,-2.5,2.5,180,0,200)
num_m_fine = TH2F("num_m_fine","Passed Muons Fine",180,-2.5,2.5,180,0,200)
den_m = TH2F("den_m","Gen Level Muons",4,array('d',m_eta_bin),6,array('d',m_pt_bin))
num_m = TH2F("num_m","Passed Muons",4,array('d',m_eta_bin),6,array('d',m_pt_bin))

num_m_loose = TH2F("num_m_loose","Passed Muons Loose Id",4,array('d',m_eta_bin),6,array('d',m_pt_bin))
num_m_medium = TH2F("num_m_medium","Passed Muons Medium Id",4,array('d',m_eta_bin),6,array('d',m_pt_bin))
num_m_tight = TH2F("num_m_tight","Passed Muons Tight Id",4,array('d',m_eta_bin),6,array('d',m_pt_bin))

num_e_WWZ_loose = TH2F("num_e_WWZ_loose","Passed Electrons Loose",10,array('d',e_eta_bin),6,array('d',e_pt_bin))
num_e_Zlep = TH2F("num_e_Zlep","Passed Electrons Z Cand",10,array('d',e_eta_bin),6,array('d',e_pt_bin))
num_e_Wlep = TH2F("num_e_Wlep","Passed Electrons W Cand",10,array('d',e_eta_bin),6,array('d',e_pt_bin))
num_m_WWZ_loose = TH2F("num_m_WWZ_loose","Passed Muons Loose",4,array('d',m_eta_bin),6,array('d',m_pt_bin))
num_m_Zlep = TH2F("num_m_Zlep","Passed Muons Z Cand",4,array('d',m_eta_bin),6,array('d',m_pt_bin))
num_m_Wlep = TH2F("num_m_Wlep","Passed Muons W Cand",4,array('d',m_eta_bin),6,array('d',m_pt_bin))
#-------------------------------------------------------Grabbing the Root Files -------------------------------------------------------------------

f = r.TFile.Open("/home/matthew.dittrich/mc/RunIISummer20UL18NanoAODv9/WWZJetsTo4L2Nu_4F_TuneCP5_13TeV-amcatnlo-pythia8/NANOAODSIM/106X_upgrade2018_realistic_v16_L1v1-v2/2520000/E1EF60E7-36AB-9044-B0BC-31222489C489.root","Read")
t = f.Get("Events")

f_lep = 0
not_f_lep = 0
Total_events = 0
No_z_events = 0
Z_ee_events = 0
Z_mumu_events = 0
Z_tautau_events = 0
leps_p = 0
leps_w = 0
z_alone = 0
w_alone = 0
p_alone = 0
zw_alone = 0
zp_alone = 0
wp_alone = 0
zwp_alone = 0
#----------------------------------------------------------A Function or Two-----------------------------------------------------------------------
def delta_r_match(pt_gen,pt_r,eta_gen,eta_r,phi_gen,phi_r,loose_r,medium_r,tight_r):
    match_pt = []
    match_eta = []
    is_loose = []
    is_medium = []
    is_tight = []
    for h in range(len(pt_gen)):
        best_delta_r = 999999
        best_pt_r = 0
        loose = False
        medium = False
        tight = False
        for f in range(len(pt_r)):
            delta_r = math.sqrt(((eta_r[f] - eta_gen[h])**2)+((phi_r[f] - phi_gen[h])**2))
            if (delta_r < best_delta_r):
                best_delta_r = delta_r
                best_pt_r = pt_r[f]

                if (loose_r[f] == 1):
                    loose = True
                else:
                    loose = False

                if (medium_r[f] == 1):
                    medium = True
                else:
                    medium = False

                if (tight_r[f] == 1):
                    tight = True
                else:
                    tight = False

        if (best_delta_r <= 0.4) and (best_pt_r >= 10):
            match_pt.append(pt_gen[h])
            match_eta.append(eta_gen[h])
            is_loose.append(loose)
            is_medium.append(medium)
            is_tight.append(tight)
    return match_pt, match_eta, is_loose, is_medium, is_tight


#--------------------------------------------------------------------------------------------------------------------------------------------------
for i in range (t.GetEntries()):

    Total_events +=1 

    t.GetEntry(i)
    
    status = t.GenPart_status
    pdgId = t.GenPart_pdgId
    mother = t.GenPart_genPartIdxMother
    pt = t.GenPart_pt
    mass = t.GenPart_mass
    eta = t.GenPart_eta
    phi = t.GenPart_phi
    pdgId_e_r = t.Electron_pdgId
    eta_e_r = t.Electron_eta
    phi_e_r = t.Electron_phi
    pt_e_r = t.Electron_pt 
    pdgId_m_r = t.Muon_pdgId
    eta_m_r = t.Muon_eta
    phi_m_r = t.Muon_phi
    pt_m_r = t.Muon_pt 
    mu_tight = t.Muon_tightId
    mu_medium = t.Muon_mediumId
    mu_loose = t.Muon_looseId
    mu_pfiso = t.Muon_pfIsoId

    w_index = -3
    aw_index = -3
    z_index = -3

    w_status = -3
    aw_status = -3
    z_status = -3

    w_lep_list = []
    z_lep_list = []
    part_lep_list = []
    lepton_counter = 0
    
    gen_e_pt = []
    gen_e_phi = []
    gen_e_eta = []
    gen_m_pt = []
    gen_m_phi = []
    gen_m_eta = []

#----------------------------------------------Grabbing the Correct W/Z Boson indices ----------------------------------------------------------------------------------------
    for j in range (len(t.GenPart_pdgId)):

        if (mother[j] == 1):
            print("We have a child of 1. This is an issue.")
        if (pdgId[j] == 23) and (status[j] > z_status):
            z_index = j
            z_status = status[j]
        if (pdgId[j] == 24) and (status[j] > w_status):
            w_index = j
            w_status = status[j]
        if (pdgId[j] == -24) and (status[j] > aw_status):
            aw_index = j
            aw_status = status[j]

    if z_index == -3:
        No_z_events += 1

#------------------------------------------------------------Grabbing 4 Leps including the Taus-----------------------------------------------------------------------------

    for k in range(len(t.GenPart_pdgId)):
        if ((pdgId[k] == 11) or (pdgId[k]== -11)):
            if ((mother[k] == w_index) or (mother[k] == aw_index)):
                lepton_counter += 1
                w_lep_list.append(k)
                den_e.Fill(eta[k],pt[k])
                den_e_fine.Fill(eta[k],pt[k])
                gen_e_pt.append(pt[k])
                gen_e_phi.append(phi[k])
                gen_e_eta.append(eta[k])
            if (mother[k] == z_index):
                lepton_counter += 1
                z_lep_list.append(k)
                den_e.Fill(eta[k],pt[k])
                den_e_fine.Fill(eta[k],pt[k])
                gen_e_pt.append(pt[k])
                gen_e_phi.append(phi[k])
                gen_e_eta.append(eta[k])
            if (mother[k] == 0):
                lepton_counter += 1
                part_lep_list.append(k)
                den_e.Fill(eta[k],pt[k])
                den_e_fine.Fill(eta[k],pt[k])
                gen_e_pt.append(pt[k])
                gen_e_phi.append(phi[k])
                gen_e_eta.append(eta[k])
        if ((pdgId[k] == 13) or (pdgId[k] == -13)):
            if ((mother[k] == w_index) or (mother[k] == aw_index)):
                lepton_counter += 1
                w_lep_list.append(k)
                den_m.Fill(eta[k],pt[k])
                den_m_fine.Fill(eta[k],pt[k])
                gen_m_pt.append(pt[k])
                gen_m_phi.append(phi[k])
                gen_m_eta.append(eta[k])
            if (mother[k] == z_index):
                lepton_counter += 1
                z_lep_list.append(k)
                den_m.Fill(eta[k],pt[k])
                den_m_fine.Fill(eta[k],pt[k])
                gen_m_pt.append(pt[k])
                gen_m_phi.append(phi[k])
                gen_m_eta.append(eta[k])
            if (mother[k] == 0):
                lepton_counter += 1
                part_lep_list.append(k)
                den_m.Fill(eta[k],pt[k])
                den_m_fine.Fill(eta[k],pt[k])
                gen_m_pt.append(pt[k])
                gen_m_phi.append(phi[k])
                gen_m_eta.append(eta[k])
        if ((pdgId[k] == 15) or (pdgId[k] == -15)):
            if ((mother[k] == w_index) or (mother[k] == aw_index)):
                lepton_counter += 1
                w_lep_list.append(k)
            if (mother[k] == z_index):
                lepton_counter += 1
                z_lep_list.append(k)
            if (mother[k] == 0):
                lepton_counter += 1
                part_lep_list.append(k)
#----------------------------------------------------------------Getting the Reconstructed Leptons--------------------------------------------------------------------

    for j in range(len(gen_e_pt)):
        best_dr = 999999
        best_dr_index = 999999
        if (len(pt_e_r) > 0):
            for h in range(len(pt_e_r)):
                delta_r = math.sqrt(((eta_e_r[h] - gen_e_eta[j])**2)+((phi_e_r[h] - gen_e_phi[j])**2))
                if (delta_r < best_dr):
                    best_dr = delta_r
                    best_dr_index = h
            if ((pt_e_r[best_dr_index] >= 10) and (best_dr <= 0.4)):
                num_e.Fill(gen_e_eta[j],gen_e_pt[j])
                num_e_fine.Fill(gen_e_eta[j],gen_e_pt[j]) 
                if ((eta_e_r[best_dr_index] < 2.5) and (t.Electron_mvaFall17V2noIso_WPL[best_dr_index] == 1) and (t.Electron_pfRelIso03_all[best_dr_index] < 0.4)):
                    num_e_WWZ_loose.Fill(gen_e_eta[j],gen_e_pt[j])
                    if ((t.Electron_sip3d[best_dr_index] < 4) and (t.Electron_pfRelIso03_all[best_dr_index] < 0.2)):
                        num_e_Zlep.Fill(gen_e_eta[j],gen_e_pt[j])
                        if (t.Electron_mvaFall17V2Iso_WP90[best_dr_index] == 1):
                            num_e_Wlep.Fill(gen_e_eta[j],gen_e_pt[j])

            
    for j in range(len(gen_m_pt)):
        best_dr = 999999
        best_dr_index = 999999
        if (len(pt_m_r) > 0):
            for h in range(len(pt_m_r)):
                delta_r = math.sqrt(((eta_m_r[h] - gen_m_eta[j])**2)+((phi_m_r[h] - gen_m_phi[j])**2))
                if (delta_r < best_dr):
                    best_dr = delta_r
                    best_dr_index = h
            if ((pt_m_r[best_dr_index] >= 10) and (best_dr <= 0.4)):
                num_m.Fill(gen_m_eta[j],gen_m_pt[j])
                num_m_fine.Fill(gen_m_eta[j],gen_m_pt[j])
                if (mu_tight[best_dr_index] == 1):
                    num_m_tight.Fill(gen_m_eta[j],gen_m_pt[j])
                if (mu_medium[best_dr_index] == 1):
                    num_m_medium.Fill(gen_m_eta[j],gen_m_pt[j])
                if (mu_loose[best_dr_index] == 1):
                    num_m_loose.Fill(gen_m_eta[j],gen_m_pt[j])
                    if ((eta_m_r[best_dr_index] < 2.4) and (mu_pfiso[best_dr_index] >= 1)):
                        num_m_WWZ_loose.Fill(gen_m_eta[j],gen_m_pt[j])
                        if ((t.Muon_sip3d[best_dr_index] < 4) and (t.Muon_pfIsoId[best_dr_index] >= 2)):
                            num_m_Zlep.Fill(gen_m_eta[j],gen_m_pt[j])
                            if (t.Muon_pfIsoId[best_dr_index] >= 3):
                                num_m_Wlep.Fill(gen_m_eta[j],gen_m_pt[j])       
                

#-----------------------------------------------------------------Counting Z Events and Parton/W Events ----------------------------------------------------------------
    if len(z_lep_list) > 0:

        if len(z_lep_list) != 2:
            print(len(z_lep_list),i)
        if (pdgId[z_lep_list[1]] == 11) or (pdgId[z_lep_list[1]] == -11):
            Z_ee_events += 1
        if (pdgId[z_lep_list[1]] == 13) or (pdgId[z_lep_list[1]] == -13):
            Z_mumu_events += 1
        if (pdgId[z_lep_list[1]] == 15) or (pdgId[z_lep_list[1]] == -15):
            Z_tautau_events += 1

    if len(part_lep_list) > 0:
        leps_p += 1
    if len(w_lep_list) > 0:
        leps_w += 1

#------------------------------------------------------Where did the 4 leps come from------------------------------------------------------------------------------------
    if (len(w_lep_list) == 0) and (len(z_lep_list) == 0) and (len(part_lep_list) > 0):
        p_alone += 1
    if (len(w_lep_list) > 0) and (len(z_lep_list) == 0) and (len(part_lep_list) == 0):
        w_alone += 1
    if (len(w_lep_list) == 0) and (len(z_lep_list) > 0) and (len(part_lep_list) ==  0):
        z_alone += 1
    if (len(w_lep_list) == 0) and (len(z_lep_list) > 0) and (len(part_lep_list) > 0):
        zp_alone += 1
    if (len(w_lep_list) > 0) and (len(z_lep_list) == 0) and (len(part_lep_list) > 0):
        wp_alone += 1
    if (len(w_lep_list) > 0) and (len(z_lep_list) > 0) and (len(part_lep_list) == 0):
        zw_alone += 1
    if (len(w_lep_list) > 0) and (len(z_lep_list) > 0) and (len(part_lep_list) > 0):
        zwp_alone += 1
#--------------------------------------------------------------Make Sure we actually have 4 leptons -------------------------------------------------------------------------

    if lepton_counter != 4:
        print("!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!")
    elif lepton_counter == 4:
        f_lep += 1

print("Total Events: ",Total_events)
print("Events with 4 Leps: ", f_lep)
print("No Z Events: ", No_z_events)
print("Z to ee Events: ", Z_ee_events)
print("Z to mumu Events: ", Z_mumu_events)
print("Z to tautau Events: ", Z_tautau_events)
print("Leps from W: ", leps_w)
print("Leps from partons: ",leps_p)
print("WZP: ", zwp_alone)
print("W: ", w_alone)
print("Z: ", z_alone)
print("P: ", p_alone)
print("WZ: ", zw_alone)
print("WP: ", wp_alone)
print("ZP: ", zp_alone)
print("Events in e num: ",num_e.GetEntries())
print("Events in e den: ",den_e.GetEntries())
print("Events in mu num: ",num_m.GetEntries())
print("Events in mu den: ",den_m.GetEntries())
#---------------------------------------------------------------------------------------------------------------------------------------------------------------------------

ratio_e = num_e.Clone("ratio_e")
ratio_e.Divide(den_e)
#ratio_e_fine = num_e_fine.Clone("ratio_e_fine")
#ratio_e_fine.Divide(den_e_fine)

ratio_m = num_m.Clone("ratio_m")
ratio_m.Divide(den_m)
ratio_m_fine = num_m_fine.Clone("ratio_m_fine")
ratio_m_fine.Divide(den_m_fine)

ratio_m_loose = num_m_loose.Clone("ratio_m_loose")
ratio_m_loose.Divide(den_m)
ratio_m_medium = num_m_medium.Clone("ratio_m_medium")
ratio_m_medium.Divide(den_m)
ratio_m_tight = num_m_tight.Clone("ratio_m_tight")
ratio_m_tight.Divide(den_m)

ratio_e_WWZ_loose = num_e_WWZ_loose.Clone("ratio_e_WWZ_loose")
ratio_e_WWZ_loose.Divide(den_e)
ratio_e_Zlep = num_e_Zlep.Clone("ratio_e_Zlep")
ratio_e_Zlep.Divide(den_e)
ratio_e_Wlep = num_e_Wlep.Clone("ratio_e_Wlep")
ratio_e_Wlep.Divide(den_e)

ratio_m_WWZ_loose = num_m_WWZ_loose.Clone("ratio_m_WWZ_loose")
ratio_m_WWZ_loose.Divide(den_m)
ratio_m_Zlep = num_m_Zlep.Clone("ratio_m_Zlep")
ratio_m_Zlep.Divide(den_m)
ratio_m_Wlep = num_m_Wlep.Clone("ratio_m_Wlep")
ratio_m_Wlep.Divide(den_m)

#ratio_e_loose = num_e_loose.Clone("ratio_e_loose")
#ratio_e_loose.Divide(den_e)
#ratio_e_medium = num_e_medium.Clone("ratio_e_medium")
#ratio_e_medium.Divide(den_e)
#ratio_e_tight = num_e_tight.Clone("ratio_e_tight")
#ratio_e_tight.Divide(den_e)

o = r.TFile("lepton_id_study.root", "recreate");
o.cd()
num_e.Write()
den_e.Write()
ratio_e.Write()
#num_e_fine.Write()
#den_e_fine.Write()
#ratio_e_fine.Write()
num_m.Write()
den_m.Write()
ratio_m.Write()
#num_m_fine.Write()
#den_m_fine.Write()
#ratio_m_fine.Write()
#ratio_e_loose.Write()
#ratio_e_medium.Write()
#ratio_e_tight.Write()
ratio_m_loose.Write()
ratio_m_medium.Write()
ratio_m_tight.Write()
ratio_e_WWZ_loose.Write()
ratio_e_Zlep.Write()
ratio_e_Wlep.Write()
ratio_m_WWZ_loose.Write()
ratio_m_Zlep.Write()
ratio_m_Wlep.Write()
num_m_Wlep.Write()
num_m_Zlep.Write()
num_e_Wlep.Write()
num_e_Zlep.Write()
num_e_WWZ_loose.Write()
num_m_WWZ_loose.Write()
sys.exit()






