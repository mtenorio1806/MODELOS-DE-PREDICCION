import sys, getopt, os, errno
import numpy as np
from Bio import SeqIO
import sklearn
print(sklearn.__version__)
try:
    import joblib
except ImportError:
    print('libreria joblib no instalada')
def predecir(matriz,modelo,evaluar):
  model=joblib.load('/vault2/homehpc/dtalero/SCRIPT/alelos/PR/modelos/'+modelo+'.sav')
  #evalu=[0,0,0,0,0,0,0,0,0]
  print(evaluar)
  arreglo=np.array([evaluar])
  resultado=model.predict(arreglo)
  print(resultado[0])
def matriz(csv,salida,mod,tag):
  modelo={'TIGECICLINA':['Klebsiella_pneumoniae_KpnE_6','Klebsiella_pneumoniae_KpnF_5','Escherichia_coli_soxS_with_mutation_conferring_antibiotic_resistance_cp1_17'],'TETRACICLINA':['OXA-9_1','tet(B)_1','sul2_1','tetR_4','tetR_1','tet(A)_cp1_2','tet(A)_2','Mycobacterium_tuberculosis_rpoB_mutants_conferring_resistance_to_rifampicin_1','Escherichia_coli_soxS_with_mutation_conferring_antibiotic_resistance_cp1_17'],'AMIKACINA':["Escherichia_coli_ampC_beta-lactamase_13", "arr-2_1", "BRP(MBL)_2", "AAC(6')-Ib7_5", "adeF_21", "armA_1", "AAC(6')-Il_1", "Salmonella_serovars_gyrB_conferring_resistance_to_fluoroquinolone_1"],'CIPROFLOXACINA':["Klebsiella_pneumoniae_KpnE_13", "marA_4", "Escherichia_coli_marR_mutant_conferring_antibiotic_resistance_4", "Staphylococcus_aureus_fusA_with_mutation_conferring_resistance_to_fusidic_acid_10", "CRP_1", "CRP_3", "CRP_4", "mdtE_1", "gadW_12", "gadX_5", "TEM-1_1", "AAC(3)-IId_2", "AAC(3)-IId_cp1_2", "Escherichia_coli_mdfA_1", "Escherichia_coli_mdfA_3", "Escherichia_coli_mdfA_4", "kdpE_17", "Klebsiella_pneumoniae_OmpK37_cp2_142", "Klebsiella_pneumoniae_OmpK37_cp2_223", "Klebsiella_pneumoniae_OmpK37_cp2_252", "Klebsiella_pneumoniae_OmpK37_cp2_270", "Klebsiella_pneumoniae_OmpK37_cp2_296", "Klebsiella_pneumoniae_OmpK37_cp3_296", "Klebsiella_pneumoniae_OmpK37_cp3_358", "Klebsiella_pneumoniae_OmpK37_cp5_331", "YojI_2", "baeR_2", "baeS_26", "mdtC_13", "mdtA_102", "mdtA_28", "mdtA_54", "mdtA_79", "ugd_8", "rosB_33", "Escherichia_coli_acrR_with_mutation_conferring_multidrug_antibiotic_resistance_14", "Escherichia_coli_acrA_4", "Escherichia_coli_ampH_beta-lactamase_26", "TolC_1", "TolC_3", "AcrS_19", "AcrS_42", "AcrE_29", "AcrF_9", "Enterobacter_cloacae_acrA_1", "Enterobacter_cloacae_acrA_14", "H-NS_cp1_20", "Haemophilus_influenzae_PBP3_conferring_resistance_to_beta-lactam_antibiotics_14", "Escherichia_coli_soxS_with_mutation_conferring_antibiotic_resistance_25", "mdtM_45", "mdtM_cp1_58", "evgA_11", "evgA_4", "acrD_7", "mdtH_6", "mdtG_22", "emrB_9", "Escherichia_coli_soxR_with_mutation_conferring_antibiotic_resistance_14", "mdtN_20", "eptA_52", "PmrF_10", "PmrF_24", "cpxA_13", "Clostridioides_difficile_gyrB_conferring_resistance_to_fluoroquinolone_8", "sdiA_18", "sdiA_2", "sdiA_6", "Escherichia_coli", "emrE_1", "MdtK_13", "msbA_4", "msbA_5", "Mycobacterium_tuberculosis_rpoB_mutants_conferring_resistance_to_rifampicin_1", "dfrA12_cp1_1", "qacH_1", "sul1_2", "mphA_1", "Escherichia_coli_EF-Tu_mutants_conferring_resistance_to_Pulvomycin_cp1_2", "Clostridioides_difficile_gyrA_conferring_resistance_to_fluoroquinolones_18", "Mycobacterium_tuberculosis_thyA_with_mutation_conferring_resistance_to_para-aminosalicylic_acid_13", "Mycobacterium_tuberculosis_katG_mutations_conferring_resistance_to_isoniazid_5", "Escherichia_coli_UhpT_with_mutation_conferring_resistance_to_fosfomycin_6", "CTX-M-15_1", "CTX-M-15_cp1_1", "dfrA1_2", "SAT-2_1", "aadA_4", "aadA_7", "APH(3'')-Ib_4", "Escherichia_coli_gyrA_conferring_resistance_to_fluoroquinolones_10", "Escherichia_coli_gyrA_conferring_resistance_to_fluoroquinolones_11", "Escherichia_coli_gyrA_conferring_resistance_to_fluoroquinolones_13", "Escherichia_coli_gyrA_conferring_resistance_to_fluoroquinolones_15", "Escherichia_coli_gyrA_conferring_resistance_to_fluoroquinolones_17", "Escherichia_coli_gyrA_conferring_resistance_to_fluoroquinolones_9", "Escherichia_coli_parC_conferring_resistance_to_fluoroquinolone_10", "Escherichia_coli_parC_conferring_resistance_to_fluoroquinolone_11", "Escherichia_coli_parC_conferring_resistance_to_fluoroquinolone_6", "arr-2_1", "NDM-1_1", "QnrB1_1", "AAC(6')-Ib-cr_1", "dfrA14_1", "dfrA14_3", "dfrA14_cp1_1", "dfrA14_cp1_3", "SHV-66_1", "QnrS1_1", "AAC(6')-Ib7_5", "catI_3", "AAC(3)-IIe_1", "AAC(3)-IIe_cp1_1", "Klebsiella_pneumoniae_KpnG_11", "adeF_17", "adeF_18", "oqxA_1", "ramA_4", "Klebsiella_pneumoniae", "acrA_1", "CTX-M-55_1", "ErmB_cp1_1", "catII_from_Escherichia_coli_K-12_cp1_1", "OXA-9_1", "KPC-3_1", "AAC(6')-Ib4_1", "SHV-30_1"],'COLISTINA':["Klebsiella_pneumoniae_KpnE_1", "Klebsiella_pneumoniae_KpnE_10", "Klebsiella_pneumoniae_KpnE_2", "Klebsiella_pneumoniae_KpnE_3", "Klebsiella_pneumoniae_KpnE_4", "Klebsiella_pneumoniae_KpnE_5", "Klebsiella_pneumoniae_KpnE_6", "Klebsiella_pneumoniae_KpnE_7", "Klebsiella_pneumoniae_KpnE_8", "Klebsiella_pneumoniae_KpnE_9", "Klebsiella_pneumoniae_KpnF_1", "Klebsiella_pneumoniae_KpnF_2", "Klebsiella_pneumoniae_KpnF_3", "Klebsiella_pneumoniae_KpnF_4", "Klebsiella_pneumoniae_KpnF_5", "Klebsiella_pneumoniae_KpnF_6", "Klebsiella_pneumoniae_KpnF_7", "Klebsiella_pneumoniae_KpnF_8", "Klebsiella_pneumoniae_KpnF_9", "marA_1", "marA_2", "marA_3", "Escherichia_coli_marR_mutant_conferring_antibiotic_resistance_1", "Escherichia_coli_marR_mutant_conferring_antibiotic_resistance_2", "Escherichia_coli_marR_mutant_conferring_antibiotic_resistance_3", "Staphylococcus_aureus_fusA_with_mutation_conferring_resistance_to_fusidic_acid_1", "Staphylococcus_aureus_fusA_with_mutation_conferring_resistance_to_fusidic_acid_2", "Staphylococcus_aureus_fusA_with_mutation_conferring_resistance_to_fusidic_acid_3", "Staphylococcus_aureus_fusA_with_mutation_conferring_resistance_to_fusidic_acid_4", "Staphylococcus_aureus_fusA_with_mutation_conferring_resistance_to_fusidic_acid_cp1_5", "Staphylococcus_aureus_fusA_with_mutation_conferring_resistance_to_fusidic_acid_cp1_6", "mdtE_1", "mdtF_1", "gadW_1", "gadX_1", "gadX_2", "TEM-1_cp1_1", "TEM-1_cp2_1", "TEM-1_cp3_1", "AAC(3)-IId_1", "Escherichia_coli_mdfA_1", "Escherichia_coli_mdfA_2", "kdpE_1", "kdpE_2", "kdpE_3", "kdpE_4", "kdpE_5", "kdpE_6", "kdpE_7", "kdpE_8", "Klebsiella_pneumoniae_OmpK37_12", "Klebsiella_pneumoniae_OmpK37_15", "Klebsiella_pneumoniae_OmpK37_18", "Klebsiella_pneumoniae_OmpK37_21", "Klebsiella_pneumoniae_OmpK37_3", "Klebsiella_pneumoniae_OmpK37_32", "Klebsiella_pneumoniae_OmpK37_34", "Klebsiella_pneumoniae_OmpK37_39", "Klebsiella_pneumoniae_OmpK37_43", "Klebsiella_pneumoniae_OmpK37_5", "Klebsiella_pneumoniae_OmpK37_55", "Klebsiella_pneumoniae_OmpK37_cp1_19", "Klebsiella_pneumoniae_OmpK37_cp1_21", "Klebsiella_pneumoniae_OmpK37_cp1_22", "Klebsiella_pneumoniae_OmpK37_cp1_24", "Klebsiella_pneumoniae_OmpK37_cp1_28", "Klebsiella_pneumoniae_OmpK37_cp1_34", "Klebsiella_pneumoniae_OmpK37_cp1_37", "Klebsiella_pneumoniae_OmpK37_cp1_42", "Klebsiella_pneumoniae_OmpK37_cp1_47", "Klebsiella_pneumoniae_OmpK37_cp1_65", "Klebsiella_pneumoniae_OmpK37_cp1_69", "Klebsiella_pneumoniae_OmpK37_cp2_14", "Klebsiella_pneumoniae_OmpK37_cp2_33", "Klebsiella_pneumoniae_OmpK37_cp2_34", "Klebsiella_pneumoniae_OmpK37_cp2_42", "Klebsiella_pneumoniae_OmpK37_cp2_45", "Klebsiella_pneumoniae_OmpK37_cp2_67", "Klebsiella_pneumoniae_OmpK37_cp2_73", "Klebsiella_pneumoniae_OmpK37_cp2_75", "Klebsiella_pneumoniae_OmpK37_cp2_87", "Klebsiella_pneumoniae_OmpK37_cp3_52", "Klebsiella_pneumoniae_OmpK37_cp3_76", "Klebsiella_pneumoniae_OmpK37_cp3_77", "Klebsiella_pneumoniae_OmpK37_cp3_79", "Klebsiella_pneumoniae_OmpK37_cp3_84", "Klebsiella_pneumoniae_OmpK37_cp3_86", "Klebsiella_pneumoniae_OmpK37_cp4_56", "Klebsiella_pneumoniae_OmpK37_cp4_76", "Klebsiella_pneumoniae_OmpK37_cp4_79", "Klebsiella_pneumoniae_OmpK37_cp4_82", "Klebsiella_pneumoniae_OmpK37_cp4_87", "Klebsiella_pneumoniae_OmpK37_cp4_89", "_Klebsiella_pneumoniae_OmpK37_cp4_90", "Klebsiella_pneumoniae_OmpK37_cp4_93", "Klebsiella_pneumoniae_OmpK37_cp5_87", "Klebsiella_pneumoniae_OmpK37_cp5_89", "Klebsiella_pneumoniae_OmpK37_cp5_92", "Klebsiella_pneumoniae_OmpK37_cp6_89", "YojI_1", "YojI_2", "YojI_3", "YojI_4", "baeR_1", "baeR_3", "baeR_4", "baeR_5", "baeR_6", "baeR_7", "baeR_8", "baeR_9", "baeS_10", "baeS_11", "baeS_3", "baeS_4", "baeS_5", "baeS_6", "baeS_7", "baeS_8", "baeS_9", "mdtC_1", "mdtC_2", "mdtC_3", "mdtC_4", "mdtC_5", "mdtC_6", "mdtC_7", "mdtB_1", "mdtB_2", "mdtB_3", "mdtB_4", "mdtB_5", "mdtB_8", "mdtA_1", "mdtA_10", "mdtA_11", "mdtA_12", "mdtA_2", "mdtA_4", "mdtA_5", "mdtA_6", "mdtA_7", "mdtA_8", "mdtA_9", "vatB_2", "ugd_13", "ugd_16", "ugd_17", "ugd_19", "ugd_24", "ugd_26", "ugd_29", "ugd_6", "ugd_7", "ugd_9", "ugd_cp1_19", "ugd_cp1_6", "rosA_1", "rosA_2", "rosA_3", "rosA_4", "rosA_5", "rosB_1", "rosB_2", "rosB_3", "rosB_4", "rosB_5", "rosB_6", "Escherichia_coli_acrR_with_mutation_conferring_multidrug_antibiotic_resistance_2", "Escherichia_coli_acrR_with_mutation_conferring_multidrug_antibiotic_resistance_4", "Escherichia_coli_acrR_with_mutation_conferring_multidrug_antibiotic_resistance_5", "Escherichia_coli_acrR_with_mutation_conferring_multidrug_antibiotic_resistance_7", "Escherichia_coli_acrR_with_mutation_conferring_multidrug_antibiotic_resistance_cp1_9", "Escherichia_coli_acrA_1", "Escherichia_coli_acrA_2", "Escherichia_coli_acrA_3", "Escherichia_coli_acrA_4", "Escherichia_coli_acrA_5", "acrB_2", "acrB_3", "acrB_4", "acrB_6", "acrB_9", "acrB_cp1_7", "acrB_cp1_8", "Escherichia_coli_ampH_beta-lactamase_1", "Escherichia_coli_ampH_beta-lactamase_10", "Escherichia_coli_ampH_beta-lactamase_2", "Escherichia_coli_ampH_beta-lactamase_3", "Escherichia_coli_ampH_beta-lactamase_4", "Escherichia_coli_ampH_beta-lactamase_5", "Escherichia_coli_ampH_beta-lactamase_6", "Escherichia_coli_ampH_beta-lactamase_7", "Escherichia_coli_ampH_beta-lactamase_8", "Escherichia_coli_ampH_beta-lactamase_9", "TolC_1", "TolC_10", "TolC_11", "TolC_3", "TolC_5", "TolC_7", "TolC_8", "TolC_9", "TolC_cp1_12", "bacA_1", "bacA_2", "bacA_3", "bacA_4", "AcrS_1", "AcrS_2", "AcrS_3", "AcrS_4", "AcrE_1", "AcrE_2", "AcrE_3", "AcrE_4", "AcrE_5", "AcrF_3", "AcrF_4", "Enterobacter_cloacae_acrA_1", "Enterobacter_cloacae_acrA_2", "Enterobacter_cloacae_acrA_3", "Enterobacter_cloacae_acrA_4", "H-NS_13", "H-NS_16", "H-NS_2", "H-NS_4", "H-NS_7", "H-NS_8", "H-NS_9", "H-NS_cp1_12", "H-NS_cp1_17", "H-NS_cp1_20", "H-NS_cp2_14", "H-NS_cp2_22", "H-NS_cp2_26", "H-NS_cp2_30", "H-NS_cp3_26", "Haemophilus_influenzae_PBP3_conferring_resistance_to_beta-lactam_antibiotics_1", "Haemophilus_influenzae_PBP3_conferring_resistance_to_beta-lactam_antibiotics_2", "Haemophilus_influenzae_PBP3_conferring_resistance_to_beta-lactam_antibiotics_3", "Haemophilus_influenzae_PBP3_conferring_resistance_to_beta-lactam_antibiotics_4", "Haemophilus_influenzae_PBP3_conferring_resistance_to_beta-lactam_antibiotics_5", "Haemophilus_influenzae_PBP3_conferring_resistance_to_beta-lactam_antibiotics_6", "dfrA3_1", "dfrA3_2", "dfrA3_3", "Escherichia_coli_soxS_with_mutation_conferring_antibiotic_resistance_1", "Escherichia_coli_soxS_with_mutation_conferring_antibiotic_resistance_10", "Escherichia_coli_soxS_with_mutation_conferring_antibiotic_resistance_11", "Escherichia_coli_soxS_with_mutation_conferring_antibiotic_resistance_2", "Escherichia_coli_soxS_with_mutation_conferring_antibiotic_resistance_5", "Escherichia_coli_soxS_with_mutation_conferring_antibiotic_resistance_6", "Escherichia_coli_soxS_with_mutation_conferring_antibiotic_resistance_7", "Escherichia_coli_soxS_with_mutation_conferring_antibiotic_resistance_cp1_15", "Escherichia_coli_soxS_with_mutation_conferring_antibiotic_resistance_cp1_16", "Escherichia_coli_soxS_with_mutation_conferring_antibiotic_resistance_cp1_17", "Escherichia_coli_soxS_with_mutation_conferring_antibiotic_resistance_cp1_2", "Escherichia_coli_soxS_with_mutation_conferring_antibiotic_resistance_cp1_4", "Escherichia_coli_soxS_with_mutation_conferring_antibiotic_resistance_cp1_8mdtM_1", "mdtM_2", "mdtM_3", "mdtM_cp1_4", "emrY_1", "emrY_2", "emrK_1", "emrK_2", "evgA_1", "evgA_3", "evgA_cp1_2", "evgA_cp1_4", "evgA_cp2_4", "evgS_1", "Escherichia_coli_PtsI_with_mutation_conferring_resistance_to_fosfomycin_1", "Escherichia_coli_PtsI_with_mutation_conferring_resistance_to_fosfomycin_2", "Escherichia_coli_PtsI_with_mutation_conferring_resistance_to_fosfomycin_3", "acrD_1", "acrD_2", "mdtH_1", "mdtH_2", "mdtH_3", "mdtH_4", "mdtH_5", "mdtH_6", "mdtH_7", "mdtH_8", "mdtG_1", "mdtG_10", "mdtG_3", "mdtG_4", "mdtG_5", "mdtG_6", "mdtG_7", "mdtG_8", "mdtG_9", "mdtG_cp1_11", "mdtG_cp1_12", "mdtG_cp1_4", "emrR_1", "emrR_2", "emrR_3", "emrR_4", "emrR_5", "emrR_6", "emrR_7", "emrA_1", "emrA_2", "emrA_4", "emrA_5", "emrA_6", "emrA_cp1_3", "emrA_cp1_7", "emrB_1", "emrB_2", "emrB_3", "Escherichia_coli_soxR_with_mutation_conferring_antibiotic_resistance_1", "Escherichia_coli_soxR_with_mutation_conferring_antibiotic_resistance_2", "Escherichia_coli_soxR_with_mutation_conferring_antibiotic_resistance_4", "Escherichia_coli_soxR_with_mutation_conferring_antibiotic_resistance_5", "Escherichia_coli_soxR_with_mutation_conferring_antibiotic_resistance_6", "Escherichia_coli_soxR_with_mutation_conferring_antibiotic_resistance_8", "mdtP_1,mdtO_2", "mdtN_1", "mdtN_2", "mdtN_3", "eptA_1", "eptA_2", "eptA_3", "Escherichia_coli_ampC_beta-lactamase_1", "Escherichia_coli_ampC_beta-lactamase_cp1_2", "PmrF_2", "PmrF_3", "PmrF_4", "PmrF_5", "PmrF_6", "PmrF_7", "cpxA_1", "cpxA_2", "cpxA_3", "cpxA_4", "cpxA_5", "cpxA_6", "Clostridioides_difficile_gyrB_conferring_resistance_to_fluoroquinolone_2", "Clostridioides_difficile_gyrB_conferring_resistance_to_fluoroquinolone_3", "sdiA_1", "sdiA_2", "sdiA_3", "sdiA_4", "sdiA_5", "sdiA_6", "Escherichia_coli_emrE_1", "Escherichia_coli_emrE_2", "Escherichia_coli_emrE_3", "Escherichia_coli_emrE_4", "Escherichia_coli_emrE_5", "Escherichia_coli_emrE_6", "MdtK_3", "MdtK_4", "MdtK_5", "msbA_1", "msbA_2", "msbA_4", "msbA_5", "msbA_6", "msbA_7", "Mycobacterium_tuberculosis_rpoB_mutants_conferring_resistance_to_rifampicin_1", "Mycobacterium_tuberculosis_rpoB_mutants_conferring_resistance_to_rifampicin_2", "Mycobacterium_tuberculosis_rpoB_mutants_conferring_resistance_to_rifampicin_3", "Mycobacterium_tuberculosis_rpoB_mutants_conferring_resistance_to_rifampicin_4", "macB_1", "macB_2", "macB_3", "macB_4", "macB_6", "macB_7", "macB_8", "dfrA12_1", "aadA2_1", "aadA2_2", "aadA2_cp2_2", "qacH_1", "qacH_3", "qacH_4", "qacH_cp1_1", "qacH_cp2_1", "qacH_cp2_4", "sul1_1", "sul1_2", "sul1_cp1_2", "sul1_cp2_2", "mphA_1", "Escherichia_coli_EF-Tu_mutants_conferring_resistance_to_Pulvomycin_1", "Escherichia_coli_EF-Tu_mutants_conferring_resistance_to_Pulvomycin_2", "Escherichia_coli_EF-Tu_mutants_conferring_resistance_to_Pulvomycin_3", "Escherichia_coli_EF-Tu_mutants_conferring_resistance_to_Pulvomycin_cp1_1", "Escherichia_coli_EF-Tu_mutants_conferring_resistance_to_Pulvomycin_cp1_2", "Escherichia_coli_EF-Tu_mutants_conferring_resistance_to_Pulvomycin_cp1_3", "tet(A)_2", "tet(A)_3", "Clostridioides_difficile_gyrA_conferring_resistance_to_fluoroquinolones_2", "Clostridioides_difficile_gyrA_conferring_resistance_to_fluoroquinolones_5", "Clostridioides_difficile_gyrA_conferring_resistance_to_fluoroquinolones_6", "Mycobacterium_tuberculosis_thyA_with_mutation_conferring_resistance_to_para-aminosalicylic_acid_1", "Mycobacterium_tuberculosis_thyA_with_mutation_conferring_resistance_to_para-aminosalicylic_acid_2", "Escherichia_coli_GlpT_with_mutation_conferring_resistance_to_fosfomycin_1", "Escherichia_coli_GlpT_with_mutation_conferring_resistance_to_fosfomycin_2", "Mycobacterium_tuberculosis_katG_mutations_conferring_resistance_to_isoniazid_1", "Mycobacterium_tuberculosis_katG_mutations_conferring_resistance_to_isoniazid_2", "Mycobacterium_tuberculosis_katG_mutations_conferring_resistance_to_isoniazid_cp1_2", "Escherichia_coli_UhpT_with_mutation_conferring_resistance_to_fosfomycin_1", "Escherichia_coli_UhpT_with_mutation_conferring_resistance_to_fosfomycin_2", "Escherichia_coli_UhpT_with_mutation_conferring_resistance_to_fosfomycin_3", "tetR_1", "tetR_2", "tetR_3", "tetR_5", "tetR_6", "tetR_cp1_3", "tetR_cp1_9", "mphB_1", "mphB_2", "mphB_5", "Escherichia_coli_ampC1_beta-lactamase_1", "arnA_2", "arnA_4", "arnA_5", "arnA_6", "arnA_7", "arnA_8", "CTX-M-15_1", "CTX-M-15_cp1_1", "CTX-M-15_cp2_1", "dfrA1_cp1_1", "SAT-2_1", "aadA_1", "floR_1", "floR_2", "APH(6)-Id_2", "APH(6)-Id_cp2_1", "sul2_1", "sul2_cp1_1", "Escherichia_coli_gyrA_conferring_resistance_to_fluoroquinolones_1", "Escherichia_coli_gyrA_conferring_resistance_to_fluoroquinolones_3", "Escherichia_coli_gyrA_conferring_resistance_to_fluoroquinolones_4", "Escherichia_coli_gyrA_conferring_resistance_to_fluoroquinolones_5", "Escherichia_coli_parC_conferring_resistance_to_fluoroquinolone_1", "Escherichia_coli_parC_conferring_resistance_to_fluoroquinolone_2", "Escherichia_coli_parC_conferring_resistance_to_fluoroquinolone_4", "Escherichia_coli_parC_conferring_resistance_to_fluoroquinolone_5", "APH(3')-Ia_1", "BRP(MBL)_2", "dfrA14_1", "SHV-66_1", "SHV-66_cp1_1", "SHV-66_cp4_1", "SHV-66_cp5_1,SHV-66_cp6_1", "SHV-66_cp7_1", "Staphylococcus_aureus_GlpT_with_mutation_conferring_resistance_to_fosfomycin_1", "Staphylococcus_aureus_GlpT_with_mutation_conferring_resistance_to_fosfomycin_3", "Staphylococcus_aureus_GlpT_with_mutation_conferring_resistance_to_fosfomycin_4", "Staphylococcus_aureus_GlpT_with_mutation_conferring_resistance_to_fosfomycin_5", "QnrS1_1", "aadA5_1", "Escherichia_coli_CyaA_with_mutation_conferring_resistance_to_fosfomycin_1", "CMY-2_1", "CMY-42_1", "dfrA8_1", "AAC(6')-Ib7_1", "catI_1", "QnrB20_1", "AAC(3)-IIe_1", "aadA16_2", "dfrA27_1", "MexA_1", "basR_1", "basR_2", "basR_3", "basR_4", "Klebsiella_pneumoniae_KpnH_2", "Klebsiella_pneumoniae_KpnH_3", "Klebsiella_pneumoniae_KpnH_4", "Klebsiella_pneumoniae_KpnH_6", "Klebsiella_pneumoniae_KpnG_1", "Klebsiella_pneumoniae_KpnG_2", "Klebsiella_pneumoniae_KpnG_3", "Klebsiella_pneumoniae_KpnG_4", "Klebsiella_pneumoniae_KpnG_5", "Klebsiella_pneumoniae_KpnG_6", "adeF_10", "adeF_3", "adeF_4", "adeF_5", "adeF_6", "adeF_cp1_11", "adeF_cp1_12", "adeF_cp1_13", "adeF_cp1_7", "adeF_cp1_9", "oqxA_2", "oqxA_3", "oqxA_4", "adeB_3", "adeB_4", "adeB_5", "Acinetobacter_baumannii_AmvA_4", "Acinetobacter_baumannii_AmvA_5", "Acinetobacter_baumannii_AmvA_6", "norB_1", "norB_2", "ramA_4", "ramA_5", "Klebsiella_pneumoniae_acrA_1", "FosA6_1", "CTX-M-55_1", "ErmB_1", "catII_from_Escherichia_coli_K-12_1", "OXA-232_1", "OXA-9_2", "CTX-M-14_cp1_1", "tet(D)_1", "tet(D)_2", "tet(D)_3", "tet(D)_4", "tet(D)_cp1_2", "SHV-1_1", "ANT(3'')-IIa_1", "ANT(3'')-IIa_2", "ANT(3'')-IIa_3", "KPC-3_cp1_1", "QnrB10_1", "QnrB10_cp1_1", "Morganella_morganii_gyrB_conferring_resistance_to_fluoroquinolone_1", "Morganella_morganii_gyrB_conferring_resistance_to_fluoroquinolone_2", "MexI_1", "MexI_2", "MexI_3", "MexH_1", "MexH_2", "MexH_3", "MexH_4", "MexG_1", "MexG_2", "MexG_4", "golS_2", "golS_3", "golS_6", "golS_cp1_4", "MexB_1", "MexB_2", "MexB_3", "catA4_1", "catA4_5", "ANT(2'')-Ia_1", "SHV-7_1", "Klebsiella_pneumoniae_ramR_mutants_1", "DHA-17_cp1_1", "CMY-108_1", "FosA5_2", "FosA5_3", "FosA5_4", "Klebsiella_pneumoniae_acrR_with_mutation_conferring_multidrug_antibiotic_resistance_1", "Klebsiella_pneumoniae_acrR_with_mutation_conferring_multidrug_antibiotic_resistance_2", "Klebsiella_pneumoniae_acrR_with_mutation_conferring_multidrug_antibiotic_resistance_3", "Klebsiella_pneumoniae_acrR_with_mutation_conferring_multidrug_antibiotic_resistance_4", "Klebsiella_pneumoniae_acrR_with_mutation_conferring_multidrug_antibiotic_resistance_6", "vatF_2", "VIM-27_1", "sul3_1", "FosA2_3", "Agrobacterium_fabrum_chloramphenicol_acetyltransferase_2", "Agrobacterium_fabrum_chloramphenicol_acetyltransferase_3", "Agrobacterium_fabrum_chloramphenicol_acetyltransferase_4", "Agrobacterium_fabrum_chloramphenicol_acetyltransferase_5", "KPC-2_cp1_1", "Salmonella_serovars_gyrB_conferring_resistance_to_fluoroquinolone_1", "ACT-25_1", "cmlA1_1", "OXA-181_1", "SHV-33_1", "NmcR_1", "NmcR_4", "NmcR_5", "SRT-2_1", "SRT-2_2", "SRT-2_3", "SRT-2_5", "Bartonella_bacilliformis_gyrB_conferring_resistance_to_aminocoumarin_1", "AAC(2')-Ia_1", "AAC(2')-Ia_4", "AAC(2')-Ia_5", "aadA14_1", "catIII_2", "catIII_3", "catIII_4", "abeS_1", "abeS_2", "tet(59)_1", "tet(59)_3", "APH(4)-Ia_1", "SHV-30_1", "VIM-1_1", "TEM-135_1", "rmtF_1", "rmtF_2", "NDM-7_1", "DHA-13_1", "Staphylococcus_aureus_murA_with_mutation_conferring_resistance_to_fosfomycin_1", "mdsA_2", "dfrA30_1", "ACT-15_1", "IMP-4_1", "MIR-10_1", "Escherichia_coli_nfsA_mutations_conferring_resistance_to_nitrofurantoin_2", "FosA3_1", "tetM_1", "MCR-1.1_1", "dfrA20_1", "dfrA20_cp1_1", "Tet(X4)_1", "MCR-3.5_1", "FosA4_1", "MCR-3.2_1", "MCR-3.4_1", "vatD_1", "MCR-3.1_1", "MCR-2.1_1", "SHV-60_1", "aadA8b_1", "mef(B)_1"],'CEFEPIME':["Klebsiella_pneumoniae_KpnE_13", "TEM-1_1", "Escherichia_coli_acrR_with_mutation_conferring_multidrug_antibiotic_resistance_14", "sul1_2", "Escherichia_coli_EF-Tu_mutants_conferring_resistance_to_Pulvomycin_1", "Mycobacterium_tuberculosis_thyA_with_mutation_conferring_resistance_to_para-aminosalicylic_acid_3", "CTX-M-15_1", "OXA-1_1"],'IMIPENEM':["Klebsiella_pneumoniae_KpnE_6", "Klebsiella_pneumoniae_KpnF_6", "mdtE_1", "TEM-1_1", "ugd_21", "Escherichia_coli_acrR_with_mutation_conferring_multidrug_antibiotic_resistance_6", "AcrE_1", "Escherichia_coli_ampC_beta-lactamase_1", "qacH_1", "sul1_2", "sul1_cp1_2", "tet(A)_cp1_2", "tetR_cp1_5", "aadA_1", "APH(3'')-Ib_1", "NDM-1_1", "dfrA17_cp1_1", "OXA-9_1", "KPC-3_1", "MexI_1", "KPC-2_1"],'TRIMETROPIM-SULFAMETOXAZOL':["Staphylococcus_aureus_fusA_with_mutation_conferring_resistance_to_fusidic_acid_3", "TEM-1_1", "TEM-1_cp1_1", "Klebsiella_pneumoniae_OmpK37_cp3_89", "Klebsiella_pneumoniae_OmpK37_cp4_89", "Escherichia_coli_ampC_beta-lactamase_1", "Mycobacterium_tuberculosis_rpoB_mutants_conferring_resistance_to_rifampicin_1", "dfrA12_1", "aadA2_2", "qacH_1", "sul1_2", "mphA_1", "Escherichia_coli_EF-Tu_mutants_conferring_resistance_to_Pulvomycin_1", "Escherichia_coli_UhpT_with_mutation_conferring_resistance_to_fosfomycin_3", "dfrA1_1", "aadA_1", "floR_1", "APH(6)-Id_1", "APH(6)-Id_cp1_1", "Escherichia_coli_gyrA_conferring_resistance_to_fluoroquinolones_4", "Escherichia_coli_parC_conferring_resistance_to_fluoroquinolone_4", "NDM-1_1", "catB3_2", "dfrA14_1", "dfrA17_1", "CMY-2_1", "tet(B)_1", "catI_1", "AAC(3)-IIe_1"]}
  #print(modelo.keys())
  #print(modelo[mod])
  relevantes=[]
  lim_id=80
  lim_cov=80
  with open(csv, 'r') as alin:
    for al in alin:
      linea=(al.strip()).split(',')
      #print(linea[2],linea[10])
      if float(linea[2])>lim_id and float(linea[10])>lim_cov:
        relevantes.append(linea[1])#Nombre del query
  print('los relevantes son :')
  print(relevantes)
  presencia=(np.zeros(len(modelo[mod])))
  if modelo.get(mod,None)!=None:
    if len(relevantes)>0:
      cont=0
      mcp='-'
      for m in modelo[mod][:]:
        recp=''
        for re in relevantes:
          if ('cp' in m):
            ub=m.find('_cp')
            #print(ub)
            ub_=m[ub+1:].find('_')
            #print(ub_)
            mcp=str(m[:ub])+str(m[ub+1+ub_:])
            #print(mcp,m)
          if ('cp' in re):
            ubr=re.find('_cp')
            #print(ubr)
            ub_r=m[ubr+1:].find('_')
            #print(ub_r)
            recp=str(re[:ubr])+str(re[ubr+1+ub_r:])
            #print(recp,m)
            #sys.exit() 
          if re==m or recp==m or re==mcp or recp==mcp:
            print('se encontro una congruencia ,,,',re,m)
            presencia[cont]=1
        cont+=1          
    else:
      print('No se encontraron genes de prediccion en los resultados del RGI')
  print(presencia)
  linea=''
  with open(salida, 'w') as sal:
    linea+=tag
    for gen in presencia:
      linea+='\t'+str(gen)
    linea+='\n'
    sal.write(linea)
  return presencia
def blast(fold,ofold,modelo):
  #salida
  print('Ejecutando software Blastp contra base de datos de entrenamiento ...')
  BD='/vault2/homehpc/dtalero/SCRIPT/alelos/PR/BD/'+str(modelo)+'/BlastDB_cp/'+str(modelo)
  print('ejecutando alineamiento contra base de datos :'+BD) 
  os.system('source /vault2/soft/miniconda2/bin/activate sgig;blastp -query '+fold+' -db '+BD+' -task blastp -max_hsps 10 -max_target_seqs 50 -num_threads 20 -outfmt "10 qseqid sseqid pident length qstart qend sstart send qlen slen qcovhsp" -out ' + ofold)
def rgi(fold,ofold):
  print('Ejecutando software RGI ...')
  os.system('source /vault2/soft/miniconda2/bin/activate rgi;rgi main -i '+fold+' -t contig -a DIAMOND -n 8 --include_loose --exclude_nudge -d wgs -o '+ofold)
  print(fold)
def extraer_dat(tabla_rgi,fo,tag):
  limbs,limid,limcov=50,50,80
  general=[]
  if os.path.isfile(tabla_rgi):
    fasta_out=fo
    nombre=tag
    print('nombre')
    #sys.exit()
    arch = open(tabla_rgi,'r')
    line = arch.readlines()
    print(len(line))
    arch.close()
    c=1    
    #with open(fasta_out, 'w') as fasta:
    for l in line:
      if c != 1:
        columns = l.split('\t')
        if str(columns[5])=='Loose':
          if float(columns[7]) > limbs and float(columns[9]) > limid and float(columns[20]) > limcov:
            linea ='>'+nombre+'_'+columns[8]+'\n'+columns[18]+'\n'
            linea_cp=[nombre,columns[8],columns[18]]
            #fasta.write(str(linea))
            general.append(linea_cp)
        else:
          linea ='>'+nombre+'_'+columns[8]+'\n'+columns[18]+'\n'
          linea_cp=[nombre,columns[8],columns[18]]
          #fasta.write(str(linea))
          general.append(linea_cp)
      c+=1
    print('lineas= '+str(c))#al descomentarlo muestra la cantidad de lineas extraidas de cada archivo menos 1 del encabezado.
  else:
    print('Error con software RGI, no se generan archivos necesarios')
  print('busqueda de repeticiones para hacer genes con copias')
  #--------------------------------------------------------------# lista de genes sin repeticion
  #print('esto es general ', general)
  genes=[]
  cuenta,c=0,1
  for gene in general:
      c=1
      for ge in genes:
          if gene[1]==ge:
              c=0
      if c == 1:
          cuenta+=1
          genes.append(gene[1])
  print(str(cuenta)+' genes entre todas las sepas.')#cantidad de genes sin repeticion
  #print(len(genes))
  #---------------------------------------------creacion de lista[aislamiento,gen,repeticion]
  aisgenesrep=[]
  cuenta,cuentagen=0,0
  
  for g in genes:
      cuenta+=1
      cuentagen=0
      for gene in general:
          if nombre == gene[0] and g == gene[1]:
              cuentagen+=1
      aisgenesrep.append([nombre,g,cuentagen])
  #print('esto es genes con rep ',aisgenesrep)
  #-----------------------------------------Escribir el archivo de busqueda implementando las copias
  print('var sin rep ',len(general))
  print('var con rep ',len(aisgenesrep))
  
  with open(fasta_out, 'w') as fasta:
    escribir=''
    escritos=0
    for gens in aisgenesrep:
      copias=1
      for g in general:
        if gens[2]>1:
          if gens[1]==g[1] and gens[2]>=copias:
            #escribir+='>'+nombre+'_'+gens[1]+'_cp'+str(copias)+'\n'+g[2]+'\n'
            escribir+='>'+str(gens[1]).replace(' ','_')+'_cp'+str(copias)+'\n'+g[2]+'\n'
            escritos+=1
            copias+=1
        else:
          if gens[1]==g[1]:
            escribir+='>'+str(gens[1]).replace(' ','_')+'\n'+g[2]+'\n'
            escritos+=1
    fasta.write(escribir)
  print('lo que escribiria es ', escribir)
  print('lo que escribi jue :',escritos)
def help_message():
  print('Script para la prediccion de resistencia')
  print('Usage:')
  print('\tpython3 antibiotico_pred.py [options]')
  print('options:')
  print('\t-f direccion absoluta del archivo fasta del genoma para su prediccion \n \t  \n \t -o direccion de la carpeta de salida')
  print('\t Codigo de los modelos relacionados a un antibiotico  ')
  print('\t <1> TIGECICLINA')
  print('\t <2> TETRACICLINA')
  print('\t <3> AMIKACINA')
  print('\t <4> CIPROFLOXACINA')
  print('\t <5> COLISTINA')
  print('\t <6> CEFEPIME_no')
  print('\t <7> IMIPENEM')
  print('\t <8> TRIMETROPIM-SULFAMETOXAZOL')
  print('\t -m # del modelo a implementar ')  
def main(a):
    antibi={1:'TIGECICLINA',2:'TETRACICLINA',3:'AMIKACINA',4:'CIPROFLOXACINA',5:'COLISTINA',6:'CEFEPIME',7:'IMIPENEM',8:'TRIMETROPIM-SULFAMETOXAZOL'}
    f,o,m=0,0,0
    tag=''
    evaluar=[]
    try:
        opts, args = getopt.getopt(a, 'f:o:m:h')#'hp:k:f:o:r'  
        print(opts)
    except getopt.GetoptError:
        #help_message()
        print('entra en el primer help que de error')
        sys.exit(2)
    for opt, arg in opts:
        if opt == '-h':
            help_message()
            sys.exit()            
        elif opt == '-f': #, '--folderFastqs')
            fold= os.path.realpath(arg)#elimina los caracteres simbolicos de la direccion
            tag=(((fold.split('/'))[-1]).split('.'))[0]
            print('carpeta de entrada :  '+str(fold))
            print('nombre genoma '+tag)
            f=1
        elif opt == '-o': 
            ofold= os.path.realpath(arg)
            print('carpeta de salida :  '+str(ofold))
            o=1
        elif opt == '-m':
          if int(arg)>=1 and int(arg)<=10: 
            modelo=antibi[int(arg)]
            print('Modelo seleccionado :  '+str(modelo))
            m=1
          else:
            print('ha seleccionado un modelo no valido ..')
            sys.exit()
        else:
            help_message()
            sys.exit()
    if f==1 and o==1:
      if os.path.isfile(ofold+'/'+tag+'.rgi.txt')==False:
        rgi(fold,ofold+'/'+tag+'.rgi')#funcion para correr rgi
      else:
        print('Ejecucion RGI previa encontrada\n')
      if os.path.isfile(ofold+'/'+tag+'.prot_out.fasta')==False:
        extraer_dat(ofold+'/'+tag+'.rgi.txt',ofold+'/'+tag+'.prot_out.fasta',tag)#genera un archivo multifasta con los genes loose que superan los limites y los perfect
      else:
        print('Ejecucion de extraccion de datos RGI previa encontrada\n')
      if os.path.isfile(ofold+'/'+tag+'.blast_out.csv')==False:
        blast(ofold+'/'+tag+'.prot_out.fasta',ofold+'/'+tag+'.blast_out.csv',modelo)#realiza un alineamiento contra la base de datos de entrenamiento
      else:
        print('Ejecucion BLASTP previa encontrada\n')
      if os.path.isfile(ofold+'/'+tag+'.matriz_in.csv')==False:#si ya esta la matriz leerla
        evaluar=matriz(ofold+'/'+tag+'.blast_out.csv',ofold+'/'+tag+'.matriz_in.csv',modelo,tag)
      """if os.path.isfile(ofold+'/'+tag+'.predicc.txt')==False:
        predecir(ofold+'/'+tag+'.matriz_in.csv',modelo,evaluar)"""
    else:
      if f==0:
        help_message()
        print('no se ha seleccionado un archivo de entrada (-f ../../g.fasta)')
        sys.exit()
      if o==0:
        help_message()
        print('no se ha seleccionado un directorio de salida (-o ../../example/)')
        sys.exit()
if __name__ == '__main__':#__name__ variable de entorno que es main si se ejecuta el mismo programa desde el compilador, pero si el script es llamado de otro script tendra el valor de 'codigo.py'
  if len(sys.argv[1:]) != 0:
     main(sys.argv[1:])
  else:
      help_message()