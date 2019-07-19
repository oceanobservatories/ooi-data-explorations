function [uframe_dataset_name, data_url] = RecoveredHost_M2M_URLs(mooring_name,node,instrument_class)
%.. Explicitly construct UFrame dataset names
if strcmp(mooring_name,'CE01ISSM') && strcmp(node,'NSIF') && strcmp(instrument_class,'CTDBP')
    uframe_dataset_name = 'CE01ISSM/RID16/03-CTDBPC000/recovered_host/ctdbp_cdef_dcl_instrument_recovered';
elseif strcmp(mooring_name,'CE01ISSM') && strcmp(node,'MFN') && strcmp(instrument_class,'CTDBP')
    uframe_dataset_name = 'CE01ISSM/MFD37/03-CTDBPC000/recovered_host/ctdbp_cdef_dcl_instrument_recovered';
elseif strcmp(mooring_name,'CE01ISSM') && strcmp(node,'BUOY') && strcmp(instrument_class,'CTDBP')
    uframe_dataset_name = 'CE01ISSM/SBD17/06-CTDBPC000/recovered_host/ctdbp_cdef_dcl_instrument_recovered';
elseif strcmp(mooring_name,'CE06ISSM') && strcmp(node,'NSIF') && strcmp(instrument_class,'CTDBP')
    uframe_dataset_name = 'CE06ISSM/RID16/03-CTDBPC000/recovered_host/ctdbp_cdef_dcl_instrument_recovered';
elseif strcmp(mooring_name,'CE06ISSM') && strcmp(node,'MFN') && strcmp(instrument_class,'CTDBP')
    uframe_dataset_name = 'CE06ISSM/MFD37/03-CTDBPC000/recovered_host/ctdbp_cdef_dcl_instrument_recovered';
elseif strcmp(mooring_name,'CE06ISSM') && strcmp(node,'BUOY') && strcmp(instrument_class,'CTDBP')
    uframe_dataset_name = 'CE06ISSM/SBD17/06-CTDBPC000/recovered_host/ctdbp_cdef_dcl_instrument_recovered';
elseif strcmp(mooring_name,'CE02SHSM') && strcmp(node,'NSIF') && strcmp(instrument_class,'CTDBP')
    uframe_dataset_name = 'CE02SHSM/RID27/03-CTDBPC000/recovered_host/ctdbp_cdef_dcl_instrument_recovered';
elseif strcmp(mooring_name,'CE07SHSM') && strcmp(node,'NSIF') && strcmp(instrument_class,'CTDBP')
    uframe_dataset_name = 'CE07SHSM/RID27/03-CTDBPC000/recovered_host/ctdbp_cdef_dcl_instrument_recovered';
elseif strcmp(mooring_name,'CE04OSSM') && strcmp(node,'NSIF') && strcmp(instrument_class,'CTDBP')
    uframe_dataset_name = 'CE04OSSM/RID27/03-CTDBPC000/recovered_host/ctdbp_cdef_dcl_instrument_recovered';
elseif strcmp(mooring_name,'CE09OSSM') && strcmp(node,'NSIF') && strcmp(instrument_class,'CTDBP')
    uframe_dataset_name = 'CE09OSSM/RID27/03-CTDBPC000/recovered_host/ctdbp_cdef_dcl_instrument_recovered';
elseif strcmp(mooring_name,'CE07SHSM') && strcmp(node,'MFN') && strcmp(instrument_class,'CTDBP')
    uframe_dataset_name = 'CE07SHSM/MFD37/03-CTDBPC000/recovered_host/ctdbp_cdef_dcl_instrument_recovered';
elseif strcmp(mooring_name,'CE09OSSM') && strcmp(node,'MFN') && strcmp(instrument_class,'CTDBP')
    uframe_dataset_name = 'CE09OSSM/MFD37/03-CTDBPE000/recovered_host/ctdbp_cdef_dcl_instrument_recovered';
elseif strcmp(mooring_name,'CE02SHSM') && strcmp(node,'BUOY') && strcmp(instrument_class,'METBK')
    uframe_dataset_name = 'CE02SHSM/SBD11/06-METBKA000/recovered_host/metbk_a_dcl_instrument_recovered';
elseif strcmp(mooring_name,'CE07SHSM') && strcmp(node,'BUOY') && strcmp(instrument_class,'METBK')
    uframe_dataset_name = 'CE07SHSM/SBD11/06-METBKA000/recovered_host/metbk_a_dcl_instrument_recovered';
elseif strcmp(mooring_name,'CE04OSSM') && strcmp(node,'BUOY') && strcmp(instrument_class,'METBK')
    uframe_dataset_name = 'CE04OSSM/SBD11/06-METBKA000/recovered_host/metbk_a_dcl_instrument_recovered';
elseif strcmp(mooring_name,'CE09OSSM') && strcmp(node,'BUOY') && strcmp(instrument_class,'METBK')
    uframe_dataset_name = 'CE09OSSM/SBD11/06-METBKA000/recovered_host/metbk_a_dcl_instrument_recovered';
elseif strcmp(mooring_name,'CE02SHSM') && strcmp(node,'BUOY') && strcmp(instrument_class,'PCO2A')
    uframe_dataset_name = 'CE02SHSM/SBD12/04-PCO2AA000/recovered_host/pco2a_a_dcl_instrument_air_recovered';
elseif strcmp(mooring_name,'CE04OSSM') && strcmp(node,'BUOY') && strcmp(instrument_class,'PCO2A')
    uframe_dataset_name = 'CE04OSSM/SBD12/04-PCO2AA000/recovered_host/pco2a_a_dcl_instrument_air_recovered';
elseif strcmp(mooring_name,'CE07SHSM') && strcmp(node,'BUOY') && strcmp(instrument_class,'PCO2A')
    uframe_dataset_name = 'CE07SHSM/SBD12/04-PCO2AA000/recovered_host/pco2a_a_dcl_instrument_air_recovered';
elseif strcmp(mooring_name,'CE09OSSM') && strcmp(node,'BUOY') && strcmp(instrument_class,'PCO2A')
    uframe_dataset_name = 'CE09OSSM/SBD12/04-PCO2AA000/recovered_host/pco2a_a_dcl_instrument_air_recovered';
elseif strcmp(mooring_name,'CE01ISSM') && strcmp(node,'NSIF') && strcmp(instrument_class,'PCO2W')
    uframe_dataset_name = 'CE01ISSM/RID16/05-PCO2WB000/recovered_host/pco2w_abc_dcl_instrument_recovered';
elseif strcmp(mooring_name,'CE01ISSM') && strcmp(node,'MFN') && strcmp(instrument_class,'PCO2W')
    uframe_dataset_name = 'CE01ISSM/MFD35/05-PCO2WB000/recovered_host/pco2w_abc_dcl_instrument_recovered';
elseif strcmp(mooring_name,'CE06ISSM') && strcmp(node,'NSIF') && strcmp(instrument_class,'PCO2W')
    uframe_dataset_name = 'CE06ISSM/RID16/05-PCO2WB000/recovered_host/pco2w_abc_dcl_instrument_recovered';
elseif strcmp(mooring_name,'CE06ISSM') && strcmp(node,'MFN') && strcmp(instrument_class,'PCO2W')
    uframe_dataset_name = 'CE06ISSM/MFD35/05-PCO2WB000/recovered_host/pco2w_abc_dcl_instrument_recovered';
elseif strcmp(mooring_name,'CE07SHSM') && strcmp(node,'MFN') && strcmp(instrument_class,'PCO2W')
    uframe_dataset_name = 'CE07SHSM/MFD35/05-PCO2WB000/recovered_host/pco2w_abc_dcl_instrument_recovered';
elseif strcmp(mooring_name,'CE09OSSM') && strcmp(node,'MFN') && strcmp(instrument_class,'PCO2W')
    uframe_dataset_name = 'CE09OSSM/MFD35/05-PCO2WB000/recovered_host/pco2w_abc_dcl_instrument_recovered';
elseif strcmp(mooring_name,'CE01ISSM') && strcmp(node,'NSIF') && strcmp(instrument_class,'PHSEN')
    uframe_dataset_name = 'CE01ISSM/RID16/06-PHSEND000/recovered_host/phsen_abcdef_dcl_instrument_recovered';
elseif strcmp(mooring_name,'CE02SHSM') && strcmp(node,'NSIF') && strcmp(instrument_class,'PHSEN')
    uframe_dataset_name = 'CE02SHSM/RID26/06-PHSEND000/recovered_host/phsen_abcdef_dcl_instrument_recovered';
elseif strcmp(mooring_name,'CE04OSSM') && strcmp(node,'NSIF') && strcmp(instrument_class,'PHSEN')
    uframe_dataset_name = 'CE04OSSM/RID26/06-PHSEND000/recovered_host/phsen_abcdef_dcl_instrument_recovered';
elseif strcmp(mooring_name,'CE06ISSM') && strcmp(node,'NSIF') && strcmp(instrument_class,'PHSEN')
    uframe_dataset_name = 'CE06ISSM/RID16/06-PHSEND000/recovered_host/phsen_abcdef_dcl_instrument_recovered';
elseif strcmp(mooring_name,'CE07SHSM') && strcmp(node,'NSIF') && strcmp(instrument_class,'PHSEN')
    uframe_dataset_name = 'CE07SHSM/RID26/06-PHSEND000/recovered_host/phsen_abcdef_dcl_instrument_recovered';
elseif strcmp(mooring_name,'CE09OSSM') && strcmp(node,'NSIF') && strcmp(instrument_class,'PHSEN')
    uframe_dataset_name = 'CE09OSSM/RID26/06-PHSEND000/recovered_host/phsen_abcdef_dcl_instrument_recovered';
elseif strcmp(mooring_name,'CE01ISSM') && strcmp(node,'MFN') && strcmp(instrument_class,'PHSEN')
    uframe_dataset_name = 'CE01ISSM/MFD35/06-PHSEND000/recovered_host/phsen_abcdef_dcl_instrument_recovered';
elseif strcmp(mooring_name,'CE06ISSM') && strcmp(node,'MFN') && strcmp(instrument_class,'PHSEN')
    uframe_dataset_name = 'CE06ISSM/MFD35/06-PHSEND000/recovered_host/phsen_abcdef_dcl_instrument_recovered';
elseif strcmp(mooring_name,'CE07SHSM') && strcmp(node,'MFN') && strcmp(instrument_class,'PHSEN')
    uframe_dataset_name = 'CE07SHSM/MFD35/06-PHSEND000/recovered_host/phsen_abcdef_dcl_instrument_recovered';
elseif strcmp(mooring_name,'CE09OSSM') && strcmp(node,'MFN') && strcmp(instrument_class,'PHSEN')
    uframe_dataset_name = 'CE09OSSM/MFD35/06-PHSEND000/recovered_host/phsen_abcdef_dcl_instrument_recovered';
elseif strcmp(mooring_name,'CE01ISSM') && strcmp(node,'NSIF') && strcmp(instrument_class,'DOSTA')
    uframe_dataset_name = 'CE01ISSM/RID16/03-DOSTAD000/recovered_host/dosta_abcdjm_ctdbp_dcl_instrument_recovered';
elseif strcmp(mooring_name,'CE02SHSM') && strcmp(node,'NSIF') && strcmp(instrument_class,'DOSTA')
    uframe_dataset_name = 'CE02SHSM/RID27/04-DOSTAD000/recovered_host/dosta_abcdjm_dcl_instrument_recovered';
elseif strcmp(mooring_name,'CE04OSSM') && strcmp(node,'NSIF') && strcmp(instrument_class,'DOSTA')
    uframe_dataset_name = 'CE04OSSM/RID27/04-DOSTAD000/recovered_host/dosta_abcdjm_dcl_instrument_recovered';
elseif strcmp(mooring_name,'CE06ISSM') && strcmp(node,'NSIF') && strcmp(instrument_class,'DOSTA')
    uframe_dataset_name = 'CE06ISSM/RID16/03-DOSTAD000/recovered_host/dosta_abcdjm_ctdbp_dcl_instrument_recovered';
elseif strcmp(mooring_name,'CE07SHSM') && strcmp(node,'NSIF') && strcmp(instrument_class,'DOSTA')
    uframe_dataset_name = 'CE07SHSM/RID27/04-DOSTAD000/recovered_host/dosta_abcdjm_dcl_instrument_recovered';
elseif strcmp(mooring_name,'CE09OSSM') && strcmp(node,'NSIF') && strcmp(instrument_class,'DOSTA')
    uframe_dataset_name = 'CE09OSSM/RID27/04-DOSTAD000/recovered_host/dosta_abcdjm_dcl_instrument_recovered';
elseif strcmp(mooring_name,'CE01ISSM') && strcmp(node,'MFN') && strcmp(instrument_class,'DOSTA')
    uframe_dataset_name = 'CE01ISSM/MFD37/03-DOSTAD000/recovered_host/dosta_abcdjm_ctdbp_dcl_instrument_recovered';
elseif strcmp(mooring_name,'CE06ISSM') && strcmp(node,'MFN') && strcmp(instrument_class,'DOSTA')
    uframe_dataset_name = 'CE06ISSM/MFD37/03-DOSTAD000/recovered_host/dosta_abcdjm_ctdbp_dcl_instrument_recovered';
elseif strcmp(mooring_name,'CE07SHSM') && strcmp(node,'MFN') && strcmp(instrument_class,'DOSTA')
    uframe_dataset_name = 'CE07SHSM/MFD37/03-DOSTAD000/recovered_host/dosta_abcdjm_ctdbp_dcl_instrument_recovered';
elseif strcmp(mooring_name,'CE09OSSM') && strcmp(node,'MFN') && strcmp(instrument_class,'DOSTA')
    uframe_dataset_name = 'CE09OSSM/MFD37/03-DOSTAD000/recovered_host/dosta_abcdjm_ctdbp_dcl_instrument_recovered';
elseif strcmp(mooring_name,'CE01ISSM') && strcmp(node,'NSIF') && strcmp(instrument_class,'FLORT')
    uframe_dataset_name = 'CE01ISSM/RID16/02-FLORTD000/recovered_host/flort_sample';
elseif strcmp(mooring_name,'CE01ISSM') && strcmp(node,'BUOY') && strcmp(instrument_class,'FLORT')
    uframe_dataset_name = 'CE01ISSM/SBD17/06-FLORTD000/recovered_host/flort_sample';
elseif strcmp(mooring_name,'CE06ISSM') && strcmp(node,'NSIF') && strcmp(instrument_class,'FLORT')
    uframe_dataset_name = 'CE06ISSM/RID16/02-FLORTD000/recovered_host/flort_sample';
elseif strcmp(mooring_name,'CE06ISSM') && strcmp(node,'BUOY') && strcmp(instrument_class,'FLORT')
    uframe_dataset_name = 'CE06ISSM/SBD17/06-FLORTD000/recovered_host/flort_sample';
elseif strcmp(mooring_name,'CE02SHSM') && strcmp(node,'NSIF') && strcmp(instrument_class,'FLORT')
    uframe_dataset_name = 'CE02SHSM/RID27/02-FLORTD000/recovered_host/flort_sample';
elseif strcmp(mooring_name,'CE07SHSM') && strcmp(node,'NSIF') && strcmp(instrument_class,'FLORT')
    uframe_dataset_name = 'CE07SHSM/RID27/02-FLORTD000/recovered_host/flort_sample';
elseif strcmp(mooring_name,'CE04OSSM') && strcmp(node,'NSIF') && strcmp(instrument_class,'FLORT')
    uframe_dataset_name = 'CE04OSSM/RID27/02-FLORTD000/recovered_host/flort_sample';
elseif strcmp(mooring_name,'CE09OSSM') && strcmp(node,'NSIF') && strcmp(instrument_class,'FLORT')
    uframe_dataset_name = 'CE09OSSM/RID27/02-FLORTD000/recovered_host/flort_sample';
elseif strcmp(mooring_name,'CE09OSPM') && strcmp(node,'WFP') && strcmp(instrument_class,'FLORT')
    uframe_dataset_name = 'CE09OSPM/WFP01/04-FLORTK000/recovered_wfp/flort_sample';
elseif strcmp(mooring_name,'CE02SHSM') && strcmp(node,'BUOY') && strcmp(instrument_class,'WAVSS')
    uframe_dataset_name = 'CE02SHSM/SBD12/05-WAVSSA000/recovered_host/wavss_a_dcl_statistics_recovered';
elseif strcmp(mooring_name,'CE04OSSM') && strcmp(node,'BUOY') && strcmp(instrument_class,'WAVSS')
    uframe_dataset_name = 'CE04OSSM/SBD12/05-WAVSSA000/recovered_host/wavss_a_dcl_statistics_recovered';
elseif strcmp(mooring_name,'CE07SHSM') && strcmp(node,'BUOY') && strcmp(instrument_class,'WAVSS')
    uframe_dataset_name = 'CE07SHSM/SBD12/05-WAVSSA000/recovered_host/wavss_a_dcl_statistics_recovered';
elseif strcmp(mooring_name,'CE09OSSM') && strcmp(node,'BUOY') && strcmp(instrument_class,'WAVSS')
    uframe_dataset_name = 'CE09OSSM/SBD12/05-WAVSSA000/recovered_host/wavss_a_dcl_statistics_recovered';
elseif strcmp(mooring_name,'CE01ISSM') && strcmp(node,'MFN') && strcmp(instrument_class,'WAVSS')
    uframe_dataset_name = 'CE01ISSM/MFD35/04-ADCPTM000/recovered_inst/adcpt_m_instrument_log9_recovered';
elseif strcmp(mooring_name,'CE06ISSM') && strcmp(node,'MFN') && strcmp(instrument_class,'WAVSS')
    uframe_dataset_name = 'CE06ISSM/MFD35/04-ADCPTM000/recovered_inst/adcpt_m_instrument_log9_recovered';
elseif strcmp(mooring_name,'CE01ISSM') && strcmp(node,'NSIF') && strcmp(instrument_class,'NUTNR')
    uframe_dataset_name{1} = 'CE01ISSM/RID16/07-NUTNRB000/recovered_host/nutnr_b_dcl_conc_instrument_recovered';
    uframe_dataset_name{2} = 'CE01ISSM/RID16/07-NUTNRB000/recovered_host/suna_dcl_recovered';
elseif strcmp(mooring_name,'CE02SHSM') && strcmp(node,'NSIF') && strcmp(instrument_class,'NUTNR')
    uframe_dataset_name{1} = 'CE02SHSM/RID26/07-NUTNRB000/recovered_host/nutnr_b_dcl_conc_instrument_recovered';
    uframe_dataset_name{2} = 'CE02SHSM/RID26/07-NUTNRB000/recovered_host/suna_dcl_recovered';
elseif strcmp(mooring_name,'CE04OSSM') && strcmp(node,'NSIF') && strcmp(instrument_class,'NUTNR')
    uframe_dataset_name{1} = 'CE04OSSM/RID26/07-NUTNRB000/recovered_host/nutnr_b_dcl_conc_instrument_recovered';
    uframe_dataset_name{2} = 'CE04OSSM/RID26/07-NUTNRB000/recovered_host/suna_dcl_recovered';
elseif strcmp(mooring_name,'CE06ISSM') && strcmp(node,'NSIF') && strcmp(instrument_class,'NUTNR')
    uframe_dataset_name{1} = 'CE06ISSM/RID16/07-NUTNRB000/recovered_host/nutnr_b_dcl_conc_instrument_recovered';
    uframe_dataset_name{2} = 'CE06ISSM/RID16/07-NUTNRB000/recovered_host/suna_dcl_recovered';
elseif strcmp(mooring_name,'CE07SHSM') && strcmp(node,'NSIF') && strcmp(instrument_class,'NUTNR')
    uframe_dataset_name{1} = 'CE07SHSM/RID26/07-NUTNRB000/recovered_host/nutnr_b_dcl_conc_instrument_recovered';
    uframe_dataset_name{2} = 'CE07SHSM/RID26/07-NUTNRB000/recovered_host/nutnr_b_dcl_full_instrument_recovered';
    uframe_dataset_name{3} = 'CE07SHSM/RID26/07-NUTNRB000/recovered_host/suna_dcl_recovered';
elseif strcmp(mooring_name,'CE09OSSM') && strcmp(node,'NSIF') && strcmp(instrument_class,'NUTNR')
    uframe_dataset_name{1} = 'CE09OSSM/RID26/07-NUTNRB000/recovered_host/nutnr_b_dcl_conc_instrument_recovered';
    uframe_dataset_name{2} = 'CE09OSSM/RID26/07-NUTNRB000/recovered_host/suna_dcl_recovered';
else
    error('Illegal mooring_name or node or instrument_class or combination thereof.');
end
m2m_url = "https://ooinet.oceanobservatories.org/api/m2m/12576/sensor/inv/";
data_url = strcat(m2m_url, uframe_dataset_name);
end
