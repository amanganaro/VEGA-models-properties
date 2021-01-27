package utils;

import insilico.core.exception.GenericFailureException;
import insilico.core.exception.InitFailureException;
import insilico.core.model.InsilicoModel;
//import insilico.model.daphnia.demetra.ismDaphniaDemetra;

import java.util.ArrayList;
import java.util.Arrays;

public class ModelsList {

//    private static final String MODEL_MUTA_CAESAR = "muta_caesar";
//    private static final String MODEL_MUTA_SARPY = "muta_sarpy";
//    private static final String MODEL_MUTA_ISS = "muta_iss";
//    private static final String MODEL_MUTA_KNN = "muta_knn";
//    private static final String MODEL_CARC_CAESAR = "carc_caesar";
//    private static final String MODEL_CARC_ISS = "carc_iss";
//    private static final String MODEL_CARC_ANTARES = "carc_antares";
//    private static final String MODEL_CARC_ISSCAN = "carc_isscan";
//    private static final String MODEL_CARC_SFOCLASS = "carc_sfoclass";
//    private static final String MODEL_CARC_SFOREGR = "carc_sforegr";
//    private static final String MODEL_CARC_SFICLASS = "carc_sficlass";
//    private static final String MODEL_CARC_SFIREGR = "carc_sfiregr";
//    private static final String MODEL_DEVTOX_CAESAR = "devtox_caesar";
//    private static final String MODEL_DEVTOX_PG = "devtox_pg";
//    private static final String MODEL_DEVTOX_ZEBRA = "devtox_zebrafish";
//    private static final String MODEL_CHORMOSOMAL_CORAL = "chrom_coral";
//    private static final String MODEL_RBA_IRFMN = "rba_irfmn";
//    private static final String MODEL_RBA_CERAPP = "rba_cerapp";
//    private static final String MODEL_RBA_COMPARA = "rba_compara";
//    private static final String MODEL_THYRALPHA_NRMEA = "tralpha_nrmea";
//    private static final String MODEL_THYRBETA_NRMEA = "trbeta_nrmea";
//    private static final String MODEL_AROMATASE_IRFMN = "aroma_irfmn";
//    private static final String MODEL_SKIN_CAESAR = "skin_caesar";
//    private static final String MODEL_SKIN_IRFMN = "skin_irfmn";
//    private static final String MODEL_HEPA_IRFMN = "hepa_irfmn";
//    private static final String MODEL_TISSUE_BLOOD_INERIS = "tissueblood_ineris";
//    private static final String MODEL_TOTAL_BODY_HL_QSARINS = "totalhl_qsarins";
//    private static final String MODEL_MICRONUCLEUS_VITRO = "mn_invitro";
//    private static final String MODEL_MICRONUCLEUS_VIVO = "mn_invivo";
//    private static final String MODEL_NOAEL_CORAL = "noael_coral";
//    private static final String MODEL_CRAMER_TOXTREE = "cramer_toxtree";
//    private static final String MODEL_MOA_TEST = "moa_test";
//    private static final String MODEL_VERHAAR_TOXTREE = "verhaar_toxtree";
//    private static final String MODEL_FISH_IRFMN = "fish_irfmn";
//    private static final String MODEL_FISH_KNN = "fish_knn";
//    private static final String MODEL_FISH_NIC = "fish_nic";
//    private static final String MODEL_FISH_LC50 = "fish_lc50";
//    private static final String MODEL_FISH_COMB = "fish_comb";
//    private static final String MODEL_FISH_NOEC = "fish_noec";
//    private static final String MODEL_FATHEAD_EPA = "fathead_epa";
//    private static final String MODEL_FATHEAD_KNN = "fathead_knn";
//    private static final String MODEL_GUPPY_KNN = "guppy_knn";
//    private static final String MODEL_DAPHNIA_EPA = "daphnia_epa";
//    private static final String MODEL_DAPHNIA_DEMETRA = "daphnia_demetra";
//    private static final String MODEL_DAPHNIA_EC50 = "daphnia_ec50";
//    private static final String MODEL_DAPHNIA_COMB = "daphnia_comb";
//    private static final String MODEL_DAPHNIA_NOEC = "daphnia_noec";
//    private static final String MODEL_ALGAE_EC50 = "algae_ec50";
//    private static final String MODEL_ALGAE_COMB = "algae_comb";
//    private static final String MODEL_ALGAE_NOEC = "algae_noec";
//    private static final String MODEL_ALGAE_COMBCLASS = "algae_combclass";
//    private static final String MODEL_BEE_KNN = "bee_knn";
//    private static final String MODEL_SLUDGE_COMB = "sludge_comb";
//    private static final String MODEL_SLUDGE_COMBCLASS = "sludge_combclass";
//    private static final String MODEL_BCF_CAESAR = "bcf_caesar";
//    private static final String MODEL_BCF_MEYLAN = "bcf_meyla";
//    private static final String MODEL_BCF_KNN = "bcf_knn";
//    private static final String MODEL_BCF_ARNOT = "bcf_arnot";
//    private static final String MODEL_KM_ARNOT = "km_arnot";
//    private static final String MODEL_RB_IRFMN = "rb_irfmn";
//    private static final String MODEL_P_SEDIMENT_IRFMN = "p_sed_irfmn";
//    private static final String MODEL_PQ_SEDIMENT_IRFMN = "pq_sed_irfmn";
//    private static final String MODEL_P_SOIL_IRFMN = "p_soil_irfmn";
//    private static final String MODEL_PQ_SOIL_IRFMN = "pq_soil_irfmn";
//    private static final String MODEL_P_WATER_IRFMN = "p_water_irfmn";
//    private static final String MODEL_PQ_WATER_IRFMN = "pq_water_irfmn";
//    private static final String MODEL_P_AIR_CORAL = "p_air_coral";
//    private static final String MODEL_KOC_OPERA = "koc_opera";
//    private static final String MODEL_KOA_OPERA = "koa_opera";
//    private static final String MODEL_LOGP_MEYLAN = "logp_meylan";
//    private static final String MODEL_LOGP_MLOGP = "logp_mlogp";
//    private static final String MODEL_LOGP_ALOGP = "logp_alogp";
//    private static final String MODEL_SKIN_PERM_POTTS = "skin_perm_potts";
//    private static final String MODEL_SKIN_PERM_TENBERGE = "skin_perm_ten";
//    private static final String MODEL_WS_IRFMN = "water_sol";
//    private static final String MODEL_HENRY_OPERA = "henry_opera";
//    private static final String MODEL_HYDROLYSIS_CORAL = "hydro_coral";
//
//    private static final String CONSENSUS_MUTA = "muta_cons";
//
//    public static final String[] MODEL_TAGS = {
//            MODEL_MUTA_CAESAR,
//            MODEL_MUTA_SARPY,
//            MODEL_MUTA_ISS,
//            MODEL_MUTA_KNN,
//            MODEL_CARC_CAESAR,
//            MODEL_CARC_ISS,
//            MODEL_CARC_ANTARES,
//            MODEL_CARC_ISSCAN,
//            MODEL_CARC_SFOCLASS,
//            MODEL_CARC_SFOREGR,
//            MODEL_CARC_SFICLASS,
//            MODEL_CARC_SFIREGR,
//            MODEL_DEVTOX_CAESAR,
//            MODEL_DEVTOX_PG,
//            MODEL_DEVTOX_ZEBRA,
//            MODEL_CHORMOSOMAL_CORAL,
//            MODEL_RBA_IRFMN,
//            MODEL_RBA_CERAPP,
//            MODEL_RBA_COMPARA,
//            MODEL_THYRALPHA_NRMEA,
//            MODEL_THYRBETA_NRMEA,
//            MODEL_AROMATASE_IRFMN,
//            MODEL_SKIN_CAESAR,
//            MODEL_SKIN_IRFMN,
//            MODEL_HEPA_IRFMN,
//            MODEL_TISSUE_BLOOD_INERIS,
//            MODEL_TOTAL_BODY_HL_QSARINS,
//            MODEL_MICRONUCLEUS_VITRO,
//            MODEL_MICRONUCLEUS_VIVO,
//            MODEL_NOAEL_CORAL,
//            MODEL_CRAMER_TOXTREE,
//            MODEL_MOA_TEST,
//            MODEL_VERHAAR_TOXTREE,
//            MODEL_FISH_IRFMN,
//            MODEL_FISH_KNN,
//            MODEL_FISH_NIC,
//            MODEL_FISH_LC50,
//            MODEL_FISH_COMB,
//            MODEL_FISH_NOEC,
//            MODEL_FATHEAD_EPA,
//            MODEL_FATHEAD_KNN,
//            MODEL_GUPPY_KNN,
//            MODEL_DAPHNIA_EPA,
//            MODEL_DAPHNIA_DEMETRA,
//            MODEL_DAPHNIA_EC50,
//            MODEL_DAPHNIA_COMB,
//            MODEL_DAPHNIA_NOEC,
//            MODEL_ALGAE_EC50,
//            MODEL_ALGAE_COMB,
//            MODEL_ALGAE_NOEC,
//            MODEL_ALGAE_COMBCLASS,
//            MODEL_BEE_KNN,
//            MODEL_SLUDGE_COMB,
//            MODEL_SLUDGE_COMBCLASS,
//            MODEL_BCF_CAESAR,
//            MODEL_BCF_MEYLAN,
//            MODEL_BCF_KNN,
//            MODEL_BCF_ARNOT,
//            MODEL_KM_ARNOT,
//            MODEL_RB_IRFMN,
//            MODEL_P_SEDIMENT_IRFMN,
//            MODEL_PQ_SEDIMENT_IRFMN,
//            MODEL_P_SOIL_IRFMN,
//            MODEL_PQ_SOIL_IRFMN,
//            MODEL_P_WATER_IRFMN,
//            MODEL_PQ_WATER_IRFMN,
//            MODEL_P_AIR_CORAL,
//            MODEL_LOGP_MEYLAN,
//            MODEL_LOGP_MLOGP,
//            MODEL_LOGP_ALOGP,
//            MODEL_WS_IRFMN,
//            MODEL_HYDROLYSIS_CORAL,
//            MODEL_HENRY_OPERA,
//            MODEL_KOC_OPERA,
//            MODEL_KOA_OPERA,
//            MODEL_SKIN_PERM_POTTS,
//            MODEL_SKIN_PERM_TENBERGE
//    };
//
//    public static final String[] CONSENSUS_TAGS = {
//            CONSENSUS_MUTA
//    };
//
//    public static String[] getModelTags() {
//        return MODEL_TAGS;
//    }
//
//    public static boolean ModelTagExists(String tag){
//        return Arrays.asList(MODEL_TAGS).contains(tag);
//    }
//
//    public static ArrayList<InsilicoModel> getAllModels() throws GenericFailureException, InitFailureException {
//        ArrayList<InsilicoModel> modelsList = new ArrayList<>();
//        for(String tag: getModelTags()){
//            modelsList.add(getModel(tag));
//        }
//        return modelsList;
//    }
//
//    public static InsilicoModel getModel(String tag) throws GenericFailureException, InitFailureException {
//
//        //todo Other models
//        if (tag.equalsIgnoreCase(MODEL_DAPHNIA_DEMETRA)) return new ismDaphniaDemetra();
//
//
//        throw  new GenericFailureException("Tag not matching with any available model");
//    }
//
//    public static void PrintAllModelsList() throws Exception {
//        ArrayList<InsilicoModel> modelsList = new ArrayList<>();
//        int index = 0;
//        for(InsilicoModel model : modelsList){
//            System.out.println(++index + "\t" + model.getInfo().getKey() + "\t" + model.getInfo().getName()
//            + "\t" + model.getInfo().getSummary());
//        }
//    }

    // todo - finish method from InsilicoModel_OLD



}
