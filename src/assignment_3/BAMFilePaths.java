package assignment_3;

import java.util.ArrayList;

import io.ConfigReader;

public class BAMFilePaths {

	private static String pathToMappinginfo = ConfigReader.readConfig().get("mappinginfo_a3t2");
	private static ArrayList<String> pathList = new ArrayList<>();
	
	public static void getPathList(){
		pathList.add(E2er_rep1_contextmap);
		pathList.add(E2er_rep1_hisat);
		pathList.add(E2er_rep1_star);
		pathList.add(E2er_rep1_tophat);
		pathList.add(E2er_rep2_contextmap);
		pathList.add(E2er_rep2_hisat);
		pathList.add(E2er_rep2_star);
		pathList.add(E2er_rep2_tophat);
		pathList.add(E2er_rep3_contextmap);
		pathList.add(E2er_rep3_hisat);
		pathList.add(E2er_rep3_star);
		pathList.add(E2er_rep3_tophat);
		pathList.add(E2erchx_rep1_contextmap);
		pathList.add(E2erchx_rep1_hisat);
		pathList.add(E2erchx_rep1_star);
		pathList.add(E2erchx_rep1_tophat);
		pathList.add(E2erchx_rep2_contextmap);
		pathList.add(E2erchx_rep2_hisat);
		pathList.add(E2erchx_rep2_star);
		pathList.add(E2erchx_rep2_tophat);
		pathList.add(E2erchx_rep3_contextmap);
		pathList.add(E2erchx_rep3_hisat);
		pathList.add(E2erchx_rep3_star);
		pathList.add(E2erchx_rep3_tophat);
		pathList.add(E2noer_rep1_contextmap);
		pathList.add(E2noer_rep1_hisat);
		pathList.add(E2noer_rep1_star);
		pathList.add(E2noer_rep1_tophat);
		pathList.add(E2noer_rep2_contextmap);
		pathList.add(E2noer_rep2_hisat);
		pathList.add(E2noer_rep2_star);
		pathList.add(E2noer_rep2_tophat);
		pathList.add(E2noer_rep3_contextmap);
		pathList.add(E2noer_rep3_hisat);
		pathList.add(E2noer_rep3_star);
		pathList.add(E2noer_rep3_tophat);
		pathList.add(E2noerchx_rep1_contextmap);
		pathList.add(E2noerchx_rep1_hisat);
		pathList.add(E2noerchx_rep1_star);
		pathList.add(E2noerchx_rep1_tophat);
		pathList.add(E2noerchx_rep2_contextmap);
		pathList.add(E2noerchx_rep2_hisat);
		pathList.add(E2noerchx_rep2_star);
		pathList.add(E2noerchx_rep2_tophat);
		pathList.add(E2noerchx_rep3_contextmap);
		pathList.add(E2noerchx_rep3_hisat);
		pathList.add(E2noerchx_rep3_star);
		pathList.add(E2noerchx_rep3_tophat);
	}
	
	/*
	 * E2_er
	 */
	public static final String E2er_rep1_contextmap = pathToMappinginfo.substring(0, pathToMappinginfo.lastIndexOf("/"))+"/EBNA2/E2_er/rep1/contextmap.bam";
	public static final String E2er_rep1_hisat = pathToMappinginfo.substring(0, pathToMappinginfo.lastIndexOf("/"))+"/EBNA2/E2_er/rep1/hisat.bam";
	public static final String E2er_rep1_star = pathToMappinginfo.substring(0, pathToMappinginfo.lastIndexOf("/"))+"/EBNA2/E2_er/rep1/star.bam";
	public static final String E2er_rep1_tophat = pathToMappinginfo.substring(0, pathToMappinginfo.lastIndexOf("/"))+"/EBNA2/E2_er/rep1/tophat2.bam";

	public static final String E2er_rep2_contextmap = pathToMappinginfo.substring(0, pathToMappinginfo.lastIndexOf("/"))+"/EBNA2/E2_er/rep2/contextmap.bam";
	public static final String E2er_rep2_hisat = pathToMappinginfo.substring(0, pathToMappinginfo.lastIndexOf("/"))+"/EBNA2/E2_er/rep2/hisat.bam";
	public static final String E2er_rep2_star = pathToMappinginfo.substring(0, pathToMappinginfo.lastIndexOf("/"))+"/EBNA2/E2_er/rep2/star.bam";
	public static final String E2er_rep2_tophat = pathToMappinginfo.substring(0, pathToMappinginfo.lastIndexOf("/"))+"/EBNA2/E2_er/rep2/tophat2.bam";

	public static final String E2er_rep3_contextmap = pathToMappinginfo.substring(0, pathToMappinginfo.lastIndexOf("/"))+"/EBNA2/E2_er/rep3/contextmap.bam";
	public static final String E2er_rep3_hisat = pathToMappinginfo.substring(0, pathToMappinginfo.lastIndexOf("/"))+"/EBNA2/E2_er/rep3/hisat.bam";
	public static final String E2er_rep3_star = pathToMappinginfo.substring(0, pathToMappinginfo.lastIndexOf("/"))+"/EBNA2/E2_er/rep3/star.bam";
	public static final String E2er_rep3_tophat = pathToMappinginfo.substring(0, pathToMappinginfo.lastIndexOf("/"))+"/EBNA2/E2_er/rep3/tophat2.bam";
	
	/*
	 * E2_er_chx
	 */
	public static final String E2erchx_rep1_contextmap = pathToMappinginfo.substring(0, pathToMappinginfo.lastIndexOf("/"))+"/EBNA2/E2_er_chx/rep1/contextmap.bam";
	public static final String E2erchx_rep1_hisat = pathToMappinginfo.substring(0, pathToMappinginfo.lastIndexOf("/"))+"/EBNA2/E2_er_chx/rep1/hisat.bam";
	public static final String E2erchx_rep1_star = pathToMappinginfo.substring(0, pathToMappinginfo.lastIndexOf("/"))+"/EBNA2/E2_er_chx/rep1/star.bam";
	public static final String E2erchx_rep1_tophat = pathToMappinginfo.substring(0, pathToMappinginfo.lastIndexOf("/"))+"/EBNA2/E2_er_chx/rep1/tophat2.bam";

	public static final String E2erchx_rep2_contextmap = pathToMappinginfo.substring(0, pathToMappinginfo.lastIndexOf("/"))+"/EBNA2/E2_er_chx/rep2/contextmap.bam";
	public static final String E2erchx_rep2_hisat = pathToMappinginfo.substring(0, pathToMappinginfo.lastIndexOf("/"))+"/EBNA2/E2_er_chx/rep2/hisat.bam";
	public static final String E2erchx_rep2_star = pathToMappinginfo.substring(0, pathToMappinginfo.lastIndexOf("/"))+"/EBNA2/E2_er_chx/rep2/star.bam";
	public static final String E2erchx_rep2_tophat = pathToMappinginfo.substring(0, pathToMappinginfo.lastIndexOf("/"))+"/EBNA2/E2_er_chx/rep2/tophat2.bam";

	public static final String E2erchx_rep3_contextmap = pathToMappinginfo.substring(0, pathToMappinginfo.lastIndexOf("/"))+"/EBNA2/E2_er_chx/rep3/contextmap.bam";
	public static final String E2erchx_rep3_hisat = pathToMappinginfo.substring(0, pathToMappinginfo.lastIndexOf("/"))+"/EBNA2/E2_er_chx/rep3/hisat.bam";
	public static final String E2erchx_rep3_star = pathToMappinginfo.substring(0, pathToMappinginfo.lastIndexOf("/"))+"/EBNA2/E2_er_chx/rep3/star.bam";
	public static final String E2erchx_rep3_tophat = pathToMappinginfo.substring(0, pathToMappinginfo.lastIndexOf("/"))+"/EBNA2/E2_er_chx/rep3/tophat2.bam";
	
	/*
	 * E2_noer
	 */
	public static final String E2noer_rep1_contextmap = pathToMappinginfo.substring(0, pathToMappinginfo.lastIndexOf("/"))+"/EBNA2/E2_noer/rep1/contextmap.bam";
	public static final String E2noer_rep1_hisat = pathToMappinginfo.substring(0, pathToMappinginfo.lastIndexOf("/"))+"/EBNA2/E2_noer/rep1/hisat.bam";
	public static final String E2noer_rep1_star = pathToMappinginfo.substring(0, pathToMappinginfo.lastIndexOf("/"))+"/EBNA2/E2_noer/rep1/star.bam";
	public static final String E2noer_rep1_tophat = pathToMappinginfo.substring(0, pathToMappinginfo.lastIndexOf("/"))+"/EBNA2/E2_noer/rep1/tophat2.bam";

	public static final String E2noer_rep2_contextmap = pathToMappinginfo.substring(0, pathToMappinginfo.lastIndexOf("/"))+"/EBNA2/E2_noer/rep2/contextmap.bam";
	public static final String E2noer_rep2_hisat = pathToMappinginfo.substring(0, pathToMappinginfo.lastIndexOf("/"))+"/EBNA2/E2_noer/rep2/hisat.bam";
	public static final String E2noer_rep2_star = pathToMappinginfo.substring(0, pathToMappinginfo.lastIndexOf("/"))+"/EBNA2/E2_noer/rep2/star.bam";
	public static final String E2noer_rep2_tophat = pathToMappinginfo.substring(0, pathToMappinginfo.lastIndexOf("/"))+"/EBNA2/E2_noer/rep2/tophat2.bam";

	public static final String E2noer_rep3_contextmap = pathToMappinginfo.substring(0, pathToMappinginfo.lastIndexOf("/"))+"/EBNA2/E2_noer/rep3/contextmap.bam";
	public static final String E2noer_rep3_hisat = pathToMappinginfo.substring(0, pathToMappinginfo.lastIndexOf("/"))+"/EBNA2/E2_noer/rep3/hisat.bam";
	public static final String E2noer_rep3_star = pathToMappinginfo.substring(0, pathToMappinginfo.lastIndexOf("/"))+"/EBNA2/E2_noer/rep3/star.bam";
	public static final String E2noer_rep3_tophat = pathToMappinginfo.substring(0, pathToMappinginfo.lastIndexOf("/"))+"/EBNA2/E2_noer/rep3/tophat2.bam";
	
	/*
	 * E2_noer_chx
	 */
	public static final String E2noerchx_rep1_contextmap = pathToMappinginfo.substring(0, pathToMappinginfo.lastIndexOf("/"))+"/EBNA2/E2_noer_chx/rep1/contextmap.bam";
	public static final String E2noerchx_rep1_hisat = pathToMappinginfo.substring(0, pathToMappinginfo.lastIndexOf("/"))+"/EBNA2/E2_noer_chx/rep1/hisat.bam";
	public static final String E2noerchx_rep1_star = pathToMappinginfo.substring(0, pathToMappinginfo.lastIndexOf("/"))+"/EBNA2/E2_noer_chx/rep1/star.bam";
	public static final String E2noerchx_rep1_tophat = pathToMappinginfo.substring(0, pathToMappinginfo.lastIndexOf("/"))+"/EBNA2/E2_noer_chx/rep1/tophat2.bam";

	public static final String E2noerchx_rep2_contextmap = pathToMappinginfo.substring(0, pathToMappinginfo.lastIndexOf("/"))+"/EBNA2/E2_noer_chx/rep2/contextmap.bam";
	public static final String E2noerchx_rep2_hisat = pathToMappinginfo.substring(0, pathToMappinginfo.lastIndexOf("/"))+"/EBNA2/E2_noer_chx/rep2/hisat.bam";
	public static final String E2noerchx_rep2_star = pathToMappinginfo.substring(0, pathToMappinginfo.lastIndexOf("/"))+"/EBNA2/E2_noer_chx/rep2/star.bam";
	public static final String E2noerchx_rep2_tophat = pathToMappinginfo.substring(0, pathToMappinginfo.lastIndexOf("/"))+"/EBNA2/E2_noer_chx/rep2/tophat2.bam";

	public static final String E2noerchx_rep3_contextmap = pathToMappinginfo.substring(0, pathToMappinginfo.lastIndexOf("/"))+"/EBNA2/E2_noer_chx/rep3/contextmap.bam";
	public static final String E2noerchx_rep3_hisat = pathToMappinginfo.substring(0, pathToMappinginfo.lastIndexOf("/"))+"/EBNA2/E2_noer_chx/rep3/hisat.bam";
	public static final String E2noerchx_rep3_star = pathToMappinginfo.substring(0, pathToMappinginfo.lastIndexOf("/"))+"/EBNA2/E2_noer_chx/rep3/star.bam";
	public static final String E2noerchx_rep3_tophat = pathToMappinginfo.substring(0, pathToMappinginfo.lastIndexOf("/"))+"/EBNA2/E2_noer_chx/rep3/tophat2.bam";
	
	
}
