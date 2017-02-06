package assignment_4;

import java.io.BufferedReader;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.io.InputStreamReader;
import java.math.BigDecimal;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Map.Entry;
import java.util.TreeMap;

import debugStuff.DebugMessageFactory;
import io.ConfigHelper;
import io.ConfigReader;
import io.ExternalFileWriter;
import javafx.util.Pair;

public class DataRetrieverTask6 {

	
	public void first(String filepath){
		
		String[] geneIds = new String[3];
		double[] smallest = new double[3];
		for (int i = 0; i < smallest.length; i++) {
			smallest[i] = 2.0;
		}
		
		try {
			
			BufferedReader br = new BufferedReader(new FileReader(filepath));
			
			String line = null;
			br.readLine();
			
			while((line = br.readLine()) != null){
				
				String[] currentLine = line.split("\t");
				
				String geneID = currentLine[0];
				double transcriptomicNRP = Integer.parseInt(currentLine[4])*1.0;
				double highestTransNRP = Integer.parseInt(currentLine[6])*1.0;
				
				if(transcriptomicNRP == 0 || highestTransNRP == 0){
					continue;
				}
				
				double proportion = highestTransNRP/transcriptomicNRP;
				
				for (int i = 0; i < smallest.length; i++) {
					if(proportion < smallest[i]){
						for(int j = smallest.length-1; j > i; j--){
							smallest[j] = smallest[j-1];
							geneIds[j] = geneIds[j-1];
						}
						smallest[i] = proportion;
						geneIds[i] = geneID;
						break;
					}
				}
			}
			
			br.close();
		} catch (FileNotFoundException e) {
			e.printStackTrace();
		} catch (IOException e) {
			e.printStackTrace();
		}

		for (int i = 0; i < smallest.length; i++) {
			getRegion(geneIds[i], "/home/proj/biosoft/praktikum/genprakt-ws16/gtf/", "1");
		}
	}
	
	public void second(){
		
		
		String er_rep1 = "/home/proj/biocluster/praktikum/genprakt-ws16/KleinPost/Solution3/EBV/EBNA2/E2_er/rep1/contextmap/contextmap_trFPKM.txt";
		String er_rep2 = "/home/proj/biocluster/praktikum/genprakt-ws16/KleinPost/Solution3/EBV/EBNA2/E2_er/rep2/contextmap/contextmap_trFPKM.txt";
		String er_rep3 = "/home/proj/biocluster/praktikum/genprakt-ws16/KleinPost/Solution3/EBV/EBNA2/E2_er/rep3/contextmap/contextmap_trFPKM.txt";
		
		String noer_rep1 = "/home/proj/biocluster/praktikum/genprakt-ws16/KleinPost/Solution3/EBV/EBNA2/E2_noer/rep1/contextmap/contextmap_trFPKM.txt";
		String noer_rep2 = "/home/proj/biocluster/praktikum/genprakt-ws16/KleinPost/Solution3/EBV/EBNA2/E2_noer/rep2/contextmap/contextmap_trFPKM.txt";
		String noer_rep3 = "/home/proj/biocluster/praktikum/genprakt-ws16/KleinPost/Solution3/EBV/EBNA2/E2_noer/rep3/contextmap/contextmap_trFPKM.txt";
		
		String er_chx_rep1 = "/home/proj/biocluster/praktikum/genprakt-ws16/KleinPost/Solution3/EBV/EBNA2/E2_er_chx/rep1/contextmap/contextmap_trFPKM.txt";
		String er_chx_rep2 = "/home/proj/biocluster/praktikum/genprakt-ws16/KleinPost/Solution3/EBV/EBNA2/E2_er_chx/rep2/contextmap/contextmap_trFPKM.txt";
		String er_chx_rep3 = "/home/proj/biocluster/praktikum/genprakt-ws16/KleinPost/Solution3/EBV/EBNA2/E2_er_chx/rep3/contextmap/contextmap_trFPKM.txt";
		
		String noer_chx_rep1 = "/home/proj/biocluster/praktikum/genprakt-ws16/KleinPost/Solution3/EBV/EBNA2/E2_noer_chx/rep1/contextmap/contextmap_trFPKM.txt";
		String noer_chx_rep2 = "/home/proj/biocluster/praktikum/genprakt-ws16/KleinPost/Solution3/EBV/EBNA2/E2_noer_chx/rep2/contextmap/contextmap_trFPKM.txt";
		String noer_chx_rep3 = "/home/proj/biocluster/praktikum/genprakt-ws16/KleinPost/Solution3/EBV/EBNA2/E2_noer_chx/rep3/contextmap/contextmap_trFPKM.txt";
		
		ArrayList<String> paths = new ArrayList<>();
		paths.add(er_rep1);
		paths.add(er_rep2);
		paths.add(er_rep3);
		paths.add(noer_rep1);
		paths.add(noer_rep2);
		paths.add(noer_rep3);
		paths.add(er_chx_rep1);
		paths.add(er_chx_rep2);
		paths.add(er_chx_rep3);
		paths.add(noer_chx_rep1);
		paths.add(noer_chx_rep2);
		paths.add(noer_chx_rep3);
		
		TreeMap<String, Pair<Double, Double>> rpkmMap = new TreeMap<>();
		
		int counter = 1;
		
		for(String path : paths){
			
			DebugMessageFactory.printInfoDebugMessage(ConfigReader.DEBUG_MODE, "["+counter+"/"+paths.size()+"]");
			
			try {
				
				BufferedReader br = new BufferedReader(new FileReader(path));
				
				String line = br.readLine();
				
				while((line = br.readLine()) != null){
					
					String[] lineArray = line.split("\t");
					String geneID = lineArray[0];
					double rpkm = Double.parseDouble(lineArray[1]);
					
					if(rpkmMap.containsKey(geneID)){
						
						/* check max und min rpkm */
						Pair<Double, Double> pair = rpkmMap.get(geneID);
	
						/* is smaller */
						if(rpkm < pair.getKey()){
							rpkmMap.put(geneID, new Pair<Double, Double>(rpkm, pair.getValue()));
						}
						/* is bigger */
						else if(rpkm > pair.getValue()){
							rpkmMap.put(geneID, new Pair<Double, Double>(pair.getKey(), rpkm));
						}
						
					}else{
						rpkmMap.put(geneID, new Pair<Double, Double>(rpkm, rpkm));
					}
					
				}
				
				br.close();
			} catch (FileNotFoundException e) {
				e.printStackTrace();
			} catch (IOException e) {
				e.printStackTrace();
			}
			
			counter++;
		}
		
		double[] diffArray = new double[3];
		String[] geneArray = new String[3];
		
		for (int i = 0; i < diffArray.length; i++) {
			diffArray[i] = Double.MIN_VALUE;
		}
		
		for(Entry<String, Pair<Double, Double>> entry : rpkmMap.entrySet()){
			
			double diff = entry.getValue().getValue()-entry.getValue().getKey();
			
			for (int i = 0; i < diffArray.length; i++) {
				if(diff > diffArray[i]){
					for(int j = diffArray.length-1; j > i; j--){
						diffArray[j] = diffArray[j-1];
						geneArray[j] = geneArray[j-1];
					}
					diffArray[i] = diff;
					geneArray[i] = entry.getKey();
					break;
				}
			}
			
		}
		
		for (int i = 0; i < geneArray.length; i++) {
			getRegion(geneArray[i], "/home/proj/biosoft/praktikum/genprakt-ws16/gtf/", "2");
		}
		
	}
	
	public void third(){
		
		double[] smallest = new double[3];
		String[] geneIds = new String[3];
		
		for (int i = 0; i < smallest.length; i++) {
			smallest[i] = Double.MAX_VALUE;
		}
		
		try {
			
			BufferedReader br = new BufferedReader(new FileReader("/home/proj/biocluster/praktikum/genprakt-ws16/KleinPost/Solution4/task_6/er_chx_vs_no_er_chx_contextmap_deseq.out"));
			
			String line = br.readLine();
			
			while((line = br.readLine()) != null){
				
				String[] lineArray = line.split("\t");
				double pval = Double.valueOf(lineArray[2]).doubleValue();
				
				for (int i = 0; i < smallest.length; i++) {
					if(pval < smallest[i]){
						for(int j = smallest.length-1; j > i; j--){
							smallest[j] = smallest[j-1];
							geneIds[j] = geneIds[j-1];
						}
						smallest[i] = pval;
						geneIds[i] = lineArray[0];
						break;
					}
				}
				
			}
			
			br.close();
		} catch (FileNotFoundException e) {
			e.printStackTrace();
		} catch (IOException e) {
			e.printStackTrace();
		}
		
		for (int i = 0; i < geneIds.length; i++) {
			getRegion(geneIds[i], "/home/proj/biosoft/praktikum/genprakt-ws16/gtf/", "3");
		}
		
	}
	
	public void fourth(){
		
		String er_rep1 = "/home/proj/biocluster/praktikum/genprakt-ws16/KleinPost/Solution3/EBV/EBNA2/E2_er/rep1/contextmap/contextmap_splitinfos.tsv";
		String er_rep2 = "/home/proj/biocluster/praktikum/genprakt-ws16/KleinPost/Solution3/EBV/EBNA2/E2_er/rep2/contextmap/contextmap_splitinfos.tsv";
		String er_rep3 = "/home/proj/biocluster/praktikum/genprakt-ws16/KleinPost/Solution3/EBV/EBNA2/E2_er/rep3/contextmap/contextmap_splitinfos.tsv";
		String noer_rep1 = "/home/proj/biocluster/praktikum/genprakt-ws16/KleinPost/Solution3/EBV/EBNA2/E2_noer/rep1/contextmap/contextmap_splitinfos.tsv";
		String noer_rep2 = "/home/proj/biocluster/praktikum/genprakt-ws16/KleinPost/Solution3/EBV/EBNA2/E2_noer/rep2/contextmap/contextmap_splitinfos.tsv";
		String noer_rep3 = "/home/proj/biocluster/praktikum/genprakt-ws16/KleinPost/Solution3/EBV/EBNA2/E2_noer/rep3/contextmap/contextmap_splitinfos.tsv";
		
		GenesplitTSVReader reader = new GenesplitTSVReader();
		
		ArrayList<String> firstCond = new ArrayList<>();
		ArrayList<String> secondCond = new ArrayList<>();
		firstCond.add(er_rep1);
		firstCond.add(er_rep2);
		firstCond.add(er_rep3);
		secondCond.add(noer_rep1);
		secondCond.add(noer_rep2);
		secondCond.add(noer_rep3);
		
		reader.readFiles(firstCond, secondCond);
		reader.writeToFile(new ConfigHelper().getDefaultOutputPath()+"consIntrons.txt");
		
	}
	
	
	public void getRegion(String geneID, String filepath, String taskNr){
		
		try {
			
			Process pb = new ProcessBuilder("grep", "-n", ".*gene.*gene_id.*"+geneID, filepath+"Homo_sapiens.GRCh37.75.gtf").start();
			
			BufferedReader br = new BufferedReader(new InputStreamReader(pb.getInputStream()));
			String line = br.readLine();
			br.close();
			
			String[] lineArray = line.split("\t");
			
			String chrID = lineArray[0].split(":")[1];
			int start = Integer.parseInt(lineArray[3]);
			int stop = Integer.parseInt(lineArray[4]);
			
			ExternalFileWriter ext = new ExternalFileWriter();
			ext.openWriter(new ConfigHelper().getDefaultOutputPath()+taskNr+"_region_"+geneID+".bed");
			ext.writeToWriter(chrID+":"+start+"-"+stop);
			ext.closeWriter();
			
		} catch (IOException e) {
			e.printStackTrace();
		}
	}
	
	public static void main(String[] args) {
		
		DataRetrieverTask6 dr = new DataRetrieverTask6();
//		dr.first("/home/proj/biocluster/praktikum/genprakt-ws16/KleinPost/Solution3/EBV/EBNA2/E2_noer/rep1/contextmap/contextmap_counts.tsv");
//		dr.second();
//		dr.third();
//		dr.fourth();
		dr.getRegion("ENSG00000196459", "/home/proj/biosoft/praktikum/genprakt-ws16/gtf/", "4");
		
	}
}
