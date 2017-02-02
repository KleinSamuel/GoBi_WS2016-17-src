package assignment_4;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;

import org.apache.commons.math3.distribution.BetaDistribution;
import org.apache.commons.math3.stat.inference.TTest;

import io.ConfigHelper;
import io.ExternalFileWriter;

public class GenesplitTSVReader {
	
	ArrayList<RelevantIntron> intronList;
	ArrayList<RelevantIntron> tempIntronList;
	
	public void read(String filePath){
		
		try {
			
			BufferedReader br = new BufferedReader(new FileReader(filePath));
			String line = null;
			
			String[] header = br.readLine().split("\t");
			
			ArrayList<RelevantIntron> tmpList = new ArrayList<>();
			tempIntronList = new ArrayList<>();
			String currentGeneID = "";
			
			while((line = br.readLine()) != null){
				
				String[] currentLine = line.split("\t");
				
				String geneID = currentLine[0];
				String chrID = currentLine[1];
				boolean isOnNegativeStrand = currentLine[2].equals("+");
				int start = Integer.parseInt(currentLine[3]);
				int stop = Integer.parseInt(currentLine[4]);
				
				if(!currentGeneID.equals(geneID)){
					
					/* intron needs to overlap another intron so there must be at least 2 introns */
					if(tmpList.size() > 1){
						processList(tmpList);
					}
					tmpList = new ArrayList<>();
					currentGeneID = geneID;
				}
				
				RelevantIntron tmp_intron = new RelevantIntron(geneID, start, stop, isOnNegativeStrand);
				
				tmpList.add(tmp_intron);
			}
			
			processList(tmpList);
			
			System.out.println("TEMPLIST SIZE:\t"+tempIntronList.size());
			
			processFinal(tempIntronList);

			System.out.println("FINALLIST SIZE:\t"+intronList.size());
			
			br.close();
			
		} catch (FileNotFoundException e) {
			e.printStackTrace();
		} catch (IOException e) {
			e.printStackTrace();
		}
		
	}
	
	public void processList(ArrayList<RelevantIntron> list){
		ArrayList<RelevantIntron> out = new ArrayList<>();
		
		for(int i = 0; i < list.size(); i++){
			
			RelevantIntron i1 = list.get(i);
			
			for(int j = i+1; j < list.size(); j++){
				
				RelevantIntron i2 = list.get(j);
				
				if(i1.overlaps(i2)){
					if(!out.contains(i1)){
						i1.addOverlappingIntron(i2);
						out.add(i1);
					}
					if(!out.contains(i2)){
						i2.addOverlappingIntron(i1);
						out.add(i2);
					}
				}
			}
		}
		
		tempIntronList.addAll(out);
	}
	
	public void processFinal(ArrayList<RelevantIntron> list){
		
		if(intronList == null){
			intronList = new ArrayList<>();
			intronList.addAll(tempIntronList);
			return;
		}
		
		ArrayList<Integer> toRemoveIndices = new ArrayList<>();
		
		for(int i = 0; i < intronList.size(); i++){
			
			RelevantIntron i1 = intronList.get(i);
			
			boolean flag = false;
			for(RelevantIntron i2 : list){
				if(i1.isTheSame(i2)){
					flag = true;
				}
			}
			
			if(!flag){
				toRemoveIndices.add(i);
			}
			
		}
		
		for(Integer i : toRemoveIndices){
			intronList.remove(i);
		}
		
	}
	
	public void writeToFile(String filePath){
		
		try {
			
			BufferedWriter bw = new BufferedWriter(new FileWriter(filePath));
			
			for(RelevantIntron i : intronList){
				bw.write(i.getGeneID()+":"+i.getStart()+":"+i.getStop()+"\n");
			}
			
			bw.flush();
			bw.close();
			
		} catch (IOException e) {
			e.printStackTrace();
		}
		
	}
	
	public RelevantIntron getIntronFromIntronlist(String geneID){
		for(RelevantIntron intron : intronList){
			if(intron.getGeneID().equals(geneID)){
				return intron;
			}
		}
		return null;
	}
	
	public double[] getFoldChanges(int count1, int count2){
		BetaDistribution betaDistrib = new BetaDistribution(count1 + 1.0, count2 + 1.0);
		
		/* 80 % credible interval */
		double cred80 = 1 - 0.8;
		double low_bound_80 = betaDistrib.inverseCumulativeProbability(cred80 / 2.0);
		double upp_bound_80 = betaDistrib.inverseCumulativeProbability(1 - cred80 / 2.0);
		double log2fc_80_low = Math.log(1.0 / low_bound_80 - 1.0) / Math.log(2);
		double log2fc_80_upp = Math.log(1.0 / upp_bound_80 - 1.0) / Math.log(2);
		
		/* 85 % credible interval */
		double cred85 = 1 - 0.85;
		double low_bound_85 = betaDistrib.inverseCumulativeProbability(cred85 / 2.0);
		double upp_bound_85 = betaDistrib.inverseCumulativeProbability(1 - cred85 / 2.0);
		double log2fc_85_low = Math.log(1.0 / low_bound_85 - 1.0) / Math.log(2);
		double log2fc_85_upp = Math.log(1.0 / upp_bound_85 - 1.0) / Math.log(2);
		
		/* 90 % credible interval */
		double cred90 = 1 - 0.9;
		double low_bound_90 = betaDistrib.inverseCumulativeProbability(cred90 / 2.0);
		double upp_bound_90 = betaDistrib.inverseCumulativeProbability(1 - cred90 / 2.0);
		double log2fc_90_low = Math.log(1.0 / low_bound_90 - 1.0) / Math.log(2);
		double log2fc_90_upp = Math.log(1.0 / upp_bound_90 - 1.0) / Math.log(2);
		
		/* 95 % credible interval */
		double cred95 = 1 - 0.95;
		double low_bound_95 = betaDistrib.inverseCumulativeProbability(cred95 / 2.0);
		double upp_bound_95 = betaDistrib.inverseCumulativeProbability(1 - cred95 / 2.0);
		double log2fc_95_low = Math.log(1.0 / low_bound_95 - 1.0) / Math.log(2);
		double log2fc_95_upp = Math.log(1.0 / upp_bound_95 - 1.0) / Math.log(2);
		
		return new double[]{log2fc_80_low, log2fc_80_upp, log2fc_85_low, log2fc_85_upp, log2fc_90_low, log2fc_90_upp, log2fc_95_low, log2fc_95_upp};
	}
	
	public int get_nfragsKidney(String geneID, String kidneyPath){
		
		try {
			BufferedReader br = new BufferedReader(new FileReader(kidneyPath));
			String lineKidney = null;
			
			while((lineKidney = br.readLine()) != null){
				
				String[] currentLineKidney = lineKidney.split("\t");
				String kidenyID = currentLineKidney[0];
				
				if(kidenyID.equals(geneID)){
					Integer nfragsKidney = Integer.parseInt(currentLineKidney[6]);
					return nfragsKidney;
				}
			}
			
			br.close();
		
		} catch (FileNotFoundException e){
			e.printStackTrace();
		} catch (IOException e) {
			e.printStackTrace();
		}
		
		return 0;
	}
	
	public int get_nfragsThyroid(String geneID, String thyroidPath){
			
		try {
			BufferedReader br2 = new BufferedReader(new FileReader(thyroidPath));
			String lineThyroid = null;
			
			while((lineThyroid = br2.readLine()) != null){
				
				String[] currentLineThyroid = lineThyroid.split("\t");
				String thyroidID = currentLineThyroid[0];
				
				if(thyroidID.equals(geneID)){
					Integer nfragsThyroid = Integer.parseInt(currentLineThyroid[6]);
					return nfragsThyroid;
				}
			}
			
			br2.close();
		
		} catch (FileNotFoundException e){
			e.printStackTrace();
		} catch (IOException e) {
			e.printStackTrace();
		}
		
		return 0;
	}
	
	public void applySteps_mode_1_nfrags(ArrayList<String> kidneyPaths, ArrayList<String> thyroidPaths){
		
		ExternalFileWriter extFW = new ExternalFileWriter();
		extFW.openWriter(new ConfigHelper().getDefaultOutputPath()+"nfrags_m1_results.txt");
		
		extFW.writeToWriter("intron_id1\tintron_id2\tmean(intron_id1_log2fc)\tSD(intron_id1_log2fc)\tmean(intron_id2_log2fc)\tSD(intron_id2_log2fc)\tpval\ttstat");
		extFW.writeToWriter("\n");
		
		TTest ttest = new TTest();
		
		int counter = 0;
		for(RelevantIntron intron : intronList){
			
			int currentPositionInValues = 0;
			double[] values_intron_main = new double[11800];
			
			for(String kidney : kidneyPaths){
				
				int nfragsKidney = get_nfragsKidney(intron.getGeneID(), kidney);
					
				for(String thyroid : thyroidPaths){
					
					int nfragsThyroid = get_nfragsThyroid(intron.getGeneID(), thyroid);
					
					for(double d : getFoldChanges(nfragsKidney, nfragsThyroid)){
						values_intron_main[currentPositionInValues] = d;
						currentPositionInValues++;
					}
					
				}
					
			}
			
			/* get intron overlapping pair */
			for(RelevantIntron overlapIntron : intron.getOverlappingIntrons()){
				
				if(intronList.indexOf(overlapIntron) < counter){
					continue;
				}
				
				int currentPositionInValuesOverlap = 0;
				double[] values_intron_overlap = new double[11800];
				
				for(String kidney : kidneyPaths){
					
					int nfragsKidney = get_nfragsKidney(overlapIntron.getGeneID(), kidney);
						
					for(String thyroid : thyroidPaths){
						
						int nfragsThyroid = get_nfragsThyroid(overlapIntron.getGeneID(), thyroid);
						
						for(double d : getFoldChanges(nfragsKidney, nfragsThyroid)){
							values_intron_overlap[currentPositionInValuesOverlap] = d;
							currentPositionInValuesOverlap++;
						}
						
					}
						
				}
				
				
				
				double pval = ttest.tTest(values_intron_main, values_intron_overlap);
				double tstat = ttest.t(values_intron_main, values_intron_overlap);
				
				extFW.writeToWriter(intron.getGeneID()+"\t");
				extFW.writeToWriter(overlapIntron.getGeneID()+"\t");
				extFW.writeToWriter("mean_i1\t");
				extFW.writeToWriter("SD_i1\t");
				extFW.writeToWriter("mean_i2\t");
				extFW.writeToWriter("SD_i2\t");
				extFW.writeToWriter(pval+"\t");
				extFW.writeToWriter(tstat+"\n");
			}
			
			counter++;
		}
		
		extFW.closeWriter();
	}
	
	public static void main(String[] args) {
		
		ArrayList<String> kidneyPaths = new ArrayList<>();
		kidneyPaths.add("H:\\GOBI\\a4\\t3\\kidney\\p1_genesplit.tsv");
		kidneyPaths.add("H:\\GOBI\\a4\\t3\\kidney\\p2_genesplit.tsv");
		kidneyPaths.add("H:\\GOBI\\a4\\t3\\kidney\\p3_genesplit.tsv");
//		kidneyPaths.add("H:\\GOBI\\a4\\t3\\kidney\\p4_genesplit.tsv");
//		kidneyPaths.add("H:\\GOBI\\a4\\t3\\kidney\\p5_genesplit.tsv");
//		kidneyPaths.add("H:\\GOBI\\a4\\t3\\kidney\\p6_genesplit.tsv");
//		kidneyPaths.add("H:\\GOBI\\a4\\t3\\kidney\\p7_genesplit.tsv");
//		kidneyPaths.add("H:\\GOBI\\a4\\t3\\kidney\\p8_genesplit.tsv");
//		kidneyPaths.add("H:\\GOBI\\a4\\t3\\kidney\\p9_genesplit.tsv");
//		kidneyPaths.add("H:\\GOBI\\a4\\t3\\kidney\\p10_genesplit.tsv");
//		kidneyPaths.add("H:\\GOBI\\a4\\t3\\kidney\\p11_genesplit.tsv");
//		kidneyPaths.add("H:\\GOBI\\a4\\t3\\kidney\\p12_genesplit.tsv");
//		kidneyPaths.add("H:\\GOBI\\a4\\t3\\kidney\\p13_genesplit.tsv");
//		kidneyPaths.add("H:\\GOBI\\a4\\t3\\kidney\\p14_genesplit.tsv");
//		kidneyPaths.add("H:\\GOBI\\a4\\t3\\kidney\\p15_genesplit.tsv");
//		kidneyPaths.add("H:\\GOBI\\a4\\t3\\kidney\\p16_genesplit.tsv");
//		kidneyPaths.add("H:\\GOBI\\a4\\t3\\kidney\\p17_genesplit.tsv");
//		kidneyPaths.add("H:\\GOBI\\a4\\t3\\kidney\\p18_genesplit.tsv");
//		kidneyPaths.add("H:\\GOBI\\a4\\t3\\kidney\\p19_genesplit.tsv");
//		kidneyPaths.add("H:\\GOBI\\a4\\t3\\kidney\\p20_genesplit.tsv");
//		kidneyPaths.add("H:\\GOBI\\a4\\t3\\kidney\\p21_genesplit.tsv");
//		kidneyPaths.add("H:\\GOBI\\a4\\t3\\kidney\\p22_genesplit.tsv");
//		kidneyPaths.add("H:\\GOBI\\a4\\t3\\kidney\\p23_genesplit.tsv");
//		kidneyPaths.add("H:\\GOBI\\a4\\t3\\kidney\\p24_genesplit.tsv");
//		kidneyPaths.add("H:\\GOBI\\a4\\t3\\kidney\\p25_genesplit.tsv");
		
		ArrayList<String> thyroidPaths = new ArrayList<>();
		thyroidPaths.add("H:\\GOBI\\a4\\t3\\thyroid\\p1_genesplit.tsv");
		thyroidPaths.add("H:\\GOBI\\a4\\t3\\thyroid\\p2_genesplit.tsv");
		thyroidPaths.add("H:\\GOBI\\a4\\t3\\thyroid\\p3_genesplit.tsv");
		
		GenesplitTSVReader r = new GenesplitTSVReader();
		r.read("H:\\GOBI\\a4\\t3\\kidney\\p1_genesplit.tsv");
		r.read("H:\\GOBI\\a4\\t3\\kidney\\p2_genesplit.tsv");
		r.read("H:\\GOBI\\a4\\t3\\kidney\\p3_genesplit.tsv");
//		r.read("H:\\GOBI\\a4\\t3\\kidney\\p4_genesplit.tsv");
//		r.read("H:\\GOBI\\a4\\t3\\kidney\\p5_genesplit.tsv");
//		r.read("H:\\GOBI\\a4\\t3\\kidney\\p6_genesplit.tsv");
//		r.read("H:\\GOBI\\a4\\t3\\kidney\\p7_genesplit.tsv");
//		r.read("H:\\GOBI\\a4\\t3\\kidney\\p8_genesplit.tsv");
//		r.read("H:\\GOBI\\a4\\t3\\kidney\\p9_genesplit.tsv");
//		r.read("H:\\GOBI\\a4\\t3\\kidney\\p10_genesplit.tsv");
//		r.read("H:\\GOBI\\a4\\t3\\kidney\\p11_genesplit.tsv");
//		r.read("H:\\GOBI\\a4\\t3\\kidney\\p12_genesplit.tsv");
//		r.read("H:\\GOBI\\a4\\t3\\kidney\\p13_genesplit.tsv");
//		r.read("H:\\GOBI\\a4\\t3\\kidney\\p14_genesplit.tsv");
//		r.read("H:\\GOBI\\a4\\t3\\kidney\\p15_genesplit.tsv");
//		r.read("H:\\GOBI\\a4\\t3\\kidney\\p16_genesplit.tsv");
//		r.read("H:\\GOBI\\a4\\t3\\kidney\\p17_genesplit.tsv");
//		r.read("H:\\GOBI\\a4\\t3\\kidney\\p18_genesplit.tsv");
//		r.read("H:\\GOBI\\a4\\t3\\kidney\\p19_genesplit.tsv");
//		r.read("H:\\GOBI\\a4\\t3\\kidney\\p20_genesplit.tsv");
//		r.read("H:\\GOBI\\a4\\t3\\kidney\\p21_genesplit.tsv");
//		r.read("H:\\GOBI\\a4\\t3\\kidney\\p22_genesplit.tsv");
//		r.read("H:\\GOBI\\a4\\t3\\kidney\\p23_genesplit.tsv");
//		r.read("H:\\GOBI\\a4\\t3\\kidney\\p24_genesplit.tsv");
//		r.read("H:\\GOBI\\a4\\t3\\kidney\\p25_genesplit.tsv");
		
		ArrayList<RelevantIntron> firstList = r.intronList;
		
		System.out.println("SIZE LIST 1: "+r.intronList.size());
		
		r = new GenesplitTSVReader();
		r.read("H:\\GOBI\\a4\\t3\\thyroid\\p1_genesplit.tsv");
		r.read("H:\\GOBI\\a4\\t3\\thyroid\\p2_genesplit.tsv");
		r.read("H:\\GOBI\\a4\\t3\\thyroid\\p3_genesplit.tsv");
		
		System.out.println("SIZE LIST 2: "+r.intronList.size());
		
		r.intronList.addAll(firstList);
		
		r.applySteps_mode_1_nfrags(kidneyPaths, thyroidPaths);
//		r.writeToFile("H:\\GOBI\\a4\\t3\\kidney\\consistent_introns.txt");
		
	}
	
}
