package assignment_4;

import java.io.BufferedReader;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.Vector;

import org.apache.commons.math3.distribution.BetaDistribution;
import org.apache.commons.math3.stat.inference.TTest;

import javafx.util.Pair;
import plotting.LinePlot;

public class RelevantIntronReader {

	ArrayList<RelevantIntron> intronList;
	
	public RelevantIntronReader(){
		
	}
	
	public RelevantIntronReader(String path){
		intronList = new ArrayList<>();
		
		try {
			
			BufferedReader br = new BufferedReader(new FileReader(path));
			
			String line = null;
			
			while((line = br.readLine()) != null){
				String[] lineArray = line.split(":");
				intronList.add(new RelevantIntron(lineArray[0], Integer.parseInt(lineArray[1]), Integer.parseInt(lineArray[2]), false));
			}
			
			br.close();
		} catch (FileNotFoundException e) {
			e.printStackTrace();
		} catch (IOException e) {
			e.printStackTrace();
		}
	}
	
	public HashMap<Pair<Integer, Integer>, double[]> deriveFoldchanges(RelevantIntron intron){
		
		HashMap<Pair<Integer, Integer>, double[]> foldChanges = new HashMap<>();
		
		for (int i = 0; i < 25; i++) {
			for (int j = 0; j < 59; j++) {
		
				BetaDistribution betaDistrib = new BetaDistribution(i + 1.0, j + 1.0);
				
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
				
				foldChanges.put(new Pair<Integer, Integer>(i, j), new double[]{log2fc_80_low, log2fc_80_upp, log2fc_85_low, log2fc_85_upp, log2fc_90_low, log2fc_90_upp, log2fc_95_low, log2fc_95_upp});
		
			}
		}
		
		return foldChanges;
	}
	
	public void executeDETest(RelevantIntron i1, RelevantIntron i2){
		
		TTest test = new TTest();
		
		HashMap<Pair<Integer, Integer>, double[]> foldChanges_i1 = deriveFoldchanges(i1);
		HashMap<Pair<Integer, Integer>, double[]> foldChanges_i2 = deriveFoldchanges(i2);
		
		double[] bigArray_i1 = new double[11800];
		
		int cnt = 0;
		for(double[] dArray : foldChanges_i1.values()){
			for(double d : dArray){
				bigArray_i1[cnt] = d;
				cnt++;
			}
		}
		
		double[] bigArray_i2 = new double[11800];
		
		cnt = 0;
		for(double[] dArray : foldChanges_i2.values()){
			for(double d : dArray){
				bigArray_i2[cnt] = d;
				cnt++;
			}
		}
		
		double result = test.tTest(bigArray_i1, bigArray_i2);
		
		System.out.println("RESULT: "+result);
	}
	
	public void generateCumulativePlot(){
		
		HashMap<String, Integer> countPerGene = new HashMap<>();
		ArrayList<RelevantIntron> tmpList = new ArrayList<>();
		String currentGeneID = "";
		
		for(RelevantIntron i : intronList){
			
			if(!currentGeneID.equals(i.getGeneID())){
				
				for(int ii = 0; ii < tmpList.size(); ii++){
					for(int ij = ii+1; ij < tmpList.size(); ij++){
						if(tmpList.get(ii).overlaps(tmpList.get(ij))){
							if(countPerGene.containsKey(tmpList.get(ii).getGeneID())){
								countPerGene.put(tmpList.get(ii).getGeneID(), countPerGene.get(tmpList.get(ii).getGeneID()) + 1);
							}else{
								countPerGene.put(tmpList.get(ii).getGeneID(), 1);
							}
						}
					}
				}
				
				tmpList = new ArrayList<>();
				currentGeneID = i.getGeneID();
			}
			
			tmpList.add(i);
			
		}
		
		for(int ii = 0; ii < tmpList.size(); ii++){
			for(int ij = ii+1; ij < tmpList.size(); ij++){
				if(tmpList.get(ii).overlaps(tmpList.get(ij))){
					if(countPerGene.containsKey(tmpList.get(ii).getGeneID())){
						countPerGene.put(tmpList.get(ii).getGeneID(), countPerGene.get(tmpList.get(ii).getGeneID()) + 1);
					}else{
						countPerGene.put(tmpList.get(ii).getGeneID(), 1);
					}
				}
			}
		}
		
		HashMap<Integer, Integer> cumList = new HashMap<>();
		
		for(Integer i : countPerGene.values()){
			if(cumList.containsKey(i)){
				cumList.put(i, cumList.get(i) + 1);
			}else{
				cumList.put(i, 1);
			}
		}
		
		Vector<Object> v1 = new Vector<>();
		Vector<Object> v2 = new Vector<>();
		
		v1.addAll(cumList.keySet());
		v2.addAll(cumList.values());
		
		Vector<Vector<Object>> key = new Vector<>();
		Vector<Vector<Object>> value = new Vector<>();
		
		key.add(v1);
		value.add(v2);
		
		Pair<Vector<Vector<Object>>, Vector<Vector<Object>>> pair = new Pair<Vector<Vector<Object>>, Vector<Vector<Object>>>(key, value);
		
		LinePlot lp = new LinePlot(pair, "consistent_introns", "Amount genes", "Amount overlapping introns", 0, 0, false, false);
		lp.filename = "consistent_introns.png";
		lp.showLegend = false;
		lp.plot();
	}
	
	public static void main(String[] args) {
		
		RelevantIntronReader r = new RelevantIntronReader();
		r.executeDETest(new RelevantIntron("", 0, 0, false), new RelevantIntron("", 0, 0, false));
		
	}
}
