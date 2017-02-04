package assignment_4;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Map.Entry;
import java.util.TreeMap;

import org.apache.commons.math3.distribution.BetaDistribution;
import org.apache.commons.math3.stat.inference.TTest;

import augmentedTree.IntervalTree;
import debugStuff.DebugMessageFactory;
import io.ConfigHelper;
import io.ConfigReader;
import io.ExternalFileWriter;
import javafx.util.Pair;

public class GenesplitTSVReader {
	
	ArrayList<Pair<RelevantIntron, RelevantIntron>> overlapPairs;
	TreeMap<String, TreeMap<Boolean, TreeMap<String, IntervalTree<RelevantIntron>>>> kidneyMap;
	TreeMap<String, TreeMap<Boolean, TreeMap<String, IntervalTree<RelevantIntron>>>> thyroidMap;
	TreeMap<String, TreeMap<Boolean, TreeMap<String, IntervalTree<RelevantIntron>>>> intronMap;
	
	/* additional tasks */
	ArrayList<Pair<RelevantIntron, RelevantIntron>> mostSignificantPairs;
	double[] smallesPVals;
	
	public GenesplitTSVReader(){
		this.intronMap = new TreeMap<>();
		this.overlapPairs = new ArrayList<>();
		this.mostSignificantPairs = new ArrayList<>();
		this.smallesPVals = new double[10];
		for (int i = 0; i < smallesPVals.length; i++) {
			smallesPVals[i] = Integer.MAX_VALUE;
			this.mostSignificantPairs.add(null);
		}
	}
	
	public void readFiles(ArrayList<String> kidneyPaths, ArrayList<String> thyroidPaths){
		
		DebugMessageFactory.printInfoDebugMessage(ConfigReader.DEBUG_MODE, "READING IN KIDNEY FILES...");
		readFirst(kidneyPaths.get(0));
		
		DebugMessageFactory.printNormalDebugMessage(ConfigReader.DEBUG_MODE, "[1/"+kidneyPaths.size()+"]");
		
		for(int i = 1; i < kidneyPaths.size(); i++){
			readOther(kidneyPaths.get(i));
			
			DebugMessageFactory.printNormalDebugMessage(ConfigReader.DEBUG_MODE, "["+(i+1)+"/"+kidneyPaths.size()+"]");
			
		}
		
		removeNonOverlappingIntrons();
		DebugMessageFactory.printInfoDebugMessage(ConfigReader.DEBUG_MODE, "FINAL SIZE: "+countIntronsInMap());
		DebugMessageFactory.printInfoDebugMessage(ConfigReader.DEBUG_MODE, "FINISHED READING IN KIDNEY FILES.");
		
		System.out.println("SIZE OF PAIR LIST: "+overlapPairs.size());
		
		this.kidneyMap = intronMap;
		this.intronMap = new TreeMap<>();
		
		DebugMessageFactory.printInfoDebugMessage(ConfigReader.DEBUG_MODE, "READING IN THYROID FILES...");
		readFirst(thyroidPaths.get(0));
		
		DebugMessageFactory.printNormalDebugMessage(ConfigReader.DEBUG_MODE, "[1/"+thyroidPaths.size()+"]");
		
		for(int i = 1; i < thyroidPaths.size(); i++){
			readOther(thyroidPaths.get(i));
			
			DebugMessageFactory.printNormalDebugMessage(ConfigReader.DEBUG_MODE, "["+(i+1)+"/"+thyroidPaths.size()+"]");
			
		}
		
		removeNonOverlappingIntrons();
		DebugMessageFactory.printInfoDebugMessage(ConfigReader.DEBUG_MODE, "FINAL SIZE: "+countIntronsInMap());
		DebugMessageFactory.printInfoDebugMessage(ConfigReader.DEBUG_MODE, "FINISHED READING IN THYROID FILES.");
		
		System.out.println("SIZE OF PAIR LIST: "+overlapPairs.size());
		
		this.thyroidMap = intronMap;
	}
	
	public void readFirst(String filePath){
		
		try {
			
			BufferedReader br = new BufferedReader(new FileReader(filePath));
			
			String line = null;
			String[] header = br.readLine().split("\t");
			
			while((line = br.readLine()) != null){
				
				String[] currentLine = line.split("\t");
				
				String geneID = currentLine[0];
				String chrID = currentLine[1];
				boolean isOnNegativeStrand = currentLine[2].equals("-");
				int start = Integer.parseInt(currentLine[3]);
				int stop = Integer.parseInt(currentLine[4]);
				
				RelevantIntron currentIntron = new RelevantIntron(geneID, chrID, start, stop, isOnNegativeStrand);
				
				/* contains chromosome */
				if(intronMap.containsKey(chrID)){
					
					TreeMap<Boolean, TreeMap<String, IntervalTree<RelevantIntron>>> strandMap = intronMap.get(chrID);
					
					/* contains strand */
					if(strandMap.containsKey(isOnNegativeStrand)){
						
						TreeMap<String, IntervalTree<RelevantIntron>> geneMap = strandMap.get(isOnNegativeStrand);
						
						/* contains gene */
						if(geneMap.containsKey(geneID)){
							geneMap.get(geneID).add(currentIntron);
						}
						/* does not contain gene */
						else{
							IntervalTree<RelevantIntron> tmpTree = new IntervalTree<>();
							tmpTree.add(currentIntron);
							geneMap.put(geneID, tmpTree);
						}
						
					}
					/* does not contain strand */
					else{
						IntervalTree<RelevantIntron> tmpTree = new IntervalTree<>();
						tmpTree.add(currentIntron);
						TreeMap<String, IntervalTree<RelevantIntron>> tmpGeneMap = new TreeMap<>();
						tmpGeneMap.put(geneID, tmpTree);
						strandMap.put(isOnNegativeStrand, tmpGeneMap);
					}
					
				}
				/* does not contain chromosome */
				else{
					IntervalTree<RelevantIntron> tmpTree = new IntervalTree<>();
					tmpTree.add(currentIntron);
					TreeMap<String, IntervalTree<RelevantIntron>> tmpGeneMap = new TreeMap<>();
					tmpGeneMap.put(geneID, tmpTree);
					TreeMap<Boolean, TreeMap<String, IntervalTree<RelevantIntron>>> tmpStrandMap = new TreeMap<>();
					tmpStrandMap.put(isOnNegativeStrand, tmpGeneMap);
					intronMap.put(chrID, tmpStrandMap);
				}
				
			}
			
			br.close();
			
		} catch (FileNotFoundException e) {
			e.printStackTrace();
		} catch (IOException e) {
			e.printStackTrace();
		}
		
		ArrayList<Pair<String, Pair<Integer, Pair<String, RelevantIntron>>>> toRemove = new ArrayList<>();
		
	}
	
	public void resetToRemove(){
		for(TreeMap<Boolean, TreeMap<String, IntervalTree<RelevantIntron>>> strandMap : intronMap.values()){
			for(TreeMap<String, IntervalTree<RelevantIntron>> geneMap : strandMap.values()){
				for(IntervalTree<RelevantIntron> intervalTree : geneMap.values()){
					for(RelevantIntron intron : intervalTree){
						intron.toRemove = true;
					}
				}
			}
		}
	}
	
	public void addOverlapPairToList(Pair<RelevantIntron, RelevantIntron> pair){
		for(Pair<RelevantIntron, RelevantIntron> p : overlapPairs){
			if(p.getKey().isTheSame(pair.getKey()) && p.getValue().isTheSame(pair.getValue())){
				return;
			}
			if(p.getKey().isTheSame(pair.getValue()) && p.getValue().isTheSame(pair.getKey())){
				return;
			}
		}
		overlapPairs.add(pair);
	}
	
	public void removeUnrelevantIntrons(){
		
		TreeMap<String, TreeMap<Boolean, TreeMap<String, IntervalTree<RelevantIntron>>>> newMap = new TreeMap<>();
		
		for(Entry<String, TreeMap<Boolean, TreeMap<String, IntervalTree<RelevantIntron>>>> chrEntry : intronMap.entrySet()){
			String chr = chrEntry.getKey();
			TreeMap<Boolean, TreeMap<String, IntervalTree<RelevantIntron>>> chrMap = chrEntry.getValue();
			
			for(Entry<Boolean, TreeMap<String, IntervalTree<RelevantIntron>>> strandEntry : chrMap.entrySet()){
				boolean isOnNegative = strandEntry.getKey();
				TreeMap<String, IntervalTree<RelevantIntron>> strandMap = strandEntry.getValue();
				
				for(Entry<String, IntervalTree<RelevantIntron>> geneEntry : strandMap.entrySet()){
					String gene = geneEntry.getKey();
					IntervalTree<RelevantIntron> intervalTree = geneEntry.getValue();
					
					for(RelevantIntron intron : intervalTree){
						if(!intron.toRemove){
							
							/* contains chromosome */
							if(newMap.containsKey(chr)){
								
								TreeMap<Boolean, TreeMap<String, IntervalTree<RelevantIntron>>> strandMap2 = newMap.get(chr);
								
								/* contains strand */
								if(strandMap2.containsKey(isOnNegative)){
									
									TreeMap<String, IntervalTree<RelevantIntron>> geneMap2 = strandMap2.get(isOnNegative);
									
									/* contains gene */
									if(geneMap2.containsKey(gene)){
										geneMap2.get(gene).add(intron);
									}
									/* does not contain gene */
									else{
										IntervalTree<RelevantIntron> tmpTree = new IntervalTree<>();
										tmpTree.add(intron);
										geneMap2.put(gene, tmpTree);
									}
									
								}
								/* does not contain strand */
								else{
									IntervalTree<RelevantIntron> tmpTree = new IntervalTree<>();
									tmpTree.add(intron);
									TreeMap<String, IntervalTree<RelevantIntron>> tmpGeneMap = new TreeMap<>();
									tmpGeneMap.put(gene, tmpTree);
									strandMap2.put(isOnNegative, tmpGeneMap);
								}
								
							}
							/* does not contain chromosome */
							else{
								IntervalTree<RelevantIntron> tmpTree = new IntervalTree<>();
								tmpTree.add(intron);
								TreeMap<String, IntervalTree<RelevantIntron>> tmpGeneMap = new TreeMap<>();
								tmpGeneMap.put(gene, tmpTree);
								TreeMap<Boolean, TreeMap<String, IntervalTree<RelevantIntron>>> tmpStrandMap = new TreeMap<>();
								tmpStrandMap.put(isOnNegative, tmpGeneMap);
								newMap.put(chr, tmpStrandMap);
							}
							
						}
					}
				}
			}
		}
		
		intronMap = newMap;
	}
	
	public void removeNonOverlappingIntrons(){
		
		TreeMap<String, TreeMap<Boolean, TreeMap<String, IntervalTree<RelevantIntron>>>> newMap = new TreeMap<>();
		
		for(Entry<String, TreeMap<Boolean, TreeMap<String, IntervalTree<RelevantIntron>>>> chrEntry : intronMap.entrySet()){
			String chr = chrEntry.getKey();
			TreeMap<Boolean, TreeMap<String, IntervalTree<RelevantIntron>>> chrMap = chrEntry.getValue();
			
			for(Entry<Boolean, TreeMap<String, IntervalTree<RelevantIntron>>> strandEntry : chrMap.entrySet()){
				boolean isOnNegative = strandEntry.getKey();
				TreeMap<String, IntervalTree<RelevantIntron>> strandMap = strandEntry.getValue();
				
				for(Entry<String, IntervalTree<RelevantIntron>> geneEntry : strandMap.entrySet()){
					String gene = geneEntry.getKey();
					IntervalTree<RelevantIntron> intervalTree = geneEntry.getValue();
					
					for(RelevantIntron intron : intervalTree){
						
						intervalTree.getIntervalsIntersecting(intron.getStart(), intron.getStop(), intron.getOverlappingIntrons());
						
						ArrayList<RelevantIntron> newOverlaps = new ArrayList<>();
						
						for(RelevantIntron i : intron.getOverlappingIntrons()){
							if(!(intron.getStart() == i.getStart() && intron.getStop() == i.getStop())){
								newOverlaps.add(i);
								addOverlapPairToList(new Pair<RelevantIntron, RelevantIntron>(intron, i));
							}
						}
						
						intron.setOverlappingIntrons(newOverlaps);
						
						if(intron.getOverlappingIntrons().size() > 0){
							
							/* contains chromosome */
							if(newMap.containsKey(chr)){
								
								TreeMap<Boolean, TreeMap<String, IntervalTree<RelevantIntron>>> strandMap2 = newMap.get(chr);
								
								/* contains strand */
								if(strandMap2.containsKey(isOnNegative)){
									
									TreeMap<String, IntervalTree<RelevantIntron>> geneMap2 = strandMap2.get(isOnNegative);
									
									/* contains gene */
									if(geneMap2.containsKey(gene)){
										geneMap2.get(gene).add(intron);
									}
									/* does not contain gene */
									else{
										IntervalTree<RelevantIntron> tmpTree = new IntervalTree<>();
										tmpTree.add(intron);
										geneMap2.put(gene, tmpTree);
									}
									
								}
								/* does not contain strand */
								else{
									IntervalTree<RelevantIntron> tmpTree = new IntervalTree<>();
									tmpTree.add(intron);
									TreeMap<String, IntervalTree<RelevantIntron>> tmpGeneMap = new TreeMap<>();
									tmpGeneMap.put(gene, tmpTree);
									strandMap2.put(isOnNegative, tmpGeneMap);
								}
								
							}
							/* does not contain chromosome */
							else{
								IntervalTree<RelevantIntron> tmpTree = new IntervalTree<>();
								tmpTree.add(intron);
								TreeMap<String, IntervalTree<RelevantIntron>> tmpGeneMap = new TreeMap<>();
								tmpGeneMap.put(gene, tmpTree);
								TreeMap<Boolean, TreeMap<String, IntervalTree<RelevantIntron>>> tmpStrandMap = new TreeMap<>();
								tmpStrandMap.put(isOnNegative, tmpGeneMap);
								newMap.put(chr, tmpStrandMap);
							}
							
						}
					}
				}
			}
		}
		
		intronMap = newMap;
	}
	
	public void readOther(String filePath){
		
		resetToRemove();
		
		try {
			
			BufferedReader br = new BufferedReader(new FileReader(filePath));
			
			String line = null;
			
			String[] header = br.readLine().split("\t");
			
			while((line = br.readLine()) != null){
				
				String[] currentLine = line.split("\t");
				
				String geneID = currentLine[0];
				String chrID = currentLine[1];
				boolean isOnNegativeStrand = currentLine[2].equals("-");
				int start = Integer.parseInt(currentLine[3]);
				int stop = Integer.parseInt(currentLine[4]);
				
				if(intronMap.containsKey(chrID)){
					TreeMap<Boolean, TreeMap<String, IntervalTree<RelevantIntron>>> strandMap = intronMap.get(chrID);
					if(strandMap.containsKey(isOnNegativeStrand)){
						TreeMap<String, IntervalTree<RelevantIntron>> geneMap = strandMap.get(isOnNegativeStrand);
						if(geneMap.containsKey(geneID)){
							IntervalTree<RelevantIntron> intervalTree = geneMap.get(geneID);
							for(RelevantIntron intron : intervalTree){
								if(intron.sameRegion(start, stop)){
									intron.toRemove = false;
								}
							}
						}
					}
				}
				
			}
			
			br.close();
			
		} catch (FileNotFoundException e) {
			e.printStackTrace();
		} catch (IOException e) {
			e.printStackTrace();
		}
		
		removeUnrelevantIntrons();
	}
	
	public void writeToFile(String filePath){
		
		try {
			
			BufferedWriter bw = new BufferedWriter(new FileWriter(filePath));
			
			for(TreeMap<Boolean, TreeMap<String, IntervalTree<RelevantIntron>>> strandMap : intronMap.values()){
				for(TreeMap<String, IntervalTree<RelevantIntron>> geneMap : strandMap.values()){
					for(IntervalTree<RelevantIntron> intervalTree : geneMap.values()){
						for(RelevantIntron i : intervalTree){
							bw.write(i.getGeneID()+":"+i.getStart()+":"+i.getStop()+"\n");
						}
					}
				}
			}
			
			bw.flush();
			bw.close();
			
		} catch (IOException e) {
			e.printStackTrace();
		}
		
		DebugMessageFactory.printNormalDebugMessage(ConfigReader.DEBUG_MODE, "Written file to: "+filePath);
		
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
	
	public int getDesiredValue(RelevantIntron intron, String path, int mode){
		
		try {
			BufferedReader br = new BufferedReader(new FileReader(path));
			String line = null;
			br.readLine();
			
			while((line = br.readLine()) != null){
				
				String[] currentLine = line.split("\t");
				String geneID = currentLine[0];
				String chromosomeID = currentLine[1];
				boolean isOnNegativeStrand = (currentLine[2].equals("-"));
				int start = Integer.parseInt(currentLine[3]);
				int stop = Integer.parseInt(currentLine[4]);
				
				if(intron.isOnSameGene(geneID) && intron.isOnSameChromosome(chromosomeID) && intron.isOnSameStrand(isOnNegativeStrand) && intron.sameRegion(start, stop)){
					
					switch (mode) {
					case 1:
						return Integer.parseInt(currentLine[6]);
					case 2:
						return Integer.parseInt(currentLine[7]);
					case 3:
						return Integer.parseInt(currentLine[8]);
					case 4:
						return Integer.parseInt(currentLine[9]);

					default:
						break;
					}
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
	
	/**
	 * 1 nfrags
	 * 2 nreads
	 * 3 mm0_nfrags
	 * 4 mm0_nreads
	 * 
	 * 
	 * @param kidneyPaths
	 * @param thyroidPaths
	 * @param mode
	 */
	public void applySteps_mode_1(ArrayList<String> kidneyPaths, ArrayList<String> thyroidPaths, int mode){
		
		ExternalFileWriter extFW = new ExternalFileWriter();
		extFW.openWriter(new ConfigHelper().getDefaultOutputPath()+"nfrags_m1_results.txt");
		
		extFW.writeToWriter("intron_id1\tintron_id2\tmean(intron_id1_log2fc)\tSD(intron_id1_log2fc)\tmean(intron_id2_log2fc)\tSD(intron_id2_log2fc)\tpval\ttstat");
		extFW.writeToWriter("\n");
		
		TTest ttest = new TTest();
		int cnt = 0;
		
		DebugMessageFactory.printNormalDebugMessage(ConfigReader.DEBUG_MODE, "Compute log2fc values.");
		
		for(Pair<RelevantIntron, RelevantIntron> pair : overlapPairs){
			
			DebugMessageFactory.printNormalDebugMessage(ConfigReader.DEBUG_MODE, "["+cnt+"/"+overlapPairs.size()+"]");
			
			/* get values for first intron */
			double[] values = new double[11800];
			int pointerInValues = 0;
			
			if(pair.getKey().values == null){
				
				for(String kidney : kidneyPaths){
					
					int desiredValueKidney = getDesiredValue(pair.getKey(), kidney, mode);
					
					for(String thyroid : thyroidPaths){
						
						int desiredValueThyroid = getDesiredValue(pair.getKey(), thyroid, mode);
						
						for(double d : getFoldChanges(desiredValueKidney, desiredValueThyroid)){
							values[pointerInValues] = d;
							pointerInValues++;
						}
					}
				}
				
				pair.getKey().values = values;
			}
			
			/* get values for overlap intron */
			values = new double[11800];
			pointerInValues = 0;
			
			if(pair.getValue().values == null){
				
				for(String kidney : kidneyPaths){
					
					int desiredValueKidney = getDesiredValue(pair.getValue(), kidney, mode);
					
					for(String thyroid : thyroidPaths){
						
						int desiredValueThyroid = getDesiredValue(pair.getValue(), thyroid, mode);
						
						for(double d : getFoldChanges(desiredValueKidney, desiredValueThyroid)){
							values[pointerInValues] = d;
							pointerInValues++;
						}
					}
				}
				
				pair.getValue().values = values;
			}
				
			double pval = ttest.tTest(pair.getKey().values, pair.getValue().values);
			double tstat = ttest.t(pair.getKey().values, pair.getValue().values);
			
			double mean_1 = getMean(pair.getKey().values);
			double mean_2 = getMean(pair.getValue().values);
			
			double sd_1 = getSD(mean_1, pair.getKey().values);
			double sd_2 = getSD(mean_2, pair.getValue().values);
			
			extFW.writeToWriter(pair.getKey().getGeneID()+"\t");
			extFW.writeToWriter(pair.getValue().getGeneID()+"\t");
			extFW.writeToWriter(mean_1+"\t");
			extFW.writeToWriter(sd_1+"\t");
			extFW.writeToWriter(mean_2+"\t");
			extFW.writeToWriter(sd_2+"\t");
			extFW.writeToWriter(pval+"\t");
			extFW.writeToWriter(tstat+"\n");
			
			sortIntoPValArray(pval, pair);
			
			cnt++;
			
			break;
		}
		
		extFW.closeWriter();
		
		System.out.println("SMALLEST PVALS:");
		for (int i = 0; i < smallesPVals.length; i++) {
			System.out.println(i+":\t"+smallesPVals[i]);
		}
	}
	
	public double getMean(double[] array){
		double mean = 0.0;
		for (int i = 0; i < array.length-1; i+=2) {
			mean += 0.5 * (array[i]+array[i+1]);
		}
		return mean/array.length/2;
	}
	
	public double getSD(double mean, double[] array){
		double temp = 0;
		for (int i = 0; i < array.length; i++) {
			temp += (array[i]-mean)*(array[i]-mean);
		}
		return Math.sqrt(temp/array.length);
	}
	
	public boolean sortIntoPValArray(double pval, Pair<RelevantIntron, RelevantIntron> pair){
		
		boolean flag = false;
		for (int i = 0; i < smallesPVals.length; i++) {
			if(pval < smallesPVals[i]){
				flag = true;
				for(int j = smallesPVals.length-1; j > i; j--){
					smallesPVals[j] = smallesPVals[j-1];
				}
				smallesPVals[i] = pval;
				for(int j = smallesPVals.length-1; j > i; j--){
					mostSignificantPairs.set(j, mostSignificantPairs.get(j-1));
				}
				mostSignificantPairs.set(i, pair);
				break;
			}
		}
		return flag;
	}
	
	public int countIntronsInMap(){
		int amountIntrons = 0;
		
		for(TreeMap<Boolean, TreeMap<String, IntervalTree<RelevantIntron>>> strandMap : intronMap.values()){
			for(TreeMap<String, IntervalTree<RelevantIntron>> geneMap : strandMap.values()){
				for(IntervalTree<RelevantIntron> intervalTree : geneMap.values()){
					for(RelevantIntron intron : intervalTree){
						amountIntrons++;
					}
				}
			}
		}
		
		return amountIntrons;
	}
	
	public static void main(String[] args) {
		
		ArrayList<String> kidneyPaths = new ArrayList<>();
		kidneyPaths.add("H:\\GOBI\\a4\\t3\\kidney\\p1_genesplit.tsv");
		kidneyPaths.add("H:\\GOBI\\a4\\t3\\kidney\\p2_genesplit.tsv");
		kidneyPaths.add("H:\\GOBI\\a4\\t3\\kidney\\p3_genesplit.tsv");
		kidneyPaths.add("H:\\GOBI\\a4\\t3\\kidney\\p4_genesplit.tsv");
		kidneyPaths.add("H:\\GOBI\\a4\\t3\\kidney\\p5_genesplit.tsv");
		kidneyPaths.add("H:\\GOBI\\a4\\t3\\kidney\\p6_genesplit.tsv");
		kidneyPaths.add("H:\\GOBI\\a4\\t3\\kidney\\p7_genesplit.tsv");
		kidneyPaths.add("H:\\GOBI\\a4\\t3\\kidney\\p8_genesplit.tsv");
		kidneyPaths.add("H:\\GOBI\\a4\\t3\\kidney\\p9_genesplit.tsv");
		kidneyPaths.add("H:\\GOBI\\a4\\t3\\kidney\\p10_genesplit.tsv");
		kidneyPaths.add("H:\\GOBI\\a4\\t3\\kidney\\p11_genesplit.tsv");
		kidneyPaths.add("H:\\GOBI\\a4\\t3\\kidney\\p12_genesplit.tsv");
		kidneyPaths.add("H:\\GOBI\\a4\\t3\\kidney\\p13_genesplit.tsv");
		kidneyPaths.add("H:\\GOBI\\a4\\t3\\kidney\\p14_genesplit.tsv");
		kidneyPaths.add("H:\\GOBI\\a4\\t3\\kidney\\p15_genesplit.tsv");
		kidneyPaths.add("H:\\GOBI\\a4\\t3\\kidney\\p16_genesplit.tsv");
		kidneyPaths.add("H:\\GOBI\\a4\\t3\\kidney\\p17_genesplit.tsv");
		kidneyPaths.add("H:\\GOBI\\a4\\t3\\kidney\\p18_genesplit.tsv");
		kidneyPaths.add("H:\\GOBI\\a4\\t3\\kidney\\p19_genesplit.tsv");
		kidneyPaths.add("H:\\GOBI\\a4\\t3\\kidney\\p20_genesplit.tsv");
		kidneyPaths.add("H:\\GOBI\\a4\\t3\\kidney\\p21_genesplit.tsv");
		kidneyPaths.add("H:\\GOBI\\a4\\t3\\kidney\\p22_genesplit.tsv");
		kidneyPaths.add("H:\\GOBI\\a4\\t3\\kidney\\p23_genesplit.tsv");
		kidneyPaths.add("H:\\GOBI\\a4\\t3\\kidney\\p24_genesplit.tsv");
		kidneyPaths.add("H:\\GOBI\\a4\\t3\\kidney\\p25_genesplit.tsv");
		
		ArrayList<String> thyroidPaths = new ArrayList<>();
		thyroidPaths.add("H:\\GOBI\\a4\\t3\\thyroid\\p1_genesplit.tsv");
		thyroidPaths.add("H:\\GOBI\\a4\\t3\\thyroid\\p2_genesplit.tsv");
		thyroidPaths.add("H:\\GOBI\\a4\\t3\\thyroid\\p3_genesplit.tsv");
		thyroidPaths.add("H:\\GOBI\\a4\\t3\\thyroid\\p4_genesplit.tsv");
		thyroidPaths.add("H:\\GOBI\\a4\\t3\\thyroid\\p5_genesplit.tsv");
		thyroidPaths.add("H:\\GOBI\\a4\\t3\\thyroid\\p6_genesplit.tsv");
		thyroidPaths.add("H:\\GOBI\\a4\\t3\\thyroid\\p7_genesplit.tsv");
		thyroidPaths.add("H:\\GOBI\\a4\\t3\\thyroid\\p8_genesplit.tsv");
		thyroidPaths.add("H:\\GOBI\\a4\\t3\\thyroid\\p9_genesplit.tsv");
		thyroidPaths.add("H:\\GOBI\\a4\\t3\\thyroid\\p10_genesplit.tsv");
		thyroidPaths.add("H:\\GOBI\\a4\\t3\\thyroid\\p11_genesplit.tsv");
		thyroidPaths.add("H:\\GOBI\\a4\\t3\\thyroid\\p12_genesplit.tsv");
		thyroidPaths.add("H:\\GOBI\\a4\\t3\\thyroid\\p13_genesplit.tsv");
		thyroidPaths.add("H:\\GOBI\\a4\\t3\\thyroid\\p14_genesplit.tsv");
		thyroidPaths.add("H:\\GOBI\\a4\\t3\\thyroid\\p15_genesplit.tsv");
		thyroidPaths.add("H:\\GOBI\\a4\\t3\\thyroid\\p16_genesplit.tsv");
		thyroidPaths.add("H:\\GOBI\\a4\\t3\\thyroid\\p17_genesplit.tsv");
		thyroidPaths.add("H:\\GOBI\\a4\\t3\\thyroid\\p18_genesplit.tsv");
		thyroidPaths.add("H:\\GOBI\\a4\\t3\\thyroid\\p19_genesplit.tsv");
		thyroidPaths.add("H:\\GOBI\\a4\\t3\\thyroid\\p20_genesplit.tsv");
		thyroidPaths.add("H:\\GOBI\\a4\\t3\\thyroid\\p21_genesplit.tsv");
		thyroidPaths.add("H:\\GOBI\\a4\\t3\\thyroid\\p22_genesplit.tsv");
		thyroidPaths.add("H:\\GOBI\\a4\\t3\\thyroid\\p23_genesplit.tsv");
		thyroidPaths.add("H:\\GOBI\\a4\\t3\\thyroid\\p24_genesplit.tsv");
		thyroidPaths.add("H:\\GOBI\\a4\\t3\\thyroid\\p25_genesplit.tsv");
		thyroidPaths.add("H:\\GOBI\\a4\\t3\\thyroid\\p26_genesplit.tsv");
		thyroidPaths.add("H:\\GOBI\\a4\\t3\\thyroid\\p27_genesplit.tsv");
		thyroidPaths.add("H:\\GOBI\\a4\\t3\\thyroid\\p28_genesplit.tsv");
		thyroidPaths.add("H:\\GOBI\\a4\\t3\\thyroid\\p29_genesplit.tsv");
		thyroidPaths.add("H:\\GOBI\\a4\\t3\\thyroid\\p30_genesplit.tsv");
		thyroidPaths.add("H:\\GOBI\\a4\\t3\\thyroid\\p31_genesplit.tsv");
		thyroidPaths.add("H:\\GOBI\\a4\\t3\\thyroid\\p32_genesplit.tsv");
		thyroidPaths.add("H:\\GOBI\\a4\\t3\\thyroid\\p33_genesplit.tsv");
		thyroidPaths.add("H:\\GOBI\\a4\\t3\\thyroid\\p34_genesplit.tsv");
		thyroidPaths.add("H:\\GOBI\\a4\\t3\\thyroid\\p35_genesplit.tsv");
		thyroidPaths.add("H:\\GOBI\\a4\\t3\\thyroid\\p36_genesplit.tsv");
		thyroidPaths.add("H:\\GOBI\\a4\\t3\\thyroid\\p37_genesplit.tsv");
		thyroidPaths.add("H:\\GOBI\\a4\\t3\\thyroid\\p38_genesplit.tsv");
		thyroidPaths.add("H:\\GOBI\\a4\\t3\\thyroid\\p39_genesplit.tsv");
		thyroidPaths.add("H:\\GOBI\\a4\\t3\\thyroid\\p40_genesplit.tsv");
		thyroidPaths.add("H:\\GOBI\\a4\\t3\\thyroid\\p41_genesplit.tsv");
		thyroidPaths.add("H:\\GOBI\\a4\\t3\\thyroid\\p42_genesplit.tsv");
		thyroidPaths.add("H:\\GOBI\\a4\\t3\\thyroid\\p43_genesplit.tsv");
		thyroidPaths.add("H:\\GOBI\\a4\\t3\\thyroid\\p44_genesplit.tsv");
		thyroidPaths.add("H:\\GOBI\\a4\\t3\\thyroid\\p45_genesplit.tsv");
		thyroidPaths.add("H:\\GOBI\\a4\\t3\\thyroid\\p46_genesplit.tsv");
		thyroidPaths.add("H:\\GOBI\\a4\\t3\\thyroid\\p47_genesplit.tsv");
		thyroidPaths.add("H:\\GOBI\\a4\\t3\\thyroid\\p48_genesplit.tsv");
		thyroidPaths.add("H:\\GOBI\\a4\\t3\\thyroid\\p49_genesplit.tsv");
		thyroidPaths.add("H:\\GOBI\\a4\\t3\\thyroid\\p50_genesplit.tsv");
		thyroidPaths.add("H:\\GOBI\\a4\\t3\\thyroid\\p51_genesplit.tsv");
		thyroidPaths.add("H:\\GOBI\\a4\\t3\\thyroid\\p52_genesplit.tsv");
		thyroidPaths.add("H:\\GOBI\\a4\\t3\\thyroid\\p53_genesplit.tsv");
		thyroidPaths.add("H:\\GOBI\\a4\\t3\\thyroid\\p54_genesplit.tsv");
		thyroidPaths.add("H:\\GOBI\\a4\\t3\\thyroid\\p55_genesplit.tsv");
		thyroidPaths.add("H:\\GOBI\\a4\\t3\\thyroid\\p56_genesplit.tsv");
		thyroidPaths.add("H:\\GOBI\\a4\\t3\\thyroid\\p57_genesplit.tsv");
		thyroidPaths.add("H:\\GOBI\\a4\\t3\\thyroid\\p58_genesplit.tsv");
		thyroidPaths.add("H:\\GOBI\\a4\\t3\\thyroid\\p59_genesplit.tsv");
		
		GenesplitTSVReader r = new GenesplitTSVReader();
		
		r.readFiles(kidneyPaths, thyroidPaths);
		r.applySteps_mode_1(kidneyPaths, thyroidPaths, 1);
//		r.writeToFile(new ConfigHelper().getDefaultOutputPath()+"consistent_introns.txt");
		
	}
	
}
