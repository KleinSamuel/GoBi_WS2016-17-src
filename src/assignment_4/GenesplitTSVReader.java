package assignment_4;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;

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
				
				tmpList.add(new RelevantIntron(geneID, start, stop, isOnNegativeStrand));
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
	
	public static void main(String[] args) {
		
		GenesplitTSVReader r = new GenesplitTSVReader();
		r.read("H:\\GOBI\\a4\\t3\\kidney\\p1_genesplit.tsv");
		r.read("H:\\GOBI\\a4\\t3\\kidney\\p2_genesplit.tsv");
		r.read("H:\\GOBI\\a4\\t3\\kidney\\p3_genesplit.tsv");
		r.read("H:\\GOBI\\a4\\t3\\kidney\\p4_genesplit.tsv");
		r.read("H:\\GOBI\\a4\\t3\\kidney\\p5_genesplit.tsv");
		r.read("H:\\GOBI\\a4\\t3\\kidney\\p6_genesplit.tsv");
		r.read("H:\\GOBI\\a4\\t3\\kidney\\p7_genesplit.tsv");
		r.read("H:\\GOBI\\a4\\t3\\kidney\\p8_genesplit.tsv");
		r.read("H:\\GOBI\\a4\\t3\\kidney\\p9_genesplit.tsv");
		r.read("H:\\GOBI\\a4\\t3\\kidney\\p10_genesplit.tsv");
		r.read("H:\\GOBI\\a4\\t3\\kidney\\p11_genesplit.tsv");
		r.read("H:\\GOBI\\a4\\t3\\kidney\\p12_genesplit.tsv");
		r.read("H:\\GOBI\\a4\\t3\\kidney\\p13_genesplit.tsv");
		r.read("H:\\GOBI\\a4\\t3\\kidney\\p14_genesplit.tsv");
		r.read("H:\\GOBI\\a4\\t3\\kidney\\p15_genesplit.tsv");
		r.read("H:\\GOBI\\a4\\t3\\kidney\\p16_genesplit.tsv");
		r.read("H:\\GOBI\\a4\\t3\\kidney\\p17_genesplit.tsv");
		r.read("H:\\GOBI\\a4\\t3\\kidney\\p18_genesplit.tsv");
		r.read("H:\\GOBI\\a4\\t3\\kidney\\p19_genesplit.tsv");
		r.read("H:\\GOBI\\a4\\t3\\kidney\\p20_genesplit.tsv");
		r.read("H:\\GOBI\\a4\\t3\\kidney\\p21_genesplit.tsv");
		r.read("H:\\GOBI\\a4\\t3\\kidney\\p22_genesplit.tsv");
		r.read("H:\\GOBI\\a4\\t3\\kidney\\p23_genesplit.tsv");
		r.read("H:\\GOBI\\a4\\t3\\kidney\\p24_genesplit.tsv");
		r.read("H:\\GOBI\\a4\\t3\\kidney\\p25_genesplit.tsv");
		r.writeToFile("H:\\GOBI\\a4\\t3\\kidney\\consistent_introns.txt");
		
	}
	
}
