package assignment_4;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;

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
			System.out.println(geneIds[i]);
			System.out.println(i+" : "+smallest[i]);
		}
		
	}
	
	public Pair<Integer, Integer> getRegion(String geneID, String filepath){
		
		try {
			ProcessBuilder pb = new ProcessBuilder("grep -n \".*gene.*gene_id \""+geneID+"\".*\" Homo_sapiens.GRCh37.75.gtf | awk -F'\t' '{ print $1\"-\"$4\"-\"$5}'");
			pb.inheritIO();
			pb.directory(new File(filepath));
			pb.start();
		} catch (IOException e) {
			e.printStackTrace();
		}
		
		return new Pair<Integer, Integer>(0, 0);
	}
	
	public static void main(String[] args) {
		
		DataRetrieverTask6 dr = new DataRetrieverTask6();
		dr.first("H:/GOBI/a4/t6/contextmap_counts.tsv");
		
	}
}
